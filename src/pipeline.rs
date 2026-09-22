//! Producer–consumer pipeline: per-contig phasing feeds per-contig vcf2fasta.

use anyhow::{anyhow, Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use crate::chunk::plan_tiles;
use crate::cli::Args;
use crate::logger::Logger;
use crate::output_manager::{OutputManager, OutputMode};
use crate::phasing::{open_maybe_gzip, phase_contig};
use crate::resource::{ResourceBudget, RunMetrics};
use crate::scheduler::{limits, DeviceChoice, ExecutionPlan};
use crate::workunit::{WorkState, WorkUnit};

// ---------------------------------------------------------------------------
// Tuning
// ---------------------------------------------------------------------------

const IN_FLIGHT_MEMORY_FRACTION: f64 = 0.50;
const MAX_IN_FLIGHT_BYTES_HARD_CAP: u64 = 8 * 1024 * 1024 * 1024;
const BATCH_PER_WORKER_MULT: usize = 2;
const MAX_TILES_PER_BATCH: usize = 256;

// ---------------------------------------------------------------------------
// Public outcome
// ---------------------------------------------------------------------------

pub struct PipelineOutcome {
    pub contigs_completed: usize,
    pub contigs_failed: usize,
    pub total_seen: usize,
    pub total_applied: usize,
    pub total_files: usize,
    pub total_warnings: usize,
    pub warnings_by_reason: BTreeMap<&'static str, usize>,
    pub failures: Vec<(String, String)>,
}

struct ContigSummary {
    seen: usize,
    applied: usize,
    files: usize,
    warnings_total: usize,
    warnings_by_reason: BTreeMap<&'static str, usize>,
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

pub fn run_pipeline(
    work_units: Vec<WorkUnit>,
    plan: &ExecutionPlan,
    budget: &ResourceBudget,
    args: &Args,
    logger: Arc<Mutex<Logger>>,
) -> Result<PipelineOutcome> {
    let ready_depth = args.ready_queue_depth.max(1);
    let (ready_tx, ready_rx) = bounded::<WorkUnit>(ready_depth);

    let contig_order: Vec<String> = work_units.iter().map(|u| u.contig.clone()).collect();

    let producer_logger = logger.clone();
    let producer_plan = plan.clone();
    let producer_args = args.clone();
    let producer_budget = budget.clone();
    let producer_handle = thread::Builder::new()
        .name("phasing-producer".into())
        .spawn(move || -> Result<()> {
            run_producer(
                work_units,
                ready_tx,
                &producer_args,
                &producer_plan,
                &producer_budget,
                producer_logger,
            )
        })
        .context("could not spawn phasing producer thread")?;

    let effective_workers = match plan.device {
        DeviceChoice::Cpu => args.vcf2fasta_workers.max(1),
        DeviceChoice::Gpu => 1,
    };

    let metrics = RunMetrics::new();
    let failure_vec = Arc::new(Mutex::new(Vec::<(String, String)>::new()));
    let completed_count = Arc::new(AtomicUsize::new(0));
    let failed_count = Arc::new(AtomicUsize::new(0));
    let seen_total = Arc::new(AtomicUsize::new(0));
    let applied_total = Arc::new(AtomicUsize::new(0));
    let files_total = Arc::new(AtomicUsize::new(0));
    let warnings_total_shared = Arc::new(AtomicUsize::new(0));
    let warnings_by_reason_shared: Arc<Mutex<BTreeMap<&'static str, usize>>> =
        Arc::new(Mutex::new(BTreeMap::new()));

    let shared_writers = Arc::new(Mutex::new(OutputManager::new(
        if args.merged_output {
            OutputMode::Merged
        } else {
            OutputMode::PerContig
        },
        args.prefix.clone(),
        args.line_width,
        contig_order,
    )?));

    let mut consumer_handles = Vec::with_capacity(effective_workers);
    for worker_id in 0..effective_workers {
        let rx = ready_rx.clone();
        let plan = plan.clone();
        let args = args.clone();
        let logger = logger.clone();
        let metrics = metrics.clone();
        let failures = failure_vec.clone();
        let completed = completed_count.clone();
        let failed = failed_count.clone();
        let seen_acc = seen_total.clone();
        let applied_acc = applied_total.clone();
        let files_acc = files_total.clone();
        let warn_total = warnings_total_shared.clone();
        let warn_reasons = warnings_by_reason_shared.clone();
        let writers = shared_writers.clone();
        let budget = budget.clone();
        let h = thread::Builder::new()
            .name(format!("vcf2fasta-consumer-{}", worker_id))
            .spawn(move || -> Result<()> {
                run_consumer(
                    rx, &plan, &args, &budget, logger, metrics, failures, completed, failed,
                    seen_acc, applied_acc, files_acc, warn_total, warn_reasons, writers,
                )
            })
            .context("could not spawn vcf2fasta consumer thread")?;
        consumer_handles.push(h);
    }
    drop(ready_rx);

    let producer_result = producer_handle
        .join()
        .map_err(|_| anyhow!("phasing producer thread panicked"))?;

    let mut consumer_errors = Vec::new();
    for h in consumer_handles {
        match h.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => consumer_errors.push(format!("{:#}", e)),
            Err(_) => consumer_errors.push("consumer thread panicked".to_string()),
        }
    }

    let (sample_names, _) = crate::vcf::read_header_metadata(&args.input)?;
    let manager_files = {
        let mgr = Arc::try_unwrap(shared_writers)
            .map_err(|_| anyhow!("output manager still shared at end of run"))?;
        let mgr = mgr
            .into_inner()
            .map_err(|_| anyhow!("output manager mutex poisoned"))?;
        mgr.finish(&sample_names)?
    };

    let total_files = if args.merged_output {
        manager_files
    } else {
        files_total.load(Ordering::Relaxed)
    };

    if !consumer_errors.is_empty() {
        return Err(anyhow!("consumer errors:\n{}", consumer_errors.join("\n")));
    }
    producer_result?;

    let _m = metrics.lock().unwrap().clone();

    let warnings_by_reason = warnings_by_reason_shared
        .lock()
        .map(|g| g.clone())
        .unwrap_or_default();

    Ok(PipelineOutcome {
        contigs_completed: completed_count.load(Ordering::Relaxed),
        contigs_failed: failed_count.load(Ordering::Relaxed),
        total_seen: seen_total.load(Ordering::Relaxed),
        total_applied: applied_total.load(Ordering::Relaxed),
        total_files,
        total_warnings: warnings_total_shared.load(Ordering::Relaxed),
        warnings_by_reason,
        failures: failure_vec.lock().unwrap().clone(),
    })
}

// ---------------------------------------------------------------------------
// Producer
// ---------------------------------------------------------------------------

fn run_producer(
    mut work_units: Vec<WorkUnit>,
    ready_tx: Sender<WorkUnit>,
    args: &Args,
    _plan: &ExecutionPlan,
    budget: &ResourceBudget,
    logger: Arc<Mutex<Logger>>,
) -> Result<()> {
    let phasing_threads = budget.phasing_threads.max(1);

    for mut unit in work_units.drain(..) {
        unit.state = WorkState::Analyzed;
        let t0 = Instant::now();

        if !unit.needs_phasing {
            unit.state = WorkState::Ready;
            {
                let mut lg = logger.lock().unwrap();
                lg.ready(&unit.contig)?;
            }
            if ready_tx.send(unit).is_err() {
                return Ok(());
            }
            continue;
        }

        unit.state = WorkState::NeedsPhasing;
        if unit.phase_backend.is_subprocess() {
            let mut lg = logger.lock().unwrap();
            lg.phasing_start(&unit.contig, unit.phase_backend.tag(), phasing_threads)?;
        }

        let slice_dir = tempfile::Builder::new()
            .prefix(&format!("vcf2fasta_slice_{}_", sanitize(&unit.contig)))
            .tempdir()
            .context("could not create slicing tempdir")?;
        let slice_path = slice_dir.path().join("slice.vcf");
        if let Err(e) = slice_contig(args.input.as_path(), &unit.contig, &slice_path) {
            unit.state = WorkState::PhasingFailed;
            {
                let mut lg = logger.lock().unwrap();
                lg.failure("phasing", &unit.contig, &format!("slice failed: {:#}", e))?;
            }
            let _ = ready_tx.send(unit);
            continue;
        }

        unit.state = WorkState::Phasing;
        let phased = match phase_contig(
            &slice_path,
            &unit,
            &args.reference,
            args,
            phasing_threads,
            budget.available_ram_bytes,
            &logger,
        ) {
            Ok(p) => p,
            Err(e) => {
                unit.state = WorkState::PhasingFailed;
                {
                    let mut lg = logger.lock().unwrap();
                    lg.failure("phasing", &unit.contig, &format!("{:#}", e))?;
                }
                let _ = ready_tx.send(unit);
                continue;
            }
        };

        drop(slice_dir);
        unit.state = WorkState::PhaseValidation;
        if let Some(p) = phased.path {
            unit.phased_vcf = Some(p);
        }
        if let Some(td) = phased._temp_dir {
            unit.phased_tempdir = Some(Arc::new(td));
        }
        unit.state = WorkState::Ready;
        {
            let mut lg = logger.lock().unwrap();
            if unit.phase_backend.is_subprocess() {
                lg.phasing_done(&unit.contig, t0.elapsed())?;
            }
            lg.ready(&unit.contig)?;
        }
        budget.sample_rss();
        if ready_tx.send(unit).is_err() {
            return Ok(());
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Consumer
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn run_consumer(
    ready_rx: Receiver<WorkUnit>,
    plan: &ExecutionPlan,
    args: &Args,
    budget: &ResourceBudget,
    logger: Arc<Mutex<Logger>>,
    metrics: Arc<Mutex<RunMetrics>>,
    failures: Arc<Mutex<Vec<(String, String)>>>,
    completed: Arc<AtomicUsize>,
    failed: Arc<AtomicUsize>,
    seen_total: Arc<AtomicUsize>,
    applied_total: Arc<AtomicUsize>,
    files_total: Arc<AtomicUsize>,
    warnings_total: Arc<AtomicUsize>,
    warnings_by_reason: Arc<Mutex<BTreeMap<&'static str, usize>>>,
    writers: Arc<Mutex<OutputManager>>,
) -> Result<()> {
    while let Ok(unit) = ready_rx.recv() {
        if unit.state.is_terminal_failure() {
            failed.fetch_add(1, Ordering::Relaxed);
            failures
                .lock()
                .unwrap()
                .push((unit.contig.clone(), format!("terminal: {}", unit.state.tag())));
            continue;
        }

        let t0 = Instant::now();
        match process_one_contig(&unit, plan, args, budget, &logger, &writers) {
            Ok(summary) => {
                let elapsed = t0.elapsed();
                completed.fetch_add(1, Ordering::Relaxed);
                seen_total.fetch_add(summary.seen, Ordering::Relaxed);
                applied_total.fetch_add(summary.applied, Ordering::Relaxed);
                files_total.fetch_add(summary.files, Ordering::Relaxed);
                warnings_total.fetch_add(summary.warnings_total, Ordering::Relaxed);
                if !summary.warnings_by_reason.is_empty() {
                    let mut shared = warnings_by_reason.lock().unwrap();
                    for (reason, count) in &summary.warnings_by_reason {
                        *shared.entry(reason).or_insert(0) += count;
                    }
                }
                metrics.lock().unwrap().vcf2fasta_total += elapsed;
                {
                    let mut lg = logger.lock().unwrap();
                    lg.complete(
                        &unit.contig,
                        summary.seen,
                        summary.applied,
                        summary.files,
                        summary.warnings_total,
                    )?;
                    lg.flush_warnings()?;
                }
            }
            Err(e) => {
                failed.fetch_add(1, Ordering::Relaxed);
                failures
                    .lock()
                    .unwrap()
                    .push((unit.contig.clone(), format!("{:#}", e)));
                {
                    let mut lg = logger.lock().unwrap();
                    lg.failure("vcf2fasta", &unit.contig, &format!("{:#}", e))?;
                }
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Per-contig processor
// ---------------------------------------------------------------------------

fn process_one_contig(
    unit: &WorkUnit,
    plan: &ExecutionPlan,
    args: &Args,
    budget: &ResourceBudget,
    logger: &Arc<Mutex<Logger>>,
    writers: &Arc<Mutex<OutputManager>>,
) -> Result<ContigSummary> {
    let mut local_args = args.clone();
    local_args.input = unit.effective_vcf(&args.input).to_path_buf();

    let contig_lengths: HashMap<String, u64> =
        HashMap::from([(unit.contig.clone(), unit.reference_length)]);
    let contigs = vec![unit.contig.clone()];

    // Always respect MAX_TILE_OUTPUT_BYTES. The previous `.max(20_000)`
    // floor silently exceeded the memory cap for cohorts above ~838
    // haplotypes; removing it lets the scheduler keep per-tile output
    // bounded at every scale, at the cost of more (smaller) tiles.
    let haps = unit.haplotype_count.max(1) as u64;
    let bounded_base_block = (limits::MAX_TILE_OUTPUT_BYTES / haps).max(1);
    let base_block_size = match args.chunk_size {
        Some(cs) => cs.max(1),
        None => bounded_base_block,
    };

    let tiles = plan_tiles(
        &contigs,
        &contig_lengths,
        unit.haplotype_count,
        plan.haplotype_block_size,
        base_block_size,
    );
    if tiles.is_empty() {
        return Ok(ContigSummary {
            seen: 0,
            applied: 0,
            files: 0,
            warnings_total: 0,
            warnings_by_reason: BTreeMap::new(),
        });
    }

    {
        let mut lg = logger.lock().unwrap();
        lg.scheduler_dispatch(&unit.contig, plan.device.as_str(), tiles.len())?;
    }

    let (sample_names, _) = crate::vcf::read_header_metadata(&local_args.input)?;

    let mut sample_hap_offset = Vec::with_capacity(sample_names.len());
    let mut acc = 0usize;
    for &p in &unit.sample_max_ploidies {
        sample_hap_offset.push(acc);
        acc += p;
    }
    let max_ploidies = unit.sample_max_ploidies.clone();

    let chunk_pad = args
        .chunk_pad
        .unwrap_or_else(|| crate::chunk::auto_chunk_pad(&local_args.input).unwrap_or(10_000));

    let worker_threads = budget.vcf2fasta_threads.max(1);

    let mut seen_total = 0usize;
    let mut applied_total = 0usize;
    let mut warnings_total = 0usize;
    let mut warnings_by_reason: BTreeMap<&'static str, usize> = BTreeMap::new();

    match plan.device {
        DeviceChoice::Cpu => {
            let hap_block = plan.haplotype_block_size.max(1) as u64;
            let per_tile_output_bytes = hap_block.saturating_mul(base_block_size).max(1);
            let memory_ceiling =
                ((budget.max_memory_bytes as f64) * IN_FLIGHT_MEMORY_FRACTION) as u64;
            let memory_ceiling = memory_ceiling.min(MAX_IN_FLIGHT_BYTES_HARD_CAP);
            let batch_by_memory = (memory_ceiling / per_tile_output_bytes).max(1) as usize;
            let batch_by_throughput =
                worker_threads.saturating_mul(BATCH_PER_WORKER_MULT).max(1);
            let batch_size = batch_by_memory
                .min(batch_by_throughput)
                .min(MAX_TILES_PER_BATCH)
                .max(1);

            {
                let mut lg = logger.lock().unwrap();
                lg.raw(&format!(
                    "[SCHEDULER] {}: tiles={} base_block={} per_tile_out≈{} KiB \
                     mem_ceiling={} MiB batch={} (mem={} thr={})",
                    unit.contig,
                    tiles.len(),
                    base_block_size,
                    per_tile_output_bytes / 1024,
                    memory_ceiling / (1024 * 1024),
                    batch_size,
                    batch_by_memory,
                    batch_by_throughput,
                ))?;
            }

            let total_batches = tiles.len().div_ceil(batch_size);
            for (batch_idx, batch) in tiles.chunks(batch_size).enumerate() {
                let mut results = crate::run_cpu_tiles_public(
                    batch,
                    &sample_names,
                    &sample_hap_offset,
                    &max_ploidies,
                    chunk_pad,
                    &local_args,
                    worker_threads,
                );
                results.sort_by_key(|r| match r {
                    Ok(tr) => (tr.base_start, tr.hap_start),
                    Err(_) => (u64::MAX, usize::MAX),
                });
                for r in results {
                    let tr = r?;
                    write_and_count(
                        tr,
                        &sample_names,
                        &sample_hap_offset,
                        &max_ploidies,
                        writers,
                        logger,
                        &mut seen_total,
                        &mut applied_total,
                        &mut warnings_total,
                        &mut warnings_by_reason,
                    )?;
                }
                budget.sample_rss();

                if batch_idx > 0 && (batch_idx % (total_batches / 10).max(1) == 0) {
                    let mut lg = logger.lock().unwrap();
                    lg.raw(&format!(
                        "[SCHEDULER] {}: batch {}/{} done (peak RSS {} MiB)",
                        unit.contig,
                        batch_idx + 1,
                        total_batches,
                        budget.peak_rss_bytes() / (1024 * 1024),
                    ))?;
                }
            }
        }
        DeviceChoice::Gpu => {
            {
                let mut lg = logger.lock().unwrap();
                lg.raw(&format!(
                    "[SCHEDULER] {}: tiles={} base_block={} (streaming via CUDA on {} device(s))",
                    unit.contig,
                    tiles.len(),
                    base_block_size,
                    plan.gpu_device_indices.len().max(1),
                ))?;
            }

            crate::run_gpu_tiles(
                &tiles,
                &sample_names,
                &sample_hap_offset,
                &max_ploidies,
                chunk_pad,
                &local_args,
                plan.gpu_streams,
                plan.in_flight_buffers,
                &plan.gpu_device_indices,
                budget,
                logger,
                |tr| {
                    write_and_count(
                        tr,
                        &sample_names,
                        &sample_hap_offset,
                        &max_ploidies,
                        writers,
                        logger,
                        &mut seen_total,
                        &mut applied_total,
                        &mut warnings_total,
                        &mut warnings_by_reason,
                    )
                },
            )?;
            budget.sample_rss();
        }
    }

    {
        let mut lg = logger.lock().unwrap();
        lg.warning_summary(&unit.contig, warnings_total, &warnings_by_reason)?;
        lg.flush_warnings()?;
    }

    {
        let mut mgr = writers.lock().unwrap();
        mgr.finish_contig(&unit.contig)?;
    }

    let contig_file_count = if args.merged_output {
        0
    } else {
        unit.haplotype_count
    };

    Ok(ContigSummary {
        seen: seen_total,
        applied: applied_total,
        files: contig_file_count,
        warnings_total,
        warnings_by_reason,
    })
}

#[allow(clippy::too_many_arguments)]
fn write_and_count(
    tr: crate::chunk::TileResult,
    sample_names: &[String],
    sample_hap_offset: &[usize],
    max_ploidies: &[usize],
    writers: &Arc<Mutex<OutputManager>>,
    logger: &Arc<Mutex<Logger>>,
    seen_total: &mut usize,
    applied_total: &mut usize,
    warnings_total: &mut usize,
    warnings_by_reason: &mut BTreeMap<&'static str, usize>,
) -> Result<()> {
    *seen_total += tr.seen;
    *applied_total += tr.applied;
    *warnings_total += tr.warning_count;

    for (reason, count) in &tr.warnings_by_reason {
        *warnings_by_reason.entry(reason).or_insert(0) += count;
    }

    if !tr.warnings.is_empty() {
        let mut lg = logger.lock().unwrap();
        for w in &tr.warnings {
            lg.warn(w)?;
        }
    }

    {
        let mut mgr = writers.lock().unwrap();
        for (local_h, bytes) in tr.per_hap.iter().enumerate() {
            let global_h = tr.hap_start + local_h;
            let Some((s_idx, l_h)) =
                global_to_sample(global_h, sample_hap_offset, max_ploidies)
            else {
                continue;
            };
            let sample_name = sample_names.get(s_idx).cloned().unwrap_or_default();
            mgr.write_segment(&tr.contig, s_idx, &sample_name, l_h, bytes)?;
        }
    }
    Ok(())
}

fn global_to_sample(
    global_h: usize,
    offsets: &[usize],
    ploidies: &[usize],
) -> Option<(usize, usize)> {
    for (s_idx, (&off, &p)) in offsets.iter().zip(ploidies.iter()).enumerate() {
        if global_h >= off && global_h < off + p {
            return Some((s_idx, global_h - off));
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Contig slicing
// ---------------------------------------------------------------------------

fn slice_contig(input: &Path, contig: &str, output: &Path) -> Result<()> {
    use std::io::Write as _;
    let reader = open_maybe_gzip(input)?;
    let mut out = File::create(output)
        .with_context(|| format!("could not create slice {}", output.display()))?;
    for line in std::io::BufRead::lines(reader) {
        let line = line?;
        if line.starts_with('#') {
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
            continue;
        }
        let chrom = match line.split('\t').next() {
            Some(s) => s,
            None => continue,
        };
        if chrom == contig {
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn slice_contig_keeps_header_and_matching_records() {
        let dir = tempdir().unwrap();
        let inp = dir.path().join("in.vcf");
        let out = dir.path().join("chr1.vcf");
        let mut f = File::create(&inp).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "##contig=<ID=chr1,length=100>").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1").unwrap();
        writeln!(f, "chr1\t1\t.\tA\tG\t.\t.\t.\tGT\t0|1").unwrap();
        writeln!(f, "chr2\t1\t.\tA\tG\t.\t.\t.\tGT\t0|1").unwrap();
        writeln!(f, "chr1\t2\t.\tC\tT\t.\t.\t.\tGT\t1|0").unwrap();
        drop(f);

        slice_contig(&inp, "chr1", &out).unwrap();
        let content = std::fs::read_to_string(&out).unwrap();
        assert!(content.contains("chr1\t1\t"));
        assert!(content.contains("chr1\t2\t"));
        assert!(!content.contains("chr2\t1\t"));
        assert!(content.starts_with("##fileformat"));
    }

    #[test]
    fn batch_size_keeps_peak_bounded_for_large_cohorts() {
        let haps: u64 = 300;
        let base_block = (limits::MAX_TILE_OUTPUT_BYTES / haps).max(1);
        let per_tile = haps * base_block;
        assert!(per_tile <= limits::MAX_TILE_OUTPUT_BYTES * 2);

        let memory_ceiling = MAX_IN_FLIGHT_BYTES_HARD_CAP;
        let batch_by_mem = (memory_ceiling / per_tile).max(1) as usize;
        let peak = batch_by_mem as u64 * per_tile;
        assert!(peak <= MAX_IN_FLIGHT_BYTES_HARD_CAP * 2);
    }

    #[test]
    fn tile_output_is_capped_at_every_cohort_size() {
        for &haps in &[100u64, 400, 838, 2_000, 10_000, 100_000] {
            let base_block = (limits::MAX_TILE_OUTPUT_BYTES / haps).max(1);
            let per_tile = haps.saturating_mul(base_block);
            assert!(
                per_tile <= limits::MAX_TILE_OUTPUT_BYTES * 2,
                "haps={} per_tile={} exceeds cap",
                haps,
                per_tile
            );
        }
    }

    #[test]
    fn batch_size_never_zero_or_huge() {
        let per_tile: u64 = 1;
        let batch_by_mem = (MAX_IN_FLIGHT_BYTES_HARD_CAP / per_tile).max(1) as usize;
        let final_batch = batch_by_mem.min(MAX_TILES_PER_BATCH);
        assert!(final_batch <= MAX_TILES_PER_BATCH);

        let per_tile = MAX_IN_FLIGHT_BYTES_HARD_CAP * 4;
        let batch_by_mem = (MAX_IN_FLIGHT_BYTES_HARD_CAP / per_tile).max(1) as usize;
        assert_eq!(batch_by_mem, 1);
    }
}