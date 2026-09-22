//! Main library entry point.

pub mod bam_resolver;
pub mod beagle_resources;
pub mod chunk;
pub mod cli;
pub mod fasta;
pub mod genotype;
pub mod logger;
pub mod output_manager;
pub mod phasing;
pub mod pinned;
pub mod pipeline;
pub mod report;
pub mod resource;
pub mod scheduler;
pub mod tile_executor;
pub mod util;
pub mod variable_ploidy;
pub mod vcf;
pub mod workload;
pub mod workunit;

#[cfg(feature = "cuda")]
pub mod cuda_pipeline;
#[cfg(feature = "cuda")]
pub mod gpu;

use anyhow::{bail, Context, Result};
use rust_htslib::faidx;
use rust_htslib::tbx::Reader as TbxReader;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::{Duration, Instant, SystemTime};

use crate::cli::{
    effective_device, parse_gpu_list, parse_size, resolve_phasing_threads, resolve_threads,
    validate_cli, Args,
};
use crate::logger::Logger;
use crate::resource::ResourceBudget;
use crate::scheduler::{DeviceChoice, ExecutionPlan, HardwareInfo, Scheduler, WorkloadProfile};
use crate::vcf::read_header_metadata;

// ---------------------------------------------------------------------------
// Process-wide CUDA device cache
// ---------------------------------------------------------------------------
//
// Creating a `CudaDevice` compiles the kernel with NVRTC (~100 ms) and
// creates a host-side CUDA context. We cache one per device index, so a
// `--gpu-devices 0,1` run pays the cost once per device rather than once
// per contig or once per worker stream.

#[cfg(feature = "cuda")]
type CudaDeviceCache =
    Mutex<HashMap<usize, std::result::Result<Arc<crate::gpu::CudaDevice>, String>>>;

#[cfg(feature = "cuda")]
static CUDA_DEVICE_CACHE: OnceLock<CudaDeviceCache> = OnceLock::new();

#[cfg(feature = "cuda")]
fn get_cuda_device(device_index: usize) -> Result<Arc<crate::gpu::CudaDevice>> {
    let cache = CUDA_DEVICE_CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    // Hold the lock across initialization so two threads asking for the
    // same device do not both pay NVRTC. Init is ~100 ms and the number
    // of devices is small, so the serialization is acceptable.
    let mut guard = cache
        .lock()
        .map_err(|_| anyhow::anyhow!("CUDA device cache poisoned"))?;

    if let Some(entry) = guard.get(&device_index) {
        return match entry {
            Ok(dev) => Ok(dev.clone()),
            Err(msg) => anyhow::bail!(
                "CUDA device {} initialization failed: {}",
                device_index,
                msg
            ),
        };
    }

    let result = crate::gpu::CudaDevice::init_device(device_index)
        .map(Arc::new)
        .map_err(|e| format!("{:#}", e));

    match result {
        Ok(dev) => {
            guard.insert(device_index, Ok(dev.clone()));
            Ok(dev)
        }
        Err(msg) => {
            let cached = msg.clone();
            guard.insert(device_index, Err(cached));
            anyhow::bail!(
                "CUDA device {} initialization failed: {}",
                device_index,
                msg
            )
        }
    }
}

// ---------------------------------------------------------------------------
// Memory telemetry
// ---------------------------------------------------------------------------
//
// `read_rss_mib` reports this process's own RSS, which is what the tool's
// internal bounds actually govern.
//
// `read_system_memory_stats` also reports the kernel's page cache and
// dirty-page counts. These are *not* process memory: they are the OS
// buffering recently-written disk pages and are reclaimed on demand. They
// often reach 10 GB on write-heavy workloads like FASTA generation and
// are the usual reason a naive `free -h` reading looks alarming.

fn read_rss_mib() -> u64 {
    #[cfg(target_os = "linux")]
    {
        if let Ok(s) = std::fs::read_to_string("/proc/self/statm") {
            if let Some(pages_str) = s.split_whitespace().nth(1) {
                if let Ok(pages) = pages_str.parse::<u64>() {
                    return pages.saturating_mul(4096) / (1024 * 1024);
                }
            }
        }
    }
    0
}

/// Returns `(process_rss_mib, kernel_cached_mib, kernel_dirty_mib)`.
fn read_system_memory_stats() -> (u64, u64, u64) {
    let rss = read_rss_mib();
    let mut cached = 0u64;
    let mut dirty = 0u64;
    #[cfg(target_os = "linux")]
    {
        if let Ok(s) = std::fs::read_to_string("/proc/meminfo") {
            for line in s.lines() {
                if let Some(rest) = line.strip_prefix("Cached:") {
                    if let Some(v) = rest.split_whitespace().next() {
                        if let Ok(kb) = v.parse::<u64>() {
                            cached = kb / 1024;
                        }
                    }
                } else if let Some(rest) = line.strip_prefix("Dirty:") {
                    if let Some(v) = rest.split_whitespace().next() {
                        if let Ok(kb) = v.parse::<u64>() {
                            dirty = kb / 1024;
                        }
                    }
                }
            }
        }
    }
    (rss, cached, dirty)
}

// ---------------------------------------------------------------------------
// Reference / index helpers
// ---------------------------------------------------------------------------

fn read_contig_lengths(reference: &Path) -> Result<HashMap<String, u64>> {
    let reader = faidx::Reader::from_path(reference)
        .with_context(|| format!("could not open indexed reference {}", reference.display()))?;
    let names = reader.seq_names()?;
    let mut out = HashMap::with_capacity(names.len());
    for name in &names {
        out.insert(name.clone(), reader.fetch_seq_len(name));
    }
    Ok(out)
}

fn tabix_indexed_contigs(path: &Path) -> Result<HashSet<String>> {
    let reader = TbxReader::from_path(path)
        .with_context(|| format!("could not open {} as tabix-indexed VCF", path.display()))?;
    let mut out = HashSet::new();
    for name in reader.seqnames() {
        out.insert(name.to_string());
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Legacy (non-pipelined) path
// ---------------------------------------------------------------------------

mod legacy {
    use super::*;
    use crate::chunk::{plan_tiles, Tile};
    use crate::fasta::FastaWriter;
    use crate::phasing::ensure_phased_vcf;
    use crate::report::MAX_WARNINGS_PER_CONTIG;
    use std::io::Write;

    pub fn run(args: &mut Args, log_file: &mut File) -> Result<()> {
        let threads = resolve_threads(args);
        let phased = ensure_phased_vcf(&args.input, &args.reference, args, threads)
            .with_context(|| "phasing stage failed")?;
        eprintln!("Phasing strategy: {}", phased.strategy.describe());
        writeln!(log_file, "Phasing strategy: {}", phased.strategy.describe())?;
        let _temp_guard = phased.temp_dir;
        args.input = phased.path;

        let (sample_names, contigs) = read_header_metadata(&args.input)?;
        let contig_lengths = read_contig_lengths(&args.reference)?;
        let used_contigs: Vec<String> = contigs
            .iter()
            .filter(|c| contig_lengths.contains_key(*c))
            .cloned()
            .collect();
        if used_contigs.is_empty() {
            bail!("No contig in both VCF header and reference FASTA.");
        }

        let mut max_ploidies = vec![0usize; sample_names.len()];
        {
            let reader = crate::phasing::open_maybe_gzip(&args.input)?;
            for line in std::io::BufRead::lines(reader) {
                let line = line?;
                if line.starts_with('#') {
                    continue;
                }
                crate::vcf::update_max_ploidies_from_line(&line, &sample_names, &mut max_ploidies);
            }
        }
        for p in &mut max_ploidies {
            if *p == 0 {
                *p = 2;
            }
        }
        let haplotype_count: usize = max_ploidies.iter().sum();

        let mut sample_hap_offset = Vec::with_capacity(sample_names.len());
        let mut acc = 0usize;
        for &p in &max_ploidies {
            sample_hap_offset.push(acc);
            acc += p;
        }

        let base_block = args.chunk_size.unwrap_or(1_000_000).max(1);
        let tiles: Vec<Tile> = plan_tiles(
            &used_contigs,
            &contig_lengths,
            haplotype_count,
            haplotype_count.max(1),
            base_block,
        );
        if tiles.is_empty() {
            bail!("no tiles produced");
        }

        let chunk_pad = args.chunk_pad.unwrap_or(10_000);

        let mut results = crate::run_cpu_tiles_public(
            &tiles,
            &sample_names,
            &sample_hap_offset,
            &max_ploidies,
            chunk_pad,
            args,
            threads,
        );

        let contig_order: HashMap<String, usize> = used_contigs
            .iter()
            .enumerate()
            .map(|(i, c)| (c.clone(), i))
            .collect();

        results.sort_by_key(|r| match r {
            Ok(tr) => (
                contig_order.get(&tr.contig).copied().unwrap_or(usize::MAX),
                tr.base_start,
                tr.hap_start,
            ),
            Err(_) => (usize::MAX, u64::MAX, usize::MAX),
        });

        let mut writers: HashMap<(String, usize, usize), FastaWriter> = HashMap::new();
        let mut aggregate_seen = 0usize;
        let mut aggregate_applied = 0usize;
        let mut warnings_emitted: HashMap<String, usize> = HashMap::new();

        for r in results {
            let tr = r?;
            aggregate_seen += tr.seen;
            aggregate_applied += tr.applied;
            let entry = warnings_emitted.entry(tr.contig.clone()).or_insert(0);
            for w in &tr.warnings {
                if *entry < MAX_WARNINGS_PER_CONTIG {
                    writeln!(log_file, "WARNING {}", w)?;
                    *entry += 1;
                }
            }
            for (local_h, bytes) in tr.per_hap.iter().enumerate() {
                let global_h = tr.hap_start + local_h;
                let (s_idx, l_h) =
                    match global_to_sample(global_h, &sample_hap_offset, &max_ploidies) {
                        Some(x) => x,
                        None => continue,
                    };
                let key = (tr.contig.clone(), s_idx, l_h);
                if !writers.contains_key(&key) {
                    let sample = &sample_names[s_idx];
                    let seq_name = format!("{}_{}:{}", sample, tr.contig, l_h);
                    let file_name = format!("{}{}.fa", args.prefix, seq_name);
                    let w = FastaWriter::create(Path::new(&file_name), &seq_name, args.line_width)?;
                    writers.insert(key.clone(), w);
                }
                writers.get_mut(&key).unwrap().write_seq(bytes)?;
            }
        }
        for (_, w) in writers {
            w.finish()?;
        }

        writeln!(
            log_file,
            "summary: seen={}, applied={}",
            aggregate_seen, aggregate_applied
        )?;
        Ok(())
    }

    fn global_to_sample(g: usize, offs: &[usize], pl: &[usize]) -> Option<(usize, usize)> {
        for (s, (&o, &p)) in offs.iter().zip(pl.iter()).enumerate() {
            if g >= o && g < o + p {
                return Some((s, g - o));
            }
        }
        None
    }
}

// ---------------------------------------------------------------------------
// run()
// ---------------------------------------------------------------------------

pub fn run(mut args: Args) -> Result<()> {
    let start_time = SystemTime::now();

    let log_path: PathBuf = if args.prefix.is_empty() {
        PathBuf::from("vcf2fasta.log")
    } else {
        let mut p = PathBuf::from(&args.prefix);
        p.set_extension("log");
        p
    };

    let warnings_path: PathBuf = if args.prefix.is_empty() {
        PathBuf::from("vcf2fasta.warnings.log")
    } else {
        let mut p = PathBuf::from(&args.prefix);
        let stem = p
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("vcf2fasta")
            .to_string();
        p.set_file_name(format!("{}.warnings.log", stem));
        p
    };

    let logger = Arc::new(Mutex::new(Logger::create(
        &log_path,
        Some(&warnings_path),
        args.quiet,
    )?));

    validate_cli(&args)?;
    let threads = resolve_threads(&args);
    args.threads = threads;

    {
        let mut lg = logger.lock().unwrap();
        lg.raw(&format!(
            "=== vcf2fasta run started at {} ===",
            chrono::Local::now().format("%Y-%m-%d %I:%M:%S %p")
        ))?;
        lg.raw(&format!(
            "Command: {}",
            std::env::args().collect::<Vec<_>>().join(" ")
        ))?;
        lg.raw(&format!("Warnings log: {}", warnings_path.display()))?;
        if args.no_call_string.is_none() {
            lg.raw("WARNING: --no-call-string not provided; defaulting to 'N'")?;
        }
    }

    if args.no_pipeline {
        let mut log_file = File::create(&log_path)?;
        return legacy::run(&mut args, &mut log_file);
    }

    let (samples, header_contigs) = read_header_metadata(&args.input)?;
    if samples.is_empty() {
        bail!("input contains no samples");
    }
    if header_contigs.is_empty() {
        bail!("input header contains no contigs");
    }

    let contig_lengths = read_contig_lengths(&args.reference)?;

    let mut used_contigs: Vec<String> = Vec::with_capacity(header_contigs.len());
    for c in &header_contigs {
        if contig_lengths.contains_key(c) {
            used_contigs.push(c.clone());
        } else {
            logger.lock().unwrap().raw(&format!(
                "WARNING: contig '{}' not present in reference; skipping",
                c
            ))?;
        }
    }

    let indexed = tabix_indexed_contigs(&args.input)?;
    let before = used_contigs.len();
    used_contigs.retain(|c| indexed.contains(c));
    let dropped = before - used_contigs.len();
    if dropped > 0 {
        logger.lock().unwrap().raw(&format!(
            "INFO: dropped {} contig(s) present in VCF header but absent from tabix index",
            dropped
        ))?;
    }
    if used_contigs.is_empty() {
        bail!("No contig is simultaneously in VCF header, reference, and tabix index.");
    }

    let discovery = crate::workunit::discover(&args.input, &contig_lengths, &used_contigs, &args)
        .with_context(|| "discovery pass failed")?;

    {
        let mut lg = logger.lock().unwrap();
        for u in &discovery.work_units {
            lg.discovery(
                &u.contig,
                u.sample_count,
                u.haplotype_count,
                u.variant_count,
                u.phased_genotypes,
                u.unphased_genotypes,
                u.phase_backend.tag(),
            )?;
        }
    }

    let reference_bases: u64 = used_contigs
        .iter()
        .filter_map(|c| contig_lengths.get(c).copied())
        .sum();
    let global_haplotypes = discovery
        .work_units
        .iter()
        .map(|u| u.haplotype_count)
        .max()
        .unwrap_or(0);
    let global_max_ploidy = discovery
        .work_units
        .iter()
        .map(|u| u.max_ploidy)
        .max()
        .unwrap_or(0);

    let workload = WorkloadProfile {
        sample_count: samples.len(),
        haplotype_count: global_haplotypes,
        max_ploidy: global_max_ploidy,
        sample_max_ploidies: discovery
            .work_units
            .first()
            .map(|u| u.sample_max_ploidies.clone())
            .unwrap_or_else(|| vec![2; samples.len()]),
        variant_count: discovery.work_units.iter().map(|u| u.variant_count).sum(),
        reference_bases,
        phased_genotypes: discovery.work_units.iter().map(|u| u.phased_genotypes).sum(),
        unphased_genotypes: discovery
            .work_units
            .iter()
            .map(|u| u.unphased_genotypes)
            .sum(),
        malformed_records: discovery
            .work_units
            .iter()
            .map(|u| u.malformed_records)
            .sum(),
        snv_count: 0,
        indel_count: 0,
        contig_count: used_contigs.len(),
        estimated_output_bytes: 0,
    };

    let max_mem = args.max_memory.as_deref().map(parse_size).transpose()?;
    let max_vram = args.max_vram.as_deref().map(parse_size).transpose()?;
    let budget = ResourceBudget::detect(max_mem, max_vram);

    // Log if the user asked for more memory than the machine can provide.
    if let Some(original) = budget.max_memory_clamped_from {
        logger.lock().unwrap().raw(&format!(
            "WARNING: --max-memory={} MiB clamped to {} MiB \
             (safe fraction of {} MiB total machine RAM).",
            original / (1024 * 1024),
            budget.max_memory_bytes / (1024 * 1024),
            budget.total_ram_bytes / (1024 * 1024),
        ))?;
    }

    let hw = HardwareInfo::detect();
    let requested = effective_device(&args);

    {
        let mut lg = logger.lock().unwrap();
        lg.raw(&format!(
            "[SCHEDULER] hardware: cores={} total_ram={} MiB available_ram={} MiB gpus={}",
            hw.cpu_cores,
            hw.total_ram_bytes / (1024 * 1024),
            hw.available_ram_bytes / (1024 * 1024),
            hw.gpus.len(),
        ))?;
        for g in &hw.gpus {
            lg.raw(&format!(
                "[SCHEDULER]   gpu[{}]: {} — {} MiB free / {} MiB total",
                g.device_index,
                g.name,
                g.free_vram_bytes / (1024 * 1024),
                g.total_vram_bytes / (1024 * 1024),
            ))?;
        }
    }

    let plan: ExecutionPlan = {
        let plan_result = if let Some(s) = &args.gpu_devices {
            let indices = parse_gpu_list(s)?;
            Scheduler::plan_multi_gpu(&hw, &workload, indices, budget.max_vram_bytes)
        } else {
            Scheduler::plan(requested, &hw, &workload, budget.max_vram_bytes)
        };

        match plan_result {
            Ok(p) => p,
            Err(e) => {
                {
                    let mut lg = logger.lock().unwrap();
                    lg.raw(&format!("[SCHEDULER] plan failed: {:#}", e))?;

                    let gpu_summary = if hw.gpus.is_empty() {
                        if cfg!(feature = "cuda") {
                            "0 (no CUDA device visible; check `nvidia-smi` and the driver)".to_string()
                        } else {
                            "0 (binary built without `--features cuda`)".to_string()
                        }
                    } else {
                        format!(
                            "{} ({})",
                            hw.gpus.len(),
                            hw.gpus
                                .iter()
                                .map(|g| g.name.as_str())
                                .collect::<Vec<_>>()
                                .join(", ")
                        )
                    };
                    lg.raw(&format!(
                        "[SCHEDULER] requested device={} but detected {} CPU core(s) and {} GPU(s)",
                        match requested {
                            crate::scheduler::Device::Auto => "auto",
                            crate::scheduler::Device::Cpu => "cpu",
                            crate::scheduler::Device::Gpu => "gpu",
                        },
                        hw.cpu_cores,
                        gpu_summary,
                    ))?;
                    lg.flush_warnings()?;
                }
                return Err(e);
            }
        }
    };

    // Report when `--max-vram` was ignored because the safe-fraction
    // heuristic already produced a smaller number. Silent no-ops are
    // confusing; users should see that their flag was considered.
    if plan.device == DeviceChoice::Gpu
        && budget.max_vram_bytes > 0
        && plan.gpu_memory_budget_bytes < budget.max_vram_bytes
    {
        logger.lock().unwrap().raw(&format!(
            "[SCHEDULER] --max-vram={} MiB was capped to {} MiB by the safe-VRAM \
             fraction over the selected GPU(s).",
            budget.max_vram_bytes / (1024 * 1024),
            plan.gpu_memory_budget_bytes / (1024 * 1024),
        ))?;
    }

    {
        let mut lg = logger.lock().unwrap();
        lg.raw(&format!(
            "[SCHEDULER] cores={} phasing_threads={} vcf2fasta_threads={} max_mem={} MiB",
            budget.total_cores,
            budget.phasing_threads,
            budget.vcf2fasta_threads,
            budget.max_memory_bytes / (1024 * 1024),
        ))?;
        lg.raw(&format!(
            "[SCHEDULER] device={} workers={} hap_block={} base_block={} streams={} in_flight={}",
            plan.device.as_str(),
            plan.cpu_workers,
            plan.haplotype_block_size,
            plan.variant_block_size,
            plan.gpu_streams,
            plan.in_flight_buffers,
        ))?;
        lg.raw(&format!("[SCHEDULER] reason: {}", plan.reason))?;
    }

    if plan.device == DeviceChoice::Gpu && args.vcf2fasta_workers > 1 {
        logger.lock().unwrap().raw(&format!(
            "WARNING: --vcf2fasta-workers={} requested but GPU mode uses 1 consumer; ignoring",
            args.vcf2fasta_workers
        ))?;
    }

    let phasing_threads = resolve_phasing_threads(&args);
    let mut budget = budget;
    if args.phasing_threads.is_some() {
        budget.phasing_threads = phasing_threads;
    }

    let outcome = pipeline::run_pipeline(
        discovery.work_units,
        &plan,
        &budget,
        &args,
        logger.clone(),
    )
    .with_context(|| "pipeline failed")?;

    if outcome.total_warnings > 0 {
        let parts: Vec<String> = outcome
            .warnings_by_reason
            .iter()
            .map(|(r, n)| format!("{}={}", r, n))
            .collect();
        logger.lock().unwrap().tag(
            "WARNINGS",
            &format!(
                "summary: total={} {{{}}} (details in {})",
                outcome.total_warnings,
                parts.join(", "),
                warnings_path.display(),
            ),
        )?;
    }

    {
        let mut lg = logger.lock().unwrap();
        lg.summary(&[
            format!("work units completed: {}", outcome.contigs_completed),
            format!("work units failed:    {}", outcome.contigs_failed),
            format!("total seen:           {}", outcome.total_seen),
            format!("total applied:        {}", outcome.total_applied),
            format!("output files:         {}", outcome.total_files),
            format!("total warnings:       {}", outcome.total_warnings),
            format!(
                "peak RSS:             {} MiB",
                budget.peak_rss_bytes() / (1024 * 1024)
            ),
        ])?;
        for (c, r) in &outcome.failures {
            lg.failure("pipeline", c, r)?;
        }
        lg.flush_warnings()?;
    }

    if let Ok(d) = SystemTime::now().duration_since(start_time) {
        logger
            .lock()
            .unwrap()
            .raw(&format!("Duration: {}", fmt_dur(d)))?;
    }
    logger.lock().unwrap().raw("=== End of log ===")?;
    Ok(())
}

fn fmt_dur(d: Duration) -> String {
    let ms = d.as_millis();
    if ms < 1000 {
        format!("{}ms", ms)
    } else if ms < 60_000 {
        format!("{:.2}s", ms as f64 / 1000.0)
    } else {
        format!("{}m{:.0}s", ms / 60_000, (ms % 60_000) as f64 / 1000.0)
    }
}

// ---------------------------------------------------------------------------
// Public re-exports used by pipeline
// ---------------------------------------------------------------------------

pub fn run_cpu_tiles_public(
    tiles: &[crate::chunk::Tile],
    sample_names: &[String],
    sample_hap_offset: &[usize],
    max_ploidies: &[usize],
    pad: u64,
    args: &Args,
    threads: usize,
) -> Vec<Result<crate::chunk::TileResult>> {
    use rayon::prelude::*;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .expect("could not create Rayon thread pool");
    pool.install(|| {
        tiles
            .par_iter()
            .map(|t| {
                crate::tile_executor::execute_tile(
                    t,
                    sample_names,
                    sample_hap_offset,
                    max_ploidies,
                    pad,
                    args,
                )
                .with_context(|| {
                    format!(
                        "tile {} haps {}-{} bases {}-{}",
                        t.contig, t.hap_start, t.hap_end, t.base_start, t.base_end,
                    )
                })
            })
            .collect()
    })
}

/// How often the collector emits a memory/telemetry line.
const TELEMETRY_INTERVAL: Duration = Duration::from_secs(5);

/// Run a contig's tiles through the CUDA pipeline, streaming results back
/// **in tile order** as they become available.
///
/// ## Producer model
///
/// A single producer thread iterates through the tile list in **windows**.
/// Within each window, tiles are built in parallel by a Rayon pool. The
/// producer waits for the whole window to complete, then submits its
/// batches to the pipeline **in tile order**, and advances to the next
/// window.
///
/// This replaces the earlier shared permit pool, which could deadlock when
/// two producers raced ahead of the collector's `next_expected` window.
/// The windowed design has no cross-thread permit handoff, so no such
/// cycle exists: the producer submits tiles in monotonic index order, and
/// the only out-of-order results the collector sees are bounded by the
/// number of CUDA workers.
///
/// ## Bounds
///
/// * **In-flight tiles**: at most `window_size` under construction plus
///   `2 * in_flight + n_streams` in the pipeline.
/// * **Ordered window** at the collector: at most `n_streams` out-of-order
///   results.
/// * **Threads**: `n_producers = min(--threads, memory_cap)`,
///   `n_streams = min(gpu_streams, n_producers, tiles)`, plus the driver
///   and collector threads.
///
/// Peak host RSS is independent of sample count, haplotype count, variant
/// count, and chromosome count.
#[cfg(feature = "cuda")]
pub fn run_gpu_tiles<F>(
    tiles: &[crate::chunk::Tile],
    sample_names: &[String],
    sample_hap_offset: &[usize],
    max_ploidies: &[usize],
    pad: u64,
    args: &Args,
    gpu_streams: usize,
    in_flight: usize,
    gpu_device_indices: &[usize],
    budget: &ResourceBudget,
    logger: &Arc<Mutex<Logger>>,
    mut on_result: F,
) -> Result<()>
where
    F: FnMut(crate::chunk::TileResult) -> Result<()>,
{
    use crate::cuda_pipeline::{CudaPipeline, PendingBatch, PipelineStats};
    use crate::tile_executor::{build_gpu_batch, new_contig_ref_cache, HostResources};
    use rayon::prelude::*;
    use std::collections::HashMap;

    if tiles.is_empty() {
        return Ok(());
    }

    // ---- Resolve every device the plan asked for -----------------------
    //
    // The plan guarantees this list is non-empty when device == Gpu. If an
    // older caller left it empty, fall back to device 0 so we do not
    // silently do nothing.
    let indices: Vec<usize> = if gpu_device_indices.is_empty() {
        vec![0]
    } else {
        gpu_device_indices.to_vec()
    };

    let mut devices = Vec::with_capacity(indices.len());
    for &idx in &indices {
        let dev = get_cuda_device(idx)?;
        logger.lock().unwrap().raw(&format!(
            "[GPU] device {}: {} ({} MiB free)",
            idx,
            dev.name,
            dev.free_vram_bytes / (1024 * 1024)
        ))?;
        devices.push(dev);
    }

    // ---- Per-tile memory estimate --------------------------------------
    let (max_hap_count, max_base_len) = tiles.iter().fold((1usize, 1u64), |acc, t| {
        (
            acc.0.max(t.hap_count().max(1)),
            acc.1.max(t.base_len().max(1)),
        )
    });
    let est_variants_per_tile = ((max_base_len / 20) as usize).max(1);
    let per_batch_input_est = max_hap_count
        .saturating_mul(est_variants_per_tile)
        .saturating_mul(4) as u64;
    let per_tile_output_est = max_hap_count.saturating_mul(max_base_len as usize) as u64;
    let per_producer_budget = per_batch_input_est
        .saturating_add(per_tile_output_est)
        .saturating_mul(3)
        .max(1);
    let memory_cap = (budget.max_memory_bytes / per_producer_budget).max(1) as usize;

    // ---- Effective thread allocation -----------------------------------
    let requested_producers = args.threads.max(1);
    let n_producers = requested_producers.min(memory_cap);

    // Streams: bounded by `gpu_streams`, by the number of producer threads
    // that can feed them, and by the number of tiles. When multiple GPUs
    // are selected the plan already sized `gpu_streams` so each device gets
    // at least one stream, but a single producer thread cannot feed two
    // streams — so we still cap by `n_producers` (with a floor of the
    // device count so every selected device gets a stream).
    let requested_streams = gpu_streams.max(devices.len());
    let effective_streams = requested_streams
        .min(n_producers.max(devices.len()))
        .min(tiles.len().max(1));

    // ---- Window size ----------------------------------------------------
    //
    // Large enough to keep the pipeline fed across the window boundary,
    // small enough that peak memory stays bounded. The window governs how
    // many tiles can be under construction at once; the pipeline's own
    // bounded channels govern how many can be in flight after submission.
    let pipeline_capacity = in_flight.saturating_mul(2) + effective_streams;
    let window_size = (pipeline_capacity + n_producers * 4).clamp(8, 64);

    let window_mem_mib = ((window_size as u64) * per_producer_budget
        + (pipeline_capacity as u64) * per_tile_output_est)
        / (1024 * 1024);

    logger.lock().unwrap().raw(&format!(
        "[GPU] host threads: producers={} streams={} devices={} window={} \
         (--threads={}, memory_cap={}, tiles={}, ~{} MiB peak)",
        n_producers,
        effective_streams,
        devices.len(),
        window_size,
        requested_producers,
        memory_cap,
        tiles.len(),
        window_mem_mib,
    ))?;
    if n_producers < requested_producers {
        logger.lock().unwrap().raw(&format!(
            "[GPU] WARNING: --threads={} reduced to {} by the memory budget \
             ({} MiB available, ~{} MiB per producer).",
            requested_producers,
            n_producers,
            budget.max_memory_bytes / (1024 * 1024),
            per_producer_budget / (1024 * 1024),
        ))?;
    }
    if effective_streams < requested_streams {
        logger.lock().unwrap().raw(&format!(
            "[GPU] stream count reduced from {} to {} (devices={}, tiles={}, producers={})",
            requested_streams,
            effective_streams,
            devices.len(),
            tiles.len(),
            n_producers,
        ))?;
    }

    let stats = PipelineStats::new();
    let pipeline = CudaPipeline::new(
        devices,
        effective_streams,
        in_flight,
        stats.clone(),
    );
    let (in_tx, out_rx, handles) = pipeline.spawn()?;

    let tiles_v = tiles.to_vec();
    let sample_names_v = sample_names.to_vec();
    let offsets_v = sample_hap_offset.to_vec();
    let ploidies_v = max_ploidies.to_vec();
    let args_v = args.clone();
    let contig_refs = new_contig_ref_cache();

    let producer = std::thread::Builder::new()
        .name("gpu-producer".to_string())
        .spawn(move || -> Result<()> {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n_producers)
                .thread_name(|i| format!("gpu-producer-{}", i))
                .build()
                .map_err(|e| anyhow::anyhow!("could not build producer pool: {e}"))?;

            let total = tiles_v.len();
            let mut start = 0usize;

            while start < total {
                let end = (start + window_size).min(total);

                // Build this window in parallel. `map_init` gives each
                // Rayon worker its own `HostResources` (reused across the
                // whole run), so the reference and VCF readers are opened
                // once per producer thread, not per tile.
                type WindowItem = (
                    usize,
                    crate::gpu::GpuBatch,
                    Arc<Vec<u8>>,
                    crate::gpu::BatchStats,
                );
                let results: Vec<Result<WindowItem>> = pool.install(|| {
                    tiles_v[start..end]
                        .par_iter()
                        .enumerate()
                        .map_init(
                            || {
                                HostResources::new(
                                    &args_v.reference,
                                    &args_v.input,
                                    contig_refs.clone(),
                                )
                            },
                            |host_result, (i, tile)| -> Result<WindowItem> {
                                if host_result.is_err() {
                                    *host_result = HostResources::new(
                                        &args_v.reference,
                                        &args_v.input,
                                        contig_refs.clone(),
                                    );
                                }
                                let host = host_result.as_mut().map_err(|e| {
                                    anyhow::anyhow!(
                                        "could not open reference/VCF for producer: {:#}",
                                        e
                                    )
                                })?;
                                let (batch, cref, bstats) = build_gpu_batch(
                                    tile,
                                    host,
                                    &sample_names_v,
                                    &offsets_v,
                                    &ploidies_v,
                                    pad,
                                    &args_v,
                                )?;
                                Ok((start + i, batch, cref, bstats))
                            },
                        )
                        .collect()
                });

                // Submit in tile order. This is the only ordering guarantee
                // the collector needs; after this, workers may complete out
                // of order but bounded by `effective_streams`.
                for r in results {
                    let (idx, batch, cref, bstats) = r?;
                    in_tx
                        .send(PendingBatch {
                            index: idx,
                            batch,
                            contig_ref: cref,
                            stats: bstats,
                        })
                        .map_err(|_| {
                            anyhow::anyhow!("GPU pipeline closed unexpectedly")
                        })?;
                }

                start = end;
            }
            Ok(())
        })
        .map_err(|e| anyhow::anyhow!("could not spawn GPU producer: {e}"))?;

    // ---- Ordered-window collection with memory telemetry ---------------
    let mut pending: HashMap<usize, crate::chunk::TileResult> = HashMap::new();
    let mut next_expected = 0usize;
    let total = tiles.len();
    let mut last_telemetry = Instant::now();
    let mut peak_pending = 0usize;
    let mut peak_rss_mib = 0u64;

    while next_expected < total {
        let (idx, res) = match out_rx.recv() {
            Ok(x) => x,
            Err(_) => break,
        };
        match res {
            Ok((hap_seqs, _timing, batch_stats)) => {
                let tile = &tiles[idx];
                if hap_seqs.len() != tile.hap_count() {
                    return Err(anyhow::anyhow!(
                        "GPU returned {} haplotypes, expected {}",
                        hap_seqs.len(),
                        tile.hap_count()
                    ));
                }
                let tr = crate::chunk::TileResult {
                    contig: tile.contig.clone(),
                    hap_start: tile.hap_start,
                    hap_end: tile.hap_end,
                    base_start: tile.base_start,
                    base_end: tile.base_end,
                    per_hap: hap_seqs,
                    seen: batch_stats.seen,
                    applied: batch_stats.applied,
                    warning_count: batch_stats.warning_count,
                    warnings: batch_stats.warnings,
                    warnings_by_reason: batch_stats.warnings_by_reason,
                };
                pending.insert(idx, tr);
                peak_pending = peak_pending.max(pending.len());

                while let Some(tr) = pending.remove(&next_expected) {
                    on_result(tr)?;
                    next_expected += 1;
                }
            }
            Err(e) => return Err(e),
        }

        if last_telemetry.elapsed() >= TELEMETRY_INTERVAL {
            last_telemetry = Instant::now();
            let (rss, cached, dirty) = read_system_memory_stats();
            peak_rss_mib = peak_rss_mib.max(rss);
            logger.lock().unwrap().raw(&format!(
                "[GPU] progress {}/{} tiles, pending={} (peak {}), \
                 process_RSS={} MiB, kernel_cache={} MiB, kernel_dirty={} MiB",
                next_expected,
                total,
                pending.len(),
                peak_pending,
                rss,
                cached,
                dirty,
            ))?;
        }
    }

    match producer.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e.context("GPU producer failed")),
        Err(_) => return Err(anyhow::anyhow!("GPU producer panicked")),
    }
    for h in handles {
        match h.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e.context("CUDA worker failed")),
            Err(_) => return Err(anyhow::anyhow!("CUDA worker panicked")),
        }
    }

    let (final_rss, final_cache, final_dirty) = read_system_memory_stats();
    peak_rss_mib = peak_rss_mib.max(final_rss);
    logger.lock().unwrap().raw(&format!(
        "[GPU] final: peak_process_RSS={} MiB, peak_ordered_window={}, \
         current_kernel_cache={} MiB, current_kernel_dirty={} MiB, batches={}",
        peak_rss_mib, peak_pending, final_cache, final_dirty, total,
    ))?;
    logger
        .lock()
        .unwrap()
        .raw(&format!("[GPU] {}", stats.summary()))?;
    Ok(())
}

#[cfg(not(feature = "cuda"))]
pub fn run_gpu_tiles<F>(
    _tiles: &[crate::chunk::Tile],
    _sample_names: &[String],
    _sample_hap_offset: &[usize],
    _max_ploidies: &[usize],
    _pad: u64,
    _args: &Args,
    _gpu_streams: usize,
    _in_flight: usize,
    _gpu_device_indices: &[usize],
    _budget: &ResourceBudget,
    _logger: &Arc<Mutex<Logger>>,
    _on_result: F,
) -> Result<()>
where
    F: FnMut(crate::chunk::TileResult) -> Result<()>,
{
    anyhow::bail!(
        "GPU execution requested, but this binary was built without the `cuda` feature. \
         Rebuild with `--features cuda` or use `--device cpu`."
    )
}