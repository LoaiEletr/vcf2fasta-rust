//! Phasing orchestration.
//!
//! Routing is **per contig** and **automatic**. See `crate::workunit` for
//! the classification table. The user only supplies `--bam-dir` and,
//! optionally, `--beagle-ref-panel` / `--beagle-genetic-map`; the tool
//! decides which phaser (if any) to run for each contig.
//!
//! The legacy whole-file entry point `ensure_phased_vcf` is retained for
//! the `--no-pipeline` path. New code uses `phase_contig` on a sliced VCF.

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use rust_htslib::bcf::{self, Read as BcfRead};
use rust_htslib::faidx;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex};

use crate::cli::Args;
use crate::logger::Logger;
use crate::report::ContigReport;
use crate::vcf::{
    clean_vcf_line, decode_and_validate_record, parse_raw_vcf_record, read_header_metadata,
    update_max_ploidies_from_line,
};
use crate::workunit::{PhaseBackend, WorkUnit};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhasingStrategy {
    NotNeeded,
    SingleSampleCanonical,
    SmallCohortBeagle,
    LargeCohortBeagle,
    PolyploidWhatsHap,
    VariablePloidy,
}

impl PhasingStrategy {
    pub fn describe(&self) -> &'static str {
        match self {
            Self::NotNeeded => "no phasing needed (input already phased)",
            Self::SingleSampleCanonical => {
                "single-sample canonical REF|ALT ordering (no statistical phase possible)"
            }
            Self::SmallCohortBeagle => {
                "small diploid cohort (2-10): Beagle with reduced phase-its"
            }
            Self::LargeCohortBeagle => {
                "large diploid cohort (>10): Beagle with default parameters"
            }
            Self::PolyploidWhatsHap => {
                "polyploid or mixed-ploidy cohort: WhatsHap polyphase with BAM files"
            }
            Self::VariablePloidy => {
                "ploidy varies within a sample: per-run WhatsHap with merge"
            }
        }
    }
}

#[derive(Debug)]
pub struct PhasedVcf {
    pub path: PathBuf,
    pub strategy: PhasingStrategy,
    pub temp_dir: Option<tempfile::TempDir>,
}

/// Per-contig phased output.
pub struct PhasedContig {
    /// `None` means the executor should read the original input directly
    /// (already-phased or canonical case).
    pub path: Option<PathBuf>,
    pub _temp_dir: Option<tempfile::TempDir>,
}

// ---------------------------------------------------------------------------
// Per-contig entry point (pipeline use)
// ---------------------------------------------------------------------------

pub fn phase_contig(
    slice_vcf: &Path,
    unit: &WorkUnit,
    reference: &Path,
    args: &Args,
    threads: usize,
    available_ram_bytes: u64,
    logger: &Arc<Mutex<Logger>>,
) -> Result<PhasedContig> {
    match unit.phase_backend {
        PhaseBackend::None | PhaseBackend::SingleSampleCanonical => {
            Ok(PhasedContig {
                path: None,
                _temp_dir: None,
            })
        }
        PhaseBackend::Beagle => {
            phase_contig_beagle(
                slice_vcf,
                unit,
                reference,
                args,
                threads,
                available_ram_bytes,
                logger,
            )
        }
        PhaseBackend::WhatsHap => {
            phase_contig_whatshap(slice_vcf, unit, reference, args, threads)
        }
        PhaseBackend::VariablePloidy => {
            phase_contig_variable(slice_vcf, unit, reference, args, threads)
        }
    }
}

fn phase_contig_beagle(
    slice_vcf: &Path,
    unit: &WorkUnit,
    reference: &Path,
    args: &Args,
    threads: usize,
    available_ram_bytes: u64,
    logger: &Arc<Mutex<Logger>>,
) -> Result<PhasedContig> {
    let temp_dir = tempfile::Builder::new()
        .prefix(&format!("vcf2fasta_phase_beagle_{}_", sanitize(&unit.contig)))
        .tempdir()
        .context("could not create Beagle tempdir")?;

    let ref_panel = args
        .beagle_ref_panel
        .as_ref()
        .and_then(|d| crate::beagle_resources::find_chromosome_file(d, &unit.contig));
    let genetic_map = args
        .beagle_genetic_map
        .as_ref()
        .and_then(|d| crate::beagle_resources::find_chromosome_file(d, &unit.contig));

    {
        let mut lg = logger.lock().unwrap();
        if let Some(d) = &args.beagle_ref_panel {
            match &ref_panel {
                Some(p) => lg.raw(&format!(
                    "[PHASING] {} -> ref panel: {}",
                    unit.contig,
                    p.display()
                ))?,
                None => lg.raw(&format!(
                    "[PHASING] {} -> no ref panel match in {} (running without ref=)",
                    unit.contig,
                    d.display()
                ))?,
            }
        }
        if let Some(d) = &args.beagle_genetic_map {
            match &genetic_map {
                Some(p) => lg.raw(&format!(
                    "[PHASING] {} -> genetic map: {}",
                    unit.contig,
                    p.display()
                ))?,
                None => lg.raw(&format!(
                    "[PHASING] {} -> no genetic map match in {} (running without map=)",
                    unit.contig,
                    d.display()
                ))?,
            }
        }
    }

    let cleaned = pre_filter_vcf(slice_vcf, reference, args, temp_dir.path())
        .with_context(|| format!("pre_filter_vcf failed for contig {}", unit.contig))?;

    let est = (unit.sample_count as u64)
        .saturating_mul(unit.variant_count)
        .saturating_mul(64);
    let budget = (available_ram_bytes as f64 * 0.75) as u64;
    if est > budget {
        bail!(
            "phasing preflight for contig '{}': estimated Beagle peak ≈ {} GiB \
             exceeds the {} GiB budget.",
            unit.contig,
            est / (1 << 30),
            budget / (1 << 30),
        );
    }

    let prefix = temp_dir.path().join("beagle_out");
    let out = run_beagle(
        &cleaned,
        &prefix,
        unit.sample_count,
        threads,
        available_ram_bytes,
        ref_panel.as_deref(),
        genetic_map.as_deref(),
    )
    .with_context(|| format!("Beagle failed for contig {}", unit.contig))?;

    let phased_gz = temp_dir.path().join("phased.vcf.gz");
    std::fs::rename(&out, &phased_gz)?;
    tabix_index(&phased_gz)?;
    validate_phased_vcf(&phased_gz, unit.sample_count)
        .with_context(|| format!("phased VCF validation failed for contig {}", unit.contig))?;

    Ok(PhasedContig {
        path: Some(phased_gz),
        _temp_dir: Some(temp_dir),
    })
}

fn phase_contig_whatshap(
    slice_vcf: &Path,
    unit: &WorkUnit,
    reference: &Path,
    args: &Args,
    threads: usize,
) -> Result<PhasedContig> {
    let bam_dir = args.bam_dir.as_ref().with_context(|| {
        format!(
            "contig '{}' requires read-based phasing (WhatsHap), \
             but --bam-dir was not provided.",
            unit.contig
        )
    })?;
    if !bam_dir.is_dir() {
        bail!("--bam-dir {} is not a directory", bam_dir.display());
    }

    let temp_dir = tempfile::Builder::new()
        .prefix(&format!("vcf2fasta_phase_whatshap_{}_", sanitize(&unit.contig)))
        .tempdir()
        .context("could not create WhatsHap tempdir")?;

    let cleaned = pre_filter_vcf(slice_vcf, reference, args, temp_dir.path())
        .with_context(|| format!("pre_filter_vcf failed for contig {}", unit.contig))?;

    let (sample_names, _) = read_header_metadata(&cleaned)?;
    if sample_names.len() != unit.sample_count {
        bail!(
            "contig {}: slice has {} samples, expected {}",
            unit.contig,
            sample_names.len(),
            unit.sample_count
        );
    }

    let mut bams: Vec<PathBuf> = Vec::with_capacity(sample_names.len());
    let mut missing: Vec<String> = Vec::new();
    for name in &sample_names {
        match crate::bam_resolver::find_sample_bam(bam_dir, name) {
            Some(p) => bams.push(p),
            None => missing.push(name.clone()),
        }
    }
    if !missing.is_empty() {
        bail!(
            "contig {}: no BAM found in {} for sample(s): {}\n\
             Looked for files named `<SAMPLE>.bam`, `<SAMPLE>.<suffix>.bam`, \
             or `<SAMPLE>.<suffix>.cram`.",
            unit.contig,
            bam_dir.display(),
            missing.join(", ")
        );
    }

    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (idx, &p) in unit.sample_max_ploidies.iter().enumerate() {
        groups.entry(p).or_default().push(idx);
    }

    let mut phased_vcfs: Vec<PathBuf> = Vec::new();
    for (&ploidy, sample_indices) in &groups {
        if ploidy < 2 {
            continue;
        }
        let subset = temp_dir.path().join(format!("subset_p{}.vcf", ploidy));
        write_subset_vcf(&cleaned, &sample_names, sample_indices, &subset)?;

        let subset_gz = temp_dir.path().join(format!("subset_p{}.vcf.gz", ploidy));
        bgzip_file(&subset, &subset_gz)?;
        tabix_index(&subset_gz)?;

        let group_bams: Vec<PathBuf> =
            sample_indices.iter().map(|&i| bams[i].clone()).collect();

        let out_plain = temp_dir.path().join(format!("phased_p{}.vcf", ploidy));
        run_whatshap_polyphase(
            &subset_gz,
            &group_bams,
            reference,
            ploidy,
            &out_plain,
            threads,
        )
        .with_context(|| {
            format!(
                "WhatsHap failed for contig {} at ploidy {}",
                unit.contig, ploidy
            )
        })?;

        let out_gz = temp_dir.path().join(format!("phased_p{}.vcf.gz", ploidy));
        bgzip_file(&out_plain, &out_gz)?;
        tabix_index(&out_gz)?;
        phased_vcfs.push(out_gz);
    }

    if phased_vcfs.is_empty() {
        bail!(
            "contig {}: no ploidy group required phasing, but backend was WhatsHap",
            unit.contig
        );
    }

    let phased_gz = temp_dir.path().join("phased.vcf.gz");
    if phased_vcfs.len() == 1 {
        std::fs::copy(&phased_vcfs[0], &phased_gz)?;
    } else {
        merge_phased_vcfs(&phased_vcfs, &phased_gz)?;
    }
    tabix_index(&phased_gz)?;
    validate_phased_vcf(&phased_gz, unit.sample_count)
        .with_context(|| format!("phased VCF validation failed for contig {}", unit.contig))?;

    Ok(PhasedContig {
        path: Some(phased_gz),
        _temp_dir: Some(temp_dir),
    })
}

fn phase_contig_variable(
    slice_vcf: &Path,
    unit: &WorkUnit,
    reference: &Path,
    args: &Args,
    threads: usize,
) -> Result<PhasedContig> {
    let temp_dir = tempfile::Builder::new()
        .prefix(&format!("vcf2fasta_phase_var_{}_", sanitize(&unit.contig)))
        .tempdir()
        .context("could not create variable-ploidy tempdir")?;

    let slice_gz = temp_dir.path().join("slice.vcf.gz");
    bgzip_file(slice_vcf, &slice_gz)?;
    tabix_index(&slice_gz)?;

    let (sample_names, _) = read_header_metadata(&slice_gz)?;
    let merged_path = crate::variable_ploidy::run_variable_ploidy_pipeline(
        &slice_gz,
        reference,
        args,
        &sample_names,
        threads,
    )
    .with_context(|| {
        format!(
            "variable-ploidy pipeline failed for contig {}",
            unit.contig
        )
    })?;

    let phased_gz = temp_dir.path().join("phased.vcf.gz");
    std::fs::copy(&merged_path, &phased_gz)?;
    tabix_index(&phased_gz)?;
    validate_phased_vcf(&phased_gz, unit.sample_count)
        .with_context(|| format!("phased VCF validation failed for contig {}", unit.contig))?;

    Ok(PhasedContig {
        path: Some(phased_gz),
        _temp_dir: Some(temp_dir),
    })
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

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

fn validate_phased_vcf(path: &Path, expected_samples: usize) -> Result<()> {
    let meta = std::fs::metadata(path)
        .with_context(|| format!("phased VCF {} does not exist", path.display()))?;
    if meta.len() == 0 {
        anyhow::bail!("phased VCF {} is empty", path.display());
    }
    let reader = open_maybe_gzip(path)?;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 10 {
                anyhow::bail!(
                    "phased VCF {} #CHROM line has {} columns; expected ≥10",
                    path.display(),
                    cols.len()
                );
            }
            let got = cols.len() - 9;
            if got != expected_samples {
                anyhow::bail!(
                    "phased VCF {} has {} samples; expected {}",
                    path.display(),
                    got,
                    expected_samples
                );
            }
            return Ok(());
        }
    }
    anyhow::bail!("phased VCF {} has no #CHROM line", path.display())
}

fn pick_beagle_heap(available_ram_bytes: u64) -> String {
    const MIN: u64 = 1 << 30;
    const MAX: u64 = 64u64 << 30;
    let target = ((available_ram_bytes as f64) * 0.60) as u64;
    let bytes = target.clamp(MIN, MAX);
    format!("{}g", bytes / (1 << 30))
}

/// Parse a `BEAGLE_HEAP`-style spec into bytes.
///
/// Accepts `1024`, `512m`, `32g`, `2t` (case-insensitive). A bare number
/// is interpreted as bytes, matching Java's `-Xmx` syntax.
fn parse_heap_spec(spec: &str) -> Option<u64> {
    let s = spec.trim();
    if s.is_empty() {
        return None;
    }
    let last = s.chars().last()?;
    let (num_str, mult) = match last.to_ascii_lowercase() {
        'k' => (&s[..s.len() - 1], 1024u64),
        'm' => (&s[..s.len() - 1], 1024 * 1024),
        'g' => (&s[..s.len() - 1], 1024u64 * 1024 * 1024),
        't' => (&s[..s.len() - 1], 1024u64 * 1024 * 1024 * 1024),
        _ => (s, 1),
    };
    let n: u64 = num_str.trim().parse().ok()?;
    Some(n.saturating_mul(mult))
}

/// Choose the heap size for Beagle, honoring `BEAGLE_HEAP` only if it is
/// within the safe fraction of machine RAM.
///
/// Setting `BEAGLE_HEAP=64g` on a 12 GB machine would otherwise cause the
/// JVM to request more RAM than the kernel can provide, and the process
/// would be killed by the OOM killer. We clamp to the same auto-derived
/// value we would have picked, and log the reduction.
fn clamped_beagle_heap(available_ram_bytes: u64) -> String {
    let auto = pick_beagle_heap(available_ram_bytes);
    let Ok(user) = std::env::var("BEAGLE_HEAP") else {
        return auto;
    };
    let Some(user_bytes) = parse_heap_spec(&user) else {
        eprintln!(
            "WARNING: BEAGLE_HEAP={:?} is not a valid memory spec; using auto {}",
            user, auto
        );
        return auto;
    };
    let Some(auto_bytes) = parse_heap_spec(&auto) else {
        return auto;
    };
    if user_bytes > auto_bytes {
        eprintln!(
            "WARNING: BEAGLE_HEAP={} exceeds the safe heap for {} MiB available RAM; \
             clamping to {}",
            user,
            available_ram_bytes / (1024 * 1024),
            auto
        );
        auto
    } else {
        user
    }
}

// ---------------------------------------------------------------------------
// Legacy whole-file entry point
// ---------------------------------------------------------------------------

pub fn ensure_phased_vcf(
    input: &Path,
    reference: &Path,
    args: &Args,
    threads: usize,
) -> Result<PhasedVcf> {
    let n_samples = count_samples(input)?;
    if n_samples == 0 {
        bail!("input VCF contains no samples");
    }

    if !has_unphased_genotypes(input)? {
        return Ok(PhasedVcf {
            path: input.to_path_buf(),
            strategy: PhasingStrategy::NotNeeded,
            temp_dir: None,
        });
    }

    let (sample_names, _) = read_header_metadata(input)?;
    let mut max_ploidies = vec![0usize; sample_names.len()];
    {
        let reader = open_maybe_gzip(input)?;
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            update_max_ploidies_from_line(&line, &sample_names, &mut max_ploidies);
        }
    }
    for p in &mut max_ploidies {
        if *p == 0 {
            *p = 2;
        }
    }

    if n_samples == 1 {
        if let Some(bam_dir) = &args.bam_dir {
            if let Some(bam) =
                crate::bam_resolver::find_sample_bam(bam_dir, &sample_names[0])
            {
                return run_whatshap_single(
                    input,
                    &bam,
                    reference,
                    max_ploidies[0],
                    args,
                    threads,
                );
            }
        }
        return Ok(PhasedVcf {
            path: input.to_path_buf(),
            strategy: PhasingStrategy::SingleSampleCanonical,
            temp_dir: None,
        });
    }

    let all_samples_uniform_diploid = max_ploidies.iter().all(|&p| p == 2)
        && !crate::variable_ploidy::has_variable_ploidy(input, &sample_names)?;
    if all_samples_uniform_diploid {
        return run_beagle_pipeline(input, reference, args, n_samples, threads);
    }

    if crate::variable_ploidy::has_variable_ploidy(input, &sample_names)? {
        let merged_path = crate::variable_ploidy::run_variable_ploidy_pipeline(
            input,
            reference,
            args,
            &sample_names,
            threads,
        )?;
        return Ok(PhasedVcf {
            path: merged_path,
            strategy: PhasingStrategy::VariablePloidy,
            temp_dir: None,
        });
    }

    run_whatshap_pipeline(
        input,
        reference,
        args,
        &sample_names,
        &max_ploidies,
        threads,
    )
}

// ---------------------------------------------------------------------------
// Beagle pipeline (legacy)
// ---------------------------------------------------------------------------

fn run_beagle_pipeline(
    input: &Path,
    reference: &Path,
    args: &Args,
    n_samples: usize,
    threads: usize,
) -> Result<PhasedVcf> {
    let temp_dir = tempfile::Builder::new()
        .prefix("vcf2fasta_phasing_")
        .tempdir()
        .context("could not create temporary directory for phasing")?;

    let cleaned = pre_filter_vcf(input, reference, args, temp_dir.path())?;
    let output_prefix = temp_dir.path().join("phased");
    let available_ram = crate::scheduler::detect_available_ram_bytes();
    let output_vcf = run_beagle(
        &cleaned,
        &output_prefix,
        n_samples,
        threads,
        available_ram,
        None,
        None,
    )?;
    tabix_index(&output_vcf)?;

    Ok(PhasedVcf {
        path: output_vcf,
        strategy: if n_samples <= 10 {
            PhasingStrategy::SmallCohortBeagle
        } else {
            PhasingStrategy::LargeCohortBeagle
        },
        temp_dir: Some(temp_dir),
    })
}

// ---------------------------------------------------------------------------
// WhatsHap pipeline (legacy)
// ---------------------------------------------------------------------------

fn run_whatshap_single(
    input: &Path,
    bam: &Path,
    reference: &Path,
    ploidy: usize,
    args: &Args,
    threads: usize,
) -> Result<PhasedVcf> {
    let temp_dir = tempfile::Builder::new()
        .prefix("vcf2fasta_phasing_")
        .tempdir()
        .context("could not create temporary directory for phasing")?;

    let cleaned = pre_filter_vcf(input, reference, args, temp_dir.path())?;
    let out_plain = temp_dir.path().join("phased.vcf");
    run_whatshap_polyphase(
        &cleaned,
        &[bam.to_path_buf()],
        reference,
        ploidy,
        &out_plain,
        threads,
    )?;

    let out_gz = temp_dir.path().join("phased.vcf.gz");
    bgzip_file(&out_plain, &out_gz)?;
    tabix_index(&out_gz)?;

    Ok(PhasedVcf {
        path: out_gz,
        strategy: PhasingStrategy::PolyploidWhatsHap,
        temp_dir: Some(temp_dir),
    })
}

fn run_whatshap_pipeline(
    input: &Path,
    reference: &Path,
    args: &Args,
    sample_names: &[String],
    max_ploidies: &[usize],
    threads: usize,
) -> Result<PhasedVcf> {
    let bam_dir = match &args.bam_dir {
        Some(d) => d,
        None => bail!(
            "This VCF contains samples with ploidy > 2.\n\
             Statistical phasing of polyploids requires aligned reads. \
             Please supply a directory containing one BAM per sample via --bam-dir:\n\
             \tvcf2fasta -f <ref.fa> --bam-dir /path/to/bams -p <prefix> <input.vcf.gz>\n\
             Each BAM must be coordinate-sorted, indexed (.bai or .csi), and \
             named `<SAMPLE>.bam` (or `<SAMPLE>.<suffix>.bam`).\n\
             Samples: {:?}",
            sample_names
        ),
    };
    if !bam_dir.is_dir() {
        bail!(
            "--bam-dir {} does not exist or is not a directory",
            bam_dir.display()
        );
    }

    let mut bam_paths: Vec<(String, PathBuf)> = Vec::with_capacity(sample_names.len());
    let mut missing_bams: Vec<String> = Vec::new();
    for name in sample_names {
        match crate::bam_resolver::find_sample_bam(bam_dir, name) {
            Some(p) => bam_paths.push((name.clone(), p)),
            None => missing_bams.push(name.clone()),
        }
    }
    if !missing_bams.is_empty() {
        bail!(
            "No BAM file found in {} for the following samples:\n  {}\n\
             Looked for files named `<SAMPLE>.bam`, `<SAMPLE>.<suffix>.bam`, \
             or `<SAMPLE>.<suffix>.cram`.",
            bam_dir.display(),
            missing_bams.join("\n  ")
        );
    }

    let temp_dir = tempfile::Builder::new()
        .prefix("vcf2fasta_phasing_")
        .tempdir()
        .context("could not create temporary directory for phasing")?;

    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (idx, &ploidy) in max_ploidies.iter().enumerate() {
        groups.entry(ploidy).or_default().push(idx);
    }

    let mut phased_vcfs: Vec<PathBuf> = Vec::new();

    for (&ploidy, sample_indices) in &groups {
        if ploidy < 2 {
            continue;
        }
        let subset_plain = temp_dir.path().join(format!("subset_p{}.vcf", ploidy));
        write_subset_vcf(input, sample_names, sample_indices, &subset_plain)?;

        let subset_gz = temp_dir.path().join(format!("subset_p{}.vcf.gz", ploidy));
        bgzip_file(&subset_plain, &subset_gz)?;
        tabix_index(&subset_gz)?;
        let cleaned = pre_filter_vcf(&subset_gz, reference, args, temp_dir.path())?;

        let group_bams: Vec<PathBuf> = sample_indices
            .iter()
            .map(|&i| bam_paths[i].1.clone())
            .collect();

        let out_plain = temp_dir.path().join(format!("phased_p{}.vcf", ploidy));
        run_whatshap_polyphase(
            &cleaned,
            &group_bams,
            reference,
            ploidy,
            &out_plain,
            threads,
        )?;

        let out_gz = temp_dir.path().join(format!("phased_p{}.vcf.gz", ploidy));
        bgzip_file(&out_plain, &out_gz)?;
        tabix_index(&out_gz)?;
        phased_vcfs.push(out_gz);
    }

    if phased_vcfs.is_empty() {
        bail!("no ploidy group produced a phased VCF");
    }

    let merged = temp_dir.path().join("merged.vcf.gz");
    if phased_vcfs.len() == 1 {
        std::fs::copy(&phased_vcfs[0], &merged)?;
    } else {
        merge_phased_vcfs(&phased_vcfs, &merged)?;
    }
    tabix_index(&merged)?;

    Ok(PhasedVcf {
        path: merged,
        strategy: PhasingStrategy::PolyploidWhatsHap,
        temp_dir: Some(temp_dir),
    })
}

pub(crate) fn write_subset_vcf(
    input: &Path,
    all_samples: &[String],
    sample_indices: &[usize],
    output: &Path,
) -> Result<()> {
    let reader = open_maybe_gzip(input)?;
    let mut out = File::create(output)?;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("#CHROM") {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                bail!("malformed #CHROM line in input VCF");
            }
            let mut new_fields: Vec<String> =
                fields[..9].iter().map(|s| s.to_string()).collect();
            for &idx in sample_indices {
                new_fields.push(all_samples[idx].clone());
            }
            out.write_all(new_fields.join("\t").as_bytes())?;
            out.write_all(b"\n")?;
            continue;
        }
        if line.starts_with('#') {
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        let mut new_fields: Vec<String> = fields[..9].iter().map(|s| s.to_string()).collect();
        for &idx in sample_indices {
            if 9 + idx < fields.len() {
                new_fields.push(fields[9 + idx].to_string());
            } else {
                new_fields.push(".".to_string());
            }
        }
        out.write_all(new_fields.join("\t").as_bytes())?;
        out.write_all(b"\n")?;
    }
    out.flush()?;
    Ok(())
}

pub(crate) fn run_whatshap_polyphase(
    vcf: &Path,
    bams: &[PathBuf],
    reference: &Path,
    ploidy: usize,
    output: &Path,
    threads: usize,
) -> Result<()> {
    let whatshap = find_in_path("whatshap").context(
        "WhatsHap is required for polyploid phasing but was not found in PATH.\n\
         Install it from https://whatshap.readthedocs.io/.",
    )?;
    for bam in bams {
        if !bam_is_indexed(bam) {
            bail!(
                "BAM/CRAM file {} has no index.\n\
                 Expected one of: {}.bai, {}.csi, {}.crai, {}.csi\n\
                 Run `samtools index {}` first.",
                bam.display(),
                bam.display(),
                bam.display(),
                bam.display(),
                bam.display(),
                bam.display(),
            );
        }
    }
    let nthreads = threads.max(1);
    let mut cmd = Command::new(&whatshap);
    cmd.arg("polyphase")
        .arg("--ploidy")
        .arg(ploidy.to_string())
        .arg("--reference")
        .arg(reference)
        .arg("--output")
        .arg(output)
        .arg("--threads")
        .arg(nthreads.to_string())
        .arg(vcf);
    for bam in bams {
        cmd.arg(bam);
    }

    let status = cmd
        .status()
        .with_context(|| format!("failed to spawn {}", whatshap.display()))?;
    if !status.success() {
        bail!("whatshap polyphase exited with non-zero status: {}", status);
    }
    let meta = std::fs::metadata(output).with_context(|| {
        format!(
            "whatshap polyphase reported success but {} does not exist",
            output.display()
        )
    })?;
    if meta.len() == 0 {
        bail!(
            "whatshap polyphase produced an empty file at {}",
            output.display()
        );
    }
    Ok(())
}

fn bam_is_indexed(bam: &Path) -> bool {
    let append_bai = PathBuf::from(format!("{}.bai", bam.display()));
    let append_csi = PathBuf::from(format!("{}.csi", bam.display()));
    let append_crai = PathBuf::from(format!("{}.crai", bam.display()));
    if append_bai.is_file() || append_csi.is_file() || append_crai.is_file() {
        return true;
    }
    if let Some(stem) = bam.file_stem() {
        let mut alt = bam.to_path_buf();
        alt.set_file_name(stem);
        let alt_bai = alt.with_extension("bai");
        let alt_crai = alt.with_extension("crai");
        if alt_bai.is_file() || alt_crai.is_file() {
            return true;
        }
    }
    false
}

pub(crate) fn merge_phased_vcfs(vcfs: &[PathBuf], output: &Path) -> Result<()> {
    let bcftools = find_in_path("bcftools")
        .context("bcftools is required to merge per-ploidy phased VCFs.")?;
    let mut cmd = Command::new(&bcftools);
    cmd.arg("merge")
        .arg("--output-type")
        .arg("z")
        .arg("--output")
        .arg(output);
    for vcf in vcfs {
        cmd.arg(vcf);
    }
    let status = cmd
        .status()
        .with_context(|| format!("failed to spawn {}", bcftools.display()))?;
    if !status.success() {
        bail!("bcftools merge exited with non-zero status: {}", status);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// VCF inspection helpers
// ---------------------------------------------------------------------------

fn count_samples(input: &Path) -> Result<usize> {
    let reader = bcf::Reader::from_path(input)
        .with_context(|| format!("could not open {}", input.display()))?;
    Ok(reader.header().samples().len())
}

pub(crate) fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file =
        File::open(path).with_context(|| format!("could not open {}", path.display()))?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn has_unphased_genotypes(input: &Path) -> Result<bool> {
    let reader = open_maybe_gzip(input)?;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue;
        }
        let format = fields[8];
        let gt_idx = match format.split(':').position(|f| f == "GT") {
            Some(i) => i,
            None => continue,
        };
        for sample_field in fields.iter().skip(9) {
            let parts: Vec<&str> = sample_field.split(':').collect();
            if let Some(gt) = parts.get(gt_idx) {
                if gt.contains('/') {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

// ---------------------------------------------------------------------------
// Pre-filter
// ---------------------------------------------------------------------------

pub(crate) fn pre_filter_vcf(
    input: &Path,
    reference: &Path,
    args: &Args,
    temp_dir: &Path,
) -> Result<PathBuf> {
    let reference_reader = faidx::Reader::from_path(reference)
        .with_context(|| format!("could not open reference {}", reference.display()))?;
    let contig_names = reference_reader.seq_names()?;
    let mut contig_lengths: HashMap<String, usize> = HashMap::with_capacity(contig_names.len());
    for name in &contig_names {
        let len = reference_reader.fetch_seq_len(name);
        contig_lengths.insert(name.clone(), len as usize);
    }

    let (sample_names, _) = read_header_metadata(input)?;
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let cleaned_plain = temp_dir.join(format!("cleaned_{}.vcf", stamp));
    let cleaned_gz = temp_dir.join(format!("cleaned_{}.vcf.gz", stamp));

    let mut max_ploidies = vec![0usize; sample_names.len()];
    {
        let reader = open_maybe_gzip(input)?;
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            update_max_ploidies_from_line(&line, &sample_names, &mut max_ploidies);
        }
    }
    for p in &mut max_ploidies {
        if *p == 0 {
            *p = 2;
        }
    }

    let mut kept = 0usize;
    let mut dropped_contig = 0usize;
    let mut dropped_invalid = 0usize;
    let mut seen_contigs: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();

    {
        let reader = open_maybe_gzip(input)?;
        let mut out = File::create(&cleaned_plain)?;
        let mut report = ContigReport::new(String::from("__pre_filter__"));
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                if line.starts_with("##FORMAT=") && !line.starts_with("##FORMAT=<ID=GT") {
                    continue;
                }
                out.write_all(line.as_bytes())?;
                out.write_all(b"\n")?;
                continue;
            }
            let raw = match parse_raw_vcf_record(&line) {
                Ok(r) => r,
                Err(_) => {
                    dropped_invalid += 1;
                    continue;
                }
            };
            seen_contigs.insert(raw.chrom.to_string());
            let contig_len = match contig_lengths.get(raw.chrom) {
                Some(&l) if l > 0 => l,
                _ => {
                    dropped_contig += 1;
                    continue;
                }
            };
            match decode_and_validate_record(
                &raw,
                raw.chrom,
                contig_len,
                &sample_names,
                &max_ploidies,
                &reference_reader,
                args,
                &mut report,
            ) {
                Ok(_) => match clean_vcf_line(&line, &sample_names, &max_ploidies) {
                    Some(new_line) => {
                        out.write_all(new_line.as_bytes())?;
                        out.write_all(b"\n")?;
                        kept += 1;
                    }
                    None => {
                        dropped_invalid += 1;
                    }
                },
                Err(_) => {
                    dropped_invalid += 1;
                    continue;
                }
            }
        }
        out.flush()?;
    }

    if kept == 0 {
        let fasta_contigs: Vec<&String> = contig_lengths.keys().collect();
        let vcf_contigs: Vec<String> = seen_contigs.iter().cloned().collect();
        bail!(
            "pre-filter removed every record from the input VCF.\n\
             VCF contigs seen:  {:?}\nFASTA contigs:     {:?}\n\
             Records dropped (contig not in FASTA): {}\n\
             Records dropped (validation failed):  {}",
            vcf_contigs,
            fasta_contigs,
            dropped_contig,
            dropped_invalid
        );
    }

    bgzip_file(&cleaned_plain, &cleaned_gz)?;
    tabix_index(&cleaned_gz)?;
    Ok(cleaned_gz)
}

// ---------------------------------------------------------------------------
// Subprocess helpers
// ---------------------------------------------------------------------------

fn run_beagle(
    input: &Path,
    output_prefix: &Path,
    n_samples: usize,
    threads: usize,
    available_ram_bytes: u64,
    ref_panel: Option<&Path>,
    genetic_map: Option<&Path>,
) -> Result<PathBuf> {
    let java = find_in_path("java")
        .context("Java is required to run Beagle but was not found in PATH.")?;
    let beagle_jar = locate_beagle_jar().context(
        "Beagle jar not found. Set BEAGLE_JAR=/path/to/beagle.jar, or put \
         `beagle.jar` on your PATH.",
    )?;

    let heap = clamped_beagle_heap(available_ram_bytes);
    let nthreads = threads.max(1);
    let burnin = if n_samples <= 10 { 3 } else { 6 };
    let iterations = if n_samples <= 10 { 12 } else { 20 };

    let mut cmd = Command::new(&java);
    cmd.arg(format!("-Xmx{}", heap))
        .arg("-XX:+UseParallelGC")
        .arg("-jar")
        .arg(&beagle_jar)
        .arg(format!("gt={}", input.display()))
        .arg(format!("out={}", output_prefix.display()))
        .arg(format!("nthreads={}", nthreads))
        .arg(format!("iterations={}", iterations))
        .arg(format!("burnin={}", burnin));

    if let Some(p) = ref_panel {
        cmd.arg(format!("ref={}", p.display()));
    }
    if let Some(p) = genetic_map {
        cmd.arg(format!("map={}", p.display()));
    }

    let status = cmd
        .status()
        .with_context(|| format!("failed to spawn Beagle ({})", java.display()))?;
    if !status.success() {
        bail!("Beagle exited with non-zero status: {}", status);
    }
    let output_vcf = PathBuf::from(format!("{}.vcf.gz", output_prefix.display()));
    let meta = std::fs::metadata(&output_vcf).with_context(|| {
        format!(
            "Beagle reported success but did not produce {}",
            output_vcf.display()
        )
    })?;
    if meta.len() == 0 {
        bail!("Beagle produced an empty file at {}", output_vcf.display());
    }
    Ok(output_vcf)
}

pub(crate) fn bgzip_file(input: &Path, output: &Path) -> Result<()> {
    let bgzip = find_in_path("bgzip").context("bgzip is required")?;
    let out = Command::new(&bgzip)
        .arg("-c")
        .arg(input)
        .output()
        .with_context(|| "failed to run bgzip")?;
    if !out.status.success() {
        bail!("bgzip exited with non-zero status: {}", out.status);
    }
    std::fs::write(output, &out.stdout)
        .with_context(|| format!("could not write {}", output.display()))?;
    Ok(())
}

pub(crate) fn tabix_index(vcf: &Path) -> Result<()> {
    let tabix = find_in_path("tabix")
        .context("tabix is required (and should already be present)")?;
    let status = Command::new(&tabix)
        .arg("-p")
        .arg("vcf")
        .arg(vcf)
        .status()
        .with_context(|| format!("failed to run tabix on {}", vcf.display()))?;
    if !status.success() {
        bail!("tabix exited with non-zero status: {}", status);
    }
    Ok(())
}

pub(crate) fn find_in_path(name: &str) -> Option<PathBuf> {
    let path_var = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path_var) {
        let candidate = dir.join(name);
        if candidate.is_file() {
            return Some(candidate);
        }
    }
    None
}

fn locate_beagle_jar() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("BEAGLE_JAR") {
        let path = PathBuf::from(p);
        if path.is_file() {
            return Some(path);
        }
    }
    for name in [
        "beagle.jar",
        "beagle.22Jul22.46e.jar",
        "beagle.08Jun19.f15.jar",
        "beagle.27Jan18.7e1.jar",
    ] {
        if let Some(p) = find_in_path(name) {
            return Some(p);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use tempfile::tempdir;

    fn write_file(path: &Path, content: &[u8]) {
        let mut f = File::create(path).unwrap();
        f.write_all(content).unwrap();
        f.flush().unwrap();
    }

    fn write_gzip(path: &Path, content: &[u8]) {
        let f = File::create(path).unwrap();
        let mut enc = GzEncoder::new(f, Compression::default());
        enc.write_all(content).unwrap();
        enc.finish().unwrap();
    }

    fn test_args() -> crate::cli::Args {
        crate::cli::Args {
            input: PathBuf::from("dummy.vcf"),
            reference: PathBuf::from("dummy.fa"),
            prefix: String::new(),
            no_call_string: Some("N".to_string()),
            threads: 1,
            line_width: 80,
            no_validate_ref: false,
            quiet: true,
            gpu: false,
            device: crate::scheduler::Device::Auto,
            bam_dir: None,
            beagle_ref_panel: None,
            beagle_genetic_map: None,
            chunk_size: None,
            chunk_pad: None,
            phasing_threads: None,
            vcf2fasta_workers: 1,
            ready_queue_depth: 4,
            no_pipeline: false,
            merged_output: false,
            max_memory: None,
            max_vram: None,
            gpu_devices: None,
        }
    }

    const TWO_DIPLOID_PHASED: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=20>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|0
";

    const TWO_DIPLOID_UNPHASED: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=20>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/0
";

    const ONE_DIPLOID_UNPHASED: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=20>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1
";

    #[test]
    fn strategy_describe_all() {
        assert!(PhasingStrategy::NotNeeded.describe().contains("no phasing"));
        assert!(PhasingStrategy::SingleSampleCanonical
            .describe()
            .contains("single-sample"));
        assert!(PhasingStrategy::SmallCohortBeagle
            .describe()
            .contains("small diploid"));
        assert!(PhasingStrategy::LargeCohortBeagle
            .describe()
            .contains("large diploid"));
        assert!(PhasingStrategy::PolyploidWhatsHap
            .describe()
            .contains("WhatsHap"));
        assert!(PhasingStrategy::VariablePloidy.describe().contains("varies"));
    }

    #[test]
    fn open_maybe_gzip_reads_plain() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("plain.vcf");
        write_file(&p, b"hello\nworld\n");
        let mut r = open_maybe_gzip(&p).unwrap();
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        assert_eq!(s, "hello\nworld\n");
    }

    #[test]
    fn open_maybe_gzip_reads_gzipped() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("x.vcf.gz");
        write_gzip(&p, b"hello\nworld\n");
        let mut r = open_maybe_gzip(&p).unwrap();
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        assert_eq!(s, "hello\nworld\n");
    }

    #[test]
    fn unphased_detection_basic() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("u.vcf");
        write_file(&p, TWO_DIPLOID_UNPHASED.as_bytes());
        assert!(has_unphased_genotypes(&p).unwrap());

        let p2 = dir.path().join("p.vcf");
        write_file(&p2, TWO_DIPLOID_PHASED.as_bytes());
        assert!(!has_unphased_genotypes(&p2).unwrap());
    }

    #[test]
    fn count_samples_works() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("two.vcf");
        write_file(&p, TWO_DIPLOID_PHASED.as_bytes());
        assert_eq!(count_samples(&p).unwrap(), 2);
    }

    #[test]
    fn passthrough_when_already_phased() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("p.vcf");
        write_file(&p, TWO_DIPLOID_PHASED.as_bytes());
        let ref_path = dir.path().join("ref.fa");
        let args = test_args();
        let result = ensure_phased_vcf(&p, &ref_path, &args, 1).unwrap();
        assert_eq!(result.strategy, PhasingStrategy::NotNeeded);
        assert_eq!(result.path, p);
    }

    #[test]
    fn single_sample_without_bam_canonicalises() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("one.vcf");
        write_file(&p, ONE_DIPLOID_UNPHASED.as_bytes());
        let ref_path = dir.path().join("ref.fa");
        let args = test_args();
        let result = ensure_phased_vcf(&p, &ref_path, &args, 1).unwrap();
        assert_eq!(result.strategy, PhasingStrategy::SingleSampleCanonical);
        assert!(result.temp_dir.is_none());
    }

    #[test]
    fn bam_is_indexed_accepts_standard_and_legacy_indexes() {
        let dir = tempdir().unwrap();

        let std_bam = dir.path().join("S1.bam");
        File::create(&std_bam).unwrap();
        File::create(dir.path().join("S1.bam.bai")).unwrap();
        assert!(bam_is_indexed(&std_bam));

        let cram = dir.path().join("S2.cram");
        File::create(&cram).unwrap();
        File::create(dir.path().join("S2.cram.crai")).unwrap();
        assert!(bam_is_indexed(&cram));

        let legacy = dir.path().join("S3.bam");
        File::create(&legacy).unwrap();
        File::create(dir.path().join("S3.bai")).unwrap();
        assert!(bam_is_indexed(&legacy));

        let missing = dir.path().join("S4.bam");
        File::create(&missing).unwrap();
        assert!(!bam_is_indexed(&missing));
    }

    #[test]
    fn parse_heap_spec_variants() {
        assert_eq!(parse_heap_spec("1024"), Some(1024));
        assert_eq!(parse_heap_spec("1k"), Some(1024));
        assert_eq!(parse_heap_spec("1K"), Some(1024));
        assert_eq!(parse_heap_spec("2m"), Some(2 * 1024 * 1024));
        assert_eq!(parse_heap_spec("3G"), Some(3 * 1024 * 1024 * 1024));
        assert_eq!(parse_heap_spec("1t"), Some(1024u64 * 1024 * 1024 * 1024));
        assert_eq!(parse_heap_spec(""), None);
        assert_eq!(parse_heap_spec("abc"), None);
    }
}