//! Work-unit abstraction for the phasing → vcf2fasta pipeline.
//!
//! One `WorkUnit` = one contig. Phasing decisions, Beagle/WhatsHap
//! invocation, and output files are all per-contig.
//!
//! The backend for a contig is chosen automatically from the data:
//!
//! | Condition (per contig)                                      | Backend               |
//! |-------------------------------------------------------------|-----------------------|
//! | all GTs phased, or max ploidy ≤ 1                           | None                  |
//! | N = 1, no BAM                                               | SingleSampleCanonical |
//! | N = 1, BAM present, constant ploidy                         | WhatsHap              |
//! | N = 1, BAM present, ploidy varies within the sample         | VariablePloidy        |
//! | N ≥ 2, ploidy varies within any sample on this contig       | VariablePloidy        |
//! | N ≥ 2, all samples uniform diploid                          | Beagle                |
//! | N ≥ 2, polyploid, constant per sample                       | WhatsHap              |
//!
//! The user never selects the phaser; only `--bam-dir` supplies reads.
//! `SingleSampleCanonical` is not a phaser — it means the vcf2fasta
//! executor's existing canonical REF|ALT ordering policy is used, exactly
//! as before.
//!
//! ## Lifetime contract for `phased_vcf` and `phased_tempdir`
//!
//! When a contig needs phasing, `phase_contig` writes the phased VCF into
//! a fresh `TempDir`. That tempdir MUST outlive every consumer that reads
//! the file. Because the producer thread may exit before the consumer
//! finishes, ownership of the tempdir travels with the `WorkUnit`:
//! `phased_tempdir` holds an `Arc<TempDir>` that is dropped only when the
//! last owner of the unit is dropped. Do not move the tempdir into a
//! producer-local vector — that is the race this field exists to prevent.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::cli::Args;
use crate::genotype::normalize_genotype;
use crate::phasing::open_maybe_gzip;
use crate::vcf::{extract_gt, parse_gt, parse_raw_vcf_record, read_header_metadata};

// ---------------------------------------------------------------------------
// Identifiers and states
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct WorkUnitId(pub u64);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkState {
    Discovered,
    Analyzed,
    NeedsPhasing,
    Phasing,
    PhaseValidation,
    Ready,
    Running,
    Complete,
    PhasingFailed,
    Vcf2FastaFailed,
    ValidationFailed,
    ResourceLimit,
}

impl WorkState {
    pub fn is_terminal_failure(&self) -> bool {
        matches!(
            self,
            Self::PhasingFailed
                | Self::Vcf2FastaFailed
                | Self::ValidationFailed
                | Self::ResourceLimit
        )
    }
    pub fn is_terminal_success(&self) -> bool {
        matches!(self, Self::Complete)
    }
    pub fn tag(&self) -> &'static str {
        match self {
            Self::Discovered => "DISCOVERED",
            Self::Analyzed => "ANALYZED",
            Self::NeedsPhasing => "NEEDS_PHASING",
            Self::Phasing => "PHASING",
            Self::PhaseValidation => "PHASE_VALIDATION",
            Self::Ready => "READY",
            Self::Running => "RUNNING",
            Self::Complete => "COMPLETE",
            Self::PhasingFailed => "PHASING_FAILED",
            Self::Vcf2FastaFailed => "VCF2FASTA_FAILED",
            Self::ValidationFailed => "VALIDATION_FAILED",
            Self::ResourceLimit => "RESOURCE_LIMIT",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseBackend {
    /// Nothing to do: every GT is phased, or every sample is haploid.
    None,
    /// Single sample, no reads: executor applies its existing REF|ALT
    /// canonical ordering. No subprocess, no biological inference.
    SingleSampleCanonical,
    /// Uniform diploid cohort (N ≥ 2): population phasing via Beagle.
    Beagle,
    /// Read-based phasing via WhatsHap. Used for:
    ///   * N = 1 with a BAM and constant ploidy,
    ///   * N ≥ 2 polyploid with constant-per-sample ploidy.
    WhatsHap,
    /// Ploidy changes within a sample on this contig: per-run WhatsHap
    /// with merge, via `variable_ploidy::run_variable_ploidy_pipeline`.
    VariablePloidy,
}

impl PhaseBackend {
    pub fn tag(&self) -> &'static str {
        match self {
            Self::None => "none",
            Self::SingleSampleCanonical => "canonical",
            Self::Beagle => "beagle",
            Self::WhatsHap => "whatshap",
            Self::VariablePloidy => "variable-ploidy",
        }
    }
    /// Whether this backend requires a `--bam-dir`.
    pub fn requires_bam(&self) -> bool {
        matches!(self, Self::WhatsHap | Self::VariablePloidy)
    }
    /// Whether a subprocess must actually be run.
    pub fn is_subprocess(&self) -> bool {
        matches!(self, Self::Beagle | Self::WhatsHap | Self::VariablePloidy)
    }
}

// ---------------------------------------------------------------------------
// WorkUnit
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct WorkUnit {
    pub id: WorkUnitId,
    pub contig: String,
    pub contig_index: usize,
    pub reference_length: u64,

    pub sample_count: usize,
    pub sample_max_ploidies: Vec<usize>,
    pub haplotype_count: usize,
    pub max_ploidy: usize,
    pub variant_count: u64,
    pub phased_genotypes: u64,
    pub unphased_genotypes: u64,
    pub malformed_records: u64,

    pub needs_phasing: bool,
    pub phase_backend: PhaseBackend,

    /// Path to the phased VCF. `None` means "use the original input, scoped
    /// to this contig".
    pub phased_vcf: Option<PathBuf>,

    /// Ownership guard for the directory containing `phased_vcf`. Must
    /// remain alive for the entire time any consumer might read the file.
    /// See the module docs for the ownership contract.
    pub phased_tempdir: Option<Arc<tempfile::TempDir>>,

    pub state: WorkState,
}

impl WorkUnit {
    /// The vcf2fasta stage reads from this path.
    pub fn effective_vcf<'a>(&'a self, fallback: &'a Path) -> &'a Path {
        self.phased_vcf.as_deref().unwrap_or(fallback)
    }
}

// ---------------------------------------------------------------------------
// Discovery
// ---------------------------------------------------------------------------

struct ContigAcc {
    max_ploidy: Vec<usize>,
    first_ploidy: Vec<Option<usize>>,
    variant_count: u64,
    phased_genotypes: u64,
    unphased_genotypes: u64,
    variable_ploidy: bool,
}

impl ContigAcc {
    fn new(n_samples: usize) -> Self {
        Self {
            max_ploidy: vec![0; n_samples],
            first_ploidy: vec![None; n_samples],
            variant_count: 0,
            phased_genotypes: 0,
            unphased_genotypes: 0,
            variable_ploidy: false,
        }
    }
}

#[derive(Debug)]
pub struct DiscoveryResult {
    pub samples: Vec<String>,
    pub contigs: Vec<String>,
    pub contig_lengths: HashMap<String, u64>,
    pub work_units: Vec<WorkUnit>,
}

pub fn discover(
    input: &Path,
    reference_contig_lengths: &HashMap<String, u64>,
    allowed_contigs: &[String],
    args: &Args,
) -> Result<DiscoveryResult> {
    let (samples, header_contigs) = read_header_metadata(input)
        .with_context(|| format!("could not read header of {}", input.display()))?;
    if samples.is_empty() {
        anyhow::bail!("input VCF contains no samples");
    }

    // A sample "has a BAM" if the directory contains an alignment file
    // whose name identifies it — allowing the common suffixes
    // (`.sorted.bam`, `.markdup.bam`, `.cram`, ...). See `bam_resolver`.
    let bams_present: Vec<bool> = samples
        .iter()
        .map(|name| {
            args.bam_dir
                .as_ref()
                .and_then(|d| crate::bam_resolver::find_sample_bam(d, name))
                .is_some()
        })
        .collect();

    let mut stats: HashMap<String, ContigAcc> = HashMap::new();
    for c in allowed_contigs {
        stats.insert(c.clone(), ContigAcc::new(samples.len()));
    }

    let reader = open_maybe_gzip(input)?;
    for line in std::io::BufRead::lines(reader) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let raw = match parse_raw_vcf_record(&line) {
            Ok(r) => r,
            Err(_) => continue,
        };
        let Some(acc) = stats.get_mut(raw.chrom) else {
            continue;
        };

        acc.variant_count += 1;

        let n_samples = samples.len();
        let orig = raw.samples.len();
        let mut fields = raw.samples.to_vec();
        if fields.len() < n_samples {
            fields.resize(n_samples, "");
        } else if fields.len() > n_samples {
            fields.truncate(n_samples);
        }

        for (idx, field) in fields.iter().enumerate() {
            if idx >= orig {
                continue;
            }
            let gt = match extract_gt(raw.format, field) {
                Ok(g) => g,
                Err(_) => continue,
            };
            let raw_gt = match parse_gt(gt) {
                Ok(g) => g,
                Err(_) => continue,
            };
            let norm = match normalize_genotype(&raw_gt, 2) {
                Ok(n) => n,
                Err(_) => continue,
            };
            let p = norm.len();
            if p > acc.max_ploidy[idx] {
                acc.max_ploidy[idx] = p;
            }
            match acc.first_ploidy[idx] {
                None => acc.first_ploidy[idx] = Some(p),
                Some(q) if q == p => {}
                Some(_) => acc.variable_ploidy = true,
            }
            if gt.contains('|') {
                acc.phased_genotypes += 1;
            } else if gt.contains('/') {
                acc.unphased_genotypes += 1;
            }
        }
    }

    let mut work_units = Vec::with_capacity(allowed_contigs.len());
    for (idx, contig) in allowed_contigs.iter().enumerate() {
        let mut acc = stats
            .remove(contig)
            .unwrap_or_else(|| ContigAcc::new(samples.len()));
        for p in &mut acc.max_ploidy {
            if *p == 0 {
                *p = 2;
            }
        }

        let ref_len = reference_contig_lengths.get(contig).copied().unwrap_or(0);
        let haplotype_count: usize = acc.max_ploidy.iter().sum();
        let max_ploidy = acc.max_ploidy.iter().copied().max().unwrap_or(0);

        let (phase_backend, needs_phasing) = classify_backend(
            samples.len(),
            &acc.max_ploidy,
            acc.variable_ploidy,
            acc.unphased_genotypes,
            max_ploidy,
            &bams_present,
        );

        work_units.push(WorkUnit {
            id: WorkUnitId(idx as u64),
            contig: contig.clone(),
            contig_index: idx,
            reference_length: ref_len,
            sample_count: samples.len(),
            sample_max_ploidies: acc.max_ploidy,
            haplotype_count,
            max_ploidy,
            variant_count: acc.variant_count,
            phased_genotypes: acc.phased_genotypes,
            unphased_genotypes: acc.unphased_genotypes,
            malformed_records: 0,
            needs_phasing,
            phase_backend,
            phased_vcf: None,
            phased_tempdir: None,
            state: WorkState::Analyzed,
        });
    }

    Ok(DiscoveryResult {
        samples,
        contigs: header_contigs,
        contig_lengths: reference_contig_lengths.clone(),
        work_units,
    })
}

/// Implements the automatic routing table. No user input required.
///
/// Precedence:
///   1. Nothing to do (all phased, or haploid-only) → `None`.
///   2. Single sample:
///        a. with BAM + variable ploidy → `VariablePloidy`
///        b. with BAM + constant ploidy → `WhatsHap`
///        c. no BAM                     → `SingleSampleCanonical`
///   3. Multi-sample with any sample's ploidy changing within the contig
///      → `VariablePloidy` (needs per-run phasing + merge).
///   4. Multi-sample uniform diploid → `Beagle`.
///   5. Multi-sample polyploid constant per sample → `WhatsHap`.
fn classify_backend(
    n_samples: usize,
    max_ploidies: &[usize],
    variable_ploidy: bool,
    unphased_genotypes: u64,
    max_ploidy_overall: usize,
    bams_present: &[bool],
) -> (PhaseBackend, bool) {
    if unphased_genotypes == 0 || max_ploidy_overall <= 1 {
        return (PhaseBackend::None, false);
    }

    if n_samples == 1 {
        let has_bam = bams_present.first().copied().unwrap_or(false);
        if has_bam {
            if variable_ploidy {
                return (PhaseBackend::VariablePloidy, true);
            }
            return (PhaseBackend::WhatsHap, true);
        }
        return (PhaseBackend::SingleSampleCanonical, false);
    }

    if variable_ploidy {
        return (PhaseBackend::VariablePloidy, true);
    }

    if max_ploidies.iter().all(|&p| p == 2) {
        return (PhaseBackend::Beagle, true);
    }

    (PhaseBackend::WhatsHap, true)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    fn write(path: &Path, content: &str) {
        let mut f = std::fs::File::create(path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
    }

    fn args() -> crate::cli::Args {
        crate::cli::Args {
            input: PathBuf::from("x.vcf"),
            reference: PathBuf::from("x.fa"),
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

    #[test]
    fn discovery_routes_uniform_diploid_to_beagle() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("in.vcf");
        write(
            &vcf,
            "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1
",
        );
        let mut lens = HashMap::new();
        lens.insert("chr1".into(), 100u64);
        let allowed = vec!["chr1".to_string()];
        let res = discover(&vcf, &lens, &allowed, &args()).unwrap();
        assert_eq!(res.work_units[0].phase_backend, PhaseBackend::Beagle);
        assert!(res.work_units[0].needs_phasing);
    }

    #[test]
    fn discovery_routes_polyploid_to_whatshap() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("in.vcf");
        write(
            &vcf,
            "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1/1\t0/0/1
chr1\t20\t.\tC\tT\t.\t.\t.\tGT\t1/0/0\t1/1/0
",
        );
        let mut lens = HashMap::new();
        lens.insert("chr1".into(), 100u64);
        let allowed = vec!["chr1".to_string()];
        let res = discover(&vcf, &lens, &allowed, &args()).unwrap();
        assert_eq!(res.work_units[0].phase_backend, PhaseBackend::WhatsHap);
    }

    #[test]
    fn discovery_routes_variable_ploidy_to_variable_pipeline() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("in.vcf");
        write(
            &vcf,
            "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1
chr1\t20\t.\tC\tT\t.\t.\t.\tGT\t0/1/1\t0/1
",
        );
        let mut lens = HashMap::new();
        lens.insert("chr1".into(), 100u64);
        let allowed = vec!["chr1".to_string()];
        let res = discover(&vcf, &lens, &allowed, &args()).unwrap();
        assert_eq!(
            res.work_units[0].phase_backend,
            PhaseBackend::VariablePloidy
        );
        assert!(res.work_units[0].needs_phasing);
    }

    #[test]
    fn discovery_single_sample_with_bam_and_variable_ploidy_routes_to_variable_pipeline() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("in.vcf");
        let bam_dir = dir.path().join("bams");
        std::fs::create_dir(&bam_dir).unwrap();
        std::fs::File::create(bam_dir.join("S1.bam")).unwrap();

        write(
            &vcf,
            "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1
chr1\t20\t.\tC\tT\t.\t.\t.\tGT\t0/1/1
",
        );
        let mut lens = HashMap::new();
        lens.insert("chr1".into(), 100u64);
        let allowed = vec!["chr1".to_string()];
        let mut a = args();
        a.bam_dir = Some(bam_dir);
        let res = discover(&vcf, &lens, &allowed, &a).unwrap();
        assert_eq!(
            res.work_units[0].phase_backend,
            PhaseBackend::VariablePloidy
        );
        assert!(res.work_units[0].needs_phasing);
    }

    #[test]
    fn discovery_single_sample_no_bam_with_variable_ploidy_falls_back_to_canonical() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("in.vcf");
        write(
            &vcf,
            "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1
chr1\t20\t.\tC\tT\t.\t.\t.\tGT\t0/1/1
",
        );
        let mut lens = HashMap::new();
        lens.insert("chr1".into(), 100u64);
        let allowed = vec!["chr1".to_string()];
        let res = discover(&vcf, &lens, &allowed, &args()).unwrap();
        assert_eq!(
            res.work_units[0].phase_backend,
            PhaseBackend::SingleSampleCanonical
        );
        assert!(!res.work_units[0].needs_phasing);
    }

    #[test]
    fn discovery_single_sample_no_bam_is_canonical() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("in.vcf");
        write(
            &vcf,
            "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1
",
        );
        let mut lens = HashMap::new();
        lens.insert("chr1".into(), 100u64);
        let allowed = vec!["chr1".to_string()];
        let res = discover(&vcf, &lens, &allowed, &args()).unwrap();
        assert_eq!(
            res.work_units[0].phase_backend,
            PhaseBackend::SingleSampleCanonical
        );
    }

    #[test]
    fn discovery_single_sample_with_suffixed_bam_routes_to_whatshap() {
        // Regression: BAMs named `S1.sorted.bam` must be found even though
        // the file is not named exactly `S1.bam`.
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("in.vcf");
        let bam_dir = dir.path().join("bams");
        std::fs::create_dir(&bam_dir).unwrap();
        std::fs::File::create(bam_dir.join("S1.sorted.markdup.bam")).unwrap();

        write(
            &vcf,
            "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1
",
        );
        let mut lens = HashMap::new();
        lens.insert("chr1".into(), 100u64);
        let allowed = vec!["chr1".to_string()];
        let mut a = args();
        a.bam_dir = Some(bam_dir);
        let res = discover(&vcf, &lens, &allowed, &a).unwrap();
        assert_eq!(res.work_units[0].phase_backend, PhaseBackend::WhatsHap);
    }

    #[test]
    fn discovery_fully_phased_needs_nothing() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("in.vcf");
        write(
            &vcf,
            "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|0
",
        );
        let mut lens = HashMap::new();
        lens.insert("chr1".into(), 100u64);
        let allowed = vec!["chr1".to_string()];
        let res = discover(&vcf, &lens, &allowed, &args()).unwrap();
        assert_eq!(res.work_units[0].phase_backend, PhaseBackend::None);
        assert!(!res.work_units[0].needs_phasing);
    }

    #[test]
    fn discovery_haploid_is_none() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("in.vcf");
        write(
            &vcf,
            "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t1
",
        );
        let mut lens = HashMap::new();
        lens.insert("chr1".into(), 100u64);
        let allowed = vec!["chr1".to_string()];
        let res = discover(&vcf, &lens, &allowed, &args()).unwrap();
        assert_eq!(res.work_units[0].phase_backend, PhaseBackend::None);
    }
}