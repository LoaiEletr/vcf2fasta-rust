//! Command-line interface definitions for vcf2fasta.

use crate::scheduler::Device;
use anyhow::{bail, Result};
use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Clone, Parser)]
#[command(
    name = "vcf2fasta",
    version,
    author,
    about = "Generate per-sample, per-haplotype FASTA files from an indexed VCF",
    long_about = "A fast, parallel tool to reconstruct haplotype sequences from a VCF.\n\n\
                  For each sample and haplotype, a FASTA file is created containing the complete \
                  contig sequence with variants applied.",
    after_help = concat!(
        "EXAMPLES:\n",
        "  vcf2fasta --reference ref.fa input.vcf.gz --prefix test/output\n",
        "  vcf2fasta -f ref.fa -p test/output --threads 8 input.vcf.gz\n",
        "  vcf2fasta -f ref.fa --device gpu input.vcf.gz\n",
        "  vcf2fasta -f ref.fa --bam-dir /bams input.vcf.gz\n",
        "  vcf2fasta -f ref.fa --beagle-ref-panel /ref --beagle-genetic-map /maps input.vcf.gz\n",
        "  vcf2fasta -f ref.fa --merged-output --prefix out input.vcf.gz\n\n",
        "Copyright (c) 2026 ",
        env!("CARGO_PKG_AUTHORS"),
        ". All rights reserved."
    )
)]
pub struct Args {
    /// Input VCF file. Must be bgzipped and indexed
    /// with a companion .tbi or .csi file.
    #[arg(value_name = "VCF")]
    pub input: PathBuf,

    /// Indexed reference FASTA. Accepts plain (`.fa`, `.fasta`) or
	/// BGZF-compressed (`.fa.gz`, `.fasta.gz`) input.
	///
	/// * Plain FASTA requires a `.fai` index (`samtools faidx ref.fa`).
	/// * BGZF-compressed FASTA requires **both** a `.fai` and a `.gzi`
	///   index. `samtools faidx ref.fa.gz` produces both when the file is
	///   BGZF. A plain-gzip file (produced by `gzip`) is not BGZF and will
	///   be rejected by the underlying htslib reader.
	///
	/// Contig names must match those in the VCF header.
    #[arg(short = 'f', long = "reference", value_name = "FASTA")]
    pub reference: PathBuf,

    /// Output path prefix. Each generated FASTA is written as
    /// `<PREFIX><sample>_<contig>:<hap>.fa`. The log file is
    /// `<PREFIX>.log` and the warnings log is `<PREFIX>.warnings.log`.
    /// If empty (the default), outputs are written into the current
    /// directory and the log is `vcf2fasta.log`.
    #[arg(short = 'p', long = "prefix", default_value = "")]
    pub prefix: String,

    /// Character used for missing or invalid alleles. Single-character
    /// strings are recommended; longer strings are permitted but will
    /// extend the output sequence by their length.
    #[arg(short = 'n', long = "no-call-string", default_value = "N")]
    pub no_call_string: Option<String>,

    /// Number of CPU threads available to the pipeline. Sets the Rayon
    /// pool size on the CPU path and the producer thread count on the
    /// GPU path. Streams on the GPU path are further capped by this
    /// value, since each stream is a host-side thread.
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    pub threads: usize,

    /// Bases per FASTA line. Zero is rejected. Typical values are 60
    /// and 80.
    #[arg(short = 'w', long = "line-width", default_value_t = 80)]
    pub line_width: usize,

    /// Skip REF-vs-FASTA validation. Disabling this check speeds up the
    /// parse pass but will silently apply variants at positions where
    /// the VCF REF does not match the reference base.
    #[arg(short = 'v', long = "no-validate-ref", default_value_t = false)]
    pub no_validate_ref: bool,

    /// Suppress stage log lines on stderr. The log file is still written.
    #[arg(short = 'q', long = "quiet", default_value_t = false)]
    pub quiet: bool,

    /// Execution device: `auto`, `cpu`, or `gpu`.
    /// `auto` picks GPU only when the workload is large enough to
    /// amortise the transfer overhead; otherwise CPU.
    #[arg(long = "device", value_enum, default_value_t = Device::Auto)]
    pub device: Device,

    /// Deprecated alias for `--device gpu`. Hidden from --help.
    #[arg(short = 'g', long = "gpu", default_value_t = false, hide = true)]
    pub gpu: bool,

    /// Directory containing one coordinate-sorted, indexed BAM (or CRAM)
    /// per sample. Only consulted when read-based phasing (WhatsHap or
    /// the variable-ploidy pipeline) is required for some contig.
    ///
    /// File naming: for a sample named `SAMPLE001`, any file whose name
    /// starts with `SAMPLE001` as a whole token is accepted. So
    /// `SAMPLE001.bam`, `SAMPLE001.sorted.bam`, `SAMPLE001.markdup.bam`,
    /// and `SAMPLE001.cram` all work. `SAMPLE0010.bam` is not matched
    /// against `SAMPLE001`. If multiple files match, the shortest name
    /// wins, with an alphabetical tiebreak.
    #[arg(long = "bam-dir", value_name = "DIR")]
    pub bam_dir: Option<PathBuf>,

    /// Directory containing Beagle reference-panel files.
    ///
    /// For each chromosome phased with Beagle, a file is looked up whose
    /// name contains the chromosome as a whole token (e.g. `chr1` or `1`
    /// appearing in `reference_chr1_asasa`, `1000G.chr1.ref.vcf.gz`,
    /// `chr1.ref`, or the double-prefix artifact `plink.chrchr1.GRCh38.map`).
    /// Matching is case-insensitive. The match is passed to Beagle as
    /// `ref=<file>`. If no file matches, `ref=` is omitted.
    #[arg(long = "beagle-ref-panel", value_name = "DIR")]
    pub beagle_ref_panel: Option<PathBuf>,

    /// Directory containing Beagle genetic-map files.
    ///
    /// Same matching rule as `--beagle-ref-panel`; the match is passed to
    /// Beagle as `map=<file>`. If no file matches for a chromosome,
    /// `map=` is omitted. Independent of `--beagle-ref-panel`.
    #[arg(long = "beagle-genetic-map", value_name = "DIR")]
    pub beagle_genetic_map: Option<PathBuf>,

    /// Override the automatic chunk size for a work unit, in base pairs.
    /// Larger chunks reduce per-chunk overhead but increase peak memory;
    /// smaller chunks give the scheduler more parallelism at the cost of
    /// more per-chunk setup. By default the chunk size is derived from
    /// `MAX_TILE_OUTPUT_BYTES` and the cohort size.
    #[arg(long = "chunk-size", value_name = "BP")]
    pub chunk_size: Option<u64>,

    /// Look-back distance (bp) used to detect variants from a previous
    /// chunk whose REF allele extends into the current chunk. Must be at
    /// least as large as the longest REF allele in the VCF. By default
    /// this is auto-detected by scanning the VCF.
    #[arg(long = "chunk-pad", value_name = "BP")]
    pub chunk_pad: Option<u64>,

    /// Threads used by each phasing subprocess (Beagle / WhatsHap).
	/// When unset, the budget's auto-derived value is used: roughly half
	/// of the machine's physical cores (the other half is reserved for
	/// the vcf2fasta stage).
    #[arg(long = "phasing-threads", value_name = "N")]
    pub phasing_threads: Option<usize>,

    /// Number of vcf2fasta consumer threads. 1 (default) processes one
    /// contig at a time; higher values parallelise across contigs on CPU.
    /// GPU mode forces 1.
    #[arg(long = "vcf2fasta-workers", value_name = "N", default_value_t = 1)]
    pub vcf2fasta_workers: usize,

    /// Depth of the bounded READY queue between phasing and vcf2fasta.
    /// Larger values give the phasing producer more slack; smaller values
    /// apply backpressure sooner.
    #[arg(long = "ready-queue-depth", value_name = "N", default_value_t = 4)]
    pub ready_queue_depth: usize,

    /// Disable the pipelined scheduler. Phasing runs to completion for
    /// the whole file before vcf2fasta starts, and GPU tiling is
    /// disabled. Useful for debugging or for bisecting regressions
    /// against the pre-pipeline implementation.
    #[arg(long = "no-pipeline", default_value_t = false)]
    pub no_pipeline: bool,

    /// Produce one merged FASTA per (sample, haplotype) containing all
    /// contigs in reference order, instead of one file per contig.
    /// Output name: `<PREFIX><sample>_<hap>.fa`.
    #[arg(long = "merged-output", default_value_t = false)]
    pub merged_output: bool,

    /// Soft ceiling on process RAM. Accepts K/M/G/T suffixes (e.g. `32G`,
    /// `512M`). The value is clamped to a safe fraction of the machine's
    /// total RAM; if the requested value is larger, the effective value
    /// and the reduction are logged.
    #[arg(long = "max-memory", value_name = "SIZE")]
    pub max_memory: Option<String>,

    /// Soft ceiling on GPU VRAM. Accepts K/M/G/T suffixes. Used by the
    /// scheduler when choosing tile sizes and stream counts. If the
    /// requested value exceeds available VRAM, the effective value is
    /// reduced.
    #[arg(long = "max-vram", value_name = "SIZE")]
    pub max_vram: Option<String>,

    /// Comma-separated CUDA device indices to use, e.g. `0,1`. Defaults
    /// to device 0 when `--device gpu` is set. When `--device auto` is
    /// used, this selects the candidate devices for the scheduler.
    #[arg(long = "gpu-devices", value_name = "LIST")]
    pub gpu_devices: Option<String>,
}

pub fn resolve_threads(args: &Args) -> usize {
    args.threads.max(1)
}

pub fn resolve_phasing_threads(args: &Args) -> usize {
    args.phasing_threads
        .unwrap_or_else(|| resolve_threads(args))
        .max(1)
}

/// Merges the deprecated `--gpu` flag into the new `--device` selection.
pub fn effective_device(args: &Args) -> Device {
    if args.gpu {
        Device::Gpu
    } else {
        args.device
    }
}

/// Parse "32G", "512M", "1024K", "1048576" into bytes.
pub fn parse_size(s: &str) -> Result<u64> {
    let s = s.trim();
    if s.is_empty() {
        bail!("empty size string");
    }
    let last = s.chars().last().unwrap();
    let (num, mult) = match last.to_ascii_uppercase() {
        'K' => (&s[..s.len() - 1], 1024u64),
        'M' => (&s[..s.len() - 1], 1024 * 1024),
        'G' => (&s[..s.len() - 1], 1024 * 1024 * 1024),
        'T' => (&s[..s.len() - 1], 1024u64 * 1024 * 1024 * 1024),
        _ => (s, 1),
    };
    let n: u64 = num
        .trim()
        .parse()
        .map_err(|e| anyhow::anyhow!("invalid size '{}': {}", s, e))?;
    Ok(n.saturating_mul(mult))
}

/// Parse "0,1,2" into a list of GPU indices.
pub fn parse_gpu_list(s: &str) -> Result<Vec<usize>> {
    let mut out = Vec::new();
    for part in s.split(',') {
        let p = part.trim();
        if p.is_empty() {
            continue;
        }
        let idx: usize = p
            .parse()
            .map_err(|e| anyhow::anyhow!("invalid GPU index '{}': {}", p, e))?;
        out.push(idx);
    }
    if out.is_empty() {
        bail!("--gpu-devices is empty");
    }
    Ok(out)
}

pub fn validate_cli(args: &Args) -> Result<()> {
    if args.threads == 0 {
        bail!("--threads must be greater than 0");
    }
    if args.line_width == 0 {
        bail!("--line-width must be greater than 0");
    }
    if let Some(0) = args.chunk_size {
        bail!("--chunk-size must be greater than 0");
    }
    if let Some(0) = args.chunk_pad {
        bail!("--chunk-pad must be greater than 0");
    }
    if args.ready_queue_depth == 0 {
        bail!("--ready-queue-depth must be greater than 0");
    }
    if args.vcf2fasta_workers == 0 {
        bail!("--vcf2fasta-workers must be greater than 0");
    }
    if let Some(s) = &args.max_memory {
        let _ = parse_size(s)?;
    }
    if let Some(s) = &args.max_vram {
        let _ = parse_size(s)?;
    }
    if let Some(s) = &args.gpu_devices {
        let _ = parse_gpu_list(s)?;
    }
    if let Some(d) = &args.beagle_ref_panel {
        if !d.is_dir() {
            bail!(
                "--beagle-ref-panel {} does not exist or is not a directory",
                d.display()
            );
        }
    }
    if let Some(d) = &args.beagle_genetic_map {
        if !d.is_dir() {
            bail!(
                "--beagle-genetic-map {} does not exist or is not a directory",
                d.display()
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_args() -> Args {
        Args {
            input: PathBuf::from("test.vcf"),
            reference: PathBuf::from("test.fa"),
            prefix: String::new(),
            no_call_string: None,
            threads: 1,
            line_width: 80,
            no_validate_ref: false,
            quiet: false,
            gpu: false,
            device: Device::Auto,
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
    fn validate_cli_ok_with_defaults() {
        assert!(validate_cli(&base_args()).is_ok());
    }
    #[test]
    fn validate_cli_fails_threads_zero() {
        let mut a = base_args();
        a.threads = 0;
        assert!(validate_cli(&a)
            .unwrap_err()
            .to_string()
            .contains("--threads"));
    }
    #[test]
    fn resolve_threads_returns_user_value() {
        let mut a = base_args();
        a.threads = 7;
        assert_eq!(resolve_threads(&a), 7);
    }
    #[test]
    fn resolve_phasing_threads_falls_back() {
        let mut a = base_args();
        a.threads = 4;
        assert_eq!(resolve_phasing_threads(&a), 4);
        a.phasing_threads = Some(2);
        assert_eq!(resolve_phasing_threads(&a), 2);
    }
    #[test]
    fn parse_size_variants() {
        assert_eq!(parse_size("1024").unwrap(), 1024);
        assert_eq!(parse_size("2K").unwrap(), 2048);
        assert_eq!(parse_size("3M").unwrap(), 3 * 1024 * 1024);
        assert_eq!(parse_size("4G").unwrap(), 4 * 1024 * 1024 * 1024);
    }
    #[test]
    fn parse_size_rejects_garbage() {
        assert!(parse_size("abc").is_err());
        assert!(parse_size("").is_err());
    }
    #[test]
    fn parse_gpu_list_ok() {
        assert_eq!(parse_gpu_list("0,1,2").unwrap(), vec![0, 1, 2]);
        assert_eq!(parse_gpu_list("0").unwrap(), vec![0]);
    }
    #[test]
    fn parse_gpu_list_rejects_empty() {
        assert!(parse_gpu_list("").is_err());
        assert!(parse_gpu_list(",,").is_err());
    }
    #[test]
    fn gpu_flag_is_alias_for_device_gpu() {
        let mut a = base_args();
        a.gpu = true;
        assert_eq!(effective_device(&a), Device::Gpu);
        a.gpu = false;
        a.device = Device::Cpu;
        assert_eq!(effective_device(&a), Device::Cpu);
    }
    #[test]
    fn validate_cli_rejects_missing_beagle_ref_panel_dir() {
        let mut a = base_args();
        a.beagle_ref_panel = Some(PathBuf::from("/definitely/not/a/real/path/xyz"));
        let e = validate_cli(&a).unwrap_err().to_string();
        assert!(e.contains("--beagle-ref-panel"), "got: {}", e);
    }
    #[test]
    fn validate_cli_rejects_missing_beagle_genetic_map_dir() {
        let mut a = base_args();
        a.beagle_genetic_map = Some(PathBuf::from("/definitely/not/a/real/path/xyz"));
        let e = validate_cli(&a).unwrap_err().to_string();
        assert!(e.contains("--beagle-genetic-map"), "got: {}", e);
    }
    #[test]
    fn validate_cli_accepts_existing_bam_dir() {
        let dir = tempfile::tempdir().unwrap();
        let mut a = base_args();
        a.bam_dir = Some(dir.path().to_path_buf());
        assert!(validate_cli(&a).is_ok());
    }
}