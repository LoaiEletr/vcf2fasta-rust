//! Integration tests for pipeline-level VCF2FASTA policies.
//!
//! These policies cannot be tested at the `decode_and_validate_record`
//! level because they live in the tile executor and depend on the
//! streamed, tabix-fetched view of the VCF:
//!
//! * **Policy 14** — Unsorted VCF must be rejected, not silently processed.
//!   Enforced by the requirement that the input be tabix-indexed, and
//!   tabix refuses to index unsorted files.
//!
//! * **Policy 26** — Duplicate position: keep the first, skip the second,
//!   emit a warning.
//!
//! * **Policy 27** — Overlapping variant (REF interval intersects a
//!   previously-applied one): keep the first, skip the second, emit a
//!   warning.
//!
//! ## Dependencies
//!
//! `bgzip` and `tabix` must be on `PATH`. Tests print a skip notice on
//! stderr and return early if either is missing, so the test binary does
//! not fail on systems without htslib tools.
//!
//! ## Fixture
//!
//! A 100-bp reference `ACGT` repeated 25 times:
//!
//! ```text
//! POS:   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 ...
//! idx:   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 ...
//! base:  A  C  G  T  A  C  G  T  A  C  G  T  A  C  G  T ...
//! ```
//!
//! Every test uses POS=10 (0-based 9, base `C`) or POS=10..11
//! (0-based 9..11, bases `CG`), so REF fields can be written to match
//! the reference exactly and pass validation.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use vcf2fasta::chunk::Tile;
use vcf2fasta::cli::Args;
use vcf2fasta::scheduler::Device;
use vcf2fasta::tile_executor::execute_tile;

// ---------------------------------------------------------------------------
// Fixture helpers
// ---------------------------------------------------------------------------

const TEST_SEQ: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

fn tools_available() -> bool {
    which("bgzip").is_some() && which("tabix").is_some()
}

fn which(name: &str) -> Option<PathBuf> {
    let path_var = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path_var) {
        let candidate = dir.join(name);
        if candidate.is_file() {
            return Some(candidate);
        }
    }
    None
}

/// Write a 100-bp reference FASTA and its `.fai` index. Returns the FASTA
/// path.
fn write_reference(dir: &Path) -> PathBuf {
    let fasta = dir.join("ref.fa");
    let fai = dir.join("ref.fa.fai");

    let mut f = fs::File::create(&fasta).unwrap();
    writeln!(f, ">chr1").unwrap();
    writeln!(f, "{}", TEST_SEQ).unwrap();
    drop(f);

    // .fai format: name \t length \t offset \t linebases \t linewidth
    // offset = bytes before the first base = len(">chr1\n") = 6
    // linebases = 100 (all bases on one line)
    // linewidth = 101 (100 bases + 1 newline)
    let offset = ">chr1\n".len();
    let mut f = fs::File::create(&fai).unwrap();
    writeln!(
        f,
        "chr1\t{}\t{}\t{}\t{}",
        TEST_SEQ.len(),
        offset,
        TEST_SEQ.len(),
        TEST_SEQ.len() + 1
    )
    .unwrap();

    fasta
}

/// Write a VCF, bgzip it, and tabix-index it. Returns the `.vcf.gz` path.
fn write_bgzipped_vcf(dir: &Path, body: &str) -> PathBuf {
    let plain = dir.join("input.vcf");
    fs::write(&plain, body).unwrap();

    let gz = dir.join("input.vcf.gz");
    let out = Command::new("bgzip")
        .arg("-c")
        .arg(&plain)
        .output()
        .expect("bgzip failed to spawn");
    assert!(
        out.status.success(),
        "bgzip failed: status={:?} stderr={:?}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    fs::write(&gz, &out.stdout).unwrap();

    let status = Command::new("tabix")
        .arg("-p")
        .arg("vcf")
        .arg(&gz)
        .status()
        .expect("tabix failed to spawn");
    assert!(status.success(), "tabix failed on {}", gz.display());

    gz
}

/// Build a minimal `Args` suitable for a single-sample CPU tile run.
fn make_args(input: PathBuf, reference: PathBuf) -> Args {
    Args {
        input,
        reference,
        prefix: String::new(),
        no_call_string: Some("N".to_string()),
        threads: 1,
        line_width: 80,
        no_validate_ref: false,
        quiet: true,
        gpu: false,
        device: Device::Cpu,
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

/// A tile covering the whole 100-bp chr1 contig for a single diploid sample.
fn whole_contig_tile() -> Tile {
    Tile {
        contig: "chr1".to_string(),
        hap_start: 0,
        hap_end: 2,
        base_start: 0,
        base_end: TEST_SEQ.len() as u64,
    }
}

/// Run `execute_tile` on a single-sample diploid tile.
fn run_tile(vcf: PathBuf, reference: PathBuf) -> anyhow::Result<vcf2fasta::chunk::TileResult> {
    let args = make_args(vcf, reference);
    let sample_names = vec!["S1".to_string()];
    let sample_hap_offset = vec![0usize];
    let max_ploidies = vec![2usize];
    execute_tile(
        &whole_contig_tile(),
        &sample_names,
        &sample_hap_offset,
        &max_ploidies,
        10_000,
        &args,
    )
}

// ---------------------------------------------------------------------------
// Policy 14: Unsorted VCF must be rejected
// ---------------------------------------------------------------------------

/// Tabix refuses to index a VCF whose records are out of position order.
/// That is what enforces policy 14 at the tool level: the tool always goes
/// through a tabix index, so an unsorted VCF never reaches the executor.
#[test]
fn unsorted_vcf_cannot_be_indexed() {
    if !tools_available() {
        eprintln!("skipping unsorted_vcf_cannot_be_indexed: bgzip/tabix not on PATH");
        return;
    }
    let dir = tempfile::tempdir().unwrap();

    // Two records deliberately out of position order.
    let vcf_body = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t50\t.\tA\tG\t.\t.\t.\tGT\t0|1
chr1\t10\t.\tC\tG\t.\t.\t.\tGT\t0|1
";
    let plain = dir.path().join("unsorted.vcf");
    fs::write(&plain, vcf_body).unwrap();

    // bgzip succeeds: it does not care about record order.
    let gz = dir.path().join("unsorted.vcf.gz");
    let out = Command::new("bgzip").arg("-c").arg(&plain).output().unwrap();
    assert!(out.status.success());
    fs::write(&gz, &out.stdout).unwrap();

    // tabix MUST fail: it verifies that positions are monotonically
    // non-decreasing within each contig.
    let status = Command::new("tabix")
        .arg("-p")
        .arg("vcf")
        .arg(&gz)
        .status()
        .expect("tabix failed to spawn");
    assert!(
        !status.success(),
        "tabix must reject an unsorted VCF; policy 14 depends on this"
    );
}

/// The tool must also reject a VCF that has no tabix index at all. This
/// covers the case where a user passes a plain `.vcf` or a `.vcf.gz`
/// without a `.tbi`.
#[test]
fn unindexed_vcf_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let ref_path = write_reference(dir.path());

    // A plain (non-bgzipped) VCF that looks syntactically valid.
    let vcf = dir.path().join("plain.vcf");
    fs::write(
        &vcf,
        "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tC\tG\t.\t.\t.\tGT\t0|1
",
    )
    .unwrap();

    // No index exists; the tile executor's tabix reader must fail.
    let result = run_tile(vcf, ref_path);
    assert!(
        result.is_err(),
        "tool must reject a VCF without a tabix index"
    );
    let msg = format!("{:#}", result.err().unwrap());
    assert!(
        msg.contains("tabix") || msg.contains("index") || msg.contains("could not open"),
        "unexpected error message: {msg}"
    );
}

// ---------------------------------------------------------------------------
// Policy 26: Duplicate position — keep the first, skip the second
// ---------------------------------------------------------------------------

#[test]
fn duplicate_position_keeps_first_and_warns() {
    if !tools_available() {
        eprintln!("skipping duplicate_position_keeps_first_and_warns: bgzip/tabix not on PATH");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let ref_path = write_reference(dir.path());

    // Two records at POS=10 (0-based 9, base 'C').
    //
    //   first  record: REF=C  ALT=G  GT=0|1  → hap0=C, hap1=G
    //   second record: REF=C  ALT=T  GT=1|0  → skipped
    //
    // Expected per-hap byte at 0-based 9:
    //   hap0: C  (from first record's REF)
    //   hap1: G  (from first record's ALT)
    let vcf_body = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tC\tG\t.\t.\t.\tGT\t0|1
chr1\t10\t.\tC\tT\t.\t.\t.\tGT\t1|0
";
    let vcf = write_bgzipped_vcf(dir.path(), vcf_body);

    let result = run_tile(vcf, ref_path).expect("execute_tile failed");

    assert_eq!(result.seen, 2, "both records must be counted as seen");
    assert_eq!(result.applied, 1, "only the first record must be applied");

    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("overlap") || w.contains("out-of-order")),
        "expected an overlap warning; got: {:?}",
        result.warnings
    );

    assert_eq!(result.per_hap[0].len(), TEST_SEQ.len());
    assert_eq!(result.per_hap[1].len(), TEST_SEQ.len());
    assert_eq!(result.per_hap[0][9], b'C', "hap0 must keep first record's REF");
    assert_eq!(result.per_hap[1][9], b'G', "hap1 must have first record's ALT");
}

// ---------------------------------------------------------------------------
// Policy 27: Overlapping variant — keep the first, skip the second
// ---------------------------------------------------------------------------

#[test]
fn overlapping_variant_keeps_first_and_warns() {
    if !tools_available() {
        eprintln!("skipping overlapping_variant_keeps_first_and_warns: bgzip/tabix not on PATH");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let ref_path = write_reference(dir.path());

    // First record: POS=10 (0-based 9), REF="CG" (bases 9..11), ALT="A".
    // Second record: POS=11 (0-based 10), REF="G". Its interval [10, 11)
    // overlaps the first record's applied interval [9, 11) and must be
    // skipped.
    //
    // Expected per-hap output:
    //   hap0 (GT allele 0 = REF "CG"): unchanged → 100 bytes
    //   hap1 (GT allele 1 = ALT "A"):  deletes one base → 99 bytes
    let vcf_body = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tCG\tA\t.\t.\t.\tGT\t0|1
chr1\t11\t.\tG\tT\t.\t.\t.\tGT\t1|0
";
    let vcf = write_bgzipped_vcf(dir.path(), vcf_body);

    let result = run_tile(vcf, ref_path).expect("execute_tile failed");

    assert_eq!(result.seen, 2);
    assert_eq!(result.applied, 1);

    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("overlap") || w.contains("out-of-order")),
        "expected an overlap warning; got: {:?}",
        result.warnings
    );

    // hap0 unchanged, hap1 lost one base to the deletion.
    assert_eq!(result.per_hap[0].len(), 100);
    assert_eq!(result.per_hap[1].len(), 99);

    // Bases at 0-based 9 and 10 in hap0 must match the reference.
    assert_eq!(result.per_hap[0][9], b'C');
    assert_eq!(result.per_hap[0][10], b'G');

    // Byte at 0-based 9 in hap1 must be the ALT (single 'A').
    assert_eq!(result.per_hap[1][9], b'A');

    // The trailing 89 bytes of hap1 must match reference[11..100].
    assert_eq!(&result.per_hap[1][10..], &TEST_SEQ.as_bytes()[11..]);
}

// ---------------------------------------------------------------------------
// Sanity: a well-formed sorted VCF is accepted
// ---------------------------------------------------------------------------

/// Positive control for policies 26 and 27: a sorted VCF with no
/// duplicates and no overlaps must be applied in full.
///
/// Reference layout reminder (0-based index → base):
///   8 → A, 9 → C, 10 → G, 11 → T, 12 → A, 13 → C, 14 → G, 15 → T
#[test]
fn baseline_sorted_vcf_applies_all_records() {
    if !tools_available() {
        eprintln!("skipping baseline_sorted_vcf_applies_all_records: bgzip/tabix not on PATH");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let ref_path = write_reference(dir.path());

    // Record 1: POS=10 (0-based 9, base 'C'), REF='C' ✓, ALT='G',
    //           GT=0|1 → hap0='C', hap1='G'.
    // Record 2: POS=14 (0-based 13, base 'C'), REF='C' ✓, ALT='T',
    //           GT=1|0 → hap0='T', hap1='C'.
    let vcf_body = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tC\tG\t.\t.\t.\tGT\t0|1
chr1\t14\t.\tC\tT\t.\t.\t.\tGT\t1|0
";
    let vcf = write_bgzipped_vcf(dir.path(), vcf_body);

    let result = run_tile(vcf, ref_path).expect("execute_tile failed");

    assert_eq!(result.seen, 2);
    assert_eq!(result.applied, 2);
    assert_eq!(
        result.warnings.len(),
        0,
        "no warnings expected on a clean VCF; got: {:?}",
        result.warnings
    );

    // hap0: pos9=C (REF from record 1), pos13=T (ALT from record 2)
    // hap1: pos9=G (ALT from record 1), pos13=C (REF from record 2)
    assert_eq!(result.per_hap[0][9], b'C');
    assert_eq!(result.per_hap[0][13], b'T');
    assert_eq!(result.per_hap[1][9], b'G');
    assert_eq!(result.per_hap[1][13], b'C');
}