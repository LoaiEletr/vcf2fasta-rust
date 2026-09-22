//! End-to-end tests for the pipeline.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
use vcf2fasta::workunit::{discover, PhaseBackend, WorkState};

fn write_vcf(path: &Path, body: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    f.write_all(body.as_bytes()).unwrap();
}

fn touch_reference(
    dir: &Path,
    contigs: &[(&str, u64)],
) -> (std::path::PathBuf, HashMap<String, u64>) {
    let mut lens = HashMap::new();
    for (c, l) in contigs {
        lens.insert(c.to_string(), *l);
    }
    let p = dir.join("ref.fa");
    write_vcf(&p, "");
    (p, lens)
}

fn make_args() -> vcf2fasta::cli::Args {
    vcf2fasta::cli::Args {
        input: std::path::PathBuf::from("x.vcf"),
        reference: std::path::PathBuf::from("x.fa"),
        prefix: String::new(),
        no_call_string: Some("N".into()),
        threads: 1,
        line_width: 80,
        no_validate_ref: false,
        quiet: true,
        gpu: false,
        device: vcf2fasta::scheduler::Device::Cpu,
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
fn discovery_produces_one_unit_per_contig() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    write_vcf(
        &vcf,
        "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##contig=<ID=chr3,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t0|1
chr2\t20\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1
chr3\t30\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1
",
    );
    let (_ref, lens) = touch_reference(
        dir.path(),
        &[("chr1", 1000), ("chr2", 1000), ("chr3", 1000)],
    );
    let allowed = vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()];
    let a = make_args();
    let res = discover(&vcf, &lens, &allowed, &a).unwrap();
    assert_eq!(res.work_units.len(), 3);
    assert!(!res.work_units[0].needs_phasing);
    assert!(res.work_units[1].needs_phasing);
    assert!(res.work_units[2].needs_phasing);
    assert_eq!(res.work_units[0].phase_backend, PhaseBackend::None);
    assert_eq!(res.work_units[1].phase_backend, PhaseBackend::Beagle);
}

#[test]
fn discovery_skips_contigs_not_in_allowed_list() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    write_vcf(
        &vcf,
        "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1
chr2\t20\t.\tA\tG\t.\t.\t.\tGT\t0/1
",
    );
    let (_ref, lens) = touch_reference(dir.path(), &[("chr1", 1000), ("chr2", 1000)]);
    let allowed = vec!["chr1".to_string()];
    let res = discover(&vcf, &lens, &allowed, &make_args()).unwrap();
    assert_eq!(res.work_units.len(), 1);
    assert_eq!(res.work_units[0].contig, "chr1");
}

#[test]
fn discovery_handles_haploid_contig() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    write_vcf(
        &vcf,
        "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t1
",
    );
    let (_ref, lens) = touch_reference(dir.path(), &[("chr1", 1000)]);
    let allowed = vec!["chr1".to_string()];
    let res = discover(&vcf, &lens, &allowed, &make_args()).unwrap();
    assert!(!res.work_units[0].needs_phasing);
}

#[test]
fn discovery_handles_polyploid_contig() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    write_vcf(
        &vcf,
        "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1/1\t0/1
",
    );
    let (_ref, lens) = touch_reference(dir.path(), &[("chr1", 1000)]);
    let allowed = vec!["chr1".to_string()];
    let res = discover(&vcf, &lens, &allowed, &make_args()).unwrap();
    assert!(res.work_units[0].needs_phasing);
    assert_eq!(res.work_units[0].phase_backend, PhaseBackend::WhatsHap);
    assert_eq!(res.work_units[0].sample_max_ploidies, vec![3, 2]);
    assert_eq!(res.work_units[0].haplotype_count, 5);
}

#[test]
fn workstate_failure_classification() {
    assert!(WorkState::PhasingFailed.is_terminal_failure());
    assert!(WorkState::Vcf2FastaFailed.is_terminal_failure());
    assert!(!WorkState::Ready.is_terminal_failure());
    assert!(WorkState::Complete.is_terminal_success());
}