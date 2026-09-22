//! Streaming workload profiling.

use crate::phasing::open_maybe_gzip;
use crate::scheduler::WorkloadProfile;
use crate::vcf::parse_raw_vcf_record;
use anyhow::{Context, Result};
use std::io::BufRead;
use std::path::Path;

pub fn profile_vcf(
    path: &Path,
    sample_names: &[String],
    reference_bases: u64,
    contig_count: usize,
    line_width: usize,
) -> Result<WorkloadProfile> {
    let mut profile = WorkloadProfile {
        sample_count: sample_names.len(),
        reference_bases,
        contig_count,
        ..Default::default()
    };

    let mut max_ploidies = vec![0usize; sample_names.len()];
    let reader = open_maybe_gzip(path)
        .with_context(|| format!("could not open {} for profiling", path.display()))?;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        let raw = match parse_raw_vcf_record(&line) {
            Ok(r) => r,
            Err(_) => {
                profile.malformed_records += 1;
                continue;
            }
        };

        profile.variant_count += 1;

        let snv = raw.ref_allele.len() == 1
            && !raw.alt_field.is_empty()
            && raw.alt_field != "."
            && raw.alt_field.split(',').all(|a| a.len() == 1);
        if snv {
            profile.snv_count += 1;
        } else {
            profile.indel_count += 1;
        }

        for (idx, _) in sample_names.iter().enumerate() {
            let sample_field = match raw.samples.get(idx) {
                Some(s) => *s,
                None => continue,
            };
            let gt_text = match extract_gt_field(raw.format, sample_field) {
                Some(g) => g,
                None => continue,
            };
            let ploidy = count_alleles(gt_text);
            if ploidy > max_ploidies[idx] {
                max_ploidies[idx] = ploidy;
            }
            if gt_text.contains('|') {
                profile.phased_genotypes += 1;
            } else if gt_text.contains('/') {
                profile.unphased_genotypes += 1;
            }
        }
    }

    for p in &mut max_ploidies {
        if *p == 0 {
            *p = 2;
        }
    }
    profile.sample_max_ploidies = max_ploidies.clone();
    profile.haplotype_count = max_ploidies.iter().sum();
    profile.max_ploidy = max_ploidies.iter().copied().max().unwrap_or(0);
    profile.estimated_output_bytes = WorkloadProfile::estimate_output_bytes(
        profile.haplotype_count,
        reference_bases,
        line_width,
    );
    Ok(profile)
}

fn extract_gt_field<'a>(format: &str, sample: &'a str) -> Option<&'a str> {
    let fmt: Vec<&str> = format.split(':').collect();
    let gi = fmt.iter().position(|f| *f == "GT")?;
    let smp: Vec<&str> = sample.split(':').collect();
    smp.get(gi).copied()
}

fn count_alleles(gt: &str) -> usize {
    if gt.is_empty() {
        return 0;
    }
    let mut n = 1usize;
    for b in gt.bytes() {
        if b == b'/' || b == b'|' {
            n += 1;
        }
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    fn write(path: &Path, content: &str) {
        let mut f = File::create(path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
    }

    const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t0/1
chr1\t20\t.\tC\tT\t.\t.\t.\tGT\t1|0\t0/0
chr1\t30\t.\tG\tGA\t.\t.\t.\tGT\t0|1\t0/1
";

    #[test]
    fn counts_variants_and_phasing() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("v.vcf");
        write(&p, VCF);
        let samples = vec!["S1".to_string(), "S2".to_string()];
        let prof = profile_vcf(&p, &samples, 1000, 1, 80).unwrap();
        assert_eq!(prof.variant_count, 3);
        assert_eq!(prof.snv_count, 2);
        assert_eq!(prof.indel_count, 1);
        assert_eq!(prof.haplotype_count, 4);
        assert_eq!(prof.sample_max_ploidies, vec![2, 2]);
    }

    #[test]
    fn counts_malformed_records() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("m.vcf");
        let content = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\tnot_a_number\t.\tA\tG\t.\t.\t.\tGT\t0/1
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1
";
        write(&p, content);
        let samples = vec!["S1".to_string()];
        let prof = profile_vcf(&p, &samples, 1000, 1, 80).unwrap();
        assert_eq!(prof.malformed_records, 1);
        assert_eq!(prof.variant_count, 1);
    }

    #[test]
    fn ploidy_summed_per_sample_not_times_max() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("p.vcf");
        let content = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1/1
";
        write(&p, content);
        let samples = vec!["S1".to_string(), "S2".to_string()];
        let prof = profile_vcf(&p, &samples, 1000, 1, 80).unwrap();
        assert_eq!(prof.haplotype_count, 5);
        assert_eq!(prof.max_ploidy, 3);
        assert_eq!(prof.sample_max_ploidies, vec![2, 3]);
    }
}