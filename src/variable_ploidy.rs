//! Variable-ploidy phasing.
//!
//! When a single sample has different ploidies at different positions
//! (e.g. a triploid region followed by a diploid region), no phaser can
//! process the sample in one pass. This module partitions each sample's
//! genome into maximal runs of constant ploidy, extracts a single-sample
//! mini-VCF per run, phases each with WhatsHap at the run's ploidy, and
//! merges the results back into a single VCF.
//!
//! The merge is a k-way streaming merge: memory is O(number of runs),
//! NOT O(variants × samples). This is what lets whole-genome variable-
//! ploidy inputs finish without OOM.

use crate::cli::Args;
use crate::phasing::{
    bgzip_file, open_maybe_gzip, pre_filter_vcf, run_whatshap_polyphase, tabix_index,
};
use crate::vcf::parse_raw_vcf_record;
use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct PloidyRun {
    pub sample_idx: usize,
    pub sample_name: String,
    pub contig: String,
    pub start_pos: usize,
    pub end_pos: usize,
    pub ploidy: usize,
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

pub fn has_variable_ploidy(input: &Path, sample_names: &[String]) -> Result<bool> {
    let mut first_ploidy: Vec<Option<usize>> = vec![None; sample_names.len()];
    let reader = open_maybe_gzip(input)?;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let raw = match parse_raw_vcf_record(&line) {
            Ok(r) => r,
            Err(_) => continue,
        };

        for (idx, _) in sample_names.iter().enumerate() {
            if idx >= raw.samples.len() {
                continue;
            }
            let ploidy = match ploidy_of_sample_field(raw.format, raw.samples[idx]) {
                Some(p) if p > 0 => p,
                _ => continue,
            };
            match first_ploidy[idx] {
                None => first_ploidy[idx] = Some(ploidy),
                Some(p) if p == ploidy => {}
                Some(_) => return Ok(true),
            }
        }
    }
    Ok(false)
}

pub fn run_variable_ploidy_pipeline(
    input: &Path,
    reference: &Path,
    args: &Args,
    sample_names: &[String],
    threads: usize,
) -> Result<PathBuf> {
    let bam_dir = match &args.bam_dir {
        Some(d) => d,
        None => bail!(
            "This VCF contains at least one sample whose ploidy changes \
             across the contig.\n\
             No phaser can process such a sample in a single run. Please \
             supply a directory of BAMs via --bam-dir.\n\
             Each BAM must be named `<SAMPLE>.bam` (or `<SAMPLE>.<suffix>.bam`) \
             and indexed (.bai or .csi).\n\
             Samples: {:?}",
            sample_names
        ),
    };
    if !bam_dir.is_dir() {
        bail!("--bam-dir {} is not a directory", bam_dir.display());
    }

    // Resolve one BAM per sample. Common suffixes such as `.sorted.bam`,
    // `.markdup.bam`, and `.cram` are accepted; see `bam_resolver`.
    let mut bam_paths: Vec<PathBuf> = Vec::with_capacity(sample_names.len());
    let mut missing: Vec<String> = Vec::new();
    for name in sample_names {
        match crate::bam_resolver::find_sample_bam(bam_dir, name) {
            Some(p) => bam_paths.push(p),
            None => missing.push(name.clone()),
        }
    }
    if !missing.is_empty() {
        bail!(
            "No BAM file found in {} for: {:?}\n\
             Looked for files named `<SAMPLE>.bam`, `<SAMPLE>.<suffix>.bam`, \
             or `<SAMPLE>.<suffix>.cram`.",
            bam_dir.display(),
            missing
        );
    }

    let runs = detect_ploidy_runs(input, sample_names)?;
    if runs.is_empty() {
        bail!("no ploidy runs detected in input VCF");
    }

    let temp_dir = tempfile::Builder::new()
        .prefix("vcf2fasta_var_ploidy_")
        .tempdir()
        .context("could not create temporary directory for variable-ploidy phasing")?;

    let mut phased_mini: Vec<PathBuf> = Vec::new();

    for (run_idx, run) in runs.iter().enumerate() {
        if run.end_pos == run.start_pos || run.ploidy < 2 {
            continue;
        }

        let mini_plain = temp_dir.path().join(format!("run_{}.vcf", run_idx));
        extract_run_vcf(input, run, &mini_plain)?;
        verify_single_sample(&mini_plain, &run.sample_name)?;

        let mini_gz = temp_dir.path().join(format!("run_{}.vcf.gz", run_idx));
        bgzip_file(&mini_plain, &mini_gz)?;
        tabix_index(&mini_gz)?;

        let cleaned = pre_filter_vcf(&mini_gz, reference, args, temp_dir.path())?;

        let phased_plain = temp_dir.path().join(format!("run_{}_phased.vcf", run_idx));
        let bam = bam_paths[run.sample_idx].clone();
        run_whatshap_polyphase(
            &cleaned,
            &[bam],
            reference,
            run.ploidy,
            &phased_plain,
            threads,
        )?;

        let phased_gz = temp_dir
            .path()
            .join(format!("run_{}_phased.vcf.gz", run_idx));
        bgzip_file(&phased_plain, &phased_gz)?;
        tabix_index(&phased_gz)?;
        phased_mini.push(phased_gz);
    }

    if phased_mini.is_empty() {
        bail!(
            "no ploidy run had enough positions to phase. \
             At least two adjacent variants with the same ploidy are required \
             for a run to produce a haplotype block."
        );
    }

    let merged_plain = temp_dir.path().join("merged.vcf");
    merge_runs_into_vcf(input, &phased_mini, sample_names, &merged_plain)?;

    let merged_gz = temp_dir.path().join("merged.vcf.gz");
    bgzip_file(&merged_plain, &merged_gz)?;
    tabix_index(&merged_gz)?;

    let persistent_dir = std::env::temp_dir().join(format!(
        "vcf2fasta_var_ploidy_out_{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&persistent_dir)?;
    let persistent = persistent_dir.join("variable_ploidy_phased.vcf.gz");
    std::fs::copy(&merged_gz, &persistent)?;
    tabix_index(&persistent)?;

    drop(temp_dir);
    Ok(persistent)
}

// ---------------------------------------------------------------------------
// Ploidy detection
// ---------------------------------------------------------------------------

fn ploidy_of_gt_string(gt: &str) -> usize {
    if gt.is_empty() {
        return 0;
    }
    let mut count = 1usize;
    for b in gt.bytes() {
        if b == b'/' || b == b'|' {
            count += 1;
        }
    }
    count
}

fn ploidy_of_sample_field(format: &str, sample: &str) -> Option<usize> {
    let fmt_parts: Vec<&str> = format.split(':').collect();
    let gt_pos = fmt_parts.iter().position(|f| *f == "GT")?;
    let smp_parts: Vec<&str> = sample.split(':').collect();
    let gt = smp_parts.get(gt_pos)?;
    Some(ploidy_of_gt_string(gt))
}

fn detect_ploidy_runs(input: &Path, sample_names: &[String]) -> Result<Vec<PloidyRun>> {
    let mut current: Vec<Option<PloidyRun>> = vec![None; sample_names.len()];
    let mut completed: Vec<PloidyRun> = Vec::new();
    let reader = open_maybe_gzip(input)?;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let raw = match parse_raw_vcf_record(&line) {
            Ok(r) => r,
            Err(_) => continue,
        };
        let contig = raw.chrom.to_string();
        let pos = raw.pos;

        for (idx, name) in sample_names.iter().enumerate() {
            if idx >= raw.samples.len() {
                continue;
            }
            let ploidy = match ploidy_of_sample_field(raw.format, raw.samples[idx]) {
                Some(p) if p > 0 => p,
                _ => continue,
            };
            match current[idx].as_mut() {
                Some(run) if run.contig == contig && run.ploidy == ploidy => {
                    run.end_pos = pos;
                }
                Some(run) => {
                    completed.push(run.clone());
                    current[idx] = Some(PloidyRun {
                        sample_idx: idx,
                        sample_name: name.clone(),
                        contig: contig.clone(),
                        start_pos: pos,
                        end_pos: pos,
                        ploidy,
                    });
                }
                None => {
                    current[idx] = Some(PloidyRun {
                        sample_idx: idx,
                        sample_name: name.clone(),
                        contig: contig.clone(),
                        start_pos: pos,
                        end_pos: pos,
                        ploidy,
                    });
                }
            }
        }
    }
    for opt in current.into_iter().flatten() {
        completed.push(opt);
    }
    Ok(completed)
}

// ---------------------------------------------------------------------------
// Extraction
// ---------------------------------------------------------------------------

fn extract_run_vcf(input: &Path, run: &PloidyRun, output: &Path) -> Result<()> {
    let reader = open_maybe_gzip(input)?;
    let mut out = File::create(output)
        .with_context(|| format!("could not create {}", output.display()))?;

    for line in reader.lines() {
        let line = line?;

        if line.starts_with("#CHROM") {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                bail!("malformed #CHROM line");
            }
            let mut new_fields: Vec<String> =
                fields[..9].iter().map(|s| s.to_string()).collect();
            new_fields.push(run.sample_name.clone());
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
        if fields.len() < 10 {
            continue;
        }
        let contig = fields[0];
        let pos: usize = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        if contig != run.contig {
            continue;
        }
        if pos < run.start_pos || pos > run.end_pos {
            continue;
        }
        let sample_col = 9 + run.sample_idx;
        if sample_col >= fields.len() {
            continue;
        }
        let mut new_fields: Vec<String> = fields[..9].iter().map(|s| s.to_string()).collect();
        new_fields.push(fields[sample_col].to_string());
        out.write_all(new_fields.join("\t").as_bytes())?;
        out.write_all(b"\n")?;
    }
    out.flush()?;
    Ok(())
}

fn verify_single_sample(vcf: &Path, expected_sample: &str) -> Result<()> {
    let reader = open_maybe_gzip(vcf)?;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() != 10 {
                bail!(
                    "mini-VCF {} has {} columns; expected exactly 10 (one sample)",
                    vcf.display(),
                    cols.len()
                );
            }
            if cols[9] != expected_sample {
                bail!(
                    "mini-VCF {} has sample '{}' but expected '{}'",
                    vcf.display(),
                    cols[9],
                    expected_sample
                );
            }
            return Ok(());
        }
    }
    bail!("mini-VCF {} has no #CHROM line", vcf.display())
}

// ---------------------------------------------------------------------------
// Streaming k-way merge
// ---------------------------------------------------------------------------

fn merge_runs_into_vcf(
    original: &Path,
    phased_mini: &[PathBuf],
    sample_names: &[String],
    output: &Path,
) -> Result<()> {
    let contig_order = read_contig_order(original)?;

    struct MiniStream {
        reader: Box<dyn BufRead>,
        sample: String,
        peeked: Option<(usize, usize, String)>,
    }

    let mut streams: Vec<MiniStream> = Vec::with_capacity(phased_mini.len());
    for mini in phased_mini {
        let mut reader = open_maybe_gzip(mini)?;
        let sample = read_mini_sample_name(&mut reader)?;
        let mut s = MiniStream {
            reader,
            sample,
            peeked: None,
        };
        s.peeked = peek_next_mini_record(&mut s.reader, &contig_order)?;
        streams.push(s);
    }

    let orig = open_maybe_gzip(original)?;
    let mut out = File::create(output)
        .with_context(|| format!("could not create {}", output.display()))?;

    let mut phased_here: HashMap<String, String> = HashMap::new();

    for line in orig.lines() {
        let line = line?;
        if line.starts_with('#') {
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 10 {
            continue;
        }
        let contig = cols[0];
        let pos: usize = match cols[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let cidx = match contig_order.get(contig) {
            Some(&c) => c,
            None => {
                out.write_all(line.as_bytes())?;
                out.write_all(b"\n")?;
                continue;
            }
        };

        phased_here.clear();
        for s in streams.iter_mut() {
            loop {
                match s.peeked.as_ref() {
                    Some((c, p, _)) if *c == cidx && *p == pos => {
                        let (_, _, gt) = s.peeked.take().unwrap();
                        phased_here.insert(s.sample.clone(), gt);
                        s.peeked = peek_next_mini_record(&mut s.reader, &contig_order)?;
                    }
                    Some((c, p, _)) if (*c, *p) < (cidx, pos) => {
                        s.peeked = peek_next_mini_record(&mut s.reader, &contig_order)?;
                    }
                    _ => break,
                }
            }
        }

        let mut new_fields: Vec<String> = cols[..9].iter().map(|s| s.to_string()).collect();
        new_fields[8] = "GT".to_string();
        for (idx, name) in sample_names.iter().enumerate() {
            let fallback = if 9 + idx < cols.len() {
                extract_gt_from_field(cols[8], cols[9 + idx])
            } else {
                ".".to_string()
            };
            let gt = phased_here.get(name).cloned().unwrap_or(fallback);
            new_fields.push(gt);
        }
        out.write_all(new_fields.join("\t").as_bytes())?;
        out.write_all(b"\n")?;
    }

    out.flush()?;
    Ok(())
}

fn read_mini_sample_name(reader: &mut Box<dyn BufRead>) -> Result<String> {
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            bail!("mini-VCF has no #CHROM line");
        }
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.trim_end().split('\t').collect();
            if cols.len() < 10 {
                bail!("mini-VCF #CHROM has <10 columns");
            }
            return Ok(cols[9].to_string());
        }
    }
}

fn peek_next_mini_record(
    reader: &mut Box<dyn BufRead>,
    contig_order: &HashMap<String, usize>,
) -> Result<Option<(usize, usize, String)>> {
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Ok(None);
        }
        let l = line.trim_end_matches('\n').trim_end_matches('\r');
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = l.split('\t').collect();
        if cols.len() < 10 {
            continue;
        }
        let contig = cols[0];
        let pos: usize = match cols[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let Some(&cidx) = contig_order.get(contig) else {
            continue;
        };
        let gt = extract_gt_from_field(cols[8], cols[9]);
        return Ok(Some((cidx, pos, gt)));
    }
}

fn read_contig_order(original: &Path) -> Result<HashMap<String, usize>> {
    if let Ok(reader) = rust_htslib::tbx::Reader::from_path(original) {
        let mut map = HashMap::new();
        for (i, name) in reader.seqnames().iter().enumerate() {
            map.insert(name.to_string(), i);
        }
        if !map.is_empty() {
            return Ok(map);
        }
    }
    let mut map: HashMap<String, usize> = HashMap::new();
    let reader = open_maybe_gzip(original)?;
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix("##contig=<ID=") {
            let id = rest.split(|c| c == ',' || c == '>').next().unwrap_or("");
            if !id.is_empty() {
                let n = map.len();
                map.entry(id.to_string()).or_insert(n);
            }
        } else if !line.starts_with('#') {
            break;
        }
    }
    Ok(map)
}

fn extract_gt_from_field(format: &str, sample: &str) -> String {
    let fmt_parts: Vec<&str> = format.split(':').collect();
    let gt_pos = match fmt_parts.iter().position(|f| *f == "GT") {
        Some(p) => p,
        None => return ".".to_string(),
    };
    let smp_parts: Vec<&str> = sample.split(':').collect();
    smp_parts.get(gt_pos).copied().unwrap_or(".").to_string()
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

    const MIXED_PLOIDY_VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=1,length=200>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
1\t1\tvar1\tG\tC\t60\tPASS\t.\tGT\t0/1/1\t1/0
1\t2\tvar2\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/1
1\t3\tvar3\tT\tA\t60\tPASS\t.\tGT\t0/1\t0/1
1\t4\tvar4\tC\tT\t60\tPASS\t.\tGT\t1/1/0\t1/0
1\t5\tvar5\tT\tC\t60\tPASS\t.\tGT\t0/1/0\t1/0
";

    #[test]
    fn ploidy_of_gt_string_works() {
        assert_eq!(ploidy_of_gt_string("0"), 1);
        assert_eq!(ploidy_of_gt_string("0/1"), 2);
        assert_eq!(ploidy_of_gt_string("0|1"), 2);
        assert_eq!(ploidy_of_gt_string("0/1/1"), 3);
        assert_eq!(ploidy_of_gt_string("1|1|0|0"), 4);
        assert_eq!(ploidy_of_gt_string("./."), 2);
        assert_eq!(ploidy_of_gt_string(""), 0);
    }

    #[test]
    fn ploidy_of_sample_field_works() {
        assert_eq!(ploidy_of_sample_field("GT", "0/1"), Some(2));
        assert_eq!(ploidy_of_sample_field("DP:GT", "10:0/1/1"), Some(3));
        assert_eq!(ploidy_of_sample_field("AD", "10,5"), None);
    }

    #[test]
    fn detects_variable_ploidy() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("v.vcf");
        write_file(&p, MIXED_PLOIDY_VCF.as_bytes());
        let samples = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];
        assert!(has_variable_ploidy(&p, &samples).unwrap());
    }

    #[test]
    fn constant_ploidy_is_not_flagged() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("c.vcf");
        let content = "\
##fileformat=VCFv4.2
##contig=<ID=1,length=200>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t1\tvar1\tG\tC\t60\tPASS\t.\tGT\t0/1
1\t2\tvar2\tA\tT\t60\tPASS\t.\tGT\t1/0
";
        write_file(&p, content.as_bytes());
        let samples = vec!["SAMPLE1".to_string()];
        assert!(!has_variable_ploidy(&p, &samples).unwrap());
    }

    #[test]
    fn detect_ploidy_runs_works() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("v.vcf");
        write_file(&p, MIXED_PLOIDY_VCF.as_bytes());
        let samples = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];
        let mut runs = detect_ploidy_runs(&p, &samples).unwrap();
        runs.sort_by(|a, b| {
            a.sample_idx
                .cmp(&b.sample_idx)
                .then(a.start_pos.cmp(&b.start_pos))
        });
        assert_eq!(runs.len(), 4);
        let s1: Vec<&PloidyRun> = runs.iter().filter(|r| r.sample_idx == 0).collect();
        assert_eq!(s1.len(), 3);
        assert_eq!((s1[0].start_pos, s1[0].end_pos, s1[0].ploidy), (1, 1, 3));
        assert_eq!((s1[1].start_pos, s1[1].end_pos, s1[1].ploidy), (2, 3, 2));
        assert_eq!((s1[2].start_pos, s1[2].end_pos, s1[2].ploidy), (4, 5, 3));
    }

    #[test]
    fn extract_run_vcf_produces_single_sample_vcf() {
        let dir = tempdir().unwrap();
        let p = dir.path().join("v.vcf");
        write_file(&p, MIXED_PLOIDY_VCF.as_bytes());
        let run = PloidyRun {
            sample_idx: 0,
            sample_name: "SAMPLE1".to_string(),
            contig: "1".to_string(),
            start_pos: 2,
            end_pos: 3,
            ploidy: 2,
        };
        let out = dir.path().join("mini.vcf");
        extract_run_vcf(&p, &run, &out).unwrap();
        let content = std::fs::read_to_string(&out).unwrap();
        let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 2);
        for l in &data_lines {
            let cols: Vec<&str> = l.split('\t').collect();
            assert_eq!(cols.len(), 10, "expected single-sample line: {}", l);
            assert_eq!(cols[9], "0/1");
        }
    }

    #[test]
    fn streaming_merge_preserves_uncovered_positions() {
        let dir = tempdir().unwrap();
        let original = dir.path().join("orig.vcf");
        write_file(&original, MIXED_PLOIDY_VCF.as_bytes());

        let mini = dir.path().join("mini.vcf");
        let mini_content = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t2\tvar2\tA\tT\t60\tPASS\t.\tGT\t0|1
1\t3\tvar3\tT\tA\t60\tPASS\t.\tGT\t1|0
";
        write_file(&mini, mini_content.as_bytes());

        let samples = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];
        let out = dir.path().join("merged.vcf");
        merge_runs_into_vcf(&original, &[mini], &samples, &out).unwrap();

        let content = std::fs::read_to_string(&out).unwrap();
        let data: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data.len(), 5);

        let p2: Vec<&str> = data[1].split('\t').collect();
        assert_eq!(p2[9], "0|1");
        assert_eq!(p2[10], "0/1");

        let p1: Vec<&str> = data[0].split('\t').collect();
        assert_eq!(p1[9], "0/1/1");
    }

    #[test]
    fn streaming_merge_handles_two_samples_and_multiple_runs() {
        let dir = tempdir().unwrap();
        let original = dir.path().join("orig.vcf");
        write_file(&original, MIXED_PLOIDY_VCF.as_bytes());

        let s1a = dir.path().join("s1a.vcf");
        let s1b = dir.path().join("s1b.vcf");
        let s2 = dir.path().join("s2.vcf");
        write_file(&s1a, b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n1\t1\t.\tG\tC\t60\tPASS\t.\tGT\t0|1|1\n");
        write_file(&s1b, b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n1\t2\t.\tA\tT\t60\tPASS\t.\tGT\t0|1\n1\t3\t.\tT\tA\t60\tPASS\t.\tGT\t1|0\n1\t4\t.\tC\tT\t60\tPASS\t.\tGT\t1|1|0\n1\t5\t.\tT\tC\t60\tPASS\t.\tGT\t0|1|0\n");
        write_file(&s2, b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE2\n1\t2\t.\tA\tT\t60\tPASS\t.\tGT\t0|1\n1\t4\t.\tC\tT\t60\tPASS\t.\tGT\t1|0\n");

        let samples = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];
        let out = dir.path().join("merged.vcf");
        merge_runs_into_vcf(&original, &[s1a, s1b, s2], &samples, &out).unwrap();

        let content = std::fs::read_to_string(&out).unwrap();
        let data: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data.len(), 5);

        let p1: Vec<&str> = data[0].split('\t').collect();
        assert_eq!(p1[9], "0|1|1");
        assert_eq!(p1[10], "1/0");

        let p2: Vec<&str> = data[1].split('\t').collect();
        assert_eq!(p2[9], "0|1");
        assert_eq!(p2[10], "0|1");

        let p3: Vec<&str> = data[2].split('\t').collect();
        assert_eq!(p3[9], "1|0");
        assert_eq!(p3[10], "0/1");
    }

    #[test]
    fn streaming_merge_handles_gzipped_inputs() {
        let dir = tempdir().unwrap();
        let original = dir.path().join("orig.vcf.gz");
        write_gzip(&original, MIXED_PLOIDY_VCF.as_bytes());

        let mini = dir.path().join("mini.vcf.gz");
        let mini_content = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n1\t2\t.\tA\tT\t60\tPASS\t.\tGT\t0|1\n";
        write_gzip(&mini, mini_content);

        let samples = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];
        let out = dir.path().join("merged.vcf");
        merge_runs_into_vcf(&original, &[mini], &samples, &out).unwrap();
        let content = std::fs::read_to_string(&out).unwrap();
        assert!(content.contains("0|1"));
    }
}