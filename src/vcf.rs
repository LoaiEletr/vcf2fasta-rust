//! VCF parsing and processing logic.

use crate::{
    chunk::{ChunkResult, GenomicChunk},
    cli::Args,
    genotype::{is_effectively_phased, AlleleCall, DecodedRecord, normalize_genotype},
    report::ContigReport,
    util::{fetch_half_open, is_symbolic_or_breakend, is_valid_nucleotide_sequence},
};
use anyhow::{anyhow, bail, Context, Result};
use rust_htslib::{
    bcf::{self, Read as BcfRead},
    faidx,
    tbx::{self, Read as TbxRead},
};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};
use flate2::read::GzDecoder;

// ============================================================================
// HEADER VALIDATION HELPERS
// ============================================================================

fn read_chrom_line(path: &Path) -> Result<String> {
    let file = File::open(path)
        .with_context(|| format!("cannot open {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        let decoder = GzDecoder::new(file);
        Box::new(BufReader::new(decoder))
    } else {
        Box::new(BufReader::new(file))
    };
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("#CHROM") {
            return Ok(line);
        }
    }
    bail!("no '#CHROM' line found in VCF header");
}

fn validate_vcf_header(path: &Path) -> Result<()> {
    let chrom_line = read_chrom_line(path)?;
    let parts: Vec<&str> = chrom_line.split('\t').collect();
    if parts.len() < 10 {
        bail!(
            "VCF header has only {} columns; expected at least 10",
            parts.len()
        );
    }
    Ok(())
}

pub fn read_header_metadata(path: &Path) -> Result<(Vec<String>, Vec<String>)> {
    validate_vcf_header(path)?;
    let reader = bcf::Reader::from_path(path)
        .with_context(|| format!("could not open VCF/BCF {}", path.display()))?;
    let header = reader.header();
    let samples = header
        .samples()
        .into_iter()
        .map(|s| {
            std::str::from_utf8(s)
                .map(str::to_owned)
                .context("sample name is not valid UTF-8")
        })
        .collect::<Result<Vec<_>>>()?;
    let mut contigs = Vec::with_capacity(header.contig_count() as usize);
    for rid in 0..header.contig_count() {
        let name = header
            .rid2name(rid)
            .with_context(|| format!("could not resolve VCF contig RID {rid}"))?;
        contigs.push(
            std::str::from_utf8(name)
                .context("contig name is not valid UTF-8")?
                .to_owned(),
        );
    }
    Ok((samples, contigs))
}

// ============================================================================
// RAW RECORD STRUCT
// ============================================================================

#[derive(Debug)]
pub(crate) struct RawVcfRecord<'a> {
    pub(crate) chrom: &'a str,
    pub(crate) pos: usize,
    pub(crate) ref_allele: &'a str,
    pub(crate) alt_field: &'a str,
    pub(crate) format: &'a str,
    pub(crate) samples: Vec<&'a str>,
}

pub(crate) fn parse_raw_vcf_record(line: &str) -> Result<RawVcfRecord<'_>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        bail!("too few columns: expected at least 9, got {}", fields.len());
    }
    let chrom = fields[0];
    if chrom.is_empty() {
        bail!("CHROM field is empty");
    }
    let pos_1_based = fields[1]
        .parse::<usize>()
        .with_context(|| format!("invalid POS '{}'", fields[1]))?;
    if pos_1_based == 0 {
        bail!("POS must be >= 1");
    }
    let ref_allele = fields[3];
    if ref_allele.is_empty() || ref_allele == "." {
        bail!("REF is missing");
    }
    let alt_field = fields[4];
    let format = fields[8];
    let mut samples = if fields.len() > 9 {
        fields[9..].to_vec()
    } else {
        Vec::new()
    };
    while let Some(last) = samples.last() {
        if last.is_empty() {
            samples.pop();
        } else {
            break;
        }
    }
    Ok(RawVcfRecord {
        chrom,
        pos: pos_1_based,
        ref_allele,
        alt_field,
        format,
        samples,
    })
}

// ============================================================================
// GENOTYPE PARSING
// ============================================================================
//
// NOTE: `parse_gt` and `extract_gt` are `pub(crate)` because the discovery
// pass in `crate::workunit` reuses them. They remain internal to the crate.

pub(crate) fn parse_gt(gt: &str) -> Result<Vec<bcf::record::GenotypeAllele>> {
    use bcf::record::GenotypeAllele as GA;
    if gt.is_empty() {
        bail!("GT field is empty");
    }
    let mut result = Vec::new();
    let mut start = 0;
    let bytes = gt.as_bytes();
    let mut current_is_phased = false;
    for i in 0..=bytes.len() {
        let at_end = i == bytes.len();
        let is_separator = !at_end && (bytes[i] == b'/' || bytes[i] == b'|');
        if !at_end && !is_separator {
            continue;
        }
        let allele_text = &gt[start..i];
        let allele = if allele_text == "." {
            if current_is_phased {
                GA::PhasedMissing
            } else {
                GA::UnphasedMissing
            }
        } else {
            let idx = allele_text
                .parse::<i32>()
                .with_context(|| format!("invalid allele index '{}'", allele_text))?;
            if idx < 0 {
                if current_is_phased {
                    GA::PhasedMissing
                } else {
                    GA::UnphasedMissing
                }
            } else if current_is_phased {
                GA::Phased(idx)
            } else {
                GA::Unphased(idx)
            }
        };
        result.push(allele);
        if at_end {
            break;
        }
        current_is_phased = bytes[i] == b'|';
        start = i + 1;
    }
    if result.is_empty() {
        bail!("GT contains no alleles");
    }
    Ok(result)
}

pub(crate) fn extract_gt<'a>(format: &str, sample: &'a str) -> Result<&'a str> {
    let format_fields: Vec<&str> = format.split(':').collect();
    let gt_pos = format_fields
        .iter()
        .position(|f| *f == "GT")
        .ok_or_else(|| anyhow!("FORMAT does not contain GT"))?;
    let sample_fields: Vec<&str> = sample.split(':').collect();
    sample_fields
        .get(gt_pos)
        .copied()
        .ok_or_else(|| {
            anyhow!(
                "sample has {} fields, but GT is at position {}",
                sample_fields.len(),
                gt_pos + 1
            )
        })
}

// ============================================================================
// MAX PLOIDY TRACKING
// ============================================================================

pub(crate) fn update_max_ploidies_from_line(
    line: &str,
    sample_names: &[String],
    max_ploidies: &mut [usize],
) {
    if line.starts_with('#') {
        return;
    }
    let raw = match parse_raw_vcf_record(line) {
        Ok(r) => r,
        Err(_) => return,
    };
    let sample_count = sample_names.len();
    let original_sample_count = raw.samples.len();
    let mut sample_fields = raw.samples.to_vec();
    if sample_fields.len() < sample_count {
        sample_fields.resize(sample_count, "");
    } else if sample_fields.len() > sample_count {
        sample_fields.truncate(sample_count);
    }
    for (sample_idx, sample_field) in sample_fields.iter().enumerate() {
        if sample_idx >= original_sample_count {
            continue;
        }
        let gt_text = match extract_gt(raw.format, sample_field) {
            Ok(g) => g,
            Err(_) => continue,
        };
        let raw_gt = match parse_gt(gt_text) {
            Ok(g) => g,
            Err(_) => continue,
        };
        let normalized = match normalize_genotype(&raw_gt, 2) {
            Ok(n) => n,
            Err(_) => continue,
        };
        let ploidy = normalized.len();
        if ploidy > max_ploidies[sample_idx] {
            max_ploidies[sample_idx] = ploidy;
        }
    }
}

// ============================================================================
// LINE REWRITING (used by phasing pre-filter)
// ============================================================================

pub(crate) fn clean_vcf_line(
    line: &str,
    sample_names: &[String],
    max_ploidies: &[usize],
) -> Option<String> {
    let raw = parse_raw_vcf_record(line).ok()?;
    let mut allele_count = 1usize;
    if raw.alt_field != "." {
        allele_count += raw.alt_field.split(',').filter(|s| !s.is_empty()).count();
    }
    let original_fields: Vec<&str> = line.split('\t').collect();
    if original_fields.len() < 9 {
        return None;
    }
    let mut out_fields: Vec<String> = original_fields[..8]
        .iter()
        .map(|s| s.to_string())
        .collect();
    out_fields.push("GT".to_string());
    for (sample_idx, _) in sample_names.iter().enumerate() {
        out_fields.push(clean_gt_for_sample(&raw, sample_idx, max_ploidies, allele_count));
    }
    Some(out_fields.join("\t"))
}

fn clean_gt_for_sample(
    raw: &RawVcfRecord<'_>,
    sample_idx: usize,
    max_ploidies: &[usize],
    allele_count: usize,
) -> String {
    let fallback = max_ploidies.get(sample_idx).copied().unwrap_or(2).max(1);
    let all_ref = || vec!["0"; fallback].join("|");
    if sample_idx >= raw.samples.len() {
        return all_ref();
    }
    let gt_text = match extract_gt(raw.format, raw.samples[sample_idx]) {
        Ok(g) => g,
        Err(_) => return all_ref(),
    };
    let raw_gt = match parse_gt(gt_text) {
        Ok(g) => g,
        Err(_) => return all_ref(),
    };
    let normalized = match normalize_genotype(&raw_gt, 2) {
        Ok(n) => n,
        Err(_) => return all_ref(),
    };
    let mut parts = Vec::with_capacity(normalized.len());
    for allele_opt in normalized {
        match allele_opt {
            Some(idx) if idx < allele_count => parts.push(idx.to_string()),
            _ => parts.push(".".to_string()),
        }
    }
    parts.join("|")
}

// ============================================================================
// DECODE & VALIDATE ONE RECORD
// ============================================================================

pub(crate) fn decode_and_validate_record(
    raw: &RawVcfRecord<'_>,
    contig: &str,
    contig_len: usize,
    sample_names: &[String],
    max_ploidies: &[usize],
    reference: &faidx::Reader,
    args: &Args,
    report: &mut ContigReport,
) -> Result<DecodedRecord> {
    if raw.chrom != contig {
        bail!("record belongs to '{}' not '{}'", raw.chrom, contig);
    }
    let start = raw.pos.checked_sub(1).ok_or_else(|| anyhow!("POS must be >= 1"))?;
    if start >= contig_len {
        bail!("POS {} is outside contig '{}' (length {})", raw.pos, contig, contig_len);
    }
    let mut alleles: Vec<Vec<u8>> = Vec::new();
    alleles.push(raw.ref_allele.as_bytes().to_vec());
    if raw.alt_field != "." {
        for alt in raw.alt_field.split(',') {
            if alt.is_empty() {
                bail!("ALT contains empty allele");
            }
            alleles.push(alt.as_bytes().to_vec());
        }
    }
    if alleles.is_empty() || alleles[0].is_empty() {
        bail!("record has no valid REF allele");
    }
    let mut allele_valid = vec![true; alleles.len()];
    for (idx, allele) in alleles.iter().enumerate() {
        if idx == 0 {
            if !is_valid_nucleotide_sequence(allele) || is_symbolic_or_breakend(allele) {
                bail!("REF allele is invalid (non-nucleotide or symbolic); skipping variant");
            }
            continue;
        }
        if is_symbolic_or_breakend(allele) {
            allele_valid[idx] = false;
        } else if !is_valid_nucleotide_sequence(allele) {
            allele_valid[idx] = false;
        }
    }
    let end = start
        .checked_add(alleles[0].len())
        .ok_or_else(|| anyhow!("REF interval overflow"))?;
    if end > contig_len {
        bail!(
            "REF interval [{}, {}) exceeds contig length {}",
            start, end, contig_len
        );
    }
    if !args.no_validate_ref {
        let observed_ref = fetch_half_open(reference, contig, start, end)?;
        if !observed_ref.eq_ignore_ascii_case(&alleles[0]) {
            bail!(
                "REF mismatch: VCF '{}' vs FASTA '{}'",
                String::from_utf8_lossy(&alleles[0]),
                String::from_utf8_lossy(&observed_ref)
            );
        }
    }
    let sample_count = sample_names.len();
    let original_sample_count = raw.samples.len();
    let mut sample_fields = raw.samples.to_vec();
    if sample_fields.len() < sample_count {
        sample_fields.resize(sample_count, "");
    } else if sample_fields.len() > sample_count {
        sample_fields.truncate(sample_count);
    }
    let mut decoded_genotypes = Vec::with_capacity(sample_count);
    let mut ploidies = Vec::with_capacity(sample_count);
    for (sample_idx, sample_field) in sample_fields.iter().enumerate() {
        let sample_name = &sample_names[sample_idx];
        let fallback_ploidy = max_ploidies.get(sample_idx).copied().unwrap_or(2).max(1);
        if sample_idx >= original_sample_count {
            let calls = vec![AlleleCall::Reference; fallback_ploidy];
            decoded_genotypes.push(calls);
            ploidies.push(fallback_ploidy);
            report.warn(format!(
                "[{}:{}] sample '{}' missing from record, using {} reference call(s)",
                contig, raw.pos, sample_name, fallback_ploidy
            ));
            continue;
        }
        let gt_text = match extract_gt(raw.format, sample_field) {
            Ok(g) => g,
            Err(e) => {
                let calls = vec![AlleleCall::Reference; fallback_ploidy];
                decoded_genotypes.push(calls);
                ploidies.push(fallback_ploidy);
                report.warn(format!(
                    "[{}:{}] sample '{}' cannot extract GT ({}), using {} reference call(s)",
                    contig, raw.pos, sample_name, e, fallback_ploidy
                ));
                continue;
            }
        };
        let raw_gt = match parse_gt(gt_text) {
            Ok(g) => g,
            Err(e) => {
                let calls = vec![AlleleCall::Reference; fallback_ploidy];
                decoded_genotypes.push(calls);
                ploidies.push(fallback_ploidy);
                report.warn(format!(
                    "[{}:{}] sample '{}' malformed GT '{}' ({}), using {} reference call(s)",
                    contig, raw.pos, sample_name, gt_text, e, fallback_ploidy
                ));
                continue;
            }
        };
        let normalized = match normalize_genotype(&raw_gt, 2) {
            Ok(n) => n,
            Err(e) => {
                let calls = vec![AlleleCall::Reference; fallback_ploidy];
                decoded_genotypes.push(calls);
                ploidies.push(fallback_ploidy);
                report.warn(format!(
                    "[{}:{}] sample '{}' normalisation failed ({}), using {} reference call(s)",
                    contig, raw.pos, sample_name, e, fallback_ploidy
                ));
                continue;
            }
        };
        let ploidy = normalized.len();
        let mut calls = Vec::with_capacity(ploidy);
        let placeholder = args.no_call_string.as_deref().unwrap_or("N");
        let allele_strings: Vec<&str> = gt_text.split(|c| c == '/' || c == '|').collect();
        for (allele_pos, allele_index) in normalized.into_iter().enumerate() {
            let call = match allele_index {
                Some(idx) => {
                    if idx >= alleles.len() || !allele_valid[idx] {
                        report.warn(format!(
                            "[{}:{}] sample '{}' has invalid allele {} in genotype '{}', using placeholder '{}'",
                            contig, raw.pos, sample_name, idx, gt_text, placeholder
                        ));
                        AlleleCall::Missing
                    } else {
                        AlleleCall::Index(idx)
                    }
                }
                None => {
                    let allele_str = allele_strings.get(allele_pos).unwrap_or(&"");
                    if allele_str.starts_with('-') {
                        report.warn(format!(
                            "[{}:{}] sample '{}' has negative allele '{}', using placeholder '{}'",
                            contig, raw.pos, sample_name, allele_str, placeholder
                        ));
                        AlleleCall::Missing
                    } else {
                        report.warn(format!(
                            "[{}:{}] sample '{}' missing allele in genotype '{}', using placeholder '{}'",
                            contig, raw.pos, sample_name, gt_text, placeholder
                        ));
                        AlleleCall::Missing
                    }
                }
            };
            calls.push(call);
        }
        if sample_count == 1
            && calls.len() == 2
            && !is_effectively_phased(&raw_gt)
        {
            if let (AlleleCall::Index(a), AlleleCall::Index(b)) = (&calls[0], &calls[1]) {
                if a > b {
                    calls.swap(0, 1);
                }
            }
        }
        decoded_genotypes.push(calls);
        ploidies.push(ploidy);
    }
    Ok(DecodedRecord {
        start,
        end,
        alleles,
        genotypes: decoded_genotypes,
        ploidies,
    })
}

// ============================================================================
// REGION PROCESSING (chunk and contig)
// ============================================================================

/// Processes a single genomic chunk and returns its partial output.
///
/// `pad` is the look-back distance (in bp) used to detect a variant from a
/// previous chunk whose REF allele extends into this chunk. It must be at
/// least as large as the longest REF allele in the VCF.
pub fn process_chunk(
    chunk: &GenomicChunk,
    sample_names: &[String],
    max_ploidies: &[usize],
    pad: u64,
    args: &Args,
) -> Result<ChunkResult> {
    let mut report = ContigReport::new(chunk.contig.clone());

    let reference = faidx::Reader::from_path(&args.reference)
        .with_context(|| format!("could not open reference {}", args.reference.display()))?;
    let contig_len_u64 = reference.fetch_seq_len(&chunk.contig);
    let contig_len = usize::try_from(contig_len_u64)
        .map_err(|_| anyhow!("reference contig {} too long for usize", chunk.contig))?;

    let chunk_start = chunk.start as usize;
    let chunk_end = chunk.end as usize;
    let query_start = (chunk.start.saturating_sub(pad)) as usize;

    let mut per_sample: Vec<Vec<Vec<u8>>> = (0..sample_names.len())
        .map(|s| (0..max_ploidies[s]).map(|_| Vec::new()).collect())
        .collect();

    let mut vcf = tbx::Reader::from_path(&args.input).with_context(|| {
        format!("could not open {} as tabix-indexed text VCF", args.input.display())
    })?;
    let tid = vcf
        .tid(&chunk.contig)
        .with_context(|| format!("contig {} not found in tabix index", chunk.contig))?;
    let query_end_u64 = u64::try_from(chunk_end)
        .map_err(|_| anyhow!("chunk end does not fit u64"))?;
    let query_start_u64 = u64::try_from(query_start)
        .map_err(|_| anyhow!("query start does not fit u64"))?;
    vcf.fetch(tid, query_start_u64, query_end_u64)
        .with_context(|| "could not fetch contig region".to_string())?;

    let mut ref_cursor = chunk_start;
    let mut buffer = Vec::<u8>::with_capacity(1024);
    let placeholder = args.no_call_string.as_deref().unwrap_or("N");

    while let Ok(has_record) = vcf.read(&mut buffer) {
        if !has_record {
            break;
        }
        let line = match std::str::from_utf8(&buffer) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let line = line.trim_end_matches('\r');
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[0] != chunk.contig {
            continue;
        }

        let raw = match parse_raw_vcf_record(line) {
            Ok(r) => r,
            Err(e) => {
                report.warn(format!("[{}:{}] malformed record: {}", chunk.contig, fields[1], e));
                continue;
            }
        };

        let var_start = match raw.pos.checked_sub(1) {
            Some(p) => p,
            None => continue,
        };
        let var_ref_len = raw.ref_allele.len();
        let var_end = var_start + var_ref_len;

        if var_start < chunk_start {
            if var_end > chunk_start && var_end > ref_cursor {
                ref_cursor = var_end;
            }
            continue;
        }

        if var_start >= chunk_end {
            continue;
        }

        report.seen += 1;

        let decoded = match decode_and_validate_record(
            &raw,
            &chunk.contig,
            contig_len,
            sample_names,
            max_ploidies,
            &reference,
            args,
            &mut report,
        ) {
            Ok(d) => d,
            Err(e) => {
                report.warn(format!("[{}:{}] {}", chunk.contig, raw.pos, e));
                continue;
            }
        };

        if decoded.start < ref_cursor {
            report.warn(format!(
                "[{}:{}] overlap/out-of-order (start={} < prev_end={})",
                chunk.contig, raw.pos, decoded.start, ref_cursor
            ));
            continue;
        }

        if decoded.start > ref_cursor {
            let between = fetch_half_open(&reference, &chunk.contig, ref_cursor, decoded.start)?;
            for s in 0..sample_names.len() {
                for h in 0..max_ploidies[s] {
                    per_sample[s][h].extend_from_slice(&between);
                }
            }
        }

        for (s_idx, sample_gt) in decoded.genotypes.iter().enumerate() {
            let mut padded = sample_gt.clone();
            while padded.len() < max_ploidies[s_idx] {
                padded.push(AlleleCall::Reference);
            }
            if padded.len() > max_ploidies[s_idx] {
                padded.truncate(max_ploidies[s_idx]);
            }
            for h_idx in 0..max_ploidies[s_idx] {
                match padded[h_idx] {
                    AlleleCall::Index(idx) => {
                        per_sample[s_idx][h_idx].extend_from_slice(&decoded.alleles[idx]);
                    }
                    AlleleCall::Missing => {
                        per_sample[s_idx][h_idx].extend_from_slice(placeholder.as_bytes());
                    }
                    AlleleCall::Reference => {
                        let base = fetch_half_open(
                            &reference,
                            &chunk.contig,
                            decoded.start,
                            decoded.start + 1,
                        )?;
                        per_sample[s_idx][h_idx].extend_from_slice(&base);
                    }
                }
            }
        }

        ref_cursor = decoded.end;
        report.applied += 1;
    }

    if ref_cursor < chunk_end {
        let tail = fetch_half_open(&reference, &chunk.contig, ref_cursor, chunk_end)?;
        for s in 0..sample_names.len() {
            for h in 0..max_ploidies[s] {
                per_sample[s][h].extend_from_slice(&tail);
            }
        }
    }

    let output_files: usize = max_ploidies.iter().sum();

    Ok(ChunkResult {
        contig: chunk.contig.clone(),
        start: chunk.start,
        end: chunk.end,
        per_sample,
        seen: report.seen,
        applied: report.applied,
        skipped: report.skipped,
        output_files,
        warning_count: report.warning_count,
        warnings: report.warnings,
    })
}

/// Processes an entire contig as one chunk. Retained for backwards
/// compatibility with tests and non-chunked callers.
pub fn process_contig(
    contig: &str,
    sample_names: &[String],
    args: &Args,
) -> Result<ChunkResult> {
    // Compute max ploidies by scanning the VCF.
    let mut max_ploidies = vec![0usize; sample_names.len()];
    {
        let reader = crate::phasing::open_maybe_gzip(&args.input)?;
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            update_max_ploidies_from_line(&line, sample_names, &mut max_ploidies);
        }
    }
    for p in &mut max_ploidies {
        if *p == 0 {
            *p = 2;
        }
    }

    let pad = args.chunk_pad.unwrap_or(10_000);

    let reference = faidx::Reader::from_path(&args.reference)
        .with_context(|| format!("could not open reference {}", args.reference.display()))?;
    let contig_len = reference.fetch_seq_len(contig);
    let chunk = GenomicChunk {
        contig: contig.to_string(),
        start: 0,
        end: contig_len,
    };
    process_chunk(&chunk, sample_names, &max_ploidies, pad, args)
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
pub mod test_helpers {
    use super::*;
    use crate::cli::Args;
    use std::fs::File;
    use std::io::Write;
    use std::path::PathBuf;
    use tempfile::tempdir;
    use rust_htslib::faidx;

    pub fn minimal_vcf_content() -> (String, Vec<String>) {
        let header = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=20>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
"#;
        let records = vec![
            "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|0".to_string(),
            "chr1\t10\t.\tC\tT\t.\t.\t.\tGT\t0|0\t1|1".to_string(),
        ];
        (header.to_string(), records)
    }

    pub fn create_test_reference() -> (tempfile::TempDir, PathBuf, usize) {
        let dir = tempdir().unwrap();
        let ref_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        let fasta_content = ">chr1\nACGTACGTACGTACGTACGT\n";
        File::create(&ref_path).unwrap().write_all(fasta_content.as_bytes()).unwrap();
        let fai_content = "chr1\t20\t6\t20\t21\n";
        File::create(&fai_path).unwrap().write_all(fai_content.as_bytes()).unwrap();
        (dir, ref_path, 20)
    }

    pub fn test_args() -> Args {
        Args {
            input: PathBuf::from("dummy.vcf"),
            reference: PathBuf::from("dummy.fa"),
            prefix: "test_".to_string(),
            no_call_string: Some("N".to_string()),
            threads: 1,
            line_width: 80,
            no_validate_ref: false,
            quiet: false,
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

    fn open_reference(ref_path: &PathBuf) -> faidx::Reader {
        faidx::Reader::from_path(ref_path).unwrap()
    }

    fn decode_line(
        line: &str,
        contig: &str,
        contig_len: usize,
        sample_names: &[String],
        reference: &faidx::Reader,
        args: &Args,
    ) -> Result<DecodedRecord> {
        let raw = parse_raw_vcf_record(line)?;
        let mut max_ploidies = vec![0usize; sample_names.len()];
        update_max_ploidies_from_line(line, sample_names, &mut max_ploidies);
        for p in &mut max_ploidies { if *p == 0 { *p = 2; } }
        let mut report = ContigReport::new(contig.to_owned());
        decode_and_validate_record(
            &raw, contig, contig_len, sample_names, &max_ploidies,
            reference, args, &mut report,
        )
    }

    #[test]
    fn read_header_metadata_ok() {
        let dir = tempdir().unwrap();
        let vcf_path = dir.path().join("test.vcf");
        let (header, records) = minimal_vcf_content();
        let mut content = header;
        for rec in records {
            content.push_str(&rec);
            content.push('\n');
        }
        File::create(&vcf_path).unwrap().write_all(content.as_bytes()).unwrap();
        let (samples, contigs) = read_header_metadata(&vcf_path).unwrap();
        assert_eq!(samples, vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()]);
        assert_eq!(contigs, vec!["chr1".to_string()]);
    }

    #[test]
    fn decode_valid_record() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|0";
        let decoded = decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(decoded.start, 4);
        assert_eq!(decoded.end, 5);
        assert_eq!(decoded.genotypes[0], vec![AlleleCall::Index(0), AlleleCall::Index(1)]);
        assert_eq!(decoded.genotypes[1], vec![AlleleCall::Index(1), AlleleCall::Index(0)]);
    }

    #[test]
    fn decode_record_with_malformed_gt() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["SAMPLE1".to_string()];
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\tABC";
        let decoded = decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(decoded.ploidies, vec![2]);
        assert_eq!(decoded.genotypes[0], vec![AlleleCall::Reference, AlleleCall::Reference]);
    }

    #[test]
    fn decode_record_with_unphased_genotype() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["SAMPLE1".to_string()];
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1";
        let decoded = decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(decoded.ploidies, vec![2]);
        assert_eq!(decoded.genotypes[0], vec![AlleleCall::Index(0), AlleleCall::Index(1)]);
    }

    #[test]
    fn decode_record_with_negative_allele() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["SAMPLE1".to_string()];
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t-1|1";
        let decoded = decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(decoded.genotypes[0], vec![AlleleCall::Missing, AlleleCall::Index(1)]);
    }

    #[test]
    fn decode_record_with_ref_interval_overflow() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["SAMPLE1".to_string()];
        let line = "chr1\t20\t.\tAA\tG\t.\t.\t.\tGT\t0|1";
        let result = decode_line(line, "chr1", contig_len, &sample_names, &reference, &args);
        assert!(result.is_err());
    }

    // ========================================================================
    // Policy tests
    // ========================================================================

    // ---- Invalid POS (negative, zero, non-numeric) ----

    #[test]
    fn invalid_pos_negative() {
        let result = parse_raw_vcf_record("chr1\t-1\t.\tA\tG\t.\t.\t.\tGT\t0|1");
        assert!(result.is_err(), "POS=-1 must be rejected");
    }

    #[test]
    fn invalid_pos_zero() {
        let result = parse_raw_vcf_record("chr1\t0\t.\tA\tG\t.\t.\t.\tGT\t0|1");
        assert!(result.is_err(), "POS=0 must be rejected");
    }

    #[test]
    fn invalid_pos_nan() {
        let result = parse_raw_vcf_record("chr1\tABC\t.\tA\tG\t.\t.\t.\tGT\t0|1");
        assert!(result.is_err(), "POS=ABC must be rejected");
    }

    // ---- POS > FASTA contig length ----

    #[test]
    fn pos_beyond_contig() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // Contig is 20 bp; POS=21 is out of range.
        let line = "chr1\t21\t.\tA\tG\t.\t.\t.\tGT\t0|1";
        let result = decode_line(line, "chr1", contig_len, &sample_names, &reference, &args);
        assert!(result.is_err(), "POS beyond contig length must be rejected");
    }

    // ---- Invalid allele index (0|3) ----

    #[test]
    fn invalid_allele_index_becomes_n() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // REF=A, ALT=G. Only 2 alleles (indices 0 and 1). GT=0|3 selects
        // index 3, which does not exist.
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0|3";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Missing]
        );
    }

    // ---- Negative allele index (-1|0) ----

    #[test]
    fn negative_allele_index_becomes_n() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t-1|0";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Missing, AlleleCall::Index(0)]
        );
    }

    // ---- Too few mandatory VCF columns ----

    #[test]
    fn too_few_columns() {
        // 8 columns; needs at least 9.
        let result = parse_raw_vcf_record("chr1\t5\t.\tA\tG\t.\t.\t.");
        assert!(result.is_err(), "record with <9 columns must be rejected");
    }

    // ---- Invalid REF nucleotide ----

    #[test]
    fn invalid_ref_nucleotide() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // REF='X' is not a valid nucleotide.
        let line = "chr1\t5\t.\tX\tG\t.\t.\t.\tGT\t0|1";
        let result = decode_line(line, "chr1", contig_len, &sample_names, &reference, &args);
        assert!(result.is_err(), "non-nucleotide REF must be rejected");
    }

    // ---- REF does not match FASTA ----

    #[test]
    fn ref_mismatch_fasta() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // FASTA has 'A' at pos 5; VCF claims REF='C'.
        let line = "chr1\t5\t.\tC\tG\t.\t.\t.\tGT\t0|1";
        let result = decode_line(line, "chr1", contig_len, &sample_names, &reference, &args);
        assert!(result.is_err(), "REF mismatch must be rejected when validation is on");
    }

    // ---- Invalid selected ALT nucleotide ----

    #[test]
    fn invalid_alt_selected_becomes_n() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // ALT='1' is not a valid nucleotide. Genotype 0|1 selects it.
        let line = "chr1\t5\t.\tA\t1\t.\t.\t.\tGT\t0|1";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Missing]
        );
    }

    // ---- Invalid unused ALT nucleotide is ignored ----

    #[test]
    fn invalid_unused_alt_ignored() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // ALT field contains '1' (invalid) followed by 'G' (valid).
        // Genotype 0|2 selects the valid ALT at index 2.
        // (The invalid ALT at index 1 does not affect the result.)
        let line = "chr1\t5\t.\tA\t1,G\t.\t.\t.\tGT\t0|2";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Index(2)]
        );
    }

    // ---- Record has fewer samples than the header declares ----

    #[test]
    fn record_has_fewer_samples() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string(), "S2".to_string()];
        // Header declares 2 samples; the record has only 1 sample column.
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0|1";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        // S1 gets the called genotype.
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Index(1)]
        );
        // S2 falls back to Reference calls.
        assert_eq!(
            decoded.genotypes[1],
            vec![AlleleCall::Reference, AlleleCall::Reference]
        );
    }

    // ---- Record has more samples than the header declares ----

    #[test]
    fn record_has_more_samples() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // Header declares 1 sample; the record has 2 sample columns.
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|0";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        // Only S1's genotype is decoded; the second column is ignored.
        assert_eq!(decoded.genotypes.len(), 1);
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Index(1)]
        );
    }

    // ---- Invalid genotype separator (0&1) ----

    #[test]
    fn invalid_separator() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // '&' is not a valid GT separator. Parser treats "0&1" as one
        // allele token, fails to parse, and falls back to Reference calls.
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0&1";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(decoded.ploidies, vec![2]);
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Reference, AlleleCall::Reference]
        );
    }

    // ---- Missing genotype allele (.|.) ----

    #[test]
    fn double_missing_becomes_n_n() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t.|.";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Missing, AlleleCall::Missing]
        );
    }

    // ---- Partially missing (.|0) ----

    #[test]
    fn missing_then_ref() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t.|0";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Missing, AlleleCall::Index(0)]
        );
    }

    // ---- Partially missing (0|.) ----

    #[test]
    fn ref_then_missing() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        let line = "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0|.";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Missing]
        );
    }

    // ---- ALT=. selected by genotype ----

    #[test]
    fn alt_dot_selected_becomes_n() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // ALT is '.', so the record has only the REF allele.
        // GT 0|1 selects allele index 1, which does not exist.
        let line = "chr1\t5\t.\tA\t.\t.\t.\t.\tGT\t0|1";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Missing]
        );
    }

    // ---- ALT=. but not selected ----

    #[test]
    fn alt_dot_not_selected_ok() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        let line = "chr1\t5\t.\tA\t.\t.\t.\t.\tGT\t0|0";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Index(0)]
        );
    }

    // ---- REF=. ----

    #[test]
    fn ref_dot() {
        let result = parse_raw_vcf_record("chr1\t5\t.\t.\tG\t.\t.\t.\tGT\t0|1");
        assert!(result.is_err(), "REF='.' must be rejected");
    }

    // ---- REF=. and ALT=. ----

    #[test]
    fn ref_and_alt_dot() {
        let result = parse_raw_vcf_record("chr1\t5\t.\t.\t.\t.\t.\t.\tGT\t0|1");
        assert!(result.is_err(), "REF='.' ALT='.' must be rejected");
    }

    // ---- Symbolic ALT selected ----

    #[test]
    fn symbolic_alt_selected_becomes_n() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        // <DEL> is symbolic. GT 0|1 selects it.
        let line = "chr1\t5\t.\tA\t<DEL>\t.\t.\t.\tGT\t0|1";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Missing]
        );
    }

    // ---- Symbolic ALT not selected ----

    #[test]
    fn symbolic_alt_not_selected_ok() {
        let (_ref_dir, ref_path, contig_len) = create_test_reference();
        let reference = open_reference(&ref_path);
        let args = test_args();
        let sample_names = vec!["S1".to_string()];
        let line = "chr1\t5\t.\tA\t<DEL>\t.\t.\t.\tGT\t0|0";
        let decoded =
            decode_line(line, "chr1", contig_len, &sample_names, &reference, &args).unwrap();
        assert_eq!(
            decoded.genotypes[0],
            vec![AlleleCall::Index(0), AlleleCall::Index(0)]
        );
    }
}