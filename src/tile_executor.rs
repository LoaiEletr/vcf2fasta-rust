//! CPU tile executor and GPU batch builder.
//!
//! The GPU batch builder caches three things across calls:
//!  * the faidx reference reader,
//!  * the tabix VCF reader,
//!  * the full contig sequence for each contig seen — shared across all
//!    producer threads via an `Arc<Mutex<HashMap<...>>>`, so N producers
//!    do not each hold N copies of a 50 MB contig.
//!
//! Without these caches, each tile re-reads the ~50 MB chr22 reference
//! and re-opens two index files — costs that dominate GPU wall time.

use crate::chunk::{Tile, TileResult};
use crate::cli::Args;
use crate::genotype::{AlleleCall, DecodedRecord};
use crate::report::ContigReport;
use crate::util::fetch_half_open;
use crate::vcf::{decode_and_validate_record, parse_raw_vcf_record};
use anyhow::{anyhow, Context, Result};
use rust_htslib::tbx::Read;
use rust_htslib::{faidx, tbx};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

// ---------------------------------------------------------------------------
// CPU tile executor
// ---------------------------------------------------------------------------

pub fn execute_tile(
    tile: &Tile,
    sample_names: &[String],
    sample_hap_offset: &[usize],
    max_ploidies: &[usize],
    pad: u64,
    args: &Args,
) -> Result<TileResult> {
    let mut report = ContigReport::new(tile.contig.clone());

    let reference = faidx::Reader::from_path(&args.reference)
        .with_context(|| format!("could not open reference {}", args.reference.display()))?;
    let contig_len = usize::try_from(reference.fetch_seq_len(&tile.contig))
        .map_err(|_| anyhow!("contig {} too long for usize", tile.contig))?;

    let base_start = tile.base_start as usize;
    let base_end = tile.base_end as usize;
    let query_start = (tile.base_start.saturating_sub(pad)) as usize;
    let hap_count = tile.hap_end - tile.hap_start;

    let mut per_hap: Vec<Vec<u8>> = (0..hap_count).map(|_| Vec::new()).collect();

    let mut vcf = tbx::Reader::from_path(&args.input)
        .with_context(|| format!("could not open {} as tabix-indexed VCF", args.input.display()))?;
    let tid = vcf
        .tid(&tile.contig)
        .with_context(|| format!("contig {} not in tabix index", tile.contig))?;
    vcf.fetch(tid, query_start as u64, base_end as u64)
        .with_context(|| "tabix fetch failed".to_string())?;

    let mut ref_cursor = base_start;
    let mut buf = Vec::<u8>::with_capacity(1024);
    let placeholder = args.no_call_string.as_deref().unwrap_or("N");

    while let Ok(has) = vcf.read(&mut buf) {
        if !has {
            break;
        }
        let line = match std::str::from_utf8(&buf) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let line = line.trim_end_matches('\r');
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[0] != tile.contig {
            continue;
        }

        let raw = match parse_raw_vcf_record(line) {
            Ok(r) => r,
            Err(e) => {
                report.warn(format!("[{}:{}] malformed: {}", tile.contig, fields[1], e));
                continue;
            }
        };
        let var_start = match raw.pos.checked_sub(1) {
            Some(p) => p,
            None => continue,
        };
        let var_end = var_start + raw.ref_allele.len();

        if var_start < base_start {
            if var_end > base_start && var_end > ref_cursor {
                ref_cursor = var_end;
            }
            continue;
        }
        if var_start >= base_end {
            continue;
        }
        report.seen += 1;

        let decoded: DecodedRecord = match decode_and_validate_record(
            &raw,
            &tile.contig,
            contig_len,
            sample_names,
            max_ploidies,
            &reference,
            args,
            &mut report,
        ) {
            Ok(d) => d,
            Err(e) => {
                report.warn(format!("[{}:{}] {}", tile.contig, raw.pos, e));
                continue;
            }
        };
        if decoded.start < ref_cursor {
            report.warn(format!(
                "[{}:{}] overlap/out-of-order (start={} < prev_end={})",
                tile.contig, raw.pos, decoded.start, ref_cursor
            ));
            continue;
        }
        if decoded.start > ref_cursor {
            let between = fetch_half_open(&reference, &tile.contig, ref_cursor, decoded.start)?;
            for h in 0..hap_count {
                per_hap[h].extend_from_slice(&between);
            }
        }
        for (s_idx, sample_gt) in decoded.genotypes.iter().enumerate() {
            let sample_off = sample_hap_offset[s_idx];
            let ploidy = max_ploidies[s_idx];
            let mut padded = sample_gt.clone();
            while padded.len() < ploidy {
                padded.push(AlleleCall::Reference);
            }
            if padded.len() > ploidy {
                padded.truncate(ploidy);
            }
            for global_h in sample_off..(sample_off + ploidy) {
                if global_h < tile.hap_start || global_h >= tile.hap_end {
                    continue;
                }
                let local_h = global_h - tile.hap_start;
                let local_idx = global_h - sample_off;
                match padded[local_idx] {
                    AlleleCall::Index(i) => {
                        per_hap[local_h].extend_from_slice(&decoded.alleles[i])
                    }
                    AlleleCall::Missing => {
                        per_hap[local_h].extend_from_slice(placeholder.as_bytes())
                    }
                    AlleleCall::Reference => {
                        let base = fetch_half_open(
                            &reference,
                            &tile.contig,
                            decoded.start,
                            decoded.start + 1,
                        )?;
                        per_hap[local_h].extend_from_slice(&base);
                    }
                }
            }
        }
        ref_cursor = decoded.end;
        report.applied += 1;
    }

    if ref_cursor < base_end {
        let tail = fetch_half_open(&reference, &tile.contig, ref_cursor, base_end)?;
        for h in 0..hap_count {
            per_hap[h].extend_from_slice(&tail);
        }
    }

    Ok(TileResult {
        contig: tile.contig.clone(),
        hap_start: tile.hap_start,
        hap_end: tile.hap_end,
        base_start: tile.base_start,
        base_end: tile.base_end,
        per_hap,
        seen: report.seen,
        applied: report.applied,
        warning_count: report.warning_count,
        warnings: report.warnings,
        warnings_by_reason: report.warnings_by_reason,
    })
}

// ---------------------------------------------------------------------------
// Cached host resources for the GPU producer
// ---------------------------------------------------------------------------

/// Shared cache of contig sequences, filled lazily. Keyed by contig name.
pub type ContigRefCache = Arc<Mutex<HashMap<String, Arc<Vec<u8>>>>>;

pub fn new_contig_ref_cache() -> ContigRefCache {
    Arc::new(Mutex::new(HashMap::new()))
}

pub struct HostResources {
    pub reference: faidx::Reader,
    pub vcf_path: PathBuf,
    pub vcf: tbx::Reader,
    pub contig_refs: ContigRefCache,
}

impl HostResources {
    pub fn new(
        reference_path: &std::path::Path,
        vcf_path: &std::path::Path,
        contig_refs: ContigRefCache,
    ) -> Result<Self> {
        let reference = faidx::Reader::from_path(reference_path)
            .with_context(|| format!("could not open reference {}", reference_path.display()))?;
        let vcf = tbx::Reader::from_path(vcf_path)
            .with_context(|| format!("could not open {} as tabix VCF", vcf_path.display()))?;
        Ok(Self {
            reference,
            vcf_path: vcf_path.to_path_buf(),
            vcf,
            contig_refs,
        })
    }

    pub fn contig_ref(&mut self, contig: &str) -> Result<Arc<Vec<u8>>> {
        {
            let guard = self
                .contig_refs
                .lock()
                .map_err(|_| anyhow!("contig_refs mutex poisoned"))?;
            if let Some(v) = guard.get(contig) {
                return Ok(v.clone());
            }
        }
        let len = self.reference.fetch_seq_len(contig);
        let seq = if len == 0 {
            Vec::new()
        } else {
            fetch_half_open(&self.reference, contig, 0, len as usize)?
        };
        let arc = Arc::new(seq);
        let mut guard = self
            .contig_refs
            .lock()
            .map_err(|_| anyhow!("contig_refs mutex poisoned"))?;
        let entry = guard
            .entry(contig.to_string())
            .or_insert_with(|| arc.clone());
        Ok(entry.clone())
    }

    pub fn reset_vcf(&mut self) -> Result<()> {
        self.vcf = tbx::Reader::from_path(&self.vcf_path)
            .with_context(|| format!("could not reopen {}", self.vcf_path.display()))?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// GPU batch builder
// ---------------------------------------------------------------------------

#[cfg(feature = "cuda")]
pub fn build_gpu_batch(
    tile: &Tile,
    host: &mut HostResources,
    sample_names: &[String],
    sample_hap_offset: &[usize],
    max_ploidies: &[usize],
    pad: u64,
    args: &Args,
) -> Result<(crate::gpu::GpuBatch, Arc<Vec<u8>>, crate::gpu::BatchStats)> {
    use crate::gpu::{BatchStats, GpuBatch};

    let contig_ref = host.contig_ref(&tile.contig)?;
    let contig_len = contig_ref.len();

    let base_start = tile.base_start as usize;
    let base_end = tile.base_end as usize;
    let query_start = (tile.base_start.saturating_sub(pad)) as usize;
    let hap_count = tile.hap_end - tile.hap_start;

    let mut variant_positions = Vec::<i32>::new();
    let mut variant_alleles = Vec::<u8>::new();
    let mut allele_lengths = Vec::<i32>::new();
    let mut allele_byte_offsets = Vec::<i32>::new();
    let mut allele_start_idx = Vec::<i32>::new();
    let mut num_alleles = Vec::<i32>::new();
    let placeholder = args.no_call_string.as_deref().unwrap_or("N").as_bytes();

    let mut per_variant_codes: Vec<Vec<i32>> = Vec::new();

    let tid = host
        .vcf
        .tid(&tile.contig)
        .with_context(|| format!("contig {} not in tabix index", tile.contig))?;
    host.vcf
        .fetch(tid, query_start as u64, base_end as u64)
        .with_context(|| "tabix fetch failed")?;

    let mut ref_cursor = base_start;
    let mut buf = Vec::<u8>::with_capacity(1024);
    let mut report = ContigReport::new(tile.contig.clone());
    let mut pending: Vec<(usize, Vec<Vec<u8>>, Vec<i32>)> = Vec::new();

    while let Ok(has) = host.vcf.read(&mut buf) {
        if !has {
            break;
        }
        let line = match std::str::from_utf8(&buf) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let line = line.trim_end_matches('\r');
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[0] != tile.contig {
            continue;
        }

        let raw = match parse_raw_vcf_record(line) {
            Ok(r) => r,
            Err(e) => {
                report.warn(format!("[{}:{}] malformed: {}", tile.contig, fields[1], e));
                continue;
            }
        };
        let var_start = match raw.pos.checked_sub(1) {
            Some(p) => p,
            None => continue,
        };
        let var_end = var_start + raw.ref_allele.len();

        if var_start < base_start {
            if var_end > base_start && var_end > ref_cursor {
                ref_cursor = var_end;
            }
            continue;
        }
        if var_start >= base_end {
            continue;
        }
        report.seen += 1;

        let decoded = match decode_and_validate_record(
            &raw,
            &tile.contig,
            contig_len,
            sample_names,
            max_ploidies,
            &host.reference,
            args,
            &mut report,
        ) {
            Ok(d) => d,
            Err(e) => {
                report.warn(format!("[{}:{}] {}", tile.contig, raw.pos, e));
                continue;
            }
        };
        if decoded.start < ref_cursor {
            report.warn(format!(
                "[{}:{}] overlap/out-of-order (start={} < prev_end={})",
                tile.contig, raw.pos, decoded.start, ref_cursor
            ));
            continue;
        }

        let mut codes = vec![-1i32; hap_count];
        for (s_idx, sample_gt) in decoded.genotypes.iter().enumerate() {
            let sample_off = sample_hap_offset[s_idx];
            let ploidy = max_ploidies[s_idx];
            let mut padded = sample_gt.clone();
            while padded.len() < ploidy {
                padded.push(AlleleCall::Reference);
            }
            if padded.len() > ploidy {
                padded.truncate(ploidy);
            }
            for global_h in sample_off..(sample_off + ploidy) {
                if global_h < tile.hap_start || global_h >= tile.hap_end {
                    continue;
                }
                let local_h = global_h - tile.hap_start;
                let local_idx = global_h - sample_off;
                codes[local_h] = match padded[local_idx] {
                    AlleleCall::Index(i) => i as i32,
                    AlleleCall::Missing => -2,
                    AlleleCall::Reference => -1,
                };
            }
        }
        pending.push((decoded.start, decoded.alleles.clone(), codes));
        ref_cursor = decoded.end;
        report.applied += 1;
    }

    let ref_start = {
        let mut cursor = base_start;
        if let Some((p, _, _)) = pending.first() {
            if *p < base_start {
                cursor = *p;
            }
        }
        cursor
    };

    for (pos, alleles, codes) in pending {
        if pos < ref_start {
            continue;
        }
        let start_idx = allele_lengths.len() as i32;
        allele_start_idx.push(start_idx);
        variant_positions.push(pos as i32);
        for a in &alleles {
            allele_lengths.push(a.len() as i32);
            allele_byte_offsets.push(variant_alleles.len() as i32);
            variant_alleles.extend_from_slice(a);
        }
        allele_lengths.push(placeholder.len() as i32);
        allele_byte_offsets.push(variant_alleles.len() as i32);
        variant_alleles.extend_from_slice(placeholder);
        num_alleles.push(alleles.len() as i32 + 1);
        per_variant_codes.push(codes);
    }

    // ---- Flatten to [variant][hap] layout ----
    //
    // Use `into_iter` so each intermediate row's allocation is freed
    // immediately after being copied into `genotype_indices`. This halves
    // the peak transient memory compared to `iter()` + per-element push,
    // which keeps the entire `per_variant_codes` vector alive while
    // `genotype_indices` grows to full size alongside it.
    let num_variants = variant_positions.len();
    let mut genotype_indices = Vec::<i32>::with_capacity(hap_count * num_variants);
    for codes in per_variant_codes.into_iter() {
        genotype_indices.extend(codes);
    }
    debug_assert_eq!(genotype_indices.len(), hap_count * num_variants);

    let batch = GpuBatch {
        contig: tile.contig.clone(),
        ref_start,
        ref_end: base_end,
        variant_positions,
        variant_alleles,
        allele_lengths,
        allele_byte_offsets,
        allele_start_idx,
        num_alleles,
        genotype_indices,
        num_haps: hap_count,
    };

    let stats = BatchStats {
        seen: report.seen,
        applied: report.applied,
        warning_count: report.warning_count,
        warnings: report.warnings,
        warnings_by_reason: report.warnings_by_reason,
    };
    Ok((batch, contig_ref, stats))
}