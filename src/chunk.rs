//! Genomic chunk planning, two-dimensional tiling, and per-tile results.
//!
//! 1-D [`GenomicChunk`] planning is retained for backwards compatibility.
//! New code should use [`Tile`] / [`plan_tiles`], which split both the
//! haplotype dimension and the base-range dimension.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

pub const MIN_CHUNK_SIZE: u64 = 100_000;
pub const MAX_CHUNK_SIZE: u64 = 50_000_000;
pub const TARGET_CHUNKS_PER_THREAD: u64 = 4;
pub const MIN_TOTAL_CHUNKS: u64 = 8;
pub const MIN_CHUNK_PAD: u64 = 1_000;

// ---------------------------------------------------------------------------
// 1-D chunk (kept for back-compat)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct GenomicChunk {
    pub contig: String,
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone)]
pub struct ChunkPlanInfo {
    pub total_bases: u64,
    pub chunk_size: u64,
    pub auto_chunk_size: bool,
    pub target_total_chunks: u64,
    pub actual_total_chunks: usize,
    pub per_contig_chunks: HashMap<String, usize>,
}

#[derive(Debug)]
pub struct ChunkPlan {
    pub chunks: Vec<GenomicChunk>,
    pub info: ChunkPlanInfo,
}

#[derive(Debug)]
pub struct ChunkResult {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub per_sample: Vec<Vec<Vec<u8>>>,
    pub seen: usize,
    pub applied: usize,
    pub skipped: usize,
    pub output_files: usize,
    pub warning_count: usize,
    pub warnings: Vec<String>,
}

pub fn plan_chunks(
    contigs: &[String],
    contig_lengths: &HashMap<String, u64>,
    threads: usize,
    user_chunk_size: Option<u64>,
) -> ChunkPlan {
    let total_bases: u64 = contigs
        .iter()
        .filter_map(|c| contig_lengths.get(c).copied())
        .sum();

    let threads_u64 = threads.max(1) as u64;
    let target_total_chunks =
        (threads_u64 * TARGET_CHUNKS_PER_THREAD).max(MIN_TOTAL_CHUNKS);

    let (chunk_size, auto_chunk_size) = match user_chunk_size {
        Some(cs) => (cs.max(1), false),
        None => {
            let raw = if total_bases > 0 {
                total_bases / target_total_chunks
            } else {
                10_000_000
            };
            (raw.clamp(MIN_CHUNK_SIZE, MAX_CHUNK_SIZE), true)
        }
    };

    let (chunks, per_contig_chunks) =
        plan_chunks_with_size(contigs, contig_lengths, chunk_size);

    ChunkPlan {
        chunks,
        info: ChunkPlanInfo {
            total_bases,
            chunk_size,
            auto_chunk_size,
            target_total_chunks,
            actual_total_chunks: per_contig_chunks.values().sum(),
            per_contig_chunks,
        },
    }
}

pub fn auto_chunk_pad(input: &Path) -> Result<u64> {
    let reader = open_maybe_gzip(input)?;
    let mut max_ref_len: usize = 0;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let mut it = line.split('\t');
        let _chrom = it.next();
        let _pos = it.next();
        let _id = it.next();
        if let Some(ref_allele) = it.next() {
            if ref_allele.len() > max_ref_len {
                max_ref_len = ref_allele.len();
            }
        }
    }
    Ok(((max_ref_len as u64) * 2).max(MIN_CHUNK_PAD))
}

fn plan_chunks_with_size(
    contigs: &[String],
    contig_lengths: &HashMap<String, u64>,
    chunk_size: u64,
) -> (Vec<GenomicChunk>, HashMap<String, usize>) {
    let chunk_size = chunk_size.max(1);
    let mut chunks = Vec::new();
    let mut per_contig: HashMap<String, usize> = HashMap::new();

    for contig in contigs {
        let len = match contig_lengths.get(contig) {
            Some(&l) => l,
            None => continue,
        };
        if len == 0 {
            continue;
        }
        let mut start: u64 = 0;
        let mut count = 0usize;
        while start < len {
            let end = (start + chunk_size).min(len);
            chunks.push(GenomicChunk {
                contig: contig.clone(),
                start,
                end,
            });
            count += 1;
            start = end;
        }
        per_contig.insert(contig.clone(), count);
    }

    (chunks, per_contig)
}

// ---------------------------------------------------------------------------
// Two-dimensional tiles
// ---------------------------------------------------------------------------

/// A 2-D work unit: a haplotype block × a base range within one contig.
///
/// * `hap_start..hap_end` are indices into the *global* haplotype list
///   (`sum of per-sample observed ploidy`).
/// * `base_start..base_end` are 0-based half-open reference coordinates.
#[derive(Debug, Clone)]
pub struct Tile {
    pub contig: String,
    pub hap_start: usize,
    pub hap_end: usize,
    pub base_start: u64,
    pub base_end: u64,
}

impl Tile {
    pub fn hap_count(&self) -> usize {
        self.hap_end.saturating_sub(self.hap_start)
    }
    pub fn base_len(&self) -> u64 {
        self.base_end.saturating_sub(self.base_start)
    }
}

/// Plans 2-D tiles across all contigs.
pub fn plan_tiles(
    contigs: &[String],
    contig_lengths: &HashMap<String, u64>,
    haplotype_count: usize,
    hap_block_size: usize,
    base_block_size: u64,
) -> Vec<Tile> {
    let hap_block_size = hap_block_size.max(1);
    let base_block_size = base_block_size.max(1);
    let mut tiles = Vec::new();
    if haplotype_count == 0 {
        return tiles;
    }

    for contig in contigs {
        let len = match contig_lengths.get(contig) {
            Some(&l) if l > 0 => l,
            _ => continue,
        };

        let mut h = 0usize;
        while h < haplotype_count {
            let h_end = (h + hap_block_size).min(haplotype_count);
            let mut b = 0u64;
            while b < len {
                let b_end = (b + base_block_size).min(len);
                tiles.push(Tile {
                    contig: contig.clone(),
                    hap_start: h,
                    hap_end: h_end,
                    base_start: b,
                    base_end: b_end,
                });
                b = b_end;
            }
            h = h_end;
        }
    }
    tiles
}

/// Result of executing one tile (or a merged batch of adjacent tiles).
#[derive(Debug)]
pub struct TileResult {
    pub contig: String,
    pub hap_start: usize,
    pub hap_end: usize,
    pub base_start: u64,
    pub base_end: u64,
    /// `per_hap[local_hap_idx]`, `local_hap_idx = global_hap - hap_start`.
    pub per_hap: Vec<Vec<u8>>,
    pub seen: usize,
    pub applied: usize,
    /// Total number of warnings emitted while processing this tile.
    /// Unbounded. Use this for the "total" number.
    pub warning_count: usize,
    /// Up to `MAX_WARNINGS_PER_CONTIG` stored warning strings.
    pub warnings: Vec<String>,
    /// Per-reason histogram. Unbounded — every warning is counted.
    /// Sum of values equals `warning_count`.
    pub warnings_by_reason: BTreeMap<&'static str, usize>,
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead>> {
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

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn lengths(pairs: &[(&str, u64)]) -> HashMap<String, u64> {
        pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    #[test]
    fn single_contig_auto_chunks() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 100_000_000)]);
        let plan = plan_chunks(&contigs, &lens, 8, None);
        assert!(plan.info.auto_chunk_size);
        assert_eq!(plan.info.chunk_size, 3_125_000);
        assert_eq!(plan.chunks.len(), 32);
    }
    #[test]
    fn many_chunks_per_worker_clamped_by_max() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 1_000_000_000)]);
        let plan = plan_chunks(&contigs, &lens, 4, None);
        assert_eq!(plan.info.chunk_size, MAX_CHUNK_SIZE);
        assert_eq!(plan.chunks.len(), 20);
    }
    #[test]
    fn tiny_genome_clamped_by_min() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 10_000)]);
        let plan = plan_chunks(&contigs, &lens, 8, None);
        assert_eq!(plan.info.chunk_size, MIN_CHUNK_SIZE);
        assert_eq!(plan.chunks.len(), 1);
    }
    #[test]
    fn user_chunk_size_overrides_auto() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 100_000_000)]);
        let plan = plan_chunks(&contigs, &lens, 8, Some(1_000_000));
        assert!(!plan.info.auto_chunk_size);
        assert_eq!(plan.chunks.len(), 100);
    }
    #[test]
    fn chunk_bounds_are_half_open() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 25)]);
        let plan = plan_chunks(&contigs, &lens, 1, Some(10));
        assert_eq!(plan.chunks.len(), 3);
        assert_eq!((plan.chunks[0].start, plan.chunks[0].end), (0, 10));
        assert_eq!((plan.chunks[2].start, plan.chunks[2].end), (20, 25));
    }

    #[test]
    fn plan_tiles_single_tile() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 10_000)]);
        let tiles = plan_tiles(&contigs, &lens, 2, 100, 100_000);
        assert_eq!(tiles.len(), 1);
        assert_eq!((tiles[0].hap_start, tiles[0].hap_end), (0, 2));
        assert_eq!((tiles[0].base_start, tiles[0].base_end), (0, 10_000));
    }
    #[test]
    fn plan_tiles_splits_haplotypes() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 1_000_000)]);
        let tiles = plan_tiles(&contigs, &lens, 1000, 250, 1_000_000);
        assert_eq!(tiles.len(), 4);
    }
    #[test]
    fn plan_tiles_splits_base_ranges() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 25)]);
        let tiles = plan_tiles(&contigs, &lens, 1, 100, 10);
        assert_eq!(tiles.len(), 3);
        assert_eq!((tiles[2].base_start, tiles[2].base_end), (20, 25));
    }
    #[test]
    fn plan_tiles_splits_both_dimensions() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 100)]);
        let tiles = plan_tiles(&contigs, &lens, 1000, 250, 20);
        assert_eq!(tiles.len(), 4 * 5);
    }
    #[test]
    fn plan_tiles_no_giant_allocation() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 10_000_000)]);
        let tiles = plan_tiles(&contigs, &lens, 2000, 500, 1_000_000);
        assert!(tiles.len() >= 4 * 10);
    }
    #[test]
    fn plan_tiles_empty_haps_is_empty() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 1_000)]);
        assert!(plan_tiles(&contigs, &lens, 0, 10, 10).is_empty());
    }
    #[test]
    fn plan_tiles_reproduces_1d_chunks_for_single_hap_block() {
        let contigs = vec!["chr1".to_string()];
        let lens = lengths(&[("chr1", 100_000_000)]);
        let plan = plan_chunks(&contigs, &lens, 8, Some(1_000_000));
        let tiles = plan_tiles(&contigs, &lens, 4, 4, 1_000_000);
        assert_eq!(tiles.len(), plan.chunks.len());
        for (t, c) in tiles.iter().zip(plan.chunks.iter()) {
            assert_eq!(t.contig, c.contig);
            assert_eq!(t.base_start, c.start);
            assert_eq!(t.base_end, c.end);
        }
    }
}