//! Output manager.
//!
//! Two modes:
//!
//! * **Per-contig (default):** each `(sample, contig, hap)` gets its own
//!   FASTA file, `<prefix><sample>_<contig>:<hap>.fa`.
//!
//! * **Merged (`--merged-output`):** one **multi-record** FASTA per
//!   `(sample, hap)`, `<prefix><sample>_<hap>.fa`, containing one
//!   `><sample>_<contig>:<hap>` record per contig in canonical order.
//!   Each contig's contribution is written to a `.part` file and the
//!   parts are concatenated at the end. The concatenation is a strict
//!   byte-for-byte join: `cat out_S_chr1:0.fa out_S_chr2:0.fa` in
//!   VCF-header order produces exactly `out_S_0.fa`.
//!
//! ## Bounded file handles
//!
//! Both modes stream results: each contig calls `write_segment` for
//! every `(sample, hap)` on every tile, so the number of distinct open
//! writers at any moment is `samples × haplotypes`, which is unbounded in
//! the input.
//!
//! To keep peak `open(2)` file descriptors flat, we cap the number of
//! simultaneously-open writers at `max_open` (default
//! [`DEFAULT_MAX_OPEN_WRITERS`]) and evict the least-recently-used one
//! when the cap is exceeded. Evicted writers are flushed, and their
//! column state is remembered so the next write can reopen them at the
//! correct position without duplicating the FASTA header.
//!
//! ## Memory profile of merged concatenation
//!
//! `finish()` streams each `.part` file with `std::io::copy` from a
//! `File` reader directly into the destination `File`. Peak RSS is
//! independent of contig length (a previous implementation used
//! `fs::read`, which materialized an entire contig in RAM per part).
//! Each part is removed immediately after concatenation, so scratch
//! usage shrinks as `finish()` progresses.
//!
//! ## Finalization and the double-newline bug
//!
//! There are two kinds of writer state at finalization time:
//!
//! * **Open writers** — currently in `self.open`. They hold the true
//!   column state internally, so calling `FastaWriter::finish()` on them
//!   produces the correct trailing newline (or none, if the writer
//!   happens to end on a perfect line boundary).
//! * **Evicted writers** — on disk, not in `self.open`, but tracked by
//!   an entry in `self.identities`. Their `column` field holds the value
//!   recorded at the moment of eviction. Finalizing them requires
//!   reopening the file at that column and calling `finish()`.
//!
//! A writer can be **both** — it was evicted at some point and then
//! reopened, so it is currently in `self.open` *and* still has an
//! identity entry. Its identity's `column` is **stale**: it was set at
//! eviction time and never updated while the writer was open. If we ran
//! both loops naively, we would finalize the file twice — the second
//! time at the stale column — appending an extra `\n`.
//!
//! The fix is to remove each key from `self.identities` as we finalize
//! it via `self.open`. Then the identity loop only sees writers that
//! are genuinely in the evicted-and-not-yet-finalized state.
//!
//! ## Reuse of a prefix across runs
//!
//! Output files are truncated on the **first write within a run** for
//! each `(contig, sample, hap)`, regardless of whether the file existed
//! on disk. This is what `WriterIdentity` tracks: its presence in
//! `self.identities` means "we have already opened this key at least
//! once during this process, so the file on disk contains our own
//! partial content and should be appended to." A file that exists but
//! has no corresponding identity was left over from a previous run and
//! is truncated on first touch.

use anyhow::{Context, Result};
use std::collections::{HashMap, VecDeque};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

use crate::fasta::FastaWriter;

/// Default cap on simultaneously-open output file handles.
///
/// Chosen to stay safely under typical `ulimit -n` values (commonly 1024)
/// while leaving room for input readers, reference FASTA, tabix indexes,
/// CUDA runtime handles, and log files.
pub const DEFAULT_MAX_OPEN_WRITERS: usize = 512;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputMode {
    PerContig,
    Merged,
}

type Key = (String, usize, usize); // (contig, sample_idx, hap)

/// Per-key state that survives eviction and reopen.
struct WriterIdentity {
    path: PathBuf,
    seq_name: String,
    /// Last known line-wrap column. Meaningful for both modes: merged
    /// parts are wrapped FastaWriter output, so eviction may leave them
    /// mid-line just like per-contig files.
    column: usize,
}

pub struct OutputManager {
    mode: OutputMode,
    prefix: String,
    line_width: usize,
    scratch: Option<tempfile::TempDir>,
    contig_order: Vec<String>,

    /// Ordered parts per `(sample, hap)`. Merged mode only.
    parts: HashMap<(usize, usize), Vec<(String, PathBuf)>>,

    /// Number of distinct output files produced so far (PerContig mode).
    files_created: usize,

    /// Identity for every key that has ever been written **during this
    /// run**. Presence in this map is the "we already own this file"
    /// signal; a file on disk with no matching identity is stale from a
    /// previous run and will be truncated on first write.
    ///
    /// Entries are removed when their writer is finalized, so the map
    /// only ever contains keys that still need action at finalization
    /// time.
    identities: HashMap<Key, WriterIdentity>,

    /// Currently-open writers. Bounded by `max_open`. Both modes use
    /// the same wrapping `FastaWriter` so merged parts are proper
    /// multi-record FASTA fragments.
    open: HashMap<Key, FastaWriter>,

    /// LRU order; front is least recently used, back is most recently.
    lru: VecDeque<Key>,

    max_open: usize,
}

impl OutputManager {
    pub fn new(
        mode: OutputMode,
        prefix: String,
        line_width: usize,
        contig_order: Vec<String>,
    ) -> Result<Self> {
        Self::with_max_open(
            mode,
            prefix,
            line_width,
            contig_order,
            DEFAULT_MAX_OPEN_WRITERS,
        )
    }

    pub fn with_max_open(
        mode: OutputMode,
        prefix: String,
        line_width: usize,
        contig_order: Vec<String>,
        max_open: usize,
    ) -> Result<Self> {
        let scratch = if mode == OutputMode::Merged {
            Some(
                tempfile::Builder::new()
                    .prefix("vcf2fasta_parts_")
                    .tempdir()
                    .context("could not create output scratch dir")?,
            )
        } else {
            None
        };
        Ok(Self {
            mode,
            prefix,
            line_width,
            scratch,
            contig_order,
            parts: HashMap::new(),
            files_created: 0,
            identities: HashMap::new(),
            open: HashMap::new(),
            lru: VecDeque::new(),
            max_open: max_open.max(1),
        })
    }

    /// Append one tile's bytes to the correct output stream.
    ///
    /// Safe to call many times for the same `(contig, sample, hap)`; the
    /// underlying writer is opened once, may be evicted under memory
    /// pressure, and is transparently reopened if needed.
    pub fn write_segment(
        &mut self,
        contig: &str,
        sample_idx: usize,
        sample_name: &str,
        hap: usize,
        bytes: &[u8],
    ) -> Result<()> {
        let key: Key = (contig.to_string(), sample_idx, hap);

        if !self.open.contains_key(&key) {
            self.evict_lru_if_needed()?;
            self.open_key(&key, contig, sample_idx, sample_name, hap)?;
        }

        // Move key to most-recently-used.
        if let Some(pos) = self.lru.iter().position(|k| k == &key) {
            self.lru.remove(pos);
        }
        self.lru.push_back(key.clone());

        self.open
            .get_mut(&key)
            .expect("writer was just opened")
            .write_seq(bytes)?;
        Ok(())
    }

    /// Open (or reopen) the writer for `key`.
    ///
    /// Both modes use the same header string `<sample>_<contig>:<hap>`;
    /// only the destination path differs. Merged parts are written as
    /// wrapped FASTA records, so the merged file is a strict
    /// concatenation of per-contig files.
    fn open_key(
        &mut self,
        key: &Key,
        contig: &str,
        sample_idx: usize,
        sample_name: &str,
        hap: usize,
    ) -> Result<()> {
        // `is_first_write` means "this is the first time in this process
        // we are writing to this key". Determined from our own bookkeeping,
        // not from the filesystem.
        let is_first_write = !self.identities.contains_key(key);

        let (path, seq_name, column) = if let Some(id) = self.identities.get(key) {
            (id.path.clone(), id.seq_name.clone(), id.column)
        } else {
            // Same header in both modes: <sample>_<contig>:<hap>.
            let seq_name = format!("{}_{}:{}", sample_name, contig, hap);
            let path = match self.mode {
                OutputMode::PerContig => {
                    PathBuf::from(format!("{}{}.fa", self.prefix, seq_name))
                }
                OutputMode::Merged => {
                    let scratch = self.scratch.as_ref().expect("merged needs scratch");
                    // One part per (sample, hap, contig). Sanitize the
                    // contig name for filesystem safety.
                    scratch.path().join(format!(
                        "p_s{}_h{}_c{}.part",
                        sample_idx,
                        hap,
                        sanitize(contig)
                    ))
                }
            };
            (path, seq_name, 0)
        };

        // Build the writer. Both modes use FastaWriter::create on first
        // write (truncate + write header) and FastaWriter::reopen on
        // subsequent writes (append at the recorded column).
        let writer = if is_first_write {
            FastaWriter::create(&path, &seq_name, self.line_width)?
        } else {
            FastaWriter::reopen(&path, self.line_width, column)?
        };

        if is_first_write {
            // Register the part in canonical contig order (once).
            if self.mode == OutputMode::Merged {
                let entry = self.parts.entry((sample_idx, hap)).or_default();
                if !entry.iter().any(|(c, _)| c == contig) {
                    entry.push((contig.to_string(), path.clone()));
                    let order: HashMap<&str, usize> = self
                        .contig_order
                        .iter()
                        .enumerate()
                        .map(|(i, c)| (c.as_str(), i))
                        .collect();
                    entry.sort_by_key(|(c, _)| {
                        order.get(c.as_str()).copied().unwrap_or(usize::MAX)
                    });
                }
            }

            self.identities.insert(
                key.clone(),
                WriterIdentity {
                    path,
                    seq_name,
                    column: 0,
                },
            );
            if self.mode == OutputMode::PerContig {
                self.files_created += 1;
            }
        }

        self.open.insert(key.clone(), writer);
        Ok(())
    }

    fn evict_lru_if_needed(&mut self) -> Result<()> {
        while self.open.len() >= self.max_open {
            let Some(old) = self.lru.pop_front() else {
                break;
            };
            if let Some(mut w) = self.open.remove(&old) {
                let col = w.column();
                w.flush()?;
                if let Some(id) = self.identities.get_mut(&old) {
                    id.column = col;
                }
            }
        }
        Ok(())
    }

    /// Close all writers belonging to `contig`. Safe to call multiple
    /// times; safe to call when nothing is open for `contig`.
    ///
    /// In Merged mode, the `.part` files for this contig remain on disk
    /// and are concatenated into the per-haplotype output at the final
    /// [`finish`].
    pub fn finish_contig(&mut self, contig: &str) -> Result<()> {
        // Finalize currently-open writers, and remove their identities
        // so the loop below does not re-finalize them.
        //
        // A writer that was evicted at some point and then reopened has
        // a stale `column` in its identity (only updated on eviction,
        // not on every write). If we left the identity in place, the
        // second loop would reopen and call `finish()` on a file that
        // has already been finalized, appending a spurious trailing
        // newline.
        let open_keys: Vec<Key> = self
            .open
            .keys()
            .filter(|(c, _, _)| c == contig)
            .cloned()
            .collect();
        for key in &open_keys {
            if let Some(w) = self.open.remove(key) {
                w.finish()?;
            }
            // The writer has been finalized, so its identity must not
            // survive into the second loop.
            self.identities.remove(key);
            self.lru.retain(|k| k != key);
        }

        // Finalize evicted writers that still have a non-zero column.
        // Reopening in append mode and calling `finish()` adds only the
        // trailing newline; it does not rewrite the file. This applies
        // to Merged parts too — a part left mid-line by an eviction
        // needs its final newline before concatenation.
        //
        // Writers that were evicted at a perfect line boundary
        // (`column == 0`) do not need re-finalization: `write_seq`
        // already wrote the newline when the column reached
        // `line_width`.
        let evicted_keys: Vec<Key> = self
            .identities
            .keys()
            .filter(|(c, _, _)| c == contig)
            .cloned()
            .collect();
        for key in evicted_keys {
            if let Some(id) = self.identities.remove(&key) {
                if id.column != 0 {
                    FastaWriter::reopen(&id.path, self.line_width, id.column)?.finish()?;
                }
            }
        }
        Ok(())
    }

    /// Finalize all output. Returns the number of output files produced.
    ///
    /// In PerContig mode this is the number of distinct `(sample, contig,
    /// hap)` writers ever opened during this run. In Merged mode it is
    /// the number of concatenated `(sample, hap)` files.
    ///
    /// **Memory:** the merged concatenation uses `std::io::copy` from
    /// `File` to `File`, so peak RSS is independent of contig length.
    pub fn finish(mut self, sample_names: &[String]) -> Result<usize> {
        // Close any writers still open, and remove their identities so
        // the loop below does not re-finalize them. Same reasoning as
        // in `finish_contig`: a writer that was evicted and reopened
        // has a stale `column` in its identity, and re-finalizing it
        // would append a spurious trailing newline.
        let open_keys: Vec<Key> = self.open.keys().cloned().collect();
        for key in &open_keys {
            if let Some(w) = self.open.remove(key) {
                w.finish()?;
            }
            self.identities.remove(key);
        }

        // Finalize any evicted writers.
        let identity_keys: Vec<Key> = self.identities.keys().cloned().collect();
        for key in identity_keys {
            if let Some(id) = self.identities.remove(&key) {
                if id.column != 0 {
                    FastaWriter::reopen(&id.path, self.line_width, id.column)?.finish()?;
                }
            }
        }

        match self.mode {
            OutputMode::PerContig => Ok(self.files_created),
            OutputMode::Merged => {
                let mut written = 0usize;
                for ((sample_idx, hap), parts) in self.parts {
                    let sample = sample_names
                        .get(sample_idx)
                        .map(|s| s.as_str())
                        .unwrap_or("unknown");
                    let file_name = format!("{}{}_{}.fa", self.prefix, sample, hap);

                    // Destination is a bare `File`, not `BufWriter<File>`:
                    // on Linux `std::io::copy` between two files issues
                    // `copy_file_range`, which never enters user space.
                    let mut dst = File::create(&file_name)
                        .with_context(|| format!("could not create {}", file_name))?;

                    // Each part is a complete FASTA record
                    // (`><sample>_<contig>:<hap>` + wrapped body).
                    // Concatenating them verbatim yields a valid
                    // multi-record FASTA in canonical contig order.
                    for (_, part) in &parts {
                        let mut src = File::open(part).with_context(|| {
                            format!("could not open scratch part {}", part.display())
                        })?;
                        std::io::copy(&mut src, &mut dst).with_context(|| {
                            format!(
                                "could not concatenate {} into {}",
                                part.display(),
                                file_name
                            )
                        })?;
                    }
                    dst.flush()?;
                    drop(dst);

                    // Free scratch as we go. Non-fatal on failure: the
                    // scratch dir is a TempDir that will be cleaned up
                    // on drop anyway.
                    for (_, part) in &parts {
                        let _ = fs::remove_file(part);
                    }

                    written += 1;
                }
                Ok(written)
            }
        }
    }
}

/// Replace characters that are not safe in a filename with `_`.
///
/// Contig names can contain `:` (e.g. `chr1:1000-2000`) or `/` in
/// unusual references; sanitizing keeps the scratch part filenames
/// portable across filesystems.
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

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn per_contig_mode_writes_one_file_per_segment() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let mut mgr = OutputManager::new(
            OutputMode::PerContig,
            prefix.clone(),
            80,
            vec!["chr1".to_string()],
        )
        .unwrap();
        mgr.write_segment("chr1", 0, "S1", 0, b"ACGT").unwrap();
        mgr.write_segment("chr1", 0, "S1", 1, b"TGCA").unwrap();
        mgr.finish_contig("chr1").unwrap();
        let n = mgr.finish(&["S1".to_string()]).unwrap();
        assert_eq!(n, 2);
        assert!(dir.path().join("S1_chr1:0.fa").exists());
        assert!(dir.path().join("S1_chr1:1.fa").exists());
    }

    #[test]
    fn per_contig_mode_appends_multiple_tiles_for_same_key() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let mut mgr = OutputManager::new(
            OutputMode::PerContig,
            prefix.clone(),
            80,
            vec!["chr1".to_string()],
        )
        .unwrap();
        mgr.write_segment("chr1", 0, "S1", 0, b"AAAA").unwrap();
        mgr.write_segment("chr1", 0, "S1", 0, b"CCCC").unwrap();
        mgr.write_segment("chr1", 0, "S1", 0, b"GGGG").unwrap();
        mgr.finish_contig("chr1").unwrap();
        mgr.finish(&["S1".to_string()]).unwrap();

        let content = fs::read_to_string(dir.path().join("S1_chr1:0.fa")).unwrap();
        assert_eq!(content, ">S1_chr1:0\nAAAACCCCGGGG\n");
    }

    /// Regression test for the doubling bug: a second run with the same
    /// prefix must truncate stale files, not append to them.
    #[test]
    fn second_run_truncates_stale_files() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());

        {
            let mut mgr = OutputManager::new(
                OutputMode::PerContig,
                prefix.clone(),
                80,
                vec!["chr1".to_string()],
            )
            .unwrap();
            mgr.write_segment("chr1", 0, "S1", 0, b"AAAA").unwrap();
            mgr.finish_contig("chr1").unwrap();
            mgr.finish(&["S1".to_string()]).unwrap();
        }

        let path = dir.path().join("S1_chr1:0.fa");
        assert_eq!(fs::read_to_string(&path).unwrap(), ">S1_chr1:0\nAAAA\n");

        {
            let mut mgr = OutputManager::new(
                OutputMode::PerContig,
                prefix,
                80,
                vec!["chr1".to_string()],
            )
            .unwrap();
            mgr.write_segment("chr1", 0, "S1", 0, b"CCCC").unwrap();
            mgr.finish_contig("chr1").unwrap();
            mgr.finish(&["S1".to_string()]).unwrap();
        }

        assert_eq!(fs::read_to_string(&path).unwrap(), ">S1_chr1:0\nCCCC\n");
    }

    #[test]
    fn per_contig_file_count_is_stable_after_finish_contig() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let mut mgr = OutputManager::new(
            OutputMode::PerContig,
            prefix,
            80,
            vec!["chr1".to_string(), "chr2".to_string()],
        )
        .unwrap();

        mgr.write_segment("chr1", 0, "S1", 0, b"AAAA").unwrap();
        mgr.write_segment("chr1", 0, "S1", 1, b"CCCC").unwrap();
        mgr.finish_contig("chr1").unwrap();

        mgr.write_segment("chr2", 0, "S1", 0, b"GGGG").unwrap();
        mgr.finish_contig("chr2").unwrap();

        let n = mgr.finish(&["S1".to_string()]).unwrap();
        assert_eq!(n, 3);
    }

    /// Regression test for the double-finalization bug.
    ///
    /// Force evictions by setting `max_open` below the number of
    /// distinct keys. Each eviction records a column in the identity
    /// map. Subsequent writes to an evicted key reopen it and leave it
    /// in `self.open`. At finalization time, both the open-writer loop
    /// and the identity loop would otherwise try to finalize the same
    /// file, producing an extra trailing newline.
    #[test]
    fn lru_eviction_preserves_all_bytes() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let mut mgr = OutputManager::with_max_open(
            OutputMode::PerContig,
            prefix.clone(),
            80,
            vec!["chr1".to_string()],
            2,
        )
        .unwrap();

        for tile in 0..5 {
            for s in 0..5 {
                mgr.write_segment(
                    "chr1",
                    s,
                    &format!("S{}", s),
                    0,
                    format!("t{}s{}_", tile, s).as_bytes(),
                )
                .unwrap();
            }
        }
        mgr.finish_contig("chr1").unwrap();

        let names: Vec<String> = (0..5).map(|s| format!("S{}", s)).collect();
        let n = mgr.finish(&names).unwrap();
        assert_eq!(n, 5);

        for s in 0..5 {
            let content =
                fs::read_to_string(dir.path().join(format!("S{}_chr1:0.fa", s))).unwrap();
            let expected: String = (0..5).map(|t| format!("t{}s{}_", t, s)).collect();
            assert_eq!(
                content,
                format!(">S{}_chr1:0\n{}\n", s, expected),
                "S{} has a spurious or missing trailing newline",
                s
            );
        }
    }

    /// Merged mode produces one record per contig, in canonical order.
    #[test]
    fn merged_mode_produces_multi_record_fasta_in_canonical_order() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let mut mgr = OutputManager::new(
            OutputMode::Merged,
            prefix.clone(),
            80,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();
        // Write out of order — the final file must still be in canonical order.
        mgr.write_segment("chr3", 0, "S1", 0, b"CCC").unwrap();
        mgr.finish_contig("chr3").unwrap();
        mgr.write_segment("chr1", 0, "S1", 0, b"AAA").unwrap();
        mgr.finish_contig("chr1").unwrap();
        mgr.write_segment("chr2", 0, "S1", 0, b"BBB").unwrap();
        mgr.finish_contig("chr2").unwrap();

        let n = mgr.finish(&["S1".to_string()]).unwrap();
        assert_eq!(n, 1);

        let content = fs::read_to_string(dir.path().join("S1_0.fa")).unwrap();
        assert_eq!(
            content,
            ">S1_chr1:0\nAAA\n>S1_chr2:0\nBBB\n>S1_chr3:0\nCCC\n"
        );
    }

    #[test]
    fn merged_mode_appends_multiple_tiles_per_contig() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let mut mgr = OutputManager::new(
            OutputMode::Merged,
            prefix.clone(),
            80,
            vec!["chr1".to_string()],
        )
        .unwrap();
        mgr.write_segment("chr1", 0, "S1", 0, b"AAAA").unwrap();
        mgr.write_segment("chr1", 0, "S1", 0, b"CCCC").unwrap();
        mgr.finish_contig("chr1").unwrap();
        let n = mgr.finish(&["S1".to_string()]).unwrap();
        assert_eq!(n, 1);
        let content = fs::read_to_string(dir.path().join("S1_0.fa")).unwrap();
        assert_eq!(content, ">S1_chr1:0\nAAAACCCC\n");
    }

    #[test]
    fn merged_mode_wraps_long_sequences() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let mut mgr = OutputManager::new(
            OutputMode::Merged,
            prefix.clone(),
            4,
            vec!["chr1".to_string()],
        )
        .unwrap();
        mgr.write_segment("chr1", 0, "S1", 0, b"ACGTACGTACGT").unwrap();
        mgr.finish_contig("chr1").unwrap();
        mgr.finish(&["S1".to_string()]).unwrap();
        let content = fs::read_to_string(dir.path().join("S1_0.fa")).unwrap();
        assert_eq!(content, ">S1_chr1:0\nACGT\nACGT\nACGT\n");
    }

    /// The property that `run_merged_output_check.sh` relies on:
    /// the merged file is byte-identical to the concatenation of the
    /// per-contig files, in canonical contig order.
    #[test]
    fn merged_mode_concatenation_equals_cat_of_per_contig() {
        let dir = tempdir().unwrap();
        let base = dir.path();
        let contigs = vec!["chr1".to_string(), "chr2".to_string()];

        {
            let mut mgr = OutputManager::new(
                OutputMode::PerContig,
                format!("{}/pc_", base.display()),
                80,
                contigs.clone(),
            )
            .unwrap();
            mgr.write_segment("chr1", 0, "S1", 0, b"AAAA").unwrap();
            mgr.finish_contig("chr1").unwrap();
            mgr.write_segment("chr2", 0, "S1", 0, b"CCCC").unwrap();
            mgr.finish_contig("chr2").unwrap();
            mgr.finish(&["S1".to_string()]).unwrap();
        }

        {
            let mut mgr = OutputManager::new(
                OutputMode::Merged,
                format!("{}/mg_", base.display()),
                80,
                contigs.clone(),
            )
            .unwrap();
            mgr.write_segment("chr1", 0, "S1", 0, b"AAAA").unwrap();
            mgr.finish_contig("chr1").unwrap();
            mgr.write_segment("chr2", 0, "S1", 0, b"CCCC").unwrap();
            mgr.finish_contig("chr2").unwrap();
            mgr.finish(&["S1".to_string()]).unwrap();
        }

        let pc1 = fs::read_to_string(base.join("pc_S1_chr1:0.fa")).unwrap();
        let pc2 = fs::read_to_string(base.join("pc_S1_chr2:0.fa")).unwrap();
        let merged = fs::read_to_string(base.join("mg_S1_0.fa")).unwrap();

        assert_eq!(merged, format!("{}{}", pc1, pc2));
    }

    /// Merged mode with forced evictions. Exercises the same
    /// double-finalization path as the per-contig version, but through
    /// the merged part files.
    #[test]
    fn merged_mode_lru_eviction_preserves_all_bytes() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let mut mgr = OutputManager::with_max_open(
            OutputMode::Merged,
            prefix.clone(),
            80,
            vec!["chr1".to_string()],
            2,
        )
        .unwrap();

        for tile in 0..5 {
            for s in 0..5 {
                mgr.write_segment(
                    "chr1",
                    s,
                    &format!("S{}", s),
                    0,
                    format!("t{}s{}_", tile, s).as_bytes(),
                )
                .unwrap();
            }
        }
        mgr.finish_contig("chr1").unwrap();

        let names: Vec<String> = (0..5).map(|s| format!("S{}", s)).collect();
        let n = mgr.finish(&names).unwrap();
        assert_eq!(n, 5);

        for s in 0..5 {
            let content = fs::read_to_string(dir.path().join(format!("S{}_0.fa", s))).unwrap();
            let expected: String = (0..5).map(|t| format!("t{}s{}_", t, s)).collect();
            assert_eq!(
                content,
                format!(">S{}_chr1:0\n{}\n", s, expected),
                "S{} has a spurious or missing trailing newline",
                s
            );
        }
    }

    #[test]
    fn finish_contig_is_idempotent_and_handles_unknown_contig() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let mut mgr = OutputManager::new(
            OutputMode::PerContig,
            prefix,
            80,
            vec!["chr1".to_string()],
        )
        .unwrap();
        mgr.finish_contig("chr1").unwrap();
        mgr.finish_contig("chrX").unwrap();
        mgr.write_segment("chr1", 0, "S1", 0, b"ACGT").unwrap();
        mgr.finish_contig("chr1").unwrap();
        mgr.finish_contig("chr1").unwrap();
        mgr.finish(&["S1".to_string()]).unwrap();
    }

    #[test]
    fn sanitize_replaces_unsafe_chars() {
        assert_eq!(sanitize("chr1"), "chr1");
        assert_eq!(sanitize("chr1:1000-2000"), "chr1_1000-2000");
        assert_eq!(sanitize("chr/1"), "chr_1");
    }
}