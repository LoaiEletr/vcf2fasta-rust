//! Locate per-sample alignment files inside a `--bam-dir`.
//!
//! Sequencing pipelines rarely produce files named exactly `<SAMPLE>.bam`.
//! Common variants include:
//!
//! * `SAMPLE001.sorted.bam`
//! * `SAMPLE001.markdup.bam`
//! * `SAMPLE001.sorted.markdup.bam`
//! * `SAMPLE001.bwa-mem.bam`
//! * `SAMPLE001.cram` (with `.crai` instead of `.bai`)
//!
//! The resolver accepts any file whose name starts with the sample name as
//! a **leading token** — i.e. followed by a non-alphanumeric character or
//! the file extension — and whose extension is `.bam` or `.cram`.
//!
//! ## Matching rules
//!
//! * The sample name must match at the very beginning of the filename's
//!   basename (case-sensitive, because VCF sample names are case-sensitive).
//!   A file named `NOVEL_SAMPLE001.bam` will not match `SAMPLE001`.
//! * The character immediately after the sample name must **not** be
//!   alphanumeric. So `SAMPLE001` matches `SAMPLE001.sorted.bam` (dot after)
//!   but **not** `SAMPLE0010.bam` (`0` after).
//! * The extension (after the last dot) must be `.bam` or `.cram`,
//!   case-insensitively. Index files like `.bam.bai` are ignored.
//! * If several files match, the shortest filename wins, with an
//!   alphabetical tiebreak. So `SAMPLE001.bam` beats
//!   `SAMPLE001.sorted.markdup.bam`, and the choice is deterministic
//!   across runs and machines.

use std::path::{Path, PathBuf};

/// Find the alignment file for `sample_name` inside `dir`.
///
/// Returns `None` if the directory is unreadable or no file matches.
pub fn find_sample_bam(dir: &Path, sample_name: &str) -> Option<PathBuf> {
    let mut candidates: Vec<PathBuf> = Vec::new();

    let entries = std::fs::read_dir(dir).ok()?;
    for entry in entries.flatten() {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if file_name_matches_sample(name, sample_name) {
            candidates.push(path);
        }
    }

    // Prefer shorter names (so `SAMPLE001.bam` beats
    // `SAMPLE001.sorted.markdup.bam`), then alphabetical for determinism.
    candidates.sort_by(|a, b| {
        let la = a
            .file_name()
            .and_then(|n| n.to_str())
            .map(str::len)
            .unwrap_or(0);
        let lb = b
            .file_name()
            .and_then(|n| n.to_str())
            .map(str::len)
            .unwrap_or(0);
        la.cmp(&lb).then_with(|| a.cmp(b))
    });

    candidates.into_iter().next()
}

/// True if `file_name` names the alignment file of `sample_name`.
pub fn file_name_matches_sample(file_name: &str, sample_name: &str) -> bool {
    // Extension filter first (cheapest rejection).
    let lower = file_name.to_ascii_lowercase();
    if !(lower.ends_with(".bam") || lower.ends_with(".cram")) {
        return false;
    }

    // The sample name must appear at the very beginning. `NOVEL_SAMPLE001.bam`
    // does not identify `SAMPLE001`.
    let Some(rest) = file_name.strip_prefix(sample_name) else {
        return false;
    };

    // The character immediately after the sample name must not be
    // alphanumeric, so `SAMPLE001` does not match `SAMPLE0010.bam`.
    // Any separator (`.`, `_`, `-`) is accepted, as is the file extension
    // itself.
    match rest.chars().next() {
        None => true,                                     // exact match with no suffix — rejected by extension check above
        Some(c) if c.is_ascii_alphanumeric() => false,    // SAMPLE0010.bam, SAMPLE001a.bam
        Some(_) => true,                                  // SAMPLE001.bam, SAMPLE001_x.bam
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_exact_name() {
        assert!(file_name_matches_sample("SAMPLE001.bam", "SAMPLE001"));
        assert!(file_name_matches_sample("SAMPLE001.cram", "SAMPLE001"));
    }

    #[test]
    fn matches_common_suffixes() {
        assert!(file_name_matches_sample("SAMPLE001.sorted.bam", "SAMPLE001"));
        assert!(file_name_matches_sample("SAMPLE001.markdup.bam", "SAMPLE001"));
        assert!(file_name_matches_sample(
            "SAMPLE001.sorted.markdup.bam",
            "SAMPLE001"
        ));
        assert!(file_name_matches_sample("SAMPLE001_bwa.bam", "SAMPLE001"));
        assert!(file_name_matches_sample("SAMPLE001-aln.bam", "SAMPLE001"));
    }

    #[test]
    fn rejects_different_sample() {
        // The prefix rule must not confuse SAMPLE001 with SAMPLE0010.
        assert!(!file_name_matches_sample("SAMPLE0010.bam", "SAMPLE001"));
        assert!(!file_name_matches_sample("SAMPLE001_a.bam", "SAMPLE001a"));
        // And it must not match a file that only *contains* the name.
        assert!(!file_name_matches_sample(
            "NOVEL_SAMPLE001.bam",
            "SAMPLE001"
        ));
        assert!(!file_name_matches_sample(
            "something_SAMPLE001.sorted.bam",
            "SAMPLE001"
        ));
    }

    #[test]
    fn rejects_non_alignment_extensions() {
        assert!(!file_name_matches_sample("SAMPLE001.vcf.gz", "SAMPLE001"));
        assert!(!file_name_matches_sample("SAMPLE001.bam.bai", "SAMPLE001"));
        assert!(!file_name_matches_sample("SAMPLE001.crai", "SAMPLE001"));
        assert!(!file_name_matches_sample("SAMPLE001.sorted", "SAMPLE001"));
    }

    #[test]
    fn extension_is_case_insensitive() {
        assert!(file_name_matches_sample("SAMPLE001.BAM", "SAMPLE001"));
        assert!(file_name_matches_sample("SAMPLE001.sorted.BaM", "SAMPLE001"));
    }

    #[test]
    fn sample_name_is_case_sensitive() {
        // VCF sample names are case-sensitive; so is the resolver.
        assert!(!file_name_matches_sample("sample001.bam", "SAMPLE001"));
        assert!(file_name_matches_sample("sample001.bam", "sample001"));
    }

    #[test]
    fn finds_file_in_directory() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("SAMPLE002.sorted.bam"), "").unwrap();
        std::fs::write(dir.path().join("SAMPLE001.sorted.markdup.bam"), "").unwrap();
        std::fs::write(dir.path().join("SAMPLE0010.bam"), "").unwrap();

        let found = find_sample_bam(dir.path(), "SAMPLE001").unwrap();
        assert_eq!(
            found.file_name().unwrap(),
            "SAMPLE001.sorted.markdup.bam"
        );
    }

    #[test]
    fn prefers_shorter_name_when_multiple_match() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("SAMPLE001.sorted.markdup.bam"), "").unwrap();
        std::fs::write(dir.path().join("SAMPLE001.bam"), "").unwrap();
        std::fs::write(dir.path().join("SAMPLE001.sorted.bam"), "").unwrap();

        let found = find_sample_bam(dir.path(), "SAMPLE001").unwrap();
        assert_eq!(found.file_name().unwrap(), "SAMPLE001.bam");
    }

    #[test]
    fn returns_none_when_no_match() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("SAMPLE002.bam"), "").unwrap();
        assert!(find_sample_bam(dir.path(), "SAMPLE001").is_none());
    }

    #[test]
    fn returns_none_when_dir_missing() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("nope");
        assert!(find_sample_bam(&missing, "SAMPLE001").is_none());
    }

    #[test]
    fn deterministic_tiebreak_when_same_length() {
        // Both names are same length; alphabetical order decides.
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("SAMPLE001.zzz.bam"), "").unwrap();
        std::fs::write(dir.path().join("SAMPLE001.aaa.bam"), "").unwrap();
        let found = find_sample_bam(dir.path(), "SAMPLE001").unwrap();
        assert_eq!(found.file_name().unwrap(), "SAMPLE001.aaa.bam");
    }
}