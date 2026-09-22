//! Locate Beagle resource files (reference panel, genetic map) for a
//! chromosome by scanning a directory.
//!
//! Beagle takes optional `ref=<file>` and `map=<file>` arguments. Both are
//! per-chromosome. Rather than force the user to specify them per run, we
//! ask for the directories once and look up the right file for each
//! chromosome as it is phased.
//!
//! ## Matching rules
//!
//! * Case-insensitive.
//! * The chromosome name must appear as a **whole token**: bounded by
//!   non-alphanumeric characters (or the start/end of the file name). Both
//!   `_`, `-`, and `.` count as delimiters, because in real filenames they
//!   separate fields (`reference_chr1_panel.vcf.gz`,
//!   `1000G.chr1.ref.vcf.gz`, `chr1.map`). This prevents `chr1` from
//!   matching `chr10`, `chr100`, or `chr1abc`.
//! * Three spellings are accepted:
//!     * `chr1` (canonical),
//!     * `1` (bare),
//!     * `chrchr1` (double-prefix artifact — some pipelines prepend `chr`
//!       to a name that already carries it, producing
//!       `plink.chrchr1.GRCh38.map`).
//! * If multiple files match, the alphabetically first is returned — so
//!   the choice is deterministic across runs and machines.
//! * No match is not an error. The caller simply runs Beagle without that
//!   argument.

use std::path::{Path, PathBuf};

/// Scan `dir` for a file whose name identifies `chromosome`.
///
/// Returns `None` if the directory is unreadable or nothing matches.
pub fn find_chromosome_file(dir: &Path, chromosome: &str) -> Option<PathBuf> {
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
        if file_name_matches_chromosome(name, chromosome) {
            candidates.push(path);
        }
    }

    candidates.sort();
    candidates.into_iter().next()
}

/// True if `file_name` contains `chromosome` as a delimited token.
pub fn file_name_matches_chromosome(file_name: &str, chromosome: &str) -> bool {
    let haystack = file_name.to_lowercase();
    chromosome_forms(chromosome)
        .iter()
        .any(|needle| contains_token(&haystack, needle))
}

/// The equivalent spellings of a chromosome name.
///
/// * `chr1` -> `["chr1", "1", "chrchr1"]`
/// * `1`    -> `["chr1", "1", "chrchr1"]`
/// * `chrX` -> `["chrx", "x", "chrchrx"]`
///
/// The third form handles a common artifact where a pipeline prepends
/// `chr` to a chromosome name that already carries the prefix, producing
/// files like `plink.chrchr1.GRCh38.map`.
///
/// `lower.clone()` is required because `stripped` borrows from `lower`.
fn chromosome_forms(chromosome: &str) -> Vec<String> {
    let lower = chromosome.to_lowercase();
    let stripped = lower.strip_prefix("chr").unwrap_or(&lower).to_string();

    vec![
        format!("chr{}", stripped),    // chr1
        stripped.clone(),              // 1
        format!("chrchr{}", stripped), // chrchr1
    ]
}

/// True if `needle` appears in `haystack` bounded by non-alphanumeric
/// characters on both sides.
fn contains_token(haystack: &str, needle: &str) -> bool {
    if needle.is_empty() || needle.len() > haystack.len() {
        return false;
    }
    let bytes = haystack.as_bytes();
    let needle_bytes = needle.as_bytes();
    let n = needle_bytes.len();

    let mut i = 0;
    while i + n <= bytes.len() {
        if &bytes[i..i + n] == needle_bytes {
            let before_ok = i == 0 || !is_token_char(bytes[i - 1]);
            let after_ok = i + n == bytes.len() || !is_token_char(bytes[i + n]);
            if before_ok && after_ok {
                return true;
            }
        }
        i += 1;
    }
    false
}

/// True if `b` can be part of a chromosome-name token.
///
/// **Only alphanumeric characters are token characters.** `_`, `-`, and
/// `.` are treated as delimiters, matching the way real filenames separate
/// fields (`reference_chr1_panel.vcf.gz`). Treating `_` as a token
/// character would break the common convention of prefixing the sample or
/// suffixing the source in the name.
fn is_token_char(b: u8) -> bool {
    b.is_ascii_alphanumeric()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_exact_chr_form() {
        assert!(file_name_matches_chromosome("chr1.map", "chr1"));
        assert!(file_name_matches_chromosome("chr1.map", "1"));
    }

    #[test]
    fn matches_embedded_chr_form() {
        assert!(file_name_matches_chromosome("reference_chr1_asasa", "chr1"));
        assert!(file_name_matches_chromosome("reference_chr1_asasa", "1"));
        assert!(file_name_matches_chromosome("1000G.chr1.ref.vcf.gz", "chr1"));
        assert!(file_name_matches_chromosome("genetic_map_chr1.txt", "1"));
    }

    #[test]
    fn does_not_match_different_chromosome() {
        assert!(!file_name_matches_chromosome("chr10.map", "chr1"));
        assert!(!file_name_matches_chromosome("chr10.map", "1"));
        assert!(!file_name_matches_chromosome("chr11.ref", "chr1"));
        assert!(!file_name_matches_chromosome("chr100.map", "chr1"));
        assert!(!file_name_matches_chromosome("chr2.map", "chr22"));
    }

    #[test]
    fn matches_are_case_insensitive() {
        assert!(file_name_matches_chromosome("CHR1.MAP", "chr1"));
        assert!(file_name_matches_chromosome("Chr1.map", "CHR1"));
        assert!(file_name_matches_chromosome("Reference_CHR1_Asasa", "chr1"));
    }

    #[test]
    fn does_not_match_substring_of_number() {
        assert!(!file_name_matches_chromosome("file_100.txt", "chr1"));
        assert!(!file_name_matches_chromosome("file_10.txt", "chr1"));
    }

    #[test]
    fn matches_x_chromosome() {
        assert!(file_name_matches_chromosome("chrX.map", "chrX"));
        assert!(file_name_matches_chromosome("chrX.map", "X"));
        assert!(file_name_matches_chromosome("genetic_map_chrX.txt", "chrX"));
    }

    #[test]
    fn finds_matching_file_in_dir() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("chr2.map"), "x").unwrap();
        std::fs::write(dir.path().join("genetic_map_chr1.txt"), "y").unwrap();
        std::fs::write(dir.path().join("chr10.map"), "z").unwrap();

        let found = find_chromosome_file(dir.path(), "chr1").unwrap();
        assert_eq!(found.file_name().unwrap(), "genetic_map_chr1.txt");
    }

    #[test]
    fn returns_none_when_no_match() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("chr2.map"), "x").unwrap();
        assert!(find_chromosome_file(dir.path(), "chr1").is_none());
    }

    #[test]
    fn returns_none_when_dir_missing() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("nope");
        assert!(find_chromosome_file(&missing, "chr1").is_none());
    }

    #[test]
    fn deterministic_when_multiple_matches() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("zzz_chr1.map"), "z").unwrap();
        std::fs::write(dir.path().join("aaa_chr1.map"), "a").unwrap();

        let found = find_chromosome_file(dir.path(), "chr1").unwrap();
        assert_eq!(found.file_name().unwrap(), "aaa_chr1.map");
    }

    #[test]
    fn matches_double_chr_prefix_artifact() {
        // Real artifact: a pipeline prepends `chr` to a name that already
        // has it, producing `plink.chrchr1.GRCh38.map`.
        assert!(file_name_matches_chromosome(
            "plink.chrchr1.GRCh38.map",
            "chr1"
        ));
        assert!(file_name_matches_chromosome(
            "plink.chrchr1.GRCh38.map",
            "1"
        ));
        assert!(file_name_matches_chromosome("chrchr22.map", "chr22"));
    }

    #[test]
    fn does_not_confuse_double_prefix_with_neighbours() {
        // `chrchr1` must still not match `chrchr10`, `chrchr11`, `chrchr100`.
        assert!(!file_name_matches_chromosome("plink.chrchr10.map", "chr1"));
        assert!(!file_name_matches_chromosome("plink.chrchr11.map", "chr1"));
        assert!(!file_name_matches_chromosome("plink.chrchr100.map", "chr1"));
    }

    #[test]
    fn finds_double_prefix_file_in_dir() {
        // End-to-end check: the exact filename from the benchmark dataset
        // must be located by the resolver.
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("plink.chrchr1.GRCh38.map"), "x").unwrap();

        let found = find_chromosome_file(dir.path(), "chr1").unwrap();
        assert_eq!(found.file_name().unwrap(), "plink.chrchr1.GRCh38.map");
    }
}