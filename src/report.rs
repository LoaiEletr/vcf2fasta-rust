//! Per‑contig reporting and warning aggregation.
//!
//! Warnings are counted twice:
//!
//! * `warning_count` — total, unbounded.
//! * `warnings_by_reason` — a per-reason histogram, also unbounded.
//!
//! Only the individual message strings are capped, at
//! `MAX_WARNINGS_PER_CONTIG` per `ContigReport`. This keeps the tile's
//! memory footprint bounded while letting the caller report accurate
//! by-reason totals.

use std::collections::BTreeMap;

/// Maximum number of distinct warning messages stored per contig.
///
/// To avoid excessive memory usage, we only keep the first
/// `MAX_WARNINGS_PER_CONTIG` warning strings. Additional warnings are still
/// counted in `warning_count` and in `warnings_by_reason`, just not stored
/// as text.
pub const MAX_WARNINGS_PER_CONTIG: usize = 100;

/// Report for a single contig, containing metrics and a limited list of
/// warnings.
#[derive(Debug)]
pub struct ContigReport {
    pub contig: String,
    pub seen: usize,
    pub applied: usize,
    pub skipped: usize,
    pub output_files: usize,
    pub warning_count: usize,
    pub warnings: Vec<String>,
    /// Histogram of warning reasons, keyed by a short, closed-set tag.
    /// Unlike `warnings`, this is **not** capped: every call to `warn()`
    /// increments exactly one entry.
    pub warnings_by_reason: BTreeMap<&'static str, usize>,
    pub gpu_status: Option<String>,
}

impl ContigReport {
    pub fn new(contig: String) -> Self {
        Self {
            contig,
            seen: 0,
            applied: 0,
            skipped: 0,
            output_files: 0,
            warning_count: 0,
            warnings: Vec::new(),
            warnings_by_reason: BTreeMap::new(),
            gpu_status: None,
        }
    }

    /// Records a warning message.
    ///
    /// Increments `warning_count` and the appropriate
    /// `warnings_by_reason` entry unconditionally. The message itself is
    /// stored only if fewer than `MAX_WARNINGS_PER_CONTIG` messages have
    /// already been recorded.
    pub fn warn(&mut self, message: String) {
        self.warning_count += 1;
        let reason = classify_warning(&message);
        *self.warnings_by_reason.entry(reason).or_insert(0) += 1;
        if self.warnings.len() < MAX_WARNINGS_PER_CONTIG {
            self.warnings.push(message);
        }
    }
}

/// Reduce a warning message to a short, closed-set reason tag.
///
/// The order of the checks matters: a message like
/// `"REF mismatch: ... outside contig"` should be classified as
/// `ref_mismatch`, not `pos_out_of_range`, so the most specific patterns
/// come first.
///
/// Callers rely on the *stability* of these tags across runs — the
/// benchmark harness diffs them. Do not rename without updating the
/// validator.
pub fn classify_warning(w: &str) -> &'static str {
    if w.contains("overlap/out-of-order") {
        "overlap"
    } else if w.contains("REF mismatch") {
        "ref_mismatch"
    } else if w.contains("exceeds contig length") {
        "ref_interval_overflow"
    } else if w.contains("outside contig") {
        "pos_out_of_range"
    } else if w.contains("POS must be >= 1") {
        "invalid_pos"
    } else if w.contains("malformed record") || w.contains("malformed:") {
        "malformed_record"
    } else if w.contains("malformed GT") {
        "malformed_gt"
    } else if w.contains("cannot extract GT") {
        "cannot_extract_gt"
    } else if w.contains("normalisation failed") {
        "normalisation_failed"
    } else if w.contains("invalid allele") {
        "invalid_allele_index"
    } else if w.contains("negative allele") {
        "negative_allele"
    } else if w.contains("missing allele") {
        "missing_allele"
    } else if w.contains("missing from record") {
        "sample_missing"
    } else if w.contains("record belongs to") {
        "wrong_contig"
    } else if w.contains("ALT contains empty") {
        "empty_alt"
    } else if w.contains("no valid REF") || w.contains("REF allele is invalid") {
        "invalid_ref"
    } else {
        "other"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn warn_counts_and_classifies_every_message() {
        let mut r = ContigReport::new("chr1".into());

        // Push 250 overlap warnings and 10 ref-mismatch warnings.
        // Only the first 100 strings are stored, but every call is counted.
        for i in 0..250 {
            r.warn(format!(
                "[chr1:{}] overlap/out-of-order (start={} < prev_end={})",
                i + 1,
                i,
                i + 1
            ));
        }
        for i in 0..10 {
            r.warn(format!(
                "[chr1:{}] REF mismatch: VCF 'A' vs FASTA 'G'",
                i + 1
            ));
        }

        assert_eq!(r.warning_count, 260);
        assert_eq!(r.warnings.len(), MAX_WARNINGS_PER_CONTIG);

        // Both counts are uncapped.
        assert_eq!(*r.warnings_by_reason.get("overlap").unwrap(), 250);
        assert_eq!(*r.warnings_by_reason.get("ref_mismatch").unwrap(), 10);

        let sum: usize = r.warnings_by_reason.values().sum();
        assert_eq!(sum, r.warning_count);
    }

    #[test]
    fn classify_covers_common_cases() {
        assert_eq!(
            classify_warning("[chr1:100] overlap/out-of-order (start=99 < prev_end=100)"),
            "overlap"
        );
        assert_eq!(
            classify_warning("[chr1:100] REF mismatch: VCF 'A' vs FASTA 'G'"),
            "ref_mismatch"
        );
        assert_eq!(
            classify_warning("[chr1:100] REF interval [99, 105) exceeds contig length 100"),
            "ref_interval_overflow"
        );
        assert_eq!(
            classify_warning("[chr1:100] POS 100 is outside contig 'chr1' (length 50)"),
            "pos_out_of_range"
        );
        assert_eq!(
            classify_warning("[chr1:100] sample 'S1' malformed GT 'ABC' (bad), using 2 reference call(s)"),
            "malformed_gt"
        );
        assert_eq!(
            classify_warning("[chr1:100] sample 'S1' has invalid allele 5 in genotype '0|5', using placeholder 'N'"),
            "invalid_allele_index"
        );
        assert_eq!(classify_warning("[chr1:100] something unusual"), "other");
    }
}