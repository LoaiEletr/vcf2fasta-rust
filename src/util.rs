//! Utility functions for allele handling and reference sequence fetching.
//!
//! This module provides helper functions used throughout the tool:
//! - Detecting symbolic or breakend alleles that cannot be written directly.
//! - Fetching reference sequences using half‑open intervals.
//! - Validating that a sequence consists only of IUPAC nucleotide codes.

use anyhow::{bail, Context, Result};
use rust_htslib::faidx;

/// Determines whether an allele is symbolic (e.g., `<DEL>`) or a breakend.
///
/// Such alleles cannot be directly written as nucleotide sequences.
/// Returns `true` if the allele:
/// - is empty,
/// - is the `*` allele (deletion),
/// - contains any of the characters `<`, `>`, `[`, or `]`.
///
/// # Arguments
/// - `allele`: the allele bytes.
///
/// # Returns
/// `true` if the allele is symbolic or a breakend, `false` otherwise.
pub fn is_symbolic_or_breakend(allele: &[u8]) -> bool {
    allele.is_empty()
        || allele == b"*"
        || allele.iter().any(|b| matches!(b, b'<' | b'>' | b'[' | b']'))
}

/// Fetches a reference sequence interval using half‑open coordinates [start, end).
///
/// `rust_htslib` uses inclusive end coordinates, so this function converts
/// the half‑open end to an inclusive coordinate (`end - 1`). If the interval
/// is empty (`start == end`), returns an empty `Vec<u8>`.
///
/// # Arguments
/// - `reference`: the indexed reference reader.
/// - `contig`: the contig name.
/// - `start`: 0‑based start position (inclusive).
/// - `end_exclusive`: 0‑based end position (exclusive).
///
/// # Errors
/// Returns an error if:
/// - `start > end_exclusive` (invalid interval).
/// - The reference fetch fails (e.g., contig missing or out of range).
///
/// # Returns
/// A `Vec<u8>` containing the sequence bytes.
pub fn fetch_half_open(
    reference: &faidx::Reader,
    contig: &str,
    start: usize,
    end_exclusive: usize,
) -> Result<Vec<u8>> {
    if start > end_exclusive {
        bail!("invalid reference interval [{start}, {end_exclusive})");
    }
    if start == end_exclusive {
        return Ok(Vec::new());
    }

    // rust-htslib faidx fetch_seq uses an inclusive end coordinate; the rest of this
    // program uses the much safer half-open convention [start, end).
    reference
        .fetch_seq(contig, start, end_exclusive - 1)
        .with_context(|| format!("could not fetch reference {contig}:[{start},{end_exclusive})"))
}

/// Checks if a byte sequence consists solely of valid IUPAC nucleotide characters.
///
/// Valid characters (case‑insensitive): `A`, `C`, `G`, `T`, `N`.
/// Returns `false` for empty sequences or any other byte.
///
/// # Arguments
/// - `seq`: the sequence bytes.
///
/// # Returns
/// `true` if the sequence is non‑empty and all characters are valid nucleotides.
pub fn is_valid_nucleotide_sequence(seq: &[u8]) -> bool {
    if seq.is_empty() {
        return false;
    }
    seq.iter().all(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_symbolic_or_breakend_accepts_normal_alleles() {
        // Normal nucleotide alleles should not be flagged.
        assert!(!is_symbolic_or_breakend(b"A"));
        assert!(!is_symbolic_or_breakend(b"ACGT"));
        assert!(!is_symbolic_or_breakend(b"G"));
        assert!(!is_symbolic_or_breakend(b"T"));
        assert!(!is_symbolic_or_breakend(b"N")); // N is a valid IUPAC code
    }

    #[test]
    fn is_symbolic_or_breakend_rejects_empty_allele() {
        // Empty allele is treated as symbolic (cannot be written).
        assert!(is_symbolic_or_breakend(b""));
    }

    #[test]
    fn is_symbolic_or_breakend_rejects_star_allele() {
        // '*' allele (deletion) is symbolic.
        assert!(is_symbolic_or_breakend(b"*"));
    }

    #[test]
    fn is_symbolic_or_breakend_rejects_symbolic_alleles() {
        // Standard symbolic alleles (e.g., <DEL>) are rejected.
        assert!(is_symbolic_or_breakend(b"<DEL>"));
        assert!(is_symbolic_or_breakend(b"<DUP>"));
        assert!(is_symbolic_or_breakend(b"<CNV>"));
        assert!(is_symbolic_or_breakend(b"<NON_REF>"));
        // Even single special characters are enough.
        assert!(is_symbolic_or_breakend(b"<"));
        assert!(is_symbolic_or_breakend(b">"));
    }

    #[test]
    fn is_symbolic_or_breakend_rejects_breakend_alleles() {
        // Breakend alleles (with brackets) are rejected.
        assert!(is_symbolic_or_breakend(b"A]chr2:123]"));
        assert!(is_symbolic_or_breakend(b"]chr1:456]A"));
        assert!(is_symbolic_or_breakend(b"[chr3:789[A"));
        assert!(is_symbolic_or_breakend(b"A[chr4:101["));
        // Even single brackets are enough.
        assert!(is_symbolic_or_breakend(b"["));
        assert!(is_symbolic_or_breakend(b"]"));
    }
}