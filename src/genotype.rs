//! Genotype and allele handling for VCF records.
//!
//! This module defines types and functions to decode VCF genotypes into a
//! simplified representation, handling ploidy, missing alleles, and phasing
//! information. It converts the low-level `rust_htslib` genotype structures
//! into an internal `AlleleCall` enum that the rest of the tool can use.

use anyhow::{Context, Result};
use rust_htslib::bcf::record::GenotypeAllele;

/// Represents how a single allele in a genotype should be written.
///
/// For each position in a genotype (haplotype), the call can be:
/// - A specific allele index (e.g., `0` for REF, `1` for first ALT, etc.).
/// - A missing/unknown allele (`.` or invalid) – replaced by a placeholder.
/// - A request to write the actual reference base from the FASTA file.
///
/// This type is used internally after parsing and validation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlleleCall {
    /// Write the allele at this index (0 = REF, 1 = first ALT, etc.).
    Index(usize),
    /// Missing/invalid allele; will be replaced by the `--no-call-string`.
    Missing,
    /// Write the reference base from the FASTA file at this position.
    Reference,
}

/// A fully decoded variant record with genotypes ready for sequence construction.
///
/// Contains all necessary information to apply the variant to a reference sequence:
/// - The start and end positions (half‑open interval) of the REF allele.
/// - All allele sequences (REF and ALTs) as byte vectors.
/// - For each sample, a vector of `AlleleCall` per haplotype (ploidy).
/// - The ploidy for each sample.
#[derive(Debug)]
pub struct DecodedRecord {
    /// Start position (0‑based, half‑open) of the REF allele.
    pub start: usize,
    /// End position (0‑based, exclusive) of the REF allele.
    pub end: usize,
    /// All allele sequences: index 0 is REF, followed by ALTs.
    pub alleles: Vec<Vec<u8>>,
    /// Genotype calls per sample. Outer vector = samples, inner = haplotypes.
    pub genotypes: Vec<Vec<AlleleCall>>,
    /// Ploidy (number of haplotypes) for each sample.
    pub ploidies: Vec<usize>,
}

/// Normalises a raw `GenotypeAllele` slice into a vector of optional allele indices.
///
/// Handles the following cases:
/// - Empty genotype → fill with `None` using the `default_ploidy`.
/// - Single missing allele (`.`) → haploid with one `None`.
/// - Explicit missing alleles (`.|.` or `.`) → `None` for each position.
/// - Normal alleles → `Some(index)`.
///
/// # Arguments
/// - `raw`: slice of `GenotypeAllele` from `rust_htslib`.
/// - `default_ploidy`: ploidy to use when the genotype has no alleles at all.
///
/// # Returns
/// A `Vec<Option<usize>>` where `None` means missing, `Some(i)` means allele index `i`.
pub fn normalize_genotype(raw: &[GenotypeAllele], default_ploidy: usize) -> Result<Vec<Option<usize>>> {
    // No alleles at all: fallback to default ploidy.
    if raw.is_empty() {
        return Ok(vec![None; default_ploidy]);
    }

    // Single '.' means haploid (one missing allele), not default ploidy.
    if raw.len() == 1 && raw[0].index().is_none() {
        return Ok(vec![None; 1]);
    }

    // Otherwise, iterate over the raw alleles.
    raw.iter()
        .map(|allele| {
            allele
                .index()
                .map(|idx| usize::try_from(idx).context("allele index does not fit usize"))
                .transpose()
        })
        .collect()
}

/// Determines whether a genotype is effectively phased.
///
/// In `rust_htslib`, the first allele of a phased genotype is stored as
/// `Unphased`, and separators before subsequent alleles are marked as `Phased`.
/// For example, `0|1` becomes `[Unphased(0), Phased(1)]`.
///
/// This function returns `true` if the genotype appears phased, with the
/// exception that haploid (length ≤ 1) and all‑missing genotypes are also
/// considered phased (since phase is irrelevant).
///
/// # Arguments
/// - `raw`: slice of `GenotypeAllele`.
///
/// # Returns
/// `true` if the genotype is considered phased, `false` otherwise.
pub fn is_effectively_phased(raw: &[GenotypeAllele]) -> bool {
    if raw.len() <= 1 || raw.iter().all(|a| a.index().is_none()) {
        return true;
    }

    raw.iter()
        .skip(1)
        .all(|a| matches!(*a, GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing))
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Tests for normalize_genotype ---

    #[test]
    fn test_normalize_empty_raw_uses_default_ploidy() {
        // Empty genotype → use default ploidy (2 → two missing alleles).
        let raw: [GenotypeAllele; 0] = [];
        let result = normalize_genotype(&raw, 2).unwrap();
        assert_eq!(result, vec![None, None]);
    }

    #[test]
    fn test_normalize_single_missing_is_haploid() {
        // Single '.' → haploid missing.
        let raw = [GenotypeAllele::UnphasedMissing];
        let result = normalize_genotype(&raw, 3).unwrap();
        assert_eq!(result, vec![None]);
    }

    #[test]
    fn test_normalize_explicit_missing_diploid() {
        // '.|.' → diploid missing.
        let raw = [GenotypeAllele::UnphasedMissing, GenotypeAllele::PhasedMissing];
        let result = normalize_genotype(&raw, 4).unwrap();
        assert_eq!(result, vec![None, None]);
    }

    #[test]
    fn test_normalize_valid_diploid_and_haploid() {
        // 0|1 → [Some(0), Some(1)]
        let raw_diploid = [GenotypeAllele::Unphased(0), GenotypeAllele::Phased(1)];
        let result_diploid = normalize_genotype(&raw_diploid, 2).unwrap();
        assert_eq!(result_diploid, vec![Some(0), Some(1)]);

        // 2 → haploid index 2
        let raw_haploid = [GenotypeAllele::Unphased(2)];
        let result_haploid = normalize_genotype(&raw_haploid, 2).unwrap();
        assert_eq!(result_haploid, vec![Some(2)]);
    }

    #[test]
    fn test_normalize_mixed_missing_and_present() {
        // 0|. → [Some(0), None]
        let raw = [GenotypeAllele::Unphased(0), GenotypeAllele::PhasedMissing];
        let result = normalize_genotype(&raw, 3).unwrap();
        assert_eq!(result, vec![Some(0), None]);
    }

    // --- Tests for phased detection ---

    #[test]
    fn phased_diploid_is_detected() {
        // 0|1 → effectively phased.
        let gt = [GenotypeAllele::Unphased(0), GenotypeAllele::Phased(1)];
        assert!(is_effectively_phased(&gt));
    }

    #[test]
    fn unphased_diploid_is_rejected() {
        // 0/1 → not phased.
        let gt = [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)];
        assert!(!is_effectively_phased(&gt));
    }

    #[test]
    fn haploid_does_not_require_phase() {
        // Haploid genotype always considered phased.
        let gt = [GenotypeAllele::Unphased(1)];
        assert!(is_effectively_phased(&gt));
    }
}