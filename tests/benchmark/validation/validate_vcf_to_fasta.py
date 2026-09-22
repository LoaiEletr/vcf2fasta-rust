#!/usr/bin/env python3
"""
validate_vcf_to_fasta.py — exhaustive VCF-to-FASTA validator.

Builds the expected sequence per haplotype from the VCF using the tool's
own policy (skip invalid POS / missing REF / symbolic REF / overlap /
REF-mismatch; write REF for allele 0, ALT for valid ALT indices, the
no-call string for missing or invalid alleles), then compares byte-by-byte
against the tool's FASTA output.

If vcflib FASTA outputs are supplied and exist, the same expected
sequence is compared against vcflib as an independent cross-check.

Ploidy:
    Pass a single haplotype FASTA (or omit the second positional) for
    haploid chromosomes (male chrX/chrY, chrM, …). The validator will
    only build and compare hap0 in that case.

Exit code:
    0  all requested comparisons passed
    1  at least one comparison failed
    2  setup error (missing file, bad VCF, missing reference)

Usage:
    # diploid
    validate_vcf_to_fasta.py <vcf> <my_hap0> <my_hap1> \\
        --reference <reference.fa> [options]

    # haploid
    validate_vcf_to_fasta.py <vcf> <my_hap0> \\
        --reference <reference.fa> [options]
"""

import argparse
import gzip
import os
import re
import sys


# ---------------------------------------------------------------------------
# FASTA loading
# ---------------------------------------------------------------------------

def load_fasta(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing file: {path}")
    parts = []
    with open(path) as f:
        for line in f:
            if not line.startswith('>'):
                parts.append(line.strip().upper())
    return "".join(parts)


class ReferenceLookup:
    def __init__(self, path):
        self.available = False
        self._pysam = None
        self._cache = {}
        if not path or not os.path.exists(path):
            return
        try:
            import pysam
            self._pysam = pysam.FastaFile(path)
            self.available = True
            return
        except Exception:
            pass
        try:
            seqs = {}
            current = None
            buf = []
            with open(path) as f:
                for line in f:
                    if line.startswith('>'):
                        if current is not None:
                            seqs[current] = "".join(buf).upper()
                        current = line[1:].split()[0]
                        buf = []
                    else:
                        buf.append(line.strip())
                if current is not None:
                    seqs[current] = "".join(buf).upper()
            self._cache = seqs
            self.available = True
        except Exception:
            self.available = False

    def _candidates(self, chrom):
        out = [chrom]
        if chrom.startswith("chr"):
            out.append(chrom[3:])
        else:
            out.append("chr" + chrom)
        return out

    def fetch(self, chrom, pos_1based, length):
        if not self.available or length <= 0:
            return None
        for name in self._candidates(chrom):
            if self._pysam is not None:
                try:
                    return self._pysam.fetch(
                        name, pos_1based - 1, pos_1based - 1 + length
                    ).upper()
                except Exception:
                    continue
            else:
                seq = self._cache.get(name)
                if seq is not None:
                    i = pos_1based - 1
                    return seq[i:i + length]
        return None

    def fetch_all(self, chrom):
        if not self.available:
            return None
        for name in self._candidates(chrom):
            if self._pysam is not None:
                try:
                    return self._pysam.fetch(name).upper()
                except Exception:
                    continue
            else:
                seq = self._cache.get(name)
                if seq is not None:
                    return seq
        return None


# ---------------------------------------------------------------------------
# VCF helpers
# ---------------------------------------------------------------------------

def parse_gt(gt_str):
    if not gt_str:
        return True, []
    is_phased = '|' in gt_str
    parts = re.split(r'[/|]', gt_str)
    out = []
    for p in parts:
        if p == '.' or p.startswith('-'):
            out.append(None)
        else:
            try:
                out.append(int(p))
            except ValueError:
                out.append(None)
    return is_phased, out


def allele_is_valid(alt):
    if not alt or alt == '.':
        return False
    if any(c in alt for c in '<>[]'):
        return False
    return all(c in 'ACGTNacgtn' for c in alt)


def read_variants(vcf_path, sample_name=None):
    """Read variants; if sample_name is given, extract GT from that column."""
    out = []
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as f:
        sample_col = None
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                if len(header) <= 9:
                    raise ValueError("VCF has no sample columns")
                samples = header[9:]
                if sample_name is not None:
                    if sample_name not in samples:
                        raise ValueError(
                            f"Sample '{sample_name}' not found in VCF. "
                            f"Available: {samples}"
                        )
                    sample_col = 9 + samples.index(sample_name)
                else:
                    sample_col = 9
                continue
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if sample_col is None or len(parts) <= sample_col:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            ref = parts[3].upper()
            alt_field = parts[4].upper()
            alts = [] if alt_field == '.' else alt_field.split(',')
            fmt_fields = parts[8].split(':')
            if 'GT' not in fmt_fields:
                continue
            gt_idx = fmt_fields.index('GT')
            sample_fields = parts[sample_col].split(':')
            if gt_idx >= len(sample_fields):
                continue
            out.append({
                'chrom': parts[0],
                'pos': pos,
                'ref': ref,
                'alts': alts,
                'alt_field': alt_field,
                'gt': sample_fields[gt_idx],
            })
    out.sort(key=lambda v: v['pos'])
    return out


# ---------------------------------------------------------------------------
# Policy — keep in sync with src/vcf.rs::process_chunk
# ---------------------------------------------------------------------------

def should_skip(v, last_end, no_validate_ref, ref_lookup, chrom):
    pos_0based = v['pos'] - 1
    ref = v['ref']
    if pos_0based < 0:
        return True, "invalid POS"
    if not ref or ref == '.':
        return True, "missing REF"
    if any(c in ref for c in '<>[]'):
        return True, "symbolic REF"
    if pos_0based < last_end:
        return True, "overlap"
    if not no_validate_ref and ref_lookup is not None and chrom is not None:
        fasta_ref = ref_lookup.fetch(chrom, v['pos'], len(ref))
        if fasta_ref is not None and fasta_ref != ref:
            return True, "REF mismatch"
    return False, None


def allele_for_hap(v, hap_idx, no_call='N'):
    _, alleles = parse_gt(v['gt'])
    if ('|' not in v['gt']
            and len(alleles) == 2
            and alleles[0] is not None
            and alleles[1] is not None
            and alleles[0] > alleles[1]):
        alleles = [alleles[1], alleles[0]]
    a = alleles[hap_idx] if hap_idx < len(alleles) else None
    if a is None:
        return no_call[:1]
    if a == 0:
        return v['ref']
    if 1 <= a <= len(v['alts']):
        alt = v['alts'][a - 1]
        return alt if allele_is_valid(alt) else no_call[:1]
    return no_call[:1]


# ---------------------------------------------------------------------------
# Expected-sequence builder
# ---------------------------------------------------------------------------

def build_expected(variants, hap_idx, ref_seq, ref_lookup, chrom,
                   no_validate_ref, no_call='N'):
    out_parts = []
    segments = []
    out_pos = 0
    last_end = 0
    stats = {'applied': 0, 'skipped': 0}
    skipped_details = []

    def emit(chunk, kind, variant):
        nonlocal out_pos
        if not chunk:
            return
        out_parts.append(chunk)
        segments.append((out_pos, out_pos + len(chunk), kind, variant))
        out_pos += len(chunk)

    for v in variants:
        skip, reason = should_skip(
            v, last_end, no_validate_ref, ref_lookup, chrom
        )
        if skip:
            stats['skipped'] += 1
            skipped_details.append((v, reason))
            continue

        pos_0 = v['pos'] - 1
        if pos_0 > last_end:
            emit(ref_seq[last_end:pos_0], 'ref', None)

        chunk = allele_for_hap(v, hap_idx, no_call)
        emit(chunk, 'variant', v)
        last_end = pos_0 + len(v['ref'])
        stats['applied'] += 1

    if last_end < len(ref_seq):
        emit(ref_seq[last_end:], 'ref', None)

    stats['skipped_details'] = skipped_details
    return ''.join(out_parts), segments, stats


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def compare_one(expected, segments, actual, label, max_errors, start_reported,
                strict=False):
    ok = True
    total = 0
    reported = start_reported
    first_errors = []

    length_ok = (len(expected) == len(actual))
    if not length_ok:
        first_errors.append(
            f"⚠️ {label}: length mismatch — "
            f"expected {len(expected):,}, actual {len(actual):,} "
            f"(diff {len(actual) - len(expected):+,})"
        )
        ok = False
        if strict:
            return ok, total, reported, first_errors

    n = min(len(expected), len(actual))
    for seg_start, seg_end, kind, variant in segments:
        seg_end_c = min(seg_end, n)
        if seg_start >= seg_end_c:
            continue
        e_chunk = expected[seg_start:seg_end_c]
        a_chunk = actual[seg_start:seg_end_c]
        if e_chunk == a_chunk:
            continue
        for j in range(len(e_chunk)):
            if e_chunk[j] != a_chunk[j]:
                total += 1
                ok = False
                if reported < max_errors:
                    reported += 1
                    abs_pos = seg_start + j
                    if kind == 'variant' and variant is not None:
                        origin = (f"variant {variant['chrom']}:{variant['pos']:,} "
                                  f"(REF={variant['ref']}, ALT={variant['alts']}, "
                                  f"GT={variant['gt']})")
                    else:
                        origin = f"reference interval [{seg_start:,}, {seg_end:,})"
                    ws = max(0, abs_pos - 20)
                    we = min(len(actual), abs_pos + 20)
                    first_errors.append(
                        f"\n🚨 {label} mismatch at output index {abs_pos:,}\n"
                        f"   Origin: {origin}\n"
                        f"   Expected: ...{expected[ws:we]}...\n"
                        f"   Actual:   ...{actual[ws:we]}..."
                    )
                break
        if reported >= max_errors:
            break

    return ok, total, reported, first_errors


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vcf')
    ap.add_argument('my_hap0',
                    help="First haplotype FASTA (always required).")
    ap.add_argument('my_hap1', nargs='?', default=None,
                    help="Second haplotype FASTA. Omit for haploid "
                         "chromosomes (male chrX/chrY, chrM, …).")
    ap.add_argument('--reference', required=True)
    ap.add_argument('--sample', default=None,
                    help="Sample name to extract GT from. Defaults to the "
                         "first sample column in the VCF.")
    ap.add_argument('--no-validate-ref', action='store_true')
    ap.add_argument('--vcflib-hap0', default=None)
    ap.add_argument('--vcflib-hap1', default=None)
    ap.add_argument('--max-errors', type=int, default=20)
    ap.add_argument('--list-skipped', type=int, default=20,
                    help="Print up to N skipped-variant detail lines. "
                         "0 disables the list.")
    ap.add_argument('--strict', action='store_true',
                    help="Fail on length mismatch without attempting byte comparison.")
    args = ap.parse_args()

    # ---- Load tool hap0 (required) ----------------------------------------
    try:
        tool_h0 = load_fasta(args.my_hap0)
    except Exception as e:
        print(f"❌ Tool hap0 FASTA load error: {e}")
        return 2

    # ---- Load tool hap1 (optional) ----------------------------------------
    # If the path wasn't given, or the file doesn't exist, treat as haploid.
    haploid = False
    tool_h1 = None
    if not args.my_hap1 or args.my_hap1 in ('-', ''):
        haploid = True
    elif not os.path.exists(args.my_hap1):
        haploid = True
    else:
        try:
            tool_h1 = load_fasta(args.my_hap1)
        except Exception as e:
            print(f"❌ Tool hap1 FASTA load error: {e}")
            return 2

    # ---- Load vcflib FASTA (optional) -------------------------------------
    vcf_h0 = vcf_h1 = None
    vcflib_available = False
    vcflib_haploid = False

    if args.vcflib_hap0 and os.path.exists(args.vcflib_hap0):
        try:
            vcf_h0 = load_fasta(args.vcflib_hap0)
            vcflib_available = True
        except Exception:
            vcf_h0 = None
            vcflib_available = False

    if vcflib_available:
        if args.vcflib_hap1 and os.path.exists(args.vcflib_hap1):
            try:
                vcf_h1 = load_fasta(args.vcflib_hap1)
            except Exception:
                vcf_h1 = None
        if vcf_h1 is None:
            vcflib_haploid = True

    # ---- Reference ---------------------------------------------------------
    ref_lookup = ReferenceLookup(args.reference)
    if not ref_lookup.available:
        print(f"❌ Could not open reference: {args.reference}")
        return 2

    # ---- VCF ---------------------------------------------------------------
    try:
        variants = read_variants(args.vcf, sample_name=args.sample)
    except ValueError as e:
        print(f"❌ {e}")
        return 2

    if not variants:
        print("⚠️ No variants in VCF")
        return 0

    chroms_seen = sorted({v['chrom'] for v in variants})
    if len(chroms_seen) > 1:
        print(f"❌ VCF contains multiple contigs: {chroms_seen}. "
              f"This validator expects one chromosome per call.")
        return 2
    chrom = chroms_seen[0]

    ref_seq = ref_lookup.fetch_all(chrom)
    if ref_seq is None:
        print(f"❌ Could not fetch contig {chrom} from reference")
        return 2

    # ---- Info block -------------------------------------------------------
    print(f"ℹ️  Sample:     {args.sample if args.sample else '(first in VCF)'}")
    print(f"ℹ️  Chromosome: {chrom}")
    print(f"ℹ️  Ploidy:     {'1 (haploid)' if haploid else '2 (diploid)'}")
    print(f"ℹ️  Variants:   {len(variants):,}")
    print(f"ℹ️  Reference:  {len(ref_seq):,} bases")
    print(f"ℹ️  Tool hap0:  {len(tool_h0):,} bases")
    if not haploid:
        print(f"ℹ️  Tool hap1:  {len(tool_h1):,} bases")
    if vcflib_available:
        print(f"ℹ️  vcflib hap0: {len(vcf_h0):,} bases")
        if not vcflib_haploid:
            print(f"ℹ️  vcflib hap1: {len(vcf_h1):,} bases")
    else:
        print("ℹ️  vcflib not available for this chromosome")

    # ---- Build expected ---------------------------------------------------
    exp0, seg0, stats = build_expected(
        variants, 0, ref_seq, ref_lookup, chrom, args.no_validate_ref
    )
    if haploid:
        exp1, seg1 = None, None
    else:
        exp1, seg1, _ = build_expected(
            variants, 1, ref_seq, ref_lookup, chrom, args.no_validate_ref
        )
    print(f"ℹ️  Policy: {stats['applied']:,} applied, {stats['skipped']:,} skipped")

    # ---- List skipped variants (with reasons) -----------------------------
    if args.list_skipped > 0 and stats['skipped_details']:
        shown = stats['skipped_details'][:args.list_skipped]
        print(f"ℹ️  Skipped variants (showing {len(shown)} of {stats['skipped']}):")
        for v, reason in shown:
            print(
                f"     SKIPPED {v['chrom']}:{v['pos']} "
                f"REF={v['ref']} ALT={v['alt_field']} GT={v['gt']} "
                f"reason={reason}"
            )
        if stats['skipped'] > len(shown):
            print(f"     ... {stats['skipped'] - len(shown)} more not shown "
                  f"(raise --list-skipped to see them)")

    # ---- Compare tool -----------------------------------------------------
    reported = 0
    all_errors = []

    tool_ok_0, _, reported, errs = compare_one(
        exp0, seg0, tool_h0, "TOOL hap0", args.max_errors, reported, args.strict
    )
    all_errors.extend(errs)
    tool_ok_1 = True
    if not haploid:
        tool_ok_1, _, reported, errs = compare_one(
            exp1, seg1, tool_h1, "TOOL hap1", args.max_errors, reported, args.strict
        )
        all_errors.extend(errs)
    tool_ok = tool_ok_0 and tool_ok_1

    # ---- Compare vcflib (optional) ---------------------------------------
    vcf_ok = None
    if vcflib_available:
        vcf_ok_0, _, reported, errs = compare_one(
            exp0, seg0, vcf_h0, "VCFLIB hap0", args.max_errors, reported, args.strict
        )
        all_errors.extend(errs)
        vcf_ok_1 = True
        if not vcflib_haploid and not haploid:
            vcf_ok_1, _, reported, errs = compare_one(
                exp1, seg1, vcf_h1, "VCFLIB hap1", args.max_errors, reported, args.strict
            )
            all_errors.extend(errs)
        vcf_ok = vcf_ok_0 and vcf_ok_1

    # ---- Report -----------------------------------------------------------
    print()
    if all_errors:
        for e in all_errors:
            print(e)
        print()

    print("📊 Results:")
    print(f"   TOOL   vs expected: {'✅ PASS' if tool_ok else '❌ FAIL'}")
    if vcflib_available:
        print(f"   VCFLIB vs expected: {'✅ PASS' if vcf_ok else '❌ FAIL'}")
        if tool_ok and vcf_ok:
            print("   ✅ Tool and vcflib both agree with expected (cross-checked)")
        elif tool_ok and not vcf_ok:
            print("   ⚠️  Tool agrees with expected; vcflib diverges")
        elif not tool_ok and vcf_ok:
            print("   ⚠️  vcflib agrees with expected; tool diverges")
        else:
            print("   ❌ Both tool and vcflib diverge from expected")
    else:
        print("   ℹ️  vcflib not available; tool-only check")

    if tool_ok and (vcf_ok is None or vcf_ok):
        return 0
    return 1


if __name__ == '__main__':
    sys.exit(main())