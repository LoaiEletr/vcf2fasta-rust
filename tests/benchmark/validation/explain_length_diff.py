#!/usr/bin/env python3
"""
explain_length_diff.py — explain a length mismatch between the Rust tool and
vcflib by attributing it to variants the Rust tool skipped per its overlap
policy, and print the divergence in aligned bracket form showing BOTH the
first (applied) emission and the skipped (duplicate) emission:

    Rust:    ...CTCCAGCCTGGGTGGCAGAGCGAGACTCTG[T][∅]TAAATAAATAA...
    vcflib:  ...CTCCAGCCTGGGTGGCAGAGCGAGACTCTG[T][T]TAAATAAATAA...
                                              ^

Lengths are reported in BASES.

TSV output (with --tsv) has 9 tab-separated fields:
    COMPARE <hap_label> <len_a> <len_b> <delta> <nd> <first> <disp_a> <disp_b>
where hap_label is "<sample>_<chromosome>_H<hap>" when --chromosome is given,
or "<sample>_H<hap>" otherwise.

Exit codes:
    0  lengths match and content identical, OR Δ is fully explained
    1  Δ not fully explained, or content differs at same length
    2  setup error
"""

import argparse
import gzip
import os
import re
import sys


LEFT_CTX = 30
RIGHT_CTX = 30
DISPLAY_PAD = 10


# ---------------------------------------------------------------------------
# FASTA / reference loading
# ---------------------------------------------------------------------------

def load_fasta_seq(path):
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
            current, buf = None, []
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
        return [chrom, chrom[3:]] if chrom.startswith("chr") else [chrom, "chr" + chrom]

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
# VCF parsing
# ---------------------------------------------------------------------------

def parse_gt(gt_str):
    if not gt_str:
        return []
    out = []
    for p in re.split(r'[/|]', gt_str):
        if p == '.' or p.startswith('-'):
            out.append(None)
        else:
            try:
                out.append(int(p))
            except ValueError:
                out.append(None)
    return out


def allele_is_valid(alt):
    if not alt or alt == '.':
        return False
    if any(c in alt for c in '<>[]'):
        return False
    return all(c in 'ACGTNacgtn' for c in alt)


def read_variants(vcf_path, sample_name):
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
                if sample_name not in samples:
                    raise ValueError(
                        f"Sample '{sample_name}' not in VCF. Available: {samples}"
                    )
                sample_col = 9 + samples.index(sample_name)
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
            fmt = parts[8].split(':')
            if 'GT' not in fmt:
                continue
            gt_idx = fmt.index('GT')
            sf = parts[sample_col].split(':')
            if gt_idx >= len(sf):
                continue
            out.append({
                'chrom': parts[0], 'pos': pos, 'ref': ref,
                'alts': alts, 'alt_field': alt_field, 'gt': sf[gt_idx],
            })
    out.sort(key=lambda v: v['pos'])
    return out


def allele_for_hap(v, hap_idx):
    alleles = parse_gt(v['gt'])
    if ('|' not in v['gt'] and len(alleles) == 2
            and alleles[0] is not None and alleles[1] is not None
            and alleles[0] > alleles[1]):
        alleles = [alleles[1], alleles[0]]
    if hap_idx < len(alleles):
        return alleles[hap_idx]
    return None


def emission_bytes(v, hap_idx, no_call='N'):
    a = allele_for_hap(v, hap_idx)
    if a is None:
        return no_call
    if a == 0:
        return v['ref']
    if 1 <= a <= len(v['alts']):
        alt = v['alts'][a - 1]
        return alt if allele_is_valid(alt) else no_call
    return no_call


# ---------------------------------------------------------------------------
# Rust overlap policy
# ---------------------------------------------------------------------------

def find_skipped_variants(variants):
    skipped = []
    last_end = 0
    applied_spans = []
    for i, v in enumerate(variants):
        pos_0 = v['pos'] - 1
        if pos_0 < last_end:
            orig_idx = -1
            orig_v = None
            for (s, e, j, ov) in reversed(applied_spans):
                if s <= pos_0 < e:
                    orig_idx = j
                    orig_v = ov
                    break
            skipped.append({
                'skip_idx': i,
                'skip_v': v,
                'orig_idx': orig_idx,
                'orig_v': orig_v,
            })
            continue
        last_end = pos_0 + len(v['ref'])
        applied_spans.append((pos_0, last_end, i, v))
    return skipped


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def fmt_num(n):
    return f"{n:,}"


def fmt_signed(n):
    if n == 0:
        return "0"
    if n > 0:
        return f"+{n:,}"
    return f"-{abs(n):,}"


def plural(n, singular, plural_form=None):
    if n == 1:
        return singular
    return plural_form if plural_form else singular + 's'


def first_divergent_prefix(a, b):
    n = min(len(a), len(b))
    lo, hi = 0, 1
    while hi <= n and a[:hi] == b[:hi]:
        lo = hi
        hi *= 2
    if lo == n and (len(a) == len(b) or a[:n] == b[:n]):
        return n
    hi = min(hi, n)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if a[:mid] == b[:mid]:
            lo = mid
        else:
            hi = mid
    return lo


def find_anchor(seq, ref_seq, ref_pos):
    if ref_pos < 0 or ref_pos >= len(ref_seq):
        return -1
    avail = len(ref_seq) - ref_pos
    for a_len in (128, 96, 64, 48, 32, 24, 16, 12, 8):
        if a_len > avail:
            continue
        anchor = ref_seq[ref_pos:ref_pos + a_len]
        if not anchor:
            continue
        hit = seq.find(anchor)
        if hit >= 0:
            return hit
    return -1


def render_emission(emitted):
    return emitted if emitted else '∅'


def make_display_line_single(label, left_ctx, emission, right_ctx, pad=DISPLAY_PAD):
    prefix = f"{label + ':':<{pad}}"
    return f"{prefix}...{left_ctx}[{render_emission(emission)}]{right_ctx}..."


def make_display_line_dual(label, left_ctx, first_emission, skip_emission,
                           right_ctx, pad=DISPLAY_PAD):
    prefix = f"{label + ':':<{pad}}"
    return (f"{prefix}...{left_ctx}"
            f"[{render_emission(first_emission)}]"
            f"[{render_emission(skip_emission)}]"
            f"{right_ctx}...")


def make_display_pair_single(label_a, label_b, left_ctx, emitted_a, emitted_b,
                             right_ctx):
    disp_a = make_display_line_single(label_a, left_ctx, emitted_a, right_ctx)
    disp_b = make_display_line_single(label_b, left_ctx, emitted_b, right_ctx)
    return disp_a, disp_b


def make_display_pair_dual(label_a, label_b, left_ctx,
                           first_a, skip_a,
                           first_b, skip_b,
                           right_ctx):
    disp_a = make_display_line_dual(label_a, left_ctx, first_a, skip_a, right_ctx)
    disp_b = make_display_line_dual(label_b, left_ctx, first_b, skip_b, right_ctx)
    return disp_a, disp_b


def caret_column_single(disp_line):
    return disp_line.index('[') if '[' in disp_line else -1


def caret_column_dual(disp_line):
    first = disp_line.find('[')
    if first < 0:
        return -1
    second = disp_line.find('[', first + 1)
    return second if second >= 0 else first


def emit_tsv(hap_label, la, lb, delta, nd, first, disp_a, disp_b):
    delta_str = fmt_signed(delta)
    print(f"COMPARE\t{hap_label}\t{la}\t{lb}\t{delta_str}"
          f"\t{nd}\t{first}\t{disp_a or '-'}\t{disp_b or '-'}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vcf')
    ap.add_argument('rust_fa')
    ap.add_argument('vcflib_fa')
    ap.add_argument('--reference', required=True)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--haplotype', required=True, type=int, choices=[0, 1])
    ap.add_argument('--chromosome', default=None,
                    help="Chromosome name to include in the haplotype label "
                         "(e.g. chr12). If omitted, the label is "
                         "'<sample>_H<hap>'.")
    ap.add_argument('--label-a', default='Rust')
    ap.add_argument('--label-b', default='vcflib')
    ap.add_argument('--max-diffs', type=int, default=20,
                    help="Max skipped variants to print in detail (0 = all).")
    ap.add_argument('--tsv', action='store_true')
    args = ap.parse_args()

    if not os.path.exists(args.rust_fa):
        print(f"❌ Missing {args.label_a} FASTA: {args.rust_fa}"); return 2
    if not os.path.exists(args.vcflib_fa):
        print(f"❌ Missing {args.label_b} FASTA: {args.vcflib_fa}"); return 2

    seq_a = load_fasta_seq(args.rust_fa)
    seq_b = load_fasta_seq(args.vcflib_fa)
    len_a, len_b = len(seq_a), len(seq_b)
    delta = len_b - len_a

    # Build the label. Include the chromosome when provided so results for
    # different chromosomes don't collide in the summary table.
    if args.chromosome:
        hap_label = f"{args.sample}_{args.chromosome}_H{args.haplotype}"
    else:
        hap_label = f"{args.sample}_H{args.haplotype}"

    # ---------------- header ----------------
    print("=" * 78)
    print(f"Length-diff explainer: {hap_label}")
    print(f"  {args.label_a:<10} length: {fmt_num(len_a)} bases")
    print(f"  {args.label_b:<10} length: {fmt_num(len_b)} bases")
    print(f"  Δ ({args.label_b} − {args.label_a}):    "
          f"{fmt_signed(delta)} bases")

    # ---------------- same length ----------------
    if delta == 0:
        if seq_a == seq_b:
            print()
            print("✅ Sequences are identical (every base matches, "
                  "position for position).")
            if args.tsv:
                emit_tsv(hap_label, len_a, len_b, 0, "0", "-", "-", "-")
            return 0

        p = first_divergent_prefix(seq_a, seq_b)
        print()
        print(f"⚠️  Lengths match, but content differs at output position "
              f"{fmt_num(p)}:")

        left = seq_a[max(0, p - LEFT_CTX):p]
        right_a = seq_a[p + 1:p + 1 + RIGHT_CTX]
        right_b = seq_b[p + 1:p + 1 + RIGHT_CTX]
        disp_a, disp_b = make_display_pair_single(
            args.label_a, args.label_b,
            left, seq_a[p:p + 1], seq_b[p:p + 1], right_a,
        )
        print()
        print(f"  {disp_a}")
        print(f"  {disp_b}")
        cc = caret_column_single(disp_a)
        if cc >= 0:
            print(' ' * (2 + cc) + '^')
        if right_a != right_b:
            print("       (right context also differs)")

        if args.tsv:
            emit_tsv(hap_label, len_a, len_b, 0, "CONTENT_DIFF",
                     str(p), disp_a, disp_b)
        return 1

    # ---------------- lengths differ ----------------
    ref = ReferenceLookup(args.reference)
    if not ref.available:
        print(f"❌ Could not open reference: {args.reference}"); return 2

    try:
        variants = read_variants(args.vcf, sample_name=args.sample)
    except ValueError as e:
        print(f"❌ {e}"); return 2

    if not variants:
        print("⚠️ No variants in VCF"); return 1

    chroms = sorted({v['chrom'] for v in variants})
    if len(chroms) > 1:
        print(f"❌ Multi-contig VCF: {chroms}"); return 2
    chrom = chroms[0]
    ref_seq = ref.fetch_all(chrom)
    if ref_seq is None:
        print(f"❌ Could not fetch contig {chrom}"); return 2

    skipped = find_skipped_variants(variants)

    print()
    print(f"Variants in VCF:                                "
          f"{fmt_num(len(variants))}")
    print(f"Skipped by {args.label_a} overlap policy:              "
          f"{fmt_num(len(skipped))}")

    if not skipped:
        print()
        print(f"⚠️  No overlap-skips found. The {fmt_signed(delta)}-base "
              f"difference is NOT explained by the overlap policy.")
        p = first_divergent_prefix(seq_a, seq_b)
        left = seq_a[max(0, p - LEFT_CTX):p]
        right_a = seq_a[p + 1:p + 1 + RIGHT_CTX]
        disp_a, disp_b = make_display_pair_single(
            args.label_a, args.label_b,
            left, seq_a[p:p + 1], seq_b[p:p + 1], right_a,
        )
        print()
        print(f"🔎 First base divergence at output position {fmt_num(p)}:")
        print()
        print(f"  {disp_a}")
        print(f"  {disp_b}")
        cc = caret_column_single(disp_a)
        if cc >= 0:
            print(' ' * (2 + cc) + '^')
        if args.tsv:
            emit_tsv(hap_label, len_a, len_b, delta, "LENGTH_MISMATCH",
                     str(p), disp_a, disp_b)
        return 1

    limit = len(skipped) if args.max_diffs == 0 else min(args.max_diffs,
                                                          len(skipped))
    print()
    print(f"── Skipped variants (each falls inside a previous variant's "
          f"REF span) ──")

    predicted_total = 0
    confirmed_count = 0
    first_confirmed_pos = None
    first_disp_a = first_disp_b = '-'

    for idx_in_list, item in enumerate(skipped):
        v = item['skip_v']
        orig_v = item['orig_v']

        emitted_b = emission_bytes(v, args.haplotype)
        pred_extra = len(emitted_b)
        predicted_total += pred_extra

        first_emitted = emission_bytes(orig_v, args.haplotype) if orig_v else ''

        if orig_v is not None:
            display_left_ref = orig_v['pos'] - 1
            anchor_start_ref = display_left_ref + len(orig_v['ref'])
        else:
            display_left_ref = v['pos'] - 1
            anchor_start_ref = display_left_ref + len(v['ref'])

        hit_a = find_anchor(seq_a, ref_seq, anchor_start_ref)
        hit_b = find_anchor(seq_b, ref_seq, anchor_start_ref)

        left_lo = max(0, display_left_ref - LEFT_CTX)
        left_ctx = ref_seq[left_lo:display_left_ref]
        right_hi = min(len(ref_seq), anchor_start_ref + RIGHT_CTX)
        right_ctx = ref_seq[anchor_start_ref:right_hi]

        show = idx_in_list < limit

        disp_a = disp_b = None
        if show:
            print()
            print(f"   [{idx_in_list + 1}] {v['chrom']}:{fmt_num(v['pos'])}  "
                  f"REF={v['ref']}  ALT={v['alt_field']}  GT={v['gt']}")
            if orig_v is not None:
                print(f"        applies first at  "
                      f"{orig_v['chrom']}:{fmt_num(orig_v['pos'])}  "
                      f"REF={orig_v['ref']}  ALT={orig_v['alt_field']}  "
                      f"GT={orig_v['gt']}")
            print(f"       {args.label_a:<10} "
                  f"first={render_emission(first_emitted)}, "
                  f"duplicate={render_emission('')} "
                  f"(skipped)")
            print(f"       {args.label_b:<10} "
                  f"first={render_emission(first_emitted)}, "
                  f"duplicate={render_emission(emitted_b)} "
                  f"(applied again)")

            disp_a, disp_b = make_display_pair_dual(
                args.label_a, args.label_b,
                left_ctx,
                first_emitted, '',
                first_emitted, emitted_b,
                right_ctx,
            )
            print()
            print(f"  {disp_a}")
            print(f"  {disp_b}")
            cc = caret_column_dual(disp_a)
            if cc >= 0:
                print(' ' * (2 + cc) + '^')

        if hit_a >= 0 and hit_b >= 0:
            observed_extra = hit_b - hit_a
            if observed_extra == pred_extra:
                confirmed_count += 1
                if first_confirmed_pos is None:
                    first_confirmed_pos = v['pos']
                    if disp_a is None:
                        disp_a, disp_b = make_display_pair_dual(
                            args.label_a, args.label_b,
                            left_ctx,
                            first_emitted, '',
                            first_emitted, emitted_b,
                            right_ctx,
                        )
                    first_disp_a = disp_a
                    first_disp_b = disp_b
                if show:
                    print()
                    print(f"       ✅ confirmed: {args.label_b} anchor at "
                          f"output position {fmt_num(hit_b)} vs "
                          f"{args.label_a} at {fmt_num(hit_a)} = "
                          f"{observed_extra:+,} "
                          f"{plural(observed_extra, 'base')}")
            else:
                if show:
                    print()
                    print(f"       ⚠️  output check: {args.label_a} anchor at "
                          f"{fmt_num(hit_a)}, {args.label_b} anchor at "
                          f"{fmt_num(hit_b)}, offset = "
                          f"{observed_extra:+,} "
                          f"(expected {pred_extra:+,})")
        else:
            if show:
                print()
                if hit_a < 0:
                    print(f"       ℹ️  could not locate post-locus anchor in "
                          f"{args.label_a}")
                if hit_b < 0:
                    print(f"       ℹ️  could not locate post-locus anchor in "
                          f"{args.label_b}")

        if first_disp_a == '-' and show and disp_a is not None:
            first_disp_a = disp_a
            first_disp_b = disp_b

    if limit < len(skipped):
        print()
        print(f"   ... {fmt_num(len(skipped) - limit)} more skipped variants "
              f"not shown (raise --max-diffs to see them)")

    print()
    print("── Summary ──")
    print(f"   Predicted Δ from skipped variants:  "
          f"{fmt_signed(predicted_total)} bases")
    print(f"   Actual Δ:                            "
          f"{fmt_signed(delta)} bases")

    explained = (predicted_total == delta)
    if explained:
        print(f"   ✅ The base-count difference is FULLY EXPLAINED by "
              f"{args.label_b} applying")
        print(f"      the {fmt_num(len(skipped))} variant(s) that "
              f"{args.label_a} skips per its overlap policy.")
        if confirmed_count:
            print(f"      ({fmt_num(confirmed_count)} of "
                  f"{fmt_num(len(skipped))} confirmed by anchor lookup "
                  f"in both outputs)")
    else:
        print(f"   ⚠️  NOT fully explained by skipped variants.")
        p = first_divergent_prefix(seq_a, seq_b)
        left = seq_a[max(0, p - LEFT_CTX):p]
        right_a = seq_a[p + 1:p + 1 + RIGHT_CTX]
        disp_a, disp_b = make_display_pair_single(
            args.label_a, args.label_b,
            left, seq_a[p:p + 1], seq_b[p:p + 1], right_a,
        )
        print()
        print(f"🔎 First base divergence at output position {fmt_num(p)}:")
        print()
        print(f"  {disp_a}")
        print(f"  {disp_b}")
        cc = caret_column_single(disp_a)
        if cc >= 0:
            print(' ' * (2 + cc) + '^')

    if args.tsv:
        if explained:
            first_pos = str(first_confirmed_pos or skipped[0]['skip_v']['pos'])
            emit_tsv(hap_label, len_a, len_b, delta,
                     str(len(skipped)), first_pos,
                     first_disp_a, first_disp_b)
        else:
            emit_tsv(hap_label, len_a, len_b, delta, "LENGTH_MISMATCH",
                     str(first_divergent_prefix(seq_a, seq_b)),
                     disp_a, disp_b)

    return 0 if explained else 1


if __name__ == '__main__':
    sys.exit(main())