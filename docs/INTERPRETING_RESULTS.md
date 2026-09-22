# Interpreting Benchmark Results for vcf2fasta-rust

This document explains the log formats produced by each benchmark script
and how to read the validation verdicts.

For the figures themselves and the methodology behind them, see
[Benchmark Results](BENCHMARK_RESULTS.md). For how to regenerate the
figures, see [Benchmarking Guide §4.6](BENCHMARKING.md#46-plot_benchmark_resultsr--regenerate-all-figures).

---

## 1. Where results are written

| Script | Log file | Output directories |
|--------|----------|--------------------|
| `run_benchmark.sh` | `benchmark_results/benchmark_results.log` | `benchmark_results/new_<engine>_t<N>/<sample>/` |
| `run_benchmark.sh` (+ `USE_VCFLIB=true`) | `results/compare_summary.tsv` | (aggregate, no per-run directory) |
| `scaling_benchmark.sh` | `scaling_results/scaling_benchmark.log` | `scaling_results/samples_<N>/<chrom>/` |
| `unphased_scaling_benchmark.sh` | `unphased_scaling_results/unphased_scaling_benchmark.log` | `unphased_scaling_results/<dataset>/` |
| `plot_benchmark_results.R` | stdout (range checks + preflight) | `docs/figures/` |

Each script also writes stderr/stdout capture files per run under an
`errors/` subdirectory.

After all runs are complete, the raw TSV logs are aggregated by hand into
three Excel workbooks committed at `tests/benchmark/results/`. Each
workbook has five sheets (`RunTime1` … `RunTime4`, `Combined`). See
[Benchmarking Guide §3](BENCHMARKING.md#3-benchmark-results) for the
workbook layout, and
[Benchmark Results](BENCHMARK_RESULTS.md) for the figures built from them.

---

## 2. `run_benchmark.sh` log

Tab-separated. One row per (chromosome, sample).

| Column | Meaning |
|--------|---------|
| `Chromosome` | Contig tested (e.g., `chr22`). |
| `Sample` | Sample name. |
| `Ploidy` | Inferred from output files: `2` if both `:0.fa` and `:1.fa` exist, `1` if only one. |
| `vcflib` | Wall-clock time of the official `vcf2fasta` (only if `USE_VCFLIB=true`). |
| `cpu_t<N>_time` / `gpu_t<N>_time` | Wall-clock time per configuration, in the form `MMm SS.sss`. |
| `speedup_<engine>_t<N>_vs_t<M>` | Baseline time / this time. Higher is better. |
| `speedup_<engine>_t<N>_vs_vcflib` | vcflib time / this time. |
| `Notes` | `Success`, `Missing tool FASTA output`, `vcf2fasta-rust Validation Failed`, etc. |

### Validation verdicts

`run_benchmark.sh` captures the full output of
`validate_vcf_to_fasta.py` into a shell variable and prints its own
one-line summary — it does **not** forward the validator's output
verbatim. A user running the benchmark normally sees:

```
        ✅ vcf2fasta-rust Validation: PASSED
        ✅ vcflib Validation: PASSED
```

or, on failure:

```
        ❌ vcf2fasta-rust Validation: FAILED
        ── validator output (first 20 lines) ──
           (first 20 lines of the validator's captured output)
```

The validator's own verdict block (`validate_vcf_to_fasta.py` prints
this when invoked directly) reads:

```
📊 Results:
   TOOL   vs expected: ✅ PASS
   VCFLIB vs expected: ✅ PASS
   ✅ Tool and vcflib both agree with expected (cross-checked)
```

- `TOOL vs expected: ✅ PASS` — byte-for-byte match between the Rust
  tool's FASTA and the expected sequence built from the VCF using the
  tool's own policy.
- `TOOL vs expected: ❌ FAIL` — mismatch. The validator prints the first
  20 mismatching byte offsets, their originating variants, and 20-byte
  context windows. Those diagnostics appear *before* the block.
- `VCFLIB vs expected:` line — only present when `--vcflib-hap0`
  (and optionally `--vcflib-hap1`) were supplied.
- The trailing line relates the two verdicts:

  | Line | Meaning |
  |------|---------|
  | `✅ Tool and vcflib both agree with expected (cross-checked)` | Both tools produce the expected sequence. |
  | `⚠️  Tool agrees with expected; vcflib diverges` | The Rust tool matches the expected sequence; vcflib does not. |
  | `⚠️  vcflib agrees with expected; tool diverges` | The reverse. |
  | `❌ Both tool and vcflib diverge from expected` | Neither matches. |
  | `ℹ️  vcflib not available; tool-only check` | `--vcflib-hap*` flags were not supplied. |

To see the block in full, invoke the validator directly:

```bash
python3 tests/benchmark/validation/validate_vcf_to_fasta.py \
    <vcf> <my_hap0> [<my_hap1>] \
    --reference <reference.fa> \
    --sample SAMPLE \
    --no-validate-ref \
    [--vcflib-hap0 VC0.fa] [--vcflib-hap1 VC1.fa]
```

Pass only `my_hap0` for haploid chromosomes.

Exit codes: `0` all requested comparisons passed, `1` at least one
failed, `2` setup error (missing file, bad VCF, missing reference).

### Length-diff attribution

When `USE_VCFLIB=true` and the validator reports a `VCFLIB vs expected:
❌ FAIL`, `run_benchmark.sh` additionally invokes
`explain_length_diff.py` for each haplotype. Its output is printed
inline (indented under the sample) and appended to
`compare_summary.tsv`.

The inline output looks like:

```
     └─ Length-diff explainer H0: NA12877_chr12_H0
        ==============================================================================
        Length-diff explainer: NA12877_chr12_H0
          Rust       length: 133,273,707 bases
          vcflib     length: 133,273,708 bases
          Δ (vcflib − Rust):    +1 bases
        
        Variants in VCF:                                193,899
        Skipped by Rust overlap policy:              1
        
        ── Skipped variants (each falls inside a previous variant's REF span) ──
        
           [1] chr12:32,192,430  REF=T  ALT=TTAAA  GT=0|1
                applies first at  chr12:32,192,430  REF=T  ALT=TTAAA  GT=0|1
               Rust       first=T, duplicate=∅ (skipped)
               vcflib     first=T, duplicate=T (applied again)
        
          Rust:     ...CTCCAGCCTGGGTGGCAGAGCGAGACTCTG[T][∅]TAAATAAATAAATAAATAAATAAATAGATA...
          vcflib:   ...CTCCAGCCTGGGTGGCAGAGCGAGACTCTG[T][T]TAAATAAATAAATAAATAAATAAATAGATA...
                                                        ^
        
               ✅ confirmed: vcflib anchor at output position 32,191,763 vs Rust at 32,191,762 = +1 base
        
        ── Summary ──
           Predicted Δ from skipped variants:  +1 bases
           Actual Δ:                            +1 bases
           ✅ The base-count difference is FULLY EXPLAINED by vcflib applying
              the 1 variant(s) that Rust skips per its overlap policy.
              (1 of 1 confirmed by anchor lookup in both outputs)
```

Reading the bracket-form display:

- Two brackets are shown per line. The first contains the **first
  emission** of the variant, which both tools apply identically. The
  second contains either `∅` (Rust: nothing, because it skipped the
  duplicate) or the duplicate's emitted bytes (vcflib: re-applied the
  duplicate).
- The caret points at the second bracket — the point of divergence.
- Both sequences are anchored to the same reference context
  (`...CTCCAGCCTGGGTGGCAGAGCGAGACTCTG` on the left,
  `TAAATAAATAA...` on the right), so you can read the two lines
  side by side and see the single extra base in vcflib's output.

The two tools will diverge **only** on variants that the Rust tool skips
under its overlap policy (policies 25 and 26 in the
[User Guide](./USER_GUIDE.md#edge-case-handling-policy-table)). Every
variant that is not part of an overlap produces identical bytes in both
outputs.

---

## 3. `scaling_benchmark.sh` log

Tab-separated. One row per (chromosome, sample-count).

| Column | Meaning |
|--------|---------|
| `Chromosome` | Contig tested. |
| `Samples` | Number of samples selected from the 1000G VCF. |
| `Mode` | `cpu`, `gpu`, or `both`. |
| `Extraction_Time` | Time for `bcftools view -s ...`. **Not** part of the CPU/GPU timings. |
| `gpu_dev<L>_t<N>_time` | GPU run with device prefix `<L>` (commas → underscores) and `<N>` host threads. |
| `cpu_t<N>_time` | CPU run with `<N>` threads. |
| `Notes` | `Success` or per-configuration error. |

The `<L>` prefix is the CUDA device list passed via `GPU_DEVICES`
(commas replaced by underscores). For example, `gpu_dev0_1_t2_time`
represents a two-device run with 2 host threads; `gpu_dev0_t2_time` is
single-device. Runs collected before `--gpu-devices` was made binding
will show the prefix the script passed but represent single-GPU
execution regardless.

---

## 4. `unphased_scaling_benchmark.sh` log

Tab-separated. One row per dataset (one dataset per invocation).

| Column | Meaning |
|--------|---------|
| `Dataset` | `tiny`, `small`, `medium`, `larger`, `ploidy_edge`. |
| `Samples` | Number of samples in the dataset. |
| `Chromosomes` | Number of contigs. |
| `Mode` | `cpu`, `gpu`, or `both`. |
| `gpu_dev<L>_t<N>_time` | GPU run with device prefix `<L>` and `<N>` host threads. |
| `cpu_t<N>_time` | CPU run with `<N>` threads. |
| `Notes` | `Success` or error. |

Because phasing dominates on unphased input, expect smaller CPU-vs-GPU
differences here than on the real 1000G data. The `<L>` prefix carries
the same meaning as in `scaling_benchmark.sh`.

---

## 5. Reading speedup columns

Speedup is `baseline_time / this_time`.

- `1.00x` = no change.
- `2.00x` = half the time.
- `< 1.00x` = slower (typical for GPU on small datasets, or CPU with
  more threads than cores).

---

## 6. Reading warnings

Every `<prefix>.warnings.log` line has the form:

```
WARNING [chr22:12345] <reason>: <details>
```

The `<prefix>.log` file also contains one `[WARNINGS]` summary per contig:

```
[WARNINGS] chr22: total=1423 {overlap=1200, ref_mismatch=200, invalid_allele_index=23}
```

Reason tags (from `report.rs::classify_warning`):

| Tag | Meaning |
|-----|---------|
| `overlap` | Variant overlaps one already applied (policies 25/26). |
| `ref_mismatch` | VCF REF disagrees with FASTA (policy 10). |
| `ref_interval_overflow` | REF extends past contig end (policy 2). |
| `pos_out_of_range` | POS outside contig. |
| `invalid_pos` | POS ≤ 0 or non-numeric. |
| `malformed_record` / `malformed_gt` | Record/GT could not be parsed. |
| `invalid_allele_index` | Allele index exceeds ALT count. |
| `negative_allele` | Allele index is negative. |
| `missing_allele` | GT contains `.`. |
| `sample_missing` | Sample column absent from the record. |
| `invalid_ref` | REF is not a valid nucleotide sequence. |

If you see a large `ref_mismatch` count, the reference build likely does
not match the VCF's.

---

## 7. Reading `compare_summary.tsv`

Produced by `run_benchmark.sh` when `USE_VCFLIB=true`. One `COMPARE` row
per `(sample, chromosome, haplotype)`.

### 7.1. Column reference

| Column | Meaning |
|--------|---------|
| `Row_type` | Always `COMPARE`. |
| `Haplotype_id` | `<sample>_<chromosome>_H<hap>`. |
| `Rust_bases` | Output length in bases for the Rust tool. |
| `vcflib_bases` | Output length in bases for vcflib. |
| `Delta_bases` | `vcflib_bases − Rust_bases`, signed with `+`/`-`. |
| `Skipped_variant_count` | Number of variants the Rust tool skipped under its overlap policy. |
| `First_skipped_vcf_pos` | VCF POS of the first such skipped variant, or `-`. |
| `Rust_display` | Bracket-form alignment snippet for the Rust output. |
| `vcflib_display` | Bracket-form alignment snippet for the vcflib output. |

### 7.2. Common patterns

| Scenario | What the row looks like | Meaning |
|----------|-------------------------|---------|
| **Identical outputs** | `Delta_bases = +0`, `Skipped_variant_count = 0`, displays both `-` | The two tools agreed byte-for-byte. |
| **Δ fully explained** | `Delta_bases = +N`, `Skipped_variant_count = k`, `First_skipped_vcf_pos` set, displays populated | The Rust tool skipped `k` overlapping variants; vcflib applied them, producing `N` extra bases. This is expected and correct. |
| **Δ not explained** | `Delta_bases = ±N`, `Skipped_variant_count = 0`, displays show a `LENGTH_MISMATCH` snippet | A bug or an untested policy difference. Investigate. |
| **Content differs at same length** | `Delta_bases = +0`, `Skipped_variant_count = 0`, `First_skipped_vcf_pos` is a numeric position, displays show the mismatching base | The two tools produced the same length but different bytes at some position. Investigate. |

### 7.3. Interpreting the display snippets

The `Rust_display` and `vcflib_display` columns contain the same
bracket-form snippet printed in the terminal output described in [§2](#length-diff-attribution):

```
Rust:    ...CTCCAGCCTGGGTGGCAGAGCGAGACTCTG[T][∅]TAAATAAATAAATAAATAAATAAATAGATA...
vcflib:  ...CTCCAGCCTGGGTGGCAGAGCGAGACTCTG[T][T]TAAATAAATAAATAAATAAATAAATAGATA...
                                              ^
```

The display is UTF-8; `∅` is the empty-set symbol (`U+2205`) and marks
where a duplicate was **not** applied. `[T]` is a literal base. The
specific base shown depends on the variant in the input VCF — `[T]` here
is the first base of the REF allele `T` (ALT `TTAAA`) at the overlapping
position, and `[T]` on the vcflib line is the same base re-emitted by
vcflib when it applies the duplicate a second time.

---