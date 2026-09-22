# Benchmarking Guide for vcf2fasta-rust

Everything you need to run the benchmark suite, validate output, and
regenerate figures.

---

## 1. Prerequisites

Complete the [Installation Guide](INSTALLATION.md). In particular:

- `cargo build --release` (or `--features cuda`) has been run.
- `samtools`, `bcftools`, `bgzip`, `tabix`, `bwa`, `art_illumina`,
  `whatshap` are on `PATH`.
- `BEAGLE_JAR` is exported if you plan to use Beagle phasing.
- R with the plotting packages installed (only needed if you want to
  regenerate figures).
- For GPU benchmarks, ensure your CUDA runtime is on `LD_LIBRARY_PATH`:

  ```bash
  export LD_LIBRARY_PATH=/usr/local/cuda-12.8/targets/x86_64-linux/lib:/usr/local/cuda-12.8/lib64:$LD_LIBRARY_PATH
  ```

---

## 2. Directory layout

```
vcf2fasta-rust/
├── docs/
│   ├── BENCHMARKING.md
│   ├── BENCHMARK_RESULTS.md
│   ├── INTERPRETING_RESULTS.md
│   └── figures/                          
│       ├── fig1_runtime_chromosomes.png
│       ├── fig2_scaling.png
│       ├── fig3_speedup.png
│       ├── fig4_unphased.png
│       ├── fig4_workload.png
│       └── fig_combined.png
└── tests/
    └── benchmark/
        ├── scripts/
        │   ├── download_data.sh
        │   ├── generate_synthetic_datasets.py
        │   ├── generate_synthetic_datasets.sh
        │   ├── run_benchmark.sh
        │   ├── run_merged_output_check.sh
        │   ├── scaling_benchmark.sh
        │   ├── unphased_scaling_benchmark.sh
        │   └── plot_benchmark_results.R
        ├── validation/
        │   ├── explain_length_diff.py       # called by run_benchmark.sh
        │   └── validate_vcf_to_fasta.py     # called by run_benchmark.sh
        ├── results/                         
        │   ├── results_benchmark.xlsx
        │   ├── results_scaling.xlsx
        │   ├── results_unphased.xlsx
        │   └── compare_summary.tsv          # per-haplotype divergence summary
        └── datasets/                        
```

The three `validation/` scripts and `generate_synthetic_datasets.py`
are **called automatically** by the other scripts and are not intended
to be invoked by hand:

| Script | Called by |
|--------|-----------|
| `validate_vcf_to_fasta.py` | `run_benchmark.sh` |
| `explain_length_diff.py` | `run_benchmark.sh` (when `USE_VCFLIB=true`) |
| `generate_synthetic_datasets.py` | `generate_synthetic_datasets.sh` |

---

## 3. Benchmark results

The benchmark scripts produce raw TSV logs and one aggregate TSV. These
are committed under `tests/benchmark/results/` (the `compare_summary.tsv`
file) and are the input to the plotting script and the authoritative
record of each benchmark run.

### 3.1. The three workbooks

| File | Benchmark | Source TSV log |
|------|-----------|----------------|
| `results_benchmark.xlsx` | Platinum per-sample, per-chromosome (with vcflib cross-check) | `benchmark_results/benchmark_results.log` |
| `results_scaling.xlsx` | 1000G chr22, scaling with sample count | `scaling_results/scaling_benchmark.log` |
| `results_unphased.xlsx` | Synthetic unphased workloads (`tiny` … `ploidy_edge`) | `unphased_scaling_results/unphased_scaling_benchmark.log` |

### 3.2. Workbook structure

Each workbook has **five sheets**:

| Sheet | Contents |
|-------|----------|
| `RunTime1` | Raw per-configuration timings from run 1 |
| `RunTime2` | Raw per-configuration timings from run 2 |
| `RunTime3` | Raw per-configuration timings from run 3 |
| `RunTime4` | Raw per-configuration timings from run 4 |
| `Combined` | Mean ± SD across the four runs, in the form `MMm SS.sss ± MMm SS.sss` |

The per-run sheets preserve the individual measurements; `Combined` is
the aggregation that `plot_benchmark_results.R` reads. Keeping both
means you can:

- Audit the aggregation if a `Combined` value looks wrong.
- Recompute means and SDs with a different method (e.g., median, IQR)
  without re-running the benchmark.
- Inspect a single anomalous run without re-running the benchmark.

### 3.3. `compare_summary.tsv`

Produced by `run_benchmark.sh` when `USE_VCFLIB=true`. One row per
`(sample, chromosome, haplotype)`, written by `explain_length_diff.py`.

Columns (tab-separated):

| Column | Meaning |
|--------|---------|
| `Row_type` | Always `COMPARE` |
| `Haplotype_id` | `<sample>_<chromosome>_H<hap>` |
| `Rust_bases` | Output length in bases for the Rust tool |
| `vcflib_bases` | Output length in bases for vcflib |
| `Delta_bases` | `vcflib_bases − Rust_bases`, signed |
| `Skipped_variant_count` | Number of variants the Rust tool skipped per its overlap policy |
| `First_skipped_vcf_pos` | VCF POS of the first such skipped variant, or `-` |
| `Rust_display` | Bracket-form alignment snippet for the Rust output |
| `vcflib_display` | Bracket-form alignment snippet for the vcflib output |

Interpretation guidance is in [Interpreting Results §7](INTERPRETING_RESULTS.md#7-reading-compare_summarytsv).

---

## 4. Script reference

### 4.1. `download_data.sh` — fetch reference and public datasets

Downloads GRCh38, Platinum Genomes VCFs (NA12877, NA12878), 1000 Genomes
chr22, and Beagle GRCh38 genetic maps.

```bash
cd tests/benchmark/scripts
./download_data.sh
```

Run this **once**. Produces everything under `tests/benchmark/datasets/`.

---

### 4.2. `generate_synthetic_datasets.sh` — build synthetic datasets

Generates a self-contained dataset for the integrated phasing pipeline.
This shell script wraps `generate_synthetic_datasets.py`; run the shell
script, not the Python file.

```bash
cd tests/benchmark/scripts
DATASET=tiny ./generate_synthetic_datasets.sh
```

Environment variables:

| Variable | Default | Meaning |
|----------|---------|---------|
| `DATASET` | `tiny` | `tiny`, `small`, `medium`, `larger`, `ploidy_edge`. |
| `SEED` | `12345` | RNG seed for reproducible genotypes. |
| `THREADS` | `4` | Threads for BWA. |
| `COVERAGE` | `5` | Per-haplotype coverage for ART. |
| `SOURCE_VCF` | `datasets/platinum/NA12877/NA12877.vcf.gz` | Source variant pool. |
| `SOURCE_REFERENCE` | `datasets/reference/GRCh38.primary_assembly.genome.fa` | Source reference. |
| `OUTPUT_ROOT` | `datasets/synthetic_unphased` | Where the dataset is written. |

Output: `datasets/synthetic_unphased/<DATASET>/` with `reference/`,
`vcf/`, `truth_haplotypes/`, `reads/`, `bam/`, and `metadata.json`.

Build every dataset:

```bash
for D in tiny small medium larger ploidy_edge; do
    DATASET=$D ./generate_synthetic_datasets.sh
done
```

---

### 4.3. `run_benchmark.sh` — per-sample benchmark with optional vcflib cross-check

Runs the Rust tool on each sample/chromosome **in isolation**. When
`USE_VCFLIB=true`, it also runs the official vcflib `vcf2fasta` and
compares the two outputs, then invokes `explain_length_diff.py` to
attribute any length differences to the Rust tool's overlap policy.

```bash
cd tests/benchmark/scripts
MODE=cpu USE_VCFLIB=true SAMPLES="NA12877 NA12878" CHROMS="chr22" ./run_benchmark.sh
```

Environment variables:

| Variable | Default | Meaning |
|----------|---------|---------|
| `MODE` | `cpu` | `cpu`, `gpu`, or `both`. |
| `USE_VCFLIB` | `false` | Also run the official `vcf2fasta` (must be on `PATH`). |
| `DATASET` | `platinum` | `platinum`. |
| `SAMPLES` | `NA12877 NA12878` | Space-separated sample names. |
| `CHROMS` | `chr22` | Space-separated contigs. |
| `CPU_THREADS_LIST` | `1 2 4 8` | CPU thread sweep. |
| `GPU_THREADS_LIST` | `1 2` | GPU host-thread sweep. |
| `GPU_DEVICES` | `0` | Comma-separated GPU ids. |
| `CHUNK_FILES` | `true` | Slice the reference to the target contig before running. |
| `LIST_SKIPPED` | `20` | Max skipped variants to print per sample. |
| `COMPARE_HAPS` | `0 1` | Haplotypes to feed to `explain_length_diff.py`. |
| `COMPARE_MAX_DIFFS` | `20` | Max skipped variants to list per haplotype. |
| `COMPARE_TSV` | `results/compare_summary.tsv` | Aggregation TSV path. |
| `COMPARE_SCRIPT` | auto | Path to `explain_length_diff.py`. |

Outputs:

- `tests/benchmark/benchmark_results/benchmark_results.log` — the main log.
- `tests/benchmark/benchmark_results/new_<engine>_t<N>/<sample>/` — FASTA outputs.
- `tests/benchmark/benchmark_results/vcflib_out/<sample>/` — vcflib FASTA outputs (if `USE_VCFLIB=true`).
- `tests/benchmark/results/compare_summary.tsv` — per-haplotype divergence summary (if `USE_VCFLIB=true`).

At the end it prints a summary table of Rust-vs-vcflib length differences
and, for each divergent haplotype, the bracket-form alignment snippet
produced by `explain_length_diff.py`.

This is the **only** place where the official vcflib `vcf2fasta` is
compared against — the correctness regression check lives here, not in
`cargo test`.

---

### 4.4. `scaling_benchmark.sh` — combinatorial CPU/GPU scaling

Sweeps sample counts × thread counts × GPU device prefixes on a fixed
chromosome (default chr22 from the 1000 Genomes VCF).

```bash
cd tests/benchmark/scripts
MODE=auto GPU_DEVICES=0,1 ./scaling_benchmark.sh
```

Environment variables:

| Variable | Default | Meaning |
|----------|---------|---------|
| `MODE` | `auto` | `auto`, `cpu`, `gpu`, `both`. |
| `SAMPLES_LIST` | `10 50 100` | Sample-count tiers. |
| `CHROMS` | `chr22` | Space-separated contigs. |
| `CPU_THREADS_LIST` | `1 2 4 8` | CPU thread sweep. |
| `GPU_THREADS_LIST` | `1 2 4` | GPU host-thread sweep. |
| `GPU_DEVICES` | `0,1` | Pool of GPU indices; cumulative prefixes are swept. |
| `OUTPUT_BASE` | `scaling_results` | Where FASTA outputs go (use `/dev/shm` for speed). |

Result: `tests/benchmark/scaling_results/scaling_benchmark.log`.

**Note:** the tool honours `--gpu-devices` as of the fix that made it
binding. Benchmarks run before that fix used device 0 only; results
generated with `GPU_DEVICES=0,1` on the old binary would have been
single-GPU runs labelled as multi-GPU. Multi-GPU scaling only helps
when a single contig yields enough tiles to keep every stream busy.

---

### 4.5. `unphased_scaling_benchmark.sh` — integrated phasing + vcf2fasta

Runs the whole pipeline (phasing → vcf2fasta) on one synthetic dataset,
sweeping CPU threads and GPU prefixes.

**Prerequisite**:  the dataset named by `DATASET` must already exist
under `datasets/synthetic_unphased/<DATASET>/`. If it does not, run
`DATASET=<name> ./generate_synthetic_datasets.sh` first. This script
does **not** invoke the generator itself — it only reads what the
generator produced.

```bash
cd tests/benchmark/scripts
DATASET=medium MODE=both GPU_DEVICES=0,1 ./unphased_scaling_benchmark.sh
```

Environment variables:

| Variable | Default | Meaning |
|----------|---------|---------|
| `DATASET` | `tiny` | Synthetic dataset name. |
| `MODE` | `auto` | `auto`, `cpu`, `gpu`, `both`. |
| `CPU_THREADS_LIST` | `1 2 4 8` | CPU thread sweep. |
| `GPU_THREADS_LIST` | `1 2 4` | GPU host-thread sweep. |
| `GPU_DEVICES` | `0,1` | Pool of GPU indices; cumulative prefixes are swept. |

Result: `tests/benchmark/unphased_scaling_results/unphased_scaling_benchmark.log`.

**Note:** as with `scaling_benchmark.sh`, the tool honours
`--gpu-devices` only as of the fix that made it binding. Because
phasing dominates on unphased input, expect smaller CPU-vs-GPU
differences here than on already-phased data, and small multi-GPU gains
unless the cohort is large.

---

### 4.6. `plot_benchmark_results.R` — regenerate all figures

Reads the three committed workbooks from `tests/benchmark/results/` and
regenerates every figure in `docs/figures/`.

**Prerequisites:** R ≥ 4.3 with `readxl`, `dplyr`, `tidyr`, `stringr`,
`ggplot2`, `patchwork`, `scales`, and `forcats`. See
[Installation Guide §5](INSTALLATION.md#5-r-and-r-packages-for-figure-regeneration)
for setup.

**Inputs.** The script reads the `Combined` sheet from each workbook:

| File | Sheet read | Produced by |
|------|-----------|-------------|
| `results_benchmark.xlsx` | `Combined` | `run_benchmark.sh` |
| `results_scaling.xlsx` | `Combined` | `scaling_benchmark.sh` |
| `results_unphased.xlsx` | `Combined` | `unphased_scaling_benchmark.sh` |

The other four sheets in each workbook (`RunTime1` … `RunTime4`) are
ignored by the script but are kept for auditability.

**Usage:**

```bash
cd tests/benchmark/scripts
Rscript plot_benchmark_results.R
```

By default the script reads from `../../results/` and writes to
`../../../docs/figures/`. Override either with environment variables if
needed:

| Environment variable | Default | Meaning |
|----------------------|---------|---------|
| `RESULTS_DIR` | `../../results` | Directory containing the three `.xlsx` files. |
| `FIGURES_DIR` | `../../../docs/figures` | Where PDF and PNG files are written. |

The script:

1. Creates `FIGURES_DIR` if it does not exist.
2. Parses the `MMm SS.sss ± MMm SS.sss` time strings into mean and SD in
   seconds.
3. Emits every figure as **both PDF and PNG** (PNG at 300 DPI).
4. Prints a per-figure **range check** (min / max / ratio of runtimes) to
   help spot parsing mistakes.
5. Runs a **preflight** report listing any `NA` values and any cells
   where SD > mean.

**Figures produced:**

| File | Description |
|------|-------------|
| `fig1_runtime_chromosomes.pdf` / `.png` | Platinum runtime per chromosome, faceted by sample |
| `fig1_alt.pdf` / `.png` | Same data, chromosome on facet, sample on X |
| `fig2_scaling.pdf` / `.png` | 1000G chr22 scaling with sample count (log-log) |
| `fig3_speedup.pdf` / `.png` | Speedup vs vcflib per chromosome |
| `fig4_unphased.pdf` / `.png` | Synthetic unphased workloads, categorical X |
| `fig4_workload.pdf` / `.png` | Synthetic unphased workloads, workload size on X |
| `fig_combined.pdf` / `.png` | Multi-panel composite (A=Fig1, B=Fig3, C=Fig2, D=Fig4) |
---

### 4.7. `run_merged_output_check.sh` — verify `--merged-output` correctness

Runs the tool twice on the same sliced VCF (once in per-contig mode,
once with `--merged-output`) and verifies that every merged file is
byte-identical to the concatenation of the corresponding per-contig
files, in VCF-header contig order.

```bash
cd tests/benchmark/scripts
./run_merged_output_check.sh
```

Environment variables:

| Variable | Default | Meaning |
|----------|---------|---------|
| `INPUT_VCF` | NA12877 Platinum VCF | Input VCF. |
| `REFERENCE` | GRCh38 primary assembly | Reference FASTA. |
| `MODE` | `cpu` | `cpu`, `gpu`, or `both`. |
| `GPU_DEVICES` | `0` | Comma-separated GPU ids. |
| `THREADS` | `2` | Host threads. |
| `SAMPLES` | first sample in VCF | Space-separated sample names. |
| `CONTIGS` | first 2 in VCF header | Space-separated contigs. |
| `MAX_CONTIGS` | `2` | Cap when `CONTIGS` is unset. |
| `NO_VALIDATE_REF` | `false` | Skip REF-vs-FASTA validation. |
| `WORKDIR` | `tests/benchmark/merged_check` | Scratch directory. |

The script checks, for each `(sample, haplotype)`:

1. The merged file exists.
2. `diff` between the merged file and the concatenation of per-contig
   files (headers and newlines included) reports no differences.
3. The merged file has exactly one header per present contig.

Result: `tests/benchmark/merged_check/merged_check.log`. The script exits
non-zero on any failure. This is the only script that exercises
`--merged-output`. The other benchmark scripts use per-contig mode.

---

### 4.8. `validate_vcf_to_fasta.py` — expected-sequence validator

Called automatically by `run_benchmark.sh`. Builds the expected sequence
per haplotype from the VCF using the Rust tool's own policy, then
compares byte-for-byte against the tool's FASTA output. Optionally also
compares vcflib's FASTA output against the same expected sequence as an
independent cross-check.

Not intended to be invoked by hand, though it can be:

```bash
python3 tests/benchmark/validation/validate_vcf_to_fasta.py \
    <vcf> <my_hap0> [<my_hap1>] \
    --reference <reference.fa> \
    [--sample SAMPLE] \
    [--no-validate-ref] \
    [--vcflib-hap0 VC0.fa] [--vcflib-hap1 VC1.fa] \
    [--max-errors N] [--list-skipped N] [--strict]
```

Pass only `my_hap0` for haploid chromosomes (male chrX/chrY, chrM).
Exit codes: `0` all comparisons passed, `1` at least one failed,
`2` setup error.

Full output interpretation is in
[Interpreting Results §2](INTERPRETING_RESULTS.md#2-run_benchmarksh-log).

---

### 4.9. `explain_length_diff.py` — attribute length differences

Called automatically by `run_benchmark.sh` when `USE_VCFLIB=true`. Given
a VCF, a Rust FASTA, and a vcflib FASTA for one `(sample, chromosome,
haplotype)`, it:

1. Identifies variants the Rust tool skipped under its overlap policy.
2. Predicts the byte difference each skipped variant should cause.
3. Confirms the prediction by locating an anchor sequence in both
   outputs and measuring the observed offset.
4. Prints a bracket-form alignment snippet showing both the first
   (applied) emission and the skipped (duplicate) emission, so the
   divergence is visible at the byte level.

```bash
python3 tests/benchmark/validation/explain_length_diff.py \
    <vcf> <rust.fa> <vcflib.fa> \
    --reference <reference.fa> \
    --sample SAMPLE \
    --haplotype {0|1} \
    [--chromosome chrN] \
    [--label-a Rust] [--label-b vcflib] \
    [--max-diffs N] [--tsv]
```

Exit codes: `0` Δ fully explained, `1` Δ not fully explained or content
differs at same length, `2` setup error.

With `--tsv`, it prints one `COMPARE` row per invocation, which
`run_benchmark.sh` appends to `compare_summary.tsv`. The format is
documented in [§3.3](#33-compare_summarytsv) and interpreted in
[Interpreting Results §7](INTERPRETING_RESULTS.md#7-reading-compare_summarytsv).

---

## 5. Recommended workflow

```bash
cd tests/benchmark/scripts

# 1. One-time setup
./download_data.sh

# 2. Generate synthetic datasets
for D in tiny small medium larger ploidy_edge; do
    DATASET=$D ./generate_synthetic_datasets.sh
done

# 3. Scaling on real 1000G data
MODE=both ./scaling_benchmark.sh

# 4. Phasing pipeline on synthetic data
for D in tiny small medium larger ploidy_edge; do
    DATASET=$D MODE=both ./unphased_scaling_benchmark.sh
done

# 5. Per-sample validation vs vcflib, plus length-diff attribution
MODE=cpu USE_VCFLIB=true ./run_benchmark.sh

# 6. Build the three workbooks from the raw TSV logs, drop them in
#    tests/benchmark/results/, then regenerate every figure:
Rscript plot_benchmark_results.R
```

---

## 6. Notes

- **Do not call `validate_vcf_to_fasta.py` directly.** It is invoked by
  `run_benchmark.sh`.
- **Do not call `explain_length_diff.py` directly.** It is invoked by
  `run_benchmark.sh` when `USE_VCFLIB=true`.
- **Do not call `generate_synthetic_datasets.py` directly.** Use the
  shell wrapper.
- Set `OUTPUT_BASE=/dev/shm` on Linux to benchmark against a RAM-backed
  filesystem.
- All scripts print a summary and write a TSV log.
- The plotting script is the only step that needs R; the tool and the
  benchmark scripts do not.
- **The three `.xlsx` workbooks and `compare_summary.tsv` are the
  committed record of the benchmark.** The `.log` files they were built
  from are not committed; they can be regenerated by re-running the
  corresponding script.

---