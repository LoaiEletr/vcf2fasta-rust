# vcf2fasta-rust

A high-performance, parallel tool to reconstruct per-sample, per-haplotype FASTA sequences from an indexed VCF file.  
This is a Rust re-implementation of the popular `vcf2fasta` tool from [vcflib](https://github.com/vcflib/vcflib), offering:

- **Blazing speed** – multi-threaded CPU, single- or multi-GPU CUDA
  acceleration (`--gpu-devices 0,1`), with a `--max-vram` cap.
- **Automatic phasing** – chooses Beagle, WhatsHap, or variable-ploidy phasing per contig based on input data.
- **Memory efficiency** – streaming CPU path with bounded memory; GPU path uses reusable pinned buffers.
- **Full correctness** – output is byte-for-byte identical to the official tool (verified by the integration tests in `tests/`).

---

## 📚 Documentation

- **[Installation Guide](docs/INSTALLATION.md)** – full setup (Rust, conda/mamba, Beagle, GPU, R for plots).
- **[User Guide](docs/USER_GUIDE.md)** – usage, VCF/FASTA basics, CLI options, edge cases.
- **[Architecture Guide](docs/ARCHITECTURE.md)** – internal design, data flow, scheduling, GPU pipeline.
- **[Benchmarking Guide](docs/BENCHMARKING.md)** – how to run every benchmark script and regenerate figures.
- **[Benchmark Results](docs/BENCHMARK_RESULTS.md)** – methodology and figures from the Colab benchmark run.
- **[Interpreting Results](docs/INTERPRETING_RESULTS.md)** – how to read logs, validation output, and plotting diagnostics.

---

## 🚀 Quick Start

### 1. Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
```

### 2. Install conda / mamba

If you do not already have conda or mamba, install **Miniforge**, which
ships with `mamba` preconfigured for the `conda-forge` channel:

```bash
# Linux x86_64
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
```

For macOS, replace `Linux-x86_64` with `MacOSX-arm64` (Apple Silicon) or
`MacOSX-x86_64` (Intel).

Add conda to your shell so future terminals see `conda activate`:

```bash
conda init bash    # or zsh, fish, etc.
exec "$SHELL" -l   # reload the shell
```

### 3. Create a dedicated environment

```bash
mamba create -n vcf2fasta python=3.12 -y
mamba activate vcf2fasta
```

Do **not** install into `base`. A dedicated environment keeps this
project's tool versions isolated from anything else on the machine.

### 4. Install bioinformatics tools

With the `vcf2fasta` environment activated:

```bash
mamba install -y -c conda-forge -c bioconda \
    samtools \
    bcftools \
    vcflib \
    whatshap \
    beagle \
    bwa \
    pysam \
    art
```

| Tool | Used for |
|------|----------|
| `samtools` | FASTA indexing, BAM handling |
| `bcftools` | VCF slicing, merging, validation |
| `bgzip` + `tabix` | VCF compression/indexing |
| `vcflib` | Reference `vcf2fasta` for benchmark cross-checks |
| `whatshap` | Read-based phasing |
| `beagle` | Population phasing (Java) |
| `bwa` | Read alignment (synthetic dataset generation) |
| `art` (`art_illumina`) | Read simulation (synthetic dataset generation) |
| `pysam` | Python VCF/FASTA parsing (validation scripts) |

### 5. Install system build dependencies

The Rust build links against `libclang` via `rust-htslib`.

**Ubuntu / Debian:**

```bash
sudo apt-get update
sudo apt-get install -y build-essential pkg-config libclang-dev
```

**macOS:**

```bash
brew install llvm pkg-config
export LIBCLANG_PATH="$(brew --prefix llvm)/lib"
```

### 6. Install Beagle (Java jar)

Beagle ships as a `.jar` distributed separately from conda. Point the tool
at it with the `BEAGLE_JAR` environment variable:

```bash
export BEAGLE_JAR=/path/to/beagle.jar
```

Add this to your shell profile (`~/.bashrc`, `~/.zshrc`) to make it
persistent. Verify:

```bash
test -f "$BEAGLE_JAR" && echo "OK" || echo "MISSING"
java -version
```

### 7. Build

```bash
git clone https://github.com/LoaiEletr/vcf2fasta-rust
cd vcf2fasta-rust

# CPU-only build (default)
cargo build --release

# GPU-enabled build (requires NVIDIA CUDA Toolkit ≥ 12.0)
cargo build --release --features cuda
```

The binary is at `target/release/vcf2fasta`

### 8. Verify

```bash
target/release/vcf2fasta --help
```

### 9. Run

```bash
# CPU, 4 threads
target/release/vcf2fasta \
    --device cpu \
    --threads 4 \
    --reference ref.fa \
    --prefix out_ \
    input.vcf.gz

# GPU
target/release/vcf2fasta \
    --device gpu \
    --reference ref.fa \
    --prefix out_ \
    input.vcf.gz

# Multi-GPU (both devices) with a VRAM cap
target/release/vcf2fasta \
    --device gpu \
    --gpu-devices 0,1 \
    --max-vram 16G \
    --reference ref.fa \
    --prefix out_ \
    input.vcf.gz

# Auto: the scheduler picks CPU or GPU from workload size and hardware.
# This is the default, so --device can be omitted.
target/release/vcf2fasta \
    --threads 4 \
    --reference ref.fa \
    --prefix out_ \
    input.vcf.gz

# Automatic phasing (supply BAMs if needed)
target/release/vcf2fasta \
    --reference ref.fa \
    --bam-dir /path/to/bams \
    --prefix out_ \
    input.vcf.gz
```

`--device auto` is the default. It chooses GPU only when a CUDA device
is visible to the process **and** the workload is large enough to
amortise the transfer overhead. On small inputs, or on
machines without a GPU, or on binaries built without `--features cuda`,
it falls back to CPU.

`--gpu-devices` and `--max-vram` are only consulted when the plan
resolves to GPU; they are ignored on the CPU path.

Output files:

- **Per-contig (default):** `<prefix><sample>_<contig>:<hap>.fa`.
- **Merged (`--merged-output`):** `<prefix><sample>_<hap>.fa`, a
  multi-record FASTA containing every contig in VCF-header order.

---

## 🧪 Testing

The test suite lives in `tests/` and consists of two integration suites
plus the unit tests embedded in each `src/*.rs` module:

- `tests/pipeline_integration.rs` – discovery and work-unit state classification.
- `tests/policy_integration.rs` – pipeline-level VCF2FASTA policies (unsorted VCF rejection, duplicate/overlap handling). Requires `bgzip` and `tabix` on `PATH`; otherwise those tests print a skip notice and return early.

Some tests are **CUDA-gated** with `#[cfg(feature = "cuda")]` (in
`src/gpu.rs`, `src/cuda_pipeline.rs`, and the CUDA branch of `src/lib.rs`).
They are compiled and run only when the `cuda` feature is enabled.

### CPU-only build (default)

```bash
cargo test --release
```

### GPU build

```bash
cargo test --release --features cuda
```

### Notes

- `policy_integration.rs` skips gracefully if `bgzip` / `tabix` are missing — you will see a `skipping ...` line on stderr, not a failure.
- There is no separate `--test regression` suite. Regression checking against the official vcflib `vcf2fasta` is done through the benchmark scripts (see [Benchmarking Guide](docs/BENCHMARKING.md)).

---

## 📊 Benchmarking

See **[Benchmarking Guide](docs/BENCHMARKING.md)** for the full workflow.

The benchmark scripts have a **dependency chain**, not a flat list:

```
download_data.sh
    └── fetch GRCh38 + Platinum + 1000G chr22
            │
            ├── run_benchmark.sh              (reads Platinum + GRCh38)
            ├── scaling_benchmark.sh          (reads 1000G chr22 + GRCh38)
            │
            └── generate_synthetic_datasets.sh
                    └── datasets/synthetic_unphased/<DATASET>/
                            └── unphased_scaling_benchmark.sh
```

`unphased_scaling_benchmark.sh` **requires** the synthetic dataset it is
asked to run on to have been generated first. Running
`DATASET=medium ./unphased_scaling_benchmark.sh` on a machine where
`medium` was never built will fail immediately.

Minimum viable sequence:

```bash
cd tests/benchmark/scripts

# 1. One-time setup: fetch reference and public datasets.
./download_data.sh

# 2. Build the synthetic datasets you intend to use.
#    `tiny` is the smallest; add more if you want a scaling curve.
DATASET=tiny   ./generate_synthetic_datasets.sh
DATASET=medium ./generate_synthetic_datasets.sh

# 3. Benchmarks that read the *real* datasets (no synthetic step needed).
MODE=both USE_VCFLIB=true ./run_benchmark.sh    # per-sample + vcflib cross-check
MODE=both ./scaling_benchmark.sh                # sample-count scaling on 1000G chr22

# 4. Benchmark that reads a *synthetic* dataset (must be built in step 2).
#    This runs the integrated pipeline (phasing → vcf2fasta) end-to-end,
#    so it is slower than the pure-vcf2fasta benchmarks above.
DATASET=medium MODE=both ./unphased_scaling_benchmark.sh
```

To run the integrated pipeline benchmark across all five synthetic
datasets (which is what the published Figure 4 shows), build all of them
first:

```bash
for D in tiny small medium larger ploidy_edge; do
    DATASET=$D ./generate_synthetic_datasets.sh
done

for D in tiny small medium larger ploidy_edge; do
    DATASET=$D MODE=both ./unphased_scaling_benchmark.sh
done
```

Result spreadsheets are committed under
[`tests/benchmark/results/`](tests/benchmark/results/) — one `.xlsx` per
benchmark script, each with per-run sheets (`RunTime1` … `RunTime4`) and a
`Combined` sheet containing mean ± SD. When `USE_VCFLIB=true`, `run_benchmark.sh` also writes
[`tests/benchmark/results/compare_summary.tsv`](tests/benchmark/results/compare_summary.tsv)
— a per-haplotype table of Rust-vs-vcflib length differences, attributed
to the Rust tool's overlap policy.

Figures are regenerated from those spreadsheets by:

```bash
cd tests/benchmark/scripts
Rscript plot_benchmark_results.R
```

The script defaults `RESULTS_DIR` to `../../results` and `FIGURES_DIR` to
`../../../docs/figures`, so no environment variables are needed in the
normal case.

Published figures live under [`docs/figures/`](docs/figures/) and are
embedded in [Benchmark Results](docs/BENCHMARK_RESULTS.md).

---

## 📄 License

This project is released under the MIT License – see the [LICENSE](LICENSE)
file for details.

---