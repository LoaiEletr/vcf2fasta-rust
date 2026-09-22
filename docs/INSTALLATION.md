# Installation Guide for vcf2fasta-rust

Full setup for running vcf2fasta-rust, the benchmark suite, and the
figure-regeneration script.

---

## 1. Rust toolchain

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
```

Add Cargo to `PATH` for the current session if needed:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
```

Verify:

```bash
cargo --version
rustc --version
```

---

## 2. System build dependencies

The project links against `libclang` (via `rust-htslib` / `bindgen`).

### Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    pkg-config \
    libclang-dev
```

### macOS

```bash
brew install llvm pkg-config
export LIBCLANG_PATH="$(brew --prefix llvm)/lib"
```

### For GPU builds

Install the NVIDIA CUDA Toolkit (≥ 12.0). Confirm `nvcc` is on `PATH`:

```bash
nvcc --version
```

Set the CUDA library path if needed:

```bash
export LD_LIBRARY_PATH=/usr/local/cuda-12.8/targets/x86_64-linux/lib:/usr/local/cuda-12.8/lib64:$LD_LIBRARY_PATH
```

Multi-GPU requires no additional setup beyond a working CUDA
installation. All visible devices are enumerated at startup and
reported in the `[SCHEDULER]` log block; the `--gpu-devices` flag
selects which of them to use.

---

## 3. Conda / mamba

### 3.1. Install Miniforge

If you do not already have conda or mamba, install **Miniforge**, which
ships with `mamba` preconfigured for `conda-forge`.

**Linux x86_64:**

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
```

**macOS** — replace the URL with the appropriate installer from
[the Miniforge releases page](https://github.com/conda-forge/miniforge/releases):

- Apple Silicon: `Miniforge3-MacOSX-arm64.sh`
- Intel: `Miniforge3-MacOSX-x86_64.sh`

### 3.2. Initialise the shell

```bash
conda init bash    # or zsh, fish, etc.
exec "$SHELL" -l   # reload the shell
```

This is a one-time step. After it, opening a new terminal shows
`(base)` in the prompt and `conda activate` works.

### 3.3. Create a dedicated environment

```bash
mamba create -n vcf2fasta python=3.12 -y
mamba activate vcf2fasta
```

Do **not** install into `base`. A dedicated environment isolates this
project's tool versions from everything else on the machine.

### 3.4. Install bioinformatics tools

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

To reactivate this environment later:

```bash
conda activate vcf2fasta
```

To deactivate:

```bash
conda deactivate
```

---

## 4. Beagle setup

Beagle ships as a Java `.jar`. Point the tool at it with the `BEAGLE_JAR`
environment variable:

```bash
export BEAGLE_JAR=/usr/local/share/beagle-5.5_27Feb25.75f-0/beagle.jar
```

Add this to your shell profile (`~/.bashrc`, `~/.zshrc`) to make it
persistent across sessions.

Confirm:

```bash
echo "BEAGLE_JAR = $BEAGLE_JAR"
test -f "$BEAGLE_JAR" && echo "OK" || echo "MISSING"
```

Also make sure `java` is on `PATH`:

```bash
java -version
```

If `java` is missing, install OpenJDK via conda:

```bash
mamba install -y -c conda-forge openjdk
```

---

## 5. R and R packages (for figure regeneration)

The benchmark figures are produced by
`tests/benchmark/scripts/plot_benchmark_results.R`. Install R (≥ 4.3)
and the required packages once.

### Ubuntu / Debian

```bash
sudo apt-get install -y r-base
```

### macOS

```bash
brew install r
```

### Install the R packages

In an R session:

```r
install.packages(c(
  "readxl",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "patchwork",
  "scales",
  "forcats"
))
```

These are only needed if you plan to regenerate figures. The tool itself
does not depend on R.

---

## 6. Build the binary

```bash
git clone https://github.com/LoaiEletr/vcf2fasta-rust
cd vcf2fasta-rust

# CPU-only
cargo build --release

# GPU (CUDA)
cargo build --release --features cuda
```

The binary is placed at `target/release/vcf2fasta`

---

## 7. Verify installation

```bash
target/release/vcf2fasta --help
```

If CUDA was enabled:

```bash
target/release/vcf2fasta --device gpu --help
```

If R was installed:

```bash
Rscript -e 'library(ggplot2); cat("R plotting OK\n")'
```

---

## 8. Troubleshooting

- **`mamba: command not found`** — conda/mamba is not on `PATH`. Run
  `source "$HOME/miniforge3/etc/profile.d/conda.sh"` (Linux/macOS).
- **`libclang` not found during `cargo build`** — the system build
  dependencies (Section 2) are missing. On macOS also set
  `LIBCLANG_PATH`.
- **`--device gpu` fails with "no compatible GPU"** — either the binary
  was built without `--features cuda`, or no NVIDIA driver is loaded.
  Check `nvidia-smi` and rebuild with `--features cuda`.
- **`--gpu-devices` fails with "not visible to this process"** — one or
  more requested CUDA indices are not enumerated by the driver. Check
  the `[SCHEDULER] gpu[N]:` lines in the log for the indices that are
  actually visible, and adjust the list.
- **Beagle jar not found** — set `BEAGLE_JAR` ([Section 4](#4-beagle-setup)). The tool
  looks for it via the environment variable only; it does not search
  `PATH`.

---