#!/usr/bin/env bash
# ============================================================================
# Unphased Dataset Scaling Benchmark
#
# Purpose:
#   Benchmark the integrated unphased VCF -> phasing -> vcf2fasta pipeline
#   across the predefined synthetic datasets.
#
# Unlike scaling_benchmark.sh:
#   - NO sample-count selection
#   - NO bcftools sample extraction
#   - The selected dataset is used in its entirety
#   - Dataset is selected with DATASET=
#
# Supported datasets:
#   tiny
#   small
#   medium
#   larger
#   ploidy_edge
#
# MODE (default: auto):
#   auto   → GPU present ? both : cpu
#   cpu    → CPU sweep only
#   gpu    → GPU sweep only (falls back to cpu if no GPU)
#   both   → GPU sweep AND CPU sweep
#
# GPU sweep:
#   GPU_DEVICES holds a comma-separated pool of CUDA device indices, e.g.
#   "0" or "0,1,2,3". The script sweeps all cumulative prefixes of this
#   pool, one run per (prefix, host-thread) pair:
#
#       --gpu-devices 0
#       --gpu-devices 0,1
#       --gpu-devices 0,1,2
#       --gpu-devices 0,1,2,3
#
#   Each prefix is combined with every value in GPU_THREADS_LIST.
#
# Example:
#   DATASET=tiny     MODE=cpu  ./unphased_scaling_benchmark.sh
#   DATASET=medium   MODE=gpu  GPU_DEVICES=0,1 GPU_THREADS_LIST="1 2" ./unphased_scaling_benchmark.sh
#   DATASET=ploidy_edge MODE=both GPU_DEVICES=0,1 ./unphased_scaling_benchmark.sh
#
# Beagle:
#   --beagle-genetic-map <DIR> is passed automatically. Expected location:
#       datasets/genetic_map/
# ============================================================================

set -uo pipefail

# ============================================================================
# CUDA
# ============================================================================

if [ -d /usr/local/cuda-12.8/targets/x86_64-linux/lib ]; then
    export LD_LIBRARY_PATH="/usr/local/cuda-12.8/targets/x86_64-linux/lib:/usr/local/cuda-12.8/lib64:/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
fi

# ============================================================================
# PATHS
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(dirname "$SCRIPT_DIR")"

RUST_BIN="${BENCH_DIR}/../../target/release/vcf2fasta"

DATASETS_DIR="${BENCH_DIR}/datasets/synthetic_unphased"

BEAGLE_GENETIC_MAP_DIR="${BENCH_DIR}/datasets/genetic_map"

RESULTS_DIR="${BENCH_DIR}/unphased_scaling_results"
LOG_FILE="${RESULTS_DIR}/unphased_scaling_benchmark.log"

mkdir -p "$RESULTS_DIR"

# ============================================================================
# USER-CONFIGURABLE PARAMETERS
# ============================================================================

# Select exactly ONE synthetic dataset.
# Supported: tiny | small | medium | larger | ploidy_edge
DATASET="${DATASET:-tiny}"

CPU_THREADS_LIST="${CPU_THREADS_LIST:-1 2 4 8}"

# Pool of CUDA device indices. Cumulative prefixes are swept.
#   GPU_DEVICES=0             → --gpu-devices 0
#   GPU_DEVICES=0,1           → --gpu-devices 0 ; --gpu-devices 0,1
#   GPU_DEVICES=0,1,2,3       → 0 ; 0,1 ; 0,1,2 ; 0,1,2,3
#   GPU_DEVICES=2,3           → 2 ; 2,3
GPU_DEVICES="${GPU_DEVICES:-0,1}"

GPU_THREADS_LIST="${GPU_THREADS_LIST:-1 2 4}"

REQUESTED_MODE="${MODE:-auto}"   # auto | cpu | gpu | both

# ============================================================================
# DATASET CONFIGURATION
# ============================================================================

DATASET_DIR="${DATASETS_DIR}/${DATASET}"

INPUT_VCF="${DATASET_DIR}/vcf/input_mixed_phase.vcf.gz"
INPUT_VCF_INDEX="${INPUT_VCF}.tbi"

TRUTH_VCF="${DATASET_DIR}/vcf/truth_phased.vcf.gz"
TRUTH_VCF_INDEX="${TRUTH_VCF}.tbi"

REFERENCE="${DATASET_DIR}/reference/synthetic_reference.fa"
REFERENCE_FAI="${REFERENCE}.fai"

BAM_DIR="${DATASET_DIR}/bam"

METADATA="${DATASET_DIR}/metadata.json"

# ============================================================================
# TIME FORMATTER
# ============================================================================

format_readable_time() {
    local total_sec="$1"

    local hours minutes seconds millis remainder

    hours=$(echo "$total_sec / 3600" | bc)
    remainder=$(echo "$total_sec - ($hours * 3600)" | bc)
    minutes=$(echo "$remainder / 60" | bc)
    remainder=$(echo "$remainder - ($minutes * 60)" | bc)
    seconds=$(echo "$remainder / 1" | bc)
    millis=$(echo "($remainder - $seconds) * 1000" | bc | awk '{print int($1+0.5)}')

    local out=""
    [ "$hours" -gt 0 ] && out="${hours}h "
    { [ "$minutes" -gt 0 ] || [ "$hours" -gt 0 ]; } && out="${out}${minutes}m "
    out="${out}${seconds}s ${millis}ms"
    echo "$out"
}

# ============================================================================
# GPU PREFIX SWEEP HELPERS
# ============================================================================
# gpu_device_prefixes "0,1,2,3" prints:
#     0
#     0,1
#     0,1,2
#     0,1,2,3
# Each line is passed verbatim to --gpu-devices.

gpu_device_prefixes() {
    local list="$1"
    local -a devs
    IFS=',' read -ra devs <<< "$list"

    local prefix="" d
    for d in "${devs[@]}"; do
        d="${d//[[:space:]]/}"
        [ -z "$d" ] && continue
        if [ -z "$prefix" ]; then
            prefix="$d"
        else
            prefix="${prefix},${d}"
        fi
        echo "$prefix"
    done
}

# Safe label for filenames / column names: replace commas with underscores.
gpu_label_safe() {
    echo "$1" | tr ',' '_'
}

# ============================================================================
# DEPENDENCY CHECK
# ============================================================================

echo
echo "=========================================================================="
echo "🔧 Checking benchmark dependencies..."
echo "=========================================================================="

if [ ! -x "$RUST_BIN" ]; then
    echo "❌ Rust binary not found:"
    echo "   $RUST_BIN"
    echo
    echo "Build with:"
    echo "   cargo build --release --features cuda"
    exit 1
fi

if ! command -v bc >/dev/null 2>&1; then
    echo "❌ bc is not installed."
    exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "❌ bcftools is not installed."
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "❌ samtools is not installed."
    exit 1
fi

if ! command -v whatshap >/dev/null 2>&1; then
    echo "⚠️  whatshap was not found in PATH."
    echo "   This may be required by the integrated phasing pipeline."
fi

echo "✅ Rust binary found."
echo "✅ Required benchmark dependencies checked."

# ============================================================================
# VALIDATE BEAGLE GENETIC MAP DIRECTORY
# ============================================================================

echo
echo "=========================================================================="
echo "🧬 Checking Beagle genetic-map directory"
echo "=========================================================================="

if [ ! -d "$BEAGLE_GENETIC_MAP_DIR" ]; then
    echo "❌ Beagle genetic-map directory not found:"
    echo "   $BEAGLE_GENETIC_MAP_DIR"
    exit 1
fi

if ! find "$BEAGLE_GENETIC_MAP_DIR" -type f -print -quit | grep -q .; then
    echo "❌ Beagle genetic-map directory is empty:"
    echo "   $BEAGLE_GENETIC_MAP_DIR"
    exit 1
fi

echo "✅ Beagle genetic-map directory found:"
echo "   $BEAGLE_GENETIC_MAP_DIR"

echo
echo "Genetic-map files:"
find "$BEAGLE_GENETIC_MAP_DIR" -maxdepth 1 -type f -printf '   %f\n' | sort

# ============================================================================
# VALID DATASET CHECK
# ============================================================================

case "$DATASET" in
    tiny|small|medium|larger|ploidy_edge) ;;
    *)
        echo
        echo "❌ Unknown dataset: $DATASET"
        echo
        echo "Supported datasets:"
        echo "   tiny"
        echo "   small"
        echo "   medium"
        echo "   larger"
        echo "   ploidy_edge"
        exit 1
        ;;
esac

# ============================================================================
# CHECK DATASET
# ============================================================================

echo
echo "=========================================================================="
echo "📦 Checking selected dataset"
echo "=========================================================================="

echo "Dataset:"
echo "   $DATASET"

echo "Dataset directory:"
echo "   $DATASET_DIR"

if [ ! -d "$DATASET_DIR" ]; then
    echo "❌ Dataset directory does not exist:"
    echo "   $DATASET_DIR"
    echo
    echo "Generate it first with:"
    echo
    echo "   DATASET=$DATASET ./generate_synthetic_datasets.sh"
    exit 1
fi

if [ ! -s "$INPUT_VCF" ]; then
    echo "❌ Input VCF not found or empty:"
    echo "   $INPUT_VCF"
    exit 1
fi

if [ ! -s "$INPUT_VCF_INDEX" ]; then
    echo "❌ Input VCF index not found:"
    echo "   $INPUT_VCF_INDEX"
    exit 1
fi

if [ ! -s "$TRUTH_VCF" ]; then
    echo "❌ Truth VCF not found or empty:"
    echo "   $TRUTH_VCF"
    exit 1
fi

if [ ! -s "$TRUTH_VCF_INDEX" ]; then
    echo "⚠️  Truth VCF index not found:"
    echo "   $TRUTH_VCF_INDEX"
    echo "   Continuing because the benchmark does not directly require it."
fi

if [ ! -s "$REFERENCE" ]; then
    echo "❌ Synthetic reference not found:"
    echo "   $REFERENCE"
    exit 1
fi

if [ ! -s "$REFERENCE_FAI" ]; then
    echo "❌ Reference index not found:"
    echo "   $REFERENCE_FAI"
    exit 1
fi

if [ ! -d "$BAM_DIR" ]; then
    echo "❌ BAM directory not found:"
    echo "   $BAM_DIR"
    exit 1
fi

if [ ! -s "$METADATA" ]; then
    echo "⚠️  Metadata file not found:"
    echo "   $METADATA"
fi

echo "✅ Dataset files found."

# ============================================================================
# READ SAMPLE INFORMATION
# ============================================================================

echo
echo "=========================================================================="
echo "📊 Reading dataset information"
echo "=========================================================================="

if ! bcftools view -h "$INPUT_VCF" >/dev/null 2>"${RESULTS_DIR}/input_vcf_check.stderr"; then
    echo "❌ Cannot read input VCF."
    [ -s "${RESULTS_DIR}/input_vcf_check.stderr" ] && cat "${RESULTS_DIR}/input_vcf_check.stderr"
    exit 1
fi

SAMPLE_NAMES_FILE="${RESULTS_DIR}/${DATASET}_samples.txt"

if ! bcftools query -l "$INPUT_VCF" > "$SAMPLE_NAMES_FILE" 2>"${RESULTS_DIR}/${DATASET}_sample_query.stderr"; then
    echo "❌ Failed to read sample names from input VCF."
    [ -s "${RESULTS_DIR}/${DATASET}_sample_query.stderr" ] && cat "${RESULTS_DIR}/${DATASET}_sample_query.stderr"
    exit 1
fi

TOTAL_SAMPLES=$(wc -l < "$SAMPLE_NAMES_FILE")
TOTAL_SAMPLES=$(echo "$TOTAL_SAMPLES" | tr -d '[:space:]')

if [ "$TOTAL_SAMPLES" -le 0 ]; then
    echo "❌ Input VCF contains no samples."
    exit 1
fi

echo "Samples:"
echo "   $TOTAL_SAMPLES"

# ============================================================================
# DISCOVER BAM FILES
# ============================================================================

BAM_COUNT=$(find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | wc -l)
BAM_COUNT=$(echo "$BAM_COUNT" | tr -d '[:space:]')

echo "BAM files:"
echo "   $BAM_COUNT"

if [ "$BAM_COUNT" -ne "$TOTAL_SAMPLES" ]; then
    echo
    echo "❌ BAM/sample count mismatch."
    echo "   VCF samples: $TOTAL_SAMPLES"
    echo "   BAM files:   $BAM_COUNT"
    echo
    echo "BAM files found:"
    find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" -print | sort
    exit 1
fi

# ============================================================================
# VALIDATE BAM FILES
# ============================================================================

echo
echo "=========================================================================="
echo "🔍 Validating BAM files..."
echo "=========================================================================="

while IFS= read -r SAMPLE; do
    BAM="${BAM_DIR}/${SAMPLE}.sorted.bam"
    BAI="${BAM}.bai"

    if [ ! -s "$BAM" ]; then
        echo "❌ BAM missing for sample:"
        echo "   $SAMPLE"
        exit 1
    fi

    if [ ! -s "$BAI" ]; then
        echo "❌ BAM index missing for sample:"
        echo "   $SAMPLE"
        exit 1
    fi
done < "$SAMPLE_NAMES_FILE"

echo "✅ All BAM files and indexes found."

# ============================================================================
# DISCOVER CHROMOSOMES
# ============================================================================

CHROMS_FILE="${RESULTS_DIR}/${DATASET}_chromosomes.txt"

if ! bcftools query -f '%CHROM\n' "$INPUT_VCF" 2>"${RESULTS_DIR}/${DATASET}_chrom_query.stderr" \
    | sort -u > "$CHROMS_FILE"; then
    echo "❌ Failed to read chromosomes from input VCF."
    exit 1
fi

CHROM_COUNT=$(wc -l < "$CHROMS_FILE")
CHROM_COUNT=$(echo "$CHROM_COUNT" | tr -d '[:space:]')

echo "Chromosomes:"
sed 's/^/   /' "$CHROMS_FILE"

echo "Chromosome count: $CHROM_COUNT"

# ============================================================================
# GPU DETECTION + MODE RESOLUTION
# ============================================================================

HAS_GPU=false
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    HAS_GPU=true
fi

case "$REQUESTED_MODE" in
    auto)
        if [ "$HAS_GPU" = true ]; then
            RUN_CPU=true; RUN_GPU=true
            echo "ℹ️  MODE=auto and NVIDIA GPU detected → running GPU and CPU sweeps."
        else
            RUN_CPU=true; RUN_GPU=false
            echo "ℹ️  MODE=auto and no GPU detected → running CPU-only sweep."
        fi
        ;;
    cpu)
        RUN_CPU=true; RUN_GPU=false
        echo "ℹ️  MODE=cpu → running CPU-only sweep."
        ;;
    gpu)
        if [ "$HAS_GPU" = true ]; then
            RUN_CPU=false; RUN_GPU=true
            echo "ℹ️  MODE=gpu and NVIDIA GPU detected → running GPU-only sweep."
        else
            RUN_CPU=true; RUN_GPU=false
            echo "⚠️  MODE=gpu requested but no NVIDIA GPU detected → falling back to CPU-only."
        fi
        ;;
    both)
        if [ "$HAS_GPU" = true ]; then
            RUN_CPU=true; RUN_GPU=true
            echo "ℹ️  MODE=both and NVIDIA GPU detected → running GPU and CPU sweeps."
        else
            RUN_CPU=true; RUN_GPU=false
            echo "⚠️  MODE=both requested but no NVIDIA GPU detected → running CPU-only."
        fi
        ;;
    *)
        echo "❌ Invalid MODE='$REQUESTED_MODE' (expected: auto, cpu, gpu, both)"
        exit 1
        ;;
esac

# Derive EXEC_MODE for the log column (matches scaling_benchmark.sh)
if [ "$RUN_CPU" = true ] && [ "$RUN_GPU" = true ]; then
    EXEC_MODE="both"
elif [ "$RUN_GPU" = true ]; then
    EXEC_MODE="gpu"
else
    EXEC_MODE="cpu"
fi

# ============================================================================
# CONFIGURATION SUMMARY
# ============================================================================

echo
echo "=========================================================================="
echo "⚙️  Benchmark Configuration"
echo "=========================================================================="

echo "Dataset:          $DATASET"
echo "Input VCF:        $INPUT_VCF"
echo "Truth VCF:        $TRUTH_VCF"
echo "Reference:        $REFERENCE"
echo "BAM directory:    $BAM_DIR"
echo "Beagle map dir:   $BEAGLE_GENETIC_MAP_DIR"
echo "Samples:          $TOTAL_SAMPLES"
echo "Chromosomes:      $CHROM_COUNT"
echo "Execution mode:   $EXEC_MODE   (RUN_CPU=$RUN_CPU, RUN_GPU=$RUN_GPU)"

if [ "$RUN_CPU" = true ]; then
    echo "CPU threads:      $CPU_THREADS_LIST"
fi

if [ "$RUN_GPU" = true ]; then
    echo "GPU device pool:  $GPU_DEVICES"
    echo "GPU host threads: $GPU_THREADS_LIST"
    echo "GPU runs to be performed:"
    while IFS= read -r PREFIX; do
        NDEVS=$(awk -F, '{print NF}' <<< "$PREFIX")
        echo "  ${NDEVS} GPU(s):  --gpu-devices ${PREFIX}   × host threads {$(echo "$GPU_THREADS_LIST" | tr ' ' ',')}"
    done < <(gpu_device_prefixes "$GPU_DEVICES")
fi

echo "=========================================================================="

# ============================================================================
# CREATE LOG HEADER
# ============================================================================

HEADER="Dataset\tSamples\tChromosomes\tMode"

if [ "$RUN_GPU" = true ]; then
    while IFS= read -r PREFIX; do
        SAFE=$(gpu_label_safe "$PREFIX")
        for GT in $GPU_THREADS_LIST; do
            HEADER="${HEADER}\tgpu_dev${SAFE}_t${GT}_time"
        done
    done < <(gpu_device_prefixes "$GPU_DEVICES")
fi

if [ "$RUN_CPU" = true ]; then
    for T in $CPU_THREADS_LIST; do
        HEADER="${HEADER}\tcpu_t${T}_time"
    done
fi

HEADER="${HEADER}\tNotes"

printf "%b\n" "$HEADER" > "$LOG_FILE"

# ============================================================================
# RESULT DIRECTORY
# ============================================================================

DATASET_RESULTS_DIR="${RESULTS_DIR}/${DATASET}"

mkdir -p \
    "${DATASET_RESULTS_DIR}/errors" \
    "${DATASET_RESULTS_DIR}/outputs"

# ============================================================================
# TIMING STORAGE
# ============================================================================

declare -A RUN_TIMES
NOTE_TEXT="Success"

# ============================================================================
# GPU BENCHMARK (prefix sweep × host-thread sweep)
# ============================================================================

if [ "$RUN_GPU" = true ]; then
    while IFS= read -r DEV_LIST; do
        SAFE=$(gpu_label_safe "$DEV_LIST")
        NDEVS=$(awk -F, '{print NF}' <<< "$DEV_LIST")

        for GT in $GPU_THREADS_LIST; do
            CONFIG_DIR="${DATASET_RESULTS_DIR}/gpu_dev${SAFE}_t${GT}"
            mkdir -p "$CONFIG_DIR"

            GPU_STDOUT="${DATASET_RESULTS_DIR}/errors/gpu_dev${SAFE}_t${GT}.stdout"
            GPU_STDERR="${DATASET_RESULTS_DIR}/errors/gpu_dev${SAFE}_t${GT}.stderr"

            echo
            echo "=========================================================================="
            echo "🚀 GPU Benchmark"
            echo "   Dataset:       $DATASET"
            echo "   --gpu-devices: $DEV_LIST  (${NDEVS} device(s))"
            echo "   --threads:     $GT"
            echo "=========================================================================="

            START=$(date +%s.%N)
            set +e

            "$RUST_BIN" \
                --device gpu \
                --gpu-devices "$DEV_LIST" \
                --threads "$GT" \
                --reference "$REFERENCE" \
                "$INPUT_VCF" \
                --bam-dir "$BAM_DIR" \
                --beagle-genetic-map "$BEAGLE_GENETIC_MAP_DIR" \
                --no-validate-ref \
                --prefix "${CONFIG_DIR}/result_" \
                >"$GPU_STDOUT" \
                2>"$GPU_STDERR"

            EXIT=$?
            set -e

            END=$(date +%s.%N)
            ELAPSED=$(echo "$END - $START" | bc)
            RUN_TIMES["gpu_dev${SAFE}_t${GT}"]="$ELAPSED"
            READABLE_T=$(format_readable_time "$ELAPSED")

            if [ "$EXIT" -eq 0 ]; then
                echo "✅ GPU dev${SAFE} / t${GT} completed successfully."
                echo "   Time: $READABLE_T"
            else
                ERR_MSG=$(head -n 1 "$GPU_STDERR" 2>/dev/null || true)
                echo "❌ GPU dev${SAFE} / t${GT} failed."
                echo "   Exit code: $EXIT"
                echo "   Error: ${ERR_MSG:-unknown}"

                if [ "$NOTE_TEXT" = "Success" ]; then
                    NOTE_TEXT="GPU dev${SAFE}_t${GT} Error"
                else
                    NOTE_TEXT="${NOTE_TEXT}; GPU dev${SAFE}_t${GT} Error"
                fi
            fi
        done
    done < <(gpu_device_prefixes "$GPU_DEVICES")
fi

# ============================================================================
# CPU BENCHMARK
# ============================================================================

if [ "$RUN_CPU" = true ]; then
    for T in $CPU_THREADS_LIST; do
        CONFIG_DIR="${DATASET_RESULTS_DIR}/cpu_t${T}"
        mkdir -p "$CONFIG_DIR"

        CPU_STDOUT="${DATASET_RESULTS_DIR}/errors/cpu_t${T}.stdout"
        CPU_STDERR="${DATASET_RESULTS_DIR}/errors/cpu_t${T}.stderr"

        echo
        echo "=========================================================================="
        echo "🖥️  CPU Benchmark"
        echo "   Dataset: $DATASET"
        echo "   Threads: $T"
        echo "=========================================================================="

        START=$(date +%s.%N)
        set +e

        "$RUST_BIN" \
            --device cpu \
            --threads "$T" \
            --reference "$REFERENCE" \
            "$INPUT_VCF" \
            --bam-dir "$BAM_DIR" \
            --beagle-genetic-map "$BEAGLE_GENETIC_MAP_DIR" \
            --no-validate-ref \
            --prefix "${CONFIG_DIR}/result_" \
            >"$CPU_STDOUT" \
            2>"$CPU_STDERR"

        EXIT=$?
        set -e

        END=$(date +%s.%N)
        ELAPSED=$(echo "$END - $START" | bc)
        RUN_TIMES["cpu_t${T}"]="$ELAPSED"
        READABLE_T=$(format_readable_time "$ELAPSED")

        if [ "$EXIT" -eq 0 ]; then
            echo "✅ CPU t${T} completed successfully."
            echo "   Time: $READABLE_T"
        else
            ERR_MSG=$(head -n 1 "$CPU_STDERR" 2>/dev/null || true)
            echo "❌ CPU t${T} failed."
            echo "   Exit code: $EXIT"
            echo "   Error: ${ERR_MSG:-unknown}"

            if [ "$NOTE_TEXT" = "Success" ]; then
                NOTE_TEXT="CPU t${T} Error"
            else
                NOTE_TEXT="${NOTE_TEXT}; CPU t${T} Error"
            fi
        fi
    done
fi

# ============================================================================
# WRITE RESULTS
# ============================================================================

ROW_DATA="${DATASET}\t${TOTAL_SAMPLES}\t${CHROM_COUNT}\t${EXEC_MODE}"

if [ "$RUN_GPU" = true ]; then
    while IFS= read -r DEV_LIST; do
        SAFE=$(gpu_label_safe "$DEV_LIST")
        for GT in $GPU_THREADS_LIST; do
            G_TIME="${RUN_TIMES["gpu_dev${SAFE}_t${GT}"]:-0}"
            ROW_DATA="${ROW_DATA}\t$(format_readable_time "$G_TIME")"
        done
    done < <(gpu_device_prefixes "$GPU_DEVICES")
fi

if [ "$RUN_CPU" = true ]; then
    for T in $CPU_THREADS_LIST; do
        T_TIME="${RUN_TIMES["cpu_t${T}"]:-0}"
        ROW_DATA="${ROW_DATA}\t$(format_readable_time "$T_TIME")"
    done
fi

ROW_DATA="${ROW_DATA}\t${NOTE_TEXT}"

printf "%b\n" "$ROW_DATA" >> "$LOG_FILE"

# ============================================================================
# FINAL DATASET SUMMARY
# ============================================================================

echo
echo "=========================================================================="
echo "📊 Dataset Benchmark Complete"
echo "=========================================================================="

echo "Dataset:     $DATASET"
echo "Samples:     $TOTAL_SAMPLES"
echo "Chromosomes: $CHROM_COUNT"
echo "Beagle map:  $BEAGLE_GENETIC_MAP_DIR"

if [ "$RUN_GPU" = true ]; then
    while IFS= read -r DEV_LIST; do
        SAFE=$(gpu_label_safe "$DEV_LIST")
        NDEVS=$(awk -F, '{print NF}' <<< "$DEV_LIST")
        for GT in $GPU_THREADS_LIST; do
            G_TIME="${RUN_TIMES["gpu_dev${SAFE}_t${GT}"]:-0}"
            echo "GPU ${NDEVS}g (dev=${DEV_LIST}) / t${GT}: $(format_readable_time "$G_TIME")"
        done
    done < <(gpu_device_prefixes "$GPU_DEVICES")
fi

if [ "$RUN_CPU" = true ]; then
    for T in $CPU_THREADS_LIST; do
        T_TIME="${RUN_TIMES["cpu_t${T}"]:-0}"
        echo "CPU t${T}:     $(format_readable_time "$T_TIME")"
    done
fi

echo
echo "Notes: $NOTE_TEXT"

echo
echo "Results:"
echo "   $LOG_FILE"

echo
echo "Outputs:"
echo "   $DATASET_RESULTS_DIR"

echo
echo "Beagle genetic-map directory passed to Rust:"
echo "   $BEAGLE_GENETIC_MAP_DIR"

echo "=========================================================================="
echo "✅ Unphased scaling benchmark complete."
echo "=========================================================================="