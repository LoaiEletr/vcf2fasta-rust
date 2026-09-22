#!/usr/bin/env bash
# ============================================================================
# Scaling Benchmark Script: Combinatorial GPU & CPU Scaling
#
# Env vars:
#   MODE              auto | cpu | gpu | both   (default: auto)
#                       auto  → GPU present ? both : cpu
#                       cpu   → CPU only
#                       gpu   → GPU only (falls back to cpu if no GPU)
#                       both  → run GPU sweep AND CPU sweep
#
#   OUTPUT_BASE       directory under which FASTA outputs are written.
#                     Default: $RESULTS_DIR. Set to /tmp or /dev/shm to
#                     benchmark against a RAM-backed filesystem.
#
#   GPU_DEVICES       comma-separated pool of CUDA device indices, e.g.
#                     "0" or "0,1,2,3". The script sweeps all *cumulative
#                     prefixes* of this pool:
#                       0
#                       0,1
#                       0,1,2
#                       0,1,2,3
#                     Each prefix is combined with every value in
#                     GPU_THREADS_LIST. The literal string passed to the
#                     tool's --gpu-devices matches the prefix exactly.
#
#   GPU_THREADS_LIST  space-separated GPU host-thread counts (default "1 2 4")
#
#   CHROMS            space-separated contigs (default: chr22)
#   SAMPLES_LIST      space-separated sample tiers (default: "10 50 100")
#   CPU_THREADS_LIST  space-separated CPU thread counts (default: "1 2 4 8")
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

RESULTS_DIR="${BENCH_DIR}/scaling_results"
INPUTS_DIR="${RESULTS_DIR}/inputs"
LOG_FILE="${RESULTS_DIR}/scaling_benchmark.log"

OUTPUT_BASE="${OUTPUT_BASE:-$RESULTS_DIR}"

mkdir -p "$RESULTS_DIR"
mkdir -p "$INPUTS_DIR"
mkdir -p "$OUTPUT_BASE"

# ============================================================================
# USER-CONFIGURABLE PARAMETERS
# ============================================================================

CHROMS="${CHROMS:-chr22}"

SAMPLES_LIST="${SAMPLES_LIST:-10 50 100}"

CPU_THREADS_LIST="${CPU_THREADS_LIST:-1 2 4 8}"

# Pool of GPU device indices. The script sweeps cumulative prefixes of this
# pool. Examples:
#   GPU_DEVICES=0              → one GPU run: --gpu-devices 0
#   GPU_DEVICES=0,1            → --gpu-devices 0 ; --gpu-devices 0,1
#   GPU_DEVICES=0,1,2,3        → 0 ; 0,1 ; 0,1,2 ; 0,1,2,3
#   GPU_DEVICES=2,3            → 2 ; 2,3
GPU_DEVICES="${GPU_DEVICES:-0,1}"

GPU_THREADS_LIST="${GPU_THREADS_LIST:-1 2 4}"

REQUESTED_MODE="${MODE:-auto}"   # auto | cpu | gpu | both

BASE_VCF="${BASE_VCF:-${BENCH_DIR}/datasets/1000G/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz}"

CHUNK_FA="${CHUNK_FA:-${BENCH_DIR}/datasets/reference/GRCh38.primary_assembly.genome.fa}"

# ============================================================================
# CHECK REQUIRED FILES / PROGRAMS
# ============================================================================

echo "=========================================================================="
echo "🔧 Checking benchmark dependencies..."
echo "=========================================================================="

if [ ! -x "$RUST_BIN" ]; then
    echo "❌ Rust binary not found at:"
    echo "   $RUST_BIN"
    echo
    echo "Build with:"
    echo "   cargo build --release --features cuda"
    exit 1
fi

if [ ! -f "$BASE_VCF" ]; then
    echo "❌ 1000G VCF dataset not found at:"
    echo "   $BASE_VCF"
    exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "❌ bcftools is not installed."
    exit 1
fi

if ! command -v tabix >/dev/null 2>&1; then
    echo "❌ tabix is not installed."
    exit 1
fi

if ! command -v bc >/dev/null 2>&1; then
    echo "❌ bc is not installed."
    exit 1
fi

# ============================================================================
# TIME FORMATTER
# ============================================================================

format_readable_time() {
    local total_sec="$1"

    local hours
    local minutes
    local seconds
    local millis
    local remainder

    hours=$(echo "$total_sec / 3600" | bc)
    remainder=$(echo "$total_sec - ($hours * 3600)" | bc)
    minutes=$(echo "$remainder / 60" | bc)
    remainder=$(echo "$remainder - ($minutes * 60)" | bc)
    seconds=$(echo "$remainder / 1" | bc)
    millis=$(echo "($remainder - $seconds) * 1000" | bc | awk '{print int($1+0.5)}')

    local out=""

    if [ "$hours" -gt 0 ]; then
        out="${hours}h "
    fi

    if [ "$minutes" -gt 0 ] || [ "$hours" -gt 0 ]; then
        out="${out}${minutes}m "
    fi

    out="${out}${seconds}s ${millis}ms"

    echo "$out"
}

# ============================================================================
# GPU PREFIX SWEEP HELPERS
# ============================================================================
# Given a comma-separated pool "0,1,2,3", emit each cumulative prefix on its
# own line:
#
#   0
#   0,1
#   0,1,2
#   0,1,2,3
#
# The literal string is passed verbatim to the tool's --gpu-devices.
#
# A "safe" label is derived from each prefix by replacing commas with
# underscores, for use in filenames and column names:
#
#   0        → 0
#   0,1      → 0_1
#   0,1,2    → 0_1_2
#
# This keeps the tool invocation and the labels in lockstep: whatever you
# ask for in --gpu-devices is what shows up in the log column header.

gpu_device_prefixes() {
    local list="$1"
    local -a devs
    IFS=',' read -ra devs <<< "$list"

    local prefix=""
    local d
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

gpu_label_safe() {
    echo "$1" | tr ',' '_'
}

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
            EXEC_MODE="both"
            echo "ℹ️  MODE=auto and NVIDIA GPU detected → running both CPU and GPU sweeps."
        else
            EXEC_MODE="cpu"
            echo "ℹ️  MODE=auto and no GPU detected → running CPU-only sweep."
        fi
        ;;
    cpu)
        EXEC_MODE="cpu"
        echo "ℹ️  MODE=cpu → running CPU-only sweep."
        ;;
    gpu)
        if [ "$HAS_GPU" = true ]; then
            EXEC_MODE="gpu"
            echo "ℹ️  MODE=gpu and NVIDIA GPU detected → running GPU-only sweep."
        else
            EXEC_MODE="cpu"
            echo "⚠️  MODE=gpu requested but no NVIDIA GPU detected → falling back to CPU-only."
        fi
        ;;
    both)
        if [ "$HAS_GPU" = true ]; then
            EXEC_MODE="both"
            echo "ℹ️  MODE=both and NVIDIA GPU detected → running GPU and CPU sweeps."
        else
            EXEC_MODE="cpu"
            echo "⚠️  MODE=both requested but no NVIDIA GPU detected → running CPU-only."
        fi
        ;;
    *)
        echo "❌ Invalid MODE='$REQUESTED_MODE' (expected: auto, cpu, gpu, both)"
        exit 1
        ;;
esac

RUN_CPU=false
RUN_GPU=false
case "$EXEC_MODE" in
    cpu)  RUN_CPU=true ;;
    gpu)  RUN_GPU=true ;;
    both) RUN_CPU=true; RUN_GPU=true ;;
esac

# ============================================================================
# CONFIGURATION SUMMARY
# ============================================================================

echo
echo "=========================================================================="
echo "⚙️  Benchmark Configuration"
echo "=========================================================================="

echo "VCF:              $BASE_VCF"
echo "Reference:        $CHUNK_FA"
echo "Chromosomes:      $CHROMS"
echo "Sample tiers:     $SAMPLES_LIST"
echo "Execution mode:   $EXEC_MODE   (RUN_CPU=$RUN_CPU, RUN_GPU=$RUN_GPU)"
echo "Output base:      $OUTPUT_BASE"
echo "Log file:         $LOG_FILE"

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
echo

# ============================================================================
# CHECK ORIGINAL VCF
# ============================================================================

echo "🔍 Checking original VCF..."

BASE_VCF_CHECK_ERR="${RESULTS_DIR}/base_vcf_check.stderr"

if ! bcftools view -h "$BASE_VCF" >/dev/null 2>"$BASE_VCF_CHECK_ERR"; then
    echo "❌ Cannot read original VCF."

    if [ -s "$BASE_VCF_CHECK_ERR" ]; then
        cat "$BASE_VCF_CHECK_ERR"
    fi

    exit 1
fi

echo "✅ Original VCF can be read."

# ============================================================================
# REFRESH ORIGINAL VCF INDEX
# ============================================================================

echo "ℹ️  Refreshing original VCF index..."

if ! tabix -f -p vcf "$BASE_VCF"; then
    echo "❌ Failed to create/refresh original VCF index."
    exit 1
fi

echo "✅ Original VCF index refreshed."

# ============================================================================
# LOAD COMPLETE SAMPLE LIST
# ============================================================================

echo
echo "📊 Reading sample list from original VCF..."

ALL_SAMPLE_NAMES_FILE="${RESULTS_DIR}/all_samples.txt"

if ! bcftools query -l "$BASE_VCF" > "$ALL_SAMPLE_NAMES_FILE" 2>"${RESULTS_DIR}/sample_query.stderr"; then
    echo "❌ Failed to read sample names from VCF."

    if [ -s "${RESULTS_DIR}/sample_query.stderr" ]; then
        cat "${RESULTS_DIR}/sample_query.stderr"
    fi

    exit 1
fi

TOTAL_SAMPLES=$(wc -l < "$ALL_SAMPLE_NAMES_FILE")
TOTAL_SAMPLES=$(echo "$TOTAL_SAMPLES" | tr -d '[:space:]')

echo "📊 Total samples available in 1000G VCF: $TOTAL_SAMPLES"

# ============================================================================
# CREATE RESULT LOG HEADER
# ============================================================================

HEADER="Chromosome\tSamples\tMode\tExtraction_Time"

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
# MAIN BENCHMARK
# ============================================================================

for CHROM in $CHROMS; do

    for NUM_SAMPLES in $SAMPLES_LIST; do

        echo
        echo "=========================================================================="
        echo "📦 Target Chromosome: $CHROM | Samples: $NUM_SAMPLES (1000G dataset)"
        echo "=========================================================================="

        if ! [[ "$NUM_SAMPLES" =~ ^[0-9]+$ ]]; then
            echo "❌ Invalid sample count: $NUM_SAMPLES"
            exit 1
        fi

        if [ "$NUM_SAMPLES" -le 0 ]; then
            echo "❌ Sample count must be greater than zero."
            exit 1
        fi

        if [ "$NUM_SAMPLES" -gt "$TOTAL_SAMPLES" ]; then
            echo "❌ Requested $NUM_SAMPLES samples."
            echo "   Only $TOTAL_SAMPLES samples are available."
            exit 1
        fi

        SAMPLE_RES_DIR="${OUTPUT_BASE}/samples_${NUM_SAMPLES}/${CHROM}"

        mkdir -p \
            "${SAMPLE_RES_DIR}/errors" \
            "${SAMPLE_RES_DIR}/outputs"

        SAMPLE_VCF="${INPUTS_DIR}/${CHROM}_${NUM_SAMPLES}samples.vcf.gz"

        NEED_EXTRACTION=true

        if [ -f "$SAMPLE_VCF" ] && [ -f "${SAMPLE_VCF}.tbi" ]; then

            echo "ℹ️  Extracted VCF already exists:"
            echo "   $SAMPLE_VCF"

            EXISTING_SAMPLE_COUNT=$(
                bcftools query -l "$SAMPLE_VCF" 2>/dev/null | wc -l
            )
            EXISTING_SAMPLE_COUNT=$(echo "$EXISTING_SAMPLE_COUNT" | tr -d '[:space:]')

            if [ "$EXISTING_SAMPLE_COUNT" -eq "$NUM_SAMPLES" ]; then

                echo "✅ Existing VCF contains exactly $NUM_SAMPLES samples."
                NEED_EXTRACTION=false
                EXTRACTION_TIME="0"

            else

                echo "⚠️  Existing VCF contains $EXISTING_SAMPLE_COUNT samples."
                echo "   Expected $NUM_SAMPLES."
                echo "   Recreating sample subset."

                rm -f "$SAMPLE_VCF"
                rm -f "${SAMPLE_VCF}.tbi"

            fi

        fi

        if [ "$NEED_EXTRACTION" = true ]; then

            echo
            echo "└─ Extracting first $NUM_SAMPLES samples from 1000G VCF..."

            SAMPLE_NAMES_FILE="${RESULTS_DIR}/selected_${CHROM}_${NUM_SAMPLES}samples.txt"

            sed -n "1,${NUM_SAMPLES}p" \
                "$ALL_SAMPLE_NAMES_FILE" \
                > "$SAMPLE_NAMES_FILE"

            SELECTED_SAMPLE_COUNT=$(wc -l < "$SAMPLE_NAMES_FILE")
            SELECTED_SAMPLE_COUNT=$(echo "$SELECTED_SAMPLE_COUNT" | tr -d '[:space:]')

            if [ "$SELECTED_SAMPLE_COUNT" -ne "$NUM_SAMPLES" ]; then

                echo "❌ Failed to select requested samples."
                echo "   Requested: $NUM_SAMPLES"
                echo "   Selected:  $SELECTED_SAMPLE_COUNT"

                exit 1

            fi

            echo "   ✅ Selected $SELECTED_SAMPLE_COUNT sample names."

            ALL_SAMPLES=$(paste -sd, "$SAMPLE_NAMES_FILE")

            if [ -z "$ALL_SAMPLES" ]; then

                echo "❌ Generated sample list is empty."

                exit 1

            fi

            TEMP_VCF="${SAMPLE_VCF}.tmp"

            rm -f "$TEMP_VCF"
            rm -f "${TEMP_VCF}.tbi"

            BCFT_STDERR="${RESULTS_DIR}/errors/bcftools_extract.stderr"
            BCFT_STDOUT="${RESULTS_DIR}/errors/bcftools_extract.stdout"

            mkdir -p "${RESULTS_DIR}/errors"
            : > "$BCFT_STDERR"
            : > "$BCFT_STDOUT"

            echo "   └─ Running bcftools view..."

            START_EXTRACT=$(date +%s.%N)

            set +e

            bcftools view \
                -s "$ALL_SAMPLES" \
                -r "$CHROM" \
                "$BASE_VCF" \
                -Oz \
                -o "$TEMP_VCF" \
                >"$BCFT_STDOUT" \
                2>"$BCFT_STDERR"

            BCFT_EXIT=$?

            END_EXTRACT=$(date +%s.%N)

            set -e

            EXTRACTION_TIME=$(echo "$END_EXTRACT - $START_EXTRACT" | bc)

            if [ "$BCFT_EXIT" -ne 0 ]; then

                echo
                echo "❌ bcftools extraction FAILED."
                echo "   Exit code: $BCFT_EXIT"
                echo
                echo "   stderr:"
                cat "$BCFT_STDERR"

                rm -f "$TEMP_VCF"
                rm -f "${TEMP_VCF}.tbi"

                exit 1

            fi

            if [ ! -s "$TEMP_VCF" ]; then

                echo "❌ bcftools returned success but produced no VCF."

                echo
                echo "   stderr:"
                cat "$BCFT_STDERR"

                rm -f "$TEMP_VCF"
                rm -f "${TEMP_VCF}.tbi"

                exit 1

            fi

            echo "   ✅ bcftools extraction finished."
            echo "   Extraction time: $(format_readable_time "$EXTRACTION_TIME")"

            echo "   └─ Indexing extracted VCF..."

            if ! tabix -f -p vcf "$TEMP_VCF"; then

                echo "❌ Failed to index extracted VCF."

                rm -f "$TEMP_VCF"
                rm -f "${TEMP_VCF}.tbi"

                exit 1

            fi

            mv "$TEMP_VCF" "$SAMPLE_VCF"
            mv "${TEMP_VCF}.tbi" "${SAMPLE_VCF}.tbi"

            EXTRACTED_SAMPLE_COUNT=$(
                bcftools query -l "$SAMPLE_VCF" 2>/dev/null | wc -l
            )
            EXTRACTED_SAMPLE_COUNT=$(echo "$EXTRACTED_SAMPLE_COUNT" | tr -d '[:space:]')

            if [ "$EXTRACTED_SAMPLE_COUNT" -ne "$NUM_SAMPLES" ]; then

                echo "❌ Extracted VCF has incorrect sample count."
                echo "   Expected: $NUM_SAMPLES"
                echo "   Found:    $EXTRACTED_SAMPLE_COUNT"

                exit 1

            fi

            echo "   ✅ Extracted VCF verified: $EXTRACTED_SAMPLE_COUNT samples."

        fi

        if [ ! -s "$SAMPLE_VCF" ]; then

            echo "❌ Sample VCF does not exist or is empty:"
            echo "   $SAMPLE_VCF"

            exit 1

        fi

        FINAL_SAMPLE_COUNT=$(
            bcftools query -l "$SAMPLE_VCF" 2>/dev/null | wc -l
        )
        FINAL_SAMPLE_COUNT=$(echo "$FINAL_SAMPLE_COUNT" | tr -d '[:space:]')

        if [ "$FINAL_SAMPLE_COUNT" -ne "$NUM_SAMPLES" ]; then

            echo "❌ Final VCF sample count mismatch."
            echo "   Expected: $NUM_SAMPLES"
            echo "   Found:    $FINAL_SAMPLE_COUNT"

            exit 1

        fi

        echo
        echo "📄 Benchmark input:"
        echo "   $SAMPLE_VCF"
        echo "   Samples: $FINAL_SAMPLE_COUNT"

        declare -A RUN_TIMES

        NOTE_TEXT="Success"

        # ====================================================================
        # GPU BENCHMARK (prefix sweep × thread sweep)
        # ====================================================================

        if [ "$RUN_GPU" = true ]; then

            while IFS= read -r DEV_LIST; do

                SAFE=$(gpu_label_safe "$DEV_LIST")
                NDEVS=$(awk -F, '{print NF}' <<< "$DEV_LIST")

                for GT in $GPU_THREADS_LIST; do

                    CONFIG_DIR="${SAMPLE_RES_DIR}/gpu_dev${SAFE}_t${GT}"
                    mkdir -p "$CONFIG_DIR"

                    GPU_STDOUT="${SAMPLE_RES_DIR}/errors/gpu_dev${SAFE}_t${GT}.stdout"
                    GPU_STDERR="${SAMPLE_RES_DIR}/errors/gpu_dev${SAFE}_t${GT}.stderr"

                    echo
                    echo "  └─ Running GPU Engine [${NDEVS} GPU(s): --gpu-devices ${DEV_LIST}, Host threads: ${GT}]..."

                    START=$(date +%s.%N)

                    set +e

                    "$RUST_BIN" \
                        --device gpu \
                        --gpu-devices "$DEV_LIST" \
                        --threads "$GT" \
                        --reference "$CHUNK_FA" \
                        "$SAMPLE_VCF" \
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

                        echo "     ✅ GPU dev${SAFE} / t${GT} finished successfully in $READABLE_T"

                    else

                        ERR_MSG=$(head -n 1 "$GPU_STDERR" 2>/dev/null || true)
                        echo "     ❌ GPU dev${SAFE} / t${GT} failed (Exit $EXIT): ${ERR_MSG:-unknown}"

                        if [ "$NOTE_TEXT" = "Success" ]; then
                            NOTE_TEXT="GPU dev${SAFE}_t${GT} Error"
                        else
                            NOTE_TEXT="${NOTE_TEXT}; GPU dev${SAFE}_t${GT} Error"
                        fi

                    fi

                done

            done < <(gpu_device_prefixes "$GPU_DEVICES")

        fi

        # ====================================================================
        # CPU BENCHMARK
        # ====================================================================

        if [ "$RUN_CPU" = true ]; then

            for T in $CPU_THREADS_LIST; do

                CONFIG_DIR="${SAMPLE_RES_DIR}/cpu_t${T}"
                mkdir -p "$CONFIG_DIR"

                CPU_STDOUT="${SAMPLE_RES_DIR}/errors/cpu_t${T}.stdout"
                CPU_STDERR="${SAMPLE_RES_DIR}/errors/cpu_t${T}.stderr"

                echo
                echo "  └─ Running CPU Engine [Threads: $T]..."

                START=$(date +%s.%N)

                set +e

                "$RUST_BIN" \
                    --device cpu \
                    --threads "$T" \
                    --reference "$CHUNK_FA" \
                    "$SAMPLE_VCF" \
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

                    echo "     ✅ CPU t${T} finished successfully in $READABLE_T"

                else

                    ERR_MSG=$(head -n 1 "$CPU_STDERR" 2>/dev/null || true)
                    echo "     ❌ CPU t${T} failed (Exit $EXIT): ${ERR_MSG:-unknown}"

                    if [ "$NOTE_TEXT" = "Success" ]; then
                        NOTE_TEXT="CPU t${T} Error"
                    else
                        NOTE_TEXT="${NOTE_TEXT}; CPU t${T} Error"
                    fi

                fi

            done

        fi

        # ====================================================================
        # WRITE RESULTS
        # ====================================================================

        ROW_DATA="${CHROM}\t${NUM_SAMPLES}\t${EXEC_MODE}\t$(format_readable_time "$EXTRACTION_TIME")"

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

        # ====================================================================
        # SAMPLE-TIER SUMMARY
        # ====================================================================

        echo
        echo "──────────────────────────────────────────────────────────────────────────"
        echo "📊 Completed sample tier: $NUM_SAMPLES"
        echo "   Extraction: $(format_readable_time "$EXTRACTION_TIME")"

        if [ "$RUN_GPU" = true ]; then
            while IFS= read -r DEV_LIST; do
                SAFE=$(gpu_label_safe "$DEV_LIST")
                NDEVS=$(awk -F, '{print NF}' <<< "$DEV_LIST")
                for GT in $GPU_THREADS_LIST; do
                    G_TIME="${RUN_TIMES["gpu_dev${SAFE}_t${GT}"]:-0}"
                    echo "   GPU ${NDEVS}g (dev=${DEV_LIST}) / t${GT}:  $(format_readable_time "$G_TIME")"
                done
            done < <(gpu_device_prefixes "$GPU_DEVICES")
        fi

        if [ "$RUN_CPU" = true ]; then
            for T in $CPU_THREADS_LIST; do
                T_TIME="${RUN_TIMES["cpu_t${T}"]:-0}"
                echo "   CPU t${T}:   $(format_readable_time "$T_TIME")"
            done
        fi

        echo "──────────────────────────────────────────────────────────────────────────"

    done

done

# ============================================================================
# FINAL
# ============================================================================

echo
echo "=========================================================================="
echo "✅ Hierarchical combinatorial scaling benchmark complete."
echo "=========================================================================="
echo
echo "Results saved to:"
echo "   $LOG_FILE"
echo
echo "Input VCFs saved to:"
echo "   $INPUTS_DIR"
echo
echo "Output FASTAs written under:"
echo "   $OUTPUT_BASE"
echo
echo "Extraction time is recorded separately and is NOT included"
echo "in CPU/GPU benchmark timings."
echo "=========================================================================="