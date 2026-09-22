#!/usr/bin/env bash
# ============================================================================
# vcf2fasta-rust benchmark — per-sample, per-chromosome, CPU/GPU/both, optional
# vcflib cross-check. No merging: every sample is processed in isolation and
# writes to its own sample-named subdirectory.
#
# Ploidy is auto-detected from the FASTA files the tool actually produces:
#   - two files (":0.fa" and ":1.fa") → diploid
#   - one file  (":0.fa" only)        → haploid (e.g. male chrX/chrY)
#   - zero files                       → error
#
# Environment variables (all optional):
#   MODE              cpu | gpu | both         (default: cpu)
#   USE_VCFLIB        true | false             (default: false)
#   DATASET           platinum | hybrid        (default: platinum)
#   SAMPLES           space-separated names    (default: "NA12877 NA12878")
#   CHROMS            space-separated contigs  (default: chr22)
#   CPU_THREADS_LIST  CPU thread sweep         (default: "1 2 4 8")
#   GPU_THREADS_LIST  GPU host-thread sweep    (default: "1 2")
#   THREADS_LIST      legacy fallback for both (default: unset)
#   GPU_DEVICES       comma-separated GPU ids  (default: 0)
#   CHUNK_FILES       true | false             (default: true)
#   LIST_SKIPPED      max skipped variants to print (default: 20)
#
# Length-diff explainer (used when USE_VCFLIB=true):
#   COMPARE_HAPS      haplotypes to compare               (default: "0 1")
#   COMPARE_MAX_DIFFS max skipped variants to list        (default: 20)
#   COMPARE_TSV       aggregation TSV path                (default: results/compare_summary.tsv)
#   COMPARE_SCRIPT    path to explain_length_diff.py      (auto)
# ============================================================================

set -uo pipefail

if [ -d /usr/local/cuda-12.8/targets/x86_64-linux/lib ]; then
    export LD_LIBRARY_PATH="/usr/local/cuda-12.8/targets/x86_64-linux/lib:/usr/local/cuda-12.8/lib64:/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(dirname "$SCRIPT_DIR")"
RUST_BIN="${BENCH_DIR}/../../target/release/vcf2fasta"

RESULTS_DIR="${BENCH_DIR}/benchmark_results"
LOG_FILE="${RESULTS_DIR}/benchmark_results.log"

# --- Config ------------------------------------------------------------------
MODE="${MODE:-cpu}"
USE_VCFLIB="${USE_VCFLIB:-false}"
DATASET="${DATASET:-platinum}"
SAMPLES="${SAMPLES:-NA12877 NA12878}"
CHROMS="${CHROMS:-chr22}"
CHUNK_FILES="${CHUNK_FILES:-true}"
GPU_DEVICES="${GPU_DEVICES:-0}"
LIST_SKIPPED="${LIST_SKIPPED:-20}"

CPU_THREADS_LIST="${CPU_THREADS_LIST:-${THREADS_LIST:-1 2 4 8}}"
GPU_THREADS_LIST="${GPU_THREADS_LIST:-${THREADS_LIST:-1 2}}"

# --- Length-diff explainer config -------------------------------------------
COMPARE_HAPS="${COMPARE_HAPS:-0 1}"
COMPARE_MAX_DIFFS="${COMPARE_MAX_DIFFS:-20}"
COMPARE_TSV="${COMPARE_TSV:-${RESULTS_DIR}/compare_summary.tsv}"
COMPARE_SCRIPT="${COMPARE_SCRIPT:-${SCRIPT_DIR}/../validation/explain_length_diff.py}"

case "$MODE" in
    cpu|gpu|both) ;;
    *) echo "❌ Invalid MODE='$MODE' (expected: cpu, gpu, or both)"; exit 1 ;;
esac

# --- Helpers -----------------------------------------------------------------
threads_for_engine() {
    case "$1" in
        cpu) echo "$CPU_THREADS_LIST" ;;
        gpu) echo "$GPU_THREADS_LIST" ;;
    esac
}
first_thread_of() { echo "$1" | awk '{print $1}'; }

format_num() {
    awk -v n="$1" 'BEGIN{
        s = sprintf("%d", n)
        out = ""
        while (length(s) > 3) {
            out = "," substr(s, length(s)-2) out
            s = substr(s, 1, length(s)-3)
        }
        print s out
    }'
}

declare -A ENGINE_FIRST_T
ENGINE_FIRST_T[cpu]=$(first_thread_of "$CPU_THREADS_LIST")
ENGINE_FIRST_T[gpu]=$(first_thread_of "$GPU_THREADS_LIST")

# --- Directories -------------------------------------------------------------
mkdir -p "${RESULTS_DIR}/errors" \
         "${RESULTS_DIR}/chunked_inputs"

printf 'Row_type\tHaplotype_id\tRust_bases\tvcflib_bases\tDelta_bases\tSkipped_variant_count\tFirst_skipped_vcf_pos\tRust_display\tvcflib_display\n' \
    > "$COMPARE_TSV"

case "$MODE" in
    cpu|both)
        for T in $CPU_THREADS_LIST; do
            for S in $SAMPLES; do
                mkdir -p "${RESULTS_DIR}/new_cpu_t${T}/${S}"
            done
        done
        ;;
esac
case "$MODE" in
    gpu|both)
        for T in $GPU_THREADS_LIST; do
            for S in $SAMPLES; do
                mkdir -p "${RESULTS_DIR}/new_gpu_t${T}/${S}"
            done
        done
        ;;
esac

if [ "$USE_VCFLIB" = "true" ]; then
    for S in $SAMPLES; do
        mkdir -p "${RESULTS_DIR}/vcflib_out/${S}"
    done
fi

if [ ! -x "$RUST_BIN" ]; then
    echo "❌ Rust binary not found at $RUST_BIN. Build with: cargo build --release --features cuda"
    exit 1
fi

# --- vcflib lookup -----------------------------------------------------------
VCFLIB_BIN=""
if [ "$USE_VCFLIB" = "true" ]; then
    OFFICIAL=$(command -v vcf2fasta 2>/dev/null || true)
    if [ -z "$OFFICIAL" ]; then
        echo "⚠️  USE_VCFLIB=true but no 'vcf2fasta' on PATH. Disabling vcflib."
        USE_VCFLIB=false
    else
        OUR_REAL=$(readlink -f "$RUST_BIN" 2>/dev/null || echo "")
        OFFICIAL_REAL=$(readlink -f "$OFFICIAL" 2>/dev/null || echo "")
        if [ -n "$OUR_REAL" ] && [ "$OFFICIAL_REAL" = "$OUR_REAL" ]; then
            echo "⚠️  'vcf2fasta' on PATH resolves to our binary. Disabling vcflib."
            USE_VCFLIB=false
        else
            VCFLIB_BIN="$OFFICIAL"
            echo "ℹ️  vcflib binary: $VCFLIB_BIN"
        fi
    fi
fi

if [ "$USE_VCFLIB" = "true" ] && [ ! -f "$COMPARE_SCRIPT" ]; then
    echo "⚠️  explain_length_diff.py not found at $COMPARE_SCRIPT"
    echo "    Direct Rust-vs-vcflib comparison will be skipped."
fi

format_readable_time() {
    local total_sec=$1
    if [ "$(echo "$total_sec == 0" | bc)" -eq 1 ]; then
        echo "0s 0ms"; return
    fi
    local hours=$(echo "$total_sec / 3600" | bc)
    local remainder=$(echo "$total_sec % 3600" | bc)
    local minutes=$(echo "$remainder / 60" | bc)
    remainder=$(echo "$remainder % 60" | bc)
    local seconds=$(echo "$remainder / 1" | bc)
    local millis=$(echo "($remainder - $seconds) * 1000" | bc | awk '{print int($1+0.5)}')
    local out=""
    [ "$hours" -gt 0 ] && out="${hours}h "
    { [ "$minutes" -gt 0 ] || [ "$hours" -gt 0 ]; } && out="${out}${minutes}m "
    out="${out}${seconds}s ${millis}ms"
    echo "$out"
}

# --- Engine list -------------------------------------------------------------
ENGINES=()
case "$MODE" in
    cpu)  ENGINES=("cpu") ;;
    gpu)  ENGINES=("gpu") ;;
    both) ENGINES=("cpu" "gpu") ;;
esac

# --- Header for the main log --------------------------------------------------
HEADER="Chromosome\tSample\tPloidy"
[ "$USE_VCFLIB" = "true" ] && HEADER="${HEADER}\tvcflib"

for E in "${ENGINES[@]}"; do
    for T in $(threads_for_engine "$E"); do
        HEADER="${HEADER}\t${E}_t${T}_time"
    done
done

for E in "${ENGINES[@]}"; do
    FIRST=$(first_thread_of "$(threads_for_engine "$E")")
    for T in $(threads_for_engine "$E"); do
        [ "$T" != "$FIRST" ] && HEADER="${HEADER}\tspeedup_${E}_t${T}_vs_t${FIRST}"
    done
done

if [ "$USE_VCFLIB" = "true" ]; then
    for E in "${ENGINES[@]}"; do
        for T in $(threads_for_engine "$E"); do
            HEADER="${HEADER}\tspeedup_${E}_t${T}_vs_vcflib"
        done
    done
fi

HEADER="${HEADER}\tNotes"
echo -e "$HEADER" > "$LOG_FILE"

# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------
prepare_reference() {
    local chrom="$1"
    if [ "$CHUNK_FILES" = "true" ]; then
        local chunk_fa="${RESULTS_DIR}/chunked_inputs/${chrom}_ref.fa"
        if [ ! -f "$chunk_fa" ]; then
            samtools faidx "${BENCH_DIR}/datasets/reference/GRCh38.primary_assembly.genome.fa" "$chrom" > "$chunk_fa"
            samtools faidx "$chunk_fa"
        fi
        echo "$chunk_fa"
    else
        echo "${BENCH_DIR}/datasets/reference/GRCh38.primary_assembly.genome.fa"
    fi
}

prepare_sample_vcf() {
    local sample="$1" chrom="$2" src

    if [ "$DATASET" = "platinum" ]; then
        src="${BENCH_DIR}/datasets/platinum/${sample}/${sample}.vcf.gz"
    else
        src="${BENCH_DIR}/datasets/hybrid/hg38.hybrid.vcf.gz"
    fi

    if [ ! -f "$src" ]; then
        echo "❌ Missing source VCF: $src" >&2
        return 1
    fi

    local out="${RESULTS_DIR}/chunked_inputs/${sample}_${chrom}.vcf.gz"
    if [ ! -f "$out" ]; then
        if [ "$DATASET" = "platinum" ]; then
            bcftools view -r "$chrom" "$src" -O z -o "$out"
        else
            bcftools view -r "$chrom" -s "$sample" "$src" -O z -o "$out"
        fi
        bcftools index -t -f "$out"
    fi
    echo "$out"
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
for CHROM in $CHROMS; do
    echo "=========================================================================="
    echo "📦 Target Chromosome: $CHROM"

    CHUNK_FA=$(prepare_reference "$CHROM")

    for S in $SAMPLES; do
        echo "  ── Sample: $S"
        if ! SAMPLE_VCF=$(prepare_sample_vcf "$S" "$CHROM"); then
            echo "     ❌ Failed to prepare VCF for sample $S / $CHROM — skipping"
            continue
        fi
        echo "     ℹ️  VCF: $SAMPLE_VCF"

        declare -A TIMES
        declare -A EXITS

        # ---- Tool runs ----------------------------------------------------
        for E in "${ENGINES[@]}"; do
            for T in $(threads_for_engine "$E"); do
                KEY="${E}_t${T}"
                EXTRA=""
                [ "$E" = "gpu" ] && EXTRA=" [GPUs: $GPU_DEVICES]"
                echo "     └─ Running ${E^^} [Host threads: $T]${EXTRA}..."

                OUT_DIR="${RESULTS_DIR}/new_${E}_t${T}/${S}"
                mkdir -p "$OUT_DIR"
                PREFIX="${OUT_DIR}/loai_${CHROM}"

                START=$(date +%s.%N)
                set +e
                if [ "$E" = "gpu" ]; then
                    "$RUST_BIN" \
                        --device gpu --gpu-devices "$GPU_DEVICES" --threads "$T" \
                        --reference "$CHUNK_FA" "$SAMPLE_VCF" --no-validate-ref \
                        --prefix "$PREFIX" \
                        > "${RESULTS_DIR}/errors/${CHROM}_${S}_${E}_t${T}.stdout" \
                        2> "${RESULTS_DIR}/errors/${CHROM}_${S}_${E}_t${T}.stderr"
                else
                    "$RUST_BIN" \
                        --device cpu --threads "$T" \
                        --reference "$CHUNK_FA" "$SAMPLE_VCF" --no-validate-ref \
                        --prefix "$PREFIX" \
                        > "${RESULTS_DIR}/errors/${CHROM}_${S}_${E}_t${T}.stdout" \
                        2> "${RESULTS_DIR}/errors/${CHROM}_${S}_${E}_t${T}.stderr"
                fi
                EXIT=$?
                set -e
                ELAPSED=$(echo "$(date +%s.%N) - $START" | bc)
                TIMES["$KEY"]="$ELAPSED"
                EXITS["$KEY"]="$EXIT"

                if [ $EXIT -eq 0 ]; then
                    echo "        ✅ ${E^^} t${T} finished in $(format_readable_time "$ELAPSED")"
                else
                    ERR=$(head -n 1 "${RESULTS_DIR}/errors/${CHROM}_${S}_${E}_t${T}.stderr" 2>/dev/null || echo "unknown")
                    echo "        ❌ ${E^^} t${T} failed (exit $EXIT): ${ERR}"
                fi
            done
        done

        # ---- vcflib (per sample) -----------------------------------------
        VCFLIB_TIME="0.00"; VCFLIB_EXISTS=false
        VCFLIB_FILES_OK=false
        VCF_PREFIX=""
        if [ "$USE_VCFLIB" = "true" ]; then
            echo "     └─ Running Legacy vcflib..."
            VCF_DIR="${RESULTS_DIR}/vcflib_out/${S}"
            mkdir -p "$VCF_DIR"
            VCF_PREFIX="${VCF_DIR}/vcflib_${CHROM}"

            START=$(date +%s.%N)
            set +e
            "$VCFLIB_BIN" \
                --reference "$CHUNK_FA" \
                --prefix "$VCF_PREFIX" \
                "$SAMPLE_VCF" \
                > "${RESULTS_DIR}/errors/${CHROM}_${S}_vcflib.stdout" \
                2> "${RESULTS_DIR}/errors/${CHROM}_${S}_vcflib.stderr"
            VCF_EXIT=$?
            set -e
            VCFLIB_TIME=$(echo "$(date +%s.%N) - $START" | bc)
            READABLE_VCF_T=$(format_readable_time "$VCFLIB_TIME")
            if [ $VCF_EXIT -eq 0 ]; then
                VCFLIB_EXISTS=true
                echo "        ✅ vcflib finished in $READABLE_VCF_T"
            else
                ERR=$(head -n 1 "${RESULTS_DIR}/errors/${CHROM}_${S}_vcflib.stderr" 2>/dev/null || echo "unknown")
                echo "        ❌ vcflib failed (exit $VCF_EXIT) in $READABLE_VCF_T: ${ERR}"
            fi
        fi

        # ---- Ploidy detection --------------------------------------------
        PRIMARY_ENGINE="${ENGINES[0]}"
        PRIMARY_FIRST_T="${ENGINE_FIRST_T[$PRIMARY_ENGINE]}"
        PRIMARY_PREFIX="${RESULTS_DIR}/new_${PRIMARY_ENGINE}_t${PRIMARY_FIRST_T}/${S}/loai_${CHROM}"

        OUT_0="${PRIMARY_PREFIX}${S}_${CHROM}:0.fa"
        OUT_1="${PRIMARY_PREFIX}${S}_${CHROM}:1.fa"

        HAS_H0=false; HAS_H1=false
        [ -f "$OUT_0" ] && HAS_H0=true
        [ -f "$OUT_1" ] && HAS_H1=true

        if [ "$HAS_H0" = "true" ] && [ "$HAS_H1" = "true" ]; then
            PLOIDY=2
        elif [ "$HAS_H0" = "true" ]; then
            PLOIDY=1
        elif [ "$HAS_H1" = "true" ]; then
            PLOIDY=1
        else
            PLOIDY=0
        fi

        echo "     └─ Validating sample $S (against ${PRIMARY_ENGINE^^} t${PRIMARY_FIRST_T} output)..."
        echo "        ℹ️  Ploidy detected: $PLOIDY (hap0=$HAS_H0, hap1=$HAS_H1)"

        NOTE_TEXT=""
        RES=""

        if [ "$PLOIDY" = "0" ]; then
            echo "        ❌ Tool produced no FASTA for $S / $CHROM"
            echo "           expected: $OUT_0"
            echo "           expected: $OUT_1"
            NOTE_TEXT="Missing tool FASTA output"
        else
            VALIDATE_ARGS=("$SAMPLE_VCF" "$OUT_0")
            if [ "$PLOIDY" = "2" ]; then
                VALIDATE_ARGS+=("$OUT_1")
            fi
            VALIDATE_ARGS+=(--reference "$CHUNK_FA"
                            --no-validate-ref
                            --sample "$S"
                            --list-skipped "$LIST_SKIPPED")

            if [ "$VCFLIB_EXISTS" = "true" ]; then
                VCF_OUT_0="${VCF_PREFIX}${S}_${CHROM}:0.fa"
                VCF_OUT_1="${VCF_PREFIX}${S}_${CHROM}:1.fa"
                if [ -f "$VCF_OUT_0" ]; then
                    if [ "$PLOIDY" = "2" ] && [ -f "$VCF_OUT_1" ]; then
                        VALIDATE_ARGS+=(--vcflib-hap0 "$VCF_OUT_0" --vcflib-hap1 "$VCF_OUT_1")
                        VCFLIB_FILES_OK=true
                    elif [ "$PLOIDY" = "1" ]; then
                        VALIDATE_ARGS+=(--vcflib-hap0 "$VCF_OUT_0")
                        VCFLIB_FILES_OK=true
                    else
                        echo "        ⚠️  vcflib FASTA incomplete for $S — running tool-only check"
                    fi
                else
                    echo "        ⚠️  vcflib FASTA not found for $S — running tool-only check"
                fi
            fi

            RES=$(python3 "${SCRIPT_DIR}/../validation/validate_vcf_to_fasta.py" "${VALIDATE_ARGS[@]}" 2>&1 || true)

            NUM_VARS=$(grep -iE "Variants:"  <<< "$RES" | head -n1 | grep -oE '[0-9][0-9,]*'        | head -n1 || true)
            REF_B=$(   grep -iE "Reference:" <<< "$RES" | head -n1 | grep -oE '[0-9][0-9,]* bases' | head -n1 | awk '{print $1}' || true)
            H0_B=$(    grep -iE "Tool hap0:" <<< "$RES" | head -n1 | grep -oE '[0-9][0-9,]* bases' | head -n1 | awk '{print $1}' || true)
            H1_B=$(    grep -iE "Tool hap1:" <<< "$RES" | head -n1 | grep -oE '[0-9][0-9,]* bases' | head -n1 | awk '{print $1}' || true)
            POLICY=$(  grep -iE "Policy:"    <<< "$RES" | head -n1 | sed 's/^.*Policy:[[:space:]]*//' || true)
            [ -z "$NUM_VARS" ] && NUM_VARS="N/A"
            [ -z "$REF_B"   ] && REF_B="N/A"
            [ -z "$H0_B"    ] && H0_B="N/A"
            [ -z "$H1_B"    ] && H1_B="-"
            [ -z "$POLICY"  ] && POLICY="N/A"

            if [ "$PLOIDY" = "2" ]; then
                echo "        ℹ️  Variants: ${NUM_VARS}   Ref: ${REF_B} b   hap0: ${H0_B} b   hap1: ${H1_B} b"
            else
                echo "        ℹ️  Variants: ${NUM_VARS}   Ref: ${REF_B} b   hap0: ${H0_B} b   (haploid)"
            fi
            echo "        ℹ️  Policy:   ${POLICY}"

            SKIPPED_LINES=$(grep -E "^[[:space:]]*SKIPPED " <<< "$RES" || true)
            if [ -n "$SKIPPED_LINES" ]; then
                echo "        ── skipped variants ──"
                sed 's/^/        /' <<< "$SKIPPED_LINES"
            fi

            if echo "$RES" | grep -q "TOOL   vs expected: ✅ PASS"; then
                TOOL_VERDICT=pass
            elif echo "$RES" | grep -q "TOOL   vs expected: ❌ FAIL"; then
                TOOL_VERDICT=fail
            else
                TOOL_VERDICT=error
            fi

            case "$TOOL_VERDICT" in
                pass)
                    echo "        ✅ vcf2fasta-rust Validation: PASSED"
                    NOTE_TEXT="vcf2fasta-rust Validation Passed"
                    ;;
                fail)
                    echo "        ❌ vcf2fasta-rust Validation: FAILED"
                    echo "        ── validator output (first 20 lines) ──"
                    sed 's/^/           /' <<< "$RES" | head -n 20
                    NOTE_TEXT="vcf2fasta-rust Validation Failed"
                    ;;
                error)
                    echo "        ⚠️  vcf2fasta-rust Validation: VALIDATOR DID NOT RUN"
                    echo "        ── validator output (first 20 lines) ──"
                    sed 's/^/           /' <<< "$RES" | head -n 20
                    NOTE_TEXT="Validator error (no verdict)"
                    ;;
            esac

            if [ "$VCFLIB_FILES_OK" = true ]; then
                if echo "$RES" | grep -q "VCFLIB vs expected: ✅ PASS"; then
                    echo "        ✅ vcflib Validation: PASSED"
                    NOTE_TEXT="${NOTE_TEXT}; vcflib Validation Passed"
                elif echo "$RES" | grep -q "VCFLIB vs expected: ❌ FAIL"; then
                    echo "        ❌ vcflib Validation: FAILED"
                    NOTE_TEXT="${NOTE_TEXT}; vcflib Validation Failed"
                else
                    echo "        ⚠️  vcflib Validation: no verdict in validator output"
                    NOTE_TEXT="${NOTE_TEXT}; vcflib validator produced no verdict"
                fi
            fi
        fi

        # ---- Length-diff explainer (per haplotype) ----------------------
        if [ "$USE_VCFLIB" = "true" ] \
           && [ "$VCFLIB_EXISTS" = "true" ] \
           && [ -n "$VCF_PREFIX" ] \
           && [ "$PLOIDY" != "0" ] \
           && [ -f "$COMPARE_SCRIPT" ]; then
            for H in $COMPARE_HAPS; do
                RUST_H="${PRIMARY_PREFIX}${S}_${CHROM}:${H}.fa"
                VCF_H="${VCF_PREFIX}${S}_${CHROM}:${H}.fa"

                [ -f "$RUST_H" ] || continue

                if [ ! -f "$VCF_H" ]; then
                    echo "     └─ Length-diff H${H}: skipped (vcflib hap${H} FASTA missing)"
                    continue
                fi

                echo "     └─ Length-diff explainer H${H}: ${S}_${CHROM}_H${H}"
                CMP_OUT=$(python3 "$COMPARE_SCRIPT" \
                    "$SAMPLE_VCF" "$RUST_H" "$VCF_H" \
                    --reference "$CHUNK_FA" \
                    --sample "$S" --haplotype "$H" \
                    --chromosome "$CHROM" \
                    --label-a Rust --label-b vcflib \
                    --max-diffs "$COMPARE_MAX_DIFFS" \
                    --tsv 2>&1 || true)

                grep -v $'^COMPARE\t' <<< "$CMP_OUT" | sed 's/^/        /'
                TSV_LINE=$(grep $'^COMPARE\t' <<< "$CMP_OUT" | head -n1 || true)
                if [ -n "$TSV_LINE" ]; then
                    echo -e "$TSV_LINE" >> "$COMPARE_TSV"
                fi
            done
        fi

        # ---- Log row ------------------------------------------------------
        VCFLIB_STR="N/A"
        [ "$VCFLIB_EXISTS" = true ] && VCFLIB_STR=$(format_readable_time "$VCFLIB_TIME")

        ROW="${CHROM}\t${S}\t${PLOIDY}"
        [ "$USE_VCFLIB" = "true" ] && ROW="${ROW}\t${VCFLIB_STR}"

        for E in "${ENGINES[@]}"; do
            for T in $(threads_for_engine "$E"); do
                ROW="${ROW}\t$(format_readable_time "${TIMES["${E}_t${T}"]:-0}")"
            done
        done

        for E in "${ENGINES[@]}"; do
            FIRST="${ENGINE_FIRST_T[$E]}"
            BASE_TIME="${TIMES["${E}_t${FIRST}"]:-0}"
            for T in $(threads_for_engine "$E"); do
                [ "$T" = "$FIRST" ] && continue
                T_TIME="${TIMES["${E}_t${T}"]:-0}"
                if [ "$(echo "$T_TIME > 0 && $BASE_TIME > 0" | bc)" -eq 1 ]; then
                    ROW="${ROW}\t$(echo "scale=2; $BASE_TIME / $T_TIME" | bc)x"
                else
                    ROW="${ROW}\tN/A"
                fi
            done
        done

        if [ "$USE_VCFLIB" = "true" ]; then
            for E in "${ENGINES[@]}"; do
                for T in $(threads_for_engine "$E"); do
                    T_TIME="${TIMES["${E}_t${T}"]:-0}"
                    if [ "$VCFLIB_EXISTS" = true ] && [ "$(echo "$VCFLIB_TIME > 0 && $T_TIME > 0" | bc)" -eq 1 ]; then
                        ROW="${ROW}\t$(echo "scale=2; $VCFLIB_TIME / $T_TIME" | bc)x"
                    else
                        ROW="${ROW}\tN/A"
                    fi
                done
            done
        fi

        ROW="${ROW}\t${NOTE_TEXT:-Success}"
        echo -e "$ROW" >> "$LOG_FILE"

        unset TIMES EXITS
    done
done

# ---------------------------------------------------------------------------
# Summary: per-haplotype length-diff table + divergence detail
# ---------------------------------------------------------------------------
if [ -s "$COMPARE_TSV" ]; then
    echo ""
    echo "=========================================================================="
    echo "📊 Per-haplotype Rust vs vcflib DIVERGENCE SUMMARY"
    echo "=========================================================================="
    printf "%-30s %15s %15s %10s %13s %17s\n" \
        "Haplotype" "Rust (bases)" "vcflib (bases)" "Δ bases" "Divergences" "First divergence"
    printf "%-30s %15s %15s %10s %13s %17s\n" \
        "---------" "-----------" "-------------" "-------" "-----------" "----------------"

    while IFS=$'\t' read -r tag hap la lb delta nd first disp_a disp_b; do
        [ "$tag" = "COMPARE" ] || continue
        la_f=$(format_num "$la")
        lb_f=$(format_num "$lb")
        if [ "$first" = "-" ]; then
            first_disp="—"
        else
            first_disp=$(format_num "$first")
        fi
        printf "%-30s %15s %15s %10s %13s %17s\n" \
            "$hap" "$la_f" "$lb_f" "$delta" "$nd" "$first_disp"
    done < "$COMPARE_TSV"
    echo "=========================================================================="

    # ---- Per-haplotype divergence detail --------------------------------
    echo ""
    echo "── Per-haplotype divergence detail ──"
    any_detail=0
    while IFS=$'\t' read -r tag hap la lb delta nd first disp_a disp_b; do
        [ "$tag" = "COMPARE" ] || continue
        if [ -z "$disp_a" ] || [ "$disp_a" = "-" ]; then
            continue
        fi
        any_detail=1
        echo ""
        if [ "$first" = "-" ]; then
            echo "  ${hap}  (Δ = ${delta} bases; ${nd})"
        else
            echo "  ${hap}  (Δ = ${delta} bases; ${nd} divergence(s); first at ${first})"
        fi
        echo "    ${disp_a}"
        echo "    ${disp_b}"

        # Caret: point at the SECOND '[' in disp_a (the divergent bracket
        # in the dual-bracket form).
        #
        # NOTE: `sub` is a reserved awk function — do not use it as a
        # variable name. `rest` is safe.
        caret=$(awk -v s="$disp_a" 'BEGIN{
            i = index(s, "[");
            if (i == 0) { print -1; exit }
            rest = substr(s, i + 1);
            j = index(rest, "[");
            if (j == 0) { print i - 1 } else { print i + j - 1 }
        }')
        if [ "$caret" -ge 0 ]; then
            printf "    %*s^\n" "$caret" ""
        fi
    done < "$COMPARE_TSV"
    if [ "$any_detail" -eq 0 ]; then
        echo "  (no divergent content — all haplotypes identical)"
    fi
    echo "=========================================================================="

    echo ""
    echo "ℹ️  Raw rows: $COMPARE_TSV"
    echo "ℹ️  Tip: COMPARE_HAPS=\"0\" to focus on haplotype 0 only,"
    echo "        COMPARE_MAX_DIFFS=0 to list every skipped variant."
fi

echo "✅ Benchmark complete. Logs saved to: $LOG_FILE"