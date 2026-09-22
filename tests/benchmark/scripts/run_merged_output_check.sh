#!/usr/bin/env bash
# ============================================================================
# run_merged_output_check.sh
#
# Verifies that `--merged-output` produces a byte-correct concatenation of
# the per-contig output, in VCF-header contig order.
#
# Since the fix, `--merged-output` produces a **multi-record FASTA** per
# `(sample, hap)`: one header line + one wrapped sequence per contig, in
# canonical order. The merged file is therefore byte-identical to
# `cat out_<sample>_<contig1>:<hap>.fa out_<sample>_<contig2>:<hap>.fa ...`
# in VCF-header order.
#
# Strategy:
#   1. Slice the input VCF to the selected contigs (default: first 2 present
#      in both the VCF header and the reference .fai).
#   2. Run the tool in per-contig mode (default) → baseline FASTA files.
#   3. Run the tool with --merged-output on the same sliced VCF.
#   4. For every (sample, haplotype) pair:
#        a. Verify the merged file exists and is named <prefix><sample>_<hap>.fa.
#        b. Byte-compare the merged file to the concatenation of the
#           per-contig files, in VCF-header contig order, using `diff -q`.
#        c. Verify the merged file has exactly one header per present contig.
#   5. Report PASS/FAIL per (sample, haplotype) and an overall verdict.
#
# Exits non-zero if any check fails.
#
# Env vars (all optional):
#   INPUT_VCF         input VCF                  (default: NA12877 Platinum VCF)
#   REFERENCE         reference FASTA            (default: GRCh38 primary assembly)
#   MODE              cpu | gpu | both           (default: cpu)
#   GPU_DEVICES       comma-separated GPU ids    (default: 0)
#   THREADS           host threads               (default: 2)
#   SAMPLES           space-separated samples    (default: first sample in VCF)
#   CONTIGS           space-separated contigs    (default: first 2 in VCF header)
#   MAX_CONTIGS       cap when CONTIGS unset     (default: 2)
#   NO_VALIDATE_REF   true | false               (default: false)
#   WORKDIR           scratch directory          (default: tests/benchmark/merged_check)
# ============================================================================

set -uo pipefail

# --- CUDA library path (same convention as other benchmark scripts) ---
if [ -d /usr/local/cuda-12.8/targets/x86_64-linux/lib ]; then
    export LD_LIBRARY_PATH="/usr/local/cuda-12.8/targets/x86_64-linux/lib:/usr/local/cuda-12.8/lib64:/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
fi

# --- Paths ------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(cd "$BENCH_DIR/../.." && pwd)"
RUST_BIN="${PROJECT_ROOT}/target/release/vcf2fasta"

# --- Config -----------------------------------------------------------------
INPUT_VCF="${INPUT_VCF:-${BENCH_DIR}/datasets/platinum/NA12877/NA12877.vcf.gz}"
REFERENCE="${REFERENCE:-${BENCH_DIR}/datasets/reference/GRCh38.primary_assembly.genome.fa}"
MODE="${MODE:-cpu}"
GPU_DEVICES="${GPU_DEVICES:-0}"
THREADS="${THREADS:-2}"
NO_VALIDATE_REF="${NO_VALIDATE_REF:-false}"
WORKDIR="${WORKDIR:-${BENCH_DIR}/merged_check}"
MAX_CONTIGS="${MAX_CONTIGS:-2}"

# --- Parse MODE -------------------------------------------------------------
case "$MODE" in
    cpu)  ENGINES=("cpu") ;;
    gpu)  ENGINES=("gpu") ;;
    both) ENGINES=("cpu" "gpu") ;;
    *) echo "❌ Invalid MODE='$MODE' (expected: cpu, gpu, both)" >&2; exit 1 ;;
esac

# --- Preflight --------------------------------------------------------------
if [ ! -x "$RUST_BIN" ]; then
    echo "❌ Rust binary not found: $RUST_BIN" >&2
    echo "   Build with: cargo build --release --features cuda" >&2
    exit 1
fi

for tool in bcftools tabix; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "❌ Required tool not found on PATH: $tool" >&2
        exit 1
    fi
done

if [ ! -s "$INPUT_VCF" ]; then
    echo "❌ Input VCF not found or empty: $INPUT_VCF" >&2
    exit 1
fi
if [ ! -s "${INPUT_VCF}.tbi" ] && [ ! -s "${INPUT_VCF}.csi" ]; then
    echo "❌ Input VCF has no tabix index: $INPUT_VCF" >&2
    echo "   Run: tabix -p vcf $INPUT_VCF" >&2
    exit 1
fi
if [ ! -s "$REFERENCE" ]; then
    echo "❌ Reference not found or empty: $REFERENCE" >&2
    exit 1
fi
if [ ! -s "${REFERENCE}.fai" ]; then
    echo "❌ Reference index not found: ${REFERENCE}.fai" >&2
    echo "   Run: samtools faidx $REFERENCE" >&2
    exit 1
fi

mkdir -p "$WORKDIR"
LOG_FILE="${WORKDIR}/merged_check.log"
: > "$LOG_FILE"

# --- Banner -----------------------------------------------------------------
{
    echo "=========================================================================="
    echo "Merged Output Verification"
    echo "=========================================================================="
    echo "Binary:      $RUST_BIN"
    echo "Input VCF:   $INPUT_VCF"
    echo "Reference:   $REFERENCE"
    echo "Mode:        $MODE"
    echo "Threads:     $THREADS"
    echo "GPU devices: $GPU_DEVICES"
    echo "Workdir:     $WORKDIR"
    echo
} | tee -a "$LOG_FILE"

# --- Determine samples ------------------------------------------------------
if [ -z "${SAMPLES:-}" ]; then
    SAMPLES=$(bcftools query -l "$INPUT_VCF" | head -n1)
fi
if [ -z "$SAMPLES" ]; then
    echo "❌ No samples found in VCF header." | tee -a "$LOG_FILE"
    exit 1
fi
echo "Samples:     $SAMPLES" | tee -a "$LOG_FILE"

# --- Determine contigs ------------------------------------------------------
VCF_CONTIGS_RAW=$(bcftools view -h "$INPUT_VCF" \
    | grep '^##contig=' \
    | sed 's/^##contig=<ID=//' \
    | sed 's/[,>].*//' \
    | tr -d '"')

REF_CONTIGS=$(awk '{print $1}' "${REFERENCE}.fai")

AVAILABLE_CONTIGS=""
for C in $VCF_CONTIGS_RAW; do
    if echo "$REF_CONTIGS" | grep -qFx "$C"; then
        AVAILABLE_CONTIGS="$AVAILABLE_CONTIGS $C"
    fi
done
AVAILABLE_CONTIGS=$(echo "$AVAILABLE_CONTIGS" | tr -s ' ' | sed 's/^ //;s/ $//')

if [ -z "${CONTIGS:-}" ]; then
    CONTIGS=$(echo "$AVAILABLE_CONTIGS" | tr ' ' '\n' | head -n "$MAX_CONTIGS" | tr '\n' ' ')
fi
CONTIGS=$(echo "$CONTIGS" | tr -s ' ' | sed 's/^ //;s/ $//')

if [ -z "$CONTIGS" ]; then
    echo "❌ No usable contigs (must be present in both VCF header and reference)." \
        | tee -a "$LOG_FILE"
    exit 1
fi

echo "Contigs:     $CONTIGS" | tee -a "$LOG_FILE"

N_CONTIGS=$(echo "$CONTIGS" | wc -w)
if [ "$N_CONTIGS" -lt 2 ]; then
    {
        echo
        echo "⚠️  Only one contig selected. The merged output will be identical"
        echo "   to the per-contig output for that contig. Testing with ≥2"
        echo "   contigs is strongly recommended to exercise concatenation."
    } | tee -a "$LOG_FILE"
fi
echo | tee -a "$LOG_FILE"

# --- Slice the VCF ----------------------------------------------------------
REGIONS=$(echo "$CONTIGS" | tr ' ' ',')
SLICED_VCF="${WORKDIR}/sliced.vcf.gz"

echo "→ Slicing VCF to: $REGIONS" | tee -a "$LOG_FILE"
if ! bcftools view -r "$REGIONS" "$INPUT_VCF" -Oz -o "$SLICED_VCF" 2>>"$LOG_FILE"; then
    echo "❌ bcftools view failed (see $LOG_FILE)" | tee -a "$LOG_FILE"
    exit 1
fi
if ! tabix -f -p vcf "$SLICED_VCF" 2>>"$LOG_FILE"; then
    echo "❌ tabix failed to index $SLICED_VCF" | tee -a "$LOG_FILE"
    exit 1
fi
echo | tee -a "$LOG_FILE"

# --- Arrays for iteration ---------------------------------------------------
read -ra CONTIGS_ARR <<< "$CONTIGS"
read -ra SAMPLES_ARR <<< "$SAMPLES"

# ============================================================================
# Helpers
# ============================================================================

# Detect ploidy for (sample, first_contig) from per-contig output.
# Returns the number of consecutive haps starting at 0 whose files exist.
detect_ploidy() {
    local prefix="$1"
    local sample="$2"
    local first_contig="$3"
    local h=0
    while [ -f "${prefix}${sample}_${first_contig}:${h}.fa" ]; do
        h=$((h + 1))
    done
    echo "$h"
}

# Compare the merged file to the byte-for-byte concatenation of per-contig
# files, in canonical contig order, using `diff`. Headers are included, so
# this catches any header, ordering, or wrapping mismatch.
#
# Prints a one-line verdict to stdout; returns 0 on PASS, 1 on FAIL.
# Args: <percontig_prefix> <merged_prefix> <sample> <hap> <contigs...>
check_pair() {
    local percontig_prefix="$1"
    local merged_prefix="$2"
    local sample="$3"
    local hap="$4"
    shift 4
    local -a contigs=("$@")

    local merged_file="${merged_prefix}${sample}_${hap}.fa"
    if [ ! -f "$merged_file" ]; then
        printf 'MISSING_MERGED_FILE (looked for %s)\n' "$merged_file"
        return 1
    fi

    # Build the expected merged file by concatenating the per-contig files
    # in canonical contig order. The headers are included — the merged file
    # should equal `cat out_S_chr1:0.fa out_S_chr2:0.fa ...`.
    local tmp_expected
    tmp_expected=$(mktemp)
    local n_present=0
    local c
    for c in "${contigs[@]}"; do
        local f="${percontig_prefix}${sample}_${c}:${hap}.fa"
        if [ -f "$f" ]; then
            cat "$f" >> "$tmp_expected"
            n_present=$((n_present + 1))
        fi
    done

    local exp_size act_size
    exp_size=$(wc -c < "$tmp_expected" | tr -d ' ')
    act_size=$(wc -c < "$merged_file"   | tr -d ' ')

    if diff -q "$tmp_expected" "$merged_file" >/dev/null 2>&1; then
        rm -f "$tmp_expected"
        printf 'PASS  contigs=%d  bytes=%d\n' "$n_present" "$act_size"
        return 0
    fi

    # On mismatch, produce a small diagnostic showing the first few
    # differing lines to help debug.
    local diff_summary
    diff_summary=$(diff "$tmp_expected" "$merged_file" 2>&1 | head -n 8 || true)
    rm -f "$tmp_expected"

    printf 'MISMATCH  expected_bytes=%d  actual_bytes=%d  delta=%+d\n' \
        "$exp_size" "$act_size" "$(( act_size - exp_size ))"
    printf '          first differences:\n'
    while IFS= read -r line; do
        printf '            %s\n' "$line"
    done <<< "$diff_summary"
    return 1
}

# ============================================================================
# Main loop
# ============================================================================

OVERALL_PASS=true

for ENGINE in "${ENGINES[@]}"; do
    {
        echo "=========================================================================="
        echo "Engine: $ENGINE"
        echo "=========================================================================="
    } | tee -a "$LOG_FILE"

    PERCONTIG_DIR="${WORKDIR}/${ENGINE}_percontig"
    MERGED_DIR="${WORKDIR}/${ENGINE}_merged"
    rm -rf "$PERCONTIG_DIR" "$MERGED_DIR"
    mkdir -p "$PERCONTIG_DIR" "$MERGED_DIR"

    EXTRA=()
    if [ "$NO_VALIDATE_REF" = "true" ]; then
        EXTRA+=(--no-validate-ref)
    fi
    if [ "$ENGINE" = "gpu" ]; then
        ENGINE_ARGS=(--device gpu --gpu-devices "$GPU_DEVICES")
    else
        ENGINE_ARGS=(--device cpu)
    fi

    # --- Per-contig run ----------------------------------------------------
    echo | tee -a "$LOG_FILE"
    echo "→ Running per-contig mode ($ENGINE)..." | tee -a "$LOG_FILE"
    PERCONTIG_START=$(date +%s.%N)
    if ! "$RUST_BIN" \
            "${ENGINE_ARGS[@]}" \
            --threads "$THREADS" \
            --reference "$REFERENCE" \
            --prefix "${PERCONTIG_DIR}/run_" \
            "${EXTRA[@]}" \
            "$SLICED_VCF" \
            > "${WORKDIR}/${ENGINE}_percontig.stdout" \
            2> "${WORKDIR}/${ENGINE}_percontig.stderr"; then
        echo "  ❌ Per-contig run failed. See ${WORKDIR}/${ENGINE}_percontig.stderr" \
            | tee -a "$LOG_FILE"
        OVERALL_PASS=false
        continue
    fi
    PERCONTIG_END=$(date +%s.%N)
    PERCONTIG_TIME=$(echo "$PERCONTIG_END - $PERCONTIG_START" | bc)
    echo "  Time:   ${PERCONTIG_TIME}s" | tee -a "$LOG_FILE"
    echo "  Output: $PERCONTIG_DIR/" | tee -a "$LOG_FILE"

    # --- Merged run --------------------------------------------------------
    echo | tee -a "$LOG_FILE"
    echo "→ Running merged mode ($ENGINE)..." | tee -a "$LOG_FILE"
    MERGED_START=$(date +%s.%N)
    if ! "$RUST_BIN" \
            --merged-output \
            "${ENGINE_ARGS[@]}" \
            --threads "$THREADS" \
            --reference "$REFERENCE" \
            --prefix "${MERGED_DIR}/run_" \
            "${EXTRA[@]}" \
            "$SLICED_VCF" \
            > "${WORKDIR}/${ENGINE}_merged.stdout" \
            2> "${WORKDIR}/${ENGINE}_merged.stderr"; then
        echo "  ❌ Merged run failed. See ${WORKDIR}/${ENGINE}_merged.stderr" \
            | tee -a "$LOG_FILE"
        OVERALL_PASS=false
        continue
    fi
    MERGED_END=$(date +%s.%N)
    MERGED_TIME=$(echo "$MERGED_END - $MERGED_START" | bc)
    echo "  Time:   ${MERGED_TIME}s" | tee -a "$LOG_FILE"
    echo "  Output: $MERGED_DIR/" | tee -a "$LOG_FILE"

    # --- Compare -----------------------------------------------------------
    echo | tee -a "$LOG_FILE"
    echo "→ Comparing per-contig vs merged ($ENGINE)..." | tee -a "$LOG_FILE"
    echo | tee -a "$LOG_FILE"

    N_PAIRS=0
    N_PASS=0
    N_FAIL=0

    FIRST_CONTIG="${CONTIGS_ARR[0]}"

    for SAMPLE in "${SAMPLES_ARR[@]}"; do
        PLOIDY=$(detect_ploidy "${PERCONTIG_DIR}/run_" "$SAMPLE" "$FIRST_CONTIG")

        if [ "$PLOIDY" -eq 0 ]; then
            echo "  $SAMPLE: no per-contig output — skipping" | tee -a "$LOG_FILE"
            continue
        fi

        echo "  $SAMPLE (ploidy=$PLOIDY):" | tee -a "$LOG_FILE"

        for HAP in $(seq 0 $((PLOIDY - 1))); do
            N_PAIRS=$((N_PAIRS + 1))
            if RESULT=$(check_pair "${PERCONTIG_DIR}/run_" "${MERGED_DIR}/run_" \
                          "$SAMPLE" "$HAP" "${CONTIGS_ARR[@]}"); then
                N_PASS=$((N_PASS + 1))
                echo "    hap $HAP: ✅ $RESULT" | tee -a "$LOG_FILE"
            else
                N_FAIL=$((N_FAIL + 1))
                OVERALL_PASS=false
                echo "    hap $HAP: ❌ $RESULT" | tee -a "$LOG_FILE"
            fi
        done
    done

    # --- Structural checks -------------------------------------------------
    echo | tee -a "$LOG_FILE"
    echo "→ Structural checks ($ENGINE)..." | tee -a "$LOG_FILE"

    # Count merged files
    N_MERGED_ACTUAL=$(find "$MERGED_DIR" -maxdepth 1 -type f -name 'run_*.fa' 2>/dev/null | wc -l | tr -d ' ')
    N_MERGED_EXPECTED=0
    for SAMPLE in "${SAMPLES_ARR[@]}"; do
        P=$(detect_ploidy "${PERCONTIG_DIR}/run_" "$SAMPLE" "$FIRST_CONTIG")
        N_MERGED_EXPECTED=$((N_MERGED_EXPECTED + P))
    done

    if [ "$N_MERGED_ACTUAL" -eq "$N_MERGED_EXPECTED" ]; then
        echo "  ✅ Merged file count: $N_MERGED_ACTUAL (expected $N_MERGED_EXPECTED)" \
            | tee -a "$LOG_FILE"
    else
        echo "  ❌ Merged file count mismatch: got $N_MERGED_ACTUAL, expected $N_MERGED_EXPECTED" \
            | tee -a "$LOG_FILE"
        OVERALL_PASS=false
    fi

    # Every merged file should have exactly one header line per present
    # contig — one record per contig in a multi-record FASTA.
    N_BAD_HEADER=0
    while IFS= read -r f; do
        N_HEADERS=$(grep -c '^>' "$f" || true)
        if [ "$N_HEADERS" -ne "$N_CONTIGS" ]; then
            echo "  ❌ $(basename "$f") has $N_HEADERS header line(s) (expected $N_CONTIGS)" \
                | tee -a "$LOG_FILE"
            N_BAD_HEADER=$((N_BAD_HEADER + 1))
        fi
    done < <(find "$MERGED_DIR" -maxdepth 1 -type f -name 'run_*.fa' 2>/dev/null)
    if [ "$N_BAD_HEADER" -eq 0 ]; then
        echo "  ✅ Every merged file has exactly $N_CONTIGS header line(s)" \
            | tee -a "$LOG_FILE"
    else
        OVERALL_PASS=false
    fi

    # --- Engine summary ----------------------------------------------------
    echo | tee -a "$LOG_FILE"
    echo "  Engine $ENGINE: $N_PASS / $N_PAIRS pairs PASS, $N_FAIL FAIL" \
        | tee -a "$LOG_FILE"
    echo | tee -a "$LOG_FILE"
done

# ============================================================================
# Overall summary
# ============================================================================

{
    echo "=========================================================================="
    echo "Summary"
    echo "=========================================================================="
    echo "Log: $LOG_FILE"
} | tee -a "$LOG_FILE"

if [ "$OVERALL_PASS" = true ]; then
    echo "Result: ✅ ALL CHECKS PASSED" | tee -a "$LOG_FILE"
    exit 0
else
    echo "Result: ❌ SOME CHECKS FAILED" | tee -a "$LOG_FILE"
    exit 1
fi