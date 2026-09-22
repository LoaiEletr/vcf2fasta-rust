#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd -- "$BENCHMARK_ROOT/../.." && pwd)"

DATASET_ROOT="$BENCHMARK_ROOT/datasets"

SOURCE_VCF="${SOURCE_VCF:-$DATASET_ROOT/platinum/NA12877/NA12877.vcf.gz}"
SOURCE_REFERENCE="${SOURCE_REFERENCE:-$DATASET_ROOT/reference/GRCh38.primary_assembly.genome.fa}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$DATASET_ROOT/synthetic_unphased}"

DATASET="${DATASET:-tiny}"
SEED="${SEED:-12345}"
THREADS="${THREADS:-4}"
COVERAGE="${COVERAGE:-5}"

echo "============================================================"
echo " Synthetic vcf2fasta Benchmark Dataset Generator"
echo "============================================================"
echo "Project root:     $PROJECT_ROOT"
echo "Benchmark root:   $BENCHMARK_ROOT"
echo "Script directory: $SCRIPT_DIR"
echo "Source VCF:       $SOURCE_VCF"
echo "Reference:        $SOURCE_REFERENCE"
echo "Output root:      $OUTPUT_ROOT"
echo "Dataset:          $DATASET"
echo "Seed:             $SEED"
echo "Threads:          $THREADS"
echo "Coverage/hap:     ${COVERAGE}x"
echo "============================================================"

if [[ ! -f "$SOURCE_VCF" ]]; then
    echo "ERROR: Source VCF does not exist:"
    echo "  $SOURCE_VCF"
    exit 1
fi

if [[ ! -f "$SOURCE_REFERENCE" ]]; then
    echo "ERROR: Source reference does not exist:"
    echo "  $SOURCE_REFERENCE"
    exit 1
fi

if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: python not found."
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found."
    exit 1
fi

if ! command -v bwa >/dev/null 2>&1; then
    echo "ERROR: bwa not found."
    exit 1
fi

if ! command -v art_illumina >/dev/null 2>&1; then
    echo "ERROR: art_illumina not found."
    exit 1
fi

python "$SCRIPT_DIR/generate_synthetic_datasets.py" \
    --dataset "$DATASET" \
    --source-vcf "$SOURCE_VCF" \
    --source-reference "$SOURCE_REFERENCE" \
    --output-root "$OUTPUT_ROOT" \
    --seed "$SEED" \
    --threads "$THREADS" \
    --coverage "$COVERAGE"

echo
echo "============================================================"
echo " Dataset generation completed successfully."
echo "============================================================"