#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Locate benchmark directory
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(dirname "$SCRIPT_DIR")"

###############################################################################
# Create directory structure
###############################################################################

mkdir -p \
    "$BENCH_DIR/datasets/reference" \
    "$BENCH_DIR/datasets/platinum/NA12877" \
    "$BENCH_DIR/datasets/platinum/NA12878" \
    "$BENCH_DIR/datasets/1000G" \
    "$BENCH_DIR/datasets/genetic_map" \
    "$BENCH_DIR/chunked_inputs" \
    "$BENCH_DIR/errors" \
    "$BENCH_DIR/benchmark_results"

###############################################################################
# Reference genome
###############################################################################

echo
echo "============================================================"
echo "Downloading GENCODE GRCh38 reference"
echo "============================================================"

curl -L --fail \
    -o "$BENCH_DIR/datasets/reference/GRCh38.primary_assembly.genome.fa.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz"

echo "Decompressing reference FASTA..."

gunzip -fk \
    "$BENCH_DIR/datasets/reference/GRCh38.primary_assembly.genome.fa.gz"

echo "Indexing reference FASTA..."

samtools faidx \
    "$BENCH_DIR/datasets/reference/GRCh38.primary_assembly.genome.fa"

###############################################################################
# Platinum Genomes - NA12877
###############################################################################

echo
echo "============================================================"
echo "Downloading Platinum Genomes NA12877"
echo "============================================================"

curl -L --fail \
    -o "$BENCH_DIR/datasets/platinum/NA12877/NA12877.vcf.gz" \
    "https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12877/NA12877.vcf.gz"

curl -L --fail \
    -o "$BENCH_DIR/datasets/platinum/NA12877/NA12877.vcf.gz.tbi" \
    "https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12877/NA12877.vcf.gz.tbi"

###############################################################################
# Platinum Genomes - NA12878
###############################################################################

echo
echo "============================================================"
echo "Downloading Platinum Genomes NA12878"
echo "============================================================"

curl -L --fail \
    -o "$BENCH_DIR/datasets/platinum/NA12878/NA12878.vcf.gz" \
    "https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz"

curl -L --fail \
    -o "$BENCH_DIR/datasets/platinum/NA12878/NA12878.vcf.gz.tbi" \
    "https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz.tbi"

###############################################################################
# 1000 Genomes - Chr22
###############################################################################

echo
echo "============================================================"
echo "Downloading 1000 Genomes Chr22"
echo "============================================================"

curl -L --fail \
    -o "$BENCH_DIR/datasets/1000G/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

curl -L --fail \
    -o "$BENCH_DIR/datasets/1000G/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"

###############################################################################
# Genetic Maps (Beagle GRCh38)
###############################################################################

echo
echo "============================================================"
echo "Downloading Genetic Maps (GRCh38)"
echo "============================================================"

curl -L --fail \
    -o "$BENCH_DIR/datasets/genetic_map/plink.GRCh38.map.zip" \
    "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"

echo "Unzipping genetic maps into datasets/genetic_map/ ..."

unzip -o "$BENCH_DIR/datasets/genetic_map/plink.GRCh38.map.zip" \
    -d "$BENCH_DIR/datasets/genetic_map/"

# Remove the zip archive after extraction
rm -f "$BENCH_DIR/datasets/genetic_map/plink.GRCh38.map.zip"

###############################################################################
# Final checks
###############################################################################

echo
echo "============================================================"
echo "Checking downloaded files"
echo "============================================================"

FILES=(
    "$BENCH_DIR/datasets/reference/GRCh38.primary_assembly.genome.fa"
    "$BENCH_DIR/datasets/reference/GRCh38.primary_assembly.genome.fa.fai"

    "$BENCH_DIR/datasets/platinum/NA12877/NA12877.vcf.gz"
    "$BENCH_DIR/datasets/platinum/NA12877/NA12877.vcf.gz.tbi"

    "$BENCH_DIR/datasets/platinum/NA12878/NA12878.vcf.gz"
    "$BENCH_DIR/datasets/platinum/NA12878/NA12878.vcf.gz.tbi"

    "$BENCH_DIR/datasets/1000G/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    "$BENCH_DIR/datasets/1000G/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"

    # Genetic map: check for a representative file (chr1)
    "$BENCH_DIR/datasets/genetic_map/plink.chr1.GRCh38.map"
)

for file in "${FILES[@]}"; do

    if [[ ! -f "$file" ]]; then
        echo "ERROR: missing:"
        echo "  $file"
        exit 1
    fi

    echo "OK: $file"

done

###############################################################################
# Final message
###############################################################################

echo
echo "============================================================"
echo "Data download and indexing complete."
echo "============================================================"
echo
echo "Benchmark directory:"
echo "  $BENCH_DIR"
echo
echo "Platinum:"
echo "  $BENCH_DIR/datasets/platinum/NA12877/NA12877.vcf.gz"
echo "  $BENCH_DIR/datasets/platinum/NA12878/NA12878.vcf.gz"
echo
echo "Reference:"
echo "  $BENCH_DIR/datasets/reference/GRCh38.primary_assembly.genome.fa"
echo
echo "1000 Genomes:"
echo "  $BENCH_DIR/datasets/1000G/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
echo
echo "Genetic Maps:"
echo "  $BENCH_DIR/datasets/genetic_map/"
echo
echo "Results:"
echo "  $BENCH_DIR/benchmark_results/"
echo