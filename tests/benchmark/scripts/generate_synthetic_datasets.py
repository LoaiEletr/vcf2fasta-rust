#!/usr/bin/env python3

"""
Synthetic vcf2fasta benchmark dataset generator.

Generates:

    synthetic_unphased/
        tiny/
            reference/
                synthetic_reference.fa
                synthetic_reference.fa.fai
                synthetic_reference.fa.amb
                synthetic_reference.fa.ann
                synthetic_reference.fa.bwt
                synthetic_reference.fa.pac
                synthetic_reference.fa.sa

            vcf/
                truth_phased.vcf.gz
                truth_phased.vcf.gz.tbi
                input_mixed_phase.vcf.gz
                input_mixed_phase.vcf.gz.tbi

            truth_haplotypes/
                SAMPLE001.hap1.fa
                SAMPLE001.hap2.fa
                ...

            reads/
                ...

            bam/
                SAMPLE001.sorted.bam
                SAMPLE001.sorted.bam.bai
                ...

            metadata.json

Important:
    - Source VCF is never modified.
    - Source GRCh38 reference is never modified.
    - Generated VCF uses original GRCh38 coordinates.
    - Generated reference starts at chromosome position 1.
    - Synthetic genotypes are deterministic.
    - Truth VCF is fully phased.
    - Input VCF contains intentionally mixed phased/unphased GTs.
    - BAM reads are generated from truth haplotypes.
    - BAMs are aligned BACK to the synthetic chromosome reference,
      making them suitable for WhatsHap.
    - Both samtools and BWA indexes are created for the synthetic
      reference.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import shutil
import subprocess
import sys
import traceback

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


# ============================================================
# PATHS
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent

BENCHMARK_ROOT = SCRIPT_DIR.parent

PROJECT_ROOT = BENCHMARK_ROOT.parent.parent

DEFAULT_DATASET_ROOT = (
    BENCHMARK_ROOT / "datasets"
)

DEFAULT_SOURCE_VCF = (
    DEFAULT_DATASET_ROOT
    / "hybrid"
    / "hg38.hybrid.vcf.gz"
)

# YOUR REFERENCE IS .fa, NOT .fa.gz
DEFAULT_SOURCE_REFERENCE = (
    DEFAULT_DATASET_ROOT
    / "reference"
    / "GRCh38.primary_assembly.genome.fa"
)

DEFAULT_OUTPUT_ROOT = (
    DEFAULT_DATASET_ROOT
    / "synthetic_unphased"
)


# ============================================================
# DATA STRUCTURES
# ============================================================

@dataclass(frozen=True)
class Variant:
    chromosome: str
    pos: int
    ref: str
    alt: str
    source_id: str = "."


@dataclass(frozen=True)
class Region:
    chromosome: str
    start: int
    end: int
    variant_count: int


# ============================================================
# DATASET CONFIGURATION
# ============================================================

DATASET_CONFIG = {

    "tiny": {
        "samples": 4,
        "chromosomes": [
            "chr1",
            "chr2",
            "chr3",
        ],
        "variants_per_chromosome": 1000,
    },

    "small": {
        "samples": 8,
        "chromosomes": [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
        ],
        "variants_per_chromosome": 1500,
    },

    "medium": {
        "samples": 50,
        "chromosomes": [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
        ],
        "variants_per_chromosome": 2200,
    },

    "larger": {
        "samples": 100,
        "chromosomes": [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
        ],
        "variants_per_chromosome": 3000,
    },

    "ploidy_edge": {
        "samples": 4,
        "chromosomes": [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
        ],
        "variants_per_chromosome": 1200,
    },
}


# ============================================================
# LOGGING
# ============================================================

def log(message: str = "") -> None:
    print(message, flush=True)


# ============================================================
# COMMAND HELPERS
# ============================================================

def require_command(command: str) -> None:

    if shutil.which(command) is None:
        raise RuntimeError(
            f"Required command not found: {command}"
        )


def run_command(
    command: List[str],
    *,
    cwd: Path | None = None,
    capture_stdout: bool = True,
) -> subprocess.CompletedProcess:

    log(
        "$ " + " ".join(map(str, command))
    )

    result = subprocess.run(
        command,
        cwd=str(cwd) if cwd else None,
        stdout=(
            subprocess.PIPE
            if capture_stdout
            else None
        ),
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:

        if result.stdout:
            log(result.stdout)

        if result.stderr:
            log(result.stderr)

        raise RuntimeError(
            f"Command failed with exit code "
            f"{result.returncode}: "
            f"{' '.join(command)}"
        )

    return result


# ============================================================
# FILE HELPERS
# ============================================================

def ensure_parent(path: Path) -> None:

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


def remove_if_exists(path: Path) -> None:

    if path.exists():

        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def open_text(path: Path, mode: str = "rt"):

    if str(path).endswith(".gz"):

        return gzip.open(
            path,
            mode,
            encoding="utf-8",
        )

    return open(
        path,
        mode,
        encoding="utf-8",
    )


# ============================================================
# REFERENCE
# ============================================================

def validate_reference_index(
    reference: Path,
) -> None:

    fai = Path(
        str(reference) + ".fai"
    )

    if not fai.exists():

        log(
            f"Reference index missing: "
            f"{fai}"
        )

        run_command(
            [
                "samtools",
                "faidx",
                str(reference),
            ]
        )


def load_reference_lengths(
    reference: Path,
) -> Dict[str, int]:

    validate_reference_index(
        reference
    )

    fai = Path(
        str(reference) + ".fai"
    )

    lengths = {}

    with open(
        fai,
        "r",
        encoding="utf-8",
    ) as handle:

        for line in handle:

            fields = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(fields) < 2:
                continue

            lengths[fields[0]] = int(
                fields[1]
            )

    if not lengths:
        raise RuntimeError(
            f"No contigs found in {fai}"
        )

    return lengths


# ============================================================
# SOURCE VCF
# ============================================================

def parse_source_vcf(
    source_vcf: Path,
    requested_chromosomes: List[str],
) -> Dict[str, List[Variant]]:

    log(
        "Loading and validating source variants..."
    )

    chromosome_set = set(
        requested_chromosomes
    )

    variants = {
        chromosome: []
        for chromosome in requested_chromosomes
    }

    previous_position = {}

    result = run_command(
        [
            "bcftools",
            "view",
            "-H",
            str(source_vcf),
        ]
    )

    for line in result.stdout.splitlines():

        fields = line.split("\t")

        if len(fields) < 5:
            continue

        chromosome = fields[0]

        if chromosome not in chromosome_set:
            continue

        try:
            pos = int(fields[1])
        except ValueError:
            continue

        ref = fields[3]
        alt_field = fields[4]

        if not ref or ref == ".":
            continue

        if not alt_field or alt_field == ".":
            continue

        # Only biallelic variants.
        if "," in alt_field:
            continue

        alt = alt_field

        # Skip symbolic variants.
        if (
            alt.startswith("<")
            or alt.endswith(">")
            or "[" in alt
            or "]" in alt
            or "*" in alt
        ):
            continue

        valid_bases = set(
            "ACGTN"
        )

        if any(
            base not in valid_bases
            for base in ref.upper()
        ):
            continue

        if any(
            base not in valid_bases
            for base in alt.upper()
        ):
            continue

        if pos <= 0:
            continue

        previous = previous_position.get(
            chromosome,
            0,
        )

        if pos <= previous:
            continue

        source_id = fields[2]

        variants[chromosome].append(
            Variant(
                chromosome=chromosome,
                pos=pos,
                ref=ref.upper(),
                alt=alt.upper(),
                source_id=source_id,
            )
        )

        previous_position[
            chromosome
        ] = pos

    for chromosome in requested_chromosomes:

        count = len(
            variants[chromosome]
        )

        if count == 0:
            raise RuntimeError(
                f"No valid variants found for "
                f"{chromosome}"
            )

        log(
            f"{chromosome}: "
            f"{count:,} valid "
            f"non-overlapping variants available"
        )

    return variants


# ============================================================
# REGION SELECTION
# ============================================================

def select_variant_region(
    chromosome: str,
    variants: List[Variant],
    target_variants: int,
) -> Tuple[
    Region,
    List[Variant],
]:

    if len(variants) < target_variants:

        raise RuntimeError(
            f"{chromosome} contains only "
            f"{len(variants):,} valid variants; "
            f"{target_variants:,} required."
        )

    selected = variants[
        :target_variants
    ]

    # --------------------------------------------------------
    # BUG FIX:
    #
    # The region must extend to cover the full REF allele of
    # every selected variant, not merely the POS of the last
    # variant. An indel at the last position can have REF
    # longer than 1 base, and even an earlier variant can
    # reach further than the last variant's POS.
    #
    # Compute the true maximum END coordinate over all
    # selected variants:
    #
    #     end = POS + len(REF) - 1
    #
    # This is the last 1-based genomic base that any selected
    # variant touches, and is exactly the coordinate the
    # synthetic reference must be truncated to.
    # --------------------------------------------------------

    region_start = selected[0].pos

    region_end = max(
        v.pos + len(v.ref) - 1
        for v in selected
    )

    region = Region(
        chromosome=chromosome,
        start=region_start,
        end=region_end,
        variant_count=len(selected),
    )

    return region, selected


# ============================================================
# PLOIDY
# ============================================================

def create_ploidy_segments(
    dataset: str,
    chromosome: str,
    samples: List[str],
) -> Dict[
    str,
    List[Tuple[int, int, int]]
]:

    del chromosome

    segments = {}

    for sample in samples:

        segments[sample] = [
            (
                1,
                10**12,
                2,
            )
        ]

    if dataset != "ploidy_edge":
        return segments

    if "SAMPLE001" in segments:

        segments["SAMPLE001"] = [
            (
                1,
                500_000,
                2,
            ),
            (
                500_001,
                10**12,
                4,
            ),
        ]

    if "SAMPLE002" in segments:

        segments["SAMPLE002"] = [
            (
                1,
                500_000,
                4,
            ),
            (
                500_001,
                10**12,
                2,
            ),
        ]

    return segments


def ploidy_at(
    segments: List[
        Tuple[int, int, int]
    ],
    position: int,
) -> int:

    for start, end, ploidy in segments:

        if start <= position <= end:
            return ploidy

    raise RuntimeError(
        f"No ploidy state covers "
        f"position {position}"
    )


def create_dataset_ploidy(
    dataset: str,
    chromosomes: List[str],
    samples: List[str],
) -> Dict[
    str,
    Dict[
        str,
        List[Tuple[int, int, int]]
    ]
]:

    result = {}

    for chromosome in chromosomes:

        result[chromosome] = (
            create_ploidy_segments(
                dataset,
                chromosome,
                samples,
            )
        )

    if dataset == "ploidy_edge":

        if "chr1" in result:

            result["chr1"][
                "SAMPLE001"
            ] = [
                (
                    1,
                    500_000,
                    2,
                ),
                (
                    500_001,
                    10**12,
                    4,
                ),
            ]

            result["chr1"][
                "SAMPLE002"
            ] = [
                (
                    1,
                    500_000,
                    4,
                ),
                (
                    500_001,
                    10**12,
                    2,
                ),
            ]

        if "chr2" in result:

            result["chr2"][
                "SAMPLE003"
            ] = [
                (
                    1,
                    500_000,
                    2,
                ),
                (
                    500_001,
                    10**12,
                    3,
                ),
            ]

            result["chr2"][
                "SAMPLE004"
            ] = [
                (
                    1,
                    500_000,
                    3,
                ),
                (
                    500_001,
                    10**12,
                    4,
                ),
            ]

    return result


def validate_variants_against_ploidy(
    variants_by_chrom,
    ploidy_map,
) -> None:

    for chromosome, variants in (
        variants_by_chrom.items()
    ):

        for sample, segments in (
            ploidy_map[
                chromosome
            ].items()
        ):

            for variant in variants:

                start_ploidy = ploidy_at(
                    segments,
                    variant.pos,
                )

                end_position = (
                    variant.pos
                    + len(variant.ref)
                    - 1
                )

                end_ploidy = ploidy_at(
                    segments,
                    end_position,
                )

                if (
                    start_ploidy
                    != end_ploidy
                ):

                    raise RuntimeError(
                        "Variant crosses "
                        "ploidy boundary: "
                        f"{sample} "
                        f"{chromosome}:"
                        f"{variant.pos}"
                    )


# ============================================================
# GENOTYPES
# ============================================================

def deterministic_genotype(
    rng: random.Random,
    ploidy: int,
) -> Tuple[int, ...]:

    if ploidy not in (
        2,
        3,
        4,
    ):
        raise RuntimeError(
            f"Unsupported ploidy: "
            f"{ploidy}"
        )

    alleles = [
        1 if rng.random() < 0.30
        else 0
        for _ in range(ploidy)
    ]

    # Ensure heterozygosity occurs.
    if all(
        allele == 0
        for allele in alleles
    ):
        alleles[0] = 1

    if all(
        allele == 1
        for allele in alleles
    ):
        alleles[-1] = 0

    return tuple(alleles)


def generate_genotypes(
    dataset: str,
    samples: List[str],
    variants_by_chrom,
    ploidy_map,
    seed: int,
):

    genotype_map = {}
    phase_map = {}

    rng = random.Random(seed)

    for chromosome, variants in (
        variants_by_chrom.items()
    ):

        for variant_index, variant in enumerate(
            variants
        ):

            for sample_index, sample in enumerate(
                samples,
                start=1,
            ):

                ploidy = ploidy_at(
                    ploidy_map[
                        chromosome
                    ][sample],
                    variant.pos,
                )

                genotype = (
                    deterministic_genotype(
                        rng,
                        ploidy,
                    )
                )

                key = (
                    sample_index,
                    chromosome,
                    variant_index,
                )

                genotype_map[
                    key
                ] = genotype

                # ------------------------------------------------
                # Phase design
                # ------------------------------------------------

                if dataset == "tiny":

                    if chromosome == "chr2":
                        phased = True

                    elif chromosome == "chr3":
                        phased = False

                    elif chromosome == "chr1":
                        phased = (
                            sample_index
                            in (1, 2)
                        )

                    else:
                        phased = False

                else:

                    try:

                        chromosome_number = int(
                            chromosome.replace(
                                "chr",
                                "",
                            )
                        )

                    except ValueError:

                        chromosome_number = 1

                    # Even = phased
                    # Odd = unphased
                    phased = (
                        chromosome_number % 2
                        == 0
                    )

                    # chr1 sample-level variation
                    if chromosome == "chr1":

                        phased = (
                            sample_index
                            in (1, 2)
                        )

                if dataset == "ploidy_edge":

                    if chromosome == "chr2":

                        phased = (
                            sample_index % 2
                            == 0
                        )

                phase_map[
                    key
                ] = phased

    return (
        genotype_map,
        phase_map,
    )


# ============================================================
# VCF
# ============================================================

def make_vcf_header(
    chromosomes,
    chromosome_lengths,
    samples,
) -> List[str]:

    lines = [
        "##fileformat=VCFv4.3",
        "##source=synthetic_vcf2fasta_benchmark_generator",
        "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Synthetic source dataset\">",
        "##INFO=<ID=SOURCE_POS,Number=1,Type=Integer,Description=\"Original source VCF genomic position\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    ]

    for chromosome in chromosomes:

        lines.append(
            f"##contig=<ID={chromosome},"
            f"length={chromosome_lengths[chromosome]}>"
        )

    header = [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
    ]

    header.extend(samples)

    lines.append(
        "\t".join(header)
    )

    return lines


def genotype_to_string(
    genotype,
    phased,
) -> str:

    separator = (
        "|"
        if phased
        else "/"
    )

    return separator.join(
        str(allele)
        for allele in genotype
    )


def write_vcf(
    path: Path,
    chromosomes,
    chromosome_lengths,
    samples,
    variants_by_chrom,
    genotype_map,
    phase_map,
    truth: bool,
) -> None:

    ensure_parent(path)

    if not str(path).endswith(
        ".vcf.gz"
    ):

        raise RuntimeError(
            f"VCF path must end with "
            f".vcf.gz: {path}"
        )

    uncompressed = Path(
        str(path)[:-3]
    )

    remove_if_exists(
        uncompressed
    )

    remove_if_exists(path)

    remove_if_exists(
        Path(str(path) + ".tbi")
    )

    remove_if_exists(
        Path(str(path) + ".csi")
    )

    log(
        f"  Writing temporary VCF: "
        f"{uncompressed.name}"
    )

    header = make_vcf_header(
        chromosomes,
        chromosome_lengths,
        samples,
    )

    with open(
        uncompressed,
        "w",
        encoding="utf-8",
        newline="\n",
    ) as out:

        for line in header:

            out.write(
                line + "\n"
            )

        for chromosome in chromosomes:

            variants = variants_by_chrom[
                chromosome
            ]

            previous_pos = 0

            for variant_index, variant in enumerate(
                variants
            ):

                if variant.pos <= previous_pos:

                    raise RuntimeError(
                        f"Variants are not sorted "
                        f"on {chromosome}"
                    )

                previous_pos = variant.pos

                samples_out = []

                for sample_index, sample in enumerate(
                    samples,
                    start=1,
                ):

                    key = (
                        sample_index,
                        chromosome,
                        variant_index,
                    )

                    genotype = genotype_map[
                        key
                    ]

                    if truth:

                        phased = True

                    else:

                        phased = bool(
                            phase_map[key]
                        )

                    samples_out.append(
                        genotype_to_string(
                            genotype,
                            phased,
                        )
                    )

                info = (
                    "SOURCE=PG;"
                    f"SOURCE_POS={variant.pos}"
                )

                row = [
                    chromosome,
                    str(variant.pos),
                    variant.source_id,
                    variant.ref,
                    variant.alt,
                    "0",
                    "PASS",
                    info,
                    "GT",
                ]

                row.extend(
                    samples_out
                )

                out.write(
                    "\t".join(row)
                    + "\n"
                )

    log(
        f"  Compressing VCF: "
        f"{path.name}"
    )

    run_command(
        [
            "bgzip",
            "-f",
            str(uncompressed),
        ]
    )

    if not path.exists():

        raise RuntimeError(
            f"bgzip failed to create "
            f"{path}"
        )

    log(
        f"  Indexing VCF: "
        f"{path.name}"
    )

    run_command(
        [
            "tabix",
            "-f",
            "-p",
            "vcf",
            str(path),
        ]
    )

    if not Path(
        str(path) + ".tbi"
    ).exists():

        raise RuntimeError(
            f"VCF index missing for "
            f"{path}"
        )


# ============================================================
# REFERENCE SUBSET
# ============================================================

def create_reference_subset(
    source_reference: Path,
    output_reference: Path,
    regions: Dict[str, Region],
    source_reference_lengths: Dict[str, int],
) -> None:

    """
    Keep coordinates compatible with original GRCh38.

    For each chromosome we extract:

        chr:start=1 -> region.end

    where region.end is the last 1-based base touched by any
    selected variant (POS + len(REF) - 1), not merely the POS
    of the last variant. This guarantees that no variant's REF
    allele extends past the end of the synthetic reference.

    Then:

        >chr1:1-region.end

    is normalized to:

        >chr1

    After creating the reference, both:

        samtools faidx

    and:

        bwa index

    are generated.

    BWA indexing is required because the generated BAMs are
    aligned against this synthetic reference using bwa mem.
    """

    log(
        "Creating dataset-specific "
        "reference subset..."
    )

    ensure_parent(
        output_reference
    )

    remove_if_exists(
        output_reference
    )

    # Remove old indexes if they exist.
    remove_if_exists(
        Path(str(output_reference) + ".fai")
    )

    remove_if_exists(
        Path(str(output_reference) + ".amb")
    )

    remove_if_exists(
        Path(str(output_reference) + ".ann")
    )

    remove_if_exists(
        Path(str(output_reference) + ".bwt")
    )

    remove_if_exists(
        Path(str(output_reference) + ".pac")
    )

    remove_if_exists(
        Path(str(output_reference) + ".sa")
    )

    with open(
        output_reference,
        "w",
        encoding="utf-8",
        newline="\n",
    ) as out:

        for chromosome, region in (
            regions.items()
        ):

            reference_length = (
                source_reference_lengths[
                    chromosome
                ]
            )

            if region.end > reference_length:

                raise RuntimeError(
                    f"Region exceeds reference "
                    f"length for {chromosome}: "
                    f"need 1-{region.end:,}, "
                    f"but {chromosome} is only "
                    f"{reference_length:,} bp"
                )

            log(
                f"  {chromosome}: "
                f"1-{region.end:,} bp "
                f"(selected variants "
                f"{region.start:,}-"
                f"{region.end:,}; "
                f"region covers last REF end)"
            )

            result = run_command(
                [
                    "samtools",
                    "faidx",
                    str(source_reference),
                    f"{chromosome}:1-{region.end}",
                ]
            )

            lines = (
                result.stdout.splitlines()
            )

            if not lines:
                raise RuntimeError(
                    f"No FASTA returned for "
                    f"{chromosome}"
                )

            if not lines[0].startswith(">"):

                raise RuntimeError(
                    f"Invalid FASTA returned for "
                    f"{chromosome}"
                )

            sequence = "".join(
                line.strip()
                for line in lines[1:]
                if line.strip()
            ).upper()

            if not sequence:

                raise RuntimeError(
                    f"Empty reference sequence "
                    f"for {chromosome}"
                )

            if len(sequence) != region.end:

                raise RuntimeError(
                    f"Reference length mismatch "
                    f"for {chromosome}: "
                    f"expected {region.end:,}, "
                    f"got {len(sequence):,}"
                )

            # IMPORTANT:
            #
            # samtools faidx returns a header such as:
            #
            #   >chr1:1-region.end
            #
            # but our VCF uses:
            #
            #   chr1
            #
            # Therefore we explicitly normalize the header.
            out.write(
                f">{chromosome}\n"
            )

            for start in range(
                0,
                len(sequence),
                80,
            ):

                out.write(
                    sequence[
                        start:start + 80
                    ]
                    + "\n"
                )

    if not output_reference.exists():

        raise RuntimeError(
            f"Generated reference was not created: "
            f"{output_reference}"
        )

    if output_reference.stat().st_size == 0:

        raise RuntimeError(
            f"Generated reference is empty: "
            f"{output_reference}"
        )

    # --------------------------------------------------------
    # SAMTOOLS FASTA INDEX
    # --------------------------------------------------------

    log(
        "  Building samtools FASTA index..."
    )

    run_command(
        [
            "samtools",
            "faidx",
            str(output_reference),
        ]
    )

    fai = Path(
        str(output_reference) + ".fai"
    )

    if not fai.exists():

        raise RuntimeError(
            "Generated reference index "
            "was not created."
        )

    # --------------------------------------------------------
    # Validate samtools index.
    # --------------------------------------------------------

    indexed = set()

    with open(
        fai,
        "r",
        encoding="utf-8",
    ) as handle:

        for line in handle:

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if fields:

                indexed.add(
                    fields[0]
                )

    missing = [
        chromosome
        for chromosome in regions
        if chromosome not in indexed
    ]

    if missing:

        raise RuntimeError(
            "Generated reference missing "
            f"chromosomes in .fai: {missing}"
        )

    # --------------------------------------------------------
    # BWA INDEX
    # --------------------------------------------------------

    log(
        "  Building BWA index..."
    )

    run_command(
        [
            "bwa",
            "index",
            str(output_reference),
        ]
    )

    # --------------------------------------------------------
    # Validate all expected BWA index files.
    # --------------------------------------------------------

    bwa_index_files = [
        Path(
            str(output_reference)
            + ".amb"
        ),
        Path(
            str(output_reference)
            + ".ann"
        ),
        Path(
            str(output_reference)
            + ".bwt"
        ),
        Path(
            str(output_reference)
            + ".pac"
        ),
        Path(
            str(output_reference)
            + ".sa"
        ),
    ]

    missing_bwa = [
        str(path)
        for path in bwa_index_files
        if not path.exists()
    ]

    if missing_bwa:

        raise RuntimeError(
            "BWA index generation failed. "
            "Missing files:\n"
            + "\n".join(missing_bwa)
        )

    # --------------------------------------------------------
    # Validate chromosome retrieval.
    # --------------------------------------------------------

    log(
        "  Validating synthetic reference..."
    )

    for chromosome, region in (
        regions.items()
    ):

        result = run_command(
            [
                "samtools",
                "faidx",
                str(output_reference),
                chromosome,
            ]
        )

        lines = (
            result.stdout.splitlines()
        )

        if not lines:

            raise RuntimeError(
                f"Failed to retrieve "
                f"{chromosome} from generated "
                f"reference."
            )

        expected_header = (
            f">{chromosome}"
        )

        if lines[0].strip() != expected_header:

            raise RuntimeError(
                f"Unexpected FASTA header for "
                f"{chromosome}: "
                f"{lines[0]!r}; expected "
                f"{expected_header!r}"
            )

        sequence = "".join(
            line.strip()
            for line in lines[1:]
            if line.strip()
        )

        if len(sequence) != region.end:

            raise RuntimeError(
                f"{chromosome}: generated "
                f"reference validation failed: "
                f"expected {region.end:,} bp, "
                f"got {len(sequence):,} bp"
            )

    log(
        "  Synthetic reference validation: OK"
    )

    log(
        "  samtools FASTA index: OK"
    )

    log(
        "  BWA index: OK"
    )


# ============================================================
# REFERENCE READING
# ============================================================

def read_fasta_sequence(
    fasta: Path,
    chromosome: str,
) -> str:

    result = run_command(
        [
            "samtools",
            "faidx",
            str(fasta),
            chromosome,
        ]
    )

    lines = (
        result.stdout.splitlines()
    )

    if not lines:

        raise RuntimeError(
            f"No sequence returned for "
            f"{chromosome}"
        )

    sequence = "".join(
        line.strip()
        for line in lines
        if line.strip()
        and not line.startswith(">")
    ).upper()

    if not sequence:

        raise RuntimeError(
            f"Empty sequence for "
            f"{chromosome}"
        )

    return sequence


# ============================================================
# HAPLOTYPE CONSTRUCTION
# ============================================================

def build_haplotype_sequence(
    reference_sequence: str,
    variants: List[Variant],
    genotypes: List[Tuple[int, ...]],
    haplotype_index: int,
) -> str:

    operations = []

    for variant, genotype in zip(
        variants,
        genotypes,
    ):

        # This haplotype doesn't exist
        # for lower-ploidy segments.
        if (
            haplotype_index
            >= len(genotype)
        ):
            continue

        allele = genotype[
            haplotype_index
        ]

        if allele == 0:
            continue

        if allele != 1:

            raise RuntimeError(
                f"Unsupported allele "
                f"{allele}"
            )

        start = (
            variant.pos - 1
        )

        end = (
            start
            + len(variant.ref)
        )

        if start < 0:

            raise RuntimeError(
                f"Invalid variant position "
                f"{variant.pos}"
            )

        if end > len(
            reference_sequence
        ):

            # ----------------------------------------------------
            # Improved diagnostic: show exactly what was wrong so
            # a future mismatch is easy to chase down.
            # ----------------------------------------------------
            raise RuntimeError(
                f"Variant exceeds reference: "
                f"{variant.chromosome}:{variant.pos} "
                f"(REF={variant.ref}, "
                f"needs bases "
                f"{variant.pos}-"
                f"{variant.pos + len(variant.ref) - 1}, "
                f"but synthetic reference for "
                f"{variant.chromosome} is only "
                f"{len(reference_sequence):,} bp)"
            )

        observed = (
            reference_sequence[
                start:end
            ].upper()
        )

        if observed != variant.ref:

            raise RuntimeError(
                "Reference mismatch at "
                f"{variant.chromosome}:"
                f"{variant.pos}: "
                f"VCF={variant.ref}, "
                f"reference={observed}"
            )

        operations.append(
            (
                start,
                end,
                variant.alt,
            )
        )

    sequence = list(
        reference_sequence
    )

    for start, end, alt in reversed(
        operations
    ):

        sequence[
            start:end
        ] = list(alt)

    return "".join(sequence)


# ============================================================
# TRUTH HAPLOTYPES
# ============================================================

def create_truth_haplotypes(
    reference: Path,
    output_dir: Path,
    chromosomes: List[str],
    samples: List[str],
    variants_by_chrom,
    genotype_map,
) -> Dict[
    str,
    List[Path]
]:

    log(
        "Generating truth haplotype FASTA..."
    )

    haplotype_dir = (
        output_dir
        / "truth_haplotypes"
    )

    remove_if_exists(
        haplotype_dir
    )

    haplotype_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Cache reference sequences.
    reference_sequences = {}

    for chromosome in chromosomes:

        reference_sequences[
            chromosome
        ] = read_fasta_sequence(
            reference,
            chromosome,
        )

    result = {}

    for sample_index, sample in enumerate(
        samples,
        start=1,
    ):

        # Find maximum ploidy.
        max_ploidy = 2

        for chromosome in chromosomes:

            for variant_index in range(
                len(
                    variants_by_chrom[
                        chromosome
                    ]
                )
            ):

                genotype = genotype_map[
                    (
                        sample_index,
                        chromosome,
                        variant_index,
                    )
                ]

                max_ploidy = max(
                    max_ploidy,
                    len(genotype),
                )

        sample_paths = []

        for haplotype_index in range(
            max_ploidy
        ):

            path = (
                haplotype_dir
                / f"{sample}.hap"
                f"{haplotype_index + 1}.fa"
            )

            with open(
                path,
                "w",
                encoding="utf-8",
            ) as out:

                for chromosome in chromosomes:

                    variants = (
                        variants_by_chrom[
                            chromosome
                        ]
                    )

                    genotypes = []

                    for variant_index in range(
                        len(variants)
                    ):

                        genotypes.append(
                            genotype_map[
                                (
                                    sample_index,
                                    chromosome,
                                    variant_index,
                                )
                            ]
                        )

                    sequence = (
                        build_haplotype_sequence(
                            reference_sequences[
                                chromosome
                            ],
                            variants,
                            genotypes,
                            haplotype_index,
                        )
                    )

                    out.write(
                        f">{chromosome}\n"
                    )

                    for start in range(
                        0,
                        len(sequence),
                        80,
                    ):

                        out.write(
                            sequence[
                                start:start + 80
                            ]
                            + "\n"
                        )

            sample_paths.append(
                path
            )

        result[sample] = (
            sample_paths
        )

    return result


# ============================================================
# READ SIMULATION
# ============================================================

def simulate_reads(
    haplotypes: Dict[
        str,
        List[Path]
    ],
    output_dir: Path,
    coverage: float,
) -> Dict[
    str,
    Tuple[
        List[Path],
        List[Path],
    ]
]:

    """
    Generate paired-end reads independently from every
    truth haplotype.

    Reads are NOT aligned here.

    They are later combined per sample and aligned against
    the synthetic chromosome reference.
    """

    log(
        "Simulating paired-end reads..."
    )

    reads_dir = (
        output_dir / "reads"
    )

    reads_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    result = {}

    for sample, sample_haplotypes in (
        haplotypes.items()
    ):

        sample_r1 = []
        sample_r2 = []

        for hap_index, haplotype in enumerate(
            sample_haplotypes,
            start=1,
        ):

            prefix = (
                reads_dir
                / f"{sample}.hap"
                f"{hap_index}."
            )

            run_command(
                [
                    "art_illumina",
                    "-ss",
                    "HS25",
                    "-i",
                    str(haplotype),
                    "-p",
                    "-l",
                    "150",
                    "-m",
                    "350",
                    "-s",
                    "30",
                    "-f",
                    str(coverage),
                    "-o",
                    str(prefix),
                ]
            )

            r1 = Path(
                str(prefix) + "1.fq"
            )

            r2 = Path(
                str(prefix) + "2.fq"
            )

            if not r1.exists():

                raise RuntimeError(
                    f"ART did not create "
                    f"{r1}"
                )

            if not r2.exists():

                raise RuntimeError(
                    f"ART did not create "
                    f"{r2}"
                )

            sample_r1.append(r1)
            sample_r2.append(r2)

        result[sample] = (
            sample_r1,
            sample_r2,
        )

    return result


# ============================================================
# FASTQ CONCATENATION
# ============================================================

def concatenate_fastq(
    inputs: List[Path],
    output: Path,
) -> None:

    ensure_parent(output)

    with open(
        output,
        "wb",
    ) as out:

        for path in inputs:

            with open(
                path,
                "rb",
            ) as inp:

                shutil.copyfileobj(
                    inp,
                    out,
                )


# ============================================================
# BAM GENERATION
# ============================================================

def generate_bams(
    reference: Path,
    simulated_reads,
    output_dir: Path,
    threads: int,
) -> Dict:

    """
    Align reads against the synthetic chromosome reference.

    This is important:

        Reads originate from truth haplotypes.

        BUT

        BAM alignment is against:

            synthetic_reference.fa

        whose contigs are:

            chr1
            chr2
            chr3
            ...

    Therefore BAMs are suitable for WhatsHap.
    """

    log(
        "Aligning reads and generating BAM files..."
    )

    bam_dir = (
        output_dir / "bam"
    )

    bam_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    statistics = {}

    # --------------------------------------------------------
    # Safety check: BWA index must exist before any sample
    # is aligned.
    # --------------------------------------------------------

    required_bwa_index_files = [
        Path(
            str(reference) + ".amb"
        ),
        Path(
            str(reference) + ".ann"
        ),
        Path(
            str(reference) + ".bwt"
        ),
        Path(
            str(reference) + ".pac"
        ),
        Path(
            str(reference) + ".sa"
        ),
    ]

    missing_bwa_index = [
        str(path)
        for path in required_bwa_index_files
        if not path.exists()
    ]

    if missing_bwa_index:

        raise RuntimeError(
            "Cannot generate BAM files because "
            "the BWA index is missing for:\n"
            + "\n".join(missing_bwa_index)
        )

    for sample, (
        r1_files,
        r2_files,
    ) in simulated_reads.items():

        log(
            f"  Aligning {sample}..."
        )

        combined_r1 = (
            output_dir
            / "reads"
            / f"{sample}.R1.combined.fq"
        )

        combined_r2 = (
            output_dir
            / "reads"
            / f"{sample}.R2.combined.fq"
        )

        concatenate_fastq(
            r1_files,
            combined_r1,
        )

        concatenate_fastq(
            r2_files,
            combined_r2,
        )

        sam_path = (
            output_dir
            / "reads"
            / f"{sample}.sam"
        )

        sorted_bam = (
            bam_dir
            / f"{sample}.sorted.bam"
        )

        # ----------------------------------------------------
        # BWA is run EXACTLY ONCE.
        # ----------------------------------------------------

        with open(
            sam_path,
            "w",
            encoding="utf-8",
        ) as sam_out:

            bwa = subprocess.run(
                [
                    "bwa",
                    "mem",
                    "-t",
                    str(threads),
                    str(reference),
                    str(combined_r1),
                    str(combined_r2),
                ],
                stdout=sam_out,
                stderr=subprocess.PIPE,
                text=True,
            )

        if bwa.returncode != 0:

            raise RuntimeError(
                f"BWA failed for {sample}:\n"
                f"{bwa.stderr}"
            )

        # ----------------------------------------------------
        # Sort BAM.
        # ----------------------------------------------------

        run_command(
            [
                "samtools",
                "sort",
                "-@",
                str(threads),
                "-o",
                str(sorted_bam),
                str(sam_path),
            ]
        )

        # ----------------------------------------------------
        # Index BAM.
        # ----------------------------------------------------

        run_command(
            [
                "samtools",
                "index",
                "-@",
                str(threads),
                str(sorted_bam),
            ]
        )

        # ----------------------------------------------------
        # Validate BAM.
        # ----------------------------------------------------

        run_command(
            [
                "samtools",
                "quickcheck",
                "-v",
                str(sorted_bam),
            ]
        )

        flagstat = run_command(
            [
                "samtools",
                "flagstat",
                str(sorted_bam),
            ]
        )

        statistics[sample] = {
            "bam": str(sorted_bam),
            "flagstat": flagstat.stdout,
        }

        # Temporary files.
        remove_if_exists(
            sam_path
        )

        remove_if_exists(
            combined_r1
        )

        remove_if_exists(
            combined_r2
        )

    return statistics


# ============================================================
# VALIDATION
# ============================================================

def validate_vcf(
    vcf: Path,
    reference: Path,
) -> None:

    log(
        f"Validating VCF: "
        f"{vcf.name}"
    )

    # Check that it can be read.
    run_command(
        [
            "bcftools",
            "view",
            "-h",
            str(vcf),
        ]
    )

    # Validate REF alleles against reference.
    run_command(
        [
            "bcftools",
            "norm",
            "-f",
            str(reference),
            "-c",
            "e",
            str(vcf),
        ]
    )


def validate_bams(
    bam_dir: Path,
) -> None:

    log(
        "Validating BAM files..."
    )

    bam_files = sorted(
        bam_dir.glob(
            "*.sorted.bam"
        )
    )

    if not bam_files:

        raise RuntimeError(
            "No BAM files generated."
        )

    for bam in bam_files:

        bai = Path(
            str(bam) + ".bai"
        )

        if not bai.exists():

            raise RuntimeError(
                f"Missing BAM index: "
                f"{bai}"
            )

        run_command(
            [
                "samtools",
                "quickcheck",
                "-v",
                str(bam),
            ]
        )


# ============================================================
# METADATA
# ============================================================

def write_metadata(
    path: Path,
    dataset: str,
    seed: int,
    threads: int,
    coverage: float,
    samples: List[str],
    chromosomes: List[str],
    regions,
    variants_by_chrom,
    source_vcf: Path,
    source_reference: Path,
    bam_statistics,
) -> None:

    metadata = {

        "generator": {
            "name": (
                "synthetic_vcf2fasta_"
                "benchmark_generator"
            ),
            "version": "3.0.0",
        },

        "dataset": dataset,

        "seed": seed,

        "threads": threads,

        "coverage_per_haplotype": coverage,

        "samples": samples,

        "chromosomes": chromosomes,

        "total_variants": sum(
            len(variants)
            for variants in
            variants_by_chrom.values()
        ),

        "variants_per_chromosome": {
            chromosome: len(
                variants_by_chrom[
                    chromosome
                ]
            )
            for chromosome in chromosomes
        },

        "regions_original_grch38_coordinates": {
            chromosome: {
                "start": region.start,
                "end": region.end,
                "variant_count": (
                    region.variant_count
                ),
            }
            for chromosome, region in (
                regions.items()
            )
        },

        "source": {
            "vcf": str(source_vcf),
            "reference": str(source_reference),
            "reference_format": "FA",
            "source_files_modified": False,
        },

        "phase_design": {

            "truth": (
                "Every genotype is fully phased."
            ),

            "input": (
                "Intentionally mixed "
                "phased/unphased."
            ),

            "tiny": {
                "chr1": (
                    "SAMPLE001 and SAMPLE002 "
                    "phased; SAMPLE003 and "
                    "SAMPLE004 unphased."
                ),
                "chr2": "fully phased.",
                "chr3": "fully unphased.",
            },

            "other_datasets": (
                "Even chromosomes phased; "
                "odd chromosomes unphased; "
                "chr1 has sample-level "
                "phase variation."
            ),
        },

        "ploidy": {

            "supported": [
                2,
                3,
                4,
            ],

            "edge_cases": {
                "SAMPLE001_chr1": "P2->P4",
                "SAMPLE002_chr1": "P4->P2",
                "SAMPLE003_chr2": "P2->P3",
                "SAMPLE004_chr2": "P3->P4",
            },
        },

        "sequencing": {

            "simulator": (
                "ART Illumina"
            ),

            "read_length": 150,

            "insert_mean": 350,

            "insert_sd": 30,

            "coverage_per_haplotype": coverage,
        },

        "alignment": {

            "aligner": "BWA-MEM",

            "reference": (
                "synthetic chromosome "
                "reference"
            ),

            "sorting": (
                "samtools sort"
            ),

            "indexing": (
                "samtools index"
            ),

            "bwa_indexing": (
                "bwa index"
            ),

            "whats_hap_compatible": True,
        },

        "bam_statistics": (
            bam_statistics
        ),

        "notes": [

            "Synthetic genotypes are "
            "deterministic.",

            "No population frequency "
            "model is claimed.",

            "Truth VCF is fully phased.",

            "Input VCF contains mixed "
            "phased/unphased genotypes.",

            "Source VCF is never modified.",

            "Source GRCh38 reference "
            "is never modified.",

            "Generated VCF retains "
            "original GRCh38 coordinates.",

            "Generated reference starts "
            "at chromosome position 1.",

            "Reads are generated from "
            "truth haplotypes.",

            "Reads are aligned against "
            "the synthetic chromosome "
            "reference.",

            "BWA index is generated "
            "for the synthetic reference.",

            "BAMs are intended for "
            "WhatsHap-based phasing.",
        ],
    }

    ensure_parent(path)

    with open(
        path,
        "w",
        encoding="utf-8",
    ) as handle:

        json.dump(
            metadata,
            handle,
            indent=2,
        )

        handle.write("\n")


# ============================================================
# MAIN GENERATOR
# ============================================================

def generate_dataset(
    dataset: str,
    source_vcf: Path,
    source_reference: Path,
    output_root: Path,
    seed: int,
    threads: int,
    coverage: float,
) -> None:

    if dataset not in DATASET_CONFIG:

        raise RuntimeError(
            f"Unknown dataset: "
            f"{dataset}"
        )

    config = DATASET_CONFIG[
        dataset
    ]

    samples = [
        f"SAMPLE{i:03d}"
        for i in range(
            1,
            config["samples"] + 1,
        )
    ]

    chromosomes = (
        config["chromosomes"]
    )

    variants_per_chromosome = (
        config[
            "variants_per_chromosome"
        ]
    )

    output_dir = (
        output_root / dataset
    )

    if output_dir.exists():

        log(
            f"Removing existing dataset: "
            f"{output_dir}"
        )

        shutil.rmtree(
            output_dir
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not source_vcf.exists():

        raise RuntimeError(
            f"Source VCF not found: "
            f"{source_vcf}"
        )

    if not source_reference.exists():

        raise RuntimeError(
            f"Source reference not found: "
            f"{source_reference}"
        )

    # --------------------------------------------------------
    # Required tools.
    # --------------------------------------------------------

    required = [
        "bcftools",
        "samtools",
        "bgzip",
        "tabix",
        "bwa",
        "art_illumina",
    ]

    for command in required:
        require_command(command)

    # --------------------------------------------------------
    # Reference.
    # --------------------------------------------------------

    reference_lengths = (
        load_reference_lengths(
            source_reference
        )
    )

    missing = [
        chromosome
        for chromosome in chromosomes
        if chromosome
        not in reference_lengths
    ]

    if missing:

        raise RuntimeError(
            "Reference is missing: "
            + ", ".join(missing)
        )

    # --------------------------------------------------------
    # Source variants.
    # --------------------------------------------------------

    source_variants = (
        parse_source_vcf(
            source_vcf,
            chromosomes,
        )
    )

    # --------------------------------------------------------
    # Dynamic regions.
    # --------------------------------------------------------

    log()
    log(
        "Selecting dynamic "
        "variant-containing regions..."
    )

    regions = {}

    variants_by_chrom = {}

    for chromosome in chromosomes:

        region, selected = (
            select_variant_region(
                chromosome,
                source_variants[
                    chromosome
                ],
                variants_per_chromosome,
            )
        )

        regions[
            chromosome
        ] = region

        variants_by_chrom[
            chromosome
        ] = selected

        log(
            f"{chromosome}: "
            f"selected genomic span: "
            f"{region.start:,}-"
            f"{region.end:,}"
        )

        log(
            f"    selected variants: "
            f"{region.variant_count:,}"
        )

    # --------------------------------------------------------
    # Ploidy.
    # --------------------------------------------------------

    log()
    log(
        "Creating ploidy segment definitions..."
    )

    ploidy_map = (
        create_dataset_ploidy(
            dataset,
            chromosomes,
            samples,
        )
    )

    log(
        "Filtering variants that cross "
        "ploidy boundaries..."
    )

    validate_variants_against_ploidy(
        variants_by_chrom,
        ploidy_map,
    )

    # --------------------------------------------------------
    # Reference.
    # --------------------------------------------------------

    generated_reference = (
        output_dir
        / "reference"
        / "synthetic_reference.fa"
    )

    create_reference_subset(
        source_reference,
        generated_reference,
        regions,
        reference_lengths,
    )

    # --------------------------------------------------------
    # Genotypes.
    # --------------------------------------------------------

    log()

    log(
        "Generating synthetic genotypes..."
    )

    genotype_map, phase_map = (
        generate_genotypes(
            dataset,
            samples,
            variants_by_chrom,
            ploidy_map,
            seed,
        )
    )

    # --------------------------------------------------------
    # Truth VCF.
    # --------------------------------------------------------

    truth_vcf = (
        output_dir
        / "vcf"
        / "truth_phased.vcf.gz"
    )

    log()

    log(
        "Writing fully phased truth VCF..."
    )

    write_vcf(
        truth_vcf,
        chromosomes,
        reference_lengths,
        samples,
        variants_by_chrom,
        genotype_map,
        phase_map,
        True,
    )

    # --------------------------------------------------------
    # Mixed input VCF.
    # --------------------------------------------------------

    input_vcf = (
        output_dir
        / "vcf"
        / "input_mixed_phase.vcf.gz"
    )

    log()

    log(
        "Writing mixed phased/unphased "
        "input VCF..."
    )

    write_vcf(
        input_vcf,
        chromosomes,
        reference_lengths,
        samples,
        variants_by_chrom,
        genotype_map,
        phase_map,
        False,
    )

    # --------------------------------------------------------
    # Truth haplotypes.
    # --------------------------------------------------------

    haplotypes = (
        create_truth_haplotypes(
            generated_reference,
            output_dir,
            chromosomes,
            samples,
            variants_by_chrom,
            genotype_map,
        )
    )

    # --------------------------------------------------------
    # Reads.
    # --------------------------------------------------------

    simulated_reads = (
        simulate_reads(
            haplotypes,
            output_dir,
            coverage,
        )
    )

    # --------------------------------------------------------
    # BAM.
    # --------------------------------------------------------

    bam_statistics = (
        generate_bams(
            generated_reference,
            simulated_reads,
            output_dir,
            threads,
        )
    )

    # --------------------------------------------------------
    # Validation.
    # --------------------------------------------------------

    log()

    log(
        "Validating generated VCFs..."
    )

    validate_vcf(
        truth_vcf,
        generated_reference,
    )

    validate_vcf(
        input_vcf,
        generated_reference,
    )

    validate_bams(
        output_dir / "bam"
    )

    # --------------------------------------------------------
    # Metadata.
    # --------------------------------------------------------

    metadata = (
        output_dir
        / "metadata.json"
    )

    write_metadata(
        metadata,
        dataset,
        seed,
        threads,
        coverage,
        samples,
        chromosomes,
        regions,
        variants_by_chrom,
        source_vcf,
        source_reference,
        bam_statistics,
    )

    # --------------------------------------------------------
    # Final summary.
    # --------------------------------------------------------

    log()
    log("=" * 60)
    log(
        f"DATASET '{dataset}' "
        f"GENERATED SUCCESSFULLY"
    )
    log("=" * 60)

    log(
        f"Output: {output_dir}"
    )

    log(
        f"Samples: {len(samples)}"
    )

    log(
        f"Chromosomes: "
        f"{', '.join(chromosomes)}"
    )

    log(
        f"Variants/chromosome: "
        f"{variants_per_chromosome:,}"
    )

    log(
        f"Total variants: "
        f"{variants_per_chromosome * len(chromosomes):,}"
    )

    log(
        f"Reference: "
        f"{generated_reference}"
    )

    log(
        f"Truth VCF: "
        f"{truth_vcf}"
    )

    log(
        f"Input VCF: "
        f"{input_vcf}"
    )

    log(
        f"Metadata: "
        f"{metadata}"
    )

    log("=" * 60)


# ============================================================
# CLI
# ============================================================

def parse_path(
    value: str,
) -> Path:

    path = Path(
        value
    ).expanduser()

    if path.is_absolute():
        return path

    return (
        BENCHMARK_ROOT / path
    )


def build_parser():

    parser = argparse.ArgumentParser(
        description=(
            "Generate synthetic "
            "vcf2fasta benchmark datasets."
        )
    )

    parser.add_argument(
        "--dataset",
        choices=sorted(
            DATASET_CONFIG.keys()
        ),
        default="tiny",
    )

    parser.add_argument(
        "--source-vcf",
        type=Path,
        default=DEFAULT_SOURCE_VCF,
    )

    parser.add_argument(
        "--source-reference",
        type=Path,
        default=DEFAULT_SOURCE_REFERENCE,
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--coverage",
        type=float,
        default=5.0,
    )

    return parser


# ============================================================
# MAIN
# ============================================================

def main() -> int:

    parser = build_parser()

    args = parser.parse_args()

    source_vcf = parse_path(
        str(args.source_vcf)
    )

    source_reference = parse_path(
        str(args.source_reference)
    )

    output_root = parse_path(
        str(args.output_root)
    )

    log()
    log("=" * 60)
    log(
        " Synthetic vcf2fasta "
        "Benchmark Dataset Generator"
    )
    log("=" * 60)

    log(
        f"Project root:     {PROJECT_ROOT}"
    )

    log(
        f"Benchmark root:   {BENCHMARK_ROOT}"
    )

    log(
        f"Script directory: {SCRIPT_DIR}"
    )

    log(
        f"Source VCF:       {source_vcf}"
    )

    log(
        f"Reference:        {source_reference}"
    )

    log(
        "Reference format: FA"
    )

    log(
        f"Output root:      {output_root}"
    )

    log(
        f"Dataset:          {args.dataset}"
    )

    log(
        f"Seed:             {args.seed}"
    )

    log(
        f"Threads:          {args.threads}"
    )

    log(
        f"Coverage/hap:     "
        f"{args.coverage}x"
    )

    log("=" * 60)

    try:

        if args.threads < 1:

            raise RuntimeError(
                "--threads must be >= 1"
            )

        if args.coverage <= 0:

            raise RuntimeError(
                "--coverage must be > 0"
            )

        generate_dataset(
            dataset=args.dataset,
            source_vcf=source_vcf,
            source_reference=source_reference,
            output_root=output_root,
            seed=args.seed,
            threads=args.threads,
            coverage=args.coverage,
        )

        return 0

    except Exception as exc:

        log()
        log(
            "DATASET GENERATION FAILED: "
            f"{exc}"
        )

        log()
        log(
            "Full traceback:"
        )

        traceback.print_exc()

        return 1


if __name__ == "__main__":
    sys.exit(main())