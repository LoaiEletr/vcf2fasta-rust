# User Guide for vcf2fasta-rust

## What are VCF and FASTA files?

### VCF (Variant Call Format)
Standard format for genetic variation. Key columns:
- **CHROM** – chromosome/contig name (must match reference FASTA).
- **POS** – 1-based position.
- **REF** – reference allele sequence.
- **ALT** – alternative allele(s), comma-separated.
- **FORMAT** – defines fields for each sample (e.g., `GT`).
- **Sample columns** – one per sample, containing values defined by FORMAT.

#### Genotype (GT) field
- Alleles: REF = `0`, ALT alleles = `1`, `2`, ...
- Ploidy: number of alleles (e.g., `0/1` diploid, `0|1|1` triploid).
- Phased (`|`) vs unphased (`/`).
- Missing allele: `.` (e.g., `./.`).

### FASTA (Reference Sequence)
Simple format for nucleotide sequences. Each contig starts with `>name` followed by sequence lines.
Must be indexed with `samtools faidx` (creates `.fai`).

---

## Input requirements

- VCF must be **bgzipped** (`.vcf.gz`) and **tabix-indexed** (`.tbi` or `.csi`).
- Reference FASTA must have a `.fai` index.
- Contigs in VCF must exist in the reference FASTA; others are skipped with a warning.

Prepare input:

```bash
bgzip input.vcf
tabix -p vcf input.vcf.gz
samtools faidx reference.fa
```

---

## Command-line options

| Option | Description |
|--------|-------------|
| `-f, --reference <FASTA>` | Indexed reference FASTA (`.fai` required). |
| `-p, --prefix <PREFIX>` | Prefix for output files (default: `empty`). |
| `-n, --no-call-string <STRING>` | Placeholder for missing/invalid alleles (default: `N`). |
| `-t, --threads <N>` | Number of threads (default: `1`). |
| `-w, --line-width <W>` | FASTA line width (default: `80`). |
| `-v, --no-validate-ref` | Skip REF vs FASTA validation. |
| `-q, --quiet` | Suppress stderr; only log file. |
| `--device <auto\|cpu\|gpu>` | Execution device (default: `auto`). |
| `-g, --gpu` | Deprecated alias for `--device gpu`. |
| `--bam-dir <DIR>` | Directory with per-sample BAM/CRAM files for read-based phasing. |
| `--beagle-ref-panel <DIR>` | Directory with Beagle reference panel files. |
| `--beagle-genetic-map <DIR>` | Directory with Beagle genetic map files. |
| `--chunk-size <BP>` | Override automatic chunk size (base pairs). |
| `--chunk-pad <BP>` | Look-back distance for overlapping variants (auto-detected). |
| `--phasing-threads <N>` | Threads for phasing subprocesses (default: `--threads`). |
| `--vcf2fasta-workers <N>` | Number of consumer threads (default: `1`). |
| `--ready-queue-depth <N>` | Depth of ready queue between phasing and vcf2fasta (default: `4`). |
| `--no-pipeline` | Disable pipelined scheduler (legacy mode). |
| `--merged-output` | One merged FASTA per (sample, haplotype) instead of per-contig. |
| `--max-memory <SIZE>` | Soft ceiling on process RAM (e.g., `32G`, `512M`). Clamped to a safe fraction of machine RAM. |
| `--max-vram <SIZE>` | Soft ceiling on GPU VRAM. Accepts `K`/`M`/`G`/`T` suffixes (e.g., `8G`). Acts as a **cap only**: it cannot raise the VRAM budget above the safe fraction of free VRAM the scheduler derives automatically. If the requested value exceeds that safe fraction, the effective value is logged. Default: unset (auto). |
| `--gpu-devices <LIST>` | Comma-separated CUDA device indices (e.g., `0,1`). When omitted, the scheduler uses the first visible GPU. Duplicate indices are collapsed. Requesting a device that is not visible to the process is a hard error. Multi-GPU distributes GPU worker streams round-robin across the selected devices. |

---

## Output files

### Per-contig mode (default)

One file per `(sample, contig, haplotype)`:

```
<prefix><sample>_<contig>:<hap>.fa
```

Example:

```
out_NA12877_chr22:0.fa
out_NA12877_chr22:1.fa
out_NA12877_chrX:0.fa
out_NA12877_chrX:1.fa
```

Each file contains one record:

```
>NA12877_chr22:0
ACGT...(wrapped at --line-width)
```

### Merged mode (`--merged-output`)

One **multi-record** FASTA per `(sample, haplotype)`:

```
<prefix><sample>_<hap>.fa
```

Example:

```
out_NA12877_0.fa
out_NA12877_1.fa
```

Each file contains one record per contig, in VCF-header order:

```
>NA12877_chr22:0
ACGT...(wrapped at --line-width)
...
>NA12877_chrX:0
GCTA...(wrapped at --line-width)
...
```

A merged file is exactly the byte-for-byte concatenation of the
corresponding per-contig files, in VCF-header order. Extract a single
contig with `samtools faidx <merged> <sample>_<contig>:<hap>`.

### Logs

- Log: `<prefix>.log` (or `vcf2fasta.log` if no prefix).
- Warnings log: `<prefix>.warnings.log`.

---

## Phasing behaviour

- **Already phased** → no phasing.
- **Single sample, no BAM** → canonical REF|ALT ordering.
- **Uniform diploid cohort** → Beagle.
- **Polyploid or single-sample with BAM** → WhatsHap.
- **Ploidy varies within a sample** → per-run WhatsHap + merge.

Supply `--bam-dir` if read-based phasing is required. Beagle resources can be supplied via `--beagle-ref-panel` and `--beagle-genetic-map`. Beagle itself is located via the `BEAGLE_JAR` environment variable (see [Installation Guide](INSTALLATION.md)).

---

## Edge case handling (policy table)

The example below assumes a diploid sample with genotype context `0|1` unless the situation itself specifies a different GT, and a reference base `A` at the position with `ALT = G`.

| # | Situation | Policy | Effect on output FASTA |
|---|-----------|--------|------------------------|
| 1 | Invalid POS (`-1`, `0`, `ABC`) | Skip variant + warning | Reference base `A` stays at that position; variant not applied |
| 2 | POS > FASTA contig length | Skip variant + warning | Reference base `A` stays at that position; variant not applied |
| 3 | Invalid allele index, e.g. `0\|3` (only indices 0 and 1 exist) | N only for invalid haplotype | H0 = `A` (REF), H1 = `N` |
| 4 | Negative allele index, e.g. `-1\|0` | N only for invalid haplotype | H0 = `N`, H1 = `A` (REF) |
| 5 | Too few mandatory VCF columns (< 9) | Abort the run | No output files produced |
| 6 | Missing VCF header | Abort the run | No output files produced |
| 7 | Header present but VCF has no variant records | Continue | Each output file is the reference sequence verbatim |
| 8 | VCF contig absent from reference FASTA | Skip that contig + warning | Skipped contig produces no files; other contigs still complete |
| 9 | Invalid REF nucleotide (e.g. `X`) | Skip variant + warning | Reference base `A` stays at that position; variant not applied |
| 10 | REF does not match FASTA | Skip variant + warning | Reference base `A` stays at that position; variant not applied |
| 11 | Selected ALT is invalid (e.g. `1`) | N for affected haplotype + warning | H0 = `A`, H1 = `N`; other samples unaffected |
| 12 | Unselected ALT is invalid | Ignore | No effect on output |
| 13 | Unsorted VCF | Reject input | Tool refuses to start; no output produced |
| 14 | Header declares 2 samples, record has 1 | First sample gets variant; missing sample writes REF | S1 = `A`/`G` per GT; S2 = `A` at that position |
| 15 | Header declares 1 sample, record has 2 | First sample gets variant; extra column discarded | S1 = `A`/`G` per GT; second column ignored |
| 16 | Malformed GT (e.g. `ABC`) | Skip variant for that sample | Affected sample gets `A` on both haplotypes; other samples unaffected |
| 17 | Invalid GT separator (e.g. `0&1`) | Skip variant for that sample | Affected sample gets `A` on both haplotypes; other samples unaffected |
| 18 | Missing genotype `.\|.` | Write placeholder | H0 = `N`, H1 = `N` |
| 19 | Partially missing `.\|0` | Write placeholder, then REF | H0 = `N`, H1 = `A` |
| 20 | Partially missing `0\|.` | Write REF, then placeholder | H0 = `A`, H1 = `N` |
| 21 | ALT=`.` selected by GT (e.g. `0\|1`) | N for affected haplotype + warning | H0 = `A`, H1 = `N` |
| 22 | ALT=`.` not selected by GT (`0\|0`) | Ignore | H0 = `A`, H1 = `A` (no change) |
| 23 | REF=`.` | Skip variant + warning | Reference base `A` stays at that position; variant not applied |
| 24 | REF=`.` and ALT=`.` | Skip variant + warning | Reference base `A` stays at that position; variant not applied |
| 25 | Duplicate position | Keep first, skip subsequent + warning | First record's alleles written; subsequent record has no effect |
| 26 | Overlapping variant | Keep first, skip subsequent + warning | First record's alleles written; subsequent record has no effect |
| 27 | Symbolic ALT selected (e.g. `<DEL>`) | N for affected haplotype + warning | H0 = `A`, H1 = `N` |
| 28 | Symbolic ALT present but not selected | Ignore | No effect; valid selected alleles still written normally |

---

## Log file interpretation

- `<prefix>.log`: stage-tagged lines with progress and summaries.
- `<prefix>.warnings.log`: every individual warning message.
- stderr: live feed of warnings (capped at 100 lines unless `--quiet`).

---

## Performance tips

- CPU: set `--threads` to number of physical cores.
- GPU: works best with many variants and many samples; small datasets may be slower.
- Use `--no-validate-ref` to skip reference validation if you trust your VCF.
- Use `--merged-output` to reduce file count.
- Adjust `--max-memory` to constrain RAM usage.
- **Multi-GPU**: `--gpu-devices 0,1` distributes GPU work across both
  cards. This only helps when a single contig produces enough tiles to
  keep every stream busy — cohorts of a few hundred haplotypes on a
  large chromosome scale well; small cohorts or short contigs do not.
  Multi-GPU does **not** parallelise phasing, and does **not** process
  more than one contig at a time.
- **VRAM cap**: `--max-vram 8G` lowers the GPU budget used for tile
  sizing and stream count. It cannot raise the budget above the
  scheduler's safe-VRAM heuristic. If the requested value exceeds what
  the scheduler can safely grant, the effective value is logged under
  `[SCHEDULER]`.

---

## Troubleshooting

- **“could not open indexed reference”** – ensure `.fai` exists.
- **“contig not found in tabix index”** – VCF must be bgzipped and indexed with `tabix -p vcf`.
- **“GPU not available”** – either CUDA not installed, no GPU, or build without `--features cuda`.
- **`--gpu-devices` fails with "not visible to this process"** – one or
  more requested CUDA indices are not enumerated by the driver. Run
  `nvidia-smi` (or the binary's startup log block, tagged `[SCHEDULER]
  gpu[N]:`) to see which indices are actually visible, and adjust the
  list.
- **Empty output files** – check warnings; variants might have been skipped entirely.
- **Phasing fails** – ensure BAMs are indexed and named correctly (`<SAMPLE>.bam` or with suffixes). Use `--bam-dir`.
- **Beagle not found** – set `BEAGLE_JAR` (see [Installation Guide](INSTALLATION.md)).

---