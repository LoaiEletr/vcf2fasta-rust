# Architecture Guide for vcf2fasta-rust

## 1. High-level architecture

```
                    +------------------------------------------------------------------+
                    |                    Command Line Interface                        |
                    |                  (clap::Parser → Args struct)                    |
                    +------------------------------------------------------------------+
                                                  |
                                                  v
                    +------------------------------------------------------------------+
                    |                        main::run()                               |
                    |  - Parse CLI, set up logging                                     |
                    |  - Read VCF header (samples, contigs)                            |
                    |  - Verify contigs exist in reference FASTA (.fai)                |
                    |  - Discover work units (one per contig)                          |
                    |  - Build execution plan (CPU/GPU, threads, memory budget)        |
                    |  - Launch producer-consumer pipeline                             |
                    +------------------------------------------------------------------+
                                                  |
                                                  v
                    +------------------------------------------------------------------+
                    |                     Pipeline (pipeline.rs)                       |
                    |  Producer thread:                                                |
                    |    - For each contig:                                            |
                    |      - Slice VCF for contig                                      |
                    |      - Determine phasing backend (workunit.rs)                   |
                    |      - Run phasing if needed (Beagle, WhatsHap, VariablePloidy)  |
                    |      - Send WorkUnit to ready queue                              |
                    |  Consumer threads (CPU) or single GPU consumer:                  |
                    |    - Receive WorkUnit                                            |
                    |    - Plan 2-D tiles (haplotypes × base ranges)                   |
                    |    - Execute tiles on CPU (Rayon) or GPU (CUDA)                  |
                    |    - Write results via OutputManager                             |
                    +------------------------------------------------------------------+
```

---

## 2. Core data structures

### `Args` (cli.rs)
Parsed command-line flags: input, reference, prefix, no-call string, threads, line width, validation flags, quiet, device selection, BAM dir, Beagle resource dirs, chunk sizes, memory limits, GPU devices, VRAM limit, etc.

### `WorkUnit` (workunit.rs)
Represents one contig. Contains:
- `contig`, `reference_length`, `sample_count`, `sample_max_ploidies`, `haplotype_count`
- `variant_count`, `phased_genotypes`, `unphased_genotypes`, `malformed_records`
- `needs_phasing`, `phase_backend` (None, SingleSampleCanonical, Beagle, WhatsHap, VariablePloidy)
- `phased_vcf` (optional path to phased VCF), `phased_tempdir` (ownership guard)
- `state` (WorkState enum)

### `Tile` (chunk.rs)
A 2-D work unit: a haplotype block × a base range within one contig.

```rust
pub struct Tile {
    pub contig: String,
    pub hap_start: usize,
    pub hap_end: usize,
    pub base_start: u64,
    pub base_end: u64,
}
```

### `DecodedRecord` (genotype.rs)
Validated variant record with:
- `start`, `end` (0-based half-open REF interval)
- `alleles`: vector of allele sequences (bytes)
- `genotypes`: per sample, vector of `AlleleCall` per haplotype
- `ploidies`: per sample

### `AlleleCall` (genotype.rs)
- `Index(usize)` – write the allele at this index
- `Missing` – write placeholder string (default `N`)
- `Reference` – write reference base from FASTA at this position

### `ContigReport` (report.rs)
Per-contig statistics: `seen`, `applied`, `skipped`, `output_files`, `warning_count`, `warnings` (capped at 100), `warnings_by_reason` (unbounded histogram).

---

## 3. Input preparation

- VCF must be bgzipped (`.vcf.gz`) and indexed with `tabix -p vcf` (produces `.tbi` or `.csi`).
- Reference FASTA must have a `.fai` index (`samtools faidx`).
- Contigs in VCF must exist in the reference; others are skipped with a warning.

---

## 4. Genotype normalisation

```
String GT (e.g. "0|1")
        ↓
parse_gt() → Vec<GenotypeAllele> (rust-htslib)
        ↓
normalize_genotype() → Vec<Option<usize>>
        ↓
Convert to AlleleCall:
- Some(i) → Index(i) if i < alleles.len() and allele is valid
- Some(i) → Missing if i invalid (symbolic, out-of-range, etc.)
- None → Missing (for '.' alleles)
- If GT malformed or sample missing → Reference
```

Special cases:
- Unphased (`0/1`) is warned but mapped deterministically (first before `/` → H0, second after `/` → H1). For single-sample canonical mode, alleles are sorted to enforce REF|ALT order.
- Negative allele indices (`-1`) are treated as Missing.
- Symbolic ALT (`<DEL>`) selected → Missing.
- ALT `.` selected → Missing.

---

## 5. Phasing decision (workunit.rs)

The backend for a contig is chosen automatically from the data:

| Condition (per contig)                                      | Backend               |
|-------------------------------------------------------------|-----------------------|
| all GTs phased, or max ploidy ≤ 1                           | None                  |
| N = 1, no BAM                                               | SingleSampleCanonical |
| N = 1, BAM present, constant ploidy                         | WhatsHap              |
| N = 1, BAM present, ploidy varies within the sample         | VariablePloidy        |
| N ≥ 2, ploidy varies within any sample on this contig       | VariablePloidy        |
| N ≥ 2, all samples uniform diploid                          | Beagle                |
| N ≥ 2, polyploid, constant per sample                       | WhatsHap              |

- **SingleSampleCanonical** is not a phaser; it uses the executor's canonical REF|ALT ordering.
- **Beagle** runs population phasing for uniform diploid cohorts.
- **WhatsHap** uses read-based phasing (requires `--bam-dir`).
- **VariablePloidy** partitions the contig into constant-ploidy runs, phases each with WhatsHap, and merges.

---

## 6. CPU execution path

- Each WorkUnit is processed by a consumer thread.
- Tiles are planned using `plan_tiles()`: splits haplotypes into blocks and base ranges into chunks.
- Tiles are grouped into batches to bound memory.
- Each tile is executed by `execute_tile()` (tile_executor.rs):
  - Opens reference and tabix VCF.
  - Iterates variants in the tile's base range.
  - Validates and decodes records.
  - Builds per-haplotype sequences by concatenating reference segments and chosen alleles.
- Rayon is used to parallelise tiles within a batch.
- Results are written via `OutputManager`.

---

## 7. GPU execution path

- Requires `--features cuda` and `--device gpu`.
- Uses a custom CUDA kernel (`apply_variants`) to apply variants in parallel.
- One producer thread builds `GpuBatch` objects (host-side data: variant
  positions, alleles, genotype indices).
- N CUDA worker threads (streams) pull batches from a shared bounded
  channel. Each worker owns a `GpuWorker` (CUDA stream + reusable pinned
  and device scratch buffers).
- **Device selection.** The scheduler records the selected CUDA device
  indices in `ExecutionPlan::gpu_device_indices`. `CudaPipeline::spawn`
  assigns worker `i` to device `i % devices.len()`, so streams are
  distributed round-robin across the selected GPUs. Each device gets
  its own `CudaDevice` (own context, own NVRTC-compiled module) and its
  own `ReferenceCache`, because device-side reference buffers are not
  portable across CUDA devices.
- **Multi-GPU granularity.** Parallelism is tile-level *within* one
  contig. The producer submits tiles from a single contig to the shared
  channel, and whichever worker pulls a tile first processes it. There
  is no cross-contig GPU parallelism (GPU mode forces
  `--vcf2fasta-workers 1`), and a single tile is never split across
  devices. The ordered collector (`next_expected`) means throughput is
  bounded by the slowest stream for the duration of one tile.
- **VRAM budget.** `ExecutionPlan::gpu_memory_budget_bytes` is derived
  from `min(free_vram × 0.5, max_alloc_bytes)` summed across the
  selected devices, then capped by `--max-vram` if supplied. Tile size,
  stream count, and in-flight depth are all sized from this budget.
- Kernel layout:
  - 2-D grid: `(hap_blocks, variant_blocks)`.
  - Each thread handles one haplotype and a slice of variants.
  - Reference segments are copied, chosen alleles inserted, tail appended.
- Output is read back and split into per-hap byte vectors.
- Reusable pinned host buffers and device buffers minimise allocation overhead.

### Feature gating

All CUDA code is behind `#[cfg(feature = "cuda")]`:

- `src/gpu.rs`, `src/cuda_pipeline.rs` — compiled only with the feature.
- `src/lib.rs::run_gpu_tiles` — two implementations: a CUDA one behind `#[cfg(feature = "cuda")]`, and a stub that errors out behind `#[cfg(not(feature = "cuda"))]`.
- `src/scheduler.rs::detect_all_gpus` — CUDA enumeration gated by the feature; returns an empty vector otherwise.

When built without the feature, `--device gpu` fails with a clear message telling you to rebuild with `--features cuda` or use `--device cpu`.

---

## 8. Output manager

- `OutputManager` supports two modes:
  - **Per-contig (default):** one FASTA file per `(sample, contig, hap)`.
    File name: `<prefix><sample>_<contig>:<hap>.fa`. One header line,
    body wrapped at `--line-width`.
  - **Merged (`--merged-output`):** one **multi-record** FASTA per
    `(sample, hap)`. File name: `<prefix><sample>_<hap>.fa`. One
    `><sample>_<contig>:<hap>` header per contig, body wrapped at
    `--line-width`, records ordered by VCF-header contig order.
- A merged file is a **strict byte-for-byte concatenation** of the
  corresponding per-contig files, in canonical contig order:
  `cat out_S_chr1:0.fa out_S_chr2:0.fa` in VCF-header order equals
  `out_S_0.fa`. This property is verified by
  `tests/benchmark/scripts/run_merged_output_check.sh`.
- Keeps a bounded number of open file handles (default 512) using an
  LRU cache. Evicted writers are flushed and their column state
  remembered; reopened in append mode without duplicating the FASTA
  header.
- `finish_contig()` finalises all writers for a contig.
- Merged mode uses a scratch directory (a `tempfile::TempDir` by
  default, respecting `$TMPDIR`) to hold per-contig `.part` files
  between the streaming writes and the final concatenation. Each part
  is deleted immediately after concatenation, so scratch usage shrinks
  as `finish()` progresses. Concatenation uses `std::io::copy`
  (kernel-level `copy_file_range` on Linux), so peak RSS is
  independent of contig length.

---

## 9. Resource budget and scheduling

- `ResourceBudget` (resource.rs) splits CPU cores between phasing and vcf2fasta, and clamps `--max-memory` to a safe fraction of machine RAM.
- `Scheduler` (scheduler.rs) chooses CPU or GPU based on workload size, available hardware, and user request.
- Execution plan includes: device, worker counts, tile sizes, GPU
  stream count, in-flight buffers, memory budgets, and the selected
  CUDA device indices.
- `--gpu-devices 0,1` forces a multi-device plan.
  `Scheduler::plan_multi_gpu` validates each requested index against
  `HardwareInfo::gpus`, deduplicates the list preserving first-seen
  order, and errors if any index is not visible to the process. When
  the flag is not given, `Scheduler::plan` uses the first visible
  device.
- `--max-vram` caps the auto-derived GPU budget. The cap is monotone:
  it can lower the budget, never raise it above the safe fraction. The
  effective value is reported in the `[SCHEDULER]` startup log block,
  and a `[SCHEDULER]` line is emitted when the requested cap is lower
  than the auto-derived budget and therefore binding.

---

## 10. Logging and warnings

- Main log (`<prefix>.log`): stage-tagged lines (`[DISCOVERY]`, `[PHASING]`, `[READY]`, `[SCHEDULER]`, `[COMPLETE]`, `[SUMMARY]`).
- Warnings log (`<prefix>.warnings.log`): every individual warning message (capped at 1,000,000 lines).
- stderr receives a capped live feed of warnings (100 lines) unless `--quiet`.

The `[SCHEDULER]` block at startup reports:
- Detected hardware (cores, RAM, GPUs, per-GPU free/total VRAM).
- Effective device, worker counts, haplotype and variant block sizes.
- GPU stream count, in-flight depth, and (for GPU plans) the selected
  device indices and effective VRAM budget.
- A one-line `reason` string summarising the plan's derivation.

---

## 11. Dependencies

| Crate | Purpose |
|-------|---------|
| `clap` | CLI parsing |
| `anyhow` | Error handling |
| `rayon` | CPU parallelisation |
| `rust-htslib` | VCF reading, FASTA fetching |
| `crossbeam-channel` | Bounded channels for pipeline |
| `cudarc` (optional) | CUDA bindings; enabled by the `cuda` feature |
| `flate2` | Gzip decompression |
| `tempfile` | Temporary directories for phasing |
| `chrono` | Timestamps in logs |

Cargo features:

| Feature | Effect |
|---------|--------|
| (default) | CPU-only build. No CUDA dependency compiled in. |
| `cuda` | Enables `cudarc` and compiles `src/gpu.rs` and `src/cuda_pipeline.rs`. |

---

## 12. Test layout

- **Unit tests** — embedded in each `src/*.rs` module under `#[cfg(test)]`. Some are CUDA-gated.
- **Integration tests** — `tests/pipeline_integration.rs` and `tests/policy_integration.rs`. Both exercise CPU paths.
  - `pipeline_integration.rs`: `discover()` and `WorkState` classification.
  - `policy_integration.rs`: end-to-end policy behaviour (unsorted VCF rejection, duplicate/overlap handling). Skips gracefully when `bgzip` / `tabix` are missing.
- **Benchmark harness** — the shell scripts under `tests/benchmark/scripts/` and the two Python validators under `tests/benchmark/validation/`

Run CPU tests only:

```bash
cargo test --release
```

Run everything (needs a CUDA device):

```bash
cargo test --release --features cuda
```

---