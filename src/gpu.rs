//! CUDA-accelerated variant application.
//!
//! Design notes:
//! * Reference is uploaded once per contig into a shared [`ReferenceCache`].
//! * Each worker thread owns a [`GpuWorker`] which owns a [`BatchScratch`].
//!   Staging buffers (pinned host + device) are allocated once and reused
//!   across every batch.
//! * Genotypes use a **[variant][hap]** layout so that threads in a warp
//!   (which process consecutive haplotypes for the *same* variant) read
//!   consecutive addresses. This coalesces what would otherwise be 32
//!   scattered 4-byte loads into one 128-byte transaction.
//! * The kernel is **2D-tiled**: the grid is `(hap_blocks, variant_blocks)`.
//!   Each thread handles one haplotype and a fixed slice of variants. The
//!   host precomputes the output offset for each `(variant_block, hap)`
//!   pair so the block can write directly without racing.
//!
//! ## Output layout
//!
//! The flat output buffer is **hap-major**:
//!
//! ```text
//! [ hap 0: block 0 | block 1 | ... | block N-1 ]
//! [ hap 1: block 0 | block 1 | ... | block N-1 ]
//! ...
//! ```
//!
//! The `hap_block_offsets` table passed to the kernel gives the **flat**
//! position where each `(block, hap)` pair starts writing. It is built by
//! [`GpuBatch::compute_block_offsets`], which performs two passes over the
//! variants: one to size each haplotype's region, one to record the
//! starting offset of every block in flat coordinates.

#![cfg(feature = "cuda")]

use crate::pinned::PinnedBuf;
use anyhow::{anyhow, Result};
use cudarc::driver::{
    CudaContext, CudaFunction, CudaModule, CudaSlice, CudaStream, DeviceRepr, LaunchConfig,
    PushKernelArg, ValidAsZeroBits,
};
use cudarc::nvrtc::Ptx;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use std::time::Instant;

// ---------------------------------------------------------------------------
// Tuning
// ---------------------------------------------------------------------------

/// Variants per kernel block (grid.y).
pub const GPU_VAR_BLOCK: usize = 4096;

/// Maximum threads per block.
pub const GPU_MAX_BLOCK_DIM: u32 = 1024;

/// CUDA warp size on all current NVIDIA hardware.
pub const GPU_WARP_SIZE: u32 = 32;

// ---------------------------------------------------------------------------
// Kernel
// ---------------------------------------------------------------------------

pub const KERNEL_SRC: &str = r#"
extern "C" __global__ void apply_variants(
    const unsigned char* __restrict__ contig_ref,
    int ref_start,
    int ref_end,
    const int*           __restrict__ variant_positions,
    const unsigned char* __restrict__ variant_alleles,
    const int*           __restrict__ allele_lengths,
    const int*           __restrict__ allele_byte_offsets,
    const int*           __restrict__ allele_start_idx,
    const int*           __restrict__ num_alleles,
    const int*           __restrict__ genotype_indices,   // [variant][hap]
    const int*           __restrict__ hap_block_offsets,  // [block][hap] in flat output coords
    const int*           __restrict__ block_ref_starts,   // [block]
    int num_variants,
    int num_haps,
    int var_block_size,
    unsigned char*       __restrict__ output_seq)
{
    int hap_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (hap_id >= num_haps) return;

    int block_idx = blockIdx.y;
    int v_start = block_idx * var_block_size;
    int v_end = min(v_start + var_block_size, num_variants);

    int out_pos = hap_block_offsets[block_idx * num_haps + hap_id];
    int ref_pos = block_ref_starts[block_idx];

    for (int v = v_start; v < v_end; ++v) {
        int pos = variant_positions[v];
        int allele_idx = genotype_indices[v * num_haps + hap_id];

        int copy_len = pos - ref_pos;
        for (int j = 0; j < copy_len; ++j) {
            output_seq[out_pos + j] = contig_ref[ref_pos + j];
        }
        out_pos += copy_len;

        int start = allele_start_idx[v];
        int ref_len_here = allele_lengths[start];
        int n_all = num_alleles[v];

        if (allele_idx == -1) {
            output_seq[out_pos] = contig_ref[pos];
            out_pos += 1;
        } else if (allele_idx == -2) {
            int pl_idx = start + n_all - 1;
            int pl_off = allele_byte_offsets[pl_idx];
            int pl_len = allele_lengths[pl_idx];
            for (int j = 0; j < pl_len; ++j) {
                output_seq[out_pos + j] = variant_alleles[pl_off + j];
            }
            out_pos += pl_len;
        } else {
            int a_idx = start + allele_idx;
            int a_off = allele_byte_offsets[a_idx];
            int a_len = allele_lengths[a_idx];
            for (int j = 0; j < a_len; ++j) {
                output_seq[out_pos + j] = variant_alleles[a_off + j];
            }
            out_pos += a_len;
        }
        ref_pos = pos + ref_len_here;
    }

    // Only the last variant block writes the trailing reference.
    if (block_idx == gridDim.y - 1) {
        int tail_len = ref_end - ref_pos;
        if (tail_len > 0) {
            for (int j = 0; j < tail_len; ++j) {
                output_seq[out_pos + j] = contig_ref[ref_pos + j];
            }
        }
    }
}
"#;

// ---------------------------------------------------------------------------
// Device
// ---------------------------------------------------------------------------

pub struct CudaDevice {
    pub ctx: Arc<CudaContext>,
    pub module: Arc<CudaModule>,
    pub name: String,
    pub total_vram_bytes: u64,
    pub free_vram_bytes: u64,
    pub max_alloc_bytes: u64,
}

impl CudaDevice {
    pub fn init() -> Result<Self> {
        cudarc::driver::result::init()
            .map_err(|e| anyhow!("cuInit failed: {e:?}"))?;

        let ctx = CudaContext::new(0)
            .map_err(|e| anyhow!("CUDA context creation failed: {e:?}"))?;

        let ptx: Ptx = cudarc::nvrtc::compile_ptx(KERNEL_SRC)
            .map_err(|e| anyhow!("NVRTC compilation failed: {e:?}"))?;

        let module = ctx
            .load_module(ptx)
            .map_err(|e| anyhow!("CUDA module load failed: {e:?}"))?;

        module
            .load_function("apply_variants")
            .map_err(|e| anyhow!("kernel 'apply_variants' not found: {e:?}"))?;

        let (free, total) = ctx
            .mem_get_info()
            .map_err(|e| anyhow!("cudaMemGetInfo failed: {e:?}"))?;
        let total_u64 = total as u64;
        let free_u64 = free as u64;
        let name = ctx.name().unwrap_or_else(|_| "unknown".to_string());

        Ok(Self {
            ctx,
            module,
            name,
            total_vram_bytes: total_u64,
            free_vram_bytes: free_u64,
            max_alloc_bytes: total_u64 / 2,
        })
    }

    pub fn init_device(device_index: usize) -> Result<Self> {
        cudarc::driver::result::init()
            .map_err(|e| anyhow!("cuInit failed: {e:?}"))?;

        let ctx = CudaContext::new(device_index).map_err(|e| {
            anyhow!(
                "CUDA context creation failed for device {}: {e:?}",
                device_index
            )
        })?;

        let ptx: Ptx = cudarc::nvrtc::compile_ptx(KERNEL_SRC)
            .map_err(|e| anyhow!("NVRTC compilation failed: {e:?}"))?;

        let module = ctx
            .load_module(ptx)
            .map_err(|e| anyhow!("CUDA module load failed: {e:?}"))?;

        module
            .load_function("apply_variants")
            .map_err(|e| anyhow!("kernel 'apply_variants' not found: {e:?}"))?;

        let (free, total) = ctx
            .mem_get_info()
            .map_err(|e| anyhow!("cudaMemGetInfo failed: {e:?}"))?;
        let total_u64 = total as u64;
        let free_u64 = free as u64;
        let name = ctx
            .name()
            .unwrap_or_else(|_| format!("device{}", device_index));

        Ok(Self {
            ctx,
            module,
            name,
            total_vram_bytes: total_u64,
            free_vram_bytes: free_u64,
            max_alloc_bytes: total_u64 / 2,
        })
    }

    pub fn new_stream(&self) -> Result<Arc<CudaStream>> {
        self.ctx
            .new_stream()
            .map_err(|e| anyhow!("cudaStreamCreate failed: {e:?}"))
    }

    pub fn kernel(&self) -> Result<CudaFunction> {
        self.module
            .load_function("apply_variants")
            .map_err(|e| anyhow!("kernel load failed: {e:?}"))
    }
}

pub fn cuda_available() -> bool {
    CudaDevice::init().is_ok()
}

pub fn enumerate_devices() -> Result<Vec<crate::scheduler::GpuInfo>> {
    cudarc::driver::result::init().map_err(|e| anyhow!("cuInit failed: {e:?}"))?;

    let n = cudarc::driver::result::device::get_count()
        .map_err(|e| anyhow!("cuDeviceGetCount failed: {e:?}"))?;
    let mut out = Vec::with_capacity(n as usize);
    for i in 0..n {
        let dev = match CudaDevice::init_device(i as usize) {
            Ok(d) => d,
            Err(_) => continue,
        };
        out.push(crate::scheduler::GpuInfo {
            name: dev.name.clone(),
            total_vram_bytes: dev.total_vram_bytes,
            free_vram_bytes: dev.free_vram_bytes,
            max_alloc_bytes: dev.max_alloc_bytes,
            device_index: i as usize,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Growable device buffer
// ---------------------------------------------------------------------------

pub struct GrowableDev<T: DeviceRepr + ValidAsZeroBits> {
    buf: Option<CudaSlice<T>>,
    capacity: usize,
}

impl<T: DeviceRepr + ValidAsZeroBits> GrowableDev<T> {
    pub fn new() -> Self {
        Self {
            buf: None,
            capacity: 0,
        }
    }

    pub fn ensure(
        &mut self,
        stream: &Arc<CudaStream>,
        len: usize,
    ) -> Result<&mut CudaSlice<T>> {
        let want = len.max(1);
        if self.buf.is_none() || self.capacity < want {
            self.buf = Some(
                stream
                    .alloc_zeros::<T>(want)
                    .map_err(|e| anyhow!("cudaMalloc failed: {e:?}"))?,
            );
            self.capacity = want;
        }
        Ok(self.buf.as_mut().unwrap())
    }

    pub fn as_slice_ref(&self) -> Option<&CudaSlice<T>> {
        self.buf.as_ref()
    }

    pub fn len(&self) -> usize {
        self.buf.as_ref().map(|b| b.len()).unwrap_or(0)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<T: DeviceRepr + ValidAsZeroBits> Default for GrowableDev<T> {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Reference cache
// ---------------------------------------------------------------------------

pub struct ReferenceCache {
    inner: RwLock<HashMap<String, Arc<CudaSlice<u8>>>>,
}

impl ReferenceCache {
    pub fn new() -> Self {
        Self {
            inner: RwLock::new(HashMap::new()),
        }
    }

    pub fn get_or_upload(
        &self,
        contig: &str,
        host_bytes: &[u8],
        stream: &Arc<CudaStream>,
    ) -> Result<Arc<CudaSlice<u8>>> {
        {
            let guard = self.inner.read().unwrap();
            if let Some(buf) = guard.get(contig) {
                return Ok(buf.clone());
            }
        }
        let mut guard = self.inner.write().unwrap();
        if let Some(buf) = guard.get(contig) {
            return Ok(buf.clone());
        }
        let mut dev = stream
            .alloc_zeros::<u8>(host_bytes.len().max(1))
            .map_err(|e| anyhow!("cudaMalloc(contig ref) failed: {e:?}"))?;
        if !host_bytes.is_empty() {
            stream
                .memcpy_htod(host_bytes, &mut dev)
                .map_err(|e| anyhow!("H2D(contig ref) failed: {e:?}"))?;
        }
        stream
            .synchronize()
            .map_err(|e| anyhow!("stream sync after contig upload failed: {e:?}"))?;
        let buf = Arc::new(dev);
        guard.insert(contig.to_string(), buf.clone());
        Ok(buf)
    }

    pub fn len(&self) -> usize {
        self.inner.read().unwrap().len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Default for ReferenceCache {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Batch
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct GpuBatch {
    pub contig: String,
    pub ref_start: usize,
    pub ref_end: usize,
    pub variant_positions: Vec<i32>,
    pub variant_alleles: Vec<u8>,
    pub allele_lengths: Vec<i32>,
    pub allele_byte_offsets: Vec<i32>,
    pub allele_start_idx: Vec<i32>,
    pub num_alleles: Vec<i32>,
    /// **Layout: `[variant][hap]`.** Element `(v, h)` is
    /// `genotype_indices[v * num_haps + h]`.
    pub genotype_indices: Vec<i32>,
    pub num_haps: usize,
}

impl GpuBatch {
    pub fn variant_count(&self) -> usize {
        self.variant_positions.len()
    }

    pub fn n_var_blocks(&self) -> usize {
        let nv = self.variant_count();
        if nv == 0 {
            1
        } else {
            nv.div_ceil(GPU_VAR_BLOCK)
        }
    }

    fn block_ref_starts(&self) -> Vec<i32> {
        let n_blocks = self.n_var_blocks();
        let nv = self.variant_count();
        let mut out = Vec::with_capacity(n_blocks);
        out.push(self.ref_start as i32);
        for b in 1..n_blocks {
            let v_prev = (b * GPU_VAR_BLOCK).saturating_sub(1);
            if v_prev >= nv {
                out.push(self.ref_end as i32);
                continue;
            }
            let start = self.allele_start_idx[v_prev] as usize;
            let ref_len = self.allele_lengths[start] as i32;
            let pos = self.variant_positions[v_prev];
            out.push(pos + ref_len);
        }
        out
    }

    /// Per-block output offsets and total output length per haplotype.
    ///
    /// The output buffer is laid out **hap-major**:
    ///
    /// ```text
    /// [ hap 0: block 0 | block 1 | ... | block N-1 ]
    /// [ hap 1: block 0 | block 1 | ... | block N-1 ]
    /// ...
    /// ```
    ///
    /// `offsets[b * num_haps + h]` is the **flat** position where block
    /// `b` starts writing for haplotype `h`. Row `b == n_var_blocks` is
    /// the end position of every haplotype's region (after the trailing
    /// reference), which equals the start of the next hap's region.
    ///
    /// Returns `(block_offsets, lengths)` where `lengths[h]` is the total
    /// number of bytes written for haplotype `h`.
    fn compute_block_offsets(&self) -> (Vec<i32>, Vec<usize>) {
        let nh = self.num_haps;
        let nv = self.variant_count();
        let n_blocks = self.n_var_blocks();

        if nh == 0 {
            return (Vec::new(), Vec::new());
        }

        // Helper: byte contribution of one (variant, hap) pair.
        //
        //   code == -1  → 1 byte (the reference base at `pos`)
        //   code == -2  → the placeholder, whose length is the last
        //                 allele's length in this variant's slice
        //   code >= 0   → the allele at `start + code`
        let allele_len_for = |v: usize, h: usize| -> i64 {
            let start = self.allele_start_idx[v] as usize;
            let n_all = self.num_alleles[v] as usize;
            let base = v * nh;
            let code = self.genotype_indices[base + h];
            if code == -1 {
                1
            } else if code == -2 {
                let pl_idx = start + n_all - 1;
                self.allele_lengths[pl_idx] as i64
            } else {
                self.allele_lengths[start + code as usize] as i64
            }
        };

        // ---- Pass 1: total bytes per haplotype. ----
        let mut per_hap_len = vec![0i64; nh];
        {
            let mut ref_pos = self.ref_start as i64;
            for v in 0..nv {
                let pos = self.variant_positions[v] as i64;
                let gap = pos - ref_pos;
                let start = self.allele_start_idx[v] as usize;
                let ref_len = self.allele_lengths[start] as i64;
                for h in 0..nh {
                    per_hap_len[h] += gap + allele_len_for(v, h);
                }
                ref_pos = pos + ref_len;
            }
            // Trailing reference: only the last block writes it, but every
            // hap's region reserves room for it.
            let tail = (self.ref_end as i64 - ref_pos).max(0);
            for h in 0..nh {
                per_hap_len[h] += tail;
            }
        }

        // Prefix sums: flat base offset of each hap's region.
        let mut hap_base = vec![0i64; nh];
        for h in 1..nh {
            hap_base[h] = hap_base[h - 1] + per_hap_len[h - 1];
        }

        // ---- Pass 2: per-block start positions in flat coordinates. ----
        let mut offsets = vec![0i32; (n_blocks + 1) * nh];
        {
            let mut per_hap_pos = vec![0i64; nh]; // bytes written by each hap so far
            let mut ref_pos = self.ref_start as i64;

            for block in 0..n_blocks {
                // Start of this block, expressed relative to each hap's
                // region base.
                for h in 0..nh {
                    offsets[block * nh + h] = (hap_base[h] + per_hap_pos[h]) as i32;
                }

                let v_start = block * GPU_VAR_BLOCK;
                let v_end = (v_start + GPU_VAR_BLOCK).min(nv);
                for v in v_start..v_end {
                    let pos = self.variant_positions[v] as i64;
                    let gap = pos - ref_pos;
                    let start = self.allele_start_idx[v] as usize;
                    let ref_len = self.allele_lengths[start] as i64;
                    for h in 0..nh {
                        per_hap_pos[h] += gap + allele_len_for(v, h);
                    }
                    ref_pos = pos + ref_len;
                }
            }

            // End row: after the last block, which also writes the tail.
            let tail = (self.ref_end as i64 - ref_pos).max(0);
            for h in 0..nh {
                per_hap_pos[h] += tail;
                offsets[n_blocks * nh + h] = (hap_base[h] + per_hap_pos[h]) as i32;
            }
        }

        let lengths: Vec<usize> = per_hap_len.iter().map(|&l| l as usize).collect();
        (offsets, lengths)
    }
}

// ---------------------------------------------------------------------------
// Batch stats
// ---------------------------------------------------------------------------

#[derive(Debug, Default, Clone)]
pub struct BatchStats {
    pub seen: usize,
    pub applied: usize,
    pub warning_count: usize,
    pub warnings: Vec<String>,
    pub warnings_by_reason: std::collections::BTreeMap<&'static str, usize>,
}

// ---------------------------------------------------------------------------
// Timing
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default)]
pub struct BatchTiming {
    pub h2d_ns: u64,
    pub kernel_ns: u64,
    pub d2h_ns: u64,
}

// ---------------------------------------------------------------------------
// Reusable per-worker scratch space
// ---------------------------------------------------------------------------

pub struct BatchScratch {
    pub pos_host: PinnedBuf<i32>,
    pub alleles_host: PinnedBuf<u8>,
    pub lens_host: PinnedBuf<i32>,
    pub offs_host: PinnedBuf<i32>,
    pub start_host: PinnedBuf<i32>,
    pub numal_host: PinnedBuf<i32>,
    pub geno_host: PinnedBuf<i32>,
    pub blk_off_host: PinnedBuf<i32>,
    pub blk_ref_host: PinnedBuf<i32>,
    pub out_host: PinnedBuf<u8>,
    pub dev_pos: GrowableDev<i32>,
    pub dev_alleles: GrowableDev<u8>,
    pub dev_lens: GrowableDev<i32>,
    pub dev_offs: GrowableDev<i32>,
    pub dev_start: GrowableDev<i32>,
    pub dev_numal: GrowableDev<i32>,
    pub dev_geno: GrowableDev<i32>,
    pub dev_blk_off: GrowableDev<i32>,
    pub dev_blk_ref: GrowableDev<i32>,
    pub dev_out: GrowableDev<u8>,
    pub pinned_reported: bool,
}

impl BatchScratch {
    pub fn new() -> Self {
        Self {
            pos_host: PinnedBuf::new(0),
            alleles_host: PinnedBuf::new(0),
            lens_host: PinnedBuf::new(0),
            offs_host: PinnedBuf::new(0),
            start_host: PinnedBuf::new(0),
            numal_host: PinnedBuf::new(0),
            geno_host: PinnedBuf::new(0),
            blk_off_host: PinnedBuf::new(0),
            blk_ref_host: PinnedBuf::new(0),
            out_host: PinnedBuf::new(0),
            dev_pos: GrowableDev::new(),
            dev_alleles: GrowableDev::new(),
            dev_lens: GrowableDev::new(),
            dev_offs: GrowableDev::new(),
            dev_start: GrowableDev::new(),
            dev_numal: GrowableDev::new(),
            dev_geno: GrowableDev::new(),
            dev_blk_off: GrowableDev::new(),
            dev_blk_ref: GrowableDev::new(),
            dev_out: GrowableDev::new(),
            pinned_reported: false,
        }
    }
}

impl Default for BatchScratch {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// GpuWorker
// ---------------------------------------------------------------------------

pub struct GpuWorker {
    pub device: Arc<CudaDevice>,
    pub stream: Arc<CudaStream>,
    pub ref_cache: Arc<ReferenceCache>,
    pub scratch: BatchScratch,
}

impl GpuWorker {
    pub fn new(device: Arc<CudaDevice>, ref_cache: Arc<ReferenceCache>) -> Result<Self> {
        let stream = device.new_stream()?;
        Ok(Self {
            device,
            stream,
            ref_cache,
            scratch: BatchScratch::new(),
        })
    }

    pub fn run(
        &mut self,
        batch: &GpuBatch,
        contig_ref_host: &[u8],
    ) -> Result<(Vec<Vec<u8>>, BatchTiming)> {
        let num_variants = batch.variant_count();
        let num_haps = batch.num_haps;
        let n_blocks = batch.n_var_blocks();
        let mut timing = BatchTiming::default();

        if num_haps == 0 {
            return Ok((Vec::new(), timing));
        }

        let dev_ref = self
            .ref_cache
            .get_or_upload(&batch.contig, contig_ref_host, &self.stream)?;

        let (block_offsets, lengths) = batch.compute_block_offsets();
        let block_ref_starts = batch.block_ref_starts();
        let total_out: usize = lengths.iter().sum();

        self.scratch.pos_host.ensure_capacity(num_variants.max(1));
        self.scratch
            .alleles_host
            .ensure_capacity(batch.variant_alleles.len().max(1));
        self.scratch
            .lens_host
            .ensure_capacity(batch.allele_lengths.len().max(1));
        self.scratch
            .offs_host
            .ensure_capacity(batch.allele_byte_offsets.len().max(1));
        self.scratch
            .start_host
            .ensure_capacity(num_variants.max(1));
        self.scratch
            .numal_host
            .ensure_capacity(num_variants.max(1));
        self.scratch
            .geno_host
            .ensure_capacity((num_haps * num_variants).max(1));
        self.scratch
            .blk_off_host
            .ensure_capacity(block_offsets.len().max(1));
        self.scratch
            .blk_ref_host
            .ensure_capacity(block_ref_starts.len().max(1));
        self.scratch.out_host.ensure_capacity(total_out.max(1));

        if !self.scratch.pinned_reported {
            self.scratch.pinned_reported = true;
            if !self.scratch.out_host.is_pinned() {
                eprintln!(
                    "WARNING: D2H staging buffer is NOT pinned. \
                     This typically costs 5-10x in D2H throughput. \
                     Check `ulimit -l` and container memory-locking permissions."
                );
            }
        }

        if num_variants > 0 {
            self.scratch.pos_host[..num_variants]
                .copy_from_slice(&batch.variant_positions);
            self.scratch.start_host[..num_variants]
                .copy_from_slice(&batch.allele_start_idx);
            self.scratch.numal_host[..num_variants]
                .copy_from_slice(&batch.num_alleles);
        }
        if !batch.variant_alleles.is_empty() {
            let n = batch.variant_alleles.len();
            self.scratch.alleles_host[..n].copy_from_slice(&batch.variant_alleles);
        }
        if !batch.allele_lengths.is_empty() {
            let n = batch.allele_lengths.len();
            self.scratch.lens_host[..n].copy_from_slice(&batch.allele_lengths);
            self.scratch.offs_host[..n].copy_from_slice(&batch.allele_byte_offsets);
        }
        if !batch.genotype_indices.is_empty() {
            let n = batch.genotype_indices.len();
            self.scratch.geno_host[..n].copy_from_slice(&batch.genotype_indices);
        }
        self.scratch.blk_off_host[..block_offsets.len()].copy_from_slice(&block_offsets);
        self.scratch.blk_ref_host[..block_ref_starts.len()].copy_from_slice(&block_ref_starts);

        let t0 = Instant::now();

        macro_rules! upload {
            ($dev:ident, $host:ident, $n:expr) => {{
                let dst = self.scratch.$dev.ensure(&self.stream, $n.max(1))?;
                if $n > 0 {
                    self.stream
                        .memcpy_htod(&self.scratch.$host[..$n], dst)
                        .map_err(|e| anyhow!("H2D({}) failed: {e:?}", stringify!($dev)))?;
                }
            }};
        }

        upload!(dev_pos, pos_host, num_variants);
        upload!(dev_start, start_host, num_variants);
        upload!(dev_numal, numal_host, num_variants);
        upload!(dev_alleles, alleles_host, batch.variant_alleles.len());
        upload!(dev_lens, lens_host, batch.allele_lengths.len());
        upload!(dev_offs, offs_host, batch.allele_byte_offsets.len());
        upload!(dev_geno, geno_host, batch.genotype_indices.len());
        upload!(dev_blk_off, blk_off_host, block_offsets.len());
        upload!(dev_blk_ref, blk_ref_host, block_ref_starts.len());

        timing.h2d_ns = t0.elapsed().as_nanos() as u64;

        let func = self.device.kernel()?;

        let block_dim = {
            let want = (num_haps as u32).max(1);
            let rounded = ((want + GPU_WARP_SIZE - 1) / GPU_WARP_SIZE) * GPU_WARP_SIZE;
            rounded.min(GPU_MAX_BLOCK_DIM)
        };
        let grid_x = ((num_haps as u32) + block_dim - 1) / block_dim;
        let grid_y = n_blocks as u32;

        let cfg = LaunchConfig {
            grid_dim: (grid_x.max(1), grid_y.max(1), 1),
            block_dim: (block_dim, 1, 1),
            shared_mem_bytes: 0,
        };

        let ref_start_i = batch.ref_start as i32;
        let ref_end_i = batch.ref_end as i32;
        let nv_i = num_variants as i32;
        let nh_i = num_haps as i32;
        let vbs_i = GPU_VAR_BLOCK as i32;

        let t1 = Instant::now();
        unsafe {
            self.stream
                .launch_builder(&func)
                .arg(&*dev_ref)
                .arg(&ref_start_i)
                .arg(&ref_end_i)
                .arg(self.scratch.dev_pos.ensure(&self.stream, num_variants.max(1))?)
                .arg(
                    self.scratch
                        .dev_alleles
                        .ensure(&self.stream, batch.variant_alleles.len().max(1))?,
                )
                .arg(
                    self.scratch
                        .dev_lens
                        .ensure(&self.stream, batch.allele_lengths.len().max(1))?,
                )
                .arg(
                    self.scratch
                        .dev_offs
                        .ensure(&self.stream, batch.allele_byte_offsets.len().max(1))?,
                )
                .arg(self.scratch.dev_start.ensure(&self.stream, num_variants.max(1))?)
                .arg(self.scratch.dev_numal.ensure(&self.stream, num_variants.max(1))?)
                .arg(
                    self.scratch
                        .dev_geno
                        .ensure(&self.stream, (num_haps * num_variants).max(1))?,
                )
                .arg(
                    self.scratch
                        .dev_blk_off
                        .ensure(&self.stream, block_offsets.len().max(1))?,
                )
                .arg(
                    self.scratch
                        .dev_blk_ref
                        .ensure(&self.stream, block_ref_starts.len().max(1))?,
                )
                .arg(&nv_i)
                .arg(&nh_i)
                .arg(&vbs_i)
                .arg(self.scratch.dev_out.ensure(&self.stream, total_out.max(1))?)
                .launch(cfg)
                .map_err(|e| anyhow!("kernel launch failed: {e:?}"))?;
        }
        timing.kernel_ns = t1.elapsed().as_nanos() as u64;

        let t2 = Instant::now();
        {
            let src_len = self
                .scratch
                .dev_out
                .ensure(&self.stream, total_out.max(1))?
                .len();
            self.scratch.out_host.ensure_capacity(src_len);

            let dev = self
                .scratch
                .dev_out
                .as_slice_ref()
                .expect("dev_out was just ensured");
            let dst = &mut self.scratch.out_host[..src_len];
            self.stream
                .memcpy_dtoh(dev, dst)
                .map_err(|e| anyhow!("D2H(out) failed: {e:?}"))?;
        }
        self.stream
            .synchronize()
            .map_err(|e| anyhow!("stream sync failed: {e:?}"))?;
        timing.d2h_ns = t2.elapsed().as_nanos() as u64;

        // Split flat output into per-hap vectors. `lengths[h]` is the byte
        // count for hap h; `block_offsets[0 * num_haps + h]` is its
        // starting position in the hap-major layout.
        let mut hap_seqs = Vec::with_capacity(num_haps);
        for h in 0..num_haps {
            let start = block_offsets[h] as usize;
            let len = lengths[h];
            hap_seqs.push(self.scratch.out_host[start..start + len].to_vec());
        }
        Ok((hap_seqs, timing))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kernel_source_is_nonempty() {
        assert!(KERNEL_SRC.contains("apply_variants"));
        assert!(KERNEL_SRC.contains("genotype_indices[v * num_haps + hap_id]"));
    }

    #[test]
    fn n_var_blocks_is_at_least_one() {
        let batch = GpuBatch {
            contig: "chr1".into(),
            ref_start: 0,
            ref_end: 100,
            variant_positions: vec![],
            variant_alleles: vec![],
            allele_lengths: vec![],
            allele_byte_offsets: vec![],
            allele_start_idx: vec![],
            num_alleles: vec![],
            genotype_indices: vec![],
            num_haps: 2,
        };
        assert_eq!(batch.n_var_blocks(), 1);
    }

    #[test]
    fn block_offsets_no_variants() {
        let batch = GpuBatch {
            contig: "chr1".into(),
            ref_start: 0,
            ref_end: 100,
            variant_positions: vec![],
            variant_alleles: vec![],
            allele_lengths: vec![],
            allele_byte_offsets: vec![],
            allele_start_idx: vec![],
            num_alleles: vec![],
            genotype_indices: vec![],
            num_haps: 2,
        };
        let (offsets, lengths) = batch.compute_block_offsets();
        // 2 rows (n_blocks + 1 = 2) × 2 haps.
        assert_eq!(offsets.len(), 4);
        // Block 0 starts: hap 0 at flat 0, hap 1 at flat 100.
        assert_eq!(offsets[0], 0);
        assert_eq!(offsets[1], 100);
        // End row: hap 0 ends at 100, hap 1 ends at 200.
        assert_eq!(offsets[2], 100);
        assert_eq!(offsets[3], 200);
        assert_eq!(lengths, vec![100, 100]);
    }

    #[test]
    fn block_offsets_deletion_flat_layout() {
        // One variant at pos 10, REF="AC" (2 bytes), ALT="A" (1 byte).
        // Hap 0 calls ALT (code 1) → length 19.
        // Hap 1 calls REF (code 0) → length 20.
        // Flat layout: hap 0 occupies [0, 19), hap 1 occupies [19, 39).
        let batch = GpuBatch {
            contig: "chr1".into(),
            ref_start: 0,
            ref_end: 20,
            variant_positions: vec![10],
            variant_alleles: vec![b'A', b'C', b'A', b'N'],
            allele_lengths: vec![2, 1, 1],
            allele_byte_offsets: vec![0, 2, 3],
            allele_start_idx: vec![0],
            num_alleles: vec![3],
            genotype_indices: vec![1, 0],
            num_haps: 2,
        };
        let (offsets, lengths) = batch.compute_block_offsets();
        assert_eq!(lengths, vec![19, 20]);
        assert_eq!(offsets[0], 0);
        assert_eq!(offsets[1], 19);
        assert_eq!(offsets[2], 19);
        assert_eq!(offsets[3], 39);
    }

    #[test]
    fn block_offsets_reference_call() {
        let batch = GpuBatch {
            contig: "chr1".into(),
            ref_start: 0,
            ref_end: 20,
            variant_positions: vec![10],
            variant_alleles: vec![b'A', b'C', b'N'],
            allele_lengths: vec![2, 1],
            allele_byte_offsets: vec![0, 2],
            allele_start_idx: vec![0],
            num_alleles: vec![2],
            genotype_indices: vec![-1],
            num_haps: 1,
        };
        let (offsets, lengths) = batch.compute_block_offsets();
        assert_eq!(lengths, vec![19]);
        assert_eq!(offsets.len(), 2); // (n_blocks + 1) × 1
        assert_eq!(offsets[0], 0);
        assert_eq!(offsets[1], 19);
    }

    #[test]
    fn block_offsets_multi_hap_are_disjoint() {
        // Two haps, five variants, all REF calls. Each hap gets 100 bytes
        // and their regions must not overlap in flat space.
        let batch = GpuBatch {
            contig: "chr1".into(),
            ref_start: 0,
            ref_end: 100,
            variant_positions: vec![10, 20, 30, 40, 50],
            variant_alleles: vec![
                b'A', b'T', b'A', b'T', b'A', b'T', b'A', b'T', b'A', b'T',
            ],
            allele_lengths: vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            allele_byte_offsets: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            allele_start_idx: vec![0, 2, 4, 6, 8],
            num_alleles: vec![2, 2, 2, 2, 2],
            genotype_indices: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            num_haps: 2,
        };
        let (offsets, lengths) = batch.compute_block_offsets();
        assert_eq!(lengths, vec![100, 100]);
        assert_eq!(offsets.len(), (batch.n_var_blocks() + 1) * 2);

        // Block 0 starts: hap 0 → 0, hap 1 → 100.
        assert_eq!(offsets[0], 0);
        assert_eq!(offsets[1], 100);

        // End row: hap 0 → 100, hap 1 → 200.
        let last_row = batch.n_var_blocks() * 2;
        assert_eq!(offsets[last_row], 100);
        assert_eq!(offsets[last_row + 1], 200);
    }

    #[test]
    fn block_offsets_multi_block_monotonic() {
        // Verify that within each hap's region, block offsets are
        // non-decreasing. This is the invariant that guarantees no two
        // blocks write to the same byte.
        let batch = GpuBatch {
            contig: "chr1".into(),
            ref_start: 0,
            ref_end: 4 * GPU_VAR_BLOCK as usize + 1000,
            variant_positions: (0..(4 * GPU_VAR_BLOCK as i32))
                .map(|i| 10 + i * 10)
                .collect(),
            variant_alleles: (0..(4 * GPU_VAR_BLOCK))
                .flat_map(|_| vec![b'A', b'T'])
                .collect(),
            allele_lengths: (0..(4 * GPU_VAR_BLOCK))
                .flat_map(|_| vec![1, 1])
                .collect(),
            allele_byte_offsets: (0..(4 * GPU_VAR_BLOCK))
                .flat_map(|i| vec![(2 * i) as i32, (2 * i + 1) as i32])
                .collect(),
            allele_start_idx: (0..(4 * GPU_VAR_BLOCK)).map(|i| (2 * i) as i32).collect(),
            num_alleles: vec![2; 4 * GPU_VAR_BLOCK],
            genotype_indices: vec![0; 4 * GPU_VAR_BLOCK * 2],
            num_haps: 2,
        };
        let n_blocks = batch.n_var_blocks();
        assert!(n_blocks >= 4, "expected multiple blocks, got {}", n_blocks);

        let (offsets, _lengths) = batch.compute_block_offsets();

        for h in 0..2 {
            for b in 0..n_blocks {
                let cur = offsets[b * 2 + h];
                let next = offsets[(b + 1) * 2 + h];
                assert!(
                    next >= cur,
                    "hap {} block {} offset {} > block {} offset {}",
                    h, b, cur, b + 1, next
                );
            }
        }

        // And hap 1's region starts exactly where hap 0's ends.
        let hap0_end = offsets[n_blocks * 2 + 0];
        let hap1_start = offsets[0 * 2 + 1];
        assert_eq!(
            hap0_end, hap1_start,
            "hap 1 must start exactly where hap 0 ends"
        );
    }
}