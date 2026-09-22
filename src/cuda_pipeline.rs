//! Bounded, multi-stream, multi-device CUDA pipeline.
//!
//! One producer thread builds [`GpuBatch`]es and pushes them into a bounded
//! channel. Each CUDA worker thread owns a [`GpuWorker`] (stream + reusable
//! scratch) and pulls batches from the shared channel. When multiple CUDA
//! devices are configured, worker streams are assigned round-robin across
//! them, and each device gets its own [`ReferenceCache`] because device
//! buffers live in that device's memory.

#![cfg(feature = "cuda")]

use crate::gpu::{BatchStats, BatchTiming, CudaDevice, GpuBatch, GpuWorker, ReferenceCache};
use anyhow::{anyhow, Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::thread;

// ---------------------------------------------------------------------------
// Stats
// ---------------------------------------------------------------------------

pub struct PipelineStats {
    h2d_ns: AtomicU64,
    kernel_ns: AtomicU64,
    d2h_ns: AtomicU64,
    batches: AtomicU64,
}

impl PipelineStats {
    pub fn new() -> Arc<Self> {
        Arc::new(Self {
            h2d_ns: AtomicU64::new(0),
            kernel_ns: AtomicU64::new(0),
            d2h_ns: AtomicU64::new(0),
            batches: AtomicU64::new(0),
        })
    }

    pub fn observe(&self, t: BatchTiming) {
        self.h2d_ns.fetch_add(t.h2d_ns, Ordering::Relaxed);
        self.kernel_ns.fetch_add(t.kernel_ns, Ordering::Relaxed);
        self.d2h_ns.fetch_add(t.d2h_ns, Ordering::Relaxed);
        self.batches.fetch_add(1, Ordering::Relaxed);
    }

    pub fn summary(&self) -> String {
        let b = self.batches.load(Ordering::Relaxed).max(1);
        format!(
            "GPU batches={} avg_h2d={}us avg_kernel={}us avg_d2h={}us",
            b,
            self.h2d_ns.load(Ordering::Relaxed) / b / 1000,
            self.kernel_ns.load(Ordering::Relaxed) / b / 1000,
            self.d2h_ns.load(Ordering::Relaxed) / b / 1000,
        )
    }
}

// ---------------------------------------------------------------------------
// Pipeline messages
// ---------------------------------------------------------------------------

/// A batch ready to be handed to a worker, together with everything the
/// worker needs to know about its origin — including the stats the batch
/// builder already computed.
pub struct PendingBatch {
    pub index: usize,
    pub batch: GpuBatch,
    pub contig_ref: Arc<Vec<u8>>,
    pub stats: BatchStats,
}

/// What a CUDA worker sends back: the per-hap byte vectors, the timing, and
/// the stats that came with the batch.
pub type PipelineMessage = (usize, Result<(Vec<Vec<u8>>, BatchTiming, BatchStats)>);

// ---------------------------------------------------------------------------
// Pipeline
// ---------------------------------------------------------------------------

pub struct CudaPipeline {
    /// Non-empty. Streams are assigned round-robin: `stream_idx % devices.len()`.
    devices: Vec<Arc<CudaDevice>>,
    num_streams: usize,
    in_flight: usize,
    stats: Arc<PipelineStats>,
}

impl CudaPipeline {
    pub fn new(
        devices: Vec<Arc<CudaDevice>>,
        num_streams: usize,
        in_flight: usize,
        stats: Arc<PipelineStats>,
    ) -> Self {
        Self {
            devices,
            num_streams: num_streams.max(1),
            in_flight: in_flight.max(1),
            stats,
        }
    }

    pub fn spawn(
        self,
    ) -> Result<(
        Sender<PendingBatch>,
        Receiver<PipelineMessage>,
        Vec<thread::JoinHandle<Result<()>>>,
    )> {
        if self.devices.is_empty() {
            return Err(anyhow!(
                "CudaPipeline::spawn called with zero devices — the scheduler \
                 must guarantee at least one for a GPU plan."
            ));
        }

        let (in_tx, in_rx) = bounded::<PendingBatch>(self.in_flight);
        let (out_tx, out_rx) = bounded::<PipelineMessage>(self.in_flight);

        // One reference cache per device: uploaded device buffers are not
        // portable across CUDA devices.
        let ref_caches: Vec<Arc<ReferenceCache>> = (0..self.devices.len())
            .map(|_| Arc::new(ReferenceCache::new()))
            .collect();

        let n_dev = self.devices.len();
        let mut handles = Vec::with_capacity(self.num_streams);

        for stream_idx in 0..self.num_streams {
            let dev_idx = stream_idx % n_dev;
            let device = self.devices[dev_idx].clone();
            let in_rx = in_rx.clone();
            let out_tx = out_tx.clone();
            let stats = self.stats.clone();
            let ref_cache = ref_caches[dev_idx].clone();

            let name = format!("cuda-worker-{}-{}", stream_idx, device.name);
            let h = thread::Builder::new()
                .name(name)
                .spawn(move || -> Result<()> {
                    let mut worker = GpuWorker::new(device, ref_cache)
                        .with_context(|| "worker: could not create GpuWorker")?;
                    while let Ok(pending) = in_rx.recv() {
                        let res = worker
                            .run(&pending.batch, &pending.contig_ref)
                            .map(|(seqs, t)| {
                                stats.observe(t);
                                (seqs, t, pending.stats)
                            });
                        if out_tx.send((pending.index, res)).is_err() {
                            break;
                        }
                    }
                    Ok(())
                })
                .map_err(|e| anyhow!("could not spawn CUDA worker: {e}"))?;
            handles.push(h);
        }

        drop(out_tx);
        Ok((in_tx, out_rx, handles))
    }
}