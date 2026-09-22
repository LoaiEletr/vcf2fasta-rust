//! Central resource budget.
//!
//! Splits machine resources between the phasing producer and the vcf2fasta
//! consumers so that one stage cannot starve the other.
//!
//! ## Memory ceiling
//!
//! User-supplied `--max-memory` is clamped to a safe fraction of the
//! machine's RAM regardless of what the user asked for. Asking for more
//! than the machine has is a request the tool cannot safely honor, so it
//! refuses to try. The original value is preserved in
//! [`ResourceBudget::max_memory_clamped_from`] so the caller can log the
//! reduction; the effective value is always `max_memory_bytes`.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Duration;

/// Shared, cloneable budget handed to producer and consumers.
#[derive(Clone)]
pub struct ResourceBudget {
    pub total_cores: usize,
    pub phasing_threads: usize,
    pub vcf2fasta_threads: usize,
    pub total_ram_bytes: u64,
    pub available_ram_bytes: u64,
    /// Effective memory ceiling the pipeline will honor.
    pub max_memory_bytes: u64,
    pub max_vram_bytes: u64,
    /// Set when the user passed `--max-memory` larger than the safe cap.
    /// Holds the original value so the caller can log the clamp.
    pub max_memory_clamped_from: Option<u64>,
    peak_rss_bytes: Arc<AtomicU64>,
}

impl ResourceBudget {
    pub fn detect(max_memory_bytes: Option<u64>, max_vram_bytes: Option<u64>) -> Self {
        let total_cores = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        let (total_ram, avail_ram) = crate::scheduler::detect_system_ram();

        // Reserve one core for orchestration; split the rest evenly.
        let worker_cores = total_cores.saturating_sub(1).max(1);
        let phasing_threads = (worker_cores / 2).max(1);
        let vcf2fasta_threads = (worker_cores - phasing_threads).max(1);

        // Hard cap: 90% of `min(total, available)`. Using the min
        // protects the case where another process is holding RAM: the
        // user's `--max-memory` cannot push us past what is actually
        // free right now.
        let hard_cap = ((total_ram.min(avail_ram) as f64) * 0.90) as u64;

        let (effective_memory, clamped_from) = match max_memory_bytes {
            Some(user) if user > hard_cap => (hard_cap, Some(user)),
            Some(user) => (user, None),
            None => (((avail_ram as f64) * 0.80) as u64, None),
        };

        Self {
            total_cores,
            phasing_threads,
            vcf2fasta_threads,
            total_ram_bytes: total_ram,
            available_ram_bytes: avail_ram,
            max_memory_bytes: effective_memory,
            max_vram_bytes: max_vram_bytes.unwrap_or(0),
            max_memory_clamped_from: clamped_from,
            peak_rss_bytes: Arc::new(AtomicU64::new(0)),
        }
    }

    /// Sample current RSS and update the high-water mark. Cheap on Linux
    /// (`/proc/self/statm`); no-op elsewhere.
    pub fn sample_rss(&self) -> u64 {
        let now = read_rss_bytes().unwrap_or(0);
        if now == 0 {
            return self.peak_rss_bytes.load(Ordering::Relaxed);
        }
        let mut prev = self.peak_rss_bytes.load(Ordering::Relaxed);
        while now > prev {
            match self.peak_rss_bytes.compare_exchange_weak(
                prev,
                now,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => prev = actual,
            }
        }
        prev
    }

    pub fn peak_rss_bytes(&self) -> u64 {
        self.peak_rss_bytes.load(Ordering::Relaxed)
    }

    /// True if a new allocation of `extra_bytes` would exceed the budget.
    pub fn would_exceed_memory(&self, extra_bytes: u64) -> bool {
        let current = read_rss_bytes().unwrap_or(0);
        current.saturating_add(extra_bytes) > self.max_memory_bytes
    }
}

#[cfg(target_os = "linux")]
fn read_rss_bytes() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/statm").ok()?;
    let rss_pages: u64 = s.split_whitespace().nth(1)?.parse().ok()?;
    Some(rss_pages.saturating_mul(4096))
}

#[cfg(not(target_os = "linux"))]
fn read_rss_bytes() -> Option<u64> {
    None
}

/// Per-run timing measurements, logged at the end.
#[derive(Default, Clone)]
pub struct RunMetrics {
    pub phasing_total: Duration,
    pub vcf2fasta_total: Duration,
    pub queue_wait_total: Duration,
    pub contigs_phased: usize,
    pub contigs_already_phased: usize,
}

impl RunMetrics {
    pub fn new() -> Arc<std::sync::Mutex<Self>> {
        Arc::new(std::sync::Mutex::new(Self::default()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_splits_cores() {
        let b = ResourceBudget::detect(None, None);
        assert!(b.phasing_threads >= 1);
        assert!(b.vcf2fasta_threads >= 1);
        assert!(b.phasing_threads + b.vcf2fasta_threads <= b.total_cores);
    }

    #[test]
    fn budget_respects_small_max_memory_override() {
        // Any small value will be under the hard cap on any real machine,
        // so it should be honored exactly.
        let b = ResourceBudget::detect(Some(123_456_789), None);
        assert_eq!(b.max_memory_bytes, 123_456_789);
        assert!(b.max_memory_clamped_from.is_none());
    }

    #[test]
    fn budget_clamps_oversized_max_memory() {
        // Ask for 1 PiB. Must be clamped to the machine's hard cap.
        let huge = 1u64 << 50; // 1 PiB
        let b = ResourceBudget::detect(Some(huge), None);
        assert!(b.max_memory_bytes < huge);
        assert_eq!(b.max_memory_clamped_from, Some(huge));
        // Sanity: effective value should be no more than total machine RAM.
        assert!(b.max_memory_bytes <= b.total_ram_bytes);
    }

    #[test]
    fn sample_rss_is_monotonic() {
        let b = ResourceBudget::detect(None, None);
        let _ = b.sample_rss();
        let p1 = b.peak_rss_bytes();
        let _ = b.sample_rss();
        assert!(b.peak_rss_bytes() >= p1);
    }
}