//! Adaptive workload scheduler.

use anyhow::{anyhow, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Device {
    Auto,
    Cpu,
    Gpu,
}

impl Default for Device {
    fn default() -> Self {
        Device::Auto
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeviceChoice {
    Cpu,
    Gpu,
}

impl DeviceChoice {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Cpu => "cpu",
            Self::Gpu => "gpu",
        }
    }
}

#[derive(Debug, Clone)]
pub struct GpuInfo {
    pub name: String,
    pub total_vram_bytes: u64,
    pub free_vram_bytes: u64,
    pub max_alloc_bytes: u64,
    /// 0-based CUDA device index.
    pub device_index: usize,
}

#[derive(Debug, Clone)]
pub struct HardwareInfo {
    pub cpu_cores: usize,
    pub total_ram_bytes: u64,
    pub available_ram_bytes: u64,
    pub gpus: Vec<GpuInfo>,
}

impl HardwareInfo {
    pub fn detect() -> Self {
        let cpu_cores = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        let (total_ram, avail_ram) = detect_system_ram();
        Self {
            cpu_cores,
            total_ram_bytes: total_ram,
            available_ram_bytes: avail_ram,
            gpus: detect_all_gpus(),
        }
    }

    /// Back-compat accessor: first GPU if any.
    pub fn gpu(&self) -> Option<&GpuInfo> {
        self.gpus.first()
    }
}

#[derive(Debug, Clone, Default)]
pub struct WorkloadProfile {
    pub sample_count: usize,
    pub haplotype_count: usize,
    pub max_ploidy: usize,
    pub sample_max_ploidies: Vec<usize>,
    pub variant_count: u64,
    pub reference_bases: u64,
    pub phased_genotypes: u64,
    pub unphased_genotypes: u64,
    pub malformed_records: u64,
    pub snv_count: u64,
    pub indel_count: u64,
    pub contig_count: usize,
    pub estimated_output_bytes: u64,
}

impl WorkloadProfile {
    pub fn work_units(&self) -> u64 {
        let var_work = self.variant_count.saturating_mul(self.haplotype_count as u64);
        let ref_work = self.reference_bases.saturating_mul(self.haplotype_count as u64);
        var_work.saturating_add(ref_work)
    }

    pub fn estimate_output_bytes(
        haplotype_count: usize,
        reference_bases: u64,
        line_width: usize,
    ) -> u64 {
        if haplotype_count == 0 || reference_bases == 0 {
            return 0;
        }
        let wrap = if line_width == 0 { 80 } else { line_width } as u64;
        let lines = reference_bases.div_ceil(wrap);
        let per_hap = reference_bases.saturating_add(lines).saturating_add(64);
        per_hap.saturating_mul(haplotype_count as u64)
    }
}

#[derive(Debug, Clone)]
pub struct ExecutionPlan {
    pub device: DeviceChoice,
    pub cpu_workers: usize,
    pub haplotype_block_size: usize,
    pub variant_block_size: u64,
    pub gpu_streams: usize,
    pub in_flight_buffers: usize,
    /// CUDA device indices the execution layer must use, in order. The
    /// GPU pipeline distributes worker streams round-robin across these.
    /// Never empty when `device == Gpu`.
    pub gpu_device_indices: Vec<usize>,
    pub gpu_memory_budget_bytes: u64,
    pub cpu_memory_budget_bytes: u64,
    pub reason: String,
}

pub mod limits {
    pub const CPU_RAM_FRACTION: f64 = 0.50;
    pub const GPU_VRAM_FRACTION: f64 = 0.50;
    pub const MIN_FREE_VRAM_BYTES: u64 = 256 * 1024 * 1024;
    pub const GPU_MIN_WORK_UNITS: u64 = 20_000_000;
    pub const GPU_WORK_PER_STREAM: u64 = 100_000_000;
    pub const MAX_GPU_STREAMS: usize = 4;
    pub const MAX_IN_FLIGHT: usize = 4;
    pub const MIN_HAP_BLOCK: usize = 1;
    pub const MIN_VARIANT_BLOCK: u64 = 10_000;
    pub const MAX_VARIANT_BLOCK: u64 = 1_000_000;
    pub const BYTES_PER_VAR_HAP: u64 = 16;
    pub const BYTES_PER_REF_BASE: u64 = 1;
    pub const MAX_TILE_OUTPUT_BYTES: u64 = 16 * 1024 * 1024;
}

pub struct Scheduler;

impl Scheduler {
    /// Single-device plan. Uses the first visible GPU when a GPU plan is
    /// chosen; explicit multi-device plans go through
    /// [`Scheduler::plan_multi_gpu`].
    ///
    /// `max_vram_bytes` is the user's `--max-vram` budget, or `0` for
    /// "unset". When nonzero it is applied as a *cap*: the scheduler will
    /// not budget more VRAM than the user asked for, but it will still
    /// respect the safe-fraction heuristic (a user cannot raise the budget
    /// above what the auto-derived value would have been).
    pub fn plan(
        requested: Device,
        hw: &HardwareInfo,
        wl: &WorkloadProfile,
        max_vram_bytes: u64,
    ) -> Result<ExecutionPlan> {
        let cpu_workers = choose_cpu_workers(hw, wl);
        let cpu_budget = ((hw.available_ram_bytes as f64) * limits::CPU_RAM_FRACTION) as u64;

        let device = match requested {
            Device::Cpu => DeviceChoice::Cpu,
            Device::Gpu => {
                if hw.gpus.is_empty() {
                    return Err(anyhow!(
                        "--device gpu was requested but no compatible GPU is available. \
                         Causes: the binary was built without `--features cuda`, no NVIDIA \
                         driver is loaded, or no CUDA device is visible to the process."
                    ));
                }
                DeviceChoice::Gpu
            }
            Device::Auto => {
                if !hw.gpus.is_empty() && wl.work_units() >= limits::GPU_MIN_WORK_UNITS {
                    DeviceChoice::Gpu
                } else {
                    DeviceChoice::Cpu
                }
            }
        };

        match device {
            DeviceChoice::Cpu => Ok(cpu_plan(cpu_workers, cpu_budget, wl)),
            DeviceChoice::Gpu => {
                let first = hw.gpus.first().map(|g| g.device_index).unwrap_or(0);
                gpu_plan(hw, cpu_workers, cpu_budget, wl, &[first], max_vram_bytes)
            }
        }
    }

    /// Multi-GPU variant.
    ///
    /// `gpu_indices` names the CUDA device indices the plan must use. An
    /// empty list falls back to a CPU plan. Duplicate indices are
    /// collapsed preserving first-seen order. An index that is not visible
    /// to the process is a hard error, not a silent fallback to device 0.
    pub fn plan_multi_gpu(
        hw: &HardwareInfo,
        wl: &WorkloadProfile,
        gpu_indices: Vec<usize>,
        max_vram_bytes: u64,
    ) -> Result<ExecutionPlan> {
        if gpu_indices.is_empty() {
            return Scheduler::plan(Device::Cpu, hw, wl, max_vram_bytes);
        }

        // Collapse duplicates preserving first-seen order.
        let mut seen = std::collections::HashSet::with_capacity(gpu_indices.len());
        let unique: Vec<usize> = gpu_indices
            .into_iter()
            .filter(|i| seen.insert(*i))
            .collect();

        let cpu_workers = choose_cpu_workers(hw, wl);
        let cpu_budget = ((hw.available_ram_bytes as f64) * limits::CPU_RAM_FRACTION) as u64;
        gpu_plan(hw, cpu_workers, cpu_budget, wl, &unique, max_vram_bytes)
    }
}

fn choose_cpu_workers(hw: &HardwareInfo, wl: &WorkloadProfile) -> usize {
    let cores = hw.cpu_cores.max(1);
    let variant_cap = (wl.variant_count as usize).max(1);
    cores.min(variant_cap).max(1)
}

fn cpu_plan(cpu_workers: usize, cpu_budget: u64, wl: &WorkloadProfile) -> ExecutionPlan {
    ExecutionPlan {
        device: DeviceChoice::Cpu,
        cpu_workers,
        haplotype_block_size: wl.haplotype_count.max(1),
        variant_block_size: wl.variant_count.max(1),
        gpu_streams: 0,
        in_flight_buffers: 1,
        gpu_device_indices: Vec::new(),
        gpu_memory_budget_bytes: 0,
        cpu_memory_budget_bytes: cpu_budget,
        reason: format!("cpu path, {} workers", cpu_workers),
    }
}

/// Build a GPU plan over an explicit, non-empty list of CUDA device indices.
fn gpu_plan(
    hw: &HardwareInfo,
    cpu_workers: usize,
    cpu_budget: u64,
    wl: &WorkloadProfile,
    gpu_indices: &[usize],
    max_vram_bytes: u64,
) -> Result<ExecutionPlan> {
    debug_assert!(!gpu_indices.is_empty());

    // ---- Resolve requested indices to actual GPUs ----------------------
    // This is where `--gpu-devices` becomes binding: a request for a
    // device the process cannot see is a hard error, not a silent
    // fallback to device 0.
    let mut selected: Vec<&GpuInfo> = Vec::with_capacity(gpu_indices.len());
    for &idx in gpu_indices {
        match hw.gpus.iter().find(|g| g.device_index == idx) {
            Some(g) => selected.push(g),
            None => {
                let visible = if hw.gpus.is_empty() {
                    "none".to_string()
                } else {
                    hw.gpus
                        .iter()
                        .map(|g| format!("{}={}", g.device_index, g.name))
                        .collect::<Vec<_>>()
                        .join(", ")
                };
                return Err(anyhow!(
                    "--gpu-devices requested CUDA device {}, but it is not visible to this \
                     process. {} device(s) detected: {}.",
                    idx,
                    hw.gpus.len(),
                    visible,
                ));
            }
        }
    }
    if selected.is_empty() {
        return Err(anyhow!("gpu plan requested with no usable devices"));
    }

    // ---- Per-device free-VRAM floor -------------------------------------
    for g in &selected {
        if g.free_vram_bytes < limits::MIN_FREE_VRAM_BYTES {
            return Err(anyhow!(
                "GPU {} ({}) has only {} MiB free VRAM, need at least {} MiB",
                g.device_index,
                g.name,
                g.free_vram_bytes / (1024 * 1024),
                limits::MIN_FREE_VRAM_BYTES / (1024 * 1024),
            ));
        }
    }

    // ---- Auto-derived budget --------------------------------------------
    let auto_budget: u64 = selected
        .iter()
        .map(|g| {
            let frac = ((g.free_vram_bytes as f64) * limits::GPU_VRAM_FRACTION) as u64;
            frac.min(g.max_alloc_bytes.max(1))
        })
        .sum();

    // ---- User cap (`--max-vram`) ----------------------------------------
    // 0 means "unset". The user can only *lower* the auto budget; a value
    // above the auto budget is silently ignored (the caller logs the note).
    let gpu_budget = if max_vram_bytes > 0 {
        auto_budget.min(max_vram_bytes).max(1)
    } else {
        auto_budget.max(1)
    };

    // ---- Stream and in-flight counts ------------------------------------
    let work = wl.work_units().max(1);
    let streams_raw = work / limits::GPU_WORK_PER_STREAM;
    let max_streams = limits::MAX_GPU_STREAMS.saturating_mul(selected.len());
    let gpu_streams = (streams_raw as usize).clamp(selected.len(), max_streams.max(1));
    let max_in_flight = limits::MAX_IN_FLIGHT.saturating_mul(selected.len());
    let in_flight = (gpu_streams + 1).min(max_in_flight);

    // ---- Tile sizing ----------------------------------------------------
    let haps = wl.haplotype_count.max(1) as u64;
    let base_block_from_output = (limits::MAX_TILE_OUTPUT_BYTES / haps).max(20_000);

    let per_device_budget = (gpu_budget / selected.len() as u64).max(1);
    let per_device_in_flight = (in_flight as u64 / selected.len() as u64).max(1);
    let per_buffer = (per_device_budget / per_device_in_flight).max(1);
    let per_hap_per_variant = limits::BYTES_PER_VAR_HAP;
    let max_vars_per_tile = (per_buffer / (haps * per_hap_per_variant).max(1)).max(1);

    let density = if wl.reference_bases > 0 {
        wl.variant_count as f64 / wl.reference_bases as f64
    } else {
        1.0
    };
    let base_block_from_input = if density > 0.0 {
        (max_vars_per_tile as f64 / density) as u64
    } else {
        base_block_from_output
    };

    let variant_block_size = base_block_from_output
        .min(base_block_from_input.max(20_000))
        .clamp(20_000, limits::MAX_VARIANT_BLOCK);

    let haplotype_block_size = wl.haplotype_count.max(1);

    let indices_str = selected
        .iter()
        .map(|g| g.device_index.to_string())
        .collect::<Vec<_>>()
        .join(",");

    Ok(ExecutionPlan {
        device: DeviceChoice::Gpu,
        cpu_workers,
        haplotype_block_size,
        variant_block_size,
        gpu_streams,
        in_flight_buffers: in_flight,
        gpu_device_indices: selected.iter().map(|g| g.device_index).collect(),
        gpu_memory_budget_bytes: gpu_budget,
        cpu_memory_budget_bytes: cpu_budget,
        reason: format!(
            "gpu: {} device(s) [{}], {} stream(s), hap_block={}, base_block={}, budget={} MiB",
            selected.len(),
            indices_str,
            gpu_streams,
            haplotype_block_size,
            variant_block_size,
            gpu_budget / (1024 * 1024),
        ),
    })
}

pub fn detect_system_ram() -> (u64, u64) {
    if let Ok(s) = std::fs::read_to_string("/proc/meminfo") {
        let mut total_kb = 0u64;
        let mut avail_kb = 0u64;
        for line in s.lines() {
            if let Some(rest) = line.strip_prefix("MemTotal:") {
                total_kb = rest
                    .split_whitespace()
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
            } else if let Some(rest) = line.strip_prefix("MemAvailable:") {
                avail_kb = rest
                    .split_whitespace()
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
            }
        }
        if total_kb > 0 {
            return (total_kb.saturating_mul(1024), avail_kb.saturating_mul(1024));
        }
    }
    (2 * 1024 * 1024 * 1024, 1024 * 1024 * 1024)
}

pub fn detect_available_ram_bytes() -> u64 {
    detect_system_ram().1
}

/// Enumerate all CUDA devices.
///
/// When the binary was built **without** the `cuda` feature, returns an
/// empty vector silently — that is expected. When the feature **is**
/// enabled but enumeration fails (driver missing, no device, symbol
/// resolution error), the failure is reported on stderr so that a later
/// `--device gpu` refusal has a visible cause.
fn detect_all_gpus() -> Vec<GpuInfo> {
    #[cfg(feature = "cuda")]
    {
        match crate::gpu::enumerate_devices() {
            Ok(v) => v,
            Err(e) => {
                eprintln!(
                    "[SCHEDULER] CUDA device enumeration failed: {:#}. \
                     GPU execution will be unavailable.",
                    e
                );
                Vec::new()
            }
        }
    }
    #[cfg(not(feature = "cuda"))]
    {
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hw_with_gpu() -> HardwareInfo {
        HardwareInfo {
            cpu_cores: 8,
            total_ram_bytes: 32 * 1024 * 1024 * 1024,
            available_ram_bytes: 24 * 1024 * 1024 * 1024,
            gpus: vec![GpuInfo {
                name: "test".into(),
                total_vram_bytes: 8 * 1024 * 1024 * 1024,
                free_vram_bytes: 6 * 1024 * 1024 * 1024,
                max_alloc_bytes: 2 * 1024 * 1024 * 1024,
                device_index: 0,
            }],
        }
    }

    fn hw_two_gpus() -> HardwareInfo {
        let mut hw = hw_with_gpu();
        hw.gpus.push(GpuInfo {
            name: "test2".into(),
            total_vram_bytes: 8 * 1024 * 1024 * 1024,
            free_vram_bytes: 6 * 1024 * 1024 * 1024,
            max_alloc_bytes: 2 * 1024 * 1024 * 1024,
            device_index: 1,
        });
        hw
    }

    fn wl_like_chr22() -> WorkloadProfile {
        WorkloadProfile {
            sample_count: 50,
            haplotype_count: 100,
            max_ploidy: 2,
            sample_max_ploidies: vec![2; 50],
            variant_count: 1_066_557,
            reference_bases: 50_818_468,
            ..Default::default()
        }
    }

    #[test]
    fn gpu_plan_caps_tile_output() {
        let plan = Scheduler::plan(Device::Gpu, &hw_with_gpu(), &wl_like_chr22(), 0).unwrap();
        let haps = wl_like_chr22().haplotype_count as u64;
        let per_tile_output = plan.variant_block_size * haps;
        assert!(per_tile_output <= limits::MAX_TILE_OUTPUT_BYTES * 2);
        assert!(plan.variant_block_size <= limits::MAX_VARIANT_BLOCK);
    }

    #[test]
    fn auto_picks_gpu_for_large_workload() {
        let plan = Scheduler::plan(Device::Auto, &hw_with_gpu(), &wl_like_chr22(), 0).unwrap();
        assert_eq!(plan.device, DeviceChoice::Gpu);
    }

    #[test]
    fn auto_falls_back_to_cpu_without_gpus() {
        let hw = HardwareInfo {
            cpu_cores: 8,
            total_ram_bytes: 32 * 1024 * 1024 * 1024,
            available_ram_bytes: 24 * 1024 * 1024 * 1024,
            gpus: Vec::new(),
        };
        let plan = Scheduler::plan(Device::Auto, &hw, &wl_like_chr22(), 0).unwrap();
        assert_eq!(plan.device, DeviceChoice::Cpu);
    }

    #[test]
    fn explicit_gpu_without_gpus_errors() {
        let hw = HardwareInfo {
            cpu_cores: 8,
            total_ram_bytes: 32 * 1024 * 1024 * 1024,
            available_ram_bytes: 24 * 1024 * 1024 * 1024,
            gpus: Vec::new(),
        };
        let err = Scheduler::plan(Device::Gpu, &hw, &wl_like_chr22(), 0).unwrap_err();
        let msg = format!("{:#}", err);
        assert!(msg.contains("no compatible GPU"), "got: {msg}");
    }

    #[test]
    fn multi_gpu_plan_sets_device_indices() {
        let hw = hw_two_gpus();
        let plan = Scheduler::plan_multi_gpu(&hw, &wl_like_chr22(), vec![0, 1], 0).unwrap();
        assert_eq!(plan.gpu_device_indices, vec![0, 1]);
        assert!(plan.gpu_streams >= 2);
    }

    #[test]
    fn multi_gpu_plan_dedupes_indices() {
        let hw = hw_two_gpus();
        let plan = Scheduler::plan_multi_gpu(&hw, &wl_like_chr22(), vec![1, 1, 0, 1], 0).unwrap();
        assert_eq!(plan.gpu_device_indices, vec![1, 0]);
    }

    #[test]
    fn multi_gpu_plan_rejects_unknown_device() {
        let hw = hw_with_gpu();
        let err = Scheduler::plan_multi_gpu(&hw, &wl_like_chr22(), vec![7], 0).unwrap_err();
        let msg = format!("{:#}", err);
        assert!(msg.contains("device 7"), "got: {msg}");
    }

    #[test]
    fn multi_gpu_empty_indices_falls_back_to_cpu() {
        let hw = hw_with_gpu();
        let plan = Scheduler::plan_multi_gpu(&hw, &wl_like_chr22(), Vec::new(), 0).unwrap();
        assert_eq!(plan.device, DeviceChoice::Cpu);
    }

    #[test]
    fn max_vram_caps_gpu_budget() {
        let hw = hw_with_gpu();
        let unlimited = Scheduler::plan(Device::Gpu, &hw, &wl_like_chr22(), 0).unwrap();
        let capped =
            Scheduler::plan(Device::Gpu, &hw, &wl_like_chr22(), 100 * 1024 * 1024).unwrap();
        assert!(capped.gpu_memory_budget_bytes <= 100 * 1024 * 1024);
        assert!(capped.gpu_memory_budget_bytes <= unlimited.gpu_memory_budget_bytes);
    }

    #[test]
    fn max_vram_above_auto_is_ignored() {
        let hw = hw_with_gpu();
        let unlimited = Scheduler::plan(Device::Gpu, &hw, &wl_like_chr22(), 0).unwrap();
        let capped = Scheduler::plan(Device::Gpu, &hw, &wl_like_chr22(), 1u64 << 40).unwrap();
        assert_eq!(
            capped.gpu_memory_budget_bytes,
            unlimited.gpu_memory_budget_bytes
        );
    }
}