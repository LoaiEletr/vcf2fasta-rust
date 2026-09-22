//! Structured logging.
//!
//! Two log streams:
//!
//! * The main log (`<prefix>.log`) receives stage-tagged lines: `[DISCOVERY]`,
//!   `[PHASING]`, `[READY]`, `[SCHEDULER]`, `[COMPLETE]`, `[SUMMARY]`, and
//!   one `[WARNINGS]` summary per contig.
//!
//! * The warnings log (`<prefix>.warnings.log`) receives every individual
//!   warning message verbatim. It is separate so that a VCF with millions of
//!   skipped variants does not drown the stage log, while still preserving
//!   the full detail for debugging and for the PI's validation reports.
//!
//! Both files are always created. stderr receives a capped live feed of
//! warnings (100 lines) so that an interactive run is not flooded.

use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

/// Hard ceiling on lines written to the warnings log. Guards against a
/// pathological VCF filling the disk. 1 M lines is roughly 100 MB.
const MAX_WARNING_LINES_IN_FILE: usize = 1_000_000;

/// Hard ceiling on warnings echoed live to stderr. After this, stderr gets
/// one note and then goes quiet; the file keeps receiving lines.
const MAX_WARNING_LINES_TO_STDERR: usize = 100;

pub struct Logger {
    /// Main stage log. Flushed per line: low volume, high value.
    file: File,

    /// Warnings log. Buffered: potentially millions of lines.
    warnings_file: Option<BufWriter<File>>,
    warnings_path: Option<PathBuf>,

    quiet: bool,

    /// Lines written to the warnings file so far.
    warnings_written: usize,

    /// Lines echoed to stderr so far.
    warnings_printed_to_stderr: usize,
}

impl Logger {
    /// Create both logs.
    ///
    /// * `path` — main log.
    /// * `warnings_path` — warnings log, or `None` to disable it.
    /// * `quiet` — suppress stderr entirely (files are always written).
    pub fn create(
        path: &Path,
        warnings_path: Option<&Path>,
        quiet: bool,
    ) -> Result<Self> {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() && !parent.exists() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let file = File::create(path)
            .with_context(|| format!("could not create log {}", path.display()))?;

        let (warnings_file, warnings_path_buf) = if let Some(wp) = warnings_path {
            if let Some(parent) = wp.parent() {
                if !parent.as_os_str().is_empty() && !parent.exists() {
                    std::fs::create_dir_all(parent)?;
                }
            }
            let f = File::create(wp).with_context(|| {
                format!("could not create warnings log {}", wp.display())
            })?;
            (Some(BufWriter::new(f)), Some(wp.to_path_buf()))
        } else {
            (None, None)
        };

        Ok(Self {
            file,
            warnings_file,
            warnings_path: warnings_path_buf,
            quiet,
            warnings_written: 0,
            warnings_printed_to_stderr: 0,
        })
    }

    pub fn quiet(&self) -> bool {
        self.quiet
    }

    pub fn has_warnings_file(&self) -> bool {
        self.warnings_file.is_some()
    }

    pub fn warnings_path(&self) -> Option<&Path> {
        self.warnings_path.as_deref()
    }

    /// Write a raw line to the main log (and stderr unless quiet).
    pub fn raw(&mut self, msg: &str) -> Result<()> {
        if !self.quiet {
            eprintln!("{}", msg);
        }
        writeln!(self.file, "{}", msg)?;
        self.file.flush()?;
        Ok(())
    }

    /// Write a tagged line to the main log.
    pub fn tag(&mut self, tag: &str, msg: &str) -> Result<()> {
        self.raw(&format!("[{}] {}", tag, msg))
    }

    // -----------------------------------------------------------------------
    // Warnings
    // -----------------------------------------------------------------------

    /// Record one individual warning.
    ///
    /// * Writes to the warnings log (if any), up to the file cap.
    /// * Echoes to stderr (unless quiet), up to the stderr cap.
    /// * Does NOT write to the main log — the main log receives only the
    ///   per-contig summary via [`warning_summary`].
    pub fn warn(&mut self, message: &str) -> Result<()> {
        // stderr, capped.
        if !self.quiet && self.warnings_printed_to_stderr < MAX_WARNING_LINES_TO_STDERR {
            eprintln!("WARNING {}", message);
            self.warnings_printed_to_stderr += 1;
            if self.warnings_printed_to_stderr == MAX_WARNING_LINES_TO_STDERR {
                if let Some(p) = &self.warnings_path {
                    eprintln!(
                        "# further warnings suppressed on stderr; full log at {}",
                        p.display()
                    );
                } else {
                    eprintln!("# further warnings suppressed on stderr");
                }
            }
        }

        // File, capped at the file ceiling.
        if let Some(f) = self.warnings_file.as_mut() {
            if self.warnings_written < MAX_WARNING_LINES_IN_FILE {
                writeln!(f, "WARNING {}", message)?;
                self.warnings_written += 1;
            } else if self.warnings_written == MAX_WARNING_LINES_IN_FILE {
                writeln!(
                    f,
                    "# warnings log reached {} lines; further warnings suppressed",
                    MAX_WARNING_LINES_IN_FILE
                )?;
                self.warnings_written += 1;
            }
        }
        Ok(())
    }

    /// Record the per-contig warning summary. Writes one line to the main
    /// log and a matching line to the warnings log (so the warnings file is
    /// self-contained).
    pub fn warning_summary(
        &mut self,
        contig: &str,
        total: usize,
        by_reason: &BTreeMap<&'static str, usize>,
    ) -> Result<()> {
        if total == 0 {
            return Ok(());
        }
        let parts: Vec<String> = by_reason
            .iter()
            .map(|(r, n)| format!("{}={}", r, n))
            .collect();
        let line = format!(
            "{}: total={} {{{}}}",
            contig,
            total,
            parts.join(", ")
        );
        self.tag("WARNINGS", &line)?;
        if let Some(f) = self.warnings_file.as_mut() {
            writeln!(f, "[WARNINGS] {}", line)?;
        }
        Ok(())
    }

    /// Flush the warnings file. Called once at end of run and after each
    /// per-contig summary.
    pub fn flush_warnings(&mut self) -> Result<()> {
        if let Some(f) = self.warnings_file.as_mut() {
            f.flush()?;
        }
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Stage helpers
    // -----------------------------------------------------------------------

    pub fn discovery(
        &mut self,
        contig: &str,
        samples: usize,
        haplotypes: usize,
        variants: u64,
        phased: u64,
        unphased: u64,
        backend: &str,
    ) -> Result<()> {
        self.tag(
            "DISCOVERY",
            &format!(
                "{}: {} samples, {} haplotypes, {} variants, phased={} unphased={} backend={}",
                contig, samples, haplotypes, variants, phased, unphased, backend
            ),
        )
    }

    pub fn phasing_start(&mut self, contig: &str, backend: &str, threads: usize) -> Result<()> {
        self.tag(
            "PHASING",
            &format!("{} -> {} (threads={})", contig, backend, threads),
        )
    }

    pub fn phasing_done(&mut self, contig: &str, elapsed: Duration) -> Result<()> {
        self.tag(
            "PHASING",
            &format!("{} -> done in {}", contig, fmt_dur(elapsed)),
        )
    }

    pub fn ready(&mut self, contig: &str) -> Result<()> {
        self.tag("READY", &format!("{} -> ready", contig))
    }

    pub fn scheduler_dispatch(&mut self, contig: &str, device: &str, tiles: usize) -> Result<()> {
        self.tag(
            "SCHEDULER",
            &format!("{} -> {} ({} tiles)", contig, device, tiles),
        )
    }

    pub fn gpu_device(
        &mut self,
        idx: usize,
        name: &str,
        free_mib: u64,
        total_mib: u64,
    ) -> Result<()> {
        self.tag(
            "GPU",
            &format!(
                "device {}: {} ({} MiB free / {} MiB total)",
                idx, name, free_mib, total_mib
            ),
        )
    }

    pub fn gpu_batch(&mut self, msg: &str) -> Result<()> {
        self.tag("GPU", msg)
    }

    pub fn complete(
        &mut self,
        contig: &str,
        seen: usize,
        applied: usize,
        files: usize,
        warnings: usize,
    ) -> Result<()> {
        self.tag(
            "COMPLETE",
            &format!(
                "{} -> complete (seen={} applied={} files={} warnings={})",
                contig, seen, applied, files, warnings
            ),
        )
    }

    pub fn failure(&mut self, stage: &str, contig: &str, reason: &str) -> Result<()> {
        self.tag("FAILURE", &format!("{} @ {}: {}", contig, stage, reason))
    }

    pub fn summary(&mut self, lines: &[String]) -> Result<()> {
        for l in lines {
            self.tag("SUMMARY", l)?;
        }
        Ok(())
    }
}

fn fmt_dur(d: Duration) -> String {
    let ms = d.as_millis();
    if ms < 1000 {
        format!("{}ms", ms)
    } else if ms < 60_000 {
        format!("{:.2}s", ms as f64 / 1000.0)
    } else {
        format!("{}m{:.0}s", ms / 60_000, (ms % 60_000) as f64 / 1000.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn warn_writes_to_warnings_file_uncapped() {
        let dir = tempdir().unwrap();
        let log = dir.path().join("run.log");
        let warn_log = dir.path().join("run.warnings.log");
        let mut lg = Logger::create(&log, Some(&warn_log), true).unwrap();

        for i in 0..500 {
            lg.warn(&format!("test message {}", i)).unwrap();
        }
        lg.flush_warnings().unwrap();

        let content = std::fs::read_to_string(&warn_log).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 500);
        assert!(lines[0].contains("test message 0"));
        assert!(lines[499].contains("test message 499"));

        // Main log must not contain them.
        let main_content = std::fs::read_to_string(&log).unwrap();
        assert!(main_content.is_empty());
    }

    #[test]
    fn warning_summary_goes_to_both_files() {
        let dir = tempdir().unwrap();
        let log = dir.path().join("run.log");
        let warn_log = dir.path().join("run.warnings.log");
        let mut lg = Logger::create(&log, Some(&warn_log), true).unwrap();

        let mut reasons: BTreeMap<&'static str, usize> = BTreeMap::new();
        reasons.insert("overlap", 42);
        reasons.insert("ref_mismatch", 7);
        lg.warning_summary("chr1", 49, &reasons).unwrap();
        lg.flush_warnings().unwrap();

        let main = std::fs::read_to_string(&log).unwrap();
        assert!(main.contains("[WARNINGS] chr1: total=49"));
        assert!(main.contains("overlap=42"));
        assert!(main.contains("ref_mismatch=7"));

        let warn = std::fs::read_to_string(&warn_log).unwrap();
        assert!(warn.contains("[WARNINGS] chr1: total=49"));
    }

    #[test]
    fn warning_summary_is_silent_when_total_zero() {
        let dir = tempdir().unwrap();
        let log = dir.path().join("run.log");
        let warn_log = dir.path().join("run.warnings.log");
        let mut lg = Logger::create(&log, Some(&warn_log), true).unwrap();
        lg.warning_summary("chr1", 0, &BTreeMap::new()).unwrap();
        lg.flush_warnings().unwrap();

        let main = std::fs::read_to_string(&log).unwrap();
        assert!(main.is_empty());
    }

    #[test]
    fn warnings_log_optional() {
        let dir = tempdir().unwrap();
        let log = dir.path().join("run.log");
        let mut lg = Logger::create(&log, None, true).unwrap();
        assert!(!lg.has_warnings_file());
        // warn() should not panic when there is no warnings file.
        lg.warn("orphan warning").unwrap();
    }
}