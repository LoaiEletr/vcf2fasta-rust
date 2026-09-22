//! FASTA file writer with configurable line width.
//!
//! Supports two ways to open a file:
//!
//! * `FastaWriter::create` — creates a new file and writes the header.
//! * `FastaWriter::reopen` — appends to an existing file, restoring the
//!   line-wrap column so multi-write sequences stay correctly formatted.
//!
//! `reopen` exists so the output manager can keep a bounded number of
//! file handles open at any moment (see `output_manager.rs`): writers are
//! closed when evicted and reopened on demand without corrupting line
//! wrapping or duplicating the header.

use anyhow::{bail, Context, Result};
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

pub struct FastaWriter {
    writer: BufWriter<File>,
    line_width: usize,
    column: usize,
}

impl FastaWriter {
    /// Create a new FASTA file and write its header.
    pub fn create(path: &Path, sequence_name: &str, line_width: usize) -> Result<Self> {
        if line_width == 0 {
            bail!("FASTA line width must be greater than 0");
        }
        let file = File::create(path)
            .with_context(|| format!("could not create output {}", path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, ">{sequence_name}")?;
        writer.flush()?;
        Ok(Self {
            writer,
            line_width,
            column: 0,
        })
    }

    /// Reopen an existing FASTA in append mode, restoring the line-wrap
    /// column so continued writes wrap correctly.
    ///
    /// Does **not** write a header. The caller must have previously
    /// created the file via [`FastaWriter::create`] or reopened it at
    /// least once before.
    pub fn reopen(path: &Path, line_width: usize, column: usize) -> Result<Self> {
        if line_width == 0 {
            bail!("FASTA line width must be greater than 0");
        }
        if column > line_width {
            bail!(
                "FASTA reopen column {} exceeds line width {}",
                column,
                line_width
            );
        }
        let file = OpenOptions::new()
            .append(true)
            .open(path)
            .with_context(|| format!("could not reopen FASTA {}", path.display()))?;
        Ok(Self {
            writer: BufWriter::new(file),
            line_width,
            column,
        })
    }

    /// Current column position (0..line_width) in the current line.
    pub fn column(&self) -> usize {
        self.column
    }

    /// Flush buffered data to the OS. Does **not** add a trailing newline.
    /// Used when a writer is evicted from the LRU cache and will be
    /// reopened later.
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }

    /// Write a sequence segment, wrapping lines as needed.
    pub fn write_seq(&mut self, mut seq: &[u8]) -> Result<()> {
        while !seq.is_empty() {
            let room = self.line_width - self.column;
            let n = room.min(seq.len());
            self.writer.write_all(&seq[..n])?;
            self.column += n;
            seq = &seq[n..];
            if self.column == self.line_width {
                self.writer.write_all(b"\n")?;
                self.column = 0;
            }
        }
        Ok(())
    }

    /// Finalise: add a trailing newline if the last line is partial, then
    /// flush.
    pub fn finish(mut self) -> Result<()> {
        if self.column != 0 {
            self.writer.write_all(b"\n")?;
        }
        self.writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    fn read_file(path: &Path) -> String {
        let mut file = File::open(path).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        contents
    }

    #[test]
    fn create_ok() {
        let temp = NamedTempFile::new().unwrap();
        let path = temp.path();
        let writer = FastaWriter::create(path, "seq1", 80).unwrap();
        assert_eq!(writer.line_width, 80);
        assert_eq!(writer.column, 0);
        drop(writer);
        assert_eq!(read_file(path), ">seq1\n");
    }

    #[test]
    fn create_fails_on_zero_line_width() {
        let temp = NamedTempFile::new().unwrap();
        assert!(FastaWriter::create(temp.path(), "seq1", 0).is_err());
    }

    #[test]
    fn write_seq_without_wrapping() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = FastaWriter::create(temp.path(), "seq1", 80).unwrap();
        writer.write_seq(b"ACGT").unwrap();
        writer.finish().unwrap();
        assert_eq!(read_file(temp.path()), ">seq1\nACGT\n");
    }

    #[test]
    fn write_seq_multiple_lines() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = FastaWriter::create(temp.path(), "seq1", 4).unwrap();
        writer.write_seq(b"ACGTACGTACGT").unwrap();
        writer.finish().unwrap();
        assert_eq!(read_file(temp.path()), ">seq1\nACGT\nACGT\nACGT\n");
    }

    #[test]
    fn write_seq_partial_last_line_extra_newline() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = FastaWriter::create(temp.path(), "seq1", 4).unwrap();
        writer.write_seq(b"ACGTACGTA").unwrap();
        writer.finish().unwrap();
        assert_eq!(read_file(temp.path()), ">seq1\nACGT\nACGT\nA\n");
    }

    #[test]
    fn write_seq_multiple_calls_continue_line() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = FastaWriter::create(temp.path(), "seq1", 4).unwrap();
        writer.write_seq(b"AC").unwrap();
        writer.write_seq(b"GT").unwrap();
        writer.write_seq(b"AC").unwrap();
        writer.finish().unwrap();
        assert_eq!(read_file(temp.path()), ">seq1\nACGT\nAC\n");
    }

    #[test]
    fn write_seq_empty_does_nothing() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = FastaWriter::create(temp.path(), "seq1", 80).unwrap();
        writer.write_seq(b"").unwrap();
        writer.finish().unwrap();
        assert_eq!(read_file(temp.path()), ">seq1\n");
    }

    #[test]
    fn reopen_continues_mid_line() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fa");

        // First writer: create and write "AC", leaving column = 2.
        {
            let mut w = FastaWriter::create(&path, "seq1", 4).unwrap();
            w.write_seq(b"AC").unwrap();
            let col = w.column();
            assert_eq!(col, 2);
            w.flush().unwrap();
        }

        // Second writer: reopen at column 2 and write "GTAC".
        // Expected final: ACGT then AC on separate lines.
        {
            let mut w = FastaWriter::reopen(&path, 4, 2).unwrap();
            w.write_seq(b"GTAC").unwrap();
            w.finish().unwrap();
        }

        assert_eq!(read_file(&path), ">seq1\nACGT\nAC\n");
    }

    #[test]
    fn reopen_at_column_zero_adds_no_extra_newline() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fa");

        // First writer: create and write exactly one full line.
        {
            let mut w = FastaWriter::create(&path, "seq1", 4).unwrap();
            w.write_seq(b"ACGT").unwrap();
            // column is now 0 because of the implicit newline
            assert_eq!(w.column(), 0);
            w.flush().unwrap();
        }

        // Second writer: reopen at column 0 and write another line.
        {
            let mut w = FastaWriter::reopen(&path, 4, 0).unwrap();
            w.write_seq(b"TGCA").unwrap();
            w.finish().unwrap();
        }

        assert_eq!(read_file(&path), ">seq1\nACGT\nTGCA\n");
    }

    #[test]
    fn reopen_rejects_column_larger_than_line_width() {
        let temp = NamedTempFile::new().unwrap();
        assert!(FastaWriter::reopen(temp.path(), 4, 5).is_err());
    }
}