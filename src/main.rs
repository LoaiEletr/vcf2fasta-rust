//! Main binary entry point for the vcf2fasta tool.
//!
//! This executable simply parses command‑line arguments using `clap`
//! and delegates all processing to the library crate `vcf2fasta`.
//! The library handles the full workflow: reading the VCF, applying
//! variants, and writing FASTA files.

use anyhow::Result;
use clap::Parser;
use vcf2fasta::{cli::Args, run};

fn main() -> Result<()> {
    // Parse command‑line arguments into the `Args` struct defined in `cli.rs`.
    // `clap` automatically generates help, version, and error messages.
    let args = Args::parse();

    // Delegate the actual work to the library's `run` function.
    // All logging, parallel processing, and file I/O are handled there.
    run(args)
}