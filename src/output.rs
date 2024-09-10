use anyhow::Context;
use compress_io::{compress::CompressIo, compress_type::CompressType};
use serde::Serialize;

use crate::{cli::Config, process::ProcessResults};

/*
fn output_counts(kc: &KmerCounts, targets: &[Target]) -> anyhow::Result<()> {
    let mut wrt = CompressIo::new()
        .path("target_counts.gz")
        .bufwriter()
        .with_context(|| "Could not open output kmer file")?;

    for (i, p) in kc.counts().iter().enumerate() {
        let cov = p.1 as f64 / targets[i].size() as f64;
        writeln!(wrt, "{i}\t{:.2}\t{}", cov, p.0)?
    }
    Ok(())
}

 */

#[derive(Serialize)]
struct JsonReport<'a, 'b, 'c> {
    program: &'static str,
    version: &'static str,
    date: String,
    max_read_length: usize,
    #[serde(flatten)]
    cfg: &'a Config,
    #[serde(flatten)]
    res: &'b ProcessResults<'c, 'a>,
}

pub fn output_results(cfg: &Config, res: &ProcessResults, idx: usize) -> anyhow::Result<()> {
    let fname = cfg.fli()[idx].make_file_name(cfg.output_prefix());
    let wrt = CompressIo::new()
        .path(fname)
        .ctype(CompressType::Gzip)
        .bufwriter()
        .with_context(|| "Could not output output file")?;

    let jr = JsonReport {
        program: env!("CARGO_PKG_NAME"),
        version: env!("CARGO_PKG_VERSION"),
        date: cfg.date().to_rfc2822(),
        cfg,
        res,
        max_read_length: res.max_read_length(),
    };

    serde_json::to_writer_pretty(wrt, &jr)
        .with_context(|| "Error writing out JSON file with results")
}
