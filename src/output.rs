use anyhow::Context;
use compress_io::compress::CompressIo;
use serde::Serialize;

use crate::{cli::Config, process::ProcessResults};

#[derive(Serialize)]
struct JsonReport<'a, 'b, 'c> {
    program: &'static str,
    version: &'static str,
    date: String,
    #[serde(flatten)]
    cfg: &'a Config,
    max_read_length: usize,
    #[serde(flatten)]
    res: &'b ProcessResults<'c>,
}

pub fn output_results(cfg: &Config, res: &ProcessResults) -> anyhow::Result<()> {
    let wrt = CompressIo::new()
        .opt_path(cfg.output_file())
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
