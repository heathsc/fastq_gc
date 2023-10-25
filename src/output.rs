use std::io::stdout;

use anyhow::Context;
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
    let wrt = stdout();
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
