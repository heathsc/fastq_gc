use std::path::Path;

use anyhow::Context;
use compress_io::{compress::CompressIo, compress_type::CompressType};
use serde::Serialize;

use crate::{
    cli::{Config, Fli},
    control_seq::CSeq,
    process::ProcessResults,
    utils::BisulfiteType,
};

#[derive(Serialize)]
pub struct JsonConfig<'a> {
    trim: usize,
    min_qual: u8,
    threads: usize,
    bisulfite: BisulfiteType,
    control_seq: Option<&'a CSeq>,
    output_prefix: &'a str,
    input_file: Option<&'a Path>,
    fli: &'a Fli,
}

impl<'a> JsonConfig<'a> {
    fn from_config(cfg: &'a Config, ix: usize) -> Self {
        Self {
            trim: cfg.trim(),
            min_qual: cfg.min_qual(),
            threads: cfg.threads(),
            bisulfite: cfg.bisulfite_type(),
            control_seq: cfg.control_seq(),
            output_prefix: cfg.output_prefix(),
            input_file: cfg.input_files()[ix].as_deref(),
            fli: &cfg.fli()[ix],
        }
    }
}

#[derive(Serialize)]
struct JsonReport<'a, 'b, 'c> {
    program: &'static str,
    version: &'static str,
    date: String,
    max_read_length: usize,
    #[serde(flatten)]
    cfg: JsonConfig<'a>,
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
        cfg: JsonConfig::from_config(cfg, idx),
        res,
        max_read_length: res.max_read_length(),
    };

    serde_json::to_writer_pretty(wrt, &jr)
        .with_context(|| "Error writing out JSON file with results")
}
