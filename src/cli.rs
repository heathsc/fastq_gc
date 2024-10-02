use std::path::PathBuf;

use anyhow::Context;
use chrono::{DateTime, Local};
use clap::ArgMatches;
use compress_io::compress::CompressIo;
use regex::Regex;
use serde::Serialize;

mod cli_model;

use crate::{
    control_seq::{process_control_sequences, CSeq},
    kmcv::Kmcv,
    utils::BisulfiteType,
};

#[derive(Serialize)]
pub struct Config {
    trim: usize,
    min_qual: u8,
    threads: usize,
    bisulfite: BisulfiteType,
    control_seq: Option<CSeq>,
    output_prefix: String,
    input_files: Vec<Option<PathBuf>>,
    #[serde(skip_serializing)]
    date: DateTime<Local>,
    #[serde(skip_serializing)]
    kmcv: Option<Kmcv>,
    fli: Vec<Fli>,
}

#[derive(Default, Serialize)]
pub struct Fli {
    sample: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    barcode: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    library: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    flowcell: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    index: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    lane: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    read_end: Option<u8>,
}

impl Fli {
    fn from_input(p: &Option<PathBuf>, m: &ArgMatches) -> Self {
        let mut fli = Self::default();
        if let Some(fname) = p
            .as_deref()
            .and_then(|p| p.file_name().and_then(|p| p.to_str()))
        {
            let re = Regex::new(r"^([^_]+)_(\d+)_([^_]+)_?(\d*)").unwrap();
            if let Some(cap) = re.captures(fname) {
                fli.flowcell = cap.get(1).map(|m| m.as_str().to_owned());
                fli.index = cap.get(3).map(|m| m.as_str().to_owned());
                fli.lane = match cap.get(2).map(|m| m.as_str().parse::<u8>()) {
                    Some(Ok(x)) => Some(x),
                    _ => None,
                };
                fli.read_end = match cap.get(4).map(|m| m.as_str().parse::<u8>()) {
                    Some(Ok(x)) => Some(x),
                    _ => None,
                }
            }
        }

        fli.sample = m.get_one::<String>("sample").map(|s| s.to_owned());
        fli.barcode = m.get_one::<String>("barcode").map(|s| s.to_owned());
        fli.library = m.get_one::<String>("library").map(|s| s.to_owned());
        if let Some(flowcell) = m.get_one::<String>("flowcell") {
            fli.flowcell = Some(flowcell.to_owned())
        }
        if let Some(lane) = m.get_one::<u8>("lane") {
            fli.lane = Some(*lane)
        }
        if let Some(index) = m.get_one::<String>("index") {
            fli.index = Some(index.to_owned())
        }
        if let Some(end) = m.get_one::<u8>("read_end") {
            fli.read_end = Some(*end)
        }
        fli
    }
    pub fn read_end(&self) -> Option<u8> {
        self.read_end
    }
    pub fn make_file_name(&self, prefix: &str) -> PathBuf {
        let mut s = prefix.to_owned();
        if let Some(fc) = self.flowcell.as_deref() {
            s.push('_');
            s.push_str(fc);
        }
        if let Some(lane) = self.lane {
            s.push_str(format!("_{}", lane).as_str());
        }
        if let Some(index) = self.index.as_deref() {
            s.push('_');
            s.push_str(index);
        }
        if let Some(end) = self.read_end {
            s.push_str(format!("_{}", end).as_str());
        }
        let mut p = PathBuf::from(&s);
        p.set_extension("json");
        p
    }
}

impl Config {
    pub fn trim(&self) -> usize {
        self.trim
    }
    pub fn min_qual(&self) -> u8 {
        self.min_qual
    }
    pub fn threads(&self) -> usize {
        self.threads
    }
    pub fn input_files(&self) -> &[Option<PathBuf>] {
        &self.input_files
    }
    pub fn output_prefix(&self) -> &str {
        &self.output_prefix
    }
    pub fn control_seq(&self) -> Option<&CSeq> {
        self.control_seq.as_ref()
    }
    pub fn fli(&self) -> &[Fli] {
        &self.fli
    }
    pub fn bisulfite_type(&self) -> BisulfiteType {
        self.bisulfite
    }
    pub fn date(&self) -> &DateTime<Local> {
        &self.date
    }
    pub fn kmcv(&self) -> Option<&Kmcv> {
        self.kmcv.as_ref()
    }
}

pub fn handle_cli() -> anyhow::Result<Config> {
    let c = cli_model::cli_model();
    let m = c.get_matches();
    super::utils::init_log(&m);

    let input_files: Vec<Option<PathBuf>> = m
        .get_many::<PathBuf>("input")
        .map(|p| p.map(|x| Some(x.to_owned())).collect())
        .unwrap_or(vec![None]);

    let output_prefix = m
        .get_one::<String>("output_prefix")
        .map(|p| p.to_owned())
        .unwrap();

    let threads = m
        .get_one::<u64>("threads")
        .map(|x| *x as usize)
        .expect("Missing default threads option");

    let min_qual = m
        .get_one::<u8>("min_qual")
        .copied()
        .expect("Missing default min_qual option");

    let trim = m
        .get_one::<usize>("trim")
        .copied()
        .expect("Missing default trim option");

    let bisulfite = m
        .get_one::<BisulfiteType>("bisulfite_type")
        .copied()
        .unwrap_or_else(|| {
            if m.get_flag("bisulfite") {
                BisulfiteType::Forward
            } else {
                BisulfiteType::None
            }
        });

    let mut bs_error = false;
    let fli: Vec<_> = input_files
        .iter()
        .map(|p| {
            let f = Fli::from_input(p, &m);
            if f.read_end.is_none() && !matches!(bisulfite, BisulfiteType::None) {
                bs_error = true;
            }
            f
        })
        .collect();

    if bs_error {
        return Err(anyhow!(
            "Cannot infer read end so bisulfite modes cannot be used"
        ));
    }

    let control_seq = match m.get_one::<PathBuf>("control") {
        Some(p) => Some(process_control_sequences(p, !bisulfite.is_none())?),
        None => None,
    };

    let kmcv = match m.get_one::<PathBuf>("kmers") {
        Some(p) => {
            let mut rdr = CompressIo::new()
                .path(p)
                .bufreader()
                .with_context(|| "Could not open kmer file for input")?;

            debug!("Opened kmer file for input");
            Some(
                Kmcv::read(&mut rdr)
                    .with_context(|| format!("Could not read kmer file {}", p.display()))?,
            )
        }
        None => None,
    };

    Ok(Config {
        input_files,
        output_prefix,
        control_seq,
        kmcv,
        trim,
        min_qual,
        threads,
        bisulfite,
        date: Local::now(),
        fli,
    })
}
