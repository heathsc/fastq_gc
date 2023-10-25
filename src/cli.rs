use std::path::{Path, PathBuf};

use chrono::{DateTime, Local};
use regex::Regex;
use serde::Serialize;

mod cli_model;

use crate::{
    control_seq::{process_control_sequences, ControlSeq},
    utils::BisulfiteType,
};

#[derive(Serialize)]
pub struct Config {
    trim: usize,
    min_qual: u8,
    threads: usize,
    bisulfite: BisulfiteType,
    control_seq: Option<ControlSeq>,
    input_file: Option<PathBuf>,
    #[serde(skip_serializing)]
    date: DateTime<Local>,
    fli: Fli,
}

#[derive(Default, Serialize)]
pub struct Fli {
    sample: Option<String>,
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
    fn from_input(p: Option<&Path>) -> Self {
        let mut fli = Self::default();
        if let Some(fname) = p.and_then(|p| p.file_name()).and_then(|p| p.to_str()) {
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
        fli
    }
    pub fn read_end(&self) -> Option<u8> {
        self.read_end
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
    pub fn input_file(&self) -> Option<&Path> {
        self.input_file.as_deref()
    }
    pub fn control_seq(&self) -> Option<&ControlSeq> {
        self.control_seq.as_ref()
    }
    pub fn fli(&self) -> &Fli {
        &self.fli
    }
    pub fn bisulfite_type(&self) -> BisulfiteType {
        self.bisulfite
    }
    pub fn date(&self) -> &DateTime<Local> {
        &self.date
    }
}

pub fn handle_cli() -> anyhow::Result<Config> {
    let c = cli_model::cli_model();
    let m = c.get_matches();
    super::utils::init_log(&m);

    let input_file = m.get_one::<PathBuf>("input").map(|p| p.to_owned());

    let mut fli = Fli::from_input(input_file.as_deref());

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

    let control_seq = match m.get_one::<PathBuf>("control") {
        Some(p) => Some(process_control_sequences(p, !bisulfite.is_none())?),
        None => None,
    };

    fli.sample = m.get_one::<String>("sample").map(|s| s.to_owned());
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

    if matches!(bisulfite, BisulfiteType::Forward | BisulfiteType::Reverse)
        && fli.read_end.is_none()
    {
        Err(anyhow!("Can not determine read end from input file name.  Either use --read-end option to specify end, or change bisulfite type using --bisulfite-type to none or nonstranded"))
    } else {
        Ok(Config {
            input_file,
            control_seq,
            trim,
            min_qual,
            threads,
            bisulfite,
            date: Local::now(),
            fli,
        })
    }
}
