use std::io::Error;
use std::{
    io::Write,
    path::{Path, PathBuf},
};

use regex::Regex;
mod cli_model;

use crate::control_seq::{process_control_sequences, ControlSeq};

pub struct Config {
    trim: usize,
    min_qual: u8,
    threads: usize,
    bisulfite: bool,
    control_seq: Option<ControlSeq>,
    input_file: Option<PathBuf>,
    fli: Fli,
}

#[derive(Default)]
pub struct Fli {
    sample: Option<String>,
    library: Option<String>,
    flowcell: Option<String>,
    index: Option<String>,
    lane: Option<u8>,
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

    pub fn json_output<W: Write>(
        &self,
        wrt: &mut W,
        indent: usize,
        mut in_list: bool,
    ) -> Result<bool, Error> {
        let mut chk = |w: &mut W| -> Result<(), Error> {
            if !in_list {
                in_list = true;
                Ok(())
            } else {
                writeln!(w, ",")
            }
        };
        if let Some(s) = self.sample.as_deref() {
            chk(wrt)?;
            write!(wrt, "{:indent$}\"sample\": \"{}\"", " ", s)?;
        }
        if let Some(l) = self.library.as_deref() {
            chk(wrt)?;
            write!(wrt, "{:indent$}\"library\": \"{}\"", " ", l)?;
        }
        if let Some(fc) = self.flowcell.as_deref() {
            chk(wrt)?;
            write!(wrt, "{:indent$}\"flowcell\": \"{}\"", " ", fc)?;
        }
        if let Some(l) = self.lane {
            chk(wrt)?;
            write!(wrt, "{:indent$}\"lane\": \"{}\"", " ", l)?;
        }
        if let Some(ix) = self.index.as_deref() {
            chk(wrt)?;
            write!(wrt, "{:indent$}\"index\": \"{}\"", " ", ix)?;
        }
        if let Some(e) = self.read_end {
            chk(wrt)?;
            write!(wrt, "{:indent$}\"read_end\": \"{}\"", " ", e)?;
        }
        Ok(in_list)
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
    pub fn bisulfite(&self) -> bool {
        self.bisulfite
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

    let bisulfite = m.get_flag("bisulfite");

    let control_seq = match m.get_one::<PathBuf>("control") {
        Some(p) => Some(process_control_sequences(p, bisulfite)?),
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

    Ok(Config {
        input_file,
        control_seq,
        trim,
        min_qual,
        threads,
        bisulfite,
        fli,
    })
}
