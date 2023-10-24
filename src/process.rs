use anyhow::Context;
use std::{
    collections::HashMap,
    io::{Error, Write},
    ops::AddAssign,
};

use crossbeam_channel::{Receiver, Sender};

use crate::{
    cli::Config,
    reader::{Buffer, FastQRecord},
};

#[derive(Default, Debug, Copy, Clone)]
pub struct BaseCounts {
    cts: [usize; 6], // A, C, T, G, N, Other
}

impl BaseCounts {
    pub fn clear(&mut self) {
        self.cts.fill(0)
    }

    pub fn json_output<W: Write>(&self, wrt: &mut W) -> Result<(), Error> {
        write!(
            wrt,
            "{{ \"A\": {}, \"C\": {}, \"G\": {}, \"T\": {}",
            self.cts[0], self.cts[1], self.cts[3], self.cts[2]
        )?;
        if self.cts[4] > 0 {
            write!(wrt, ", \"N\": {}", self.cts[4])?
        }
        if self.cts[5] > 0 {
            write!(wrt, ", \"Other\": {}", self.cts[5])?
        }
        write!(wrt, " }}")
    }
}

impl AddAssign for BaseCounts {
    fn add_assign(&mut self, rhs: Self) {
        for (x, y) in self.cts.iter_mut().zip(rhs.cts) {
            *x += y
        }
    }
}

pub struct ProcessResults<'a> {
    gc_hash: HashMap<[u32; 2], usize>,
    cts: BaseCounts,
    temp_cts: BaseCounts,
    per_pos_cts: Vec<BaseCounts>,
    control_seq_counts: Option<Vec<usize>>,
    control_seq_ids: Option<&'a [String]>,
    base_map: [u8; 256],
}

impl<'a> ProcessResults<'a> {
    pub fn new(control_seq_ids: Option<&'a [String]>) -> Self {
        let gc_hash = HashMap::new();
        let cts = BaseCounts::default();
        let temp_cts = BaseCounts::default();
        let per_pos_cts = Vec::with_capacity(256);
        let mut base_map = [6; 256];
        for (i, c) in [b'A', b'C', b'T', b'G', b'N'].into_iter().enumerate() {
            base_map[c as usize] = i as u8; // Upper case
            base_map[(c | 32) as usize] = i as u8; // Lower case
        }
        let control_seq_counts = control_seq_ids.map(|v| vec![0; v.len()]);

        Self {
            gc_hash,
            cts,
            temp_cts,
            per_pos_cts,
            control_seq_counts,
            control_seq_ids,
            base_map,
        }
    }

    fn add_obs(&mut self, b: u8, i: usize) {
        // Resize per pos vec if necessary
        if i >= self.per_pos_cts.len() {
            self.per_pos_cts.resize_with(i + 1, BaseCounts::default)
        }
        let m = self.base_map[b as usize] as usize;
        self.cts.cts[m] += 1;
        self.per_pos_cts[i].cts[m] += 1;
    }

    fn clear_tmp_counts(&mut self) {
        self.temp_cts.clear()
    }

    fn add_to_tmp_counts(&mut self, b: u8) {
        let m = self.base_map[b as usize] as usize;
        self.temp_cts.cts[m] += 1;
    }

    fn add_gc(&mut self, a: usize, b: usize) {
        if a + b > 0 {
            let ct = [a as u32, b as u32];
            let e = self.gc_hash.entry(ct).or_insert(0);
            *e += 1
        }
    }

    pub fn json_output<W: Write>(
        &self,
        wrt: &mut W,
        trim: usize,
        indent: usize,
        in_list: bool,
    ) -> Result<(), Error> {
        if in_list {
            writeln!(wrt, ",")?;
        }
        if let Some(cs_cts) = self.control_seq_counts.as_deref() {
            writeln!(wrt, "{:indent$}\"control_seq_counts\": {{", " ")?;
            let i = indent + 3;
            let seq_ids = self.control_seq_ids.unwrap();
            let mut first = true;
            for (s, c) in seq_ids.iter().zip(cs_cts) {
                if first {
                    first = false
                } else {
                    writeln!(wrt, ",")?;
                }
                write!(wrt, "{:i$}\"{s}\": {c}", " ")?;
            }
            writeln!(wrt, "\n{:indent$}}},", " ")?;
        }
        write!(wrt, "{:indent$}\"base_counts\": ", " ")?;
        self.cts.json_output(wrt)?;
        let c0 = BaseCounts::default();
        let v = vec![c0; trim];
        let mut itr = v.iter().chain(self.per_pos_cts.iter());
        if let Some(c) = itr.next() {
            writeln!(wrt, ",\n{:indent$}\"per_position_base_counts\": [", " ")?;
            let i = indent + 3;
            write!(wrt, "{:i$}", " ")?;
            c.json_output(wrt)?;
            for c in itr {
                write!(wrt, ",\n{:i$}", " ")?;
                c.json_output(wrt)?
            }
            write!(wrt, "\n{:indent$}]", " ")?;
        }
        let mut itr = self.gc_hash.iter();
        if let Some((k, v)) = itr.next() {
            let i = indent + 3;
            write!(
                wrt,
                ",\n{:indent$}\"gc_counts\": {{\n{:i$}\"{}:{}\": {v}",
                " ", " ", k[0], k[1]
            )?;
            for (k, v) in itr {
                write!(wrt, ",\n{:i$}\"{}:{}\": {v}", " ", k[0], k[1])?
            }
            write!(wrt, "\n{:indent$}}}", " ")
        } else {
            Ok(())
        }
    }
}

impl<'a> AddAssign for ProcessResults<'a> {
    fn add_assign(&mut self, rhs: Self) {
        // Add hash
        for (key, val) in rhs.gc_hash {
            let e = self.gc_hash.entry(key).or_insert(0);
            *e += val
        }
        // And counts
        self.cts += rhs.cts;
        if rhs.per_pos_cts.len() > self.per_pos_cts.len() {
            self.per_pos_cts
                .resize_with(rhs.per_pos_cts.len(), BaseCounts::default)
        }
        for (x, y) in self.per_pos_cts.iter_mut().zip(&rhs.per_pos_cts) {
            *x += *y
        }
        // And control sequence counts
        if let Some(cs) = self.control_seq_counts.as_deref_mut() {
            // Will panic if both do not have counts, or if counts vec are different sizes
            let cs1 = rhs.control_seq_counts.as_deref().unwrap();
            assert_eq!(cs.len(), cs1.len());
            for (p, q) in cs.iter_mut().zip(cs1) {
                *p += *q
            }
        }
    }
}
fn process_record(
    rec: &FastQRecord,
    trim: usize,
    min_qual: u8,
    strand: BisulfiteStrand,
    res: &mut ProcessResults,
) {
    for (i, (s, _)) in rec.seq()[trim..]
        .iter()
        .zip(&rec.qual()[trim..])
        .enumerate()
        .filter(|(_, (_, q))| **q >= min_qual)
    {
        res.add_obs(*s, i)
    }
    let cts = &res.temp_cts.cts;
    let (a, b) = match strand {
        BisulfiteStrand::None => (cts[0] + cts[2], cts[1] + cts[3]), // A + T, C + G
        BisulfiteStrand::C2T => (cts[0], cts[3]),                    // A , G
        BisulfiteStrand::G2A => (cts[2], cts[1]),                    // T , C
    };
    res.add_gc(a, b)
}

fn base_counts_from_record(rec: &FastQRecord, trim: usize, min_qual: u8, res: &mut ProcessResults) {
    res.clear_tmp_counts();
    for (s, _) in rec.seq()[trim..]
        .iter()
        .zip(&rec.qual()[trim..])
        .filter(|(_, q)| **q >= min_qual)
    {
        res.add_to_tmp_counts(*s)
    }
}

#[derive(Copy, Clone)]
pub enum BisulfiteStrand {
    None,
    C2T,
    G2A,
}

impl BisulfiteStrand {
    fn from_counts(c: &BaseCounts) -> Self {
        let ct = &c.cts;
        // Compare f(T) - f(C) against f(A) - f(G)
        if ct[2] - ct[1] > ct[0] - ct[3] {
            Self::C2T
        } else {
            Self::G2A
        }
    }
}
fn process_buffer(cfg: &Config, b: &Buffer, res: &mut ProcessResults) -> anyhow::Result<()> {
    let trim = cfg.trim();
    let min_qual = cfg.min_qual();
    let cseq = cfg.control_seq();
    let bisulfite = cfg.bisulfite();
    for rec in b.fastq() {
        let rec = rec?;
        base_counts_from_record(&rec, trim, min_qual, res);
        let strand = if bisulfite {
            BisulfiteStrand::from_counts(&res.temp_cts)
        } else {
            BisulfiteStrand::None
        };

        if let Some(cs) = cseq {
            if let Some(ix) = cs.filter(&rec, strand) {
                res.control_seq_counts.as_mut().unwrap()[ix] += 1;
                continue;
            }
        }
        process_record(&rec, trim, min_qual, strand, res)
    }
    Ok(())
}

pub fn process_thread(
    cfg: &Config,
    ix: usize,
    rx: Receiver<Buffer>,
    sx: Sender<Buffer>,
) -> anyhow::Result<ProcessResults> {
    debug!("Starting up process thread {ix}");
    let mut res = ProcessResults::new(cfg.control_seq().map(|c| c.seq_ids()));

    while let Ok(mut b) = rx.recv() {
        trace!("Process thread {ix} received new block");
        process_buffer(cfg, &b, &mut res)
            .with_context(|| format!("Process thread {}: Error parsing input buffer", ix))?;
        trace!("Process thread {ix} finished processing block; sending empty block back to reader");
        b.clear();

        // Send empty buffer to reader to be refilled.
        // Ignore errors when sending - this will happen when the reader process has finished
        let _ = sx.send(b);
    }
    debug!("Closing down process thread {ix}");

    Ok(res)
}
