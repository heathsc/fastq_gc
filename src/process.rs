use anyhow::Context;
use std::{collections::HashMap, ops::AddAssign};

use crossbeam_channel::{Receiver, Sender};
use serde::{
    ser::{SerializeMap, Serializer},
    Serialize,
};

use crate::utils::BisulfiteType;
use crate::{
    cli::Config,
    kmcv::Kmcv,
    kmers::{KmerType, KmerWork},
    reader::{Buffer, FastQRecord},
};

#[derive(Default, Debug, Copy, Clone)]
pub struct BaseCounts {
    cts: [usize; 6], // A, C, T, G, N, Other
}

const SER_LIST: [(&str, usize); 6] = [
    ("A", 0),
    ("C", 1),
    ("G", 3),
    ("T", 2),
    ("N", 4),
    ("Other", 5),
];
impl Serialize for BaseCounts {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut n = 4;
        for c in &self.cts[4..] {
            if *c > 0 {
                n += 1
            }
        }

        let mut map = serializer.serialize_map(Some(n))?;
        for (k, v) in &SER_LIST[..4] {
            map.serialize_entry(*k, &self.cts[*v])?;
        }
        for (k, v) in &SER_LIST[4..] {
            let x = self.cts[*v];
            if x > 0 {
                map.serialize_entry(*k, &x)?;
            }
        }
        map.end()
    }
}

impl BaseCounts {
    pub fn clear(&mut self) {
        self.cts.fill(0)
    }
}

impl AddAssign for BaseCounts {
    fn add_assign(&mut self, rhs: Self) {
        for (x, y) in self.cts.iter_mut().zip(rhs.cts) {
            *x += y
        }
    }
}

#[derive(Copy, Clone, Eq, PartialOrd, PartialEq, Hash)]
pub struct GcHistKey(u32, u32);

impl Serialize for GcHistKey {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&format!("{}:{}", self.0, self.1))
    }
}

pub struct ControlSeqCounts<'a> {
    counts: Vec<usize>,
    ids: &'a [String],
}

impl<'a> ControlSeqCounts<'a> {
    fn new(ids: &'a [String]) -> Self {
        Self {
            ids,
            counts: vec![0; ids.len()],
        }
    }

    fn add(&mut self, rhs: &Self) {
        for (p, q) in self.counts.iter_mut().zip(rhs.counts.iter()) {
            *p += *q
        }
    }

    fn add_read(&mut self, ix: usize) {
        self.counts[ix] += 1
    }
}

impl<'a> Serialize for ControlSeqCounts<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut map = serializer.serialize_map(Some(self.ids.len()))?;
        for (k, v) in self.ids.iter().zip(self.counts.iter()) {
            map.serialize_entry(k, v)?;
        }
        map.end()
    }
}

struct PerPositionCounts {
    v: Vec<BaseCounts>,
    trim: usize,
}

impl Serialize for PerPositionCounts {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut map = serializer.serialize_map(Some(self.v.len()))?;
        for (ix, c) in self.v.iter().enumerate() {
            let i = self.trim + ix + 1;
            map.serialize_entry(&i, c)?;
        }
        map.end()
    }
}

impl PerPositionCounts {
    fn add_obs(&mut self, m: usize, i: usize) {
        // Resize per pos vec if necessary
        if i >= self.v.len() {
            self.v.resize_with(i + 1, BaseCounts::default)
        }
        self.v[i].cts[m] += 1;
    }

    fn add(&mut self, rhs: &Self) {
        if rhs.v.len() > self.v.len() {
            self.v.resize_with(rhs.v.len(), BaseCounts::default)
        }
        for (x, y) in self.v.iter_mut().zip(&rhs.v) {
            *x += *y
        }
    }

    fn max_read_length(&self) -> usize {
        self.v.len() + self.trim
    }
}

#[derive(Serialize)]
pub struct ProcessResults<'a, 'b> {
    #[serde(skip_serializing_if = "Option::is_none")]
    control_seq_counts: Option<ControlSeqCounts<'a>>,
    cts: BaseCounts,
    per_pos_cts: PerPositionCounts,
    gc_hash: HashMap<GcHistKey, usize>,
    #[serde(skip_serializing)]
    temp_cts: BaseCounts,
    #[serde(skip_serializing)]
    base_map: [u8; 256],
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(flatten)]
    kmer_work: Option<KmerWork<'b>>,
}

impl<'a, 'b> ProcessResults<'a, 'b> {
    pub fn new(control_seq_ids: Option<&'a [String]>, trim: usize, kmcv: Option<&'b Kmcv>) -> Self {
        let gc_hash = HashMap::new();
        let cts = BaseCounts::default();
        let temp_cts = BaseCounts::default();
        let per_pos_cts = PerPositionCounts {
            v: Vec::with_capacity(256),
            trim,
        };
        let mut base_map = [6; 256];
        for (i, c) in [b'A', b'C', b'T', b'G', b'N'].into_iter().enumerate() {
            base_map[c as usize] = i as u8; // Upper case
            base_map[(c | 32) as usize] = i as u8; // Lower case
        }
        let control_seq_counts = control_seq_ids.map(ControlSeqCounts::new);
        let kmer_work = kmcv.map(KmerWork::new);

        Self {
            gc_hash,
            cts,
            temp_cts,
            per_pos_cts,
            control_seq_counts,
            base_map,
            kmer_work,
        }
    }

    pub fn max_read_length(&self) -> usize {
        self.per_pos_cts.max_read_length()
    }
    fn add_obs(&mut self, b: u8, i: usize) {
        let m = self.base_map[b as usize] as usize;
        self.cts.cts[m] += 1;
        self.per_pos_cts.add_obs(m, i)
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
            let ct = GcHistKey(a as u32, b as u32);
            let e = self.gc_hash.entry(ct).or_insert(0);
            *e += 1
        }
    }
}

impl<'a, 'b> AddAssign for ProcessResults<'a, 'b> {
    fn add_assign(&mut self, rhs: Self) {
        // Add hash
        for (key, val) in rhs.gc_hash {
            let e = self.gc_hash.entry(key).or_insert(0);
            *e += val
        }

        // And counts
        self.cts += rhs.cts;
        self.per_pos_cts.add(&rhs.per_pos_cts);

        // plus control sequence counts
        if let Some(cs) = self.control_seq_counts.as_mut() {
            let cs1 = rhs.control_seq_counts.as_ref().unwrap();
            cs.add(cs1)
        }

        // plus kmer counts
        if let Some(kw) = self.kmer_work.as_mut() {
            if let Some(kw1) = rhs.kmer_work.as_ref() {
                kw.add(kw1)
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

fn process_kmers(rec: &FastQRecord, kw: &mut KmerWork, base_map: &[u8; 256], trim: usize) {
    let (kc, kb) = kw.counts_builder_mut();
    kc.clear();

    /*
    for v in rec.seq()[trim..].chunks_exact(kb.kmer_length()) {
        if let Some(kmer) = kb.make_from_slice(v, |b| base_map[*b as usize] as KmerType) {
            kc.add_target_hit(kmer)
        }
    }*/

    let kmer_length = kc.kmer_length();

    if kmer_length <= rec.seq().len() {
        kb.clear();
        let l = rec.seq().len();
        for (i, b) in rec.seq().iter().copied().enumerate() {
            kb.add_base(base_map[b as usize]);
            if i >= kmer_length {
                let j = i - kmer_length;
                if j < kmer_length {
                    kc.add_target_hit(j, kb.kmer())
                } else {
                    kc.update(kb.kmer())
                }
            }
        }
        kc.check_map_and_update_counts(rec.seq().len())
    }
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

#[derive(Debug, Copy, Clone)]
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
    let bisulfite_type = cfg.bisulfite_type();
    let file_index = b.file_index();
    let read_end = cfg.fli()[file_index].read_end();

    let strand = match (bisulfite_type, read_end) {
        (BisulfiteType::None, _) => Some(BisulfiteStrand::None),
        (BisulfiteType::Forward, Some(1)) | (BisulfiteType::Reverse, Some(2)) => {
            Some(BisulfiteStrand::C2T)
        }
        (BisulfiteType::Forward, Some(2)) | (BisulfiteType::Reverse, Some(1)) => {
            Some(BisulfiteStrand::G2A)
        }
        _ => None,
    };

    assert!(
        res.kmer_work.is_none() || matches!(bisulfite_type, BisulfiteType::None),
        "Cannot track kmer usage with bisulfite data"
    );

    for rec in b.fastq() {
        let rec = rec?;
        base_counts_from_record(&rec, trim, min_qual, res);
        let st = strand.unwrap_or_else(|| BisulfiteStrand::from_counts(&res.temp_cts));

        if let Some(cs) = cseq {
            if let Some(ix) = cs.filter(&rec, st) {
                res.control_seq_counts.as_mut().unwrap().add_read(ix);
                continue;
            }
        }
        process_record(&rec, trim, min_qual, st, res);
        if let Some(kw) = res.kmer_work.as_mut() {
            process_kmers(&rec, kw, &res.base_map, trim)
        }
    }

    Ok(())
}

pub fn process_thread(
    cfg: &Config,
    ix: usize,
    rx: Receiver<Buffer>,
    sx: Sender<Buffer>,
) -> anyhow::Result<Vec<(usize, ProcessResults)>> {
    debug!("Starting up process thread {ix}");
    let mut res_vec = Vec::new();

    let new_res = || {
        ProcessResults::new(
            cfg.control_seq().map(|c| c.seq_ids()),
            cfg.trim(),
            cfg.kmcv(),
        )
    };

    while let Ok(mut b) = rx.recv() {
        trace!("Process thread {ix} received new block");
        if let Some((idx, _)) = res_vec.last() {
            if b.file_index() != *idx {
                debug!(
                    "Process thread {ix} started processing blocks from file {}",
                    b.file_index()
                );
                res_vec.push((b.file_index(), new_res()));
            }
        } else {
            debug!(
                "Process thread {ix} started processing blocks from file {}",
                b.file_index()
            );
            res_vec.push((b.file_index(), new_res()));
        }
        let res = &mut res_vec.last_mut().unwrap().1;
        process_buffer(cfg, &b, res)
            .with_context(|| format!("Process thread {}: Error parsing input buffer", ix))?;
        trace!("Process thread {ix} finished processing block; sending empty block back to reader");
        b.clear();

        // Send empty buffer to reader to be refilled.
        // Ignore errors when sending - this will happen when the reader process has finished
        let _ = sx.send(b);
    }
    debug!("Closing down process thread {ix}");

    Ok(res_vec)
}
