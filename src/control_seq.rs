use anyhow::Context;
use compress_io::compress::CompressIo;
use std::{collections::HashMap, io::BufRead, ops::Deref, path::Path};

use serde::Serialize;

use crate::{process::BisulfiteStrand, reader::FastQRecord};

type CsMask = u64;

const NORMAL_KMER_LEN: u32 = 20;
const BISULFITE_KMER_LEN: u32 = 27;

#[derive(Serialize)]
pub enum CSeq {
    Regular(ControlSeq<NORMAL_KMER_LEN>),
    Bisulfite(ControlSeq<BISULFITE_KMER_LEN>),
}

impl CSeq {
    pub fn filter(&self, rec: &FastQRecord, strand: BisulfiteStrand) -> Option<usize> {
        match self {
            Self::Regular(cs) => cs.filter(rec, strand),
            Self::Bisulfite(cs) => cs.filter(rec, strand),
        }
    }

    pub fn seq_ids(&self) -> &[String] {
        match self {
            Self::Regular(cs) => cs.seq_ids(),
            Self::Bisulfite(cs) => cs.seq_ids(),
        }
    }
}
#[derive(Serialize)]
pub struct ControlSeq<const K: u32> {
    seq_ids: Vec<String>,
    #[serde(skip_serializing)]
    hash: HashMap<u64, CsMask>,
}

impl<const K: u32> ControlSeq<K> {
    const KMER_MASK: u64 = (1 << (K * 2)) - 1;
    const REV_SHIFT: u32 = (K - 1) * 2;
    pub fn new(n_seq: u32) -> Self {
        assert!(
            n_seq <= CsMask::BITS,
            "Too many control sequences - change CsMask type"
        );
        assert!(u64::BITS >= 2 * K, "Kmer length too large for u64");
        Self {
            seq_ids: Vec::with_capacity(n_seq as usize),
            hash: HashMap::new(),
        }
    }

    pub fn seq_ids(&self) -> &[String] {
        &self.seq_ids
    }

    pub fn n_seq(&self) -> usize {
        self.seq_ids.len()
    }

    fn gen_hash(&mut self, seq: &Seq, ix: usize, rev_ix: usize) {
        if seq.len() >= K as usize {
            let mut n = 0;
            let mut key: u64 = 0;
            let mut rev_key: u64 = 0;
            let msk: CsMask = 1 << ix;
            let rev_msk: CsMask = 1 << rev_ix;

            let get_x = |c: Base| {
                let x = (c as u64) & 3;
                (x, ((x + 2) & 3) << Self::REV_SHIFT)
            };

            // Initialize first window
            for c in seq.iter().take(K as usize) {
                if !c.is_gap() {
                    n += 1;
                }
                let (x, rev_x) = get_x(*c);
                key = (key << 2) | x;
                rev_key = (rev_key >> 2) | rev_x;
            }
            // Process the full sequence
            let i1 = K as usize - 1;
            for s in seq.windows(K as usize) {
                // Update hash table
                if n == K {
                    add_key(&mut self.hash, key, msk);
                    add_key(&mut self.hash, rev_key, rev_msk);
                }
                // Update keys and counts
                key = (key << 2) & Self::KMER_MASK;
                rev_key >>= 2;
                let (c, c1) = (s[0], s[i1]);
                if c1.is_gap() {
                    if !c.is_gap() {
                        n -= 1
                    }
                } else {
                    let (x, rev_x) = get_x(c1);
                    if c.is_gap() {
                        n += 1
                    }
                    key |= x;
                    rev_key |= rev_x;
                }
            }
        }
    }

    fn gen_bs_hash(&mut self, seq: &Seq, ix: usize, n_seq: usize) {
        if seq.len() >= K as usize {
            let mut n = 0;
            let mut keys = [0u64; 4];
            let msk: [CsMask; 4] = [
                1 << ix,
                1 << (ix + n_seq),
                1 << (ix + 2 * n_seq),
                1 << (ix + 3 * n_seq),
            ];

            let shift_key = |k: &mut u64| *k = (*k << 2) & Self::KMER_MASK;
            let rev_shift_key = |k: &mut u64| *k >>= 2;

            let shift_keys = |keys: &mut [u64; 4]| {
                shift_key(&mut keys[0]);
                rev_shift_key(&mut keys[1]);
                shift_key(&mut keys[2]);
                rev_shift_key(&mut keys[3]);
            };

            let get_x = |c: Base| match c {
                Base::A => [0, 2 << Self::REV_SHIFT, 0, 2 << Self::REV_SHIFT],
                Base::C => [2, 3 << Self::REV_SHIFT, 1, 0],
                Base::T => [2, 0, 2, 0],
                Base::G => [3, 2 << Self::REV_SHIFT, 0, 1 << Self::REV_SHIFT],
                _ => [0, 0, 0, 0],
            };

            let update_keys = |keys: &mut [u64; 4], c: Base| {
                let x = get_x(c);
                for (k, x) in keys.iter_mut().zip(&x) {
                    *k |= *x
                }
            };

            // Initialize first window
            for c in seq.iter().take(K as usize) {
                shift_keys(&mut keys);
                if !c.is_gap() {
                    n += 1;
                    update_keys(&mut keys, *c);
                }
            }
            // Process the full sequence
            let i1 = K as usize - 1;
            for s in seq.windows(K as usize) {
                // Update hash table
                if n == K {
                    for (k, m) in keys.iter().zip(&msk) {
                        add_key(&mut self.hash, *k, *m)
                    }
                }
                // Update keys and counts
                shift_keys(&mut keys);

                let (c, c1) = (s[0], s[i1]);
                if c1.is_gap() {
                    if !c.is_gap() {
                        n -= 1
                    }
                } else {
                    if c.is_gap() {
                        n += 1
                    }
                    update_keys(&mut keys, c1)
                }
            }
        }
    }

    pub fn filter(&self, rec: &FastQRecord, strand: BisulfiteStrand) -> Option<usize> {
        let seq = rec.seq();
        let mut i = 0;
        let mut msk: CsMask = !0;
        let mut n_match = 0;
        let mut start_match = 0;
        let mut best = (0, 0);
        let f = match strand {
            BisulfiteStrand::None => Base::from_u8,
            BisulfiteStrand::C2T => Base::from_u8_c2t,
            BisulfiteStrand::G2A => Base::from_u8_g2a,
        };
        let mut curr_key: Option<u64> = None;
        while i < seq.len() - K as usize {
            let (key, ix) = {
                // If curr_key is not None then we are just moving by one base to test the next
                // kmer.  In this case it is quicker to modify the previous key than make a new
                // key from scratch
                if let Some(mut key) = curr_key.take() {
                    // Find the new base to add to the key
                    let b = f(seq[i + K as usize - 1]);
                    if b.is_gap() {
                        // Not a base - skip ahead by KMER_LEN bases
                        (None, K as usize)
                    } else {
                        // Shift key and add new base
                        key = ((key << 2) & Self::KMER_MASK) | ((b as u64) & 3);
                        (Some(key), K as usize)
                    }
                } else {
                    // curr_key is none so we make a new key from scratch
                    Self::make_key(&seq[i..i + K as usize], f)
                }
            };
            // key exists if we have no 'N's in the KMER_LEN bases starting from seq[i]
            if let Some(key) = key {
                // Check hash
                if let Some(k) = self.hash.get(&key) {
                    msk &= *k;
                    if msk == 0 {
                        n_match = 0;
                        i = start_match + 1;
                        msk = !0;
                    } else {
                        if n_match == 0 {
                            start_match = i;
                        }
                        n_match += 1;
                        if n_match > best.0 {
                            best = (n_match, msk)
                        }
                        i += ix;
                    }
                } else {
                    // Hash not found so we will check the position one base over
                    // We will keep the previous key and modify it for the next iteration
                    curr_key = Some(key);
                    i += 1
                }
            } else {
                // An 'N' was found, so we skip to one past the N and try again
                i += ix
            }
        }
        if best.0 > 1 {
            let mut x = best.1;
            assert_ne!(x, 0);
            let mut i = 0;
            loop {
                if (x & 1) == 1 {
                    break;
                }
                x >>= 1;
                i += 1
            }
            Some(i % self.n_seq())
        } else {
            None
        }
    }
    fn make_key(seq: &[u8], f: fn(u8) -> Base) -> (Option<u64>, usize) {
        let mut key = 0;
        for (i, c) in seq.iter().enumerate() {
            let b = f(*c);
            if b.is_gap() {
                return (None, i + 1);
            }
            let x = (b as u64) & 3;
            key = ((key << 2) | x) & Self::KMER_MASK;
        }
        (Some(key), seq.len())
    }
}

#[derive(Default, Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
#[repr(u8)]
pub enum Base {
    A = 0,
    C,
    T,
    G,
    N,
    #[default]
    Other,
}

impl Base {
    pub fn from_u8(c: u8) -> Self {
        match c {
            b'A' | b'a' => Self::A,
            b'C' | b'c' => Self::C,
            b'G' | b'g' => Self::G,
            b'T' | b't' => Self::T,
            b'N' | b'n' => Self::N,
            _ => Self::Other,
        }
    }

    pub fn from_u8_c2t(c: u8) -> Self {
        match c {
            b'A' | b'a' => Self::A,
            b'C' | b'c' | b'T' | b't' => Self::T,
            b'G' | b'g' => Self::G,
            b'N' | b'n' => Self::N,
            _ => Self::Other,
        }
    }
    pub fn from_u8_g2a(c: u8) -> Self {
        match c {
            b'A' | b'a' | b'G' | b'g' => Self::A,
            b'C' | b'c' => Self::C,
            b'T' | b't' => Self::T,
            b'N' | b'n' => Self::N,
            _ => Self::Other,
        }
    }

    pub fn is_gap(&self) -> bool {
        self >= &Self::N
    }
}

#[derive(Debug, Default)]
pub struct Seq(Vec<Base>);

impl Deref for Seq {
    type Target = [Base];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Default)]
struct BuildContig {
    id: String,
    seq: Seq,
}

fn read_sequence<P: AsRef<Path>>(file: P) -> anyhow::Result<Vec<BuildContig>> {
    debug!("Reading control sequences from {}", file.as_ref().display());
    let mut rdr = CompressIo::new()
        .path(file)
        .bufreader()
        .with_context(|| "Could not open control sequence file")?;

    let mut buf = String::new();
    let mut ctg = BuildContig::default();
    let mut contigs = Vec::new();
    while rdr
        .read_line(&mut buf)
        .with_context(|| "Error reading from control sequence file")?
        > 0
    {
        let b = buf.trim();
        if let Some(s) = b.strip_prefix('>') {
            if !ctg.seq.is_empty() {
                contigs.push(ctg);
            }
            ctg = BuildContig::default();
            ctg.id.push_str(s);
            debug!("Reading control sequence {s}");
        } else {
            for c in b.as_bytes() {
                ctg.seq.0.push(Base::from_u8(*c))
            }
        }
        buf.clear();
    }
    if !ctg.seq.is_empty() {
        contigs.push(ctg);
    }
    if contigs.is_empty() {
        Err(anyhow!("No sequences read in from control sequence file"))
    } else {
        Ok(contigs)
    }
}
fn add_key(hash: &mut HashMap<u64, CsMask>, key: u64, msk: CsMask) {
    let e = hash.entry(key).or_insert(0);
    *e |= msk
}

pub fn process_control_sequences<P: AsRef<Path>>(file: P, bisulfite: bool) -> anyhow::Result<CSeq> {
    let mut contigs = read_sequence(file)?;

    if bisulfite {
        process_control_sequences_for_bisulfite(contigs).map(CSeq::Bisulfite)
    } else {
        let n_seq = contigs.len();
        // Each sequence is present in forward and reverse versions
        let n_seq2 = n_seq * 2;
        let mut cseq = ControlSeq::new(n_seq2 as u32);

        // Generate hashes
        for (ix, ctg) in contigs.drain(..).enumerate() {
            let BuildContig { id, seq } = ctg;
            debug!("Generating hashtable for {id}");
            cseq.gen_hash(&seq, ix, ix + n_seq);
            cseq.seq_ids.push(id);
        }
        Ok(CSeq::Regular(cseq))
    }
}

fn process_control_sequences_for_bisulfite(
    mut contigs: Vec<BuildContig>,
) -> anyhow::Result<ControlSeq<BISULFITE_KMER_LEN>> {
    let n_seq = contigs.len();
    // Each sequence is present in C2T and A2G and forward and reverse versions
    let n_seq2 = n_seq * 4;
    let mut cseq = ControlSeq::new(n_seq2 as u32);

    // Generate hashes
    for (ix, ctg) in contigs.drain(..).enumerate() {
        let BuildContig { id, seq } = ctg;
        debug!("Generating bisulfite hashtable for {id}");
        cseq.gen_bs_hash(&seq, ix, n_seq);
        cseq.seq_ids.push(id);
    }
    Ok(cseq)
}
