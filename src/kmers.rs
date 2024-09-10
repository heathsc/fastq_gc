use std::collections::{hash_map::Entry, HashMap, VecDeque};

use crate::kmcv::Kmcv;
use serde::Serialize;

pub type KmerType = u32;

// We use the upper 4 bits (61-64) to encode the number of hits, leaving the lower 60 bits to
// encode the index into the array of hits (region indexes)
pub struct KmerIdxHits(u64);

impl KmerIdxHits {
    const MAX_IDX: u64 = (1u64 << 60) - 1;

    #[inline]
    pub fn new(idx: usize, n_hits: u8) -> Self {
        let idx = idx as u64;
        assert!(idx <= Self::MAX_IDX && n_hits <= 8);
        Self(idx | ((n_hits as u64) << 60))
    }
    #[inline]
    pub fn hits(&self) -> u8 {
        (self.0 >> 60) as u8
    }
    #[inline]
    pub fn idx(&self) -> usize {
        (self.0 & Self::MAX_IDX) as usize
    }
}

#[derive(Default)]
struct TargetHash {
    hash: HashMap<u32, u32>, // Keep track of which targets have been tagged for a read
}

impl TargetHash {
    fn clear(&mut self) {
        self.hash.clear()
    }

    fn add_hit(&mut self, v: &[u32]) {
        for idx in v {
            *self.hash.entry(*idx).or_insert(0) += 1
        }
    }
}

#[derive(Serialize)]
pub struct KmerCounts<'a> {
    kmcv: &'a Kmcv,
    total_reads: u32,
    mapped_reads: u32,
    total_bases: u64,
    mapped_bases: u64,
    counts: Vec<(u32, u64)>,
    #[serde(skip_serializing)]
    hash_vec: VecDeque<Option<TargetHash>>,
}

impl<'a> KmerCounts<'a> {
    pub fn new(kmcv: &'a Kmcv) -> Self {
        let n_targets = kmcv.n_targets();
        let kmer_length = kmcv.kmer_length();
        let counts = vec![(0, 0); n_targets];
        let mut hash_vec = VecDeque::with_capacity(kmer_length as usize);
        for _ in 0..kmer_length {
            hash_vec.push_back(Some(TargetHash::default()));
        }

        Self {
            kmcv,
            counts,
            hash_vec,
            total_reads: 0,
            mapped_reads: 0,
            total_bases: 0,
            mapped_bases: 0,
        }
    }

    pub fn add(&mut self, rhs: &Self) {
        for (p, q) in self.counts.iter_mut().zip(rhs.counts.iter()) {
            *p = (p.0 + q.0, p.1 + q.1)
        }
        self.total_bases += rhs.total_bases;
        self.mapped_reads += rhs.mapped_reads;
        self.total_reads += rhs.total_reads;
        self.mapped_bases += rhs.mapped_bases;
    }

    pub fn clear_hash(&mut self) {
        for h in self.hash_vec.iter_mut() {
            h.as_mut().unwrap().clear()
        }
    }

    pub fn kmer_length(&self) -> usize {
        self.kmcv.kmer_length() as usize
    }

    fn add_current_contrib(&mut self, kmer: Option<KmerType>, th: &mut TargetHash) {
        if let Some(k) = kmer {
            if let Some(v) = self.kmcv.target_hits(k) {
                th.add_hit(v);
            }
        }
    }

    fn get_max(th: &mut TargetHash, th1: &TargetHash) {
        for (a, b) in th1.hash.iter() {
            match th.hash.entry(*a) {
                Entry::Vacant(mut e) => {
                    e.insert(*b);
                }
                Entry::Occupied(mut e) => {
                    let p = e.get_mut();
                    let x = *p;
                    *p = x.max(*b)
                }
            }
        }
    }

    pub fn add_target_hit(&mut self, i: usize, kmer: Option<KmerType>) {
        let mut th = self.hash_vec[i].take().unwrap();
        th.hash.clear();

        // Add contributions from current kmer
        self.add_current_contrib(kmer, &mut th);

        // Get max of previous position and current position
        if i > 0 {
            let th1 = self.hash_vec[i - 1].as_ref().unwrap();
            Self::get_max(&mut th, th1);
        }
        self.hash_vec[i] = Some(th);
    }

    pub fn update(&mut self, kmer: Option<KmerType>) {
        let mut th = self.hash_vec.pop_front().unwrap().unwrap();

        let l = self.hash_vec.len();
        // Add contributions from current kmer
        self.add_current_contrib(kmer, &mut th);

        let th1 = self.hash_vec[l - 1].as_ref().unwrap();
        Self::get_max(&mut th, th1);
        self.hash_vec.push_back(Some(th));
    }

    pub fn check_map_and_update_counts(&mut self, read_length: usize) {
        let mut hit1: Option<(u32, u32)> = None;
        let mut hit2: Option<(u32, u32)> = None;

        let th = self.hash_vec.back().unwrap().as_ref().unwrap();
        for (t, n) in th.hash.iter().filter(|(&t, &n)| t > 0 && n > 1) {
            if let Some((t1, n1)) = hit1.take() {
                if *n > n1 {
                    hit1 = Some((*t, *n));
                    hit2 = Some((t1, n1));
                } else {
                    hit1 = Some((t1, n1));
                    if *n < n1 {
                        if let Some((t2, n2)) = hit2.take() {
                            if *n > n2 {
                                hit2 = Some((*t, *n))
                            } else {
                                hit2 = Some((t2, n2))
                            }
                        } else {
                            hit2 = Some((*t, *n))
                        }
                    }
                }
            } else {
                hit1 = Some((*t, *n))
            }
        }
        if let Some(x) = match (hit1, hit2) {
            (Some((t1, n1)), Some((_, n2))) => {
                if n1 > n2 {
                    Some(t1)
                } else {
                    None
                }
            }
            (Some((t1, _)), _) => Some(t1),
            _ => None,
        } {
            assert!(x > 0);
            let p = &mut self.counts[x as usize - 1];
            *p = (p.0 + 1, p.1 + read_length as u64);
            self.mapped_reads += 1;
            self.mapped_bases += read_length as u64;
        }
        self.total_reads += 1;
        self.total_bases += read_length as u64;
    }
}

pub struct KmerBuilder {
    kmer: KmerType,
    valid: KmerType,
    mask: KmerType,
    valid_mask: KmerType,
    kmer_length: usize,
}

#[derive(Serialize)]
pub struct KmerWork<'a> {
    counts: KmerCounts<'a>,
    #[serde(skip_serializing)]
    builder: KmerBuilder,
}

impl<'a> KmerWork<'a> {
    pub fn new(k: &'a Kmcv) -> Self {
        let counts = KmerCounts::new(k);
        let builder = KmerBuilder::new(k.kmer_length());
        Self { counts, builder }
    }

    pub fn add(&mut self, rhs: &Self) {
        self.counts.add(&rhs.counts)
    }

    pub fn counts_builder_mut(&mut self) -> (&mut KmerCounts<'a>, &mut KmerBuilder) {
        (&mut self.counts, &mut self.builder)
    }
}

impl KmerBuilder {
    pub fn new(kmer_length: u8) -> Self {
        let k = kmer_length as usize;
        let nb = KmerType::BITS as usize;
        assert!(k + k <= nb, "Kmer length too large for KmerType");
        assert!(k > 0, "Kmer length cannot be zero");

        const ZERO: KmerType = 0;
        debug!("mask = {:x}", (!ZERO) >> (nb - k - k));
        Self {
            kmer: 0,
            valid: 0,
            mask: (!ZERO) >> (nb - k - k),
            valid_mask: (!ZERO) >> (nb - k),
            kmer_length: k,
        }
    }

    #[inline]
    pub fn clear(&mut self) {
        self.valid = 0;
        self.kmer = 0;
    }

    pub fn add_base(&mut self, base: u8) {
        let x = base as usize;
        let valid: KmerType = if x < 4 { 1 } else { 0 };
        self.kmer = ((self.kmer << 2) & self.mask) | ((x & 3) as KmerType);
        self.valid = ((self.valid << 1) & self.valid_mask) | valid;
    }

    pub fn make_from_slice<F>(&mut self, bases: &[u8], f: F) -> Option<KmerType>
    where
        F: Fn(&u8) -> KmerType,
    {
        self.clear();
        if bases.len() < self.kmer_length {
            None
        } else {
            for x in bases[..self.kmer_length].iter().map(f) {
                if x > 3 {
                    return None;
                }
                self.kmer = (self.kmer << 2) | x;
            }
            Some(self.kmer)
        }
    }

    /// Check if kmer is composed entirely of valid (i.e., A, C, G, T) bases
    #[inline]
    pub fn kmer(&self) -> Option<KmerType> {
        if self.valid == self.valid_mask {
            Some(self.kmer)
        } else {
            None
        }
    }

    #[inline]
    pub fn kmer_length(&self) -> usize {
        self.kmer_length
    }
}
