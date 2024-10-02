use std::collections::{HashMap, VecDeque};

use crate::kmcv::Kmcv;
use serde::Serialize;

mod builder;
mod iterators;
mod target;

use builder::*;
use iterators::*;
use target::*;

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

#[derive(Serialize)]
pub struct KmerCounts<'a> {
    kmcv: &'a Kmcv,
    total_reads: u32,
    mapped_reads: u32,
    total_bases: u64,
    mapped_bases: u64,
    counts: Vec<(u32, u64)>,
    #[serde(skip_serializing)]
    vec: VecDeque<Option<TargetVec>>,
    #[serde(skip_serializing)]
    work_vec: Option<TargetVec>,
    #[serde(skip_serializing)]
    pos_hash: Option<HashMap<Target, [u32; 2]>>,
}

impl<'a> KmerCounts<'a> {
    pub fn new(kmcv: &'a Kmcv) -> Self {
        let n_targets = kmcv.n_targets();
        let kmer_length = kmcv.kmer_length();
        let counts = vec![(0, 0); n_targets];
        let mut hash_vec = VecDeque::with_capacity(kmer_length as usize);
        for _ in 0..kmer_length {
            hash_vec.push_back(Some(TargetVec::default()));
        }

        Self {
            kmcv,
            counts,
            vec: hash_vec,
            work_vec: Some(TargetVec::default()),
            pos_hash: Some(HashMap::new()),
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

    pub fn clear(&mut self) {
        for h in self.vec.iter_mut() {
            h.as_mut().unwrap().clear()
        }
        self.pos_hash.as_mut().map(|h| h.clear());
    }

    pub fn kmer_length(&self) -> usize {
        self.kmcv.kmer_length() as usize
    }

    fn add_target_pos(pos_hash: &mut HashMap<Target, [u32; 2]>, t: &Target, ix: u32) {
        let p = pos_hash.entry(*t).or_insert([ix, ix]);
        p[1] = p[1].max(ix);
    }

    fn add_contrib_and_prev(
        &self,
        kmer: Option<KmerType>,
        th: &mut TargetVec,
        prev: Option<&TargetVec>,
        pos_hash: &mut HashMap<Target, [u32; 2]>,
        ix: u32,
    ) {
        let hits = kmer.and_then(|k| self.kmcv.target_hits(k));
        match (hits, prev) {
            (Some(v), Some(p)) => {
                let itr = v.iter().map(|id| Target::new(*id));
                for x in TVMaxVec::new(p.as_slice(), itr) {
                    th.push(x);
                    Self::add_target_pos(pos_hash, x.target(), ix);
                }
            }
            (Some(v), None) => {
                for x in v {
                    let t = Target::new(*x).into_target_info();
                    th.push(t);
                    Self::add_target_pos(pos_hash, t.target(), ix);
                }
            }
            (None, Some(p)) => {
                th.extend(p.as_slice());
            }
            (None, None) => (),
        }
    }

    fn update_contribs(
        &self,
        kmer: Option<KmerType>,
        th: &TargetVec,
        th1: &TargetVec,
        work: &mut TargetVec,
        pos_hash: &mut HashMap<Target, [u32; 2]>,
        ix: u32,
    ) {
        work.clear();
        let hits = kmer.and_then(|k| self.kmcv.target_hits(k));
        match hits {
            Some(v) => {
                let itr = v.iter().map(|id| Target::new(*id));
                for x in TVAddVecMrgMax::new(&th.as_slice(), &th1.as_slice(), itr) {
                    work.push(x);
                    Self::add_target_pos(pos_hash, x.target(), ix);
                }
            }
            None => {
                for x in TVMrgMax::new(th.as_slice(), &th1.as_slice()) {
                    work.push(x);
                    Self::add_target_pos(pos_hash, x.target(), ix);
                }
            }
        }
    }

    pub fn add_target_hit(&mut self, i: usize, kmer: Option<KmerType>, ix: u32) {
        let mut pos_hash = self.pos_hash.take().unwrap();
        let mut th = self.vec[i].take().unwrap();
        th.clear();
        let prev = if i > 0 {
            self.vec[i - 1].as_ref()
        } else {
            None
        };

        // Add contributions from current kmer
        self.add_contrib_and_prev(kmer, &mut th, prev, &mut pos_hash, ix);
        self.vec[i] = Some(th);
        self.pos_hash = Some(pos_hash);
    }

    pub fn update(&mut self, kmer: Option<KmerType>, ix: u32) {
        let mut pos_hash = self.pos_hash.take().unwrap();
        let mut work = self.work_vec.take().unwrap();
        work.clear();
        let th = self.vec.pop_front().unwrap().unwrap();

        let l = self.vec.len();
        let th1 = self.vec[l - 1].as_ref().unwrap();

        if th.is_empty() {
            self.add_contrib_and_prev(kmer, &mut work, Some(th1), &mut pos_hash, ix)
        } else {
            self.update_contribs(kmer, &th, th1, &mut work, &mut pos_hash, ix)
        }
        self.vec.push_back(Some(work));
        self.work_vec = Some(th);
        self.pos_hash = Some(pos_hash);
    }

    pub fn check_map_and_update_counts(&mut self, read_length: usize) {
        let mut hit1: Option<TargetInfo> = None;
        let mut hit2: Option<TargetInfo> = None;

        let th = self.vec.back().unwrap().as_ref().unwrap();
        for t in th.as_slice().iter().filter(|t| t.is_mapped()) {
            if let Some(t1) = hit1.take() {
                if t.count() > t1.count() {
                    hit1 = Some(*t);
                    hit2 = Some(t1);
                } else {
                    hit1 = Some(t1);
                    if let Some(t2) = hit2.take() {
                        if t.count() > t2.count() {
                            hit2 = Some(*t)
                        } else {
                            hit2 = Some(t2)
                        }
                    } else {
                        hit2 = Some(*t)
                    }
                }
            } else {
                hit1 = Some(*t)
            }
        }
        self.total_reads += 1;
        self.total_bases += read_length as u64;
        if let Some(x) = match (hit1, hit2) {
            (Some(t1), Some(t2)) => {
                if t1.count() - t2.count() > 1 {
                    Some(t1)
                } else {
                    None
                }
            }
            (Some(t1), _) => Some(t1),
            _ => None,
        } {
            assert!(x.target_id() > 0);
            let p = self.pos_hash.as_ref().unwrap().get(x.target()).unwrap();
            let mapped_length = (p[1] - p[0] + self.kmcv.kmer_length() as u32) as u64;
            if mapped_length as f64 / read_length as f64 >= 0.6 {
                // eprintln!("Ack: {:?}\t{}\t{}", x, mapped_length, read_length);
                let p = &mut self.counts[x.target_id() as usize - 1];
                *p = (p.0 + 1, p.1 + mapped_length);
                self.mapped_reads += 1;
                self.mapped_bases += mapped_length;
            }
        }
    }
}

#[derive(Serialize)]
pub struct KmerWork<'a> {
    kmer_counts: KmerCounts<'a>,
    #[serde(skip_serializing)]
    builder: KmerBuilder,
}

impl<'a> KmerWork<'a> {
    pub fn new(k: &'a Kmcv) -> Self {
        let counts = KmerCounts::new(k);
        let builder = KmerBuilder::new(k.kmer_length());
        Self {
            kmer_counts: counts,
            builder,
        }
    }

    pub fn add(&mut self, rhs: &Self) {
        self.kmer_counts.add(&rhs.kmer_counts)
    }

    pub fn counts_builder_mut(&mut self) -> (&mut KmerCounts<'a>, &mut KmerBuilder) {
        (&mut self.kmer_counts, &mut self.builder)
    }
}
