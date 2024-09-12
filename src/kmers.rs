use std::{cmp::Ordering, collections::VecDeque};

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
struct TargetVec {
    v: Vec<(u32, u32)>, // Keep track of which targets have been tagged for a read
}

impl TargetVec {
    fn clear(&mut self) {
        self.v.clear()
    }
}

struct TVMaxVec<'a, 'b> {
    left: &'a [(u32, u32)],
    right: &'b [u32],
    left_ix: usize,
    right_ix: usize,
}

impl<'a, 'b> TVMaxVec<'a, 'b> {
    fn new(left: &'a [(u32, u32)], right: &'b [u32]) -> Self {
        Self {
            left,
            right,
            left_ix: 0,
            right_ix: 0,
        }
    }
}
impl<'a, 'b> Iterator for TVMaxVec<'a, 'b> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        let l = self.left.get(self.left_ix);
        let r = self.right.get(self.right_ix);

        let which = match (l, r) {
            (Some((l, _)), Some(r)) => Some(l.cmp(r)),
            (Some((_, _)), None) => Some(Ordering::Less),
            (None, Some(_)) => Some(Ordering::Greater),
            (None, None) => None,
        };

        match which {
            Some(Ordering::Less) => {
                self.left_ix += 1;
                l.copied()
            }
            Some(Ordering::Equal) => {
                self.left_ix += 1;
                self.right_ix += 1;

                l.copied()
            }
            Some(Ordering::Greater) => {
                self.right_ix += 1;
                r.map(|a| (*a, 1))
            }
            None => None,
        }
    }
}

struct TVMrgMax<'a, 'b> {
    left: &'a [(u32, u32)],
    right: &'b [(u32, u32)],
    left_ix: usize,
    right_ix: usize,
}

impl<'a, 'b> TVMrgMax<'a, 'b> {
    fn new(left: &'a [(u32, u32)], right: &'b [(u32, u32)]) -> Self {
        Self {
            left,
            right,
            left_ix: 0,
            right_ix: 0,
        }
    }
}
impl<'a, 'b> Iterator for TVMrgMax<'a, 'b> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        let l = self.left.get(self.left_ix);
        let r = self.right.get(self.right_ix);

        let (which, n2) = match (l, r) {
            (Some((l, _)), Some((r, n2))) => (Some(l.cmp(r)), *n2),
            (Some((_, _)), None) => (Some(Ordering::Less), 0),
            (None, Some(_)) => (Some(Ordering::Greater), 0),
            (None, None) => (None, 0),
        };

        match which {
            Some(Ordering::Less) => {
                self.left_ix += 1;
                l.copied()
            }
            Some(Ordering::Equal) => {
                self.left_ix += 1;
                self.right_ix += 1;

                l.map(|(a, b)| (*a, *b.max(&n2)))
            }
            Some(Ordering::Greater) => {
                self.right_ix += 1;
                r.copied()
            }
            None => None,
        }
    }
}

struct TVAddVecMrgMax<'a, 'b, 'c> {
    left: &'a [(u32, u32)],
    right: &'b [(u32, u32)],
    v: &'c [u32],
    left_ix: usize,
    right_ix: usize,
    v_ix: usize,
}

impl<'a, 'b, 'c> TVAddVecMrgMax<'a, 'b, 'c> {
    fn new(left: &'a [(u32, u32)], right: &'b [(u32, u32)], v: &'c [u32]) -> Self {
        Self {
            left,
            right,
            v,
            left_ix: 0,
            right_ix: 0,
            v_ix: 0,
        }
    }
}

enum AVOrd {
    AllEqual(u32),
    LMin,
    RMin,
    VMin,
    LRMin(u32),
    LVMin,
    RVMin,
    None,
}

impl<'a, 'b, 'c> Iterator for TVAddVecMrgMax<'a, 'b, 'c> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        let l = self.left.get(self.left_ix);
        let r = self.right.get(self.right_ix);
        let v = self.v.get(self.v_ix);

        let which = match (l, r, v) {
            (Some((l, n1)), Some((r, n2)), Some(v)) => match l.cmp(r) {
                Ordering::Equal => match l.cmp(v) {
                    Ordering::Equal => AVOrd::AllEqual((n1 + 1).max(*n2)),
                    Ordering::Less => AVOrd::LRMin(*n1.max(n2)),
                    Ordering::Greater => AVOrd::VMin,
                },
                Ordering::Less => match l.cmp(v) {
                    Ordering::Equal => AVOrd::LVMin,
                    Ordering::Less => AVOrd::LMin,
                    Ordering::Greater => AVOrd::VMin,
                },
                Ordering::Greater => match r.cmp(v) {
                    Ordering::Equal => AVOrd::RVMin,
                    Ordering::Less => AVOrd::RMin,
                    Ordering::Greater => AVOrd::VMin,
                },
            },
            (Some((l, n1)), Some((r, n2)), None) => match l.cmp(r) {
                Ordering::Equal => AVOrd::LRMin(*n1.max(n2)),
                Ordering::Less => AVOrd::LMin,
                Ordering::Greater => AVOrd::RMin,
            },
            (Some((l, _)), None, Some(v)) => match l.cmp(v) {
                Ordering::Equal => AVOrd::LVMin,
                Ordering::Less => AVOrd::LMin,
                Ordering::Greater => AVOrd::VMin,
            },
            (None, Some((r, _)), Some(v)) => match r.cmp(v) {
                Ordering::Equal => AVOrd::RVMin,
                Ordering::Less => AVOrd::RMin,
                Ordering::Greater => AVOrd::VMin,
            },
            (Some(_), None, None) => AVOrd::LMin,
            (None, Some(_), None) => AVOrd::RMin,
            (None, None, Some(_)) => AVOrd::VMin,
            (None, None, None) => AVOrd::None,
        };

        match which {
            AVOrd::AllEqual(n) => {
                self.left_ix += 1;
                self.right_ix += 1;
                self.v_ix += 1;
                l.map(|(a, _)| (*a, n))
            }
            AVOrd::LRMin(n) => {
                self.left_ix += 1;
                self.right_ix += 1;
                l.map(|(a, _)| (*a, n))
            }
            AVOrd::LVMin => {
                self.left_ix += 1;
                self.v_ix += 1;
                l.map(|(a, b)| (*a, b + 1))
            }
            AVOrd::RVMin => {
                self.right_ix += 1;
                self.v_ix += 1;
                r.copied()
            }
            AVOrd::LMin => {
                self.left_ix += 1;
                l.copied()
            }
            AVOrd::RMin => {
                self.right_ix += 1;
                r.copied()
            }
            AVOrd::VMin => {
                self.v_ix += 1;
                v.map(|a| (*a, 1))
            }
            AVOrd::None => None,
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
    hash_vec: VecDeque<Option<TargetVec>>,
    #[serde(skip_serializing)]
    work_vec: Option<TargetVec>,
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
            hash_vec,
            work_vec: Some(TargetVec::default()),
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
        for h in self.hash_vec.iter_mut() {
            h.as_mut().unwrap().clear()
        }
    }

    pub fn kmer_length(&self) -> usize {
        self.kmcv.kmer_length() as usize
    }

    fn add_contrib_and_prev(
        &self,
        kmer: Option<KmerType>,
        th: &mut TargetVec,
        prev: Option<&TargetVec>,
    ) {
        let hits = kmer.and_then(|k| self.kmcv.target_hits(k));

        match (hits, prev) {
            (Some(v), Some(p)) => {
                for x in TVMaxVec::new(&p.v, v) {
                    th.v.push(x)
                }
            }
            (Some(v), None) => {
                for x in v {
                    th.v.push((*x, 1))
                }
            }
            (None, Some(p)) => {
                th.v.extend_from_slice(&p.v);
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
    ) {
        work.v.clear();
        let hits = kmer.and_then(|k| self.kmcv.target_hits(k));

        match hits {
            Some(v) => {
                for x in TVAddVecMrgMax::new(&th.v, &th1.v, v) {
                    work.v.push(x)
                }
            }
            None => {
                for x in TVMrgMax::new(&th.v, &th1.v) {
                    work.v.push(x)
                }
            }
        }
    }

    pub fn add_target_hit(&mut self, i: usize, kmer: Option<KmerType>) {
        let mut th = self.hash_vec[i].take().unwrap();
        th.clear();
        let prev = if i > 0 {
            self.hash_vec[i - 1].as_ref()
        } else {
            None
        };

        // Add contributions from current kmer
        self.add_contrib_and_prev(kmer, &mut th, prev);

        self.hash_vec[i] = Some(th);
    }

    pub fn update(&mut self, kmer: Option<KmerType>) {
        let mut work = self.work_vec.take().unwrap();
        work.clear();
        let th = self.hash_vec.pop_front().unwrap().unwrap();

        let l = self.hash_vec.len();
        let th1 = self.hash_vec[l - 1].as_ref().unwrap();

        if th.v.is_empty() {
            self.add_contrib_and_prev(kmer, &mut work, Some(th1))
        } else {
            self.update_contribs(kmer, &th, th1, &mut work)
        }

        self.hash_vec.push_back(Some(work));
        self.work_vec = Some(th);
    }

    pub fn check_map_and_update_counts(&mut self, read_length: usize) {
        let mut hit1: Option<(u32, u32)> = None;
        let mut hit2: Option<(u32, u32)> = None;

        let th = self.hash_vec.back().unwrap().as_ref().unwrap();
        for (t, n) in th.v.iter().filter(|(t, n)| *t > 0 && *n > 1) {
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

    /// Check if kmer is composed entirely of valid (i.e., A, C, G, T) bases
    #[inline]
    pub fn kmer(&self) -> Option<KmerType> {
        if self.valid == self.valid_mask {
            Some(self.kmer)
        } else {
            None
        }
    }
}
