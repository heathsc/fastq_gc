use std::{
    cmp::Ordering,
    iter::Peekable,
};

use super::target::*;

/// Merge sort max iterator where left is a slice of [TargetInfo] sorted by target id
/// and right is a slice of targets sorted by target ids.
/// The targets in right have an implied count of 1, so if a target is only present
/// in right we create a TargetInfo with a count of 1.
///
/// If a target id is present in both left and right, the maximum of the two counts is taken. Since the count
/// for targets in right is always 1 and in left the count is at least 1, we simply take the left
/// entry and discard the right entry.
pub(super) struct TVMaxVec<'a, I: Iterator<Item = Target>> {
    left: &'a [TargetInfo],
    right: Peekable<I>,
    left_ix: usize,
}

impl<'a, I: Iterator<Item = Target>> TVMaxVec<'a, I> {
    pub(super) fn new(left: &'a [TargetInfo], right: I) -> Self {
        let right = right.peekable();
        Self {
            left,
            right,
            left_ix: 0,
        }
    }
}
impl<'a,  I: Iterator<Item = Target>> Iterator for TVMaxVec<'a, I> {
    type Item = TargetInfo;

    fn next(&mut self) -> Option<Self::Item> {
        let l = self.left.get(self.left_ix);
        let r = self.right.peek();

        let which = match (l, r) {
            (Some(l), Some(r)) => Some(l.target_id().cmp(&r.id())),
            (Some(_), None) => Some(Ordering::Less),
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
                let _ = self.right.next();
                l.copied()
            }
            Some(Ordering::Greater) => {
                let r = self.right.next();
                r.map(|a| a.into_target_info())
            }
            None => None,
        }
    }
}

/// Merge sort max iterator where left and right are slices of [TargetInfo] sorted by target id
/// If a target id is present in both slices, the maximum of the two counts is taken.
pub(super) struct TVMrgMax<'a, 'b> {
    left: &'a [TargetInfo],
    right: &'b [TargetInfo],
    left_ix: usize,
    right_ix: usize,
}

impl<'a, 'b> TVMrgMax<'a, 'b> {
    pub(super) fn new(left: &'a [TargetInfo], right: &'b [TargetInfo]) -> Self {
        Self {
            left,
            right,
            left_ix: 0,
            right_ix: 0,
        }
    }
}
impl<'a, 'b> Iterator for TVMrgMax<'a, 'b> {
    type Item = TargetInfo;

    fn next(&mut self) -> Option<Self::Item> {
        let l = self.left.get(self.left_ix);
        let r = self.right.get(self.right_ix);

        let which = match (l, r) {
            (Some(l), Some(r)) => Some(l.cmp_target_id(r)),
            (Some(_), None) => Some(Ordering::Less),
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

                l.as_ref().unwrap().max(r.as_ref().unwrap())
            }
            Some(Ordering::Greater) => {
                self.right_ix += 1;
                r.copied()
            }
            None => None,
        }
    }
}

/// Merge sort max iterator where left and right are slices of [TargetInfo] sorted by target id
/// and v is a sorted slice of target ids with implicit count of one
/// The resulting TargetInfo has a count which is hte maximum of (left + v) and right
pub(super) struct TVAddVecMrgMax<'a, 'b, I: Iterator<Item = Target>> {
    left: &'a [TargetInfo],
    right: &'b [TargetInfo],
    v: Peekable<I>,
    left_ix: usize,
    right_ix: usize,
}

impl<'a, 'b, I: Iterator<Item = Target>> TVAddVecMrgMax<'a, 'b, I> {
    pub(super) fn new(left: &'a [TargetInfo], right: &'b [TargetInfo], v: I) -> Self {
        let v = v.peekable();
        Self {
            left,
            right,
            v,
            left_ix: 0,
            right_ix: 0,
        }
    }
}

/// Possible ordering of l, r and v. Note that Some(x) always counts as lower than None
enum AVOrd {
    AllEqual(TargetInfo), // l, r and v entries have the same target id. The target parameter is the maximum of (l+v) and r
    LMin,          // l entry has the lowest target id
    RMin,          // r entry has the lowest target id
    VMin,          // v entry has the lowest target id
    LRMin(TargetInfo), // l and r entries have the same target id which is lower than the v entry. The target parameter is the maximum of l and r
    LVMin,      // l and v entries have the same target id which is lower than the r entry
    RVMin,      // r and v entries have the same target id which is lower than the l entry
    None,       // No entry is present
}

impl<'a, 'b, I: Iterator<Item = Target>> Iterator for TVAddVecMrgMax<'a, 'b, I> {
    type Item = TargetInfo;

    fn next(&mut self) -> Option<Self::Item> {
        let l = self.left.get(self.left_ix);
        let r = self.right.get(self.right_ix);
        let v = self.v.peek();

        let which = match (l, r, v) {
            (Some(l), Some(r), Some(v)) => match l.cmp_target_id(r) {
                Ordering::Equal => match l.target_id().cmp(&v.id()) {
                    Ordering::Equal => AVOrd::AllEqual(l.max2(r, v).unwrap()),
                    Ordering::Less => AVOrd::LRMin(l.max(r).unwrap()),
                    Ordering::Greater => AVOrd::VMin,
                },
                Ordering::Less => match l.target_id().cmp(&v.id()) {
                    Ordering::Equal => AVOrd::LVMin,
                    Ordering::Less => AVOrd::LMin,
                    Ordering::Greater => AVOrd::VMin,
                },
                Ordering::Greater => match r.target_id().cmp(&v.id()) {
                    Ordering::Equal => AVOrd::RVMin,
                    Ordering::Less => AVOrd::RMin,
                    Ordering::Greater => AVOrd::VMin,
                },
            },
            (Some(l), Some(r), None) => match l.cmp_target_id(r) {
                Ordering::Equal => AVOrd::LRMin(l.max(r).unwrap()),
                Ordering::Less => AVOrd::LMin,
                Ordering::Greater => AVOrd::RMin,
            },
            (Some(l), None, Some(v)) => match l.target_id().cmp(&v.id()) {
                Ordering::Equal => AVOrd::LVMin,
                Ordering::Less => AVOrd::LMin,
                Ordering::Greater => AVOrd::VMin,
            },
            (None, Some(r), Some(v)) => match r.target_id().cmp(&v.id()) {
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
            AVOrd::AllEqual(t) => {
                self.left_ix += 1;
                self.right_ix += 1;
                let _ = self.v.next();
                Some(t)
            }
            AVOrd::LRMin(t) => {
                self.left_ix += 1;
                self.right_ix += 1;
                Some(t)
            }
            AVOrd::LVMin => {
                self.left_ix += 1;
                let v = self.v.next();
                l.and_then(|left| left.add_target(&v.unwrap()))
            }
            AVOrd::RVMin => {
                self.right_ix += 1;
                let _ = self.v.next();
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
                let v = self.v.next();
                v.map(|t| t.into_target_info())
            }
            AVOrd::None => None,
        }
    }
}

