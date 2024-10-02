use std::cmp::Ordering;

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub(super) struct Target(u32);

impl Target {
    #[inline]
    pub(super) fn new(id: u32) -> Self {
        Self(id)
    }

    #[inline]
    pub(super) fn id(&self) -> u32 {
        self.0
    }

    #[inline]
    pub(super) fn into_target_info(self) -> TargetInfo {
        TargetInfo {
            target: self,
            count: 1,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub(super) struct TargetInfo {
    target: Target,
    count: u32,
}

impl TargetInfo {
    #[inline]
    pub(super) fn target_id(&self) -> u32 {
        self.target.id()
    }

    #[inline]
    pub(super) fn target(&self) -> &Target {
        &self.target
    }

    #[inline]
    pub(super) fn cmp_target_id(&self, other: &Self) -> Ordering {
        self.target.id().cmp(&other.target.id())
    }

    #[inline]
    pub(super) fn max(&self, other: &Self) -> Option<Self> {
        if self.target.id() == other.target.id() {
            Some(match self.count.cmp(&other.count) {
                Ordering::Greater | Ordering::Equal => *self,
                Ordering::Less => *other,
            })
        } else {
            None
        }
    }

    #[inline]
    pub(super) fn max2(&self, other: &Self, t: &Target) -> Option<Self> {
        if self.target.id() == other.target.id() && self.target.id() == t.id() {
            match (self.count + 1).cmp(&other.count) {
                Ordering::Greater => self.add_target(t),
                Ordering::Less => Some(*other),
                Ordering::Equal => self.add_target(t).and_then(|n| n.max(other)),
            }
        } else {
            None
        }
    }

    #[inline]
    pub(super) fn add_target(&self, t: &Target) -> Option<Self> {
        if self.target.id() == t.id() {
            let count = self.count + 1;
            let target = self.target;
            Some(Self { target, count })
        } else {
            None
        }
    }

    #[inline]
    pub(super) fn is_mapped(&self) -> bool {
        self.target_id() > 0 && self.count > 2
    }
    #[inline]
    pub(super) fn count(&self) -> u32 {
        self.count
    }
}

#[derive(Default, Debug)]
pub(super) struct TargetVec {
    v: Vec<TargetInfo>, // Keep track of which targets have been tagged for a read
}

impl TargetVec {
    #[inline]
    pub(super) fn clear(&mut self) {
        self.v.clear();
    }
    #[inline]
    pub(super) fn as_slice(&self) -> &[TargetInfo] {
        &self.v
    }
    #[inline]
    pub(super) fn push(&mut self, t: TargetInfo) {
        self.v.push(t);
    }
    #[inline]
    pub(super) fn extend(&mut self, v: &[TargetInfo]) {
        self.v.extend_from_slice(v);
    }
    #[inline]
    pub(super) fn is_empty(&self) -> bool {
        self.v.is_empty()
    }
}

mod test {
    #[allow(unused_imports)]
    use super::*;
    #[test]
    fn target_info_basic() {
        let t = Target::new(20);
        let ti = t.into_target_info();

        assert_eq!(ti.target_id(), 20);
    }
    #[test]
    fn target_info_add_max() {
        let t = Target::new(20);
        let ti = t.into_target_info();
        let t1 = Target::new(20);
        let ti1 = ti.add_target(&t1).unwrap();
        assert_eq!(ti1.count(), 2);
        let ti2 = ti.max(&ti1);
        assert_eq!(ti2, Some(ti1));
        let t2 = Target::new(3);
        assert_eq!(ti.add_target(&t2), None);
    }
}
