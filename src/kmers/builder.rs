use super::KmerType;

pub struct KmerBuilder {
    kmer: KmerType,
    valid: KmerType,
    mask: KmerType,
    valid_mask: KmerType,
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
