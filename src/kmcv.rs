use std::{collections::HashMap, io::BufRead};

use crate::kmers::{KmerIdxHits, KmerType};
use anyhow::Context;
use log::{log_enabled, Level::Trace};
use serde::Serialize;

fn get_u16_from_slice(p: &[u8]) -> u16 {
    u16::from_le_bytes(p.try_into().expect("Slice has wrong size"))
}

fn get_u32_from_slice(p: &[u8]) -> u32 {
    u32::from_le_bytes(p.try_into().expect("Slice has wrong size"))
}

fn get_u64_from_slice(p: &[u8]) -> u64 {
    u64::from_le_bytes(p.try_into().expect("Slice has wrong size"))
}

#[derive(Debug, Serialize)]
pub struct KmcvHeader {
    #[serde(skip_serializing)]
    version: [u8; 2],
    kmer_length: u8,
    #[serde(skip_serializing)]
    max_hits: u8,
    n_contigs: u32,
    n_targets: u32,
    rnd_id: u32,
    #[serde(skip_serializing)]
    mapped: u64,
    #[serde(skip_serializing)]
    on_target: u64,
    #[serde(skip_serializing)]
    redundant: u64,
    #[serde(skip_serializing)]
    total_hits: u64,
}

impl KmcvHeader {
    pub fn read<R: BufRead>(rdr: &mut R) -> anyhow::Result<Self> {
        let mut buf = [0; 52];
        rdr.read_exact(&mut buf)
            .with_context(|| "Error reading header from kmer file")?;

        if buf[0..4] != [b'K', b'M', b'C', b'V'] {
            return Err(anyhow!(
                "Incorrect magic number from header block of  kmer file"
            ));
        }
        let version = [buf[4], buf[5]];
        if version[0] != 2 {
            return Err(anyhow!("Incorrect version for kmer file (expected V2)"));
        }
        let kmer_length = buf[6];
        let max_hits = buf[7];
        if (kmer_length as u32) << 1 > KmerType::BITS {
            return Err(anyhow!("Kmer length {kmer_length} too large for KmerType"));
        }
        let rnd_id = get_u32_from_slice(&buf[8..12]);
        let n_contigs = get_u32_from_slice(&buf[12..16]);
        let n_targets = get_u32_from_slice(&buf[16..20]);
        let mapped = get_u64_from_slice(&buf[20..28]);
        let on_target = get_u64_from_slice(&buf[28..36]);
        let redundant = get_u64_from_slice(&buf[36..44]);
        let total_hits = get_u64_from_slice(&buf[44..]);

        if mapped < on_target || on_target < redundant {
            Err(anyhow!("Kmer totals are inconsistent"))
        } else {
            Ok(Self {
                version,
                kmer_length,
                max_hits,
                rnd_id,
                n_contigs,
                n_targets,
                mapped,
                on_target,
                redundant,
                total_hits,
            })
        }
    }

    pub fn read_close<R: BufRead>(&self, rdr: &mut R) -> anyhow::Result<()> {
        let mut buf: [u8; 8] = [0; 8];

        let e = || "Error reading closing block from kmer file";

        rdr.read_exact(&mut buf).with_context(e)?;

        // Check magic number
        if buf[4..8] != [b'V', b'C', b'M', b'K'] {
            return Err(anyhow!(
                "Incorrect magic number from closing block of kmer file"
            ));
        }

        // Check random ID
        let rnd_id = get_u32_from_slice(&buf[0..4]);
        if rnd_id != self.rnd_id {
            return Err(anyhow!("Incorrect run id from closing block of kmer file"));
        }

        // Check we are at EOF
        if rdr.read(&mut buf).with_context(e)? != 0 {
            Err(anyhow!("Trailing garbage at end of kmer file"))
        } else {
            Ok(())
        }
    }

    pub fn kmer_length(&self) -> u8 {
        self.kmer_length
    }

    pub fn rnd_id(&self) -> u32 {
        self.rnd_id
    }
}

pub struct Target {
    start: u32,
    end: u32,
}

impl Target {
    fn read<R: BufRead>(rdr: &mut R, n_contigs: u32) -> anyhow::Result<(Self, u32)> {
        let mut buf = [0u8; 12];
        rdr.read_exact(&mut buf)
            .with_context(|| "Error reading target block from kmer file")?;

        let contig = Self::get_contig(&buf[..4], n_contigs)?;

        let (start, end) = Self::get_start_end(&buf[4..])?;

        Ok((Self { start, end }, contig))
    }

    fn get_contig(buf: &[u8], n_contigs: u32) -> anyhow::Result<u32> {
        let contig = get_u32_from_slice(&buf[..4]);
        if contig >= n_contigs {
            Err(anyhow!(
                "Contig id {contig}  from target definition not in range"
            ))
        } else {
            Ok(contig)
        }
    }

    fn get_start_end(buf: &[u8]) -> anyhow::Result<(u32, u32)> {
        let start = get_u32_from_slice(&buf[..4]);
        let end = get_u32_from_slice(&buf[4..]);
        if end < start {
            Err(anyhow!("End coordinate of target less than start"))
        } else {
            Ok((start, end))
        }
    }

    pub fn size(&self) -> u32 {
        self.end + 1 - self.start
    }
}

pub struct KContig {
    name: Box<str>,
    targets: Vec<u32>,
}

impl KContig {
    fn read<R: BufRead>(rdr: &mut R) -> anyhow::Result<Self> {
        let mut buf = [0u8; 2];

        rdr.read_exact(&mut buf)
            .with_context(|| "Error reading string length from kmer file")?;

        let l = get_u16_from_slice(&buf) as usize;
        if l == 0 {
            return Err(anyhow!("Contig name length is zero"));
        }

        let mut s = String::with_capacity(l);
        while s.len() < l {
            let p = rdr
                .fill_buf()
                .with_context(|| "Error while reading contig name")?;
            let m = l - s.len();
            let n = m.min(p.len());
            let s1 = std::str::from_utf8(&p[..n]).with_context(|| "Contig name not utf8")?;
            s.push_str(s1);
            rdr.consume(n);
            if n == m {
                break;
            }
        }
        trace!("Read contig {s}");
        let name = s.into_boxed_str();
        let targets = Vec::new();
        Ok(Self { name, targets })
    }
}

#[derive(Debug)]
enum KType {
    UniqueOffTarget,
    Redundant,
    Mapped(u8),
}

impl KType {
    fn from_u8(x: u8) -> anyhow::Result<KType> {
        let x = x & 0xf;
        match x {
            0..=7 => Ok(KType::Mapped(x + 1)),
            8 => Ok(KType::UniqueOffTarget),
            9 => Ok(KType::Redundant),
            _ => Err(anyhow!("Unexpected kmer type")),
        }
    }
}

fn get_type_skip<R: BufRead>(rdr: &mut R) -> anyhow::Result<(KType, usize)> {
    let mut buf = [0u8; 4];

    rdr.read_exact(&mut buf[..1])
        .with_context(|| "Error reading type byte string from kmer block")?;
    let ktype = KType::from_u8(buf[0])?;
    let mut skip = (buf[0] >> 4) as usize;
    if skip == 0xf {
        rdr.read_exact(&mut buf[..1])
            .with_context(|| "Error reading skip byte string from kmer block")?;
        if buf[0] == 0xff {
            skip += 0xff;
            rdr.read_exact(&mut buf[..2])
                .with_context(|| "Error reading skip u16from kmer block")?;
            let s = get_u16_from_slice(&buf[..2]);
            if s == 0xffff {
                skip += 0xffff;
                rdr.read_exact(&mut buf)
                    .with_context(|| "Error reading skip u16from kmer block")?;
                let s = get_u32_from_slice(&buf[..4]);
                skip += s as usize;
            } else {
                skip += s as usize;
            }
        } else {
            skip += buf[0] as usize;
        }
    }
    Ok((ktype, skip))
}

#[derive(Serialize)]
pub struct Kmcv {
    #[serde(flatten)]
    header: KmcvHeader,
    #[serde(skip_serializing)]
    kmers: HashMap<KmerType, KmerIdxHits>,
    #[serde(skip_serializing)]
    target_hits: Vec<u32>,
    #[serde(skip_serializing)]
    contigs: Vec<KContig>,
    #[serde(skip_serializing)]
    targets: Vec<Target>,
}

impl Kmcv {
    pub fn read<R: BufRead>(rdr: &mut R) -> anyhow::Result<Self> {
        debug!("Reading header from kmer file");
        let header =
            KmcvHeader::read(rdr).with_context(|| "Error reading header from kmer file")?;

        let n_hits = header.total_hits as usize;
        let kmers_to_store = (header.on_target - header.redundant) as usize;
        let n_ctgs = header.n_contigs as usize;
        let n_targets = header.n_targets as usize;

        let mut kmcv = Self {
            header,
            kmers: HashMap::with_capacity(kmers_to_store),
            target_hits: Vec::with_capacity(n_hits + 1),
            contigs: Vec::with_capacity(n_ctgs),
            targets: Vec::with_capacity(n_targets),
        };

        debug!("Reading contig blocks from kmer file");
        kmcv.read_contig_blocks(rdr)
            .with_context(|| "Error reading contig information")?;

        debug!("Reading target blocks from kmer file");
        kmcv.read_target_blocks(rdr)
            .with_context(|| "Error reading target information")?;

        debug!("Reading kmer blocks from kmer file");
        kmcv.read_kmer_blocks(rdr)
            .with_context(|| "Error reading kmer information")?;

        debug!("Reading closing blocks from kmer file");
        kmcv.header
            .read_close(rdr)
            .with_context(|| "Error reading closing block from kmer file")?;

        Ok(kmcv)
    }

    pub fn n_targets(&self) -> usize {
        self.targets.len()
    }

    pub fn kmer_length(&self) -> u8 {
        self.header.kmer_length
    }

    pub fn target_hits(&self, kmer: KmerType) -> Option<&[u32]> {
        self.kmers.get(&kmer).and_then(|k| {
            let (n, i) = (k.hits() as usize, k.idx());
            if self.target_hits[i] == 0 {
                if n > 1 {
                    Some(&self.target_hits[i + 1..i + n])
                } else {
                    None
                }
            } else {
                Some(&self.target_hits[i..i + n])
            }
        })
    }

    pub fn targets(&self) -> &[Target] {
        &self.targets
    }

    /// Private functions
    fn read_contig_blocks<R: BufRead>(&mut self, rdr: &mut R) -> anyhow::Result<()> {
        self.contigs.clear();
        for _ in 0..self.header.n_contigs {
            self.contigs.push(KContig::read(rdr)?)
        }
        Ok(())
    }

    fn read_target_blocks<R: BufRead>(&mut self, rdr: &mut R) -> anyhow::Result<()> {
        let n_contigs = self.header.n_contigs;

        for ix in 0..self.header.n_targets {
            let (target, contig) = Target::read(rdr, n_contigs)?;
            self.contigs[contig as usize].targets.push(ix);
            self.targets.push(target);
        }

        if log_enabled!(Trace) {
            for ctg in self.contigs.iter() {
                debug!(
                    "Contig {} number of targets: {}",
                    ctg.name,
                    ctg.targets.len()
                )
            }
        }

        Ok(())
    }

    fn read_kmer_blocks<R: BufRead>(&mut self, rdr: &mut R) -> anyhow::Result<()> {
        let total_kmers = 1usize << (self.header.kmer_length << 1);
        let blocks_in_file = self.header.mapped;
        self.target_hits.clear();
        self.target_hits.push(0);

        let mut kmer: usize = 0;
        for _ in 0..blocks_in_file {
            kmer = self.read_kmer_block(rdr, kmer, total_kmers)?
        }
        if self.kmers.len() != (self.header.on_target - self.header.redundant) as usize {
            Err(anyhow!("Unexpected number of kmers stored"))
        } else if self.target_hits.len() - 1 != self.header.total_hits as usize {
            Err(anyhow!("Unexpected number of target hits read"))
        } else {
            debug!(
                "Kmers stored: {}, target hits stored: {}",
                self.kmers.len(),
                self.target_hits.len() - 1,
            );
            Ok(())
        }
    }

    fn read_kmer_block<R: BufRead>(
        &mut self,
        rdr: &mut R,
        kmer: usize,
        total_kmers: usize,
    ) -> anyhow::Result<usize> {
        let (ktype, skip) = get_type_skip(rdr)?;
        let kmer = kmer + skip;

        if kmer >= total_kmers {
            Err(anyhow!("kmer is larger than maximum"))
        } else if let KType::Mapped(hits) = ktype {
            let i = self.target_hits.len();
            if i + hits as usize > (self.header.total_hits + 1) as usize {
                Err(anyhow!("Too many target hits in kmer file"))
            } else if self
                .kmers
                .insert(kmer as KmerType, KmerIdxHits::new(i, hits))
                .is_some()
            {
                Err(anyhow!("Zero skip leading to duplicate kmers"))
            } else {
                for _ in 0..hits {
                    self.process_hits(rdr).with_context(|| {
                        format!("Error processing hits for kmer {} (skip = {})", kmer, skip)
                    })?
                }
                self.target_hits[i..].sort_unstable();
                Ok(kmer)
            }
        } else {
            Ok(kmer)
        }
    }

    fn process_hits<R: BufRead>(&mut self, rdr: &mut R) -> anyhow::Result<()> {
        let mut buf = [0u8; 4];
        rdr.read_exact(&mut buf)
            .with_context(|| "Error reading contig hit from kmer block")?;
        let ix = get_u32_from_slice(&buf);
        if ix > self.header.n_targets {
            Err(anyhow!(
                "Illegal target id read from kmer block (read {ix}, maximum expected: {}",
                self.header.n_targets
            ))
        } else {
            self.target_hits.push(ix);
            Ok(())
        }
    }
}
