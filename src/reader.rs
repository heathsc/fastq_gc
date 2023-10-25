use std::io::Read;

use anyhow::Context;
use compress_io::compress::CompressIo;
use crossbeam_channel::{Receiver, Sender};

use crate::cli::Config;

const BUF_SIZE: usize = 8192;

pub struct Buffer {
    inner: Box<[u8]>,
    used: usize,
}

impl Default for Buffer {
    fn default() -> Self {
        let inner = vec![0; BUF_SIZE].into_boxed_slice();
        Self { inner, used: 0 }
    }
}

impl Buffer {
    /// Fills buffer unless EOF or an error occurs
    /// If buffer is full, returns a slice containing the possibly incomplete last entry,
    /// otherwise it returns an empty slice
    fn fill<R: Read>(&mut self, rdr: &mut R, work: &mut Vec<(usize, u8)>) -> anyhow::Result<&[u8]> {
        // read() will not always fill the buffer, even if we are not at EOF
        // so we must loop until EOF is reached or the buffer is full
        loop {
            let l = rdr
                .read(&mut self.inner[self.used..])
                .with_context(|| "Error reading from input")?;

            self.used += l;

            // If buffer is full or we are at EOF
            if self.used == self.inner.len() || l == 0 {
                break;
            }
        }

        // Return if at EOF, otherwise search backwards for beginning of last (possibly incomplete)
        // record and set used appropriately
        if self.used < self.inner.len() {
            // EOF
            Ok(&[])
        } else {
            let i = find_last_fastq_entry(&self.inner, work)
                .with_context(|| "Could not detect last FASTQ record")?;
            if i == 0 {
                Err(anyhow!("FASTQ record larger than buffer"))
            } else {
                self.used = i - 1;
                Ok(&self.inner[i..])
            }
        }
    }

    fn push_slice(&mut self, s: &[u8]) {
        let l = s.len();
        self.inner[self.used..self.used + l].copy_from_slice(s);
        self.used += s.len()
    }

    pub fn as_slice(&self) -> &[u8] {
        &self.inner[..self.used]
    }

    pub fn clear(&mut self) {
        self.used = 0
    }

    pub fn is_empty(&self) -> bool {
        self.used == 0
    }

    pub fn fastq(&self) -> FastQIter {
        FastQIter {
            inner: self.as_slice(),
        }
    }
}

pub struct Lines<'a> {
    inner: &'a [u8],
}

impl<'a> Iterator for Lines<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.inner.is_empty() {
            None
        } else if let Some((i, _)) = self.inner.iter().enumerate().find(|(_, c)| **c == b'\n') {
            let (s1, s2) = self.inner.split_at(i);
            self.inner = &s2[1..];
            Some(s1)
        } else {
            let s = self.inner;
            self.inner = &[];
            Some(s)
        }
    }
}

pub struct FastQRecord<'a> {
    id: &'a [u8],
    seq: &'a [u8],
    qual: &'a [u8],
}

impl<'a> FastQRecord<'a> {
    #[allow(dead_code)]
    pub fn id(&self) -> &[u8] {
        self.id
    }

    pub fn seq(&self) -> &[u8] {
        self.seq
    }

    pub fn qual(&self) -> &[u8] {
        self.qual
    }
}

pub struct FastQIter<'a> {
    inner: &'a [u8],
}

impl<'a> Iterator for FastQIter<'a> {
    type Item = anyhow::Result<FastQRecord<'a>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.inner.is_empty() {
            None
        } else {
            Some({
                let mut itr = Lines { inner: self.inner };
                let raw_id = itr.next();
                let seq = itr.next();
                let id2 = itr.next();
                let qual = itr.next();
                if let (Some(raw_id), Some(seq), Some(id2), Some(qual)) = (raw_id, seq, id2, qual) {
                    if raw_id.is_empty()
                        || raw_id[0] != b'@'
                        || id2.is_empty()
                        || id2[0] != b'+'
                        || seq.is_empty()
                        || seq.len() != qual.len()
                    {
                        Err(anyhow!("Invalid FASTQ record"))
                    } else {
                        let id = &raw_id[1..];
                        self.inner = itr.inner;
                        Ok(FastQRecord { id, seq, qual })
                    }
                } else {
                    Err(anyhow!("Incomplete FASTQ record"))
                }
            })
        }
    }
}

/// bit 1 set if character compatible with FASTQ line 1 (Identifier)
/// bit 2 set if character compatible with FASTQ line 2 (Sequence)
/// bit 3 set if character compatible with FASTQ line 3 (Secondary identifier)
/// bit 4 set if character compatible with FASTQ line 4 (Quality)
const MASK_TAB: [u8; 256] = [
    0, 0, 0, 0, 0, 0, 0, 0, 5, 0xf, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 5, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xf, 0xd, 0xd, 0xd, 0xd,
    0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
    0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
    0xf, 0xf, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
    0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xd, 0xd, 0xd, 0xf,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,
];

/// Returns index of start of last FASTQ entry in slice, if present
fn find_last_fastq_entry(s: &[u8], lines: &mut Vec<(usize, u8)>) -> anyhow::Result<usize> {
    lines.clear();
    let mut curr_poss = 0xf; // Initially we could be at any of the first lines
    let mut prev_poss = 0xf;
    let mut prev_c: Option<u8> = None;
    let mut i = s.len();
    while i > 0 {
        i -= 1;
        let c = *unsafe { s.get_unchecked(i) };
        if c == b'\n' {
            if let Some(c1) = prev_c.take() {
                match c1 {
                    b'@' => curr_poss &= 9,   // Either ID or Qual line
                    b'+' => curr_poss &= 0xc, // Either ID2 or Qual line
                    _ => curr_poss &= 0xa,    // Either Seq or Qual line
                }
            }
            // If possibility mask has changed, propagate info through previous lines
            let mut found = if curr_poss == 1 { Some(i + 1) } else { None };
            if curr_poss != prev_poss {
                let mut cp = curr_poss;
                for (j, st) in lines.iter_mut().rev() {
                    let st1 = *st & ((cp >> 3) | (cp << 1));
                    if st1 == *st {
                        break;
                    }
                    cp = st1;
                    if cp == 1 {
                        found = Some(*j)
                    }
                    *st = st1;
                }
            }
            if let Some(j) = found.take() {
                return Ok(j);
            }
            // Store current line position and possibility mask
            lines.push((i + 1, curr_poss));
            // Set up mask for next line (previous line in file as we are going backwards)
            prev_poss = curr_poss;
            curr_poss = ((curr_poss & 1) << 3) | (curr_poss >> 1);
        } else {
            curr_poss &= MASK_TAB[c as usize];
            if curr_poss == 0 {
                trace!("No possible states remain");
                return Err(anyhow!("Invalid FASTQ format"));
            }
            prev_c = Some(c)
        }
    }
    Err(anyhow!("Start of previous FASTQ record not found"))
}

fn get_buffer(buf_list: &mut Vec<Buffer>, rcv: &Receiver<Buffer>) -> anyhow::Result<Buffer> {
    // Collect any buffers waiting in the channel
    while let Ok(b) = rcv.try_recv() {
        buf_list.push(b)
    }

    // If buf_list not empty, pop buffer from list and return
    if let Some(b) = buf_list.pop() {
        Ok(b)
    } else {
        // Otherwise, wait until a buffer becomes available
        rcv.recv()
            .with_context(|| "Error while waiting for empty buffer")
    }
    .map(|mut b| {
        b.clear();
        b
    })
}

pub fn reader(cfg: &Config, rcv: Receiver<Buffer>, snd: Sender<Buffer>) -> anyhow::Result<()> {
    // Try to open input file/stream
    let mut rdr = CompressIo::new()
        .opt_path(cfg.input_file())
        .reader()
        .with_context(|| "Could not open input file")?;

    debug!("Opened input");
    // Create empty buffers
    let nbuffers = cfg.threads() * 2;
    let mut buffer_list = Vec::with_capacity(nbuffers);
    for _ in 0..nbuffers {
        buffer_list.push(Buffer::default())
    }

    let mut work = Vec::new();
    let mut pending = buffer_list.pop().unwrap();

    debug!("Starting main read loop");
    // Main loop - read file until empty
    loop {
        let mut b = pending;
        let rem = b
            .fill(&mut rdr, &mut work)
            .with_context(|| "Error reading FASTQ data")?;
        pending = get_buffer(&mut buffer_list, &rcv)?;
        trace!("Read buffer");
        pending.push_slice(rem);
        if b.is_empty() {
            break;
        }
        trace!("Sending full buffer for processing");
        snd.send(b).with_context(|| "Error sending full buffer")?;
    }

    debug!("Finished reading input");
    Ok(())
}

mod test {
    use super::*;
    #[allow(dead_code)]
    fn test_split(s: &str, ix: usize, wlen: usize) {
        let mut work = Vec::new();
        println!("{:?}", s);
        let i = find_last_fastq_entry(s.as_bytes(), &mut work).unwrap();
        for (i, m) in work.iter() {
            println!("{} {}", i, m);
        }
        println!("{i} {:?}", s.split_at(i));
        assert_eq!(i, ix, "Bad split");
        assert_eq!(work.len(), wlen, "Wrong size for work vec");
    }

    #[test]
    fn test_find_last_fastq_entry() {
        let s =
            "@test\nACGGC\n+\nF&D*SD\n@test 1\nACGGC\n+e+@ws\n@F+&D@*SD\n@test 2\nACGGC\n+\nF&D*SD\n";

        let (s1, _) = s.split_at(55);
        test_split(s1, 52, 1);

        let (s1, _) = s.split_at(46);
        test_split(s1, 21, 3);

        let (s1, _) = s.split_at(52);
        test_split(s1, 52, 3);
    }

    #[allow(dead_code)]
    fn get_rec<'a, I: Iterator<Item = anyhow::Result<FastQRecord<'a>>>>(
        itr: &mut I,
    ) -> Option<FastQRecord<'a>> {
        match itr.next() {
            Some(Err(e)) => panic!("FASTQ error: {e}"),
            Some(Ok(r)) => Some(r),
            None => None,
        }
    }

    #[test]
    fn test_fastq_iter() {
        let s =
            "@test\nACGGC\n+\nF&D*S\n@test 1\nACGGCGCTT\n+e+@ws\n@F+&D@*SD\n@test 2\nACGGCC\n+\nF&D*SD\n";
        let bytes = s.as_bytes();
        let used = bytes.len();
        let inner = bytes.to_vec().into_boxed_slice();
        let buf = Buffer { inner, used };

        let mut itr = buf.fastq().skip(2);
        let rec = get_rec(&mut itr).unwrap();
        println!("{:?}\n{:?}\n{:?}", rec.id(), rec.seq(), rec.qual());
        assert_eq!(rec.qual(), &[70, 38, 68, 42, 83, 68]);
    }
}
