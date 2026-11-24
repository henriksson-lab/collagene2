//! Strict one-block framing around a pure-Rust DEFLATE codec.
use crate::*;
use std::io::{Read, Write};
const EOF: [u8; 28] = [
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VirtualOffset {
    pub compressed: u64,
    pub uncompressed: u16,
}
impl VirtualOffset {
    pub fn packed(self) -> Result<u64> {
        if self.compressed >= 1 << 48 {
            return Err(Error::Limit("48-bit BGZF address"));
        }
        Ok(self.compressed << 16 | u64::from(self.uncompressed))
    }
}
pub struct BgzfReader<R> {
    inner: R,
    limits: Limits,
    buf: Vec<u8>,
    pos: usize,
    block_start: u64,
    compressed: u64,
    decoded: u64,
    blocks: u64,
    hash: blake3::Hasher,
    ended: bool,
    last_empty: bool,
    failed: bool,
}
impl<R: Read> BgzfReader<R> {
    pub fn new(inner: R, limits: Limits) -> Result<Self> {
        limits.validate()?;
        Ok(Self {
            inner,
            limits,
            buf: vec![],
            pos: 0,
            block_start: 0,
            compressed: 0,
            decoded: 0,
            blocks: 0,
            hash: blake3::Hasher::new(),
            ended: false,
            last_empty: false,
            failed: false,
        })
    }
    pub fn virtual_offset(&self) -> VirtualOffset {
        if self.pos == self.buf.len() {
            VirtualOffset {
                compressed: self.compressed,
                uncompressed: 0,
            }
        } else {
            VirtualOffset {
                compressed: self.block_start,
                uncompressed: self.pos as u16,
            }
        }
    }
    pub fn peek_virtual_offset(&mut self) -> io::Result<Option<VirtualOffset>> {
        if self.failed {
            return Err(io::Error::other(Error::Poisoned));
        }
        while self.pos == self.buf.len() {
            if self.ended {
                return Ok(None);
            }
            match self.block() {
                Ok(true) => {}
                Ok(false) => return Ok(None),
                Err(e) => {
                    self.failed = true;
                    return Err(io::Error::new(io::ErrorKind::InvalidData, e));
                }
            }
        }
        Ok(Some(self.virtual_offset()))
    }
    pub fn source_bytes(&self) -> u64 {
        self.compressed
    }
    pub fn source_blake3(&self) -> Option<String> {
        (self.ended && !self.failed).then(|| self.hash.finalize().to_hex().to_string())
    }
    fn bytes(&mut self, b: &mut [u8]) -> Result<()> {
        let next = checked_add(
            self.compressed,
            b.len() as u64,
            self.limits.max_input_bytes,
            "compressed stream bytes",
        )?;
        self.inner.read_exact(b)?;
        self.compressed = next;
        self.hash.update(b);
        Ok(())
    }
    fn block(&mut self) -> Result<bool> {
        let mut fixed = [0; 12];
        let mut first = [0];
        match self.inner.read(&mut first)? {
            0 => {
                if !self.last_empty {
                    return Err(Error::Invalid("missing final BGZF empty EOF block"));
                }
                self.ended = true;
                return Ok(false);
            }
            1 => fixed[0] = first[0],
            _ => unreachable!(),
        };
        self.block_start = self.compressed;
        self.compressed = checked_add(
            self.compressed,
            1,
            self.limits.max_input_bytes,
            "compressed stream bytes",
        )?;
        self.hash.update(&first);
        self.bytes(&mut fixed[1..])?;
        if fixed[..4] != [31, 139, 8, 4] {
            return Err(Error::Unsupported("requires BGZF gzip FEXTRA framing"));
        }
        let xlen = u16::from_le_bytes([fixed[10], fixed[11]]) as usize;
        if !(6..=65536 - 20).contains(&xlen) {
            return Err(Error::Invalid("BGZF extra header length"));
        }
        let mut extra = vec![0; xlen];
        self.bytes(&mut extra)?;
        let mut i = 0;
        let mut size = None;
        while i < extra.len() {
            if extra.len() - i < 4 {
                return Err(Error::Invalid("truncated BGZF extra field"));
            }
            let n = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
            let end = i
                .checked_add(4 + n)
                .filter(|n| *n <= extra.len())
                .ok_or(Error::Invalid("BGZF extra field length"))?;
            if &extra[i..i + 2] == b"BC" {
                if n != 2 || size.is_some() {
                    return Err(Error::Invalid("BGZF BC field"));
                }
                size = Some(u16::from_le_bytes([extra[i + 4], extra[i + 5]]) as usize + 1);
            }
            i = end;
        }
        let size = size.ok_or(Error::Invalid("missing BGZF BC field"))?;
        let head = 12 + xlen;
        if size < head + 10 || size > 65536 {
            return Err(Error::Invalid("BGZF block size"));
        }
        let mut rest = vec![0; size - head];
        self.bytes(&mut rest)?;
        let end = rest.len() - 8;
        let crc = u32::from_le_bytes(rest[end..end + 4].try_into().unwrap());
        let len = u32::from_le_bytes(rest[end + 4..].try_into().unwrap()) as usize;
        if len > 65536 {
            return Err(Error::Invalid("BGZF uncompressed block exceeds 64 KiB"));
        }
        self.blocks = checked_add(self.blocks, 1, self.limits.max_bgzf_blocks, "BGZF blocks")?;
        self.decoded = checked_add(
            self.decoded,
            len as u64,
            self.limits.max_uncompressed_bytes,
            "uncompressed stream bytes",
        )?;
        let mut output = vec![0; len + 1];
        let mut decoder = flate2::Decompress::new(false);
        let status = decoder
            .decompress(&rest[..end], &mut output, flate2::FlushDecompress::Finish)
            .map_err(|_| Error::Invalid("DEFLATE stream"))?;
        if status != flate2::Status::StreamEnd
            || decoder.total_in() != end as u64
            || decoder.total_out() != len as u64
        {
            return Err(Error::Invalid("BGZF DEFLATE size/end mismatch"));
        }
        output.truncate(len);
        if crc32fast::hash(&output) != crc {
            return Err(Error::Invalid("BGZF CRC mismatch"));
        }
        self.last_empty = output.is_empty();
        self.buf = output;
        self.pos = 0;
        Ok(true)
    }
}
impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, dst: &mut [u8]) -> io::Result<usize> {
        if dst.is_empty() {
            return Ok(0);
        }
        if self.failed {
            return Err(io::Error::other(Error::Poisoned));
        }
        while self.pos == self.buf.len() {
            if self.ended {
                return Ok(0);
            }
            match self.block() {
                Ok(true) => {}
                Ok(false) => return Ok(0),
                Err(e) => {
                    self.failed = true;
                    return Err(io::Error::new(io::ErrorKind::InvalidData, e));
                }
            }
        }
        let n = dst.len().min(self.buf.len() - self.pos);
        dst[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}
pub struct BgzfWriter<W: Write> {
    inner: W,
    buf: Vec<u8>,
    written: u64,
    finished: bool,
    failed: bool,
    output_limit: u64,
}
impl<W: Write> BgzfWriter<W> {
    pub fn new(inner: W) -> Self {
        Self::with_limit(inner, 64 << 30)
    }
    pub fn with_limit(inner: W, output_limit: u64) -> Self {
        Self {
            inner,
            buf: Vec::with_capacity(32768),
            written: 0,
            finished: false,
            failed: false,
            output_limit,
        }
    }
    pub fn virtual_offset(&self) -> VirtualOffset {
        VirtualOffset {
            compressed: self.written,
            uncompressed: self.buf.len() as u16,
        }
    }
    fn flush_block(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let mut e = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(&self.buf)?;
        let encoded = e.finish()?;
        let size = 18 + encoded.len() + 8;
        if size > 65536 {
            return Err(io::Error::other("BGZF encoded block too large"));
        }
        let next = self
            .written
            .checked_add(size as u64)
            .filter(|n| *n <= self.output_limit)
            .ok_or_else(|| io::Error::other("BGZF output bytes limit"))?;
        let mut h = EOF[..18].to_vec();
        h[16..18].copy_from_slice(&((size - 1) as u16).to_le_bytes());
        self.inner.write_all(&h)?;
        self.inner.write_all(&encoded)?;
        self.inner
            .write_all(&crc32fast::hash(&self.buf).to_le_bytes())?;
        self.inner
            .write_all(&(self.buf.len() as u32).to_le_bytes())?;
        self.written = next;
        self.buf.clear();
        Ok(())
    }
    pub fn finish(mut self) -> io::Result<W> {
        if self.failed {
            return Err(io::Error::other(Error::Poisoned));
        }
        self.flush_block()?;
        if self
            .written
            .checked_add(EOF.len() as u64)
            .is_none_or(|n| n > self.output_limit)
        {
            return Err(io::Error::other("BGZF output bytes limit"));
        }
        self.inner.write_all(&EOF)?;
        self.inner.flush()?;
        self.finished = true;
        Ok(self.inner)
    }
}
impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        if self.failed {
            return Err(io::Error::other(Error::Poisoned));
        }
        if self.finished {
            return Err(io::Error::other("writer finished"));
        }
        let mut rest = b;
        while !rest.is_empty() {
            let n = (32768 - self.buf.len()).min(rest.len());
            self.buf.extend_from_slice(&rest[..n]);
            rest = &rest[n..];
            if self.buf.len() == 32768
                && let Err(e) = self.flush_block()
            {
                self.failed = true;
                return Err(e);
            }
        }
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        if self.failed {
            return Err(io::Error::other(Error::Poisoned));
        }
        let r = self.flush_block().and_then(|_| self.inner.flush());
        if r.is_err() {
            self.failed = true;
        }
        r
    }
}
