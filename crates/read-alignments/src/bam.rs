use crate::*;
use std::io::{Read, Write};
fn i32read<R: Read>(r: &mut R) -> Result<i32> {
    let mut b = [0; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}
fn sized(n: i32, cap: usize, label: &'static str) -> Result<usize> {
    usize::try_from(n)
        .ok()
        .filter(|v| *v <= cap)
        .ok_or(Error::Limit(label))
}
struct Bytes<'a> {
    b: &'a [u8],
    p: usize,
}
impl<'a> Bytes<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .p
            .checked_add(n)
            .filter(|n| *n <= self.b.len())
            .ok_or(Error::Invalid("truncated BAM record"))?;
        let v = &self.b[self.p..end];
        self.p = end;
        Ok(v)
    }
    fn byte(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn i32(&mut self) -> Result<i32> {
        Ok(self.u32()? as i32)
    }
    fn z(&mut self) -> Result<Vec<u8>> {
        let n = self.b[self.p..]
            .iter()
            .position(|b| *b == 0)
            .ok_or(Error::Invalid("unterminated BAM string"))?;
        let v = self.take(n)?.to_vec();
        self.take(1)?;
        Ok(v)
    }
}
fn string(b: &[u8]) -> Result<String> {
    std::str::from_utf8(b)
        .map(str::to_owned)
        .map_err(|_| Error::Invalid("UTF-8 string"))
}
fn width(t: u8) -> Result<usize> {
    match t {
        b'c' | b'C' => Ok(1),
        b's' | b'S' => Ok(2),
        b'i' | b'I' | b'f' => Ok(4),
        _ => Err(Error::Unsupported("BAM auxiliary type")),
    }
}
fn bits(c: &mut Bytes<'_>, t: u8) -> Result<u32> {
    Ok(match width(t)? {
        1 => u32::from(c.byte()?),
        2 => u32::from(c.u16()?),
        _ => c.u32()?,
    })
}
fn integer(t: u8, v: u32) -> i64 {
    match t {
        b'c' => i64::from(v as u8 as i8),
        b'C' => i64::from(v as u8),
        b's' => i64::from(v as u16 as i16),
        b'S' => i64::from(v as u16),
        b'i' => i64::from(v as i32),
        _ => i64::from(v),
    }
}
fn aux(c: &mut Bytes<'_>, l: &Limits) -> Result<Vec<Aux>> {
    let mut out = Vec::new();
    while c.p < c.b.len() {
        if out.len() >= l.max_aux_tags {
            return Err(Error::Limit("auxiliary tags"));
        }
        let tag = [c.byte()?, c.byte()?];
        let t = c.byte()?;
        let value = match t {
            b'A' => {
                let v = c.byte()?;
                if !(33..=126).contains(&v) {
                    return Err(Error::Invalid("auxiliary character"));
                }
                AuxValue::Char(v)
            }
            b'c' | b'C' | b's' | b'S' | b'i' | b'I' => AuxValue::Integer(integer(t, bits(c, t)?)),
            b'f' => AuxValue::FloatBits(c.u32()?),
            b'Z' | b'H' => {
                let v = c.z()?;
                if t == b'H' {
                    if v.len() % 2 != 0
                        || !v
                            .iter()
                            .all(|b| b.is_ascii_digit() || (b'A'..=b'F').contains(b))
                    {
                        return Err(Error::Invalid("BAM hex value"));
                    }
                    AuxValue::Hex(v)
                } else {
                    string(&v)?;
                    if v.iter().any(|b| *b < 32 || *b == 127) {
                        return Err(Error::Invalid("BAM auxiliary string"));
                    }
                    AuxValue::String(v)
                }
            }
            b'B' => {
                let t = c.byte()?;
                let n = sized(c.i32()?, l.max_array_values, "BAM array count")?;
                let w = width(t)?;
                if n.checked_mul(w).is_none_or(|n| n > c.b.len() - c.p) {
                    return Err(Error::Invalid("truncated BAM array"));
                }
                let mut values = Vec::with_capacity(n);
                for _ in 0..n {
                    values.push(bits(c, t)?);
                }
                AuxValue::Array {
                    element_type: t,
                    values,
                }
            }
            _ => return Err(Error::Unsupported("BAM auxiliary extension")),
        };
        out.push(Aux {
            tag,
            wire_type: t,
            value,
        });
    }
    Ok(out)
}
fn cigar(words: impl Iterator<Item = u32>, l: &Limits) -> Result<Vec<CigarOp>> {
    let mut out = Vec::new();
    for v in words {
        if out.len() >= l.max_cigar_ops {
            return Err(Error::Limit("CIGAR operators"));
        }
        out.push(CigarOp {
            len: v >> 4,
            kind: CigarKind::from_code(v & 15)?,
        });
    }
    Ok(out)
}
fn ref_name(id: i32, h: &Header) -> Result<Option<String>> {
    if id == -1 {
        Ok(None)
    } else {
        usize::try_from(id)
            .ok()
            .and_then(|i| h.references.get(i))
            .map(|r| Some(r.name.clone()))
            .ok_or(Error::Invalid("BAM reference ID"))
    }
}
fn pos(n: i32) -> Result<Option<u32>> {
    if n == -1 {
        Ok(None)
    } else {
        u32::try_from(n)
            .map(Some)
            .map_err(|_| Error::Invalid("BAM negative position"))
    }
}
fn parse(raw: Vec<u8>, h: &Header, l: &Limits, offset: VirtualOffset) -> Result<Record> {
    let mut c = Bytes { b: &raw[4..], p: 0 };
    let reference = ref_name(c.i32()?, h)?;
    let start = pos(c.i32()?)?;
    let name_len = usize::from(c.byte()?);
    let mapq = c.byte()?;
    let _bin = c.u16()?;
    let ncigar = usize::from(c.u16()?);
    let flags = c.u16()?;
    let len = sized(c.i32()?, l.max_sequence_bases, "BAM sequence length")?;
    let mate_reference = ref_name(c.i32()?, h)?;
    let mate_start = pos(c.i32()?)?;
    let template_length = c.i32()?;
    if name_len == 0 {
        return Err(Error::Invalid("BAM read name length"));
    }
    let name = c.take(name_len)?;
    if name.last() != Some(&0) || name[..name.len() - 1].contains(&0) {
        return Err(Error::Invalid("BAM read name terminator"));
    }
    let name = if name == b"*\0" {
        None
    } else {
        let n = string(&name[..name.len() - 1])?;
        if n.is_empty() || n.len() > 254 || !n.bytes().all(|b| (33..=126).contains(&b) && b != b'@')
        {
            return Err(Error::Invalid("BAM read name"));
        }
        Some(n)
    };
    if ncigar > l.max_cigar_ops {
        return Err(Error::Limit("CIGAR operators"));
    }
    let words = c.take(ncigar * 4)?;
    let mut ops = cigar(
        words
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| u32::from_le_bytes(*b)),
        l,
    )?;
    let packed = c.take(len.div_ceil(2))?;
    let mut sequence = String::with_capacity(len);
    for i in 0..len {
        let v = if i % 2 == 0 {
            packed[i / 2] >> 4
        } else {
            packed[i / 2] & 15
        };
        sequence.push(b"=ACMGRSVTWYHKDBN"[v as usize] as char);
    }
    let q = c.take(len)?;
    let qualities = if q.iter().all(|b| *b == 255) {
        None
    } else {
        if q.iter().any(|b| *b > 93) {
            return Err(Error::Unsupported(
                "mixed missing or out-of-range BAM quality scores",
            ));
        }
        Some(q.to_vec())
    };
    let aux = aux(&mut c, l)?;
    let mut long_cigar = false;
    if let Some(a) = aux.iter().find(|a| a.tag == *b"CG") {
        let AuxValue::Array {
            element_type: b'I',
            values,
        } = &a.value
        else {
            return Err(Error::Invalid("long CIGAR CG type"));
        };
        if values.len() <= 65535
            || ops.len() != 2
            || ops[0].kind != CigarKind::SoftClip
            || ops[0].len as usize != len
            || ops[1].kind != CigarKind::Skip
        {
            return Err(Error::Unsupported("noncanonical long CIGAR transport"));
        }
        let full = cigar(values.iter().copied(), l)?;
        let ref_len = full
            .iter()
            .filter(|op| op.kind.consumes_reference())
            .try_fold(0u64, |n, op| n.checked_add(u64::from(op.len)))
            .ok_or(Error::Invalid("long CIGAR span overflow"))?;
        if ref_len != u64::from(ops[1].len) {
            return Err(Error::Invalid("long CIGAR placeholder span"));
        }
        ops = full;
        long_cigar = true;
    }
    let fields = Fields {
        name,
        flags,
        reference,
        start,
        mapping_quality: (mapq != 255).then_some(mapq),
        cigar: ops,
        mate_reference,
        mate_start,
        template_length,
        sequence: (len > 0).then_some(sequence),
        qualities: if len == 0 { None } else { qualities },
        aux,
    };
    validate_fields(&fields, l)?;
    Ok(Record {
        fields,
        format: Format::Bam,
        original: raw,
        dictionary: h.dictionary_identity,
        virtual_offset: Some(offset),
        long_cigar,
    })
}
pub struct BamReader<R> {
    inner: BgzfReader<R>,
    header: Header,
    limits: Limits,
    records: u64,
    emitted: u64,
    failed: bool,
    complete: bool,
}
impl<R: Read> BamReader<R> {
    pub fn new(inner: R, limits: Limits) -> Result<Self> {
        let mut inner = BgzfReader::new(inner, limits.clone())?;
        let mut magic = [0; 4];
        inner.read_exact(&mut magic)?;
        if &magic != b"BAM\x01" {
            detect(&magic)?;
            return Err(Error::Invalid("BAM magic"));
        }
        let n = i32read(&mut inner)?;
        let len = sized(n, limits.max_header_bytes, "BAM header length")?;
        let mut raw = Vec::from(magic);
        raw.extend(n.to_le_bytes());
        let mut text = vec![0; len];
        inner.read_exact(&mut text)?;
        raw.extend(&text);
        let meaningful = text
            .iter()
            .rposition(|b| *b != 0)
            .map(|i| i + 1)
            .unwrap_or(0);
        let mut header = Header::from_sam(&text[..meaningful], &limits)?;
        let n = i32read(&mut inner)?;
        let count = sized(n, limits.max_references, "BAM reference count")?;
        raw.extend(n.to_le_bytes());
        let mut refs = Vec::with_capacity(count);
        let mut names = std::collections::BTreeSet::new();
        for _ in 0..count {
            let n = i32read(&mut inner)?;
            let len = sized(n, limits.max_header_bytes, "BAM reference name length")?;
            if len < 2
                || raw
                    .len()
                    .checked_add(8 + len)
                    .is_none_or(|n| n > limits.max_header_bytes)
            {
                return Err(Error::Limit("BAM aggregate header bytes"));
            }
            raw.extend(n.to_le_bytes());
            let mut name = vec![0; len];
            inner.read_exact(&mut name)?;
            raw.extend(&name);
            if name.last() != Some(&0) || name[..len - 1].contains(&0) {
                return Err(Error::Invalid("BAM reference name"));
            }
            let name = string(&name[..len - 1])?;
            let length = i32read(&mut inner)?;
            raw.extend(length.to_le_bytes());
            if length <= 0 || !sam::reference_name(&name) || !names.insert(name.clone()) {
                return Err(Error::Invalid("BAM reference dictionary"));
            }
            let circular = header
                .references
                .get(refs.len())
                .is_some_and(|r| r.circular);
            refs.push(Reference {
                name,
                length: length as u32,
                circular,
            });
        }
        if !header.references.is_empty() && header.references.as_ref() != refs.as_slice() {
            return Err(Error::Invalid("BAM text/binary dictionary mismatch"));
        }
        header.sam_dictionary_matches = refs.is_empty() || !header.references.is_empty();
        header.set_references(refs);
        header.bam_original = Some(raw);
        Ok(Self {
            inner,
            header,
            limits,
            records: 0,
            emitted: 0,
            failed: false,
            complete: false,
        })
    }
    pub fn header(&self) -> &Header {
        &self.header
    }
    pub fn progress(&self) -> StreamProgress {
        StreamProgress {
            state: if self.failed {
                StreamState::Failed
            } else if self.complete {
                StreamState::Complete
            } else {
                StreamState::Reading
            },
            records_emitted: self.emitted,
            source_bytes_accounted: self.inner.source_bytes(),
        }
    }
    pub fn completion(&self) -> Option<Completion> {
        if !self.complete || self.failed {
            return None;
        }
        Some(Completion{format:Format::Bam,records:self.records,source_bytes:self.inner.source_bytes(),source_blake3:self.inner.source_blake3()?,warnings:vec!["Imported mappings/auxiliary evidence are unverified; no coverage, reference match or variant validation was performed.".into()]})
    }
    pub fn next_record(&mut self) -> Result<Option<Record>> {
        if self.failed {
            return Err(Error::Poisoned);
        }
        let r = self.next();
        if r.is_err() {
            self.failed = true;
        }
        if matches!(&r, Ok(Some(_))) {
            self.emitted += 1;
        }
        r
    }
    fn next(&mut self) -> Result<Option<Record>> {
        if self.complete {
            return Ok(None);
        }
        let Some(offset) = self.inner.peek_virtual_offset()? else {
            self.complete = true;
            return Ok(None);
        };
        let mut b = [0; 4];
        if self.inner.read(&mut b[..1])? == 0 {
            self.complete = true;
            return Ok(None);
        }
        self.inner.read_exact(&mut b[1..])?;
        let len = sized(
            i32::from_le_bytes(b),
            self.limits.max_record_bytes.saturating_sub(4),
            "BAM record block size",
        )?;
        if len < 32 {
            return Err(Error::Invalid("BAM record smaller than fixed fields"));
        }
        self.records = checked_add(self.records, 1, self.limits.max_records, "records")?;
        let mut raw = Vec::with_capacity(4 + len);
        raw.extend(b);
        raw.resize(4 + len, 0);
        self.inner.read_exact(&mut raw[4..])?;
        Ok(Some(parse(raw, &self.header, &self.limits, offset)?))
    }
}
fn put_i32(b: &mut Vec<u8>, n: i32) {
    b.extend(n.to_le_bytes());
}
fn reference_id(n: Option<&str>, h: &std::collections::BTreeMap<String, usize>) -> Result<i32> {
    match n {
        None => Ok(-1),
        Some(n) => h.get(n).map(|i| *i as i32).ok_or(Error::Unsupported(
            "BAM export requires an explicit reference dictionary entry",
        )),
    }
}
fn write_aux(b: &mut Vec<u8>, a: &Aux) -> Result<()> {
    b.extend(a.tag);
    match &a.value {
        AuxValue::Char(v) => {
            b.extend([b'A', *v]);
        }
        AuxValue::Integer(n) => {
            let t = if a.wire_type == b'i' && *n > i32::MAX as i64 {
                b'I'
            } else {
                a.wire_type
            };
            let w = width(t)?;
            b.push(t);
            b.extend(&(*n as u32).to_le_bytes()[..w]);
        }
        AuxValue::FloatBits(v) => {
            b.push(b'f');
            b.extend(v.to_le_bytes());
        }
        AuxValue::String(v) | AuxValue::Hex(v) => {
            b.push(a.wire_type);
            b.extend(v);
            b.push(0);
        }
        AuxValue::Array {
            element_type,
            values,
        } => {
            b.extend([b'B', *element_type]);
            put_i32(b, values.len() as i32);
            let w = width(*element_type)?;
            for v in values {
                b.extend(&v.to_le_bytes()[..w]);
            }
        }
    }
    Ok(())
}
fn encoded_size(r: &Record, l: &Limits, cap: usize) -> Result<usize> {
    let f = &r.fields;
    let mut size = 36usize;
    let mut add = |n: usize| -> Result<()> {
        size = size
            .checked_add(n)
            .filter(|n| *n <= cap)
            .ok_or(Error::Limit("BAM encoded record bytes"))?;
        Ok(())
    };
    add(f.name.as_deref().unwrap_or("*").len() + 1)?;
    let long = f.cigar.len() > 65535;
    for op in &f.cigar {
        if op.len > 0x0fff_ffff {
            return Err(Error::Unsupported(
                "CIGAR operation exceeds BAM 28-bit length",
            ));
        }
    }
    add(if long {
        8
    } else {
        f.cigar
            .len()
            .checked_mul(4)
            .ok_or(Error::Limit("BAM CIGAR size"))?
    })?;
    let seq = f.sequence.as_ref().map(String::len).unwrap_or(0);
    add(seq.div_ceil(2))?;
    add(seq)?;
    for a in &f.aux {
        if a.tag == *b"CG" {
            return Err(Error::Unsupported(
                "SAM CG transport tag conflicts with BAM reconstruction",
            ));
        }
        add(match &a.value {
            AuxValue::Char(_) => 4,
            AuxValue::Integer(_) => 3 + width(a.wire_type)?,
            AuxValue::FloatBits(_) => 7,
            AuxValue::String(v) | AuxValue::Hex(v) => 4 + v.len(),
            AuxValue::Array {
                element_type,
                values,
            } => {
                8 + values
                    .len()
                    .checked_mul(width(*element_type)?)
                    .ok_or(Error::Limit("BAM array size"))?
            }
        })?;
    }
    if long {
        if f.cigar.len() > l.max_array_values || f.aux.len() >= l.max_aux_tags {
            return Err(Error::Limit("long CIGAR auxiliary transport"));
        }
        add(8 + f
            .cigar
            .len()
            .checked_mul(4)
            .ok_or(Error::Limit("long CIGAR bytes"))?)?;
    }
    Ok(size)
}
fn encode(
    r: &Record,
    h: &std::collections::BTreeMap<String, usize>,
    l: &Limits,
    cap: usize,
) -> Result<Vec<u8>> {
    let f = &r.fields;
    validate_fields(f, l)?;
    let size = encoded_size(r, l, cap)?;
    let mut b = Vec::with_capacity(size);
    b.extend([0; 4]);
    put_i32(&mut b, reference_id(f.reference.as_deref(), h)?);
    put_i32(&mut b, f.start.map(|v| v as i32).unwrap_or(-1));
    let name = f.name.as_deref().unwrap_or("*");
    b.extend([(name.len() + 1) as u8, f.mapping_quality.unwrap_or(255)]);
    let ref_len = f
        .cigar
        .iter()
        .filter(|o| o.kind.consumes_reference())
        .try_fold(0u64, |v, o| v.checked_add(u64::from(o.len)))
        .ok_or(Error::Invalid("reference CIGAR overflow"))?;
    let end = u64::from(f.start.unwrap_or(0))
        .checked_add(ref_len.max(1))
        .ok_or(Error::Invalid("reference span overflow"))?;
    let bin = if let Some(start) = f.start {
        reg2bin(
            u64::from(start),
            if f.flags & 4 != 0 {
                u64::from(start) + 1
            } else {
                end
            },
        )
    } else {
        4680
    };
    b.extend(bin.to_le_bytes());
    let long = f.cigar.len() > 65535;
    let len = f.sequence.as_ref().map(String::len).unwrap_or(0);
    if long && (ref_len > 0x0fff_ffff || len > 0x0fff_ffff) {
        return Err(Error::Unsupported(
            "long CIGAR placeholder exceeds 28-bit operation",
        ));
    }
    b.extend((if long { 2 } else { f.cigar.len() as u16 }).to_le_bytes());
    b.extend(f.flags.to_le_bytes());
    put_i32(&mut b, len as i32);
    put_i32(&mut b, reference_id(f.mate_reference.as_deref(), h)?);
    put_i32(&mut b, f.mate_start.map(|v| v as i32).unwrap_or(-1));
    put_i32(&mut b, f.template_length);
    b.extend(name.as_bytes());
    b.push(0);
    let word = |op: &CigarOp| -> Result<u32> {
        if op.len > 0x0fff_ffff {
            return Err(Error::Unsupported(
                "CIGAR operation exceeds BAM 28-bit length",
            ));
        }
        Ok(op.len << 4 | op.kind.code())
    };
    if long {
        b.extend(((len as u32) << 4 | 4).to_le_bytes());
        b.extend(((ref_len as u32) << 4 | 3).to_le_bytes());
    } else {
        for op in &f.cigar {
            b.extend(word(op)?.to_le_bytes());
        }
    }
    if let Some(seq) = &f.sequence {
        let mut nibble = 0;
        for (i, v) in seq.bytes().enumerate() {
            let code = b"=ACMGRSVTWYHKDBN"
                .iter()
                .position(|b| *b == v.to_ascii_uppercase())
                .ok_or(Error::Unsupported("SAM base not representable in BAM"))?
                as u8;
            if i % 2 == 0 {
                nibble = code << 4;
            } else {
                b.push(nibble | code);
            }
        }
        if seq.len() % 2 != 0 {
            b.push(nibble);
        }
    }
    if let Some(q) = &f.qualities {
        b.extend(q);
    } else {
        b.resize(b.len() + len, 255);
    }
    for a in &f.aux {
        if a.tag == *b"CG" {
            return Err(Error::Unsupported(
                "SAM CG transport tag conflicts with BAM reconstruction",
            ));
        }
        write_aux(&mut b, a)?;
        if b.len() > l.max_record_bytes {
            return Err(Error::Limit("BAM export bytes"));
        }
    }
    if long {
        if f.cigar.len() > l.max_array_values || f.aux.len() >= l.max_aux_tags {
            return Err(Error::Limit("long CIGAR auxiliary transport"));
        }
        let values = f.cigar.iter().map(word).collect::<Result<Vec<_>>>()?;
        write_aux(
            &mut b,
            &Aux {
                tag: *b"CG",
                wire_type: b'B',
                value: AuxValue::Array {
                    element_type: b'I',
                    values,
                },
            },
        )?;
    }
    if b.len() > l.max_record_bytes {
        return Err(Error::Limit("BAM export bytes"));
    }
    if b.len() != size {
        return Err(Error::Invalid("internal BAM encoded-size disagreement"));
    }
    let n = (b.len() - 4) as i32;
    b[..4].copy_from_slice(&n.to_le_bytes());
    Ok(b)
}
fn reg2bin(begin: u64, end: u64) -> u16 {
    let end = end - 1;
    for (shift, offset) in [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if begin >> shift == end >> shift {
            return ((offset + (begin >> shift)) & 0xffff) as u16;
        }
    }
    0
}
pub struct BamWriter<W: Write> {
    inner: BgzfWriter<W>,
    header_identity: [u8; 32],
    dictionary: std::sync::Arc<std::collections::BTreeMap<String, usize>>,
    limits: Limits,
    failed: bool,
    notices: Vec<ExportNotice>,
    records: u64,
    bytes: u64,
}
impl<W: Write> BamWriter<W> {
    pub fn new(inner: W, header: &Header, limits: Limits) -> Result<Self> {
        header.validate_policy(&limits)?;
        let mut inner = BgzfWriter::with_limit(inner, limits.max_input_bytes);
        if header.references.len() > limits.max_references {
            return Err(Error::Limit("export reference dictionary"));
        }
        let header_bytes;

        if let Some(b) = &header.bam_original {
            if b.len() > limits.max_header_bytes {
                return Err(Error::Limit("BAM header"));
            }
            header_bytes = b.len() as u64;
            if header_bytes > limits.max_uncompressed_bytes {
                return Err(Error::Limit("BAM export uncompressed bytes"));
            }
            inner.write_all(b)?;
        } else {
            let mut size = 12usize
                .checked_add(header.raw.len())
                .filter(|n| *n <= limits.max_header_bytes)
                .ok_or(Error::Limit("BAM encoded header bytes"))?;
            for reference in header.references.iter() {
                size = size
                    .checked_add(9)
                    .and_then(|n| n.checked_add(reference.name.len()))
                    .filter(|n| *n <= limits.max_header_bytes)
                    .ok_or(Error::Limit("BAM encoded header bytes"))?;
            }
            if size as u64 > limits.max_uncompressed_bytes {
                return Err(Error::Limit("BAM export uncompressed bytes"));
            }
            let mut b = Vec::with_capacity(size);
            b.extend(b"BAM\x01");
            put_i32(&mut b, header.raw.len() as i32);
            b.extend(&header.raw);
            put_i32(&mut b, header.references.len() as i32);
            for r in header.references.iter() {
                put_i32(&mut b, (r.name.len() + 1) as i32);
                b.extend(r.name.as_bytes());
                b.push(0);
                put_i32(&mut b, r.length as i32);
            }
            if b.len() > limits.max_header_bytes {
                return Err(Error::Limit("BAM header bytes"));
            }
            header_bytes = b.len() as u64;
            if header_bytes > limits.max_uncompressed_bytes {
                return Err(Error::Limit("BAM export uncompressed bytes"));
            }
            inner.write_all(&b)?;
        }
        Ok(Self {
            inner,
            header_identity: header.dictionary_identity,
            dictionary: header.reference_index.clone(),
            limits,
            failed: false,
            notices: if header.bam_original.is_some() {
                vec![ExportNotice::BgzfFramingNotPreserved]
            } else {
                vec![]
            },
            records: 0,
            bytes: header_bytes,
        })
    }
    pub fn notices(&self) -> &[ExportNotice] {
        &self.notices
    }
    pub fn write_record(&mut self, r: &Record) -> Result<Vec<ExportNotice>> {
        if self.failed {
            return Err(Error::Poisoned);
        }
        let result = self.write(r);
        if result.is_err() {
            self.failed = true;
        }
        if let Ok(notices) = &result {
            for notice in notices {
                if !self.notices.contains(notice) {
                    self.notices.push(notice.clone());
                }
            }
        }
        result
    }
    fn write(&mut self, r: &Record) -> Result<Vec<ExportNotice>> {
        validate_fields(&r.fields, &self.limits)?;
        self.records = checked_add(self.records, 1, self.limits.max_records, "export records")?;

        if r.dictionary != self.header_identity {
            return Err(Error::Invalid("record/header dictionary identity"));
        }
        if r.format == Format::Bam {
            if r.original.len() > self.limits.max_record_bytes {
                return Err(Error::Limit("BAM record bytes"));
            }
            self.bytes = checked_add(
                self.bytes,
                r.original.len() as u64,
                self.limits.max_uncompressed_bytes,
                "BAM export uncompressed bytes",
            )?;
            self.inner.write_all(&r.original)?;
            return Ok(vec![ExportNotice::BgzfFramingNotPreserved]);
        }
        let cap = self
            .limits
            .max_uncompressed_bytes
            .saturating_sub(self.bytes)
            .min(self.limits.max_record_bytes as u64) as usize;
        let b = encode(r, &self.dictionary, &self.limits, cap)?;
        self.bytes = checked_add(
            self.bytes,
            b.len() as u64,
            self.limits.max_uncompressed_bytes,
            "BAM export uncompressed bytes",
        )?;
        self.inner.write_all(&b)?;
        let mut n = vec![
            ExportNotice::TextFormattingNotPreserved,
            ExportNotice::AuxIntegerWidthSelected,
        ];
        if r.fields
            .sequence
            .as_ref()
            .is_some_and(|s| s.bytes().any(|b| b.is_ascii_lowercase()))
        {
            n.push(ExportNotice::SequenceCaseNormalized);
        }
        if r.fields.cigar.len() > 65535 {
            n.push(ExportNotice::LongCigarTransportConverted);
        }
        Ok(n)
    }
    pub fn finish(self) -> Result<W> {
        if self.failed {
            return Err(Error::Poisoned);
        }
        Ok(self.inner.finish()?)
    }
}
