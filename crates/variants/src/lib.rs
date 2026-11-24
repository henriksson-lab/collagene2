//! Immutable variant-file evidence, not variant calling or reference validation.
//! Typed INFO/FORMAT/genotype access uses the read-only Noodles Record trait.
//! Original lexical bytes are retained independently of normalized conversion.
use noodles_bcf as bcf;
pub use noodles_vcf as vcf;
use std::io::{self, Cursor, Read, Write};
use vcf::variant::io::Write as VariantWrite;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Vcf,
    Bcf,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExportFormat {
    Original,
    Vcf,
    Bcf,
}
#[derive(Clone, Debug)]
pub struct Limits {
    pub source_bytes: usize,
    pub decoded_bytes: usize,
    pub header_bytes: usize,
    pub record_bytes: usize,
    pub records: usize,
    pub samples: usize,
    pub retained_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            source_bytes: 16 << 20,
            decoded_bytes: 32 << 20,
            header_bytes: 1 << 20,
            record_bytes: 256 << 10,
            records: 50_000,
            samples: 256,
            retained_bytes: 128 << 20,
        }
    }
}
impl Limits {
    fn validate(&self) -> io::Result<()> {
        if self.source_bytes == 0
            || self.source_bytes > 64 << 20
            || self.decoded_bytes == 0
            || self.decoded_bytes > 64 << 20
            || self.header_bytes == 0
            || self.header_bytes > 4 << 20
            || self.record_bytes == 0
            || self.record_bytes > 1 << 20
            || self.records > 100_000
            || self.samples > 1024
            || self.retained_bytes == 0
            || self.retained_bytes > 256 << 20
        {
            return Err(invalid("limits outside hard bounds"));
        }
        Ok(())
    }
}
fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s)
}
fn bounded(n: usize, cap: usize) -> io::Result<usize> {
    if n > cap {
        Err(invalid("variant resource limit"))
    } else {
        Ok(n)
    }
}
fn add(a: usize, b: usize, cap: usize) -> io::Result<usize> {
    bounded(
        a.checked_add(b).ok_or_else(|| invalid("size overflow"))?,
        cap,
    )
}
struct Sink {
    bytes: Vec<u8>,
    cap: usize,
}
impl Sink {
    fn new(cap: usize) -> Self {
        Self {
            bytes: Vec::new(),
            cap,
        }
    }
}
impl Write for Sink {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        add(self.bytes.len(), b.len(), self.cap)?;
        self.bytes.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
enum Stored {
    Vcf(vcf::Record),
    Bcf(bcf::Record),
}
impl Stored {
    fn record(&self) -> &dyn vcf::variant::Record {
        match self {
            Self::Vcf(r) => r,
            Self::Bcf(r) => r,
        }
    }
}
/// Private complete parser state. There is intentionally no Deserialize or mutable accessor.
pub struct Collection {
    original: Vec<u8>,
    hash: String,
    format: Format,
    header: vcf::Header,
    records: Vec<Stored>,
    limits: Limits,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Row {
    pub index: usize,
    pub reference_name: String,
    pub position_1based: Option<usize>,
    pub vcf: String,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Export {
    pub bytes: Vec<u8>,
    pub notices: Vec<String>,
}
impl Collection {
    pub fn import(bytes: &[u8], format: Format, limits: &Limits) -> io::Result<Self> {
        limits.validate()?;
        bounded(bytes.len(), limits.source_bytes)?;
        // Decoded staging + original + bounded raw records + canonical validation
        // scratch are charged conservatively before copying or decoding input.
        let estimated = bytes
            .len()
            .checked_add(
                limits
                    .decoded_bytes
                    .checked_mul(3)
                    .ok_or_else(|| invalid("size overflow"))?,
            )
            .and_then(|n| n.checked_add(limits.records.checked_mul(256)?))
            .and_then(|n| n.checked_add(limits.record_bytes))
            .ok_or_else(|| invalid("size overflow"))?;
        bounded(estimated, limits.retained_bytes)?;
        let decoded;
        let raw = if format == Format::Bcf && bytes.starts_with(&[31, 139]) {
            decoded = decode_bgzf(bytes, limits.decoded_bytes)?;
            decoded.as_slice()
        } else {
            bounded(bytes.len(), limits.decoded_bytes)?;
            bytes
        };
        let (count, samples) = match format {
            Format::Vcf => preflight_vcf(raw, limits)?,
            Format::Bcf => preflight_bcf(raw, limits)?,
        };
        bounded(count, limits.records)?;
        bounded(samples, limits.samples)?;
        let mut records = Vec::with_capacity(count);
        let header = match format {
            Format::Vcf => {
                let mut r = vcf::io::Reader::new(Cursor::new(raw));
                let h = r.read_header()?;
                for _ in 0..count {
                    let mut x = vcf::Record::default();
                    if r.read_record(&mut x)? == 0 {
                        return Err(invalid("truncated VCF"));
                    }
                    validate_record(&h, &x, limits.record_bytes)?;
                    records.push(Stored::Vcf(x));
                }
                h
            }
            Format::Bcf => {
                let mut r = bcf::io::Reader::from(Cursor::new(raw));
                let h = r.read_header()?;
                preflight_bcf_formats(raw, &h)?;
                for _ in 0..count {
                    let mut x = bcf::Record::default();
                    if r.read_record(&mut x)? == 0 {
                        return Err(invalid("truncated BCF"));
                    }
                    validate_record(&h, &x, limits.record_bytes)?;
                    records.push(Stored::Bcf(x));
                }
                h
            }
        };
        if header.sample_names().len() != samples {
            return Err(invalid("sample dictionary mismatch"));
        }
        Ok(Self {
            original: bytes.to_vec(),
            hash: blake3::hash(bytes).to_hex().to_string(),
            format,
            header,
            records,
            limits: limits.clone(),
        })
    }
    pub fn source_blake3(&self) -> &str {
        &self.hash
    }
    pub fn original_bytes(&self) -> &[u8] {
        &self.original
    }
    pub fn format(&self) -> Format {
        self.format
    }
    pub fn header(&self) -> &vcf::Header {
        &self.header
    }
    pub fn len(&self) -> usize {
        self.records.len()
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
    pub fn record(&self, index: usize) -> Option<&dyn vcf::variant::Record> {
        self.records.get(index).map(Stored::record)
    }
    pub fn page(&self, offset: usize, limit: usize) -> io::Result<Vec<Row>> {
        bounded(limit, 100)?;
        if offset > self.len() {
            return Err(invalid("page offset"));
        }
        let mut rows = Vec::with_capacity(limit.min(self.len() - offset));
        let mut page_bytes = 0;
        for index in offset..offset.saturating_add(limit).min(self.len()) {
            let r = self.records[index].record();
            let remaining = (4usize << 20) - page_bytes;
            let text = render(&self.header, r, self.limits.record_bytes.min(remaining))?;
            let reference_name = r.reference_sequence_name(&self.header)?;
            page_bytes = add(page_bytes, text.len(), 4 << 20)?;
            page_bytes = add(page_bytes, reference_name.len(), 4 << 20)?;
            rows.push(Row {
                index,
                reference_name: reference_name.into(),
                position_1based: r.variant_start().transpose()?.map(usize::from),
                vcf: String::from_utf8(text).map_err(|_| invalid("VCF UTF-8"))?,
            });
        }
        Ok(rows)
    }
    /// Conversion normalizes syntax. Original is the only byte-preserving mode.
    /// Any codec, representation or output-limit error returns no artifact.
    pub fn export(&self, format: ExportFormat, max_bytes: usize) -> io::Result<Export> {
        bounded(max_bytes, 64 << 20)?;
        let notices = if format == ExportFormat::Original {
            vec![]
        } else {
            vec!["Normalized exchange: header ordering/indices, numeric lexical spelling, BCF integer widths and missing/vector-end encodings may change. Exact original bytes remain separate. Records are declared evidence, not reference-validated calls.".into()]
        };
        let bytes = match format {
            ExportFormat::Original => {
                bounded(self.original.len(), max_bytes)?;
                self.original.clone()
            }
            ExportFormat::Vcf => {
                let mut w = vcf::io::Writer::new(Sink::new(max_bytes));
                w.write_header(&self.header)?;
                for r in &self.records {
                    w.write_variant_record(&self.header, r.record())?;
                }
                w.into_inner().bytes
            }
            ExportFormat::Bcf => {
                let mut sink = Sink::new(max_bytes);
                {
                    let mut w = bcf::io::Writer::new(&mut sink);
                    w.write_header(&self.header)?;
                    for r in &self.records {
                        // Noodles 0.85's BCF INFO encoder has an unimplemented
                        // branch for a present key with a wholly missing value.
                        // Reject the conversion rather than panic or drop it.
                        for field in r.record().info().iter(&self.header) {
                            if field?.1.is_none() {
                                return Err(invalid(
                                    "BCF conversion of wholly missing INFO values is unsupported by this codec; use Original or VCF",
                                ));
                            }
                        }
                        w.write_variant_record(&self.header, r.record())?;
                    }
                    w.try_finish()?;
                }
                sink.bytes
            }
        };
        Ok(Export { bytes, notices })
    }
}
fn render(h: &vcf::Header, r: &dyn vcf::variant::Record, cap: usize) -> io::Result<Vec<u8>> {
    let mut w = vcf::io::Writer::new(Sink::new(cap));
    w.write_variant_record(h, r)?;
    Ok(w.into_inner().bytes)
}
fn validate_record(h: &vcf::Header, r: &dyn vcf::variant::Record, cap: usize) -> io::Result<()> {
    render(h, r, cap).map(|_| ())
}
fn preflight_vcf(raw: &[u8], l: &Limits) -> io::Result<(usize, usize)> {
    if !raw.starts_with(b"##fileformat=VCFv4.") {
        return Err(invalid(
            "explicit VCF requires VCF4 header; other formats unsupported",
        ));
    }
    let mut header = 0;
    let mut header_lines = 0;
    let mut count = 0;
    let mut samples = None;
    let mut body = false;
    // PASS implicitly occupies string slot 0. Sparse IDX plus following implicit
    // entries must remain bounded even when every individual IDX is legal.
    let mut implicit = [1usize, 0];
    let mut explicit = [0usize; 2];
    for line in raw.split_inclusive(|b| *b == b'\n') {
        bounded(line.len(), l.record_bytes)?;
        if line.starts_with(b"#") {
            if body {
                return Err(invalid("header after records"));
            }
            header = add(header, line.len(), l.header_bytes)?;
            header_lines = add(header_lines, 1, 4096)?;
            if let Some((group, idx)) = header_idx(line)? {
                if let Some(idx) = idx {
                    explicit[group] = explicit[group].max(idx + 1);
                } else {
                    implicit[group] = add(implicit[group], 1, 4096)?;
                }
                add(explicit[group], implicit[group], 4096)?;
            }
            if line.starts_with(b"#CHROM\t") {
                if samples.is_some() {
                    return Err(invalid("duplicate column header"));
                }
                let cols = line.split(|b| *b == b'\t').count();
                if cols < 8 || cols == 9 {
                    return Err(invalid("sample columns"));
                }
                let n = cols.saturating_sub(9);
                bounded(n, l.samples)?;
                samples = Some(n);
            }
        } else {
            body = true;
            let n = samples.ok_or_else(|| invalid("missing column header"))?;
            let cols = line.split(|b| *b == b'\t').count();
            if cols != if n == 0 { 8 } else { 9 + n } {
                return Err(invalid("record sample column count"));
            }
            count = add(count, 1, l.records)?;
        }
    }
    Ok((
        count,
        samples.ok_or_else(|| invalid("missing column header"))?,
    ))
}
// Scan only actual structured fields, not strings such as Description="IDX=...".
// No header parser/string-map allocation is allowed before this check.
fn header_idx(line: &[u8]) -> io::Result<Option<(usize, Option<usize>)>> {
    let mut entry = None;
    for (prefix, group) in [
        (b"##INFO=<".as_slice(), 0),
        (b"##FORMAT=<".as_slice(), 0),
        (b"##FILTER=<".as_slice(), 0),
        (b"##contig=<".as_slice(), 1),
    ] {
        if let Some(rest) = line.strip_prefix(prefix) {
            entry = Some((rest, group));
            break;
        }
    }
    let Some((rest, group)) = entry else {
        return Ok(None);
    };
    let rest = rest
        .trim_ascii()
        .strip_suffix(b">")
        .ok_or_else(|| invalid("structured header terminator"))?;
    let mut quoted = false;
    let mut escaped = false;
    let mut start = 0;
    let mut idx = None;
    for i in 0..=rest.len() {
        if i < rest.len() {
            let b = rest[i];
            if escaped {
                escaped = false;
                continue;
            }
            if quoted && b == b'\\' {
                escaped = true;
                continue;
            }
            if b == b'"' {
                quoted = !quoted;
                continue;
            }
            if b != b',' || quoted {
                continue;
            }
        }
        if quoted || escaped {
            return Err(invalid("unterminated structured header string"));
        }
        let field = rest[start..i].trim_ascii();
        start = i + 1;
        if let Some(eq) = field.iter().position(|b| *b == b'=') {
            if field[..eq].trim_ascii() == b"IDX" {
                if idx.is_some() {
                    return Err(invalid("duplicate header IDX"));
                }
                let value = field[eq + 1..].trim_ascii();
                if value.is_empty() {
                    return Err(invalid("empty header IDX"));
                }
                let mut n = 0usize;
                for b in value {
                    if !b.is_ascii_digit() {
                        return Err(invalid("invalid header IDX"));
                    }
                    n = n
                        .checked_mul(10)
                        .and_then(|n| n.checked_add(usize::from(*b - b'0')))
                        .filter(|n| *n < 4096)
                        .ok_or_else(|| {
                            invalid("header IDX exceeds 4096 dictionary-slot resource limit")
                        })?;
                }
                idx = Some(n);
            }
        }
    }
    Ok(Some((group, idx)))
}
fn u32le(b: &[u8]) -> io::Result<usize> {
    Ok(u32::from_le_bytes(
        b.get(..4)
            .ok_or_else(|| invalid("truncated length"))?
            .try_into()
            .unwrap(),
    ) as usize)
}
fn preflight_bcf(raw: &[u8], l: &Limits) -> io::Result<(usize, usize)> {
    if !raw.starts_with(b"BCF\x02\x02") {
        return Err(invalid(
            "explicit BCF requires BCF2.2; CRAM/BCF1 unsupported",
        ));
    }
    let n = u32le(raw.get(5..).ok_or_else(|| invalid("BCF header"))?)?;
    bounded(n, l.header_bytes)?;
    let mut at = add(9, n, raw.len())?;
    let text = &raw[9..at];
    if text.last() != Some(&0) {
        return Err(invalid("BCF header terminator"));
    }
    let (_, samples) = preflight_vcf(&text[..text.len() - 1], l)?;
    let mut count = 0;
    while at < raw.len() {
        let shared = u32le(&raw[at..])?;
        let individual = u32le(
            raw.get(at + 4..)
                .ok_or_else(|| invalid("BCF record header"))?,
        )?;
        if shared < 24 {
            return Err(invalid("BCF shared fields"));
        }
        let size = add(shared, individual, l.record_bytes)?;
        let end = add(at, add(8, size, usize::MAX)?, raw.len())?;
        let s = u32le(&raw[at + 28..])? & 0x00ff_ffff;
        if s != samples {
            return Err(invalid("BCF record sample count"));
        }
        preflight_bcf_values(
            &raw[at + 8..at + 8 + shared],
            &raw[at + 8 + shared..end],
            samples,
        )?;
        count = add(count, 1, l.records)?;
        at = end;
    }
    Ok((count, samples))
}
fn take<'a>(src: &mut &'a [u8], n: usize) -> io::Result<&'a [u8]> {
    let part = src
        .get(..n)
        .ok_or_else(|| invalid("truncated BCF typed payload"))?;
    *src = &src[n..];
    Ok(part)
}
fn positive_integer(src: &mut &[u8]) -> io::Result<usize> {
    let descriptor = take(src, 1)?[0];
    let n = match descriptor {
        0x11 => i32::from(take(src, 1)?[0] as i8),
        0x12 => i32::from(i16::from_le_bytes(take(src, 2)?.try_into().unwrap())),
        0x13 => i32::from_le_bytes(take(src, 4)?.try_into().unwrap()),
        _ => return Err(invalid("BCF length/key must be one integer")),
    };
    usize::try_from(n).map_err(|_| invalid("negative/missing BCF length/key"))
}
fn descriptor(src: &mut &[u8]) -> io::Result<(u8, usize, usize)> {
    let byte = take(src, 1)?[0];
    let ty = byte & 15;
    let width = match ty {
        0 => 0,
        1 | 7 => 1,
        2 => 2,
        3 | 5 => 4,
        _ => return Err(invalid("unsupported BCF value type")),
    };
    let len = if byte >> 4 == 15 {
        positive_integer(src)?
    } else {
        usize::from(byte >> 4)
    };
    if ty == 0 && len != 0 {
        return Err(invalid("BCF null type with nonzero count"));
    }
    Ok((ty, len, width))
}
fn typed(src: &mut &[u8], multiplicity: usize) -> io::Result<(u8, usize)> {
    let (ty, len, width) = descriptor(src)?;
    let bytes = len
        .checked_mul(width)
        .and_then(|n| n.checked_mul(multiplicity))
        .ok_or_else(|| invalid("BCF typed count overflow"))?;
    // Limits are tied to actual encoded payload, not forged count descriptors.
    take(src, bytes)?;
    Ok((ty, len))
}
fn key(src: &mut &[u8]) -> io::Result<()> {
    bounded(positive_integer(src)?, 4095).map(|_| ())
}
fn preflight_bcf_values(shared: &[u8], individual: &[u8], samples: usize) -> io::Result<()> {
    let info = usize::from(u16::from_le_bytes(shared[16..18].try_into().unwrap()));
    let alleles = usize::from(u16::from_le_bytes(shared[18..20].try_into().unwrap()));
    let formats = usize::from(shared[23]);
    bounded(info, 4096)?;
    bounded(alleles, 4096)?;
    if alleles == 0 || (samples == 0 && formats != 0) {
        return Err(invalid("BCF allele/sample dimensions"));
    }
    let mut site = &shared[24..];
    // Noodles raw Fields::index slices these strings without a fallible bound
    // check. Validate all lengths here before its record reader is invoked.
    for _ in 0..=alleles {
        if typed(&mut site, 1)?.0 != 7 {
            return Err(invalid("BCF ID/allele string type"));
        }
    }
    let filter = typed(&mut site, 1)?.0;
    if !matches!(filter, 0..=3) {
        return Err(invalid("BCF filter integer type"));
    }
    for _ in 0..info {
        key(&mut site)?;
        typed(&mut site, 1)?;
    }
    if !site.is_empty() {
        return Err(invalid("unconsumed BCF shared payload"));
    }
    let mut values = individual;
    for _ in 0..formats {
        key(&mut values)?;
        typed(&mut values, samples)?;
    }
    if !values.is_empty() {
        return Err(invalid("unconsumed BCF FORMAT payload"));
    }
    Ok(())
}
fn preflight_bcf_formats(raw: &[u8], header: &vcf::Header) -> io::Result<()> {
    use vcf::header::record::value::map::format::{Number, Type};
    let mut at = 9 + u32le(&raw[5..])?;
    while at < raw.len() {
        let shared = u32le(&raw[at..])?;
        let individual = u32le(&raw[at + 4..])?;
        let samples = u32le(&raw[at + 28..])? & 0x00ff_ffff;
        let formats = usize::from(raw[at + 31]);
        let mut src = &raw[at + 8 + shared..at + 8 + shared + individual];
        for _ in 0..formats {
            let id = positive_integer(&mut src)?;
            let name = header
                .string_maps()
                .strings()
                .get_index(id)
                .ok_or_else(|| invalid("undefined BCF FORMAT index"))?;
            let definition = header
                .formats()
                .get(name)
                .ok_or_else(|| invalid("missing BCF FORMAT definition"))?;
            let (ty, len, width) = descriptor(&mut src)?;
            if len == 0 || matches!(definition.number(), Number::Count(0)) {
                return Err(invalid("unsupported empty BCF FORMAT value"));
            }
            let scalar = matches!(definition.number(), Number::Count(1));
            let valid = if name == "GT" {
                ty == 1
            } else {
                match definition.ty() {
                    Type::Integer => matches!(ty, 1..=3) && (!scalar || len == 1),
                    Type::Float => ty == 5 && (!scalar || len == 1),
                    Type::Character | Type::String => ty == 7,
                }
            };
            if !valid {
                return Err(invalid(
                    "unsupported BCF FORMAT header/type/width combination",
                ));
            }
            let payload = take(
                &mut src,
                len.checked_mul(width)
                    .and_then(|n| n.checked_mul(samples))
                    .ok_or_else(|| invalid("BCF sample count overflow"))?,
            )?;
            if ty == 7 {
                for sample in payload.chunks_exact(len) {
                    let sample =
                        &sample[..sample.iter().position(|b| *b == 0).unwrap_or(sample.len())];
                    std::str::from_utf8(sample).map_err(|_| invalid("invalid BCF sample UTF-8"))?;
                    if matches!(definition.ty(), Type::Character)
                        && sample.split(|b| *b == b',').any(|part| part.is_empty())
                    {
                        return Err(invalid("empty BCF Character value"));
                    }
                }
            } else if name == "GT" {
                if payload.iter().any(|b| *b >= 128 && *b != 129) {
                    return Err(invalid("reserved BCF genotype code"));
                }
            } else if scalar {
                for value in payload.chunks_exact(width) {
                    let reserved = match ty {
                        1 => value[0] > 128 && value[0] < 136,
                        2 => {
                            let n = u16::from_le_bytes(value.try_into().unwrap());
                            n > 32768 && n < 32776
                        }
                        3 => {
                            let n = u32::from_le_bytes(value.try_into().unwrap());
                            n > 0x80000000 && n < 0x80000008
                        }
                        5 => {
                            let n = u32::from_le_bytes(value.try_into().unwrap());
                            n > 0x7f800001 && n < 0x7f800008
                        }
                        _ => false,
                    };
                    if reserved {
                        return Err(invalid("reserved/vector-end scalar BCF FORMAT value"));
                    }
                }
            }
        }
        at += 8 + shared + individual;
    }
    Ok(())
}
// Decode independent <=64KiB BGZF members only after framing and declared-size
// checks; never hand a whole unbounded gzip stream to a codec.
fn decode_bgzf(bytes: &[u8], cap: usize) -> io::Result<Vec<u8>> {
    let mut out = Vec::new();
    let mut at = 0;
    let mut eof = false;
    while at < bytes.len() {
        if eof {
            return Err(invalid("data after BGZF EOF"));
        }
        let h = bytes
            .get(at..at + 12)
            .ok_or_else(|| invalid("truncated BGZF header"))?;
        if h[..4] != [31, 139, 8, 4] {
            return Err(invalid("BGZF framing"));
        }
        let xlen = u16::from_le_bytes([h[10], h[11]]) as usize;
        let extra_end = add(at + 12, xlen, bytes.len())?;
        let mut p = at + 12;
        let mut block = None;
        while p < extra_end {
            let x = bytes
                .get(p..p + 4)
                .filter(|_| p + 4 <= extra_end)
                .ok_or_else(|| invalid("BGZF extra"))?;
            let n = u16::from_le_bytes([x[2], x[3]]) as usize;
            let e = add(p + 4, n, extra_end)?;
            if &x[..2] == b"BC" {
                if n != 2 || block.is_some() {
                    return Err(invalid("BGZF BC field"));
                }
                block = Some(u16::from_le_bytes([bytes[p + 4], bytes[p + 5]]) as usize + 1);
            }
            p = e;
        }
        let end = add(
            at,
            block.ok_or_else(|| invalid("BGZF BC absent"))?,
            bytes.len(),
        )?;
        if end < extra_end + 8 {
            return Err(invalid("BGZF block length"));
        }
        let isize = u32le(&bytes[end - 4..end])?;
        bounded(isize, 65536)?;
        add(out.len(), isize, cap)?;
        let mut d = flate2::read::DeflateDecoder::new(&bytes[extra_end..end - 8]);
        let mut chunk = Vec::with_capacity(isize);
        d.by_ref().take(65537).read_to_end(&mut chunk)?;
        if chunk.len() != isize || d.total_in() as usize != end - 8 - extra_end {
            return Err(invalid("BGZF inflated size/trailing deflate bytes"));
        }
        if crc32fast::hash(&chunk) != u32le(&bytes[end - 8..end - 4])? as u32 {
            return Err(invalid("BGZF CRC"));
        }
        eof = chunk.is_empty();
        out.extend_from_slice(&chunk);
        at = end;
    }
    if !eof {
        return Err(invalid("missing BGZF EOF"));
    }
    Ok(out)
}
