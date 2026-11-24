use crate::*;
use std::io::{BufRead, Write};
pub(crate) fn line<R: BufRead>(r: &mut R, cap: usize) -> Result<Option<Vec<u8>>> {
    let mut out = Vec::new();
    loop {
        let b = r.fill_buf()?;
        if b.is_empty() {
            return Ok((!out.is_empty()).then_some(out));
        }
        let n = b
            .iter()
            .position(|v| *v == b'\n')
            .map(|i| i + 1)
            .unwrap_or(b.len());
        if out.len().checked_add(n).is_none_or(|v| v > cap) {
            return Err(Error::Limit("line bytes"));
        }
        let done = b[n - 1] == b'\n';
        out.extend_from_slice(&b[..n]);
        r.consume(n);
        if done {
            return Ok(Some(out));
        }
    }
}
fn trim_line(b: &[u8]) -> &[u8] {
    let b = b.strip_suffix(b"\n").unwrap_or(b);
    b.strip_suffix(b"\r").unwrap_or(b)
}
fn text(b: &[u8]) -> Result<&str> {
    std::str::from_utf8(b).map_err(|_| Error::Invalid("UTF-8"))
}
fn int<T: std::str::FromStr>(b: &[u8]) -> Result<T> {
    text(b)?.parse().map_err(|_| Error::Invalid("integer"))
}
pub(crate) fn reference_name(s: &str) -> bool {
    !s.is_empty()
        && !s.starts_with('*')
        && !s.starts_with('=')
        && s.bytes()
            .all(|b| (33..=126).contains(&b) && !b"\\,\"`'()[]{}<>".contains(&b))
}
pub(crate) fn parse_header(bytes: &[u8], l: &Limits) -> Result<Header> {
    l.validate()?;
    if bytes.len() > l.max_header_bytes {
        return Err(Error::Limit("header bytes"));
    }
    let mut refs = Vec::new();
    let mut names = std::collections::BTreeSet::new();
    let mut seen_hd = false;
    for (index, raw) in bytes.split_inclusive(|v| *v == b'\n').enumerate() {
        if raw.len() > l.max_record_bytes {
            return Err(Error::Limit("header line bytes"));
        }
        let b = trim_line(raw);
        if b.len() < 3
            || b[0] != b'@'
            || !b[1..3].iter().all(u8::is_ascii_uppercase)
            || b.get(3).is_some_and(|v| *v != b'\t')
        {
            return Err(Error::Invalid("header line"));
        }
        text(b)?;
        if b.iter().any(|c| *c < 32 && *c != b'\t' || *c == 127) {
            return Err(Error::Invalid("header control character"));
        }
        if &b[..3] == b"@CO" {
            continue;
        }
        if &b[..3] == b"@HD" {
            if index != 0 || seen_hd {
                return Err(Error::Invalid("header HD ordering"));
            }
            seen_hd = true;
        }
        let mut tags = std::collections::BTreeMap::new();
        for part in b.split(|v| *v == b'\t').skip(1) {
            if part.len() < 4
                || part[2] != b':'
                || !part[0].is_ascii_alphabetic()
                || !part[1].is_ascii_alphanumeric()
                || tags.insert(&part[..2], &part[3..]).is_some()
            {
                return Err(Error::Invalid("header tag"));
            }
        }
        if &b[..3] == b"@HD" {
            let vn = tags
                .get(b"VN".as_slice())
                .ok_or(Error::Invalid("HD VN required"))?;
            let mut p = vn.split(|b| *b == b'.');
            if !p
                .next()
                .is_some_and(|s| !s.is_empty() && s.iter().all(u8::is_ascii_digit))
                || !p
                    .next()
                    .is_some_and(|s| !s.is_empty() && s.iter().all(u8::is_ascii_digit))
                || p.next().is_some()
            {
                return Err(Error::Invalid("HD VN version"));
            }
        }
        if matches!(&b[..3], b"@RG" | b"@PG")
            && !tags.get(b"ID".as_slice()).is_some_and(|id| !id.is_empty())
        {
            return Err(Error::Invalid("RG/PG ID required"));
        }
        if &b[..3] == b"@SQ" {
            if refs.len() >= l.max_references {
                return Err(Error::Limit("reference dictionary"));
            }
            let name = text(
                tags.get(b"SN".as_slice())
                    .ok_or(Error::Invalid("SQ SN required"))?,
            )?
            .to_string();
            let length: u32 = int(tags
                .get(b"LN".as_slice())
                .ok_or(Error::Invalid("SQ LN required"))?)?;
            if !reference_name(&name)
                || length == 0
                || length > i32::MAX as u32
                || !names.insert(name.clone())
            {
                return Err(Error::Invalid("reference name/length/duplicate"));
            }
            if let Some(aliases) = tags.get(b"AN".as_slice()) {
                for alias in aliases.split(|b| *b == b',') {
                    let alias = text(alias)?;
                    if !reference_name(alias) || !names.insert(alias.to_string()) {
                        return Err(Error::Invalid("reference alias identity collision"));
                    }
                }
            }
            let circular = match tags.get(b"TP".as_slice()).copied() {
                Some(b"circular") => true,
                Some(b"linear") | None => false,
                Some(_) => return Err(Error::Unsupported("SQ topology")),
            };
            refs.push(Reference {
                name,
                length,
                circular,
            });
        }
    }
    Ok(Header::new(bytes.to_vec(), refs))
}
fn parse_cigar(b: &[u8], l: &Limits) -> Result<Vec<CigarOp>> {
    if b == b"*" {
        return Ok(vec![]);
    }
    let mut out = Vec::new();
    let mut n = 0u32;
    let mut digits = 0;
    for c in b {
        if c.is_ascii_digit() {
            n = n
                .checked_mul(10)
                .and_then(|v| v.checked_add(u32::from(*c - b'0')))
                .ok_or(Error::Invalid("CIGAR length overflow"))?;
            digits += 1;
        } else {
            if digits == 0 || out.len() >= l.max_cigar_ops {
                return Err(Error::Limit("CIGAR operators"));
            }
            let code = b"MIDNSHP=X"
                .iter()
                .position(|v| v == c)
                .ok_or(Error::Unsupported("CIGAR operator"))?;
            out.push(CigarOp {
                len: n,
                kind: CigarKind::from_code(code as u32)?,
            });
            n = 0;
            digits = 0;
        }
    }
    if digits != 0 || out.is_empty() {
        return Err(Error::Invalid("CIGAR"));
    }
    Ok(out)
}
fn scalar_bits(kind: u8, value: &[u8]) -> Result<u32> {
    Ok(match kind {
        b'c' => int::<i8>(value)? as u8 as u32,
        b'C' => u32::from(int::<u8>(value)?),
        b's' => int::<i16>(value)? as u16 as u32,
        b'S' => u32::from(int::<u16>(value)?),
        b'i' => int::<i32>(value)? as u32,
        b'I' => int::<u32>(value)?,
        b'f' => {
            let v: f32 = int(value)?;
            if !v.is_finite() {
                return Err(Error::Unsupported("non-finite SAM float"));
            }
            v.to_bits()
        }
        _ => return Err(Error::Unsupported("auxiliary array type")),
    })
}
fn aux(part: &[u8], l: &Limits) -> Result<Aux> {
    if part.len() < 5 || part[2] != b':' || part[4] != b':' {
        return Err(Error::Invalid("auxiliary field"));
    }
    let tag = [part[0], part[1]];
    let kind = part[3];
    let v = &part[5..];
    let value = match kind {
        b'A' => {
            if v.len() != 1 || !(33..=126).contains(&v[0]) {
                return Err(Error::Invalid("auxiliary character"));
            }
            AuxValue::Char(v[0])
        }
        b'i' => {
            let n: i64 = int(v)?;
            if n < i32::MIN as i64 || n > u32::MAX as i64 {
                return Err(Error::Invalid("SAM auxiliary integer range"));
            }
            AuxValue::Integer(n)
        }
        b'f' => AuxValue::FloatBits(scalar_bits(b'f', v)?),
        b'Z' => {
            if v.iter().any(|b| *b < 32 || *b == 127) {
                return Err(Error::Invalid("auxiliary string"));
            }
            text(v)?;
            AuxValue::String(v.to_vec())
        }
        b'H' => {
            if !v.len().is_multiple_of(2)
                || !v
                    .iter()
                    .all(|b| b.is_ascii_digit() || (b'A'..=b'F').contains(b))
            {
                return Err(Error::Invalid("hexadecimal auxiliary value"));
            }
            AuxValue::Hex(v.to_vec())
        }
        b'B' => {
            if v.is_empty() || !b"cCsSiIf".contains(&v[0]) || v.len() > 1 && v[1] != b',' {
                return Err(Error::Invalid("auxiliary array"));
            }
            let mut values = Vec::new();
            if v.len() > 1 {
                for n in v[2..].split(|b| *b == b',') {
                    if values.len() >= l.max_array_values {
                        return Err(Error::Limit("array values"));
                    }
                    values.push(scalar_bits(v[0], n)?);
                }
            }
            AuxValue::Array {
                element_type: v[0],
                values,
            }
        }
        _ => return Err(Error::Unsupported("SAM auxiliary type")),
    };
    Ok(Aux {
        tag,
        wire_type: kind,
        value,
    })
}
pub(crate) fn parse_record(raw: Vec<u8>, h: &Header, l: &Limits) -> Result<Record> {
    let b = trim_line(&raw);
    if b.starts_with(b"@") {
        return Err(Error::Invalid("header after alignment"));
    }
    let p = b
        .split(|v| *v == b'\t')
        .take(12 + l.max_aux_tags)
        .collect::<Vec<_>>();
    if p.len() < 11 || p.len() > 11 + l.max_aux_tags {
        return Err(Error::Invalid("record field count"));
    }
    let flags: u16 = int(p[1])?;
    let name = if p[0] == b"*" {
        None
    } else {
        let name = text(p[0])?;
        if name.len() > 254
            || name.is_empty()
            || !name.bytes().all(|b| (33..=126).contains(&b) && b != b'@')
        {
            return Err(Error::Invalid("QNAME"));
        }
        Some(name.to_string())
    };
    let r = if p[2] == b"*" {
        None
    } else {
        let n = text(p[2])?;
        if !reference_name(n) {
            return Err(Error::Invalid("RNAME"));
        }
        Some(n.to_string())
    };
    let position = |b: &[u8]| -> Result<Option<u32>> {
        let n: u32 = int(b)?;
        if n > i32::MAX as u32 {
            return Err(Error::Unsupported(
                "position exceeds BAM coordinate profile",
            ));
        }
        Ok(n.checked_sub(1))
    };
    let mapq: u8 = int(p[4])?;
    let mate = match p[6] {
        b"*" => None,
        b"=" => Some(
            r.clone()
                .ok_or(Error::Invalid("mate = without reference"))?,
        ),
        v => {
            let n = text(v)?;
            if !reference_name(n) {
                return Err(Error::Invalid("RNEXT"));
            }
            Some(n.into())
        }
    };
    if p[9].len() > l.max_sequence_bases && p[9] != b"*" {
        return Err(Error::Limit("sequence bases"));
    }
    let sequence = if p[9] == b"*" {
        None
    } else {
        Some(text(p[9])?.into())
    };
    let qualities = if p[10] == b"*" {
        None
    } else {
        if p[10].iter().any(|b| !(33..=126).contains(b)) {
            return Err(Error::Invalid("SAM quality encoding"));
        }
        Some(p[10].iter().map(|v| v - 33).collect())
    };
    let fields = Fields {
        name,
        flags,
        reference: r,
        start: position(p[3])?,
        mapping_quality: (mapq != 255).then_some(mapq),
        cigar: parse_cigar(p[5], l)?,
        mate_reference: mate,
        mate_start: position(p[7])?,
        template_length: int(p[8])?,
        sequence,
        qualities,
        aux: p[11..].iter().map(|p| aux(p, l)).collect::<Result<_>>()?,
    };
    if !h.references.is_empty()
        && [&fields.reference, &fields.mate_reference]
            .iter()
            .any(|name| {
                name.as_ref()
                    .is_some_and(|n| !h.references.iter().any(|r| r.name == *n))
            })
    {
        return Err(Error::Invalid("reference absent from SQ dictionary"));
    }
    validate_fields(&fields, l)?;
    Ok(Record {
        fields,
        format: Format::Sam,
        original: raw,
        dictionary: h.dictionary_identity,
        virtual_offset: None,
        long_cigar: false,
    })
}
pub struct SamReader<R> {
    inner: R,
    header: Header,
    limits: Limits,
    pending: Option<Vec<u8>>,
    bytes: u64,
    records: u64,
    emitted: u64,
    hash: blake3::Hasher,
    failed: bool,
    complete: bool,
}
impl<R: BufRead> SamReader<R> {
    pub fn new(mut inner: R, limits: Limits) -> Result<Self> {
        limits.validate()?;
        let mut raw = Vec::new();
        let mut pending = None;
        let mut hash = blake3::Hasher::new();
        let mut bytes = 0;
        while let Some(b) = line(
            &mut inner,
            (limits
                .max_input_bytes
                .min(limits.max_uncompressed_bytes)
                .saturating_sub(bytes))
            .min(limits.max_record_bytes as u64) as usize,
        )? {
            bytes = checked_add(
                bytes,
                b.len() as u64,
                limits.max_input_bytes.min(limits.max_uncompressed_bytes),
                "SAM source bytes",
            )?;
            hash.update(&b);
            if b.starts_with(b"@") {
                if raw
                    .len()
                    .checked_add(b.len())
                    .is_none_or(|v| v > limits.max_header_bytes)
                {
                    return Err(Error::Limit("SAM header bytes"));
                }
                raw.extend(b);
            } else {
                detect(&b)?;
                pending = Some(b);
                break;
            }
        }
        let header = parse_header(&raw, &limits)?;
        Ok(Self {
            inner,
            header,
            limits,
            pending,
            bytes,
            records: 0,
            emitted: 0,
            hash,
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
            source_bytes_accounted: self.bytes,
        }
    }
    pub fn completion(&self) -> Option<Completion> {
        (self.complete&&!self.failed).then(||Completion{format:Format::Sam,records:self.records,source_bytes:self.bytes,source_blake3:self.hash.finalize().to_hex().to_string(),warnings:vec!["Format integrity only; mappings, auxiliary annotations, uniqueness, variants and coverage are unverified imported assertions.".into()]})
    }
    pub fn next_record(&mut self) -> Result<Option<Record>> {
        if self.failed {
            return Err(Error::Poisoned);
        }
        let result = self.next();
        if result.is_err() {
            self.failed = true;
        }
        if matches!(&result, Ok(Some(_))) {
            self.emitted += 1;
        }
        result
    }
    fn next(&mut self) -> Result<Option<Record>> {
        if self.complete {
            return Ok(None);
        }
        let b = match self.pending.take() {
            Some(b) => Some(b),
            None => {
                let b = line(
                    &mut self.inner,
                    (self
                        .limits
                        .max_input_bytes
                        .min(self.limits.max_uncompressed_bytes)
                        .saturating_sub(self.bytes))
                    .min(self.limits.max_record_bytes as u64) as usize,
                )?;
                if let Some(b) = &b {
                    self.bytes = checked_add(
                        self.bytes,
                        b.len() as u64,
                        self.limits
                            .max_input_bytes
                            .min(self.limits.max_uncompressed_bytes),
                        "SAM source bytes",
                    )?;
                    self.hash.update(b);
                }
                b
            }
        };
        let Some(b) = b else {
            self.complete = true;
            return Ok(None);
        };
        self.records = checked_add(self.records, 1, self.limits.max_records, "records")?;
        Ok(Some(parse_record(b, &self.header, &self.limits)?))
    }
}
pub struct SamWriter<W> {
    inner: W,
    header_identity: [u8; 32],
    limits: Limits,
    line_boundary: bool,
    failed: bool,
    notices: Vec<ExportNotice>,
    records: u64,
    bytes: u64,
}
impl<W: Write> SamWriter<W> {
    pub fn new(mut inner: W, header: &Header, limits: Limits) -> Result<Self> {
        header.validate_policy(&limits)?;
        if header.raw.len() > limits.max_header_bytes
            || header.references.len() > limits.max_references
            || header.raw.len() as u64 > limits.max_input_bytes.min(limits.max_uncompressed_bytes)
        {
            return Err(Error::Limit("header bytes"));
        }
        if !header.sam_dictionary_matches {
            return Err(Error::Unsupported(
                "SAM export needs an explicit matching SQ text dictionary",
            ));
        }
        inner.write_all(&header.raw)?;
        Ok(Self {
            inner,
            header_identity: header.dictionary_identity,
            limits,
            line_boundary: header.raw.is_empty() || header.raw.ends_with(b"\n"),
            failed: false,
            notices: if header.bam_original.is_some() {
                vec![ExportNotice::TextFormattingNotPreserved]
            } else {
                vec![]
            },
            records: 0,
            bytes: header.raw.len() as u64,
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
        if !self.line_boundary {
            return Err(Error::Unsupported(
                "append after unterminated original SAM line",
            ));
        }
        validate_fields(&r.fields, &self.limits)?;
        self.records = checked_add(self.records, 1, self.limits.max_records, "export records")?;
        if r.original.len() > self.limits.max_record_bytes {
            return Err(Error::Limit("export record bytes"));
        }

        if r.dictionary != self.header_identity {
            return Err(Error::Invalid("record/header dictionary identity"));
        }
        if r.format == Format::Sam {
            self.bytes = checked_add(
                self.bytes,
                r.original.len() as u64,
                self.limits
                    .max_input_bytes
                    .min(self.limits.max_uncompressed_bytes),
                "SAM export bytes",
            )?;
            self.inner.write_all(&r.original)?;
            self.line_boundary = r.original.ends_with(b"\n");
            return Ok(vec![]);
        }
        let cap = (self
            .limits
            .max_input_bytes
            .min(self.limits.max_uncompressed_bytes)
            .saturating_sub(self.bytes))
        .min(self.limits.max_record_bytes as u64) as usize;
        let text = render_record(r, cap)?;
        if text.len() > self.limits.max_record_bytes {
            return Err(Error::Limit("SAM export record bytes"));
        }
        self.bytes = checked_add(
            self.bytes,
            text.len() as u64,
            self.limits
                .max_input_bytes
                .min(self.limits.max_uncompressed_bytes),
            "SAM export bytes",
        )?;
        self.inner.write_all(text.as_bytes())?;
        let mut n = vec![
            ExportNotice::TextFormattingNotPreserved,
            ExportNotice::AuxIntegerWidthNotPreserved,
        ];
        if r.long_cigar {
            n.push(ExportNotice::LongCigarTransportConverted);
        }
        Ok(n)
    }
    pub fn finish(mut self) -> Result<W> {
        if self.failed {
            return Err(Error::Poisoned);
        }
        self.inner.flush()?;
        Ok(self.inner)
    }
}
pub(crate) fn render_record(r: &Record, cap: usize) -> Result<String> {
    struct Counter {
        used: usize,
        cap: usize,
    }
    impl std::fmt::Write for Counter {
        fn write_str(&mut self, text: &str) -> std::fmt::Result {
            self.used = self
                .used
                .checked_add(text.len())
                .filter(|n| *n <= self.cap)
                .ok_or(std::fmt::Error)?;
            Ok(())
        }
    }
    let mut count = Counter { used: 0, cap };
    render_into(r, &mut count)?;
    let mut out = String::with_capacity(count.used);
    render_into(r, &mut out)?;
    if out.len() != count.used {
        return Err(Error::Invalid("internal SAM encoded-size disagreement"));
    }
    Ok(out)
}
fn render_into(r: &Record, s: &mut impl std::fmt::Write) -> Result<()> {
    let f = &r.fields;
    if [f.start, f.mate_start]
        .into_iter()
        .flatten()
        .any(|p| p >= i32::MAX as u32)
    {
        return Err(Error::Unsupported(
            "BAM position exceeds SAM 1-based coordinate range",
        ));
    }
    if f.qualities.as_deref() == Some([9u8].as_slice()) {
        return Err(Error::Unsupported(
            "solitary Q9 would become SAM missing-quality sentinel; explicit user conversion is required",
        ));
    }
    write!(
        s,
        "{}\t{}\t{}\t{}\t{}\t",
        f.name.as_deref().unwrap_or("*"),
        f.flags,
        f.reference.as_deref().unwrap_or("*"),
        f.start.map(|n| u64::from(n) + 1).unwrap_or(0),
        f.mapping_quality.unwrap_or(255)
    )?;
    if f.cigar.is_empty() {
        s.write_char('*')?;
    }
    for op in &f.cigar {
        write!(s, "{}{}", op.len, op.kind.symbol())?;
    }
    let mate = if f.mate_reference.is_some() && f.mate_reference == f.reference {
        "="
    } else {
        f.mate_reference.as_deref().unwrap_or("*")
    };
    write!(
        s,
        "\t{mate}\t{}\t{}\t{}\t",
        f.mate_start.map(|n| u64::from(n) + 1).unwrap_or(0),
        f.template_length,
        f.sequence.as_deref().unwrap_or("*")
    )?;
    if let Some(q) = &f.qualities {
        for c in q {
            s.write_char((c + 33) as char)?;
        }
    } else {
        s.write_char('*')?;
    }
    for a in &f.aux {
        if r.long_cigar && a.tag == *b"CG" {
            continue;
        }
        write!(s, "\t{}{}:", a.tag[0] as char, a.tag[1] as char)?;
        match &a.value {
            AuxValue::Char(c) => write!(s, "A:{}", *c as char)?,
            AuxValue::Integer(n) => write!(s, "i:{n}")?,
            AuxValue::FloatBits(v) => {
                let f = f32::from_bits(*v);
                if !f.is_finite() {
                    return Err(Error::Unsupported(
                        "non-finite BAM float cannot preserve IEEE payload in SAM",
                    ));
                }
                write!(s, "f:{f}")?;
            }
            AuxValue::String(v) => {
                s.write_str("Z:")?;
                s.write_str(text(v)?)?;
            }
            AuxValue::Hex(v) => {
                s.write_str("H:")?;
                s.write_str(text(v)?)?;
            }
            AuxValue::Array {
                element_type,
                values,
            } => {
                write!(s, "B:{}", *element_type as char)?;
                for v in values {
                    match *element_type {
                        b'c' => write!(s, ",{}", *v as u8 as i8),
                        b'C' => write!(s, ",{}", *v as u8),
                        b's' => write!(s, ",{}", *v as u16 as i16),
                        b'S' => write!(s, ",{}", *v as u16),
                        b'i' => write!(s, ",{}", *v as i32),
                        b'I' => write!(s, ",{v}"),
                        b'f' => {
                            let n = f32::from_bits(*v);
                            if !n.is_finite() {
                                return Err(Error::Unsupported("non-finite auxiliary array float"));
                            }
                            write!(s, ",{n}")
                        }
                        _ => return Err(Error::Unsupported("array type")),
                    }?;
                }
            }
        }
    }
    s.write_char('\n')?;
    Ok(())
}
