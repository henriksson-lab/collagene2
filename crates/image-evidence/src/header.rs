//! Container admission, not an alternate pixel codec. All offsets address bounded slices.
use crate::{Error, Format, Limits};
use std::collections::BTreeSet;

const MAX_RECORDS: usize = 4096;
const MAX_IFDS: usize = 32;
const MAX_VALUES: u64 = 65_536;
const MAX_METADATA: usize = 2 * 1024 * 1024;

pub(crate) struct Header {
    pub format: Format,
    pub width: u32,
    pub height: u32,
    pub orientation: u8,
    pub precision: u8,
    pub multiple: bool,
}
impl Header {
    pub fn read(bytes: &[u8], limits: Limits) -> Result<Self, Error> {
        if bytes.starts_with(b"\x89PNG\r\n\x1a\n") {
            Self::png(bytes, limits)
        } else if bytes.starts_with(&[0xff, 0xd8]) {
            Self::jpeg(bytes, limits, true)
        } else if bytes.starts_with(b"II") || bytes.starts_with(b"MM") {
            let scan = Tiff::read(bytes, limits)?;
            scan.image
                .ok_or(Error::Malformed("TIFF first IFD has no dimensions"))
        } else {
            Err(Error::Unsupported("expected JPEG, TIFF or PNG signature"))
        }
    }

    fn png(bytes: &[u8], limits: Limits) -> Result<Self, Error> {
        let mut result = Self {
            format: Format::Png,
            width: 0,
            height: 0,
            orientation: 1,
            precision: 0,
            multiple: false,
        };
        let mut offset = 8;
        let mut count = 0;
        let mut exif = false;
        let mut pixels = false;
        while offset < bytes.len() {
            count += 1;
            if count > MAX_RECORDS {
                return Err(Error::Limit("PNG chunk count"));
            }
            let length = be32(take(bytes, offset, 4)?) as usize;
            let kind = take(bytes, offset + 4, 4)?;
            let body = take(bytes, offset + 8, length)?;
            take(bytes, offset + 8 + length, 4)?; // CRC is verified by the pixel decoder.
            offset += length + 12;
            match kind {
                b"IHDR" => {
                    if count != 1 || body.len() != 13 {
                        return Err(Error::Malformed("PNG IHDR"));
                    }
                    result.width = be32(&body[..4]);
                    result.height = be32(&body[4..8]);
                    limits.dimensions(result.width, result.height)?;
                    result.precision = body[8];
                }
                b"IDAT" => {
                    pixels = true;
                }
                b"IEND" => {
                    if !body.is_empty() || !pixels || offset != bytes.len() {
                        return Err(Error::Malformed("PNG ending"));
                    }
                    return Ok(result);
                }
                b"eXIf" => {
                    if exif {
                        return Err(Error::Malformed("duplicate PNG EXIF"));
                    }
                    exif = true;
                    result.orientation = Tiff::read(body, limits)?.orientation;
                }
                b"acTL" => {
                    if body.len() != 8 || be32(&body[..4]) == 0 {
                        return Err(Error::Malformed("PNG animation control"));
                    }
                    result.multiple = be32(&body[..4]) > 1;
                }
                _ => {}
            }
        }
        Err(Error::Malformed("PNG missing IEND"))
    }

    fn jpeg(bytes: &[u8], limits: Limits, read_exif: bool) -> Result<Self, Error> {
        let mut result = None;
        let mut orientation = 1;
        let mut exif = false;
        let mut offset = 2;
        let mut records = 0;
        let mut entropy = false;
        let mut saw_scan = false;
        let mut scans = 0;
        while offset < bytes.len() {
            if entropy {
                // Scan compressed bytes without interpreting coefficients. Escaped FF and
                // restart markers remain entropy data; the next real marker ends the scan.
                while offset < bytes.len() && bytes[offset] != 0xff {
                    offset += 1;
                }
                if offset == bytes.len() {
                    break;
                }
            } else if bytes[offset] != 0xff {
                return Err(Error::Malformed("JPEG marker prefix"));
            }
            while offset < bytes.len() && bytes[offset] == 0xff {
                offset += 1;
            }
            let marker = *bytes
                .get(offset)
                .ok_or(Error::Malformed("truncated JPEG marker"))?;
            offset += 1;
            if entropy && (marker == 0 || (0xd0..=0xd7).contains(&marker)) {
                continue;
            }
            entropy = false;
            records += 1;
            if records > MAX_RECORDS {
                return Err(Error::Limit("JPEG marker count"));
            }
            if marker == 0xd9 {
                if !saw_scan || offset != bytes.len() {
                    return Err(Error::Malformed("JPEG ending"));
                }
                let mut header: Self = result.ok_or(Error::Malformed("JPEG missing frame"))?;
                header.orientation = orientation;
                return Ok(header);
            }
            if marker == 0 || marker == 0xd8 || (0xd0..=0xd7).contains(&marker) {
                return Err(Error::Malformed("unexpected JPEG marker"));
            }
            if marker == 1 {
                continue;
            } // TEM has no length.
            let length_bytes = take(bytes, offset, 2)?;
            let length = usize::from(u16::from_be_bytes([length_bytes[0], length_bytes[1]]));
            if length < 2 {
                return Err(Error::Malformed("JPEG segment length"));
            }
            let body = take(bytes, offset + 2, length - 2)?;
            offset += length;
            match marker {
                0xc0..=0xc3 | 0xc5..=0xc7 | 0xc9..=0xcb | 0xcd..=0xcf => {
                    if result.is_some() {
                        return Err(Error::Unsupported("multiple JPEG frame headers"));
                    }
                    if !matches!(marker, 0xc0..=0xc2) {
                        return Err(Error::Unsupported("JPEG coding process"));
                    }
                    if body.len() < 6 || body[0] != 8 {
                        return Err(Error::Unsupported("JPEG sample precision"));
                    }
                    let height = u32::from(u16::from_be_bytes([body[1], body[2]]));
                    let width = u32::from(u16::from_be_bytes([body[3], body[4]]));
                    limits.dimensions(width, height)?;
                    let channels = usize::from(body[5]);
                    if channels == 0 || channels > 4 || body.len() != 6 + channels * 3 {
                        return Err(Error::Unsupported("JPEG component count"));
                    }
                    let mut blocks = 0u32;
                    for component in body[6..].chunks_exact(3) {
                        let horizontal = component[1] >> 4;
                        let vertical = component[1] & 15;
                        if !(1..=4).contains(&horizontal) || !(1..=4).contains(&vertical) {
                            return Err(Error::Malformed("JPEG sampling factor"));
                        }
                        blocks += u32::from(horizontal) * u32::from(vertical);
                    }
                    if blocks > 10 {
                        return Err(Error::Unsupported("JPEG MCU component bound"));
                    }
                    result = Some(Self {
                        format: Format::Jpeg,
                        width,
                        height,
                        orientation: 1,
                        precision: 8,
                        multiple: false,
                    });
                }
                0xe1 if read_exif && body.starts_with(b"Exif\0\0") => {
                    if exif {
                        return Err(Error::Malformed("duplicate JPEG EXIF"));
                    }
                    exif = true;
                    orientation = Tiff::read(&body[6..], limits)?.orientation;
                }
                0xda => {
                    if result.is_none() {
                        return Err(Error::Malformed("JPEG scan before frame"));
                    }
                    scans += 1;
                    if scans > 64 {
                        return Err(Error::Limit("JPEG scan count"));
                    }
                    if body.is_empty()
                        || body[0] == 0
                        || body[0] > 4
                        || body.len() != 1 + usize::from(body[0]) * 2 + 3
                    {
                        return Err(Error::Malformed("JPEG scan header"));
                    }
                    entropy = true;
                    saw_scan = true;
                }
                0xdc => return Err(Error::Unsupported("JPEG delayed height")),
                _ => {}
            }
        }
        Err(Error::Malformed("JPEG missing EOI"))
    }
}

fn take(bytes: &[u8], offset: usize, length: usize) -> Result<&[u8], Error> {
    let end = offset
        .checked_add(length)
        .ok_or(Error::Malformed("container offset overflow"))?;
    bytes
        .get(offset..end)
        .ok_or(Error::Malformed("container offset outside input"))
}
fn be32(bytes: &[u8]) -> u32 {
    u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])
}

/// Classic and BigTIFF share one bounded IFD traversal. Metadata pointer edges are
/// traversed without trusting declared lengths or allocating their claimed counts.
struct Tiff<'a> {
    bytes: &'a [u8],
    little: bool,
    big: bool,
    orientation: u8,
    image: Option<Header>,
}
#[derive(Clone, Copy)]
struct Field<'a> {
    tag: u16,
    kind: u16,
    count: usize,
    bytes: &'a [u8],
}
impl<'a> Tiff<'a> {
    fn read(bytes: &'a [u8], limits: Limits) -> Result<Self, Error> {
        let signature = take(bytes, 0, 4)?;
        let little = match &signature[..2] {
            b"II" => true,
            b"MM" => false,
            _ => return Err(Error::Malformed("EXIF/TIFF byte order")),
        };
        let mut this = Self {
            bytes,
            little,
            big: false,
            orientation: 1,
            image: None,
        };
        let magic = this.number(&signature[2..4]);
        let first = match magic {
            42 => this.number(take(bytes, 4, 4)?),
            43 => {
                this.big = true;
                if this.number(take(bytes, 4, 2)?) != 8 || this.number(take(bytes, 6, 2)?) != 0 {
                    return Err(Error::Unsupported("BigTIFF offset layout"));
                }
                this.number(take(bytes, 8, 8)?)
            }
            _ => return Err(Error::Unsupported("TIFF version")),
        };
        if first == 0 {
            return Err(Error::Malformed("TIFF missing first IFD"));
        }
        let mut pending = vec![(first, true, true)];
        let mut visited = BTreeSet::new();
        let mut records = 0usize;
        let mut metadata = 0usize;
        let mut compressed_extents = 0usize;
        let mut seen_orientation = false;
        while let Some((offset, primary, is_first)) = pending.pop() {
            if visited.len() >= MAX_IFDS {
                return Err(Error::Limit("TIFF IFD count"));
            }
            if !visited.insert(offset) {
                return Err(Error::Malformed("TIFF repeated/cyclic IFD"));
            }
            let offset =
                usize::try_from(offset).map_err(|_| Error::Malformed("TIFF offset size"))?;
            let (count_width, entry_width, inline_width) =
                if this.big { (8, 20, 8) } else { (2, 12, 4) };
            let count = this.number(take(bytes, offset, count_width)?);
            if count > MAX_RECORDS as u64 {
                return Err(Error::Limit("TIFF IFD entry count"));
            }
            records += count as usize;
            if records > MAX_RECORDS {
                return Err(Error::Limit("TIFF total entry count"));
            }
            let entries_start = offset
                .checked_add(count_width)
                .ok_or(Error::Malformed("TIFF offset overflow"))?;
            let entries = take(bytes, entries_start, count as usize * entry_width)?;
            let mut fields = Vec::with_capacity(count as usize);
            let mut tags = BTreeSet::new();
            for entry in entries.chunks_exact(entry_width) {
                let tag = this.number(&entry[..2]) as u16;
                if !tags.insert(tag) {
                    return Err(Error::Malformed("duplicate TIFF tag"));
                }
                let kind = this.number(&entry[2..4]) as u16;
                let count = this.number(&entry[4..entry_width - inline_width]);
                if count > MAX_VALUES {
                    return Err(Error::Limit("TIFF tag value count"));
                }
                let type_width = match kind {
                    1 | 2 | 6 | 7 => 1,
                    3 | 8 => 2,
                    4 | 9 | 11 | 13 => 4,
                    5 | 10 | 12 | 16 | 17 | 18 => 8,
                    _ => return Err(Error::Unsupported("TIFF field type")),
                };
                let length = count as usize * type_width;
                metadata = metadata
                    .checked_add(length)
                    .ok_or(Error::Limit("TIFF metadata length"))?;
                if metadata > MAX_METADATA {
                    return Err(Error::Limit("TIFF metadata bytes"));
                }
                let inline = &entry[entry_width - inline_width..];
                let value = if length <= inline_width {
                    &inline[..length]
                } else {
                    let pointer = usize::try_from(this.number(inline))
                        .map_err(|_| Error::Malformed("TIFF value offset size"))?;
                    take(bytes, pointer, length)?
                };
                fields.push(Field {
                    tag,
                    kind,
                    count: count as usize,
                    bytes: value,
                });
            }
            let next = this.number(take(bytes, entries_start + entries.len(), inline_width)?);
            if next != 0 {
                pending.push((next, primary, false));
            }
            for field in &fields {
                match field.tag {
                    330 | 34665 | 34853 | 40965 => {
                        for index in 0..field.count {
                            let pointer = this.unsigned(*field, index)?;
                            if pointer != 0 {
                                pending.push((pointer, false, false));
                            }
                            if pending.len() > MAX_IFDS {
                                return Err(Error::Limit("TIFF IFD pointer count"));
                            }
                        }
                    }
                    274 if is_first => {
                        if seen_orientation || field.count != 1 {
                            return Err(Error::Malformed("TIFF orientation count"));
                        }
                        seen_orientation = true;
                        let orientation = this.unsigned(*field, 0)?;
                        if !(1..=8).contains(&orientation) {
                            return Err(Error::Malformed("TIFF orientation value"));
                        }
                        this.orientation = orientation as u8;
                    }
                    _ => {}
                }
            }
            // Check all primary image sizes, without decoding secondary images.
            let width = this.scalar(&fields, 256)?;
            let height = this.scalar(&fields, 257)?;
            if primary && (width.is_some() || height.is_some()) {
                let width = u32::try_from(width.ok_or(Error::Malformed("TIFF missing width"))?)
                    .map_err(|_| Error::Limit("TIFF width"))?;
                let height = u32::try_from(height.ok_or(Error::Malformed("TIFF missing height"))?)
                    .map_err(|_| Error::Limit("TIFF height"))?;
                limits.dimensions(width, height)?;
                if is_first {
                    let mut precision = 1;
                    if let Some(bits) = fields.iter().find(|f| f.tag == 258) {
                        if bits.count > 4 {
                            return Err(Error::Unsupported("TIFF sample count"));
                        }
                        for index in 0..bits.count {
                            let value = this.unsigned(*bits, index)?;
                            if value == 0 || value > 32 {
                                return Err(Error::Unsupported("TIFF sample precision"));
                            }
                            precision = precision.max(value as u8);
                        }
                    }
                    if this.scalar(&fields, 277)?.unwrap_or(1) > 4 {
                        return Err(Error::Unsupported("TIFF samples per pixel"));
                    }
                    this.image = Some(Header {
                        format: Format::Tiff,
                        width,
                        height,
                        orientation: this.orientation,
                        precision,
                        multiple: next != 0 || fields.iter().any(|f| f.tag == 330),
                    });
                }
            }
            // Strip/tile compressed extents must address this file, with equal counts.
            for (offset_tag, size_tag) in [(273, 279), (324, 325)] {
                let offsets = fields.iter().find(|f| f.tag == offset_tag);
                let sizes = fields.iter().find(|f| f.tag == size_tag);
                if offsets.is_some() != sizes.is_some() {
                    return Err(Error::Malformed("TIFF unpaired block offsets/counts"));
                }
                if let (Some(offsets), Some(sizes)) = (offsets, sizes) {
                    if offsets.count != sizes.count {
                        return Err(Error::Malformed("TIFF block count mismatch"));
                    }
                    for index in 0..offsets.count {
                        let start = usize::try_from(this.unsigned(*offsets, index)?)
                            .map_err(|_| Error::Malformed("TIFF block offset size"))?;
                        let length = usize::try_from(this.unsigned(*sizes, index)?)
                            .map_err(|_| Error::Malformed("TIFF block length size"))?;
                        let block = take(bytes, start, length)?;
                        compressed_extents = compressed_extents
                            .checked_add(length)
                            .ok_or(Error::Limit("TIFF compressed extent sum"))?;
                        if compressed_extents > limits.input_bytes * 2 {
                            return Err(Error::Limit("TIFF repeated compressed extents"));
                        }
                        // TIFF's JPEG codec allocates from the embedded SOF, independently
                        // of the outer TIFF dimensions. Admit that header before invoking it.
                        if this.scalar(&fields, 259)? == Some(7) {
                            let mut joined = Vec::new();
                            let jpeg = if let Some(tables) = fields.iter().find(|f| f.tag == 347) {
                                if tables.bytes.len() < 4
                                    || !tables.bytes.starts_with(&[0xff, 0xd8])
                                    || !tables.bytes.ends_with(&[0xff, 0xd9])
                                    || !block.starts_with(&[0xff, 0xd8])
                                {
                                    return Err(Error::Malformed("TIFF JPEG tables/block framing"));
                                }
                                // Shared tables are replayed by the codec for every block;
                                // charge that work too, not just the physical strip bytes.
                                compressed_extents = compressed_extents
                                    .checked_add(tables.bytes.len() - 2)
                                    .ok_or(Error::Limit("TIFF JPEG table replay sum"))?;
                                if compressed_extents > limits.input_bytes * 2 {
                                    return Err(Error::Limit("TIFF JPEG table replay bytes"));
                                }
                                joined.extend_from_slice(&tables.bytes[..tables.bytes.len() - 2]);
                                joined.extend_from_slice(&block[2..]);
                                joined.as_slice()
                            } else {
                                block
                            };
                            let embedded = Header::jpeg(jpeg, limits, false)?;
                            let source_width = this
                                .scalar(&fields, 256)?
                                .ok_or(Error::Malformed("TIFF JPEG missing width"))?;
                            let source_height = this
                                .scalar(&fields, 257)?
                                .ok_or(Error::Malformed("TIFF JPEG missing height"))?;
                            let width = this.scalar(&fields, 322)?.unwrap_or(source_width);
                            let height = this.scalar(&fields, 323)?.unwrap_or(
                                this.scalar(&fields, 278)?
                                    .unwrap_or(source_height)
                                    .min(source_height),
                            );
                            if u64::from(embedded.width) > width.saturating_add(32)
                                || u64::from(embedded.height) > height.saturating_add(32)
                            {
                                return Err(Error::Limit(
                                    "TIFF JPEG block/header dimension mismatch",
                                ));
                            }
                        }
                    }
                }
            }
            if let Some(tile_width) = this.scalar(&fields, 322)? {
                let tile_height = this
                    .scalar(&fields, 323)?
                    .ok_or(Error::Malformed("TIFF tile height"))?;
                if tile_width == 0
                    || tile_height == 0
                    || tile_width > u64::from(limits.dimension)
                    || tile_height > u64::from(limits.dimension)
                    || tile_width * tile_height > limits.source_pixels
                    || tile_width * tile_height * 16 > limits.decoded_bytes
                {
                    return Err(Error::Limit("TIFF tile intermediate dimensions"));
                }
            }
        }
        Ok(this)
    }
    fn number(&self, bytes: &[u8]) -> u64 {
        if self.little {
            bytes.iter().rev().fold(0, |n, b| (n << 8) | u64::from(*b))
        } else {
            bytes.iter().fold(0, |n, b| (n << 8) | u64::from(*b))
        }
    }
    fn unsigned(&self, field: Field<'_>, index: usize) -> Result<u64, Error> {
        let width = match field.kind {
            1 => 1,
            3 => 2,
            4 | 13 => 4,
            16 | 18 => 8,
            _ => return Err(Error::Unsupported("TIFF unsigned field representation")),
        };
        if index >= field.count {
            return Err(Error::Malformed("TIFF field index"));
        }
        Ok(self.number(take(field.bytes, index * width, width)?))
    }
    fn scalar(&self, fields: &[Field<'_>], tag: u16) -> Result<Option<u64>, Error> {
        fields
            .iter()
            .find(|f| f.tag == tag)
            .map(|field| {
                if field.count != 1 {
                    return Err(Error::Malformed("TIFF scalar count"));
                }
                self.unsigned(*field, 0)
            })
            .transpose()
    }
}
