use collagene_image_evidence::{Error, Evidence, Format, Limits, Warning};

fn prepare(bytes: Vec<u8>) -> Evidence {
    Evidence::prepare(bytes, Limits::default()).unwrap()
}
fn pixels(evidence: &Evidence) -> Vec<u8> {
    image::load_from_memory_with_format(evidence.preview_png(), image::ImageFormat::Png)
        .unwrap()
        .to_rgba8()
        .pixels()
        .map(|p| p[0])
        .collect()
}

// Independent, uncompressed TIFF fixture writer. It does not use the decoder's encoder.
fn tiff(big: bool, little: bool, orientation: u16, bits: u16) -> Vec<u8> {
    fn number(value: u64, size: usize, little: bool) -> Vec<u8> {
        let bytes = if little {
            value.to_le_bytes()
        } else {
            value.to_be_bytes()
        };
        if little {
            bytes[..size].to_vec()
        } else {
            bytes[8 - size..].to_vec()
        }
    }
    let (header, count_width, entry_width, inline) =
        if big { (16, 8, 20, 8) } else { (8, 2, 12, 4) };
    let entries = 11;
    let metadata_end = header + count_width + entries * entry_width + inline;
    let samples = if bits == 8 { 3 } else { 1 };
    let bits_offset = metadata_end;
    let bits_length = if bits == 8 && !big { 6 } else { 0 };
    let pixel_offset = metadata_end + bits_length;
    let values: Vec<u8> = if bits == 8 {
        [10, 20, 30, 40, 50, 60]
            .into_iter()
            .flat_map(|x| [x, 0, 0])
            .collect()
    } else {
        [0, 0x1234, 0x5678, 0x9abc, 0xdef0, 0xffff]
            .into_iter()
            .flat_map(|x| number(x, 2, little))
            .collect()
    };
    let mut out = if little {
        b"II".to_vec()
    } else {
        b"MM".to_vec()
    };
    out.extend(number(if big { 43 } else { 42 }, 2, little));
    if big {
        out.extend(number(8, 2, little));
        out.extend(number(0, 2, little));
    }
    out.extend(number(header as u64, inline, little));
    out.extend(number(entries as u64, count_width, little));
    for (tag, kind, count, value) in [
        (256, 4, 1, 2),
        (257, 4, 1, 3),
        (
            258,
            3,
            samples,
            if bits_length > 0 {
                bits_offset as u64
            } else {
                bits as u64
            },
        ),
        (259, 3, 1, 1),
        (262, 3, 1, if samples == 3 { 2 } else { 1 }),
        (273, 4, 1, pixel_offset as u64),
        (274, 3, 1, orientation as u64),
        (277, 3, 1, samples),
        (278, 4, 1, 3),
        (279, 4, 1, values.len() as u64),
        (284, 3, 1, 1),
    ] {
        out.extend(number(tag, 2, little));
        out.extend(number(kind, 2, little));
        out.extend(number(count, if big { 8 } else { 4 }, little));
        if tag == 258 && big && samples == 3 {
            for _ in 0..3 {
                out.extend(number(8, 2, little));
            }
            out.resize(out.len() + inline - 6, 0);
        } else if kind == 3 && !(tag == 258 && bits_length > 0) {
            out.extend(number(value, 2, little));
            out.resize(out.len() + inline - 2, 0);
        } else {
            out.extend(number(value, 4, little));
            out.resize(out.len() + inline - 4, 0);
        }
    }
    out.extend(number(0, inline, little));
    if bits_length > 0 {
        for _ in 0..3 {
            out.extend(number(8, 2, little));
        }
    }
    assert_eq!(out.len(), pixel_offset);
    out.extend(values);
    out
}

fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = !0u32;
    for byte in bytes {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xedb88320 & (0u32.wrapping_sub(crc & 1)));
        }
    }
    !crc
}
fn chunk(out: &mut Vec<u8>, kind: &[u8; 4], body: &[u8]) {
    out.extend((body.len() as u32).to_be_bytes());
    let start = out.len();
    out.extend(kind);
    out.extend(body);
    let crc = crc32(&out[start..]);
    out.extend(crc.to_be_bytes());
}
fn png(exif: Option<&[u8]>) -> Vec<u8> {
    let mut out = b"\x89PNG\r\n\x1a\n".to_vec();
    let mut ihdr = vec![];
    ihdr.extend(2u32.to_be_bytes());
    ihdr.extend(3u32.to_be_bytes());
    ihdr.extend([8, 2, 0, 0, 0]);
    chunk(&mut out, b"IHDR", &ihdr);
    if let Some(exif) = exif {
        chunk(&mut out, b"eXIf", exif);
    }
    let raw: Vec<u8> = [
        [0, 10, 0, 0, 20, 0, 0],
        [0, 30, 0, 0, 40, 0, 0],
        [0, 50, 0, 0, 60, 0, 0],
    ]
    .into_iter()
    .flatten()
    .collect();
    let length = raw.len() as u16;
    let mut zlib = vec![0x78, 0x01, 1];
    zlib.extend(length.to_le_bytes());
    zlib.extend((!length).to_le_bytes());
    zlib.extend(&raw);
    let (mut a, mut b) = (1u32, 0u32);
    for value in raw {
        a = (a + u32::from(value)) % 65521;
        b = (b + a) % 65521;
    }
    zlib.extend(((b << 16) | a).to_be_bytes());
    chunk(&mut out, b"IDAT", &zlib);
    chunk(&mut out, b"IEND", &[]);
    out
}
// Hand-built baseline JPEG: one DC=0/EOB block, single-symbol Huffman tables.
fn jpeg(exif: Option<&[u8]>) -> Vec<u8> {
    fn segment(out: &mut Vec<u8>, marker: u8, body: &[u8]) {
        out.extend([0xff, marker]);
        out.extend(((body.len() + 2) as u16).to_be_bytes());
        out.extend(body);
    }
    let mut out = vec![0xff, 0xd8];
    if let Some(exif) = exif {
        let mut body = b"Exif\0\0".to_vec();
        body.extend(exif);
        segment(&mut out, 0xe1, &body);
    }
    let mut quantization = vec![0];
    quantization.extend([1; 64]);
    segment(&mut out, 0xdb, &quantization);
    segment(&mut out, 0xc0, &[8, 0, 8, 0, 8, 1, 1, 0x11, 0]);
    let mut huffman = vec![];
    for table in [0, 0x10] {
        huffman.extend([table, 1]);
        huffman.extend([0; 15]);
        huffman.push(0);
    }
    segment(&mut out, 0xc4, &huffman);
    segment(&mut out, 0xda, &[1, 1, 0, 0, 63, 0]);
    out.extend([0x3f, 0xff, 0xd9]);
    out
}

#[test]
fn classic_and_big_tiff_both_byte_orders_preserve_exact_original() {
    for big in [false, true] {
        for little in [false, true] {
            let original = tiff(big, little, 1, 8);
            let result = prepare(original.clone());
            assert_eq!(result.format(), Format::Tiff);
            assert_eq!(result.source_dimensions(), (2, 3));
            assert_eq!(result.preview_dimensions(), (2, 3));
            assert_eq!(pixels(&result), [10, 20, 30, 40, 50, 60]);
            assert_eq!(result.content_hash(), blake3::hash(&original).as_bytes());
            assert_eq!(result.original_bytes(), original);
            assert_eq!(result.into_original(), original);
        }
    }
}
#[test]
fn all_eight_orientations_have_independent_pixel_oracles() {
    let oracles: [&[u8]; 8] = [
        &[10, 20, 30, 40, 50, 60],
        &[20, 10, 40, 30, 60, 50],
        &[60, 50, 40, 30, 20, 10],
        &[50, 60, 30, 40, 10, 20],
        &[10, 30, 50, 20, 40, 60],
        &[50, 30, 10, 60, 40, 20],
        &[60, 40, 20, 50, 30, 10],
        &[20, 40, 60, 10, 30, 50],
    ];
    for orientation in 1..=8 {
        let result = prepare(tiff(false, true, orientation, 8));
        assert_eq!(result.source_orientation(), orientation as u8);
        assert_eq!(
            result.preview_dimensions(),
            if orientation < 5 { (2, 3) } else { (3, 2) }
        );
        assert_eq!(pixels(&result), oracles[orientation as usize - 1]);
    }
}
#[test]
fn independent_png_and_jpeg_fixtures_decode_and_reexport() {
    for (original, format, dimensions) in [
        (png(None), Format::Png, (2, 3)),
        (jpeg(None), Format::Jpeg, (8, 8)),
    ] {
        let result = prepare(original.clone());
        assert_eq!(result.format(), format);
        assert_eq!(result.source_dimensions(), dimensions);
        assert_eq!(result.original_bytes(), original);
        assert!(result.warnings().contains(&Warning::DocumentaryPreview));
        assert!(result.warnings().contains(&Warning::MetadataNotTransferred));
        assert!(result.warnings().contains(&Warning::UncalibratedDisplay));
        if format == Format::Png {
            assert_eq!(pixels(&result), [10, 20, 30, 40, 50, 60]);
        } else {
            assert!(pixels(&result).iter().all(|p| *p == 128));
        }
    }
}
#[test]
fn png_and_jpeg_exif_is_applied_once() {
    let exif = tiff(false, true, 6, 8);
    let result = prepare(png(Some(&exif)));
    assert_eq!(result.preview_dimensions(), (3, 2));
    assert_eq!(pixels(&result), [50, 30, 10, 60, 40, 20]);
    let result = prepare(jpeg(Some(&exif)));
    assert_eq!(result.source_orientation(), 6);
    // Rendered PNG carries no EXIF orientation, so re-import cannot rotate it twice.
    let again = prepare(result.preview_png().to_vec());
    assert_eq!(again.source_orientation(), 1);
    assert_eq!(again.preview_dimensions(), result.preview_dimensions());
}
#[test]
fn high_bit_original_is_retained_and_preview_precision_is_explicit() {
    for big in [false, true] {
        for little in [false, true] {
            let original = tiff(big, little, 1, 16);
            let result = prepare(original.clone());
            assert_eq!(result.original_bytes(), original);
            assert!(
                result
                    .warnings()
                    .contains(&Warning::PrecisionReducedForPreview)
            );
            assert_eq!(pixels(&result).last(), Some(&255));
        }
    }
}
#[test]
fn multipage_keeps_original_and_warns_without_claiming_other_pages_decoded() {
    let mut original = tiff(false, true, 1, 8);
    // Append an independent second IFD and adjust its data/array offsets.
    let mut second = tiff(false, true, 1, 8);
    let shift = original.len() as u32;
    for index in [2usize, 5] {
        let location = 8 + 2 + index * 12 + 8;
        let value = u32::from_le_bytes(second[location..location + 4].try_into().unwrap()) + shift;
        second[location..location + 4].copy_from_slice(&value.to_le_bytes());
    }
    let next = 8 + 2 + 11 * 12;
    original[next..next + 4].copy_from_slice(&(shift + 8).to_le_bytes());
    original.extend(second);
    let result = prepare(original.clone());
    assert_eq!(result.original_bytes(), original);
    assert!(
        result
            .warnings()
            .contains(&Warning::AdditionalImagesNotPreviewed)
    );
    assert_eq!(pixels(&result), [10, 20, 30, 40, 50, 60]);
}
#[test]
fn limits_reject_before_pixel_decoding_and_preview_is_bounded() {
    let input = png(None);
    let mut limits = Limits::default();
    limits.input_bytes = input.len() - 1;
    assert!(matches!(
        Evidence::prepare(input.clone(), limits),
        Err(Error::Limit(_))
    ));
    let mut limits = Limits::default();
    limits.source_pixels = 5;
    assert!(matches!(
        Evidence::prepare(input.clone(), limits),
        Err(Error::Limit(_))
    ));
    let mut limits = Limits::default();
    limits.decoded_bytes = 1;
    assert!(matches!(
        Evidence::prepare(input.clone(), limits),
        Err(Error::Limit(_))
    ));
    let mut limits = Limits::default();
    limits.working_bytes = 1;
    assert!(matches!(
        Evidence::prepare(input.clone(), limits),
        Err(Error::Limit(_))
    ));
    let mut limits = Limits::default();
    limits.preview_pixels = 2;
    limits.preview_dimension = 2;
    let result = Evidence::prepare(input.clone(), limits).unwrap();
    assert!(
        u64::from(result.preview_dimensions().0) * u64::from(result.preview_dimensions().1) <= 2
    );
    assert!(result.warnings().contains(&Warning::DownsampledForPreview));
    let mut limits = Limits::default();
    limits.preview_bytes = 1;
    assert!(Evidence::prepare(input, limits).is_err());
}
#[test]
fn unsupported_limits_and_signatures_fail_explicitly() {
    let mut limits = Limits::default();
    limits.input_bytes += 1;
    assert!(matches!(
        Evidence::prepare(png(None), limits),
        Err(Error::Limit(_))
    ));
    assert!(matches!(
        Evidence::prepare(b"GIF89a".to_vec(), Limits::default()),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        Evidence::prepare(b"II\x2c\0\0\0\0\0".to_vec(), Limits::default()),
        Err(Error::Unsupported(_))
    ));
}
#[test]
fn malformed_tiff_ifds_counts_offsets_cycles_and_orientation_are_rejected() {
    let valid = tiff(false, true, 1, 8);
    let mut cases = vec![];
    let mut count = valid.clone();
    count[8..10].copy_from_slice(&u16::MAX.to_le_bytes());
    cases.push(count);
    let mut values = valid.clone();
    values[8 + 2 + 2 * 12 + 4..8 + 2 + 2 * 12 + 8].copy_from_slice(&u32::MAX.to_le_bytes());
    cases.push(values);
    let mut offset = valid.clone();
    offset[8 + 2 + 5 * 12 + 8..8 + 2 + 5 * 12 + 12].copy_from_slice(&u32::MAX.to_le_bytes());
    cases.push(offset);
    let mut cycle = valid.clone();
    let next = 8 + 2 + 11 * 12;
    cycle[next..next + 4].copy_from_slice(&8u32.to_le_bytes());
    cases.push(cycle);
    let mut duplicate = valid.clone();
    duplicate[8 + 2 + 12..8 + 2 + 14].copy_from_slice(&256u16.to_le_bytes());
    cases.push(duplicate);
    let mut giant_dimension = valid.clone();
    giant_dimension[8 + 2 + 8..8 + 2 + 12].copy_from_slice(&u32::MAX.to_le_bytes());
    cases.push(giant_dimension);
    cases.push(tiff(false, true, 0, 8));
    cases.push(tiff(false, true, 9, 8));
    let mut big_count = tiff(true, true, 1, 8);
    big_count[16..24].copy_from_slice(&u64::MAX.to_le_bytes());
    cases.push(big_count);
    for case in cases {
        assert!(Evidence::prepare(case, Limits::default()).is_err());
    }
}
#[test]
fn every_truncated_container_prefix_is_rejected_without_panicking() {
    for valid in [
        png(None),
        jpeg(None),
        tiff(false, true, 1, 8),
        tiff(true, false, 1, 8),
    ] {
        for end in 0..valid.len() {
            assert!(
                Evidence::prepare(valid[..end].to_vec(), Limits::default()).is_err(),
                "accepted prefix {end}/{}",
                valid.len()
            );
        }
    }
}
#[test]
fn huge_png_dimensions_and_jpeg_sampling_are_rejected_before_decode() {
    let mut huge = b"\x89PNG\r\n\x1a\n".to_vec();
    let mut ihdr = vec![];
    ihdr.extend(u32::MAX.to_be_bytes());
    ihdr.extend(u32::MAX.to_be_bytes());
    ihdr.extend([16, 6, 0, 0, 0]);
    chunk(&mut huge, b"IHDR", &ihdr);
    assert!(matches!(
        Evidence::prepare(huge, Limits::default()),
        Err(Error::Limit(_))
    ));
    let mut invalid = jpeg(None);
    let sof = invalid.windows(2).position(|w| w == [0xff, 0xc0]).unwrap();
    invalid[sof + 4 + 7] = 0xff;
    assert!(Evidence::prepare(invalid, Limits::default()).is_err());
}
#[test]
fn png_crc_corruption_and_jpeg_missing_end_are_rejected() {
    let mut corrupt = png(None);
    corrupt[29] ^= 1;
    assert!(Evidence::prepare(corrupt, Limits::default()).is_err());
    let mut truncated = jpeg(None);
    truncated.truncate(truncated.len() - 2);
    assert!(Evidence::prepare(truncated, Limits::default()).is_err());
}

#[test]
fn tiff_embedded_jpeg_cannot_override_outer_resource_admission() {
    let mut input = tiff(false, true, 1, 16);
    // Convert the outer single-channel fixture to JPEG compression; embed an
    // attacker-controlled SOF whose dimensions exceed the outer 2x3 image.
    let compression = 8 + 2 + 3 * 12 + 8;
    input[compression..compression + 2].copy_from_slice(&7u16.to_le_bytes());
    let bits = 8 + 2 + 2 * 12 + 8;
    input[bits..bits + 2].copy_from_slice(&8u16.to_le_bytes());
    let mut embedded = jpeg(None);
    let sof = embedded.windows(2).position(|w| w == [0xff, 0xc0]).unwrap();
    embedded[sof + 4 + 1..sof + 4 + 5].copy_from_slice(&[0xff; 4]);
    let offset = 8 + 2 + 11 * 12 + 4;
    input.truncate(offset);
    let count = 8 + 2 + 9 * 12 + 8;
    input[count..count + 4].copy_from_slice(&(embedded.len() as u32).to_le_bytes());
    input.extend(embedded);
    assert!(matches!(
        Evidence::prepare(input, Limits::default()),
        Err(Error::Limit(_))
    ));
}

#[test]
fn valid_jpeg_compressed_tiff_uses_admitted_embedded_dimensions() {
    let mut input = tiff(false, true, 1, 16);
    for (index, value, width) in [
        (0usize, 8u32, 4usize),
        (1, 8, 4),
        (2, 8, 2),
        (3, 7, 2),
        (8, 8, 4),
    ] {
        let field = 8 + 2 + index * 12 + 8;
        input[field..field + width].copy_from_slice(&value.to_le_bytes()[..width]);
    }
    let embedded = jpeg(None);
    let offset = 8 + 2 + 11 * 12 + 4;
    input.truncate(offset);
    let count = 8 + 2 + 9 * 12 + 8;
    input[count..count + 4].copy_from_slice(&(embedded.len() as u32).to_le_bytes());
    input.extend(embedded);
    let result = prepare(input.clone());
    assert_eq!(result.original_bytes(), input);
    assert_eq!(result.source_dimensions(), (8, 8));
    assert!(pixels(&result).iter().all(|p| *p == 128));
}
