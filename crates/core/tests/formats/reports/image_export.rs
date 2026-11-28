//! Native static PNG pixels are checked independently of SVG string identity.
use collagene_core::{
    output::{
        gel::{self, Calibration, LaneInput, Preparation},
        gel_export, map_export,
        raster::{Format, Image, MAX_PIXELS, PngImage, PngOptions},
    },
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Feature, Location, Molecule, Topology};
use std::io::Cursor;

fn decode(image: &PngImage) -> Vec<u8> {
    image.validate_binding(&image.context_revision).unwrap();
    let decoder = png::Decoder::new_with_limits(
        Cursor::new(&image.bytes),
        png::Limits {
            bytes: 32 * 1024 * 1024,
        },
    );
    let mut reader = decoder.read_info().unwrap();
    assert_eq!(reader.info().width, image.width);
    assert_eq!(reader.info().height, image.height);
    assert_eq!(reader.info().color_type, png::ColorType::Rgba);
    assert_eq!(reader.info().bit_depth, png::BitDepth::Eight);
    let size = reader.output_buffer_size().unwrap();
    assert_eq!(size, image.width as usize * image.height as usize * 4);
    assert!(size as u64 <= MAX_PIXELS * 4);
    let mut pixels = vec![0; size];
    reader.next_frame(&mut pixels).unwrap();
    reader.finish().unwrap();
    pixels
}
fn pixel(pixels: &[u8], width: u32, x: usize, y: usize) -> &[u8] {
    let index = (y * width as usize + x) * 4;
    &pixels[index..index + 4]
}
fn circular_map(format: Format) -> map_export::Request {
    let mut source = Molecule::dna(
        "Synthetic <reference> & control",
        "ACGTACGTACGT",
        Topology::Circular,
    )
    .unwrap();
    for (name, location) in [
        ("forward CDS", Location::range(0, 4)),
        (
            "reverse CDS",
            Location::Complement(Box::new(Location::range(8, 12))),
        ),
    ] {
        source.features.push(Feature {
            kind: "CDS".into(),
            location,
            qualifiers: vec![("label".into(), Some(name.into()))],
            needs_review: false,
        });
    }
    map_export::Request {
        source_revision: source.revision(),
        source,
        restriction_sites: None,
        format,
    }
}
#[test]
fn png_map_contains_canonical_backbone_and_bundled_font_text_with_source_identity() {
    let request = circular_map(Format::Png(PngOptions { width: 1000 }));
    let source = request.source.clone();
    let rendered = request.render().unwrap();
    rendered.validate_binding(&source.revision()).unwrap();
    let Image::Png(image) = &rendered.image else {
        panic!("PNG requested")
    };
    let pixels = decode(image);
    assert_eq!(image.width, 1000);
    assert_eq!(pixel(&pixels, image.width, 10, 10), &[255, 255, 255, 255]);
    // Canonical circular backbone is centered at (500,350), radius190.
    assert_ne!(pixel(&pixels, image.width, 690, 350), &[255, 255, 255, 255]);
    assert_eq!(pixel(&pixels, image.width, 500, 300), &[255, 255, 255, 255]);
    assert_ne!(
        pixel(&pixels, image.width, 645, 205),
        &[255, 255, 255, 255],
        "the forward feature arc must survive rasterization"
    );
    let dark_title = (20..44)
        .flat_map(|y| (250..750).map(move |x| (x, y)))
        .filter(|&(x, y)| pixel(&pixels, image.width, x, y)[0] < 180)
        .count();
    assert!(dark_title > 50, "font text must not silently disappear");
    assert!(
        image
            .bytes
            .windows(source.revision().len())
            .any(|s| s == source.revision().as_bytes())
    );
    assert_eq!(request.source, source);
    assert!(rendered.warnings.iter().any(|w| w.contains("static pixel")));
    let mut stale = request;
    stale.source.name.push_str(" changed");
    assert!(stale.render().is_err());
}
#[test]
fn predicted_gel_png_retains_linear_ladder_band_positions_and_unplaced_circle_evidence() {
    let linear = Molecule::dna("1000 bp", &"ACGT".repeat(250), Topology::Linear).unwrap();
    let circular = Molecule::dna("circle", &"ACGT".repeat(250), Topology::Circular).unwrap();
    let report = gel::predict_gel(
        &[
            LaneInput {
                label: "linear".into(),
                molecule: linear,
                preparation: Preparation::Uncut,
                mass_ng: Some(20.),
            },
            LaneInput {
                label: "unknown circle".into(),
                molecule: circular,
                preparation: Preparation::Uncut,
                mass_ng: None,
            },
        ],
        &Calibration::default(),
    )
    .unwrap();
    assert_eq!(report.lanes[0].bands[0].migration, 0.5);
    assert_eq!(report.lanes[1].unpositioned_fragments, vec![0]);
    let request =
        gel_export::Request::new(report.clone(), Format::Png(PngOptions { width: 920 })).unwrap();
    let revision = request.revision().unwrap();
    let image = request.render().unwrap();
    image.validate_binding(&revision).unwrap();
    let Image::Png(png) = image.image else {
        panic!("PNG requested")
    };
    assert_eq!((png.width, png.height), (920, 680));
    let pixels = decode(&png);
    // The illustrative ladder's1000bp point is0.5 of the110..510 migration axis.
    assert_eq!(pixel(&pixels, 920, 170, 310), &[181, 238, 214, 255]);
    assert_eq!(pixel(&pixels, 920, 170, 290), &[19, 42, 39, 255]);
    assert_eq!(
        pixel(&pixels, 920, 260, 310),
        &[19, 42, 39, 255],
        "unknown circular mobility is not a fake linear band"
    );
    assert_eq!(request.report.lanes[1].unpositioned_fragments, vec![0]);
    let mut stale = request;
    stale.report.assumptions.push("changed evidence".into());
    assert!(stale.render().is_err());
}
#[test]
fn png_bounds_crc_compressed_stream_and_context_are_validated_without_unbounded_allocation() {
    let mut request = circular_map(Format::Png(PngOptions { width: 400 }));
    let Image::Png(valid) = request.render().unwrap().image else {
        unreachable!()
    };
    let mut forged = valid.clone();
    forged.width = u32::MAX;
    assert!(forged.validate_binding(&request.source_revision).is_err());
    assert!(valid.validate_binding(&"f".repeat(64)).is_err());
    let mut oversized = valid.clone();
    oversized
        .bytes
        .resize(collagene_core::output::raster::MAX_PNG_BYTES + 1, 0);
    assert!(matches!(
        oversized.validate_binding(&request.source_revision),
        Err(collagene_model::Error::Format(_))
    ));
    let mut damaged = valid.clone();
    damaged.bytes[35] ^= 1;
    damaged.content_hash = blake3::hash(&damaged.bytes).to_hex().to_string();
    assert!(damaged.validate_binding(&request.source_revision).is_err());
    let mut damaged = valid;
    let mut offset = 8;
    loop {
        let len =
            u32::from_be_bytes(damaged.bytes[offset..offset + 4].try_into().unwrap()) as usize;
        if &damaged.bytes[offset + 4..offset + 8] == b"IDAT" {
            damaged.bytes[offset + 8] = 0; // Invalid zlib header, with all outer checksums recomputed.
            let crc = crc32fast::hash(&damaged.bytes[offset + 4..offset + 8 + len]);
            damaged.bytes[offset + 8 + len..offset + 12 + len].copy_from_slice(&crc.to_be_bytes());
            break;
        }
        offset += len + 12;
    }
    damaged.content_hash = blake3::hash(&damaged.bytes).to_hex().to_string();
    // UI admission intentionally checks bounded container/binding structure, not pixel decoding.
    damaged.validate_binding(&request.source_revision).unwrap();
    let decoder = png::Decoder::new_with_limits(
        Cursor::new(&damaged.bytes),
        png::Limits {
            bytes: 32 * 1024 * 1024,
        },
    );
    let decode_failed = match decoder.read_info() {
        Err(_) => true,
        Ok(mut reader) => {
            let mut pixels = vec![0; reader.output_buffer_size().unwrap()];
            reader.next_frame(&mut pixels).is_err()
        }
    };
    assert!(decode_failed);
    request.format = Format::Png(PngOptions { width: 4097 });
    assert!(request.render().is_err());
    request.format = Format::Png(PngOptions { width: 4096 });
    assert!(
        request.render().is_err(),
        "square-like map exceeds4million pixels at4096px width"
    );
}
#[test]
fn map_and_gel_image_jobs_share_current_stateless_session_and_reply_bounds() {
    let source = circular_map(Format::Png(PngOptions { width: 400 }));
    let map = Request {
        id: 11,
        revision: "a".repeat(64),
        job: AnalysisJob::RenderMap { request: source },
    };
    let molecule = Molecule::dna("linear", &"ACGT".repeat(250), Topology::Linear).unwrap();
    let report = gel::predict_gel(
        &[LaneInput {
            label: "linear".into(),
            molecule,
            preparation: Preparation::Uncut,
            mass_ng: None,
        }],
        &Calibration::default(),
    )
    .unwrap();
    let export = gel_export::Request::new(report, Format::Png(PngOptions { width: 460 })).unwrap();
    let gel = Request {
        id: 12,
        revision: export.revision().unwrap(),
        job: AnalysisJob::RenderGelImage { request: export },
    };
    let mut session = protocol::session::AnalysisSession::default();
    for request in [map, gel] {
        let bytes = request.encode().unwrap();
        assert!(bytes.len() < 32 * 1024 * 1024);
        let stateless = protocol::execution::process_request(&bytes).unwrap();
        let session = session.process_request(&bytes).unwrap();
        assert_eq!(stateless, session);
        // Neither new image reply may acquire an archive-sized budget. Prefix admission
        // refuses a declared oversized body before copying/decoding that body.
        let mut prefix = stateless[..stateless.len().min(256)].to_vec();
        let oversized = (32 * 1024 * 1024 + 1) as u32;
        prefix[6..10].copy_from_slice(&oversized.to_le_bytes());
        assert!(matches!(
            protocol::response::Response::admitted_frame_limit(&prefix, 42 + oversized as usize),
            Err(collagene_model::Error::TooLarge)
        ));

        let reply = protocol::response::Response::decode(&stateless).unwrap();
        assert_eq!(reply.id, request.id);
        assert_eq!(reply.revision, request.revision);
        match reply.result {
            AnalysisResult::MapExported(r) => assert!(matches!(r.image, Image::Png(_))),
            AnalysisResult::GelImageExported(r) => r.validate_binding(&request.revision).unwrap(),
            other => panic!("unexpected {other:?}"),
        }
        let mut old = bytes;
        old[4..6].copy_from_slice(&(protocol::schema::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
        assert!(protocol::execution::process_request(&old).is_err());
    }
}
