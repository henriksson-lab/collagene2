//! Source-preserving documentary images. A preview is never measurement or verification.
//!
//! Import in a worker: codecs are synchronous and bounded, but still computationally expensive.
//! This crate performs no filesystem, network, sequence, or physical-material operations.

mod header;

/// Bump when preview admission or rendering semantics change.
pub const PREVIEW_ENGINE_VERSION: &str = "documentary-image-1";

use image::{GenericImageView, ImageDecoder, ImageEncoder, ImageReader};
use std::io::{Cursor, Write};

/// Hard ceilings also apply to caller-supplied admission settings.
#[derive(Clone, Copy, Debug)]
pub struct Limits {
    pub input_bytes: usize,
    pub dimension: u32,
    pub source_pixels: u64,
    pub decoded_bytes: u64,
    /// Conservative admission estimate, plus the codec's best-effort allocation limit.
    /// This is not a process-wide allocator or a promise about codec allocator overhead.
    pub working_bytes: u64,
    pub preview_pixels: u64,
    pub preview_dimension: u32,
    pub preview_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            input_bytes: 16 * 1024 * 1024,
            dimension: 16_384,
            source_pixels: 4_194_304,
            decoded_bytes: 64 * 1024 * 1024,
            working_bytes: 256 * 1024 * 1024,
            preview_pixels: 1_000_000,
            preview_dimension: 1_600,
            preview_bytes: 5 * 1024 * 1024,
        }
    }
}
impl Limits {
    fn validate(self) -> Result<Self, Error> {
        let hard = Self::default();
        if self.input_bytes == 0
            || self.input_bytes > hard.input_bytes
            || self.dimension == 0
            || self.dimension > hard.dimension
            || self.source_pixels == 0
            || self.source_pixels > hard.source_pixels
            || self.decoded_bytes == 0
            || self.decoded_bytes > hard.decoded_bytes
            || self.working_bytes == 0
            || self.working_bytes > hard.working_bytes
            || self.preview_pixels == 0
            || self.preview_pixels > hard.preview_pixels
            || self.preview_dimension == 0
            || self.preview_dimension > hard.preview_dimension
            || self.preview_bytes == 0
            || self.preview_bytes > hard.preview_bytes
        {
            return Err(Error::Limit("invalid admission limits"));
        }
        Ok(self)
    }

    pub(crate) fn dimensions(self, width: u32, height: u32) -> Result<(), Error> {
        let pixels = u64::from(width) * u64::from(height);
        if width == 0
            || height == 0
            || width > self.dimension
            || height > self.dimension
            || pixels > self.source_pixels
        {
            return Err(Error::Limit("source dimensions/pixel count"));
        }
        // Padding covers JPEG MCU edges. 64 bytes/pixel conservatively admits source,
        // coefficient/scanline intermediates, and source conversion. Small RGBA preview,
        // its rotated copy and PNG encoder workspace are admitted separately.
        let padded = (u64::from(width) + 32) * (u64::from(height) + 32);
        let estimate =
            padded * 64 + self.preview_pixels * 24 + self.input_bytes as u64 * 2 + 4 * 1024 * 1024;
        if pixels * 16 > self.decoded_bytes || estimate > self.working_bytes {
            return Err(Error::Limit("decoded/intermediate allocation admission"));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Jpeg,
    Tiff,
    Png,
}
impl Format {
    pub fn media_type(self) -> &'static str {
        match self {
            Self::Jpeg => "image/jpeg",
            Self::Tiff => "image/tiff",
            Self::Png => "image/png",
        }
    }
    pub fn extension(self) -> &'static str {
        match self {
            Self::Jpeg => "jpg",
            Self::Tiff => "tif",
            Self::Png => "png",
        }
    }
    fn decoder_format(self) -> image::ImageFormat {
        match self {
            Self::Jpeg => image::ImageFormat::Jpeg,
            Self::Tiff => image::ImageFormat::Tiff,
            Self::Png => image::ImageFormat::Png,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Warning {
    /// Display-only: no bands, intensities, concentrations, or verification inferred.
    DocumentaryPreview,
    /// Original metadata is retained only in original bytes, not copied into the preview.
    MetadataNotTransferred,
    /// Color profiles, acquisition calibration, and physical dimensions are not applied.
    UncalibratedDisplay,
    /// Only the first/default image is previewed; other frames/pages remain in original.
    AdditionalImagesNotPreviewed,
    /// The original has >8 bits/sample; the preview has 8 bits/channel.
    PrecisionReducedForPreview,
    DownsampledForPreview,
}

impl Warning {
    pub fn message(self) -> &'static str {
        match self {
            Self::DocumentaryPreview => {
                "Documentary preview; no measurement or biological verification inferred."
            }
            Self::MetadataNotTransferred => {
                "Original metadata is retained in the source file, not copied to the preview."
            }
            Self::UncalibratedDisplay => {
                "ICC color profiles, physical units and acquisition calibration are not applied."
            }
            Self::AdditionalImagesNotPreviewed => {
                "Only the first/default image is previewed; additional images remain in the original."
            }
            Self::PrecisionReducedForPreview => {
                "The preview uses 8 bits per channel; higher-precision samples remain in the original."
            }
            Self::DownsampledForPreview => {
                "The preview is downsampled; the original image is unchanged."
            }
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("image resource limit: {0}")]
    Limit(&'static str),
    #[error("malformed image: {0}")]
    Malformed(&'static str),
    #[error("unsupported image: {0}")]
    Unsupported(&'static str),
    #[error("image codec rejected input: {0}")]
    Codec(#[from] image::ImageError),
}

/// The only construction path validates the original and renders a bounded preview.
/// No unchecked Deserialize implementation is provided.
pub struct Evidence {
    original: Vec<u8>,
    content_hash: [u8; 32],
    format: Format,
    source_dimensions: (u32, u32),
    source_orientation: u8,
    preview_dimensions: (u32, u32),
    preview: Vec<u8>,
    warnings: Vec<Warning>,
}
impl Evidence {
    pub fn prepare(original: Vec<u8>, limits: Limits) -> Result<Self, Error> {
        let limits = limits.validate()?;
        if original.len() > limits.input_bytes {
            return Err(Error::Limit("compressed input bytes"));
        }
        let header = header::Header::read(&original, limits)?;
        limits.dimensions(header.width, header.height)?;
        let mut decoder_limits = image::Limits::default();
        decoder_limits.max_image_width = Some(limits.dimension);
        decoder_limits.max_image_height = Some(limits.dimension);
        decoder_limits.max_alloc = Some(limits.working_bytes / 2);
        let mut reader = ImageReader::with_format(
            Cursor::new(original.as_slice()),
            header.format.decoder_format(),
        );
        reader.limits(decoder_limits);
        let decoder = reader.into_decoder()?;
        if decoder.dimensions() != (header.width, header.height) {
            return Err(Error::Malformed("decoder/header dimensions disagree"));
        }
        if decoder.total_bytes() > limits.decoded_bytes {
            return Err(Error::Limit("decoder output bytes"));
        }
        let precision = header.precision > 8
            || decoder.color_type().bits_per_pixel()
                > u16::from(decoder.color_type().channel_count()) * 8;
        let decoded = image::DynamicImage::from_decoder(decoder)?;
        let (width, height) = preview_dimensions(header.width, header.height, limits);
        let downsampled = (width, height) != (header.width, header.height);
        // Sample straight into the small RGBA8 raster. Avoid image::resize's
        // full-width float intermediate, especially for extremely wide sources.
        // Nearest-neighbour is documentary and makes no quantitative inference.
        let small = image::RgbaImage::from_fn(width, height, |x, y| {
            let sx =
                ((u64::from(x) * 2 + 1) * u64::from(header.width) / (u64::from(width) * 2)) as u32;
            let sy = ((u64::from(y) * 2 + 1) * u64::from(header.height) / (u64::from(height) * 2))
                as u32;
            decoded.get_pixel(sx, sy)
        });
        drop(decoded);
        let mut preview_image = image::DynamicImage::ImageRgba8(small);
        preview_image.apply_orientation(
            image::metadata::Orientation::from_exif(header.orientation)
                .ok_or(Error::Malformed("invalid orientation"))?,
        );
        let dimensions = (preview_image.width(), preview_image.height());
        let mut output = LimitedOutput {
            bytes: Vec::new(),
            maximum: limits.preview_bytes,
        };
        image::codecs::png::PngEncoder::new(&mut output).write_image(
            preview_image.as_bytes(),
            dimensions.0,
            dimensions.1,
            image::ExtendedColorType::Rgba8,
        )?;
        let mut warnings = vec![
            Warning::DocumentaryPreview,
            Warning::MetadataNotTransferred,
            Warning::UncalibratedDisplay,
        ];
        if header.multiple {
            warnings.push(Warning::AdditionalImagesNotPreviewed);
        }
        if precision {
            warnings.push(Warning::PrecisionReducedForPreview);
        }
        if downsampled {
            warnings.push(Warning::DownsampledForPreview);
        }
        Ok(Self {
            content_hash: *blake3::hash(&original).as_bytes(),
            original,
            format: header.format,
            source_dimensions: (header.width, header.height),
            source_orientation: header.orientation,
            preview_dimensions: dimensions,
            preview: output.bytes,
            warnings,
        })
    }
    pub fn original_bytes(&self) -> &[u8] {
        &self.original
    }
    pub fn into_original(self) -> Vec<u8> {
        self.original
    }
    pub fn content_hash(&self) -> &[u8; 32] {
        &self.content_hash
    }
    pub fn format(&self) -> Format {
        self.format
    }
    /// Encoded raster dimensions, before EXIF orientation.
    pub fn source_dimensions(&self) -> (u32, u32) {
        self.source_dimensions
    }
    pub fn source_orientation(&self) -> u8 {
        self.source_orientation
    }
    pub fn preview_dimensions(&self) -> (u32, u32) {
        self.preview_dimensions
    }
    pub fn preview_png(&self) -> &[u8] {
        &self.preview
    }
    pub fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

fn preview_dimensions(width: u32, height: u32, limits: Limits) -> (u32, u32) {
    // Integer-only ratio fitting; no float rounding can cross the pixel limit.
    let mut divisor = 1u32;
    loop {
        let w = (width / divisor).max(1);
        let h = (height / divisor).max(1);
        if w <= limits.preview_dimension
            && h <= limits.preview_dimension
            && u64::from(w) * u64::from(h) <= limits.preview_pixels
        {
            return (w, h);
        }
        divisor += 1;
    }
}

struct LimitedOutput {
    bytes: Vec<u8>,
    maximum: usize,
}
impl Write for LimitedOutput {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() > self.maximum.saturating_sub(self.bytes.len()) {
            return Err(std::io::Error::other("PNG preview byte limit"));
        }
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl std::fmt::Debug for Evidence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Evidence")
            .field("format", &self.format)
            .field("source_dimensions", &self.source_dimensions)
            .field("source_orientation", &self.source_orientation)
            .field("original_bytes", &self.original.len())
            .field("preview_dimensions", &self.preview_dimensions)
            .field("preview_bytes", &self.preview.len())
            .field("warnings", &self.warnings)
            .finish_non_exhaustive()
    }
}
