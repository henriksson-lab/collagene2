//! Portable static image encoding. Only our bounded map/gel SVG generators call the renderer.
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::io::{Cursor, Write};

pub const MAX_DIMENSION: u32 = 4096;
pub const MAX_PIXELS: u64 = 4_000_000;
pub const MAX_PNG_BYTES: usize = 24 * 1024 * 1024;
/// Gels can contain more vector primitives than maps; map owners retain their 2MiB cap.
pub const MAX_VECTOR_BYTES: usize = 16 * 1024 * 1024;
const FONT: &[u8] = include_bytes!("../../assets/fonts/LiberationSans-Regular.ttf");
const CONTEXT_KEY: &[u8] = b"CollageneContext\0";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PngOptions {
    /// Output width in pixels. Height follows the canonical diagram's aspect ratio.
    pub width: u32,
}
impl Default for PngOptions {
    fn default() -> Self {
        Self { width: 1600 }
    }
}
impl PngOptions {
    pub fn validate(&self) -> Result<(), Error> {
        if !(16..=MAX_DIMENSION).contains(&self.width) {
            return Err(bad("PNG width must be 16–4096 pixels"));
        }
        Ok(())
    }
    fn dimensions(&self, width: f32, height: f32) -> Result<(u32, u32), Error> {
        self.validate()?;
        if !width.is_finite() || !height.is_finite() || width <= 0. || height <= 0. {
            return Err(bad("invalid diagram dimensions"));
        }
        let height = (f64::from(height) / f64::from(width) * f64::from(self.width)).ceil();
        if !height.is_finite() || !(1.0..=f64::from(MAX_DIMENSION)).contains(&height) {
            return Err(Error::TooLarge);
        }
        let height = height as u32;
        check_dimensions(self.width, height)?;
        Ok((self.width, height))
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Format {
    Svg,
    Png(PngOptions),
}
impl Format {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Svg => Ok(()),
            Self::Png(options) => options.validate(),
        }
    }
    pub fn extension(self) -> &'static str {
        match self {
            Self::Svg => "svg",
            Self::Png(_) => "png",
        }
    }
    pub fn mime_type(self) -> &'static str {
        match self {
            Self::Svg => "image/svg+xml",
            Self::Png(_) => "image/png",
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PngImage {
    pub width: u32,
    pub height: u32,
    pub context_revision: String,
    pub content_hash: String,
    pub bytes: Vec<u8>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Image {
    Svg(String),
    Png(PngImage),
}
impl Image {
    pub fn format(&self) -> Format {
        match self {
            Self::Svg(_) => Format::Svg,
            Self::Png(png) => Format::Png(PngOptions { width: png.width }),
        }
    }
    pub fn bytes(&self) -> &[u8] {
        match self {
            Self::Svg(svg) => svg.as_bytes(),
            Self::Png(png) => &png.bytes,
        }
    }
    pub fn validate_binding(&self, revision: &str) -> Result<(), Error> {
        match self {
            Self::Svg(svg) => {
                if svg.is_empty()
                    || svg.len() > MAX_VECTOR_BYTES
                    || !svg.starts_with("<svg ")
                    || !svg.ends_with("</svg>")
                {
                    return Err(bad("SVG output bounds or framing differ"));
                }
                Ok(())
            }
            Self::Png(png) => png.validate_binding(revision),
        }
    }
}
fn bad(message: &str) -> Error {
    Error::Format(format!("image export: {message}"))
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn check_dimensions(width: u32, height: u32) -> Result<(), Error> {
    if width == 0
        || height == 0
        || width > MAX_DIMENSION
        || height > MAX_DIMENSION
        || u64::from(width) * u64::from(height) > MAX_PIXELS
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
impl PngImage {
    /// Bounded transport admission: framing, dimensions, CRCs, identity and content hash.
    /// This deliberately does not decompress pixels on a UI response thread. Pixels are
    /// documentary presentation; their structure is not physical/scientific authority.
    pub fn validate_binding(&self, revision: &str) -> Result<(), Error> {
        PngOptions { width: self.width }.validate()?;
        check_dimensions(self.width, self.height)?;
        if !hex(revision)
            || self.context_revision != revision
            || !hex(&self.content_hash)
            || self.bytes.len() > MAX_PNG_BYTES
            || self.bytes.len() < 57
            || &self.bytes[..8] != b"\x89PNG\r\n\x1a\n"
            || blake3::hash(&self.bytes).to_hex().as_str() != self.content_hash
        {
            return Err(bad("PNG content or captured context differs"));
        }
        let mut offset = 8usize;
        let mut header = false;
        let mut context = false;
        let mut data = false;
        let mut ended = false;
        let mut chunks = 0usize;
        while offset < self.bytes.len() {
            chunks += 1;
            if chunks > 4096 || self.bytes.len() - offset < 12 {
                return Err(bad("PNG chunk framing exceeds bounds"));
            }
            let len =
                u32::from_be_bytes(self.bytes[offset..offset + 4].try_into().unwrap()) as usize;
            let end = offset
                .checked_add(12)
                .and_then(|n| n.checked_add(len))
                .filter(|n| *n <= self.bytes.len())
                .ok_or_else(|| bad("truncated PNG chunk"))?;
            let kind = &self.bytes[offset + 4..offset + 8];
            let body = &self.bytes[offset + 8..end - 4];
            let crc = u32::from_be_bytes(self.bytes[end - 4..end].try_into().unwrap());
            if crc32fast::hash(&self.bytes[offset + 4..end - 4]) != crc {
                return Err(bad("PNG chunk checksum differs"));
            }
            match kind {
                b"IHDR" if !header && offset == 8 && len == 13 => {
                    if body[..4] != self.width.to_be_bytes()
                        || body[4..8] != self.height.to_be_bytes()
                        || body[8..] != [8, 6, 0, 0, 0]
                    {
                        return Err(bad("PNG dimensions or RGBA encoding differs"));
                    }
                    header = true;
                }
                b"tEXt" if header && !data && !context => {
                    if body.strip_prefix(CONTEXT_KEY) != Some(revision.as_bytes()) {
                        return Err(bad("PNG embedded source context differs"));
                    }
                    context = true;
                }
                b"IDAT" if header && context && !ended => data = true,
                b"IEND" if data && len == 0 && end == self.bytes.len() => ended = true,
                _ => return Err(bad("unexpected PNG export chunk or order")),
            }
            offset = end;
        }
        if !header || !context || !data || !ended {
            return Err(bad("incomplete PNG export"));
        }
        Ok(())
    }
    /// Worker/native renderer check of the encoded pixels. Decoding is row-bounded and
    /// never invoked by Response::decode or browser capture admission.
    fn validate_pixels(&self) -> Result<(), Error> {
        let decoder = png::Decoder::new_with_limits(
            Cursor::new(&self.bytes),
            png::Limits {
                bytes: 32 * 1024 * 1024,
            },
        );
        let mut reader = decoder.read_info().map_err(|e| bad(&e.to_string()))?;
        let mut rows = 0u32;
        while let Some(row) = reader.next_row().map_err(|e| bad(&e.to_string()))? {
            if row.data().len() != self.width as usize * 4 {
                return Err(bad("PNG row size differs"));
            }
            rows += 1;
            if rows > self.height {
                return Err(Error::TooLarge);
            }
        }
        reader.finish().map_err(|e| bad(&e.to_string()))?;
        if rows != self.height {
            return Err(bad("PNG row count differs"));
        }
        Ok(())
    }
}
struct BoundedBytes(Vec<u8>);
impl Write for BoundedBytes {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() > MAX_PNG_BYTES.saturating_sub(self.0.len()) {
            return Err(std::io::Error::other("PNG byte budget exceeded"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
/// This is deliberately crate-private: SVG is generated by our typed map or gel owners,
/// never uploaded SVG markup, user paths, embedded images or remotely loaded fonts.
pub(crate) fn render_generated_svg(
    svg: &str,
    options: PngOptions,
    revision: &str,
) -> Result<(PngImage, Vec<String>), Error> {
    if svg.len() > MAX_VECTOR_BYTES || !hex(revision) {
        return Err(Error::TooLarge);
    }
    options.validate()?;
    let face = ttf_parser::Face::parse(FONT, 0).map_err(|_| bad("bundled font is invalid"))?;
    let document = roxmltree::Document::parse_with_options(
        svg,
        roxmltree::ParsingOptions {
            allow_dtd: false,
            nodes_limit: 100_000,
            entity_resolver: None,
        },
    )
    .map_err(|e| bad(&e.to_string()))?;
    let mut replacements = Vec::new();
    let mut missing = 0usize;
    for node in document.descendants() {
        if node.is_element()
            && matches!(
                node.tag_name().name(),
                "image" | "script" | "foreignObject" | "style"
            )
        {
            return Err(bad(
                "non-generated or external SVG resources are unsupported",
            ));
        }
        for attribute in node.attributes() {
            if matches!(attribute.name(), "href" | "src") && !attribute.value().starts_with('#') {
                return Err(bad("external SVG references are unsupported"));
            }
        }
        if node.is_text()
            && node.parent_element().is_some_and(|parent| {
                matches!(parent.tag_name().name(), "text" | "tspan" | "textPath")
            })
        {
            let mut text = String::new();
            let mut changed = false;
            for c in node.text().unwrap_or_default().chars() {
                if !c.is_whitespace() && face.glyph_index(c).is_none() {
                    use std::fmt::Write;
                    let _ = write!(text, "[U+{:04X}]", c as u32);
                    missing += 1;
                    changed = true;
                } else {
                    match c {
                        '&' => text.push_str("&amp;"),
                        '<' => text.push_str("&lt;"),
                        '>' => text.push_str("&gt;"),
                        _ => text.push(c),
                    }
                }
            }
            if changed {
                replacements.push((node.range(), text));
            }
        }
    }
    let mut source = svg.to_owned();
    for (range, text) in replacements.into_iter().rev() {
        source.replace_range(range, &text);
    }
    if source.len() > MAX_VECTOR_BYTES {
        return Err(Error::TooLarge);
    }
    let mut config = resvg::usvg::Options {
        font_family: "Liberation Sans".into(),
        resources_dir: None,
        image_href_resolver: resvg::usvg::ImageHrefResolver {
            resolve_data: Box::new(|_, _, _| None),
            resolve_string: Box::new(|_, _| None),
        },
        ..Default::default()
    };
    config.fontdb_mut().load_font_data(FONT.to_vec());
    config.fontdb_mut().set_sans_serif_family("Liberation Sans");
    config.fontdb_mut().set_serif_family("Liberation Sans");
    config.fontdb_mut().set_monospace_family("Liberation Sans");
    let tree = resvg::usvg::Tree::from_str(&source, &config).map_err(|e| bad(&e.to_string()))?;
    let (width, height) = options.dimensions(tree.size().width(), tree.size().height())?;
    let mut pixmap = resvg::tiny_skia::Pixmap::new(width, height).ok_or(Error::TooLarge)?;
    let scale = width as f32 / tree.size().width();
    resvg::render(
        &tree,
        resvg::tiny_skia::Transform::from_scale(scale, scale),
        &mut pixmap.as_mut(),
    );
    let rgba = pixmap.take_demultiplied();
    let mut bytes = BoundedBytes(Vec::new());
    {
        let mut encoder = png::Encoder::new(&mut bytes, width, height);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        encoder
            .add_text_chunk("CollageneContext".into(), revision.into())
            .map_err(|e| bad(&e.to_string()))?;
        let mut writer = encoder.write_header().map_err(|e| bad(&e.to_string()))?;
        writer
            .write_image_data(&rgba)
            .map_err(|e| bad(&e.to_string()))?;
        writer.finish().map_err(|e| bad(&e.to_string()))?;
    }
    let image = PngImage {
        width,
        height,
        context_revision: revision.into(),
        content_hash: blake3::hash(&bytes.0).to_hex().to_string(),
        bytes: bytes.0,
    };
    image.validate_binding(revision)?;
    image.validate_pixels()?;
    let mut warnings = vec!["PNG is a static pixel image; selectable vector text, coordinate attributes and SVG evidence metadata are not retained. Keep the source document or gel report for complete evidence. Rendering uses the bundled Liberation Sans font, independently of installed fonts.".into()];
    if missing != 0 {
        warnings.push(format!("The bundled font cannot represent {missing} label scalars. They are visibly escaped as [U+XXXX]; original names remain in the source document or report."));
    }
    Ok((image, warnings))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn generated_images_cannot_resolve_files_urls_embedded_images_or_stylesheets() {
        for body in [
            "<image href=\"file:///unavailable.png\" width=\"16\" height=\"16\"/>",
            "<use href=\"https://example.invalid/source.svg#shape\"/>",
            "<image href=\"data:image/svg+xml;base64,AAAA\"/>",
            "<style>@import url(https://example.invalid/fonts.css);</style>",
            "<script>alert(1)</script>",
            "<foreignObject width=\"16\" height=\"16\"/>",
        ] {
            let svg = format!(
                "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"16\" height=\"16\">{body}</svg>"
            );
            assert!(render_generated_svg(&svg, PngOptions { width: 16 }, &"a".repeat(64)).is_err());
        }
    }
    #[test]
    fn missing_font_scalars_are_visible_escapes_with_explicit_warnings() {
        let svg = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"400\" height=\"40\"><rect width=\"400\" height=\"40\" fill=\"white\"/><text x=\"4\" y=\"28\" font-size=\"20\">A漢🧬&amp;B</text></svg>";
        let (image, warnings) =
            render_generated_svg(svg, PngOptions { width: 400 }, &"b".repeat(64)).unwrap();
        assert!(
            warnings
                .iter()
                .any(|w| w.contains("2 label scalars") && w.contains("[U+XXXX]"))
        );
        image.validate_binding(&"b".repeat(64)).unwrap();
    }
}
