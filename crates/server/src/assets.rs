use crate::{Error, Result};
use axum::body::Bytes;
use cap_fs_ext::{DirExt, FollowSymlinks, OpenOptionsFollowExt};
use cap_std::{
    ambient_authority,
    fs::{Dir, OpenOptions},
};
use std::{collections::BTreeMap, io::Read, path::Path};

pub const MAX_ASSET_FILES: usize = 4096;
pub const MAX_ASSET_FILE_BYTES: usize = 128 * 1024 * 1024;
pub const MAX_ASSET_TOTAL_BYTES: usize = 512 * 1024 * 1024;
const MAX_ASSET_DEPTH: usize = 32;
#[derive(Clone)]
pub(crate) struct Asset {
    pub bytes: Bytes,
    pub mime: &'static str,
    pub etag: String,
}
/// All HTTP reads are from this bounded immutable memory snapshot. No request
/// path is ever used for a filesystem operation.
pub struct AssetSnapshot {
    pub(crate) files: BTreeMap<String, Asset>,
}
impl AssetSnapshot {
    pub fn load(root: &Path) -> Result<Self> {
        if std::fs::symlink_metadata(root)?.file_type().is_symlink() {
            return Err(Error::Configuration("asset root symlink"));
        }
        let root = Dir::open_ambient_dir(root, ambient_authority())?;
        let mut files = BTreeMap::new();
        let mut total = 0;
        snapshot(&root, "", 0, &mut files, &mut total)?;
        let index = files
            .get_mut("/index.html")
            .ok_or(Error::Configuration("asset index.html missing"))?;
        let enabled = enable_index_marker(&index.bytes, "collagene-library")?;
        let enabled = enable_index_marker(&enabled, "collagene-source-gateway")?;
        let adjusted = total - index.bytes.len() + enabled.len();
        if enabled.len() > MAX_ASSET_FILE_BYTES || adjusted > MAX_ASSET_TOTAL_BYTES {
            return Err(Error::Configuration("asset byte limit"));
        }
        index.etag = format!("\"{}\"", blake3::hash(&enabled).to_hex());
        index.bytes = Bytes::from(enabled);
        Ok(Self { files })
    }
    pub fn file_count(&self) -> usize {
        self.files.len()
    }
}
fn snapshot(
    directory: &Dir,
    prefix: &str,
    depth: usize,
    files: &mut BTreeMap<String, Asset>,
    total: &mut usize,
) -> Result<()> {
    if depth > MAX_ASSET_DEPTH {
        return Err(Error::Configuration("asset directory depth"));
    }
    for entry in directory.entries()? {
        let entry = entry?;
        let name = entry
            .file_name()
            .into_string()
            .map_err(|_| Error::Configuration("asset filename encoding"))?;
        let kind = entry.file_type()?;
        if kind.is_symlink() {
            return Err(Error::Configuration("asset symlinks are forbidden"));
        }
        if name.starts_with('.') {
            continue;
        }
        if name.contains(['\\', '%', '?', '#']) || name.chars().any(char::is_control) {
            return Err(Error::Configuration("asset filename"));
        }
        let path = format!("{prefix}/{name}");
        if kind.is_dir() {
            let child = directory.open_dir_nofollow(&name)?;
            snapshot(&child, &path, depth + 1, files, total)?;
            continue;
        }
        if !kind.is_file() {
            return Err(Error::Configuration("asset must be regular file"));
        }
        if files.len() >= MAX_ASSET_FILES {
            return Err(Error::Configuration("asset file count"));
        }
        let budget = MAX_ASSET_FILE_BYTES.min(MAX_ASSET_TOTAL_BYTES.saturating_sub(*total));
        let mut options = OpenOptions::new();
        options.read(true).follow(FollowSymlinks::No);
        let mut file = directory.open_with(&name, &options)?;
        let metadata = file.metadata()?;
        if !metadata.is_file() || metadata.len() > budget as u64 {
            return Err(Error::Configuration("asset byte limit"));
        }
        let mut bytes = Vec::new();
        file.by_ref()
            .take(budget as u64 + 1)
            .read_to_end(&mut bytes)?;
        if bytes.len() > budget {
            return Err(Error::Configuration("asset byte limit"));
        }
        *total += bytes.len();
        let etag = format!("\"{}\"", blake3::hash(&bytes).to_hex());
        files.insert(
            path,
            Asset {
                bytes: Bytes::from(bytes),
                mime: mime(&name),
                etag,
            },
        );
    }
    Ok(())
}
fn mime(name: &str) -> &'static str {
    match name.rsplit('.').next().unwrap_or("") {
        "html" => "text/html; charset=utf-8",
        "js" | "mjs" => "text/javascript; charset=utf-8",
        "css" => "text/css; charset=utf-8",
        "wasm" => "application/wasm",
        "json" | "map" => "application/json",
        "svg" => "image/svg+xml",
        "png" => "image/png",
        "jpg" | "jpeg" => "image/jpeg",
        "webp" => "image/webp",
        "gif" => "image/gif",
        "ico" => "image/x-icon",
        "woff" => "font/woff",
        "woff2" => "font/woff2",
        "txt" => "text/plain; charset=utf-8",
        _ => "application/octet-stream",
    }
}
pub(crate) fn asset_path(raw: &str) -> Option<String> {
    let mut decoded = Vec::with_capacity(raw.len());
    let bytes = raw.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' {
            if i + 2 >= bytes.len() {
                return None;
            }
            let pair = std::str::from_utf8(&bytes[i + 1..i + 3]).ok()?;
            decoded.push(u8::from_str_radix(pair, 16).ok()?);
            i += 3;
        } else {
            decoded.push(bytes[i]);
            i += 1;
        }
    }
    let path = String::from_utf8(decoded).ok()?;
    if !path.starts_with('/')
        || path.contains('\\')
        || path.chars().any(char::is_control)
        || path.split('/').any(|part| part.starts_with('.'))
    {
        return None;
    }
    Some(if path == "/" {
        "/index.html".into()
    } else {
        path
    })
}

const MAX_INDEX_BYTES: usize = 2 * 1024 * 1024;
/// Transform only the trusted index snapshot. This deliberately accepts a
/// restricted HTML document structure and fails closed on ambiguous markers;
/// it does not rewrite script/link text or other asset bytes.
fn enable_index_marker(bytes: &[u8], marker_name: &str) -> Result<Vec<u8>> {
    let canonical_meta = format!("<meta name=\"{marker_name}\" content=\"enabled\">");
    let bad = || Error::Configuration("ambiguous capability marker/index HTML");
    if bytes.len() > MAX_INDEX_BYTES {
        return Err(bad());
    }
    let html = std::str::from_utf8(bytes).map_err(|_| bad())?;
    let lower = html.to_ascii_lowercase();
    let mut cursor = 0;
    let mut head = None;
    let mut in_head = false;
    let mut closed_head = false;
    let mut marker = None;
    while let Some(relative) = html[cursor..].find('<') {
        let start = cursor + relative;
        if html[start..].starts_with("<!--") {
            cursor = start + 4 + html[start + 4..].find("-->").ok_or_else(bad)? + 3;
            continue;
        }
        let end = tag_end(html, start).ok_or_else(bad)?;
        let mut inside = html[start + 1..end - 1].trim();
        let closing = inside.starts_with('/');
        if closing {
            inside = inside[1..].trim_start();
        }
        let name_end = inside
            .find(|c: char| !c.is_ascii_alphanumeric() && !matches!(c, '-' | ':'))
            .unwrap_or(inside.len());
        let name = &inside[..name_end];
        cursor = end;
        if name.is_empty() {
            continue;
        } // doctype/declaration, never an element marker
        if name.eq_ignore_ascii_case("template") {
            return Err(bad());
        }
        if name.eq_ignore_ascii_case("head") {
            if closing {
                if !in_head || closed_head {
                    return Err(bad());
                }
                in_head = false;
                closed_head = true;
            } else {
                if head.is_some() || closed_head || inside.trim_end().ends_with('/') {
                    return Err(bad());
                }
                head = Some(end);
                in_head = true;
            }
        }
        if !closing && name.eq_ignore_ascii_case("meta") {
            let attributes = attributes(&inside[name_end..]).ok_or_else(bad)?;
            let names = attributes
                .iter()
                .filter(|(key, _)| key.eq_ignore_ascii_case("name"))
                .collect::<Vec<_>>();
            let is_marker = names
                .iter()
                .any(|(_, value)| value.trim().eq_ignore_ascii_case(marker_name));
            if is_marker {
                let contents = attributes
                    .iter()
                    .filter(|(key, _)| key.eq_ignore_ascii_case("content"))
                    .collect::<Vec<_>>();
                let mut keys = std::collections::BTreeSet::new();
                if !in_head
                    || marker.is_some()
                    || names.len() != 1
                    || contents.len() != 1
                    || names[0].1 != marker_name
                    || contents[0].1 != "enabled"
                    || attributes
                        .iter()
                        .any(|(key, _)| !keys.insert(key.to_ascii_lowercase()))
                {
                    return Err(bad());
                }
                marker = Some((start, end));
            }
        }
        if !closing
            && [
                "script", "style", "title", "textarea", "xmp", "iframe", "noembed", "noframes",
            ]
            .iter()
            .any(|raw| name.eq_ignore_ascii_case(raw))
        {
            let needle = format!("</{}", name.to_ascii_lowercase());
            let mut search = cursor;
            loop {
                let close = search + lower[search..].find(&needle).ok_or_else(bad)?;
                let after = close + needle.len();
                if lower
                    .as_bytes()
                    .get(after)
                    .is_some_and(|c| c.is_ascii_whitespace() || matches!(c, b'>' | b'/'))
                {
                    cursor = tag_end(html, close).ok_or_else(bad)?;
                    break;
                }
                search = after;
            }
        }
    }
    let insertion = head.filter(|_| closed_head && !in_head).ok_or_else(bad)?;
    let mut result = Vec::with_capacity(bytes.len() + canonical_meta.len() + 2);
    if let Some((start, end)) = marker {
        result.extend_from_slice(&bytes[..start]);
        result.extend_from_slice(canonical_meta.as_bytes());
        result.extend_from_slice(&bytes[end..]);
    } else {
        result.extend_from_slice(&bytes[..insertion]);
        result.push(b'\n');
        result.extend_from_slice(canonical_meta.as_bytes());
        result.push(b'\n');
        result.extend_from_slice(&bytes[insertion..]);
    }
    Ok(result)
}
fn tag_end(html: &str, start: usize) -> Option<usize> {
    let mut quote = None;
    for (offset, byte) in html.as_bytes()[start + 1..].iter().copied().enumerate() {
        match (quote, byte) {
            (Some(active), value) if active == value => quote = None,
            (None, b'\'' | b'"') => quote = Some(byte),
            (None, b'>') => return Some(start + offset + 2),
            _ => (),
        }
    }
    None
}
fn attributes(mut text: &str) -> Option<Vec<(&str, String)>> {
    let mut values = Vec::new();
    loop {
        text = text.trim_start_matches(|c: char| c.is_ascii_whitespace());
        if text.is_empty() || text == "/" {
            return Some(values);
        }
        let end = text
            .find(|c: char| c.is_ascii_whitespace() || matches!(c, '=' | '/'))
            .unwrap_or(text.len());
        if end == 0
            || !text[..end]
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || matches!(c, b'-' | b'_' | b':'))
        {
            return None;
        }
        let key = &text[..end];
        text = text[end..].trim_start_matches(|c: char| c.is_ascii_whitespace());
        let mut value = String::new();
        if let Some(rest) = text.strip_prefix('=') {
            text = rest.trim_start_matches(|c: char| c.is_ascii_whitespace());
            if text.starts_with(['\'', '"']) {
                let quote = text.as_bytes()[0] as char;
                let end = text[1..].find(quote)? + 1;
                value = decode_attribute(&text[1..end])?;
                text = &text[end + 1..];
            } else {
                let end = text
                    .find(|c: char| c.is_ascii_whitespace())
                    .unwrap_or(text.len());
                if end == 0 {
                    return None;
                }
                value = decode_attribute(&text[..end])?;
                text = &text[end..];
            }
        }
        values.push((key, value));
    }
}
fn decode_attribute(value: &str) -> Option<String> {
    let mut result = String::new();
    let mut rest = value;
    while let Some(amp) = rest.find('&') {
        result.push_str(&rest[..amp]);
        rest = &rest[amp + 1..];
        if let Some(number) = rest.strip_prefix('#') {
            let (number, radix) = if number.starts_with(['x', 'X']) {
                (&number[1..], 16)
            } else {
                (number, 10)
            };
            let digits = number
                .find(|c: char| !c.is_digit(radix))
                .unwrap_or(number.len());
            if digits == 0 {
                return None;
            }
            result.push(char::from_u32(
                u32::from_str_radix(&number[..digits], radix).ok()?,
            )?);
            rest = &number[digits..];
            if let Some(after) = rest.strip_prefix(';') {
                rest = after;
            }
        } else {
            let entity = [
                ("amp;", '&'),
                ("quot;", '"'),
                ("apos;", '\''),
                ("lt;", '<'),
                ("gt;", '>'),
            ]
            .into_iter()
            .find(|(name, _)| rest.starts_with(name));
            if let Some((name, character)) = entity {
                result.push(character);
                rest = &rest[name.len()..];
            } else {
                result.push('&');
            }
        }
    }
    result.push_str(rest);
    Some(result)
}
