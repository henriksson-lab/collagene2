//! Build-time immutable worker packaging. No shell commands or network access.
#![forbid(unsafe_code)]
use std::{
    collections::BTreeSet,
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};
pub const SENTINEL: &str = "__COLLAGENE_WORKER_URL__";
pub const FILES: [&str; 3] = [
    "analysis-worker_loader.js",
    "analysis-worker.js",
    "analysis-worker_bg.wasm",
];
pub const MAX_HTML_BYTES: u64 = 8 * 1024 * 1024;
pub const MAX_LOADER_BYTES: u64 = 1024 * 1024;
pub const MAX_JS_BYTES: u64 = 16 * 1024 * 1024;
pub const MAX_WASM_BYTES: u64 = 64 * 1024 * 1024;
pub const MAX_TOTAL_BYTES: u64 = 80 * 1024 * 1024;
/// Raw meta content URL cap shared with the browser worker factory.
pub const MAX_WORKER_URL_BYTES: usize = 4096;
const MARKER: &str = "<!-- collagene-worker-bundle:v1:";
static SERIAL: AtomicU64 = AtomicU64::new(0);
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Asset {
    pub name: &'static str,
    pub bytes: u64,
    pub blake3: String,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Manifest {
    pub schema_version: u16,
    pub bundle_hash: String,
    pub loader_url: String,
    pub relative_loader_path: String,
    pub files: [Asset; 3],
}
#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    InvalidStaging,
    UnsafePath,
    TooLarge(&'static str),
    InvalidLoader,
    InvalidGlue,
    InvalidWasm,
    InvalidHtml,
    MissingSentinel,
    AmbiguousSentinel,
    InvalidPublicUrl,
    CorruptBundle,
    ChangedInput,
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Io(_) => "worker asset filesystem operation failed",
            Self::InvalidStaging => "staging must be an existing ordinary directory",
            Self::UnsafePath => "symlinks or nonordinary managed paths are forbidden",
            Self::TooLarge(_) => "worker build input exceeds its byte limit",
            Self::InvalidLoader => {
                "unsupported worker loader: expected only the two named relative imports"
            }
            Self::InvalidGlue => "unsupported or malformed wasm-bindgen worker glue",
            Self::InvalidWasm => "missing WebAssembly version 1 header",
            Self::InvalidHtml => "expected one quoted collagene-worker meta content attribute",
            Self::MissingSentinel => {
                "staged HTML lacks the worker sentinel or valid prior packaging marker"
            }
            Self::AmbiguousSentinel => "staged HTML contains multiple worker sentinels",
            Self::InvalidPublicUrl => {
                "public URL must be an unambiguous root/relative or HTTP(S) directory prefix"
            }
            Self::CorruptBundle => {
                "an existing immutable bundle directory does not match its complete content"
            }
            Self::ChangedInput => "staged input changed during packaging",
        })
    }
}
impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        if let Self::Io(e) = self {
            Some(e)
        } else {
            None
        }
    }
}
impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}
type Result<T> = std::result::Result<T, Error>;
fn reparse(m: &fs::Metadata) -> bool {
    #[cfg(windows)]
    {
        use std::os::windows::fs::MetadataExt;
        if m.file_attributes() & 0x400 != 0 {
            return true;
        }
    }
    m.file_type().is_symlink()
}
fn directory(path: &Path) -> Result<()> {
    let m = fs::symlink_metadata(path)?;
    if reparse(&m) || !m.is_dir() {
        return Err(Error::UnsafePath);
    }
    Ok(())
}
fn read_file(path: &Path, limit: u64, label: &'static str) -> Result<Vec<u8>> {
    let before = fs::symlink_metadata(path)?;
    if reparse(&before) || !before.is_file() {
        return Err(Error::UnsafePath);
    }
    if before.len() > limit {
        return Err(Error::TooLarge(label));
    }
    let f = File::open(path)?;
    let opened = f.metadata()?;
    if !opened.is_file() || opened.len() != before.len() {
        return Err(Error::ChangedInput);
    }
    let capacity = usize::try_from(before.len()).map_err(|_| Error::TooLarge(label))?;
    let mut bytes = Vec::with_capacity(capacity.saturating_add(1));
    f.take(limit + 1).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > limit {
        return Err(Error::TooLarge(label));
    }
    let after = fs::symlink_metadata(path)?;
    if reparse(&after)
        || !after.is_file()
        || bytes.len() as u64 != before.len()
        || after.len() != before.len()
        || before.modified().ok() != after.modified().ok()
    {
        return Err(Error::ChangedInput);
    }
    Ok(bytes)
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'a'..=b'f'))
}
/// Domain-separated BLAKE3 over ordered filename length/name, byte length, and
/// exact content. Files cannot collide through concatenation or name swapping.
pub fn worker_hash(files: [&[u8]; 3]) -> String {
    let mut h = blake3::Hasher::new();
    h.update(b"collagene-worker-bundle-v1\0");
    for (name, data) in FILES.into_iter().zip(files) {
        h.update(&(name.len() as u64).to_le_bytes());
        h.update(name.as_bytes());
        h.update(&(data.len() as u64).to_le_bytes());
        h.update(data);
    }
    h.finalize().to_hex().to_string()
}
fn loader_valid(bytes: &[u8]) -> bool {
    let Ok(text) = std::str::from_utf8(bytes) else {
        return false;
    };
    let mut rest = text;
    fn token(rest: &mut &str, s: &str) -> bool {
        *rest = rest.trim_start_matches(|c: char| c.is_ascii_whitespace());
        if let Some(t) = rest.strip_prefix(s) {
            *rest = t;
            true
        } else {
            false
        }
    }
    fn reference(rest: &mut &str, name: &str) -> bool {
        *rest = rest.trim_start_matches(|c: char| c.is_ascii_whitespace());
        let Some(q) = rest.chars().next().filter(|q| matches!(q, '\'' | '"')) else {
            return false;
        };
        *rest = &rest[1..];
        let Some(end) = rest.find(q) else {
            return false;
        };
        let value = &rest[..end];
        if value != name && value != format!("./{name}") {
            return false;
        }
        *rest = &rest[end + 1..];
        true
    }
    token(&mut rest, "importScripts")
        && token(&mut rest, "(")
        && reference(&mut rest, FILES[1])
        && token(&mut rest, ")")
        && token(&mut rest, ";")
        && token(&mut rest, "wasm_bindgen")
        && token(&mut rest, "(")
        && reference(&mut rest, FILES[2])
        && token(&mut rest, ")")
        && token(&mut rest, ";")
        && rest.trim().is_empty()
}
fn prefix(value: &str) -> Result<String> {
    if value.len() > 8192
        || value.bytes().any(|b| {
            b.is_ascii_control() || b.is_ascii_whitespace() || matches!(b, b'?' | b'#' | b'\\')
        })
    {
        return Err(Error::InvalidPublicUrl);
    }
    let value = if value.is_empty() { "/" } else { value };
    let path = if let Some(rest) = value
        .strip_prefix("https://")
        .or_else(|| value.strip_prefix("http://"))
    {
        let (authority, path) = rest.split_once('/').unwrap_or((rest, ""));
        if authority.is_empty()
            || authority.contains('@')
            || authority.contains(['\'', '"', '<', '>'])
        {
            return Err(Error::InvalidPublicUrl);
        }
        path
    } else {
        if value.starts_with("//") || value.contains(':') {
            return Err(Error::InvalidPublicUrl);
        }
        value
    };
    if path
        .split('/')
        .any(|p| matches!(p, ".." | "%2e%2e" | "%2E%2E"))
    {
        return Err(Error::InvalidPublicUrl);
    }
    Ok(if value.ends_with('/') {
        value.into()
    } else {
        format!("{value}/")
    })
}
/// Escape one HTML attribute value, including both quote styles; no JS string is emitted.
pub fn escape_attribute(value: &str) -> String {
    let mut s = String::with_capacity(value.len());
    for c in value.chars() {
        match c {
            '&' => s.push_str("&amp;"),
            '"' => s.push_str("&quot;"),
            '\'' => s.push_str("&#39;"),
            '<' => s.push_str("&lt;"),
            '>' => s.push_str("&gt;"),
            _ => s.push(c),
        }
    }
    s
}
struct Meta {
    start: usize,
    end: usize,
    quote: u8,
}
fn tag_end(text: &[u8], start: usize) -> Option<usize> {
    let mut quote = None;
    for (i, &b) in text.iter().enumerate().skip(start) {
        match quote {
            Some(q) if b == q => quote = None,
            None if matches!(b, b'\'' | b'"') => quote = Some(b),
            None if b == b'>' => return Some(i),
            _ => (),
        }
    }
    None
}
fn meta_attribute(html: &str) -> Result<Meta> {
    let lower = html.to_ascii_lowercase();
    let bytes = html.as_bytes();
    let mut cursor = 0;
    let mut found = None;
    while let Some(rel) = html[cursor..].find('<') {
        let start = cursor + rel;
        if html[start..].starts_with("<!--") {
            let end = html[start + 4..].find("-->").ok_or(Error::InvalidHtml)?;
            cursor = start + 4 + end + 3;
            continue;
        }
        let end = tag_end(bytes, start + 1).ok_or(Error::InvalidHtml)?;
        let mut p = start + 1;
        let closing = bytes.get(p) == Some(&b'/');
        if closing {
            p += 1;
        }
        let name_start = p;
        while p < end && bytes[p].is_ascii_alphanumeric() {
            p += 1;
        }
        let name = &lower[name_start..p];
        if !closing && matches!(name, "script" | "style") {
            let close = format!("</{name}");
            let at = lower[end + 1..].find(&close).ok_or(Error::InvalidHtml)? + end + 1;
            cursor = tag_end(bytes, at + close.len()).ok_or(Error::InvalidHtml)? + 1;
            continue;
        }
        if !closing && name == "meta" {
            if end - start > 64 * 1024 {
                return Err(Error::InvalidHtml);
            }
            let mut attrs = BTreeSet::new();
            let (mut worker, mut content) = (false, None);
            while p < end {
                while p < end && (bytes[p].is_ascii_whitespace() || bytes[p] == b'/') {
                    p += 1;
                }
                if p == end {
                    break;
                }
                let begin = p;
                while p < end
                    && (bytes[p].is_ascii_alphanumeric() || matches!(bytes[p], b'-' | b'_' | b':'))
                {
                    p += 1;
                }
                if p == begin || p - begin > 128 || attrs.len() >= 64 {
                    return Err(Error::InvalidHtml);
                }
                let key = &lower[begin..p];
                if !attrs.insert(key) {
                    return Err(Error::InvalidHtml);
                }
                while p < end && bytes[p].is_ascii_whitespace() {
                    p += 1;
                }
                if bytes.get(p) != Some(&b'=') {
                    continue;
                }
                p += 1;
                while p < end && bytes[p].is_ascii_whitespace() {
                    p += 1;
                }
                let q = bytes.get(p).copied().filter(|b| matches!(b, b'\'' | b'"'));
                let a;
                if let Some(q) = q {
                    p += 1;
                    a = p;
                    while p < end && bytes[p] != q {
                        p += 1;
                    }
                    if p == end {
                        return Err(Error::InvalidHtml);
                    }
                } else {
                    a = p;
                    while p < end && !bytes[p].is_ascii_whitespace() && bytes[p] != b'>' {
                        p += 1;
                    }
                }
                let value = &html[a..p];
                if key == "name" && value == "collagene-worker" {
                    worker = true;
                }
                if key == "content" {
                    let quote = q.ok_or(Error::InvalidHtml)?;
                    content = Some(Meta {
                        start: a,
                        end: p,
                        quote,
                    });
                }
                if q.is_some() {
                    p += 1;
                }
            }
            if worker {
                if found.is_some() {
                    return Err(Error::InvalidHtml);
                }
                found = Some(content.ok_or(Error::InvalidHtml)?);
            }
        }
        cursor = end + 1;
    }
    found.ok_or(Error::InvalidHtml)
}
fn prior_marker(html: &str) -> Option<&str> {
    let suffix = html.trim_end();
    let start = suffix.rfind(MARKER)?;
    if suffix[..start].contains(MARKER) {
        return None;
    }
    let marker = &suffix[start + MARKER.len()..];
    let hash = marker.strip_suffix(" -->")?;
    hex(hash).then_some(hash)
}
fn index_bytes(html: &str, url: &str, hash: &str) -> Result<Vec<u8>> {
    let meta = meta_attribute(html)?;
    let count = html.matches(SENTINEL).count();
    if count > 1 {
        return Err(Error::AmbiguousSentinel);
    }
    let escaped = escape_attribute(url);
    let mut base = html.to_owned();
    if count == 1 {
        if &html[meta.start..meta.end] != SENTINEL {
            return Err(Error::InvalidHtml);
        }
        if html.contains(MARKER) {
            return Err(Error::InvalidHtml);
        }
    } else {
        let old = prior_marker(html).ok_or(Error::MissingSentinel)?;
        let expected = format!("workers/{old}/{}", FILES[0]);
        if !html[meta.start..meta.end].ends_with(&expected) {
            return Err(Error::MissingSentinel);
        }
        let at = html.rfind(MARKER).ok_or(Error::MissingSentinel)?;
        base.truncate(at);
        base = base.trim_end().to_owned();
    }
    if !matches!(meta.quote, b'\'' | b'"') {
        return Err(Error::InvalidHtml);
    }
    base.replace_range(meta.start..meta.end, &escaped);
    let out = format!("{}\n{MARKER}{hash} -->\n", base.trim_end());
    if out.len() as u64 > MAX_HTML_BYTES {
        return Err(Error::TooLarge("index.html"));
    }
    Ok(out.into_bytes())
}
struct Temp {
    path: PathBuf,
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = if self.path.is_dir() {
            fs::remove_dir_all(&self.path)
        } else {
            fs::remove_file(&self.path)
        };
    }
}
fn unique(parent: &Path, label: &str, directory: bool) -> Result<Temp> {
    for _ in 0..32 {
        let n = SERIAL.fetch_add(1, Ordering::Relaxed);
        let path = parent.join(format!(".collagene-{label}-{}-{n}", std::process::id()));
        let result = if directory {
            fs::create_dir(&path)
        } else {
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
                .map(drop)
        };
        match result {
            Ok(()) => return Ok(Temp { path }),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e.into()),
        }
    }
    Err(Error::CorruptBundle)
}
fn verify_bundle(dir: &Path, data: &[Vec<u8>; 3]) -> Result<()> {
    directory(dir)?;
    let mut names = BTreeSet::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let name = entry
            .file_name()
            .into_string()
            .map_err(|_| Error::CorruptBundle)?;
        if !FILES.contains(&name.as_str()) || !names.insert(name) {
            return Err(Error::CorruptBundle);
        }
    }
    if names.len() != 3 {
        return Err(Error::CorruptBundle);
    }
    for (i, name) in FILES.iter().enumerate() {
        let actual = read_file(&dir.join(name), data[i].len() as u64, name)
            .map_err(|_| Error::CorruptBundle)?;
        if actual != data[i] {
            return Err(Error::CorruptBundle);
        }
    }
    Ok(())
}
/// Publish a complete immutable three-file bundle, then atomically replace staged
/// HTML. Existing root copies and older bundle directories are retained unchanged.
/// Run after Trunk staging completes, in a build-owned directory without concurrent
/// hostile filesystem mutation; std filesystem checks do not provide openat-style isolation.
pub fn package_worker(staging: &Path, public_url: &str) -> Result<Manifest> {
    let staging_meta = fs::symlink_metadata(staging).map_err(|_| Error::InvalidStaging)?;
    if reparse(&staging_meta) || !staging_meta.is_dir() {
        return Err(Error::InvalidStaging);
    }
    let root = fs::canonicalize(staging)?;
    directory(&root)?;
    let prefix = prefix(public_url)?;
    let html = read_file(&root.join("index.html"), MAX_HTML_BYTES, "index.html")?;
    let html_text = std::str::from_utf8(&html).map_err(|_| Error::InvalidHtml)?;
    let data = [
        read_file(&root.join(FILES[0]), MAX_LOADER_BYTES, FILES[0])?,
        read_file(&root.join(FILES[1]), MAX_JS_BYTES, FILES[1])?,
        read_file(&root.join(FILES[2]), MAX_WASM_BYTES, FILES[2])?,
    ];
    if data.iter().map(|d| d.len() as u64).sum::<u64>() > MAX_TOTAL_BYTES {
        return Err(Error::TooLarge("complete worker bundle"));
    }
    if !loader_valid(&data[0]) {
        return Err(Error::InvalidLoader);
    }
    let glue = std::str::from_utf8(&data[1]).map_err(|_| Error::InvalidGlue)?;
    if glue.contains('\0')
        || !glue.trim_start().starts_with("let wasm_bindgen")
        || !glue.contains("WebAssembly")
    {
        return Err(Error::InvalidGlue);
    }
    if !data[2].starts_with(b"\0asm\x01\0\0\0") {
        return Err(Error::InvalidWasm);
    }
    let hash = worker_hash([&data[0], &data[1], &data[2]]);
    let relative = format!("workers/{hash}/{}", FILES[0]);
    let url = format!("{prefix}{relative}");
    if url.len() > MAX_WORKER_URL_BYTES {
        return Err(Error::InvalidPublicUrl);
    }
    let next_html = index_bytes(html_text, &url, &hash)?;
    let parent = root.join("workers");
    match fs::create_dir(&parent) {
        Ok(()) => (),
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => directory(&parent)?,
        Err(e) => return Err(e.into()),
    }
    let dest = parent.join(&hash);
    if fs::symlink_metadata(&dest).is_ok() {
        verify_bundle(&dest, &data)?;
    } else {
        let temp = unique(&parent, "worker", true)?;
        for (i, name) in FILES.iter().enumerate() {
            let mut f = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(temp.path.join(name))?;
            f.write_all(&data[i])?;
            f.sync_all()?;
        }
        match fs::rename(&temp.path, &dest) {
            Ok(()) => (),
            Err(e) => {
                if fs::symlink_metadata(&dest).is_ok() {
                    verify_bundle(&dest, &data)?;
                } else {
                    return Err(e.into());
                }
            }
        }
        verify_bundle(&dest, &data)?;
    }
    if read_file(&root.join("index.html"), MAX_HTML_BYTES, "index.html")? != html {
        return Err(Error::ChangedInput);
    }
    if next_html != html {
        let temp = unique(&root, "index", false)?;
        {
            let mut f = OpenOptions::new()
                .write(true)
                .truncate(true)
                .open(&temp.path)?;
            f.write_all(&next_html)?;
            f.sync_all()?;
        }
        fs::rename(&temp.path, root.join("index.html"))?;
    }
    let files = std::array::from_fn(|i| Asset {
        name: FILES[i],
        bytes: data[i].len() as u64,
        blake3: blake3::hash(&data[i]).to_hex().to_string(),
    });
    Ok(Manifest {
        schema_version: 1,
        bundle_hash: hash,
        loader_url: url,
        relative_loader_path: relative,
        files,
    })
}

/// Remove only the optional, non-executable WebAssembly `name` custom section.
/// Core sections, exports, split metadata and target features remain byte-exact.
/// This keeps non-LTO development previews small without optimizing their code.
pub fn strip_debug_names(staging: &Path) -> Result<Vec<(String, u64, u64)>> {
    use base64::Engine;
    use sha2::{Digest, Sha384};
    use std::io::{Seek, SeekFrom};
    directory(staging)?;
    let root = fs::canonicalize(staging)?;
    let mut stripped = Vec::new();
    let mut integrity = Vec::new();
    for (index, entry) in fs::read_dir(&root)?.enumerate() {
        if index >= 4096 {
            return Err(Error::TooLarge("staging entries"));
        }
        let entry = entry?;
        let name = entry.file_name();
        let Some(name) = name.to_str().filter(|s| s.ends_with(".wasm")) else {
            continue;
        };
        if !name
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
        {
            return Err(Error::UnsafePath);
        }
        let path = entry.path();
        let before = fs::symlink_metadata(&path)?;
        if reparse(&before) || !before.is_file() {
            return Err(Error::UnsafePath);
        }
        if before.len() > 1024 * 1024 * 1024 {
            return Err(Error::TooLarge("compiled WebAssembly"));
        }
        let mut input = File::open(&path)?;
        let mut header = [0u8; 8];
        input.read_exact(&mut header)?;
        if &header != b"\0asm\x01\0\0\0" {
            return Err(Error::InvalidWasm);
        }
        let mut remove = Vec::new();
        let mut sections = 0usize;
        loop {
            let start = input.stream_position()?;
            if start == before.len() {
                break;
            }
            sections += 1;
            if sections > 10_000 {
                return Err(Error::InvalidWasm);
            }
            let mut id = [0u8; 1];
            input.read_exact(&mut id)?;
            let length = wasm_u32(&mut input)? as u64;
            let end = input
                .stream_position()?
                .checked_add(length)
                .filter(|n| *n <= before.len())
                .ok_or(Error::InvalidWasm)?;
            if id[0] == 0 {
                let name_length = wasm_u32(&mut input)? as u64;
                let at = input.stream_position()?;
                if name_length > end.saturating_sub(at) || at > end {
                    return Err(Error::InvalidWasm);
                }
                if name_length == 4 {
                    let mut custom = [0u8; 4];
                    input.read_exact(&mut custom)?;
                    if &custom == b"name" {
                        remove.push((start, end));
                    }
                }
            }
            input.seek(SeekFrom::Start(end))?;
        }
        if remove.is_empty() {
            continue;
        }
        let temp = unique(&root, "wasm-names", false)?;
        let mut output = OpenOptions::new()
            .write(true)
            .truncate(true)
            .open(&temp.path)?;
        input.seek(SeekFrom::Start(0))?;
        let mut at = 0;
        for (start, end) in remove {
            copy_wasm_bytes(&mut input, &mut output, start - at)?;
            input.seek(SeekFrom::Start(end))?;
            at = end;
        }
        copy_wasm_bytes(&mut input, &mut output, before.len() - at)?;
        output.sync_all()?;
        let after = output.metadata()?.len();
        drop(output);
        drop(input);
        let unchanged = fs::symlink_metadata(&path)?;
        if reparse(&unchanged)
            || unchanged.len() != before.len()
            || unchanged.modified().ok() != before.modified().ok()
        {
            return Err(Error::ChangedInput);
        }
        fs::rename(&temp.path, &path)?;
        let mut digest = Sha384::new();
        let mut file = File::open(&path)?;
        let mut buffer = [0u8; 64 * 1024];
        loop {
            let n = file.read(&mut buffer)?;
            if n == 0 {
                break;
            }
            digest.update(&buffer[..n]);
        }
        integrity.push((
            name.to_owned(),
            format!(
                "sha384-{}",
                base64::engine::general_purpose::STANDARD.encode(digest.finalize())
            ),
        ));
        stripped.push((name.to_owned(), before.len(), after));
    }
    if !integrity.is_empty() {
        let path = root.join("index.html");
        let bytes = read_file(&path, MAX_HTML_BYTES, "index.html")?;
        let html = std::str::from_utf8(&bytes).map_err(|_| Error::InvalidHtml)?;
        let mut replacements = Vec::new();
        let mut cursor = 0;
        while let Some(offset) = html[cursor..].find("<link") {
            let start = cursor + offset;
            let end = tag_end(html.as_bytes(), start + 5).ok_or(Error::InvalidHtml)?;
            let tag = &html[start..=end];
            let attributes = link_attributes(tag)?;
            if let (Some(href), Some(sri)) = (
                attributes.iter().find(|a| a.0 == "href"),
                attributes.iter().find(|a| a.0 == "integrity"),
            ) {
                if let Some((_, hash)) = integrity
                    .iter()
                    .find(|(name, _)| tag[href.1..href.2].ends_with(name))
                {
                    replacements.push((start + sri.1, start + sri.2, hash.clone()));
                }
            }
            cursor = end + 1;
        }
        let mut html = html.to_owned();
        for (start, end, value) in replacements.into_iter().rev() {
            html.replace_range(start..end, &value);
        }
        if html.as_bytes() != bytes {
            let temp = unique(&root, "wasm-integrity", false)?;
            {
                let mut file = OpenOptions::new()
                    .write(true)
                    .truncate(true)
                    .open(&temp.path)?;
                file.write_all(html.as_bytes())?;
                file.sync_all()?;
            }
            fs::rename(&temp.path, path)?;
        }
    }
    Ok(stripped)
}
fn wasm_u32(reader: &mut impl Read) -> Result<u32> {
    let mut value = 0u32;
    for shift in (0..35).step_by(7) {
        let mut b = [0u8; 1];
        reader.read_exact(&mut b)?;
        if shift == 28 && b[0] & 0xf0 != 0 {
            return Err(Error::InvalidWasm);
        }
        value |= u32::from(b[0] & 0x7f) << shift;
        if b[0] & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(Error::InvalidWasm)
}
fn copy_wasm_bytes(reader: &mut File, writer: &mut File, length: u64) -> Result<()> {
    if std::io::copy(&mut reader.take(length), writer)? != length {
        return Err(Error::InvalidWasm);
    }
    Ok(())
}
// Trunk emits quoted link attributes. Preserve every value except the matching
// WASM preload integrity value; stylesheet/JavaScript integrity stays untouched.
fn link_attributes(tag: &str) -> Result<Vec<(&str, usize, usize)>> {
    let b = tag.as_bytes();
    let mut at = 5;
    let mut out = Vec::new();
    while at < b.len() {
        while b.get(at).is_some_and(u8::is_ascii_whitespace) {
            at += 1;
        }
        if matches!(b.get(at), Some(b'>') | Some(b'/') | None) {
            break;
        }
        let start = at;
        while b
            .get(at)
            .is_some_and(|v| v.is_ascii_alphanumeric() || b"_:-".contains(v))
        {
            at += 1;
        }
        if at == start {
            return Err(Error::InvalidHtml);
        }
        let name = &tag[start..at];
        while b.get(at).is_some_and(u8::is_ascii_whitespace) {
            at += 1;
        }
        if b.get(at) != Some(&b'=') {
            continue;
        }
        at += 1;
        while b.get(at).is_some_and(u8::is_ascii_whitespace) {
            at += 1;
        }
        let quote = b
            .get(at)
            .copied()
            .filter(|v| matches!(v, b'\'' | b'"'))
            .ok_or(Error::InvalidHtml)?;
        at += 1;
        let value = at;
        while b.get(at).is_some_and(|v| *v != quote) {
            at += 1;
        }
        if at >= b.len() {
            return Err(Error::InvalidHtml);
        }
        out.push((name, value, at));
        at += 1;
    }
    Ok(out)
}
