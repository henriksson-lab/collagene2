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
