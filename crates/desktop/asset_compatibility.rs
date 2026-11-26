use std::{fs, io::Read, path::Path};
// Compatibility check for the actual Trunk bootstrap target. A marker in a
// worker, stale sibling WASM or arbitrary extra file cannot satisfy this check.
// Native IPC framing/authorization remain separate runtime boundaries.
pub fn check_native_bridge(frontend: &Path) -> std::io::Result<()> {
    let bad = || std::io::Error::other("missing/ambiguous referenced main native bridge");
    let index_path = frontend.join("index.html");
    if fs::metadata(&index_path)?.len() > 2 * 1024 * 1024 {
        return Err(bad());
    }
    let index = fs::read_to_string(index_path)?;
    let mut targets = std::collections::BTreeSet::new();
    for chunk in index.split("module_or_path:").skip(1) {
        let chunk = chunk.trim_start();
        let quote = *chunk.as_bytes().first().ok_or_else(bad)?;
        if !matches!(quote, b'\'' | b'"') {
            return Err(bad());
        }
        let end = chunk[1..].find(quote as char).ok_or_else(bad)? + 1;
        let url = &chunk[1..end];
        let name = url
            .strip_prefix('/')
            .or_else(|| url.strip_prefix("./"))
            .unwrap_or(url);
        let hash = name
            .strip_prefix("collagene-web-")
            .and_then(|s| s.strip_suffix("_bg.wasm"))
            .ok_or_else(bad)?;
        if hash.len() < 8
            || hash.len() > 128
            || !hash
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(bad());
        }
        // Require the same URL in the real preload link emitted by Trunk.
        let quoted = [format!("href=\"{url}\""), format!("href='{url}'")];
        if !index
            .split("<link")
            .skip(1)
            .filter_map(|tag| tag.split('>').next())
            .any(|tag| {
                quoted.iter().any(|href| tag.contains(href)) && tag.contains("application/wasm")
            })
        {
            return Err(bad());
        }
        targets.insert(name.to_owned());
    }
    if targets.len() != 1 {
        return Err(bad());
    }
    let target = frontend.join(targets.into_iter().next().ok_or_else(bad)?);
    let mut wasm = fs::File::open(target)?;
    const MARKER: &[u8] = b"collagene_rpc";
    let mut buffer = vec![0u8; 64 * 1024 + MARKER.len()];
    let mut carry = 0;
    loop {
        let read = wasm.read(&mut buffer[carry..])?;
        let end = carry + read;
        if buffer[..end]
            .windows(MARKER.len())
            .any(|slice| slice == MARKER)
        {
            return Ok(());
        }
        if read == 0 {
            return Err(bad());
        }
        carry = (MARKER.len() - 1).min(end);
        buffer.copy_within(end - carry..end, 0);
    }
}
