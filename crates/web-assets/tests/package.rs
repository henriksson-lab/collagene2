use collagene_web_assets::*;
use std::{
    fs,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};
static SEQ: AtomicU64 = AtomicU64::new(0);
struct Stage(PathBuf);
impl Stage {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "collagene-web-assets-test-{}-{}",
            std::process::id(),
            SEQ.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        let s = Self(p);
        fs::write(s.0.join("index.html"), html()).unwrap();
        fs::write(
            s.0.join(FILES[0]),
            b"importScripts('./analysis-worker.js'); wasm_bindgen('./analysis-worker_bg.wasm');\n",
        )
        .unwrap();
        fs::write(
            s.0.join(FILES[1]),
            b"let wasm_bindgen = function (path) { return WebAssembly.instantiate(path); };\n",
        )
        .unwrap();
        fs::write(s.0.join(FILES[2]), b"\0asm\x01\0\0\0\0\x02\x01a").unwrap();
        s
    }
    fn path(&self) -> &Path {
        &self.0
    }
}
impl Drop for Stage {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn html() -> String {
    format!(
        "<!doctype html><html><head><meta name=\"collagene-worker\" content=\"{SENTINEL}\"></head><body></body></html>\n"
    )
}
#[test]
fn complete_bundle_is_byte_identical_domain_hashed_and_idempotent() {
    let s = Stage::new();
    let a = package_worker(s.path(), "/app/").unwrap();
    let index = fs::read(s.0.join("index.html")).unwrap();
    assert_eq!(a.schema_version, 1);
    assert_eq!(a.bundle_hash.len(), 64);
    assert!(a.loader_url.starts_with("/app/workers/"));
    let files: [Vec<u8>; 3] = std::array::from_fn(|i| fs::read(s.0.join(FILES[i])).unwrap());
    assert_eq!(
        a.bundle_hash,
        worker_hash([&files[0], &files[1], &files[2]])
    );
    for (i, f) in a.files.iter().enumerate() {
        let p = s.0.join("workers").join(&a.bundle_hash).join(f.name);
        assert_eq!(fs::read(p).unwrap(), files[i]);
        assert_eq!(f.bytes, files[i].len() as u64);
        assert_eq!(f.blake3, blake3::hash(&files[i]).to_hex().as_str());
    }
    let b = package_worker(s.path(), "/app/").unwrap();
    assert_eq!(a, b);
    assert_eq!(fs::read(s.0.join("index.html")).unwrap(), index);
    assert!(!String::from_utf8(index).unwrap().contains(SENTINEL));
}
#[test]
fn each_input_mutation_changes_the_entire_bundle_uri_and_retains_old_complete_versions() {
    for (i, name) in FILES.iter().enumerate() {
        let s = Stage::new();
        let first = package_worker(s.path(), "/").unwrap();
        let path = s.0.join(name);
        let mut data = fs::read(&path).unwrap();
        if i == 0 {
            data.insert(0, b' ');
        } else if i == 1 {
            data.extend_from_slice(b"// distinct compiled glue\n");
        } else {
            *data.last_mut().unwrap() = b'b';
        }
        fs::write(path, data).unwrap();
        let next = package_worker(s.path(), "/").unwrap();
        assert_ne!(first.bundle_hash, next.bundle_hash);
        assert_ne!(first.loader_url, next.loader_url);
        for name in FILES {
            assert!(
                s.0.join("workers")
                    .join(&first.bundle_hash)
                    .join(name)
                    .is_file()
            );
            assert!(
                s.0.join("workers")
                    .join(&next.bundle_hash)
                    .join(name)
                    .is_file()
            );
        }
    }
}
#[test]
fn existing_same_hash_corruption_or_extraneous_files_is_never_repaired_in_place() {
    for extra in [false, true] {
        let s = Stage::new();
        let m = package_worker(s.path(), "/").unwrap();
        let index = fs::read(s.0.join("index.html")).unwrap();
        let dir = s.0.join("workers").join(m.bundle_hash);
        if extra {
            fs::write(dir.join("unexpected.js"), b"not in bundle").unwrap();
        } else {
            fs::write(dir.join(FILES[1]), b"changed").unwrap();
        }
        assert!(matches!(
            package_worker(s.path(), "/"),
            Err(Error::CorruptBundle)
        ));
        assert_eq!(fs::read(s.0.join("index.html")).unwrap(), index);
    }
}
#[test]
fn invalid_inputs_and_missing_ambiguous_or_nonattribute_sentinels_fail_before_publication() {
    for text in [
        html().replace(SENTINEL, "missing"),
        format!("{}<!-- {SENTINEL} -->", html()),
        format!(
            "<script>let fake='{SENTINEL}';</script><meta name='collagene-worker' content='missing'>"
        ),
        format!("<!-- <meta name='collagene-worker' content='{SENTINEL}'> -->"),
        format!("<meta name='collagene-worker' content={SENTINEL}>"),
        format!("{}<meta name='collagene-worker' content='other'>", html()),
    ] {
        let s = Stage::new();
        fs::write(s.0.join("index.html"), &text).unwrap();
        assert!(package_worker(s.path(), "/").is_err());
        assert!(!s.0.join("workers").exists());
        assert_eq!(fs::read_to_string(s.0.join("index.html")).unwrap(), text);
    }
    for name in FILES {
        let s = Stage::new();
        fs::remove_file(s.0.join(name)).unwrap();
        assert!(package_worker(s.path(), "/").is_err());
        assert!(!s.0.join("workers").exists());
    }
}
#[test]
fn public_prefix_is_escaped_in_both_quote_styles_and_urls_remain_prefix_relative() {
    let public = "/app/a&\"'<b>/";
    for quote in ['\'', '"'] {
        let s = Stage::new();
        let index =
            format!("<meta content={quote}{SENTINEL}{quote} name={quote}collagene-worker{quote}>");
        fs::write(s.0.join("index.html"), index).unwrap();
        let m = package_worker(s.path(), public).unwrap();
        assert!(m.loader_url.starts_with(public));
        let out = fs::read_to_string(s.0.join("index.html")).unwrap();
        assert!(out.contains("/app/a&amp;&quot;&#39;&lt;b&gt;/workers/"));
        assert!(!out.contains(public));
        assert_eq!(package_worker(s.path(), public).unwrap(), m);
    }
    for public in ["https://cdn.example.org/app", "./", "nested/path/", "/"] {
        let s = Stage::new();
        let m = package_worker(s.path(), public).unwrap();
        assert!(m.loader_url.contains("workers/"));
        assert!(!m.loader_url.contains("//workers/"));
    }
    for public in [
        "javascript:alert(1)",
        "//cdn.example.org/app/",
        "https://user:secret@example.org/",
        "/app/?query=x",
        "/app/#fragment",
        "/../",
        "/app\\other/",
        "https:///broken",
    ] {
        let s = Stage::new();
        assert!(matches!(
            package_worker(s.path(), public),
            Err(Error::InvalidPublicUrl)
        ));
    }
}
#[test]
fn loader_chain_and_glue_wasm_shapes_fail_closed() {
    for loader in [
        "importScripts('https://elsewhere/analysis-worker.js');wasm_bindgen('./analysis-worker_bg.wasm');",
        "importScripts('../analysis-worker.js');wasm_bindgen('./analysis-worker_bg.wasm');",
        "importScripts('./analysis-worker.js');wasm_bindgen('/analysis-worker_bg.wasm');",
        "importScripts('./analysis-worker.js');wasm_bindgen('./analysis-worker_bg.wasm');fetch('extra.js');",
    ] {
        let s = Stage::new();
        fs::write(s.0.join(FILES[0]), loader).unwrap();
        assert!(matches!(
            package_worker(s.path(), "/"),
            Err(Error::InvalidLoader)
        ));
        assert!(!s.0.join("workers").exists());
    }
    let s = Stage::new();
    fs::write(s.0.join(FILES[1]), b"<html>error page</html>").unwrap();
    assert!(matches!(
        package_worker(s.path(), "/"),
        Err(Error::InvalidGlue)
    ));
    let s = Stage::new();
    fs::write(s.0.join(FILES[2]), b"not wasm").unwrap();
    assert!(matches!(
        package_worker(s.path(), "/"),
        Err(Error::InvalidWasm)
    ));
}
#[test]
fn each_file_limit_is_checked_before_reading_unbounded_bytes() {
    for (name, limit) in [
        ("index.html", MAX_HTML_BYTES),
        (FILES[0], MAX_LOADER_BYTES),
        (FILES[1], MAX_JS_BYTES),
        (FILES[2], MAX_WASM_BYTES),
    ] {
        let s = Stage::new();
        fs::OpenOptions::new()
            .write(true)
            .truncate(true)
            .open(s.0.join(name))
            .unwrap()
            .set_len(limit + 1)
            .unwrap();
        assert!(matches!(
            package_worker(s.path(), "/"),
            Err(Error::TooLarge(_))
        ));
        assert!(!s.0.join("workers").exists());
    }
}
#[cfg(unix)]
#[test]
fn symlinked_inputs_output_parent_and_existing_hash_dir_are_never_followed() {
    use std::os::unix::fs::symlink;
    for name in FILES {
        let s = Stage::new();
        let outside = Stage::new();
        fs::remove_file(s.0.join(name)).unwrap();
        symlink(outside.0.join(name), s.0.join(name)).unwrap();
        assert!(matches!(
            package_worker(s.path(), "/"),
            Err(Error::UnsafePath)
        ));
    }
    let s = Stage::new();
    let outside = Stage::new();
    symlink(outside.path(), s.0.join("workers")).unwrap();
    assert!(matches!(
        package_worker(s.path(), "/"),
        Err(Error::UnsafePath)
    ));
    let s = Stage::new();
    let first = package_worker(s.path(), "/").unwrap();
    let dir = s.0.join("workers").join(first.bundle_hash);
    fs::remove_dir_all(&dir).unwrap();
    symlink(outside.path(), &dir).unwrap();
    assert!(package_worker(s.path(), "/").is_err());
    assert!(!outside.0.join("unexpected").exists());
}

#[test]
fn oversized_final_url_fails_before_creating_any_bundle_or_mutating_html() {
    let s = Stage::new();
    let original = fs::read(s.0.join("index.html")).unwrap();
    let prefix = format!("/{}/", "a".repeat(MAX_WORKER_URL_BYTES));
    assert!(matches!(
        package_worker(s.path(), &prefix),
        Err(Error::InvalidPublicUrl)
    ));
    assert_eq!(fs::read(s.0.join("index.html")).unwrap(), original);
    assert!(!s.0.join("workers").exists());
}

#[test]
fn filename_and_length_framing_prevent_ambiguous_concatenation_hashes() {
    assert_ne!(
        worker_hash([b"ab", b"c", b"d"]),
        worker_hash([b"a", b"bc", b"d"])
    );
    assert_ne!(
        worker_hash([b"a", b"b", b"c"]),
        worker_hash([b"b", b"a", b"c"])
    );
}

#[test]
fn cli_uses_only_staged_environment_and_rejects_bad_inputs_without_partial_index_update() {
    let s = Stage::new();
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_collagene-web-assets"))
        .env("TRUNK_STAGING_DIR", s.path())
        .env("TRUNK_PUBLIC_URL", "/nested/")
        .output()
        .unwrap();
    assert!(output.status.success());
    let url = String::from_utf8(output.stdout).unwrap();
    assert!(url.trim().starts_with("/nested/workers/"));
    assert!(
        fs::read_to_string(s.0.join("index.html"))
            .unwrap()
            .contains(url.trim())
    );
    let before = fs::read(s.0.join("index.html")).unwrap();
    let error = std::process::Command::new(env!("CARGO_BIN_EXE_collagene-web-assets"))
        .env("TRUNK_STAGING_DIR", s.path())
        .env("TRUNK_PUBLIC_URL", "javascript:bad")
        .output()
        .unwrap();
    assert!(!error.status.success());
    assert_eq!(fs::read(s.0.join("index.html")).unwrap(), before);
}
