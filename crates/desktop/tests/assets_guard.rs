//! Standalone std-only fixtures: rustc --test tests/assets_guard.rs.
//! They do not compile Tauri, instantiate WASM or start a GUI/listener.
#[path = "../asset_compatibility.rs"]
mod asset_compatibility;
use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
const MAIN: &str = "collagene-web-0123456789abcdef_bg.wasm";
const OTHER: &str = "collagene-web-fedcba9876543210_bg.wasm";
struct Assets(PathBuf);
impl Assets {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let root = std::env::temp_dir().join(format!(
            "collagene-desktop-assets-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).unwrap();
        Self(root)
    }
    fn index(&self, bootstrap: &str, preload: &str) -> String {
        let index = format!(
            r#"<!doctype html><html><head><script type="module">import init from '/app.js'; const wasm = await init({{ module_or_path: '/{bootstrap}' }});</script><link rel="preload" href="/{preload}" crossorigin="anonymous" integrity="sha384-untouched" as="fetch" type="application/wasm"></head><body></body></html>"#
        );
        fs::write(self.0.join("index.html"), &index).unwrap();
        index
    }
    fn check(&self) -> std::io::Result<()> {
        asset_compatibility::check_native_bridge(&self.0)
    }
}
impl Drop for Assets {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
#[test]
fn referenced_native_main_is_accepted_without_altering_preload_or_disk_bytes() {
    let assets = Assets::new();
    let index = assets.index(MAIN, MAIN);
    // Put the invoke literal across the streaming guard's read boundary.
    let mut main = vec![0u8; 64 * 1024 + b"collagene_rpc".len() - 5];
    main.extend_from_slice(b"collagene_rpc");
    main.extend_from_slice(b"unchanged trailing frontend data");
    fs::write(assets.0.join(MAIN), &main).unwrap();
    assets.check().unwrap();
    assert_eq!(
        fs::read_to_string(assets.0.join("index.html")).unwrap(),
        index
    );
    assert_eq!(fs::read(assets.0.join(MAIN)).unwrap(), main);
}
#[test]
fn native_literal_only_in_worker_or_unreferenced_stale_main_cannot_enable_desktop_packaging() {
    let assets = Assets::new();
    assets.index(MAIN, MAIN);
    fs::write(assets.0.join(MAIN), b"old browser-only main").unwrap();
    fs::write(assets.0.join(OTHER), b"collagene_rpc").unwrap();
    fs::create_dir(assets.0.join("workers")).unwrap();
    fs::write(
        assets.0.join("workers/analysis-worker_bg.wasm"),
        b"collagene_rpc",
    )
    .unwrap();
    assert!(assets.check().is_err());
    fs::write(assets.0.join(MAIN), b"current native main collagene_rpc").unwrap();
    assets.check().unwrap();
}
#[test]
fn actual_trunk_bootstrap_and_preload_must_reference_one_same_main_wasm() {
    let assets = Assets::new();
    fs::write(assets.0.join(MAIN), b"collagene_rpc").unwrap();
    fs::write(assets.0.join(OTHER), b"collagene_rpc").unwrap();
    assets.index(MAIN, OTHER);
    assert!(assets.check().is_err());
    let first = assets.index(MAIN, MAIN);
    let second = assets.index(OTHER, OTHER);
    fs::write(assets.0.join("index.html"), format!("{first}{second}")).unwrap();
    assert!(assets.check().is_err());
    assets.index(
        "workers/analysis-worker_bg.wasm",
        "workers/analysis-worker_bg.wasm",
    );
    assert!(assets.check().is_err());
    fs::write(
        assets.0.join("index.html"),
        format!(r#"<html><head><link href="/{MAIN}" type="application/wasm"></head></html>"#),
    )
    .unwrap();
    assert!(assets.check().is_err());
}
