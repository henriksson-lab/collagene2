use std::{fs, path::Path};
mod asset_compatibility;
fn main() {
    println!("cargo:rerun-if-changed=tauri.conf.json");
    println!("cargo:rerun-if-changed=capabilities");
    println!("cargo:rerun-if-changed=../../dist");
    // Package a fresh frontend, never the historical public-server worker archive.
    // No generated files or frontend bytes are changed by this check.
    let frontend = Path::new("../../dist");
    let mut count = 0;
    let mut total = 0;
    check_assets(frontend, 0, &mut count, &mut total)
        .expect("prepare a fresh confined ../../dist before building Collagene desktop");
    assert!(
        frontend.join("index.html").is_file(),
        "desktop frontend requires dist/index.html"
    );
    asset_compatibility::check_native_bridge(frontend).expect("desktop requires a fresh main Leptos WASM with its native collagene_rpc bridge; browser-only assets cannot provide durable storage in an ephemeral desktop window");
    if let Ok(workers) = fs::read_dir(frontend.join("workers")) {
        let versions = workers
            .collect::<Result<Vec<_>, _>>()
            .expect("read current worker directory");
        assert!(
            versions.len() <= 1,
            "desktop dist must contain only the current worker version; use a fresh Trunk output rather than the public history snapshot"
        );
    }
    tauri_build::try_build(
        tauri_build::Attributes::new()
            .app_manifest(tauri_build::AppManifest::new().commands(&["collagene_rpc"])),
    )
    .expect("build the restricted Collagene desktop manifest");
}
fn check_assets(
    path: &Path,
    depth: usize,
    count: &mut usize,
    total: &mut u64,
) -> std::io::Result<()> {
    if depth > 16 {
        return Err(std::io::Error::other("frontend depth limit"));
    }
    let metadata = fs::symlink_metadata(path)?;
    if metadata.file_type().is_symlink() {
        return Err(std::io::Error::other("frontend symlink forbidden"));
    }
    if metadata.is_dir() {
        for child in fs::read_dir(path)? {
            check_assets(&child?.path(), depth + 1, count, total)?;
        }
    } else if metadata.is_file() {
        *count += 1;
        *total = total
            .checked_add(metadata.len())
            .ok_or_else(|| std::io::Error::other("frontend size overflow"))?;
        if *count > 512 || metadata.len() > 128 * 1024 * 1024 || *total > 256 * 1024 * 1024 {
            return Err(std::io::Error::other("frontend asset limit"));
        }
    } else {
        return Err(std::io::Error::other("frontend special file forbidden"));
    }
    Ok(())
}
