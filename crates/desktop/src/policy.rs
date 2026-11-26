//! Desktop IPC is available only to the one embedded application document.
pub fn trusted_url(url: &tauri::Url) -> bool {
    #[cfg(target_os = "windows")]
    let origin = url.scheme() == "https" && url.host_str() == Some("tauri.localhost");
    #[cfg(not(target_os = "windows"))]
    let origin = url.scheme() == "tauri" && url.host_str() == Some("localhost");
    origin
        && url.username().is_empty()
        && url.password().is_none()
        && url.port().is_none()
        && url.query().is_none()
        && matches!(url.path(), "" | "/" | "/index.html")
}
/// Do not place a writable repository within a build tree that a later desktop
/// package could embed. Packaged assets themselves are immutable compiled bytes.
pub fn separate_assets(root: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    use std::path::Component;
    let absolute = if root.is_absolute() {
        root.to_path_buf()
    } else {
        std::env::current_dir()?.join(root)
    };
    if absolute
        .components()
        .any(|part| matches!(part, Component::ParentDir))
    {
        return Err("Use a storage root without parent-directory traversal".into());
    }
    let mut ancestor = absolute.as_path();
    let mut suffix = Vec::new();
    while !ancestor.exists() {
        suffix.push(
            ancestor
                .file_name()
                .ok_or("Invalid desktop storage root")?
                .to_owned(),
        );
        ancestor = ancestor.parent().ok_or("Invalid desktop storage root")?;
    }
    let mut candidate = ancestor.canonicalize()?;
    for part in suffix.into_iter().rev() {
        candidate.push(part);
    }
    if let Ok(assets) = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../dist")
        .canonicalize()
    {
        if candidate.starts_with(&assets) || assets.starts_with(&candidate) {
            return Err("Desktop storage root must be separate from frontend assets".into());
        }
    }
    Ok(())
}
