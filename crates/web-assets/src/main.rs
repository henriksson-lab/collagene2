use std::path::PathBuf;
fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("worker asset packaging failed: {e}");
            std::process::ExitCode::FAILURE
        }
    }
}
fn run() -> Result<(), Box<dyn std::error::Error>> {
    let staging = std::env::var_os("TRUNK_STAGING_DIR").ok_or("TRUNK_STAGING_DIR is required")?;
    let public = std::env::var("TRUNK_PUBLIC_URL").unwrap_or_else(|_| "/".into());
    let staging = PathBuf::from(staging);
    for (name, before, after) in collagene_web_assets::strip_debug_names(&staging)? {
        eprintln!("Stripped WebAssembly debug names: {name}: {before} -> {after} bytes");
    }
    let m = collagene_web_assets::package_worker(&staging, &public)?;
    for asset in &m.files {
        eprintln!(
            "Packaged worker asset {}: {} bytes",
            asset.name, asset.bytes
        );
    }
    println!("{}", m.loader_url);
    Ok(())
}
