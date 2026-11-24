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
    let m = collagene_web_assets::package_worker(&PathBuf::from(staging), &public)?;
    println!("{}", m.loader_url);
    Ok(())
}
