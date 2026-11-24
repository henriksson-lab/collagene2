use collagene_enzyme_catalog::{Edition, seal};
use std::io::Read;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("supply an independently curated proposal JSON path")?;
    let mut bytes = Vec::new();
    std::fs::File::open(&path)?
        .take((collagene_enzyme_catalog::MAX_JSON_BYTES + 1) as u64)
        .read_to_end(&mut bytes)?;
    if bytes.len() > collagene_enzyme_catalog::MAX_JSON_BYTES {
        return Err("proposal exceeds limit".into());
    }
    let proposal: Edition = serde_json::from_slice(&bytes)?;
    let sealed = seal(proposal)?;
    std::fs::write(path, serde_json::to_vec(&sealed)?)?;
    Ok(())
}
