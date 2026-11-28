//! Maintainer-only reference-data curation; application design loads a compact profile.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let profile = collagene_core::sequence::codon_design::hosts::HostProfile::derive_mg1655()?;
    profile.validate()?;
    serde_json::to_writer_pretty(std::io::stdout().lock(), &profile)?;
    println!();
    Ok(())
}
