//! Development smoke helper; normal applications call the same embedded SourceClient directly.
use collagene_sources::{Config, SourceClient};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = std::env::args()
        .nth(1)
        .ok_or("supply an accession.version or public HTTPS sequence URL")?;
    let client = SourceClient::new(Config::from_env()?)?;
    let report = client.import_url(&input).await?;
    for m in report.molecules {
        println!(
            "{} · {} bases · {:?} · {}",
            m.name,
            m.sequence.len(),
            m.topology,
            m.metadata
                .get("source:canonical_url")
                .map_or("", String::as_str)
        );
    }
    Ok(())
}
