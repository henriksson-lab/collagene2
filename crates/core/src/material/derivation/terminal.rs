//! Terminal processing preview and retained source/chemistry bindings.
use super::ENGINE_VERSION;
use super::FreshMaterial;
use super::SourceOrigin;
use super::TerminalAction;
use super::bad;
use super::encoded;
use super::fresh_source;
use super::hex;
use super::limitations;
use super::summary_shape;
use super::workspace_source;
use crate::chemistry;
use crate::chemistry::tracked_duplex as modified_duplex;
use crate::chemistry::tracked_duplex::ModificationSummary;
use crate::chemistry::tracked_duplex::TrackedDuplex;
use crate::documents::workspace::Workspace;
use crate::material::lineage::EventBinding;
use crate::material::run as run_material;
use crate::verification;
use collagene_model::Alphabet;
use collagene_model::Error;
use collagene_model::Molecule;
use serde::Deserialize;
use serde::Serialize;

// Documentary source/chemistry comparison; serialized summaries never restore state.
dto! {pub struct TerminalProcessingPreview{
pub engine_version:String,pub context_revision:String,pub source_revision:String,pub source_sequence_blake3:String,
pub source:EventBinding,pub origin_revision:Option<String>,pub origin_incomplete:bool,pub action:TerminalAction,
pub product:Molecule,pub source_duplex:chemistry::DuplexSummary,pub duplex:chemistry::DuplexSummary,
pub source_modifications:ModificationSummary,pub modifications:ModificationSummary,
pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
fn source_sequence_hash(source: &Molecule) -> String {
    blake3::hash(source.sequence.as_str().as_bytes())
        .to_hex()
        .to_string()
}
fn terminal_context_fields(
    source: &EventBinding,
    source_sequence_blake3: &str,
    origin_revision: &Option<String>,
    action: &TerminalAction,
) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        "terminal-processing-preview",
        verification::ENGINE_VERSION,
        source,
        source_sequence_blake3,
        origin_revision,
        action,
    ))
}
pub fn terminal_processing_preview_revision(
    w: &Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
    action: &TerminalAction,
) -> Result<String, Error> {
    let molecule = workspace_source(w, source, origin)?;
    terminal_context_fields(
        source,
        &source_sequence_hash(molecule),
        &match origin {
            SourceOrigin::OriginalReference => None,
            SourceOrigin::Saved(origin) => Some(crate::encoding::bounded_json::hash(origin)?),
        },
        action,
    )
}
pub(super) fn terminal_from_fresh(
    source: &Molecule,
    fresh: &FreshMaterial,
    action: &TerminalAction,
    name: &str,
) -> Result<(Molecule, TrackedDuplex, Vec<String>), Error> {
    let (tracked, findings) = match action {
        TerminalAction::TreatEnds { treatment } => {
            let tracked = fresh.tracked.treat(*treatment)?;
            (
                tracked,
                vec![format!(
                    "Explicit {treatment:?} changes only the two external 5′ groups on actual intact linear strands. Strand letters, 3′ groups and positional modifications remain; no enzyme yield, nick repair or independent-target success is inferred."
                )],
            )
        }
        TerminalAction::Tail3Prime { base } => {
            let tracked = fresh.tracked.tail(*base)?;
            (tracked,vec![format!("Exactly one explicitly declared {} base is added to each actual 3′ hydroxyl on the blunt linear duplex. Source calls, 5′ groups, existing nicks and modifications remain; new calls are unmodified under the existing standard-nucleotide model.",base.as_str()),"This ideal endpoint transform does not predict enzyme/kit activity, tail heterogeneity, yield or experimental completion. Symmetric single-base tails do not establish insertion orientation or unique products; later ligation must retain alternative products and the independent intended path.".into()])
        }
    };
    let mut reference = tracked.duplex.molecule(name)?;
    // Source top-reference calls remain in precisely the same coordinates. Appended
    // terminal bases are outside every original feature, never inserted into a CDS.
    reference.features = source.features.clone();
    reference.metadata.insert(
        "source:material_derivation:terminal_source_revision".into(),
        source.revision(),
    );
    reference.metadata.insert(
        "source:material_derivation:terminal_source_sequence_blake3".into(),
        source_sequence_hash(source),
    );
    reference.metadata.insert(
        "source:material_derivation:terminal_action".into(),
        serde_json::to_string(action).map_err(|_| Error::TooLarge)?,
    );
    reference.metadata.insert("source:material_derivation:terminal_scope".into(),"Native source-backed endpoint processing. Detached reference and summaries are not reusable stock; retain exact source events and replay the native recipe. No intended target was checked.".into());
    reference.metadata.insert("source:material_derivation:scope".into(),"Terminal processing reference preview. Exact native source ancestry and processing recipe are required for material reuse.".into());
    reference.validate()?;
    Ok((reference, tracked, findings))
}
pub fn preview_terminal_processing(
    w: &Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
    action: &TerminalAction,
) -> Result<TerminalProcessingPreview, Error> {
    let context_revision = terminal_processing_preview_revision(w, source, origin, action)?;
    let molecule = run_material::resolve(w, source)?;
    let fresh = fresh_source(molecule, origin)?;
    let (product, tracked, local) =
        terminal_from_fresh(molecule, &fresh, action, "Terminally processed material")?;
    let total = fresh
        .findings_total
        .checked_add(local.len())
        .ok_or(Error::TooLarge)?;
    let (findings, findings_total, findings_truncated) = limitations(
        local
            .into_iter()
            .chain(fresh.findings.iter().cloned())
            .collect(),
        total,
        fresh.findings_truncated,
    );
    let result = TerminalProcessingPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: molecule.revision(),
        source_sequence_blake3: source_sequence_hash(molecule),
        source: source.clone(),
        origin_revision: match origin {
            SourceOrigin::OriginalReference => None,
            SourceOrigin::Saved(origin) => Some(crate::encoding::bounded_json::hash(origin)?),
        },
        origin_incomplete: fresh.incomplete,
        action: *action,
        product,
        source_duplex: fresh.tracked.duplex.summary()?,
        duplex: tracked.duplex.summary()?,
        source_modifications: modified_duplex::summary(&fresh.tracked)?,
        modifications: modified_duplex::summary(&tracked)?,
        incomplete: fresh.incomplete,
        findings,
        findings_total,
        findings_truncated,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl TerminalProcessingPreview {
    /// Bounded documentary correspondence only. Never replays source material or
    /// executes endpoint processing, and cannot authorize installation on its own.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        run_material::event_shape(&self.source)?;
        collagene_model::bounds::model_preflight(&[&self.product])?;
        self.product.validate()?;
        summary_shape(&self.source_modifications)?;
        summary_shape(&self.modifications)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || self.source_revision != self.source.revision
            || !hex(&self.source_sequence_blake3)
            || self.origin_revision.as_ref().is_some_and(|s| !hex(s))
            || self.origin_revision.is_none() && self.origin_incomplete
            || self.context_revision
                != terminal_context_fields(
                    &self.source,
                    &self.source_sequence_blake3,
                    &self.origin_revision,
                    &self.action,
                )?
            || self.incomplete != self.origin_incomplete
            || self.findings.len() != self.findings_total.min(128)
            || self.findings_total > crate::material::lineage::MAX_RUN_BYTES
            || self.findings.iter().any(|s| s.chars().count() > 2048)
            || self.findings_total > self.findings.len() && !self.findings_truncated
            || self.source_duplex.top_length != self.source_modifications.top_length
            || self.source_duplex.bottom_length != self.source_modifications.bottom_length
            || self.duplex.top_length != self.modifications.top_length
            || self.duplex.bottom_length != self.modifications.bottom_length
            || self.product.sequence.alphabet() != Alphabet::Dna
            || self.product.sequence.len() != self.duplex.top_length
            || self.product.topology != self.duplex.topology
            || self.duplex.topology != collagene_model::Topology::Linear
            || self.source_duplex.topology != collagene_model::Topology::Linear
            || self.source_duplex.nicks != self.duplex.nicks
            || self.product.metadata.get("duplex_model")
                != Some(&serde_json::to_string(&self.duplex).map_err(|_| Error::TooLarge)?)
            || self
                .product
                .metadata
                .get("source:material_derivation:terminal_source_revision")
                != Some(&self.source_revision)
            || self
                .product
                .metadata
                .get("source:material_derivation:terminal_source_sequence_blake3")
                != Some(&self.source_sequence_blake3)
            || self
                .product
                .metadata
                .get("source:material_derivation:terminal_action")
                != Some(&serde_json::to_string(&self.action).map_err(|_| Error::TooLarge)?)
            || self.source_modifications.known_modified != self.modifications.known_modified
            || self.source_modifications.unknown_eligible != self.modifications.unknown_eligible
            || self.source_modifications.generation != self.modifications.generation
        {
            return Err(bad(
                "Terminal preview changed its exact source, action, chemistry, inherited limitations or bounded findings",
            ));
        }
        let source_ends = self
            .source_duplex
            .ends
            .as_ref()
            .ok_or_else(|| bad("Terminal processing source has no free ends"))?;
        let product_ends = self
            .duplex
            .ends
            .as_ref()
            .ok_or_else(|| bad("Terminal processing product has no free ends"))?;
        match self.action {
            TerminalAction::TreatEnds { treatment } => {
                if !self.duplex.nicks.is_empty()
                    || self.duplex.top_length != self.source_duplex.top_length
                    || self.duplex.bottom_length != self.source_duplex.bottom_length
                    || self.duplex.bottom_start != self.source_duplex.bottom_start
                    || source_sequence_hash(&self.product) != self.source_sequence_blake3
                    || self.modifications.known_unmodified
                        != self.source_modifications.known_unmodified
                {
                    return Err(bad(
                        "Terminal treatment changed source geometry, letters or modifications",
                    ));
                }
                for (before, after) in [
                    (&source_ends.0, &product_ends.0),
                    (&source_ends.1, &product_ends.1),
                ] {
                    let mut wanted = before.clone();
                    wanted.chemistry.five_prime = match treatment {
                        chemistry::EndTreatment::Phosphorylate => chemistry::FivePrime::Phosphate,
                        chemistry::EndTreatment::Dephosphorylate => chemistry::FivePrime::Hydroxyl,
                    };
                    if &wanted != after {
                        return Err(bad(
                            "Terminal treatment changed more than the declared external 5′ groups",
                        ));
                    }
                }
            }
            TerminalAction::Tail3Prime { base } => {
                let length = self.source_duplex.top_length;
                if self.source_duplex.bottom_start != 0
                    || self.source_duplex.bottom_length != length
                    || self.duplex.bottom_start != -1
                    || self.duplex.top_length != length.checked_add(1).ok_or(Error::TooLarge)?
                    || self.duplex.bottom_length != self.duplex.top_length
                    || self.product.sequence.as_str().as_bytes().last().copied()
                        != base.as_str().as_bytes().first().copied()
                    || blake3::hash(&self.product.sequence.as_str().as_bytes()[..length])
                        .to_hex()
                        .as_str()
                        != self.source_sequence_blake3
                    || self.modifications.known_unmodified
                        != self
                            .source_modifications
                            .known_unmodified
                            .checked_add(
                                if matches!(
                                    base,
                                    crate::cloning::ligation::tailing::TailingBase::A
                                        | crate::cloning::ligation::tailing::TailingBase::C
                                ) {
                                    2
                                } else {
                                    0
                                },
                            )
                            .ok_or(Error::TooLarge)?
                {
                    return Err(bad(
                        "Terminal tailing differs from exactly one declared newly synthesized base per strand",
                    ));
                }
                for (before, after) in [
                    (&source_ends.0, &product_ends.0),
                    (&source_ends.1, &product_ends.1),
                ] {
                    if before.polarity != chemistry::Polarity::Blunt
                        || !before.overhang.is_empty()
                        || before.chemistry.three_prime != chemistry::ThreePrime::Hydroxyl
                        || after.polarity != chemistry::Polarity::ThreePrime
                        || after.overhang != base.as_str()
                        || after.chemistry != before.chemistry
                    {
                        return Err(bad(
                            "Terminal tailing changed inherited groups or did not expose the declared 3′ tails",
                        ));
                    }
                }
            }
        }
        encoded(self, crate::material::lineage::MAX_RUN_BYTES)?;
        Ok(())
    }
}
