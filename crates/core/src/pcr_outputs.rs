//! Explicit selection of freshly reconstructed exact-PCR reference outputs.
//! These references and order tables assert neither an independent target check
//! nor physical terminal chemistry. Hosts bind them to exact source history events.
use crate::{
    analysis::{PcrProduct, Primer},
    oligo_formats::OligoTable,
    verification,
};
use collagene_model::{Alphabet, Error, Molecule};
use serde::{Deserialize, Serialize};

pub const ENGINE_VERSION: &str = "pcr-outputs-1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Options {pub name:String,pub save_product:bool,pub save_primers:bool}}
dto! {pub struct Request {pub source:Molecule,pub forward:Primer,pub reverse:Primer,pub template_strand:crate::chemistry::Strand,pub selected_index:usize,pub expected_product_revision:String,pub options:Options}}
dto! {pub struct Prepared {
    pub context_revision:String,
    pub selected_index:usize,
    #[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]
    pub geometry:PcrProduct,
    pub product:Option<Molecule>,
    #[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]
    pub primers:Option<OligoTable>,
    pub findings:Vec<String>,
}}
fn bad(message: &str) -> Error {
    Error::Format(format!("PCR outputs: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn name(value: &str) -> Result<(), Error> {
    if value.trim().is_empty() || value.len() > 256 || value.chars().any(char::is_control) {
        return Err(bad(
            "names must be nonempty, at most 256 bytes, and contain no control characters",
        ));
    }
    Ok(())
}
pub(crate) fn primer(value: &Primer) -> Result<(), Error> {
    name(&value.name)?;
    if value.annealing.alphabet() != Alphabet::Dna
        || value.annealing.is_empty()
        || !value
            .annealing
            .as_str()
            .bytes()
            .all(|b| b"ACGT".contains(&b))
    {
        return Err(bad("exact-PCR annealing domains require concrete DNA"));
    }
    // The existing geometry engine supports IUPAC tails as reference bases;
    // uncertainty in a synthetic tail does not create ambiguous template binding.
    if !value.tail.trim().is_empty() {
        collagene_model::Sequence::dna(&value.tail)?;
    }
    Ok(())
}
impl Options {
    pub fn validate(&self) -> Result<(), Error> {
        name(&self.name)?;
        if self.name.trim() != self.name {
            return Err(bad("output name must be trimmed"));
        }
        if !self.save_product && !self.save_primers {
            return Err(bad("choose at least one output"));
        }
        Ok(())
    }
}
/// Bounded identity over all source content, oligos, selected preview and choices.
/// Source history events remain the host's separate binding responsibility.
pub fn revision(request: &Request) -> Result<String, Error> {
    collagene_model::bounds::model_preflight(&[&request.source])?;
    request.source.validate()?;
    primer(&request.forward)?;
    primer(&request.reverse)?;
    request.options.validate()?;
    if request.selected_index >= 10_000 || !digest(&request.expected_product_revision) {
        return Err(bad(
            "invalid selected geometry or expected preview revision",
        ));
    }
    crate::assembly_recipe::hash(&(ENGINE_VERSION, request))
}
/// Two literal named oligos, including domain/tail decomposition and full ordered
/// bases. Reuses the standard collection constructor and its export safety rules.
pub fn primer_table(forward: &Primer, reverse: &Primer) -> Result<OligoTable, Error> {
    primer(forward)?;
    primer(reverse)?;
    crate::oligo_orders::from_primers(&[forward.clone(), reverse.clone()])
}
/// Recompute before selecting; serialized previews and user-facing indices are
/// never accepted as proof of a source/primer geometry. Run in an analysis worker.
pub fn prepare(request: &Request) -> Result<Prepared, Error> {
    let context_revision = revision(request)?;
    let previews = verification::simulate_pcr_documents(
        &request.source,
        &request.forward,
        &request.reverse,
        request.template_strand,
    )?;
    prepare_from_previews(request, context_revision, previews)
}
/// Only native callers which just reconstructed source-bound previews may use
/// this selection adapter. No serialized preview is accepted by public APIs.
pub(crate) fn prepare_from_previews(
    request: &Request,
    context_revision: String,
    previews: Vec<verification::PcrDocument>,
) -> Result<Prepared, Error> {
    let preview = previews
        .into_iter()
        .nth(request.selected_index)
        .ok_or_else(|| bad("selected PCR geometry is no longer available"))?;
    if preview.molecule.revision() != request.expected_product_revision {
        return Err(bad(
            "selected preview changed; recompute and select its geometry again",
        ));
    }
    let mut product = preview.molecule;
    product.name = request.options.name.clone();
    let primers = if request.options.save_primers {
        Some(primer_table(&request.forward, &request.reverse)?)
    } else {
        None
    };
    let mut findings = preview.findings;
    if primers
        .as_ref()
        .is_some_and(|table| table.records.iter().any(|row| row.primer.is_none()))
    {
        findings.push("Some declared oligos exceed the existing reusable-primer import profile. Their full ordered bases remain collection/order evidence; inspect row diagnostics before reuse.".into());
    }
    let prepared = Prepared {
        context_revision,
        selected_index: request.selected_index,
        geometry: preview.geometry,
        product: request.options.save_product.then_some(product),
        primers,
        findings,
    };
    prepared.validate_binding(&prepared.context_revision)?;
    Ok(prepared)
}
impl Prepared {
    /// Read-only structural guard. Fresh native reconstruction and host-owned
    /// exact event/epoch checks supply authority, not a checksum alone.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !digest(expected)
            || self.context_revision != expected
            || self.selected_index >= 10_000
            || (self.product.is_none() && self.primers.is_none())
        {
            return Err(bad("output receipt does not match its request"));
        }
        if let Some(product) = &self.product {
            collagene_model::bounds::model_preflight(&[product])?;
            product.validate()?;
            name(&product.name)?;
            if product.sequence != self.geometry.sequence
                || product.topology != collagene_model::Topology::Linear
            {
                return Err(bad("output sequence does not match selected geometry"));
            }
        }
        if let Some(table) = &self.primers {
            table.validate()?;
            if table.records.len() != 2
                || table.records.iter().any(|r| {
                    r.raw_ordered_sequence.is_none()
                        || r.raw_annealing_sequence.is_none()
                        || r.raw_tail.is_none()
                })
            {
                return Err(bad("output collection must contain both declared primers"));
            }
        }
        crate::assembly_recipe::hash(self)?;
        Ok(())
    }
}
