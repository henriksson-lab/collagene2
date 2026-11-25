//! Explicit source-backed preparation declarations over the existing TOPO engine.
//! Compiling a declaration neither activates a stock nor adds unrequested DNA.
use crate::{Error, analysis::Primer, tailing::TailingBase, topo};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io};
use topo::physical::{self as physical, projection};
pub const ENGINE_VERSION: &str = "topo-preparation-declaration-1";
macro_rules! dto{($($tokens:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($tokens)*};}
dto! {#[derive(Copy,Default,PartialEq,Eq)]pub enum PreparationMode{#[default]Ta,DirectionalCacc}}
dto! {pub struct PreparedVectorDeclaration{
    pub declaration_id:String,pub assumption:String,pub strands:Vec<topo::PreparedStrand>,pub pairing:Vec<physical::PairSpan>,
}}
dto! {pub struct ExactPcrDeclaration{pub material_id:String,pub template_role:String,pub forward:Primer,pub reverse:Primer,pub tail:Option<TailingBase>}}
dto! {pub struct PreparationRequest{pub mode:PreparationMode,pub vector:PreparedVectorDeclaration,pub insert:ExactPcrDeclaration,
pub engine_limits:physical::Limits,pub projection_options:projection::Options}}
dto! {pub struct CompiledPreparation{pub request:topo::Request,pub vector_strands:[usize;2],pub insert_strands:[usize;2],pub directional:bool,pub assumptions:Vec<String>}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("TOPO preparation declaration: {message}"))
}
fn text(value: &str, max: usize) -> bool {
    !value.trim().is_empty()
        && value.len() <= max
        && !value
            .chars()
            .any(|c| c.is_control() && !matches!(c, '\n' | '\r' | '\t'))
}
fn length(strand: &topo::PreparedStrand) -> Result<usize, Error> {
    if !text(&strand.id, 256) || strand.segments.is_empty() || strand.segments.len() > 64 {
        return Err(bad(
            "supply bounded explicit strand intervals and identities",
        ));
    }
    let mut length = 0usize;
    for interval in &strand.segments {
        if !text(&interval.source_role, 256) || interval.length == 0 {
            return Err(bad(
                "each vector segment names its actual source role and nonempty interval",
            ));
        }
        length = length
            .checked_add(interval.length)
            .filter(|value| *value <= 100_000)
            .ok_or(Error::TooLarge)?;
    }
    Ok(length)
}
/// Pure declaration conversion. Original models are not copied or inferred here;
/// call source_preflight/compile_preparation before native orientation inspection.
pub fn compile(request: &PreparationRequest) -> Result<CompiledPreparation, Error> {
    let vector = &request.vector;
    let insert = &request.insert;
    if !text(&vector.declaration_id, 256)
        || !text(&vector.assumption, 4096)
        || vector.strands.len() != 2
        || vector.pairing.len() != 1
        || !text(&insert.material_id, 256)
        || !text(&insert.template_role, 256)
    {
        return Err(bad(
            "declare two actual vector strands, one exact paired core, and named original-template PCR",
        ));
    }
    let top = length(&vector.strands[0])?;
    let bottom = length(&vector.strands[1])?;
    if vector.strands[0].id == vector.strands[1].id {
        return Err(bad("vector strand identities must be distinct"));
    }
    let mut activation_ids = BTreeSet::new();
    for strand in &vector.strands {
        if let physical::ThreePrime::Phosphotyrosyl { activation_id, .. } = &strand.three_prime {
            if !text(activation_id, 256) || !activation_ids.insert(activation_id) {
                return Err(bad(
                    "every declared activation has a bounded distinct identity",
                ));
            }
        }
    }
    let pair = &vector.pairing[0];
    let (paired, bottom_start) = match request.mode {
        PreparationMode::Ta => {
            if top < 2 || top != bottom || insert.tail != Some(TailingBase::A) {
                return Err(bad(
                    "TA requires equal vector lengths and explicit single A addition to each actual PCR strand",
                ));
            }
            (top - 1, 0)
        }
        PreparationMode::DirectionalCacc => {
            if bottom != top.checked_add(4).ok_or(Error::TooLarge)? || insert.tail.is_some() {
                return Err(bad(
                    "directional preparation requires the declared bottom-strand four-base5′ extension and blunt exact PCR with no terminal addition",
                ));
            }
            (top, 4)
        }
    };
    if pair.strand_a != 0
        || pair.start_a != 0
        || pair.strand_b != 1
        || pair.start_b != bottom_start
        || pair.length != paired
    {
        return Err(bad(
            "the explicit paired core does not match the selected TA or directional strand geometry",
        ));
    }
    for primer in [&insert.forward, &insert.reverse] {
        if !text(&primer.name, 256)
            || primer.annealing.is_empty()
            || primer
                .annealing
                .len()
                .checked_add(primer.tail.len())
                .is_none_or(|length| length > crate::primers::MAX_OLIGO)
            || primer.annealing.alphabet() != crate::Alphabet::Dna
        {
            return Err(bad("supply actual bounded named ordered DNA primers"));
        }
        // Validate literal tail notation without adding, trimming or replacing it.
        if !primer.tail.is_empty() {
            if crate::Sequence::dna(&primer.tail)?.as_str() != primer.tail {
                return Err(bad(
                    "primer tails must be explicit canonical DNA; no whitespace or sequence is silently rewritten",
                ));
            }
        }
    }
    let compiled=CompiledPreparation{
        request:topo::Request{recipe:topo::GenerationRecipe{schema_version:1,materials:vec![
            topo::MaterialRecipe::Prepared{declaration_id:vector.declaration_id.clone(),assumption:vector.assumption.clone(),strands:vector.strands.clone(),pairing:vector.pairing.clone()},
            topo::MaterialRecipe::ExactPcr{material_id:insert.material_id.clone(),template_role:insert.template_role.clone(),forward:insert.forward.clone(),reverse:insert.reverse.clone(),tail:insert.tail},
        ]},transfers:vec![],engine_limits:request.engine_limits.clone(),projection_options:request.projection_options.clone()},
        vector_strands:[0,1],insert_strands:[2,3],directional:request.mode==PreparationMode::DirectionalCacc,
        assumptions:vec![
            "Vector preparation, activation and every terminal state are explicit hypotheses. Original reference sequence or a commercial name does not establish activated stock; this builder does not simulate vector manufacture.".into(),
            "The existing native engine must reconstruct unique exact PCR, actual pairing, donor recognition and acceptor chemistry. No CACC, A tail or phosphorylation is silently added.".into(),
            "Directional displacement can retain an attached flap and withhold a paired contour. No repair, contour choice, circular material or kit efficiency is inferred from the intended target.".into(),
        ],
    };
    let mut count = Count(0);
    serde_json::to_writer(&mut count, &compiled).map_err(|_| Error::TooLarge)?;
    Ok(compiled)
}
struct Count(usize);
impl io::Write for Count {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|size| *size <= topo::MAX_REPORT_BYTES)
            .ok_or_else(|| io::Error::other("TOPO declaration exceeds wire budget"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
/// Borrowed original-source/operation guards, without simulating preparation.
pub fn source_preflight(inputs: &[topo::Input], request: &PreparationRequest) -> Result<(), Error> {
    let compiled = compile(request)?;
    topo::source_preflight(inputs, &compiled.request)
}
pub fn compile_preparation(
    inputs: &[topo::Input],
    request: &PreparationRequest,
) -> Result<CompiledPreparation, Error> {
    let compiled = compile(request)?;
    topo::source_preflight(inputs, &compiled.request)?;
    Ok(compiled)
}
