//! Source/edit-declared Golden Gate intent with selectively regenerated actual
//! PCR→catalog-digest preparations. Retained rows are never silently repaired.
use crate::{
    Error, Molecule, Sequence,
    analysis::Primer,
    assembly_intent,
    assembly_junction_edit::{self as junction, Edit, Intent},
    assembly_recipe::{self, Preparation},
    catalog_digest::{self, CatalogSelection},
    chemistry::{End, Polarity, Terminus},
    golden_gate::GoldenGateEnzyme,
    golden_gate_primer_design as original, ligation_recipe,
    verification::PlanInput,
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "golden-gate-junction-edit-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum Row {Retain(ligation_recipe::Fragment),DesignPcr{role:String,input_id:String,reverse:bool,forward_length:usize,reverse_length:usize}}}
dto! {pub struct Request {pub base:assembly_intent::Request,pub edits:Vec<Edit>,pub rows:Vec<Row>,pub enzyme:GoldenGateEnzyme,pub clamp:String,pub spacer:String,pub terminal_right_overhang:Option<String>,pub allow_nicked:bool}}
dto! {pub struct TailSegments {pub clamp:String,pub motif:String,pub spacer:String,pub fusion:String,pub inserted:String}}
impl TailSegments {
    pub fn sequence(&self) -> String {
        format!(
            "{}{}{}{}{}",
            self.clamp, self.motif, self.spacer, self.fusion, self.inserted
        )
    }
}
dto! {pub struct Design {pub role:String,pub input_id:String,pub input_revision:String,pub forward:Primer,pub reverse:Primer,pub forward_tail:TailSegments,pub reverse_tail:TailSegments,pub selection:CatalogSelection,pub fragment_index:usize,pub oriented_top:String,pub left:End,pub right:End}}
dto! {pub struct Report {pub context_revision:String,pub request:Request,pub intent:Intent,pub recipe:ligation_recipe::Recipe,pub designs:Vec<Option<Design>>,pub retained_rows:Vec<usize>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("Golden Gate junction edit: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn dna(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b"ACGT".contains(&b))
}
fn rc(s: &str) -> Result<String, Error> {
    if s.is_empty() {
        Ok(String::new())
    } else {
        Ok(Sequence::dna(s)?.reverse_complement()?.as_str().into())
    }
}
pub(crate) fn dimensions(r: &Request) -> Result<usize, Error> {
    if !r.base.material_bindings.is_empty() {
        return Err(bad(
            "primer proposals require reference ranges; material-strand primer templates are unsupported",
        ));
    }
    let p = r.enzyme.profile();
    let k = usize::try_from(p.bottom_cut - p.top_cut).map_err(|_| Error::TooLarge)?;
    let spacer = usize::try_from(p.top_cut)
        .map_err(|_| Error::TooLarge)?
        .checked_sub(p.recognition.len())
        .ok_or_else(|| bad("unsupported cut geometry"))?;
    if r.rows.len() != r.base.fragments.len()
        || r.rows.is_empty()
        || r.rows.len() > 32
        || r.clamp.len() > 12
        || !dna(&r.clamp)
        || r.spacer.len() != spacer
        || !dna(&r.spacer)
    {
        return Err(bad(
            "bounded concrete clamp/spacer and one row per base contribution are required",
        ));
    }
    let terminal_required =
        !r.base.circular && matches!(r.rows.last(), Some(Row::DesignPcr { .. }));
    match (&r.terminal_right_overhang, terminal_required) {
        (None, false) => {}
        (Some(t), true) if t.len() == k && dna(t) => {}
        _ => {
            return Err(bad(
                "only a selected final linear PCR row needs an explicit profile-sized terminal fusion; all other cases requireNone",
            ));
        }
    }
    for (row, f) in r.rows.iter().zip(&r.base.fragments) {
        let f = f.reference_fragment()?;
        let (role, input, reverse) = match row {
            Row::Retain(retained) => {
                retained_raw(retained)?;
                (
                    &retained.source.role,
                    &retained.source.input_id,
                    retained.source.reverse,
                )
            }
            Row::DesignPcr {
                role,
                input_id,
                reverse,
                forward_length,
                reverse_length,
            } => {
                if !(8..=80).contains(forward_length) || !(8..=80).contains(reverse_length) {
                    return Err(bad("selected annealing lengths must be8..80"));
                }
                (role, input_id, *reverse)
            }
        };
        if role != &f.role || input != &f.input_id || reverse != f.reverse {
            return Err(bad(
                "row role/order/input/orientation differs from the independent base intent",
            ));
        }
    }
    Ok(k)
}
fn retained_raw(row: &ligation_recipe::Fragment) -> Result<(), Error> {
    if row.source.role.len() > 128 || row.source.input_id.len() > 128 {
        return Err(Error::TooLarge);
    }
    let selection = match &row.source.preparation {
        Preparation::CatalogDigest { selection, .. } => selection,
        Preparation::PcrCatalogDigest {
            forward,
            reverse,
            selection,
            template_strand: crate::chemistry::Strand::Top,
            ..
        } => {
            for p in [forward, reverse] {
                if p.name.len() > 256 || p.annealing.len() > 1000 || p.tail.len() > 1000 {
                    return Err(Error::TooLarge);
                }
            }
            selection
        }
        _ => {
            return Err(bad(
                "retained Type IIS rows require a pinned CatalogDigest or PcrCatalogDigest preparation",
            ));
        }
    };
    selection.validate()?;
    Ok(())
}
pub(crate) fn retained_pin(
    row: &ligation_recipe::Fragment,
    catalog: &catalog_digest::catalog::Catalog,
    enzyme: GoldenGateEnzyme,
) -> Result<(), Error> {
    let selection = match &row.source.preparation {
        Preparation::CatalogDigest { selection, .. }
        | Preparation::PcrCatalogDigest {
            selection,
            template_strand: crate::chemistry::Strand::Top,
            ..
        } => selection,
        _ => return Err(bad("retained digest is not catalog pinned")),
    };
    let enzymes = selection.resolved_enzymes(catalog)?;
    if enzymes.len() != 1
        || selection.enzymes[0].id != enzyme.profile().name.to_ascii_lowercase()
        || assembly_recipe::hash(&enzymes[0])? != assembly_recipe::hash(&enzyme.profile())?
    {
        return Err(bad(
            "retained enzyme pin differs from the chosen Type IIS profile",
        ));
    }
    Ok(())
}
fn lengths(inputs: &[(&str, &Molecule)], r: &Request) -> Result<Vec<usize>, Error> {
    r.base
        .fragments
        .iter()
        .map(|f| {
            inputs
                .iter()
                .find(|(id, _)| *id == f.input_id)
                .map(|(_, m)| m.sequence.len())
                .ok_or_else(|| bad("missing original source"))
        })
        .collect()
}
fn right_fusion(
    inputs: &[(&str, &Molecule)],
    r: &Request,
    adjusted: &[crate::overlap_primer_design::Fragment],
    i: usize,
    k: usize,
) -> Result<String, Error> {
    if i + 1 < r.rows.len() || r.base.circular {
        let next = &adjusted[(i + 1) % r.rows.len()];
        if next.end - next.start < k {
            return Err(bad(
                "downstream adjusted core is shorter than the required fusion",
            ));
        }
        let source = inputs
            .iter()
            .find(|(id, _)| *id == next.input_id)
            .map(|(_, m)| *m)
            .ok_or_else(|| bad("missing downstream original"))?;
        let fusion = junction::prefix(source, next, k)?;
        if !dna(&fusion) {
            return Err(bad(
                "selected outgoing fusion requires concrete downstream source calls",
            ));
        }
        Ok(fusion)
    } else {
        r.terminal_right_overhang
            .clone()
            .ok_or_else(|| bad("missing explicit terminal fusion"))
    }
}
/// Borrowed model/catalog/work guards before source copies and scientific scans.
pub fn preflight(inputs: &[(&str, &Molecule)], r: &Request) -> Result<(), Error> {
    let k = dimensions(r)?;
    junction::intent_preflight(inputs, &r.base, &r.edits, None)?;
    let adjusted = junction::adjustment(&r.base, &r.edits, &lengths(inputs, r)?)?
        .iter()
        .map(assembly_intent::Contribution::reference_fragment)
        .collect::<Result<Vec<_>, _>>()?;
    let inserted = junction::inserts(&r.base, &r.edits)?;
    let catalog = catalog_digest::builtin_catalog()?;
    let mut retained = Vec::new();
    for row in &r.rows {
        if let Row::Retain(row) = row {
            retained_pin(row, &catalog, r.enzyme)?;
            retained.push(row.source.clone());
        }
    }
    if !retained.is_empty() {
        assembly_recipe::preparation_rows_preflight(inputs, &r.base.name, &retained)?;
    }
    let p = r.enzyme.profile();
    let tail = r.clamp.len() + p.recognition.len() + r.spacer.len();
    let mut work = 0usize;
    let mut output = 0usize;
    for (i, (row, f)) in r.rows.iter().zip(&adjusted).enumerate() {
        let Row::DesignPcr {
            forward_length,
            reverse_length,
            ..
        } = row
        else {
            continue;
        };
        let source = inputs
            .iter()
            .find(|(id, _)| *id == f.input_id)
            .map(|(_, m)| *m)
            .ok_or_else(|| bad("missing selected original"))?;
        original::physical_source(source)?;
        let core = f.end - f.start;
        let annealed=core.checked_sub(if f.reverse{k}else{0}).filter(|n|*n>=forward_length+reverse_length&&core>k).ok_or_else(||bad("selected adjusted core must fit both nonoverlapping domains and its oriented fusion"))?;
        let _ = annealed;
        if !inserted[i].bytes().all(|b| b"ACGT".contains(&b)) {
            return Err(bad("selected literal insert requires concrete DNA"));
        }
        let forward_tail = tail + if f.reverse { k + inserted[i].len() } else { 0 };
        let reverse_tail = tail + if f.reverse { k } else { k + inserted[i].len() };
        if forward_length + forward_tail > 200 || reverse_length + reverse_tail > 200 {
            return Err(bad(
                "complete selected oligo including fusion and insert exceeds200 bases",
            ));
        }
        right_fusion(inputs, r, &adjusted, i, k)?;
        work = work
            .checked_add(
                source
                    .sequence
                    .len()
                    .checked_mul((forward_length + reverse_length + 16) * 8)
                    .ok_or(Error::TooLarge)?,
            )
            .and_then(|n| {
                core.checked_add(tail * 2 + k + inserted[i].len())
                    .and_then(|amp| amp.checked_mul(p.recognition.len() * 4))
                    .and_then(|cuts| n.checked_add(cuts))
            })
            .filter(|n| *n <= 100_000_000)
            .ok_or(Error::TooLarge)?;
        output = output
            .checked_add(core + inserted[i].len() + 128)
            .filter(|n| *n <= 2 << 20)
            .ok_or(Error::TooLarge)?;
    }
    // Whole-original uniqueness is the existing Type IIS preparation profile,
    // and applies only to explicitly selected rows.
    for (row, f) in r.rows.iter().zip(&adjusted) {
        if matches!(row, Row::DesignPcr { .. }) {
            let source = inputs
                .iter()
                .find(|(id, _)| *id == f.input_id)
                .map(|(_, m)| *m)
                .ok_or_else(|| bad("missing original"))?;
            if !dna(source.sequence.as_str()) {
                return Err(bad(
                    "selected original contains ambiguous calls; unique-PCR/cleavage preparation is withheld",
                ));
            }
        }
    }
    assembly_recipe::hash(&(inputs, r))?;
    Ok(())
}
pub fn source_preflight(inputs: &[(&str, &Molecule)], r: &Request) -> Result<(), Error> {
    preflight(inputs, r)
}
fn pins(inputs: &[PlanInput]) -> Result<Vec<(&str, &str)>, Error> {
    if inputs.is_empty() || inputs.len() > 32 {
        return Err(Error::TooLarge);
    }
    for i in inputs {
        if !hex(&i.revision) || i.revision != i.molecule.revision() {
            return Err(bad("original input revision changed"));
        }
    }
    Ok(inputs
        .iter()
        .map(|i| (i.id.as_str(), i.revision.as_str()))
        .collect())
}
pub fn design_revision(inputs: &[PlanInput], r: &Request) -> Result<String, Error> {
    if inputs.len() > 32 {
        return Err(Error::TooLarge);
    }
    preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
    )?;
    assembly_recipe::hash(&(ENGINE_VERSION, pins(inputs)?, r))
}
fn segments(r: &Request, fusion: String, inserted: String) -> TailSegments {
    TailSegments {
        clamp: r.clamp.clone(),
        motif: r.enzyme.profile().recognition,
        spacer: r.spacer.clone(),
        fusion,
        inserted,
    }
}
pub fn design(inputs: &[PlanInput], r: &Request) -> Result<Report, Error> {
    let context_revision = design_revision(inputs, r)?;
    let intent = junction::build_intent(inputs, &r.base, &r.edits, None)?;
    let adjusted = intent
        .adjusted_fragments
        .iter()
        .map(assembly_intent::Contribution::reference_fragment)
        .collect::<Result<Vec<_>, _>>()?;
    let insert = junction::inserts(&r.base, &r.edits)?;
    let k = dimensions(r)?;
    let catalog = catalog_digest::builtin_catalog()?;
    let selection =
        catalog_digest::select(&catalog, &[r.enzyme.profile().name.to_ascii_lowercase()])?;
    let enzymes = selection.resolved_enzymes(&catalog)?;
    if enzymes.len() != 1
        || assembly_recipe::hash(&enzymes[0])? != assembly_recipe::hash(&r.enzyme.profile())?
    {
        return Err(bad("builtin selected enzyme geometry changed"));
    }
    let borrowed = inputs
        .iter()
        .map(|i| (i.id.as_str(), &i.molecule))
        .collect::<Vec<_>>();
    let mut specs = Vec::new();
    for (i, (row, f)) in r.rows.iter().zip(&adjusted).enumerate() {
        if let Row::DesignPcr {
            forward_length,
            reverse_length,
            ..
        } = row
        {
            specs.push(original::Fragment {
                role: f.role.clone(),
                input_id: f.input_id.clone(),
                start: f.start,
                end: f.end,
                reverse: f.reverse,
                forward_length: *forward_length,
                reverse_length: *reverse_length,
                right_overhang: right_fusion(&borrowed, r, &adjusted, i, k)?,
            });
        }
    }
    let template = original::Request {
        name: r.base.name.clone(),
        fragments: specs,
        enzyme: r.enzyme,
        clamp: r.clamp.clone(),
        spacer: r.spacer.clone(),
    };
    let mut recipe_rows = Vec::with_capacity(r.rows.len());
    let mut designs = Vec::with_capacity(r.rows.len());
    let mut retained_rows = Vec::new();
    for (i, (row, f)) in r.rows.iter().zip(&adjusted).enumerate() {
        match row {
            Row::Retain(row) => {
                recipe_rows.push(row.clone());
                designs.push(None);
                retained_rows.push(i);
            }
            Row::DesignPcr { .. } => {
                let original = inputs
                    .iter()
                    .find(|i| i.id == f.input_id)
                    .ok_or_else(|| bad("missing selected original"))?;
                let spec = template
                    .fragments
                    .iter()
                    .find(|s| s.role == f.role)
                    .ok_or_else(|| bad("missing selected specification"))?;
                let (preparation, evidence) = original::design_row(
                    original, &template, spec, &insert[i], &selection, &enzymes,
                )?;
                let forward_tail = segments(
                    r,
                    if f.reverse {
                        rc(&spec.right_overhang)?
                    } else {
                        String::new()
                    },
                    if f.reverse {
                        rc(&insert[i])?
                    } else {
                        String::new()
                    },
                );
                let reverse_tail = segments(
                    r,
                    if f.reverse {
                        evidence.oriented_top[..k].into()
                    } else {
                        rc(&spec.right_overhang)?
                    },
                    if f.reverse {
                        String::new()
                    } else {
                        rc(&insert[i])?
                    },
                );
                if evidence.forward.tail != forward_tail.sequence()
                    || evidence.reverse.tail != reverse_tail.sequence()
                {
                    return Err(bad(
                        "shared physical preparation differs from explicit tail segments",
                    ));
                }
                recipe_rows.push(ligation_recipe::Fragment {
                    source: preparation,
                    treatment: None,
                });
                designs.push(Some(Design {
                    role: evidence.role,
                    input_id: evidence.input_id,
                    input_revision: evidence.input_revision,
                    forward: evidence.forward,
                    reverse: evidence.reverse,
                    forward_tail,
                    reverse_tail,
                    selection: evidence.selection,
                    fragment_index: evidence.fragment_index,
                    oriented_top: evidence.oriented_top,
                    left: evidence.left,
                    right: evidence.right,
                }));
            }
        }
    }
    let recipe = ligation_recipe::Recipe {
        name: r.base.name.clone(),
        fragments: recipe_rows,
        circular: r.base.circular,
        allow_nicked: r.allow_nicked,
    };
    ligation_recipe::source_preflight(&borrowed, &intent.target, &recipe)?;
    let result=Report{context_revision:context_revision.clone(),request:r.clone(),intent,recipe,designs,retained_rows,warnings:vec!["Independent source trim/insert intent plus explicit selected Type IIS preparations only. Retained rows and their chemistry/treatments remain literal and can be incompatible; fresh full ordered ligation verification is required.".into(),"Selected primers have unique exact original-source geometry and exactly two eligible added sites under the existing bounded model. This establishes neither reaction yield, methylation clearance, empirical activity nor exhaustive mixture uniqueness. No internal-site domestication or automatic fusion optimization is performed.".into()]};
    result.validate_binding(&context_revision)?;
    Ok(result)
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let k = dimensions(&self.request)?;
        if !hex(expected)
            || self.context_revision != expected
            || self.designs.len() != self.request.rows.len()
            || self.recipe.fragments.len() != self.request.rows.len()
        {
            return Err(bad("report context or row count changed"));
        }
        self.intent
            .validate_binding(&self.intent.context_revision)?;
        if assembly_recipe::hash(&self.request.base)?
            != assembly_recipe::hash(&self.intent.base_request)?
            || self.request.edits != self.intent.edits
            || self.recipe.name != self.request.base.name
            || self.recipe.circular != self.request.base.circular
            || self.recipe.allow_nicked != self.request.allow_nicked
        {
            return Err(bad("independent intent or ligation options changed"));
        }
        let raw = self
            .intent
            .target
            .metadata
            .get("assembly_junction_edit:inputs")
            .filter(|s| s.len() < 16_384)
            .ok_or_else(|| bad("missing source revision pins"))?;
        let pins: Vec<(String, String)> =
            serde_json::from_str(raw).map_err(|_| bad("invalid source pins"))?;
        if assembly_recipe::hash(&(ENGINE_VERSION, &pins, &self.request))? != self.context_revision
        {
            return Err(bad("source/request context changed"));
        }
        let inserted = junction::inserts(&self.request.base, &self.request.edits)?;
        let catalog = catalog_digest::builtin_catalog()?;
        let expected_selection = catalog_digest::select(
            &catalog,
            &[self.request.enzyme.profile().name.to_ascii_lowercase()],
        )?;
        let mut retained = Vec::new();
        for (i, ((request, row), design)) in self
            .request
            .rows
            .iter()
            .zip(&self.recipe.fragments)
            .zip(&self.designs)
            .enumerate()
        {
            match request {
                Row::Retain(original) => {
                    retained_pin(original, &catalog, self.request.enzyme)?;
                    if design.is_some()
                        || assembly_recipe::hash(original)? != assembly_recipe::hash(row)?
                    {
                        return Err(bad("retained preparation, treatment or repair changed"));
                    }
                    retained.push(i);
                }
                Row::DesignPcr {
                    forward_length,
                    reverse_length,
                    ..
                } => {
                    let f = &self.intent.fragments[i];
                    let geometry = self.intent.adjusted_fragments[i].reference_fragment()?;
                    let evidence = design
                        .as_ref()
                        .ok_or_else(|| bad("selected design missing"))?;
                    let fusion = if i + 1 < self.designs.len() || self.request.base.circular {
                        let ni = (i + 1) % self.designs.len();
                        let next = &self.intent.fragments[ni];
                        let end = next
                            .output_start
                            .checked_add(k)
                            .filter(|end| *end <= next.output_end - inserted[ni].len())
                            .ok_or_else(|| bad("next original core is shorter than fusion"))?;
                        self.intent.target.sequence.as_str()[next.output_start..end].to_owned()
                    } else {
                        self.request
                            .terminal_right_overhang
                            .clone()
                            .ok_or_else(|| bad("missing terminal fusion"))?
                    };
                    let core_length = geometry.end - geometry.start;
                    if core_length <= k
                        || core_length - (if geometry.reverse { k } else { 0 })
                            < *forward_length + *reverse_length
                    {
                        return Err(bad(
                            "selected returned core cannot supply its declared domains/fusion",
                        ));
                    }
                    let wanted =
                        &self.intent.target.sequence.as_str()[f.output_start..f.output_end];
                    if !dna(wanted)
                        || evidence.oriented_top != wanted
                        || evidence.role != f.role
                        || evidence.input_id != f.input_id
                        || evidence.input_revision != f.input_revision
                        || evidence.fragment_index != 1
                        || assembly_recipe::hash(&evidence.selection)?
                            != assembly_recipe::hash(&expected_selection)?
                        || row.source.role != f.role
                        || row.source.input_id != f.input_id
                        || row.source.reverse != geometry.reverse
                        || row.treatment.is_some()
                    {
                        return Err(bad(
                            "selected source, actual retained top or enzyme pin changed",
                        ));
                    }
                    let ft = segments(
                        &self.request,
                        if geometry.reverse {
                            rc(&fusion)?
                        } else {
                            String::new()
                        },
                        if geometry.reverse {
                            rc(&inserted[i])?
                        } else {
                            String::new()
                        },
                    );
                    let rt = segments(
                        &self.request,
                        if geometry.reverse {
                            wanted[..k].into()
                        } else {
                            rc(&fusion)?
                        },
                        if geometry.reverse {
                            String::new()
                        } else {
                            rc(&inserted[i])?
                        },
                    );
                    if assembly_recipe::hash(&ft)? != assembly_recipe::hash(&evidence.forward_tail)?
                        || assembly_recipe::hash(&rt)?
                            != assembly_recipe::hash(&evidence.reverse_tail)?
                    {
                        return Err(bad("fusion/insert segmentation changed"));
                    }
                    for (p, t, length, name) in [
                        (&evidence.forward, &ft, *forward_length, "forward"),
                        (&evidence.reverse, &rt, *reverse_length, "reverse"),
                    ] {
                        if p.name != format!("{}-{name}", f.role)
                            || p.annealing.len() != length
                            || !dna(p.annealing.as_str())
                            || p.annealing.alphabet() != crate::Alphabet::Dna
                            || p.annealing.len() + p.tail.len() > 200
                            || p.tail != t.sequence()
                        {
                            return Err(bad("selected complete primer differs from its segments"));
                        }
                    }
                    if evidence.left.polarity != Polarity::FivePrime
                        || evidence.right.polarity != Polarity::FivePrime
                        || evidence.left.overhang != wanted[..k]
                        || evidence.right.overhang != rc(&fusion)?
                        || evidence.left.chemistry != Terminus::CLEAVED
                        || evidence.right.chemistry != Terminus::CLEAVED
                    {
                        return Err(bad("selected actual oriented end geometry changed"));
                    }
                    let Preparation::PcrCatalogDigest {
                        forward,
                        reverse,
                        selection,
                        fragment_index,
                        repair: None,
                        template_strand: crate::chemistry::Strand::Top,
                        ..
                    } = &row.source.preparation
                    else {
                        return Err(bad(
                            "selected preparation must remain PCR→pinned digest without implicit repair",
                        ));
                    };
                    if assembly_recipe::hash(forward)? != assembly_recipe::hash(&evidence.forward)?
                        || assembly_recipe::hash(reverse)?
                            != assembly_recipe::hash(&evidence.reverse)?
                        || assembly_recipe::hash(selection)?
                            != assembly_recipe::hash(&evidence.selection)?
                        || *fragment_index != evidence.fragment_index
                    {
                        return Err(bad(
                            "selected preparation does not match fresh physical evidence",
                        ));
                    }
                }
            }
        }
        if self.retained_rows != retained
            || self.warnings.len() > 1024
            || self.warnings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad("invalid retained-index or warning receipt"));
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
