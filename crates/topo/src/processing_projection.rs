//! Fresh paired-contour proof after explicit ordered physical graph actions. No imported
//! graph or positive report is accepted as an input to projection.
use crate::{Budget, Call, Limits, Substrates, Transfer, processing as actions, projection as p};
use serde::{Deserialize, Serialize};
pub const ENGINE: &str = "topo-actions-paired-contour-proof-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Selection{pub processing_context_revision:String,pub component:usize,pub reference_seed:Call}}
dto! {pub struct ProjectionReport{pub engine:String,pub context_revision:String,pub selection:Selection,pub options:p::Options,pub processing:actions::Report,pub outcome:p::Outcome}}
pub type Error = p::Error;
fn bad(s: &'static str) -> Error {
    Error::InvalidSelection(s)
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn context(select: &Selection, o: &p::Options) -> Result<String, Error> {
    Ok(crate::hash(&(ENGINE, select, o))?)
}
/// Borrowed bounds and combined live-workspace reservation before simulation.
/// Projection's own narrower budgets still yield an explicit withheld outcome.
pub fn preflight(
    s: &Substrates,
    t: &[Transfer],
    l: &Limits,
    r: &actions::Request,
    select: &Selection,
    o: &p::Options,
) -> Result<usize, Error> {
    p::limits(&o.limits)?;
    let n = actions::preflight(s, t, l, r)?;
    if !hex(&select.processing_context_revision)
        || select.component >= r.limits.max_product_chains
        || s.strands
            .get(select.reference_seed.strand)
            .is_none_or(|v| select.reference_seed.index >= v.sequence.len())
    {
        return Err(bad("bounded source/component/seed identity required"));
    }
    if select.processing_context_revision != actions::revision(s, t, l, r)? {
        return Err(bad("source, transfer or processing context changed"));
    }
    let steps = crate::add(t.len(), r.actions.len())?;
    let bytes = crate::add(
        crate::mul(n, 768)?,
        crate::add(
            crate::mul((n / 2).min(o.limits.max_mapping_segments), 128)?,
            crate::add(
                crate::mul(s.strands.len(), 8192)?,
                crate::add(
                    crate::mul(s.pairing.len(), 256)?,
                    crate::add(crate::mul(steps, 16384)?, 65536)?,
                )?,
            )?,
        )?,
    )?;
    crate::bounded(bytes, l.max_work_bytes, Budget::WorkBytes)?;
    crate::bounded(
        crate::add(crate::mul(n, 384)?, crate::mul(steps, 256)?)?,
        l.max_work_cells,
        Budget::WorkCells,
    )?;
    // Parent and treated products each retain n calls; the paired view adds
    // at most n further displayed calls across both reference strands.
    crate::bounded(crate::mul(n, 3)?, l.max_output_bases, Budget::OutputBases)?;
    Ok(n)
}
pub fn revision(
    s: &Substrates,
    t: &[Transfer],
    l: &Limits,
    r: &actions::Request,
    select: &Selection,
    o: &p::Options,
) -> Result<String, Error> {
    preflight(s, t, l, r, select, o)?;
    context(select, o)
}
/// Runs original transfers and processing exactly once, then borrows the freshly
/// produced child graph for the private contour algorithm. Full ancestry moves
/// intact into the result, including withheld views and other components.
pub fn project(
    s: &Substrates,
    t: &[Transfer],
    l: &Limits,
    r: &actions::Request,
    select: &Selection,
    o: &p::Options,
) -> Result<ProjectionReport, Error> {
    let context_revision = revision(s, t, l, r, select, o)?;
    let processing = actions::simulate(s, t, l, r)?;
    if select.component >= processing.product.components.len() {
        return Err(bad("component absent from freshly processed graph"));
    }
    // Empty operations are vacuously complete for policy purposes; their native
    // report flags intentionally do not assert a successful nonempty reaction.
    let all = processing
        .parent
        .steps
        .iter()
        .all(|x| x.outcome == crate::Outcome::Applied)
        && processing
            .steps
            .iter()
            .all(|x| x.outcome() == crate::Outcome::Applied);
    let outcome = match p::prove(
        s,
        &processing.product,
        processing.source_bases,
        all,
        select.component,
        select.reference_seed,
        o,
    ) {
        Ok(v) => p::Outcome::Proven(v),
        Err(w) => p::Outcome::Withheld(w),
    };
    let report = ProjectionReport {
        engine: ENGINE.into(),
        context_revision,
        selection: select.clone(),
        options: o.clone(),
        processing,
        outcome,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl ProjectionReport {
    /// Readonly envelope/geometry-shape binding, never a chemical reconstruction.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        p::limits(&self.options.limits)?;
        if self.engine != ENGINE
            || !hex(expected)
            || self.context_revision != expected
            || !hex(&self.selection.processing_context_revision)
            || self.context_revision != context(&self.selection, &self.options)?
        {
            return Err(bad("processed projection engine/context changed"));
        }
        self.processing
            .validate_binding(&self.selection.processing_context_revision)?;
        if self.selection.component >= self.processing.product.components.len()
            || self.selection.reference_seed.strand >= self.processing.parent.limits.max_strands
            || self.selection.reference_seed.index >= self.processing.source_bases
        {
            return Err(bad("processed selection shape changed"));
        }
        let all = self
            .processing
            .parent
            .steps
            .iter()
            .all(|x| x.outcome == crate::Outcome::Applied)
            && self
                .processing
                .steps
                .iter()
                .all(|x| x.outcome() == crate::Outcome::Applied);
        if let p::Outcome::Proven(v) = &self.outcome {
            let n = v.top.len();
            let cap = self.options.limits.max_pairs;
            if n == 0
                || n > cap
                || n.checked_mul(2)
                    .is_none_or(|x| x > self.processing.source_bases)
                || v.bottom.len() != n
                || v.seed_reference_position >= n
                || v.source_maps.len() > self.options.limits.max_mapping_segments
                || v.nicks.len() > self.processing.product.nicks.len()
                || v.all_requested_transfers_applied != all
                || (self.options.require_all_requested_transfers_applied && !all)
                || (self.options.require_known_terminal_chemistry
                    && !v.fully_known_terminal_chemistry)
            {
                return Err(bad("processed contour shape or policy changed"));
            }
            let mut cursor = 0usize;
            for map in &v.source_maps {
                if map.reference_start != cursor
                    || map.length == 0
                    || map.top_strand >= self.processing.parent.limits.max_strands
                    || map.bottom_strand >= self.processing.parent.limits.max_strands
                    || map
                        .top_start
                        .checked_add(map.length)
                        .is_none_or(|x| x > self.processing.source_bases)
                    || map
                        .bottom_start
                        .checked_add(map.length)
                        .is_none_or(|x| x > self.processing.source_bases)
                {
                    return Err(bad("processed paired-source map changed"));
                }
                cursor = cursor
                    .checked_add(map.length)
                    .ok_or_else(|| bad("paired-source map overflow"))?;
                if cursor > n {
                    return Err(bad("paired-source map exceeds contour"));
                }
            }
            if cursor != n {
                return Err(bad("paired-source map incomplete"));
            }
            for nick in &v.nicks {
                if nick.boundary > n
                    || !self.processing.product.nicks.iter().any(|x| {
                        x.three == nick.three
                            && x.five == nick.five
                            && x.three_prime == nick.three_prime
                            && x.five_prime == nick.five_prime
                    })
                {
                    return Err(bad("projected nick absent from processed graph"));
                }
            }
        }
        crate::report_size(self, self.options.limits.max_report_bytes).map_err(|_| {
            Error::ReportBytes {
                limit: self.options.limits.max_report_bytes,
            }
        })?;
        Ok(())
    }
}
pub fn validate_report(
    s: &Substrates,
    t: &[Transfer],
    l: &Limits,
    r: &actions::Request,
    select: &Selection,
    o: &p::Options,
    report: &ProjectionReport,
) -> Result<(), Error> {
    report.validate_binding(&revision(s, t, l, r, select, o)?)?;
    if &report.selection != select || &report.options != o {
        return Err(bad("processed selection/options changed"));
    }
    actions::validate_report(s, t, l, r, &report.processing)?;
    if let p::Outcome::Proven(v) = &report.outcome {
        for map in &v.source_maps {
            let top = s
                .strands
                .get(map.top_strand)
                .and_then(|x| {
                    x.sequence
                        .as_str()
                        .as_bytes()
                        .get(map.top_start..map.top_start.checked_add(map.length)?)
                })
                .ok_or_else(|| bad("top source map out of bounds"))?;
            let bottom = s
                .strands
                .get(map.bottom_strand)
                .and_then(|x| {
                    x.sequence
                        .as_str()
                        .as_bytes()
                        .get(map.bottom_start..map.bottom_start.checked_add(map.length)?)
                })
                .ok_or_else(|| bad("bottom source map out of bounds"))?;
            if v.top
                .as_str()
                .as_bytes()
                .get(map.reference_start..map.reference_start + map.length)
                != Some(top)
            {
                return Err(bad("projected top differs from immutable calls"));
            }
            for (k, b) in bottom.iter().rev().enumerate() {
                if crate::comp(*b) != Some(top[k]) {
                    return Err(bad("projected paired-source calls differ"));
                }
            }
            // bottom is stored actual 5′→3′, antiparallel to the complete top.
            let start = v.bottom.len() - map.reference_start - map.length;
            if v.bottom.as_str().as_bytes().get(start..start + map.length) != Some(bottom) {
                return Err(bad("projected bottom differs from immutable calls"));
            }
        }
    }
    Ok(())
}
pub fn replay(
    s: &Substrates,
    t: &[Transfer],
    l: &Limits,
    r: &actions::Request,
    select: &Selection,
    o: &p::Options,
    report: &ProjectionReport,
) -> Result<(), Error> {
    validate_report(s, t, l, r, select, o, report)?;
    if &project(s, t, l, r, select, o)? != report {
        return Err(bad(
            "processed projection differs from fresh source-based reconstruction",
        ));
    }
    Ok(())
}
