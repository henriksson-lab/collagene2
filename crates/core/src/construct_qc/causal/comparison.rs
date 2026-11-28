//! One bounded assessment over cached native history and declared CDS pairs.
use super::*;
use super::{coding::Coding, history::History};
use crate::{
    construct_qc::{Citation, Verdict},
    material::{lineage, run},
    verification::ancestry,
};
use collagene_model::{Location, annotations::biological_path::BiologicalCall};
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Clone, Copy)]
pub(super) struct ParentCall {
    pub ordinal: usize,
    pub call: BiologicalCall,
}
#[derive(Clone, Copy)]
pub(super) enum CellKind {
    Parent(ParentCall),
    Introduced,
    Foreign,
    Outside,
    Unknown,
}
#[derive(Clone, Copy)]
pub(super) struct Cell {
    pub kind: CellKind,
    pub changed: bool,
    pub causes: usize,
}
pub(super) struct CausePool {
    sets: Vec<Vec<TransitionPin>>,
    index: BTreeMap<String, usize>,
    last: usize,
    bytes: usize,
}
impl CausePool {
    fn new() -> Self {
        Self {
            sets: vec![vec![]],
            index: BTreeMap::new(),
            last: 0,
            bytes: 0,
        }
    }
    fn intern(&mut self, causes: Vec<TransitionPin>, work: &mut Work) -> Result<usize, Error> {
        if causes.len() > 64 || !work.charge(causes.len().max(1)) {
            return Err(Error::TooLarge);
        }
        if self.sets[self.last] == causes {
            return Ok(self.last);
        }
        crate::encoding::bounded_json::preflight(&causes, 64 << 10)?;
        let key = crate::encoding::bounded_json::hash(&causes)?;
        if let Some(&index) = self.index.get(&key) {
            if self.sets[index] != causes {
                return Err(Error::Unsupported(
                    "Native cause context identity collision".into(),
                ));
            }
            self.last = index;
            return Ok(index);
        }
        if self.sets.len() >= 4096 {
            return Err(Error::TooLarge);
        }
        let bytes = crate::verification::evidence::evidence_bytes(&causes, (1 << 20) - self.bytes)?;
        self.bytes += bytes;
        let index = self.sets.len();
        self.sets.push(causes);
        self.index.insert(key, index);
        self.last = index;
        Ok(index)
    }
    pub fn get(&self, index: usize) -> &[TransitionPin] {
        &self.sets[index]
    }
    pub fn combine(&self, indices: impl IntoIterator<Item = usize>) -> Vec<TransitionPin> {
        let mut result = Vec::new();
        for index in indices {
            for cause in self.get(index) {
                if !result.contains(cause) {
                    if result.len() >= 64 {
                        return vec![];
                    }
                    result.push(cause.clone());
                }
            }
        }
        result
    }
}
pub(super) struct OutputBudget {
    remaining: usize,
    pub exhausted: bool,
    pub segments: usize,
    pub stops: usize,
}
impl OutputBudget {
    fn new() -> Self {
        Self {
            remaining: MAX_REPORT_BYTES / 2,
            exhausted: false,
            segments: 0,
            stops: 0,
        }
    }
    pub fn retain(&mut self, value: &impl Serialize) -> bool {
        match crate::verification::evidence::evidence_bytes(value, self.remaining) {
            Ok(bytes) => {
                self.remaining -= bytes;
                true
            }
            Err(_) => {
                self.exhausted = true;
                false
            }
        }
    }
}
pub(super) fn finding(
    rule: Rule,
    verdict: Verdict,
    reason: impl Into<String>,
    coordinates: Vec<Location>,
    assumed: bool,
) -> Finding {
    Finding {
        rule,
        verdict,
        reason: reason.into(),
        product_coordinates: coordinates,
        citations: vec![Citation {
            source: "Collagene native causal guards".into(),
            locator: "construct_qc/causal; verification/ancestry; canonical translation/code pins"
                .into(),
        }],
        under_declared_annotation_assumptions: assumed,
    }
}
fn diagnostic(error: impl std::fmt::Display) -> String {
    let message = error.to_string();
    if message.len() <= 3000 {
        message
    } else {
        "Native owner refused this comparison; overlong diagnostic detail was not retained.".into()
    }
}
fn unresolved(row: &mut PairReport, reason: &str) {
    for rule in [
        Rule::CodingCorrespondence,
        Rule::JunctionRegister,
        Rule::NewInternalStop,
    ] {
        row.findings
            .push(finding(rule, Verdict::Indeterminate, reason, vec![], false));
    }
}
pub(super) fn assess(workspace: &Workspace, request: &Request) -> Result<Report, Error> {
    request.validate_shape()?;
    let fingerprint = lineage::workspace_fingerprint(workspace)?;
    let trace_budget = ancestry::TraceBudget::new(ancestry::Limits {
        max_lookup_work: request.limits.max_work,
        ..Default::default()
    })?;
    let history = History::prepare(workspace, request, trace_budget)?;
    let mut session = history.session();
    let mut work = Work::new(request.limits.max_work);
    let mut coding_remaining = request.limits.max_coding_calls;
    let mut output = OutputBudget::new();
    let mut rows = Vec::with_capacity(request.pairs.len());
    for pair in &request.pairs {
        let mut row = PairReport {
            id: pair.id.clone(), source: pair.source.clone(), product: pair.product.clone(),
            source_cds: pair.source_cds.clone(), product_cds: pair.product_cds.clone(),
            code: crate::construct_qc::code_pin(pair.genetic_code)?,
            annotation_correspondence: AnnotationCorrespondence { verdict: Verdict::Indeterminate, reason: "Native annotation identity has not been established; pairing is explicitly declared.".into() },
            mapping_complete: false, findings: vec![], register_segments: vec![], new_stops: vec![],
        };
        let history_evidence = session.history(&pair.product);
        row.findings.push(finding(
            Rule::NativeHistory,
            if history_evidence.available {
                Verdict::Pass
            } else {
                Verdict::Indeterminate
            },
            history_evidence.reason,
            vec![],
            false,
        ));
        let data = run::resolve(workspace, &pair.source).and_then(|source| {
            run::resolve(workspace, &pair.product).map(|product| (source, product))
        });
        let (source, product) = match data {
            Ok(data) => data,
            Err(error) => {
                unresolved(&mut row, &diagnostic(error));
                rows.push(row);
                continue;
            }
        };
        let before = Coding::load(
            source,
            &pair.source.revision,
            &pair.source_cds,
            pair.genetic_code,
            &mut coding_remaining,
        );
        let after = Coding::load(
            product,
            &pair.product.revision,
            &pair.product_cds,
            pair.genetic_code,
            &mut coding_remaining,
        );
        let (before, after) = match (before, after) {
            (Ok(before), Ok(after)) => (before, after),
            (before, after) => {
                for error in [before.err(), after.err()].into_iter().flatten() {
                    row.findings.push(finding(
                        Rule::AnnotationBasis,
                        Verdict::Indeterminate,
                        diagnostic(error),
                        vec![],
                        false,
                    ));
                }
                unresolved(
                    &mut row,
                    "Current explicit CDS/code assumptions are unavailable, stale, unreviewed or outside the bounded canonical model.",
                );
                rows.push(row);
                continue;
            }
        };
        row.findings.push(finding(Rule::AnnotationBasis, Verdict::Informational, "Both current CDS roles/coordinates are explicitly declared and not flagged for review. This is not verification of biological identity.", vec![], true));
        if !history_evidence.available || output.exhausted || work.remaining() == 0 {
            unresolved(
                &mut row,
                "Causal comparison requires available native history and remaining shared work/evidence budget.",
            );
            rows.push(row);
            continue;
        }
        let annotation = session.map_feature(
            &pair.product,
            &pair.source,
            pair.product_cds.feature_index,
            &mut work,
        );
        let native_annotation = annotation.complete
            && !annotation.parent_feature_indices.is_empty()
            && annotation
                .parent_feature_indices
                .iter()
                .all(|index| *index == pair.source_cds.feature_index);
        row.annotation_correspondence = AnnotationCorrespondence {
            verdict: if native_annotation {
                Verdict::Pass
            } else {
                Verdict::Indeterminate
            },
            reason: if native_annotation {
                "Actual native feature ordinals connect the selected annotations; biological roles remain declared assumptions.".into()
            } else {
                annotation.reason.unwrap_or_else(|| "The declared CDS pair is not established as the same inherited annotation; DNA correspondence is checked independently.".into())
            },
        };
        row.findings.push(finding(
            Rule::AnnotationCorrespondence,
            row.annotation_correspondence.verdict.clone(),
            row.annotation_correspondence.reason.clone(),
            vec![],
            true,
        ));
        let mut pool = CausePool::new();
        let mut cells = Vec::new();
        let mut complete = true;
        for ordinal in after.offset..after.path.len() {
            if !work.charge(1) {
                complete = false;
                break;
            }
            let call = after
                .path
                .call_at(ordinal)
                .expect("admitted coding ordinal");
            let mapped = session.map_call(&pair.product, &pair.source, call, &mut work);
            complete &= mapped.complete;
            let kind = if !mapped.complete {
                CellKind::Unknown
            } else if mapped.introduced && mapped.parent_calls.is_empty() && !mapped.foreign {
                CellKind::Introduced
            } else if mapped.foreign && mapped.parent_calls.is_empty() && !mapped.introduced {
                CellKind::Foreign
            } else if !mapped.foreign && !mapped.introduced && !mapped.parent_calls.is_empty() {
                let mut chosen = None;
                let mut uncertain = false;
                for parent in &mapped.parent_calls {
                    match before
                        .path
                        .ordinal_of(*parent)
                        .filter(|index| *index >= before.offset)
                    {
                        Some(index) => {
                            let candidate = ParentCall {
                                ordinal: index - before.offset,
                                call: *parent,
                            };
                            if chosen.is_some_and(|old: ParentCall| {
                                old.ordinal != candidate.ordinal || old.call != candidate.call
                            }) {
                                uncertain = true;
                            }
                            chosen = Some(candidate);
                        }
                        None => uncertain = true,
                    }
                }
                if uncertain {
                    CellKind::Outside
                } else {
                    chosen.map_or(CellKind::Unknown, CellKind::Parent)
                }
            } else {
                CellKind::Unknown
            };
            if matches!(kind, CellKind::Unknown | CellKind::Outside) {
                complete = false;
            }
            let causes = match pool.intern(mapped.causes, &mut work) {
                Ok(index) => index,
                Err(_) => {
                    complete = false;
                    break;
                }
            };
            cells.push(Cell {
                kind,
                changed: mapped.changed,
                causes,
            });
        }
        complete &= cells.len() == after.coding_len();
        row.mapping_complete = complete;
        row.findings.push(finding(Rule::CodingCorrespondence, if complete { Verdict::Informational } else { Verdict::Indeterminate },
            if complete { "Native call ancestry was traversed completely within bounds; other-source and outside-CDS calls are not silently paired by similarity." } else { "Native call ancestry or the shared comparison budget is incomplete; no exhaustive absence or unique mapping is asserted." }, vec![], true));
        let register = super::register::inspect(
            &before,
            &after,
            &cells,
            &pool,
            request,
            &mut output,
            &mut work,
        );
        row.findings.push(finding(
            Rule::JunctionRegister,
            register.verdict.clone(),
            register.reason.clone(),
            vec![],
            true,
        ));
        let stops = super::stops::inspect(
            &before,
            &after,
            &cells,
            &pool,
            &register,
            native_annotation,
            request,
            &mut output,
            &mut work,
        );
        row.register_segments = register.segments;
        row.new_stops = stops.evidence;
        row.findings.push(finding(
            Rule::NewInternalStop,
            stops.verdict,
            stops.reason,
            vec![],
            true,
        ));
        if output.exhausted {
            row.mapping_complete = false;
            row.findings.push(finding(Rule::Scope, Verdict::Indeterminate, "Shared evidence limits refused additional detail; partial evidence does not establish absence.", vec![], false));
        }
        row.findings.push(finding(Rule::Scope, Verdict::Informational, "These are reference-coordinate/codon facts under the paired declarations. Annotation correspondence, physical strand inheritance, expression and experimental outcomes are separate questions.", vec![], false));
        rows.push(row);
    }
    let mut report = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision: contract::context(&fingerprint, request)?,
        workspace_fingerprint: fingerprint,
        request_revision: contract::request_revision(request)?,
        pairs: rows,
        work: work.used(),
        coding_calls: request.limits.max_coding_calls - coding_remaining,
        history_events: history.history_events(),
        statement: STATEMENT.into(),
        report_revision: String::new(),
    };
    report.report_revision = contract::digest(&report)?;
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
