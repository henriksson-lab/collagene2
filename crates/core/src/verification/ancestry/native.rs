//! Observes the one native reconstruction without changing its result semantics.
use super::*;
use crate::verification::projection::NativeProjection;

/// A trace failure is retained separately from ordinary scientific reconstruction.
/// This owner cannot turn a failed receipt into an apparently complete causal trace.
pub(crate) struct NativeStep<'a> {
    step: Option<StepBuilder<'a>>,
    error: Option<Error>,
    occurrence: u32,
}
impl<'a> NativeStep<'a> {
    pub fn new(builder: Option<&'a mut TraceBuilder>, index: usize) -> Self {
        match builder.map(|builder| builder.begin_step(index)).transpose() {
            Ok(step) => Self {
                step,
                error: None,
                occurrence: 0,
            },
            Err(error) => Self {
                step: None,
                error: Some(error),
                occurrence: 0,
            },
        }
    }
    pub fn budget(&self) -> Option<TraceBudget> {
        if self.error.is_some() {
            return None;
        }
        self.step.as_ref().map(|step| step.builder.budget.clone())
    }
    fn apply<T>(&mut self, f: impl FnOnce(&mut StepBuilder<'a>) -> Result<T, Error>) -> Option<T> {
        if self.error.is_some() {
            return None;
        }
        let step = self.step.as_mut()?;
        match f(step) {
            Ok(value) => Some(value),
            Err(error) => {
                self.error = Some(error);
                None
            }
        }
    }
    pub fn source(&mut self, source: &Source) -> Option<SnapshotToken> {
        self.apply(|step| {
            step.builder
                .trace
                .snapshots
                .iter()
                .find(|snapshot| {
                    snapshot.identity.as_ref().is_some_and(|identity| {
                        identity.scope.is_empty() && matches!(&identity.source, SourceHandle::Plan(handle) if handle == source)
                    })
                })
                .map(|snapshot| snapshot.token)
                .ok_or_else(|| bad("actual plan source has no committed native snapshot"))
        })
    }
    pub fn destination(&mut self, molecule: &Molecule) -> Option<SnapshotToken> {
        self.apply(|step| step.destination(molecule))
    }
    pub fn seal(&mut self, token: Option<SnapshotToken>, molecule: &Molecule) {
        self.apply(|step| {
            step.seal(
                token.ok_or_else(|| bad("native snapshot token is unavailable"))?,
                molecule,
            )
        });
    }
    #[allow(clippy::too_many_arguments)]
    pub fn project(
        &mut self,
        source: Option<SnapshotToken>,
        source_molecule: &Molecule,
        destination: &mut Option<SnapshotToken>,
        destination_molecule: &Molecule,
        segments: &[MappingSegment],
        native: &NativeProjection,
        kind: ProjectionKind,
        changes: &[MappedChange],
    ) {
        let occurrence = self.occurrence;
        self.occurrence = match self.occurrence.checked_add(1) {
            Some(next) => next,
            None => {
                self.error = Some(Error::TooLarge);
                return;
            }
        };
        self.apply(|step| {
            let source = source.ok_or_else(|| bad("actual source token is unavailable"))?;
            let token = match *destination {
                Some(token) => token,
                None => {
                    let token = step.destination(destination_molecule)?;
                    *destination = Some(token);
                    token
                }
            };
            step.remap_features(token, &native.previous_feature_indices)?;
            step.project(ProjectionCapture {
                source,
                source_molecule,
                destination: token,
                destination_molecule,
                occurrence,
                kind,
                segments,
                projection: &native.projection,
                destination_feature_indices: &native.destination_feature_indices,
                changes,
            })
        });
    }
    pub fn introduce(
        &mut self,
        token: Option<SnapshotToken>,
        start: usize,
        end: usize,
        kind: IntroductionKind,
    ) {
        self.apply(|step| {
            step.introduce(
                token.ok_or_else(|| bad("introduced-call snapshot is unavailable"))?,
                start,
                end,
                kind,
            )
        });
    }
    #[allow(clippy::too_many_arguments)]
    pub fn reference(
        &mut self,
        source: Option<SnapshotToken>,
        source_molecule: &Molecule,
        destination: &mut Option<SnapshotToken>,
        destination_molecule: &Molecule,
        receipt: &crate::verification::projection::ReferenceProjection,
        kind: ProjectionKind,
        changes: &[MappedChange],
        introduced_kind: IntroductionKind,
    ) {
        self.project(
            source,
            source_molecule,
            destination,
            destination_molecule,
            &receipt.map,
            &receipt.native,
            kind,
            changes,
        );
        for &(start, end) in &receipt.introduced {
            if start < end {
                self.introduce(*destination, start, end, introduced_kind.clone());
            }
        }
    }
    /// The native caller has just cloned this exact source without changing its
    /// DNA or annotations. No source is discovered by content equality.
    pub fn unchanged_reference(
        &mut self,
        source: Option<SnapshotToken>,
        source_molecule: &Molecule,
        destination: &mut Option<SnapshotToken>,
        destination_molecule: &Molecule,
    ) {
        if self.step.is_none() || self.error.is_some() {
            return;
        }
        let map = [MappingSegment {
            source_start: 0,
            source_end: source_molecule.sequence.len(),
            destination_start: 0,
            reverse: false,
        }];
        let projection = match collagene_model::annotations::project_features_with_topology(
            source_molecule,
            &map,
            destination_molecule.sequence.len(),
            destination_molecule.topology,
        ) {
            Ok(projection) => projection,
            Err(error) => {
                self.error = Some(error);
                return;
            }
        };
        let native = NativeProjection {
            destination_feature_indices: (0..projection.features.len()).collect(),
            previous_feature_indices: vec![],
            projection,
        };
        self.project(
            source,
            source_molecule,
            destination,
            destination_molecule,
            &map,
            &native,
            ProjectionKind::ReferenceCorrespondence,
            &[],
        );
    }
    pub fn sources(
        &mut self,
        bindings: &[(&Source, &Molecule)],
        destination: &mut Option<SnapshotToken>,
        molecule: &Molecule,
        receipts: Option<&[crate::verification::projection::SourceProjection]>,
        kind: ProjectionKind,
    ) {
        if self.step.is_none() || self.error.is_some() {
            return;
        }
        let Some(receipts) = receipts else {
            self.unavailable(
                "native graph does not provide a complete current-namespace coordinate receipt",
            );
            return;
        };
        for receipt in receipts {
            let Some((source, source_molecule)) = bindings.get(receipt.input) else {
                self.unavailable(
                    "native source ordinal exceeds its actual explicit operand bindings",
                );
                return;
            };
            let source_token = self.source(source);
            self.reference(
                source_token,
                source_molecule,
                destination,
                molecule,
                &receipt.reference,
                kind.clone(),
                &[],
                IntroductionKind::NewlySynthesizedUnmappedCalls,
            );
        }
    }
    pub fn graph_sources(
        &mut self,
        bindings: &[(&Source, &Molecule)],
        destination: &mut Option<SnapshotToken>,
        molecule: &Molecule,
        receipts: Option<crate::verification::projection::GraphProjections>,
    ) {
        use crate::verification::projection::GraphProjectionSource;
        use std::borrow::Cow;
        if self.step.is_none() || self.error.is_some() {
            return;
        }
        let Some(receipts) = receipts else {
            self.unavailable("native graph coordinate receipt retention is incomplete");
            return;
        };
        if receipts.sources.len() > 64
            || receipts.retained.len() > 64
            || receipts.projections.len() > 4096
        {
            self.unavailable("native graph capture exceeds bounded root/projection counts");
            return;
        }
        let mut retained = Vec::with_capacity(receipts.retained.len());
        for replay in receipts.retained {
            let trace = match replay.trace {
                Ok(trace) => trace,
                Err(_) => {
                    self.unavailable("native graph cached material trace is incomplete");
                    return;
                }
            };
            let Some(actual) = replay.actual_reference else {
                self.unavailable("native graph cached trace lacks its exact native reference");
                return;
            };
            let token = self.import(replay.binding, trace, replay.selected_step);
            retained.push((token, actual));
        }
        let mut roots = Vec::with_capacity(receipts.sources.len());
        for source in receipts.sources {
            let (token, actual) = match source {
                GraphProjectionSource::CurrentInput { index } => match bindings.get(index) {
                    Some((source, molecule)) => (self.source(source), Cow::Borrowed(*molecule)),
                    None => {
                        self.unavailable(
                            "native root input ordinal exceeds actual caller bindings",
                        );
                        return;
                    }
                },
                GraphProjectionSource::RetainedGraph {
                    role,
                    scopes,
                    retained_index,
                    molecule,
                } => {
                    let token = if let Some(index) = retained_index {
                        let Some((token, actual)) = retained.get(index) else {
                            self.unavailable(
                                "native graph root names an absent cached replay occurrence",
                            );
                            return;
                        };
                        let mut presented = None;
                        self.unchanged_reference(*token, actual, &mut presented, &molecule);
                        self.seal(presented, &molecule);
                        self.apply(|step| {
                            step.bind_topo_source(
                                presented.ok_or_else(|| {
                                    bad("native graph presentation token is unavailable")
                                })?,
                                &role,
                                &molecule,
                                &scopes,
                            )
                        });
                        presented
                    } else {
                        self.apply(|step| step.register_topo_source(&role, &molecule, &scopes))
                    };
                    (token, Cow::Owned(molecule))
                }
            };
            roots.push((token, actual));
        }
        for receipt in receipts.projections {
            let Some((token, source)) = roots.get(receipt.input) else {
                self.unavailable("native graph receipt source exceeds actual expansion roots");
                return;
            };
            self.reference(
                *token,
                source.as_ref(),
                destination,
                molecule,
                &receipt.reference,
                ProjectionKind::ReferenceCorrespondence,
                &[],
                IntroductionKind::NewlySynthesizedUnmappedCalls,
            );
        }
    }
    pub fn import(
        &mut self,
        binding: OriginBinding,
        trace: ReconstructionTrace,
        selected_step: usize,
    ) -> Option<SnapshotToken> {
        let occurrence = self.occurrence;
        self.occurrence = match self.occurrence.checked_add(1) {
            Some(next) => next,
            None => {
                self.error = Some(Error::TooLarge);
                return None;
            }
        };
        self.apply(|step| step.import(binding, trace, selected_step, occurrence))
    }
    pub fn bind_reference_event(
        &mut self,
        token: Option<SnapshotToken>,
        event: &crate::material::lineage::EventBinding,
        molecule: &Molecule,
    ) {
        self.apply(|step| {
            step.bind_reference_event(
                token.ok_or_else(|| bad("native retrieved event snapshot is unavailable"))?,
                event,
                molecule,
            )
        });
    }
    /// Capture the actual saved reference before the enclosing owner appends its
    /// own metadata. Event identity is supplied by that validated native owner.
    pub fn present_event(
        &mut self,
        source: Option<SnapshotToken>,
        actual: &Molecule,
        event: &crate::material::lineage::EventBinding,
        destination: &mut Option<SnapshotToken>,
        presented: &Molecule,
    ) {
        let mut event_snapshot = None;
        self.unchanged_reference(source, actual, &mut event_snapshot, presented);
        self.seal(event_snapshot, presented);
        self.bind_reference_event(event_snapshot, event, presented);
        self.unchanged_reference(event_snapshot, presented, destination, presented);
    }
    pub fn unavailable(&mut self, reason: &str) {
        if self.step.is_some() && self.error.is_none() {
            self.error = Some(bad(reason));
        }
    }
    /// Finish only after all native product metadata/annotations are final.
    pub fn finish(
        mut self,
        token: Option<SnapshotToken>,
        product: Option<&Molecule>,
    ) -> Result<(), Error> {
        if let Some(error) = self.error.take() {
            return Err(error);
        }
        let Some(mut step) = self.step.take() else {
            return Ok(());
        };
        let Some(product) = product else {
            return Ok(());
        }; // Drop rolls back a failed native operation.
        let token = token
            .ok_or_else(|| bad("reconstructed product has no actual captured coordinate recipe"))?;
        step.seal(token, product)?;
        step.commit(token)?;
        Ok(())
    }
}
