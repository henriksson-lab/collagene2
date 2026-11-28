//! Receipts retained during the one native derivation, with explicit recipe roles.
use super::{SourceOrigin, SourceRole};
use crate::{
    material::run::RetainedReplay,
    verification::{
        ancestry::{
            IntroductionKind, MappedChange, ProjectionKind, ReconstructionTrace, TraceBudget,
            TraceBuilder, native::NativeStep,
        },
        projection::ReferenceProjection,
    },
};
use collagene_model::{Error, Molecule};

struct Root {
    role: SourceRole,
    molecule: Molecule,
    retained: Option<RetainedReplay>,
}
struct Applied {
    input: usize,
    projection: ReferenceProjection,
    kind: ProjectionKind,
    changes: Vec<MappedChange>,
    introduced: IntroductionKind,
}
pub(super) struct Capture {
    budget: TraceBudget,
    context_revision: String,
    roots: Vec<Root>,
    applied: Vec<Applied>,
    unavailable: Option<&'static str>,
    oe: Option<crate::cloning::overlap::extension::execution::NativeStages>,
}
impl Capture {
    pub fn new(budget: TraceBudget, context_revision: &str) -> Self {
        Self {
            budget,
            context_revision: context_revision.into(),
            roots: vec![],
            applied: vec![],
            unavailable: None,
            oe: None,
        }
    }
    pub fn source(
        &mut self,
        role: SourceRole,
        source: &Molecule,
        origin: &SourceOrigin,
        replay: &mut crate::verification::replay_budget::ReplayBudget<'_>,
    ) -> Result<super::FreshMaterial, Error> {
        let (fresh, retained) = match origin {
            SourceOrigin::OriginalReference => (super::fresh_source(source, origin)?, None),
            SourceOrigin::Saved(origin) => {
                let (fresh, retained) =
                    origin.replay_scoped(source, replay, Some(self.budget.clone()))?;
                (super::FreshMaterial::from_run_source(fresh)?, retained)
            }
        };
        self.root(role, source, retained)?;
        Ok(fresh)
    }
    pub fn root(
        &mut self,
        role: SourceRole,
        source: &Molecule,
        retained: Option<RetainedReplay>,
    ) -> Result<(), Error> {
        if self.roots.len() >= super::MAX_SOURCES || self.roots.iter().any(|root| root.role == role)
        {
            return Err(Error::TooLarge);
        }
        if self.unavailable.is_some() {
            return Ok(());
        }
        if self.budget.retain_model(source).is_err() {
            self.unavailable("Shared native receipt source-copy budget was exhausted");
            return Ok(());
        }
        self.roots.push(Root {
            role,
            molecule: source.clone(),
            retained,
        });
        Ok(())
    }
    pub fn record(
        &mut self,
        input: usize,
        projection: ReferenceProjection,
        kind: ProjectionKind,
        changes: Vec<MappedChange>,
        introduced: IntroductionKind,
    ) {
        if self.unavailable.is_some() {
            return;
        }
        self.applied.push(Applied {
            input,
            projection,
            kind,
            changes,
            introduced,
        });
    }
    /// The terminal owner copied these source features verbatim while retaining
    /// every original top-reference coordinate. This records that actual copy.
    pub fn record_top_copy(&mut self, source: &Molecule, product: &Molecule) -> Result<(), Error> {
        if self.unavailable.is_some() {
            return Ok(());
        }
        if self.budget.retain_model(source).is_err() {
            self.unavailable("Shared native receipt feature-copy budget was exhausted");
            return Ok(());
        }
        let count = source.features.len();
        let map = vec![collagene_model::annotations::MappingSegment {
            source_start: 0,
            source_end: source.sequence.len(),
            destination_start: 0,
            reverse: false,
        }];
        let projection = collagene_model::annotations::Projection {
            features: source.features.clone(),
            source_feature_indices: (0..count).collect(),
            warnings: vec![],
            unmapped_features: vec![],
            unmapped_source_feature_indices: vec![],
            dropped_features: 0,
            dropped_source_feature_indices: vec![],
        };
        self.record(
            0,
            ReferenceProjection {
                map,
                native: crate::verification::projection::NativeProjection {
                    projection,
                    destination_feature_indices: (0..count).collect(),
                    previous_feature_indices: vec![],
                },
                introduced: if source.sequence.len() < product.sequence.len() {
                    vec![(source.sequence.len(), product.sequence.len())]
                } else {
                    vec![]
                },
            },
            ProjectionKind::ReferenceCorrespondence,
            vec![],
            IntroductionKind::DeclaredLiteral,
        );
        Ok(())
    }
    pub fn record_oe(
        &mut self,
        native: Option<crate::cloning::overlap::extension::execution::NativeStages>,
    ) {
        if native.is_none() {
            self.unavailable("Overlap-extension native stage receipt retention was refused");
        }
        self.oe = native;
    }
    pub fn unavailable(&mut self, reason: &'static str) {
        self.unavailable = Some(reason);
    }
    pub fn finish(self, product: &Molecule) -> Result<ReconstructionTrace, Error> {
        if let Some(reason) = self.unavailable {
            return Err(Error::Unsupported(reason.into()));
        }
        let mut builder = TraceBuilder::for_material(self.budget, &self.context_revision)?;
        let mut tokens = Vec::with_capacity(self.roots.len());
        for root in &self.roots {
            tokens.push(if root.retained.is_none() {
                Some(builder.register_material_source(root.role, &root.molecule)?)
            } else {
                None
            });
        }
        let mut step = NativeStep::new(Some(&mut builder), 0);
        let mut roots = Vec::with_capacity(self.roots.len());
        for (index, root) in self.roots.into_iter().enumerate() {
            if let Some(retained) = root.retained {
                let child = retained.trace?;
                let actual = retained.actual_reference.ok_or_else(|| {
                    Error::Format("Native retained trace has no actual reference model".into())
                })?;
                let event = match &retained.binding {
                    crate::verification::ancestry::OriginBinding::Run { selection, .. } => {
                        selection.event.clone()
                    }
                    crate::verification::ancestry::OriginBinding::Material { product, .. } => {
                        product.clone()
                    }
                    _ => {
                        return Err(Error::Format(
                            "Material replay returned a different native origin kind".into(),
                        ));
                    }
                };
                let token = step.import(retained.binding, child, retained.selected_step);
                let mut presented = None;
                step.unchanged_reference(token, &actual, &mut presented, &root.molecule);
                step.seal(presented, &root.molecule);
                step.bind_reference_event(presented, &event, &root.molecule);
                tokens[index] = presented;
            }
            roots.push(root.molecule);
        }
        let mut destination = None;
        if let Some(oe) = self.oe {
            let source = roots
                .first()
                .ok_or_else(|| Error::Format("Missing actual OE template role".into()))?;
            let mut left = None;
            step.reference(
                tokens[0],
                source,
                &mut left,
                &oe.left,
                &oe.left_projection,
                ProjectionKind::TemplateSynthesis,
                &[],
                IntroductionKind::PrimerTail,
            );
            step.seal(left, &oe.left);
            let mut right = None;
            step.reference(
                tokens[0],
                source,
                &mut right,
                &oe.right,
                &oe.right_projection,
                ProjectionKind::TemplateSynthesis,
                &[],
                IntroductionKind::PrimerTail,
            );
            step.seal(right, &oe.right);
            let mut extended = None;
            let mut projections = oe.extension_projections.into_iter();
            let left_projection = projections
                .next()
                .ok_or_else(|| Error::Format("Missing native left OE projection".into()))?;
            let right_projection = projections
                .next()
                .ok_or_else(|| Error::Format("Missing native right OE projection".into()))?;
            step.reference(
                left,
                &oe.left,
                &mut extended,
                &oe.extended,
                &left_projection,
                ProjectionKind::HomologyCollapse,
                &[],
                IntroductionKind::NewlySynthesizedUnmappedCalls,
            );
            step.reference(
                right,
                &oe.right,
                &mut extended,
                &oe.extended,
                &right_projection,
                ProjectionKind::HomologyCollapse,
                &[],
                IntroductionKind::NewlySynthesizedUnmappedCalls,
            );
            step.seal(extended, &oe.extended);
            step.reference(
                extended,
                &oe.extended,
                &mut destination,
                product,
                &oe.outer_projection,
                ProjectionKind::TemplateSynthesis,
                &[],
                IntroductionKind::PrimerTail,
            );
        }
        for applied in self.applied {
            let source = roots.get(applied.input).ok_or_else(|| {
                Error::Format("Native derivation receipt source ordinal differs".into())
            })?;
            step.reference(
                tokens[applied.input],
                source,
                &mut destination,
                product,
                &applied.projection,
                applied.kind,
                &applied.changes,
                applied.introduced,
            );
        }
        step.finish(destination, Some(product))?;
        Ok(builder.finish())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cloning::ligation::tailing::TailingBase,
        verification::ancestry::{Limits, SourceHandle},
    };
    use collagene_model::{Feature, Location, Topology};
    #[test]
    fn actual_terminal_copy_retains_material_role_and_marks_only_new_tail_calls() {
        let mut source = Molecule::dna("literal source", "ACGT", Topology::Linear).unwrap();
        source.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(1, 3),
            qualifiers: vec![],
            needs_review: false,
        });
        let fresh = super::super::FreshMaterial {
            reference: source.clone(),
            tracked: crate::chemistry::tracked_duplex::TrackedDuplex::exact_pcr(&source.sequence)
                .unwrap(),
            incomplete: false,
            findings: vec![],
            findings_total: 0,
            findings_truncated: false,
            reaction_conditions: vec![],
        };
        let (product, _, _) = super::super::terminal::terminal_from_fresh(
            &source,
            &fresh,
            &super::super::TerminalAction::Tail3Prime {
                base: TailingBase::A,
            },
            "actual tail product",
        )
        .unwrap();
        assert_eq!(product.sequence.as_str(), "ACGTA");
        let mut capture = Capture::new(
            TraceBudget::new(Limits::default()).unwrap(),
            &"a".repeat(64),
        );
        capture.root(SourceRole::Template, &source, None).unwrap();
        capture.record_top_copy(&source, &product).unwrap();
        let trace = capture.finish(&product).unwrap();
        assert!(trace.step(0).unwrap().local_coverage_complete);
        let mut lookup = trace.lookup_session();
        let existing = lookup.at(0, 2, None).unwrap();
        assert!(existing.complete);
        assert_eq!(
            existing.origins[0].source.source,
            SourceHandle::Material(SourceRole::Template)
        );
        assert_eq!(existing.origins[0].position, 2);
        let added = lookup.at(0, 4, None).unwrap();
        assert!(added.complete && added.origins.is_empty());
        assert_eq!(added.introductions.len(), 1);
        let annotation = lookup.feature_at(0, 0, None).unwrap();
        assert!(annotation.complete);
        assert_eq!(annotation.origins[0].feature_index, 0);
    }
}
