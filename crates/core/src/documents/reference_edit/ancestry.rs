//! Observes one ordered native reference transaction; failures refuse only trace.
use super::{Operation, Receipt};
use crate::{
    material::lineage::EventBinding,
    sequence::editing::AnnotationIdentity,
    verification::{
        ancestry::{
            IntroductionKind, ProjectionKind, ReconstructionTrace, SnapshotToken, TraceBudget,
            TraceBuilder, native::NativeStep,
        },
        projection::{NativeProjection, ReferenceProjection},
    },
};
use collagene_model::{
    Error, Molecule,
    annotations::{MappingSegment, Projection, ReplacementProjection},
};

pub(super) struct Recorder {
    builder: Option<TraceBuilder>,
    error: Option<Error>,
    budget: TraceBudget,
    current: Option<SnapshotToken>,
    donors: Vec<Option<SnapshotToken>>,
}
impl Recorder {
    pub fn new(
        source: &Molecule,
        event: &EventBinding,
        actions: &[Operation],
        budget: TraceBudget,
    ) -> Self {
        let prepared = (|| {
            let revision = Operation::actions_revision(event, actions)?;
            let mut builder =
                TraceBuilder::for_reference_edits(budget.clone(), event, &revision, actions.len())?;
            let source = builder.register_reference_source(event, source)?;
            let mut donors = vec![None; actions.len()];
            for (index, action) in actions.iter().enumerate() {
                if let Operation::PasteReference {
                    source,
                    source_event,
                    ..
                } = action
                {
                    donors[index] = Some(builder.register_reference_donor(source_event, source)?);
                }
            }
            Ok::<_, Error>((builder, source, donors))
        })();
        match prepared {
            Ok((builder, current, donors)) => Self {
                builder: Some(builder),
                error: None,
                budget,
                current: Some(current),
                donors,
            },
            Err(error) => Self {
                builder: None,
                error: Some(error),
                budget,
                current: None,
                donors: vec![],
            },
        }
    }
    fn refuse(&mut self, error: Error) {
        self.builder = None;
        self.current = None;
        if self.error.is_none() {
            self.error = Some(error);
        }
    }
    /// Charge actual evidence copying before native helpers retain projections.
    pub fn capture(&mut self, source: &Molecule) -> bool {
        if self.builder.is_none() {
            return false;
        }
        match self.budget.retain_model(source) {
            Ok(()) => true,
            Err(error) => {
                self.refuse(error);
                false
            }
        }
    }
    pub fn record(
        &mut self,
        index: usize,
        action: &Operation,
        source: &Molecule,
        product: &Molecule,
        receipt: Option<Receipt>,
    ) {
        if self.builder.is_none() {
            return;
        }
        if let Err(error) = self.budget.retain_model(product) {
            self.refuse(error);
            return;
        }
        let Some(receipt) = receipt else {
            self.refuse(Error::Unsupported(
                "native reference owner did not retain a causal receipt".into(),
            ));
            return;
        };
        let prepared = (|| {
            Ok::<_, Error>(match receipt {
                Receipt::Retained => Observed::Single(projected(
                    source,
                    product,
                    crate::sequence::annotations::retained_reference_projection(source),
                )?),
                Receipt::Projection(receipt) => {
                    Observed::Single(projected(source, product, receipt)?)
                }
                Receipt::Replacement(receipt) => {
                    Observed::Single(replacement(product, receipt, true)?)
                }
                Receipt::Annotations(receipt) => {
                    Observed::Single(annotation(source, product, receipt)?)
                }
                Receipt::Paste { destination, donor } => Observed::Paste {
                    destination: replacement(product, destination, false)?,
                    donor: match action {
                        Operation::PasteReference { source, .. } => {
                            projected(source, product, donor)?
                        }
                        _ => {
                            return Err(Error::Format(
                                "paste receipt supplied by another native action".into(),
                            ));
                        }
                    },
                },
            })
        })();
        let receipt = match prepared {
            Ok(receipt) => receipt,
            Err(error) => {
                self.refuse(error);
                return;
            }
        };
        let Some(builder) = self.builder.as_mut() else {
            return;
        };
        let mut native = NativeStep::new(Some(builder), index);
        let mut destination = None;
        match receipt {
            Observed::Single(receipt) => native.reference(
                self.current,
                source,
                &mut destination,
                product,
                &receipt,
                ProjectionKind::ReferenceCorrespondence,
                &[],
                IntroductionKind::DeclaredLiteral,
            ),
            Observed::Paste {
                destination: receipt,
                donor,
            } => {
                native.reference(
                    self.current,
                    source,
                    &mut destination,
                    product,
                    &receipt,
                    ProjectionKind::ReferenceCorrespondence,
                    &[],
                    IntroductionKind::DeclaredLiteral,
                );
                if let Operation::PasteReference { source, .. } = action {
                    native.reference(
                        self.donors.get(index).copied().flatten(),
                        source,
                        &mut destination,
                        product,
                        &donor,
                        ProjectionKind::ReferenceCorrespondence,
                        &[],
                        IntroductionKind::DeclaredLiteral,
                    );
                }
            }
        }
        match native.finish(destination, Some(product)) {
            Ok(()) => self.current = destination,
            Err(error) => self.refuse(error),
        }
    }
    pub fn finish(self) -> Result<ReconstructionTrace, Error> {
        match (self.error, self.builder) {
            (Some(error), _) => Err(error),
            (None, Some(builder)) => Ok(builder.finish()),
            (None, None) => Err(Error::Unsupported(
                "native reference trace unavailable".into(),
            )),
        }
    }
}
enum Observed {
    Single(ReferenceProjection),
    Paste {
        destination: ReferenceProjection,
        donor: ReferenceProjection,
    },
}
fn replacement(
    product: &Molecule,
    receipt: ReplacementProjection,
    introduce_literal: bool,
) -> Result<ReferenceProjection, Error> {
    let count = receipt.source_feature_indices.len();
    let features = product
        .features
        .get(..count)
        .ok_or_else(|| {
            Error::Format("native replacement ordinals exceed retained features".into())
        })?
        .to_vec();
    let introduced = if introduce_literal && !receipt.unchanged && receipt.inserted_length > 0 {
        vec![(receipt.start, receipt.start + receipt.inserted_length)]
    } else {
        vec![]
    };
    Ok(ReferenceProjection {
        map: receipt.map,
        native: NativeProjection {
            projection: Projection {
                features,
                source_feature_indices: receipt.source_feature_indices,
                warnings: vec![],
                unmapped_features: vec![],
                unmapped_source_feature_indices: vec![],
                dropped_features: receipt.dropped_source_feature_indices.len(),
                dropped_source_feature_indices: receipt.dropped_source_feature_indices,
            },
            destination_feature_indices: (0..count).collect(),
            previous_feature_indices: vec![],
        },
        introduced,
    })
}
fn annotation(
    source: &Molecule,
    product: &Molecule,
    identity: AnnotationIdentity,
) -> Result<ReferenceProjection, Error> {
    let mut features = Vec::with_capacity(identity.retained.len());
    let mut source_indices = Vec::with_capacity(identity.retained.len());
    let mut destination_indices = Vec::with_capacity(identity.retained.len());
    for (source_index, destination_index) in identity.retained {
        features.push(
            product
                .features
                .get(destination_index)
                .ok_or_else(|| {
                    Error::Format("native annotation destination ordinal differs".into())
                })?
                .clone(),
        );
        source_indices.push(source_index);
        destination_indices.push(destination_index);
    }
    let mut dropped = identity.reassigned;
    dropped.extend(identity.removed);
    let map = if source.sequence.is_empty() {
        vec![]
    } else {
        vec![MappingSegment {
            source_start: 0,
            source_end: source.sequence.len(),
            destination_start: 0,
            reverse: false,
        }]
    };
    Ok(ReferenceProjection {
        map,
        native: NativeProjection {
            projection: Projection {
                features,
                source_feature_indices: source_indices,
                warnings: vec![],
                unmapped_features: vec![],
                unmapped_source_feature_indices: vec![],
                dropped_features: dropped.len(),
                dropped_source_feature_indices: dropped,
            },
            destination_feature_indices: destination_indices,
            previous_feature_indices: vec![],
        },
        introduced: vec![],
    })
}

fn projected(
    source: &Molecule,
    product: &Molecule,
    receipt: crate::sequence::annotations::ReferenceProjectionReceipt,
) -> Result<ReferenceProjection, Error> {
    let mut features = Vec::with_capacity(receipt.destination_feature_indices.len());
    for index in &receipt.destination_feature_indices {
        features.push(
            product
                .features
                .get(*index)
                .ok_or_else(|| {
                    Error::Format("native transform destination ordinal differs".into())
                })?
                .clone(),
        );
    }
    let mut unmapped = Vec::with_capacity(receipt.unmapped_source_feature_indices.len());
    if receipt.unmapped_source_feature_indices.len() != receipt.unmapped_review.len() {
        return Err(Error::Format(
            "native unmapped annotation review dispositions differ".into(),
        ));
    }
    for (index, needs_review) in receipt
        .unmapped_source_feature_indices
        .iter()
        .zip(&receipt.unmapped_review)
    {
        let mut evidence = source
            .features
            .get(*index)
            .ok_or_else(|| Error::Format("native unmapped source ordinal differs".into()))?
            .clone();
        evidence.needs_review = *needs_review;
        unmapped.push(evidence);
    }
    Ok(ReferenceProjection {
        map: receipt.map,
        native: NativeProjection {
            projection: Projection {
                features,
                source_feature_indices: receipt.source_feature_indices,
                warnings: vec![],
                unmapped_features: unmapped,
                unmapped_source_feature_indices: receipt.unmapped_source_feature_indices,
                dropped_features: receipt.dropped_source_feature_indices.len(),
                dropped_source_feature_indices: receipt.dropped_source_feature_indices,
            },
            destination_feature_indices: receipt.destination_feature_indices,
            previous_feature_indices: receipt.previous_feature_indices,
        },
        introduced: vec![],
    })
}
