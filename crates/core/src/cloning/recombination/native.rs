//! Private receipts from actual annotation append operations, never replayed maps.
use super::{AnnotationLimits, Selection, Work};
use crate::verification::projection::{NativeProjection, ReferenceProjection, SourceProjection};
use collagene_model::{
    Error,
    annotations::{MappingSegment, Projection},
};

pub(super) struct Capture {
    pub selection: Selection,
    pub complete: bool,
    pub projections: Vec<SourceProjection>,
    bytes: usize,
    work: usize,
}
impl Capture {
    pub fn new(selection: Selection) -> Self {
        Self {
            selection,
            complete: true,
            projections: vec![],
            bytes: 0,
            work: 0,
        }
    }
    pub fn append(
        &mut self,
        input: usize,
        map: &[MappingSegment],
        projection: &Projection,
        previous: usize,
        limits: &AnnotationLimits,
        _work: &mut Work,
    ) -> Result<(), Error> {
        // Receipt admission has its own bounded ledger. Refusal must not alter
        // the scientific batch's annotation counters or reconstructed products.
        if !self.complete {
            return Ok(());
        }
        let admitted = (|| {
            let bytes = super::size(
                &(
                    &projection.features,
                    &projection.unmapped_features,
                    &projection.source_feature_indices,
                    &projection.unmapped_source_feature_indices,
                    &projection.dropped_source_feature_indices,
                    &projection.warnings,
                    map,
                ),
                limits.max_copy_bytes,
            )?;
            let ordinals = previous
                .checked_add(projection.features.len())
                .ok_or(Error::TooLarge)?;
            let bytes = bytes
                .checked_add(
                    ordinals
                        .checked_mul(std::mem::size_of::<usize>())
                        .ok_or(Error::TooLarge)?,
                )
                .ok_or(Error::TooLarge)?;
            let bytes = self
                .bytes
                .checked_add(bytes)
                .filter(|n| *n <= limits.max_copy_bytes.min(8 << 20))
                .ok_or(Error::TooLarge)?;
            let work = self
                .work
                .checked_add(map.len() + projection.features.len())
                .filter(|n| *n <= limits.max_work)
                .ok_or(Error::TooLarge)?;
            Ok::<_, Error>((bytes, work, ordinals))
        })();
        let (bytes, units, ordinals) = match admitted {
            Ok(value) => value,
            Err(_) => {
                self.complete = false;
                self.projections.clear();
                return Ok(());
            }
        };
        self.bytes = bytes;
        self.work = units;
        self.projections.push(SourceProjection {
            input,
            reference: ReferenceProjection {
                map: map.to_vec(),
                introduced: vec![],
                native: NativeProjection {
                    destination_feature_indices: (previous..ordinals).collect(),
                    previous_feature_indices: (0..previous).collect(),
                    projection: projection.clone(),
                },
            },
        });
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::*;
    #[test]
    fn capture_keeps_whole_batch_identical_and_records_selected_literal_feature_ordinals() {
        let molecules = crate::formats::sequence::import_genbank(include_str!(
            "../../../../../fixtures/recombination-bp-demo.gb"
        ))
        .unwrap()
        .molecules;
        let inputs = vec![
            reference::Input {
                role: "insert".into(),
                molecule: molecules[0].clone(),
            },
            reference::Input {
                role: "donor".into(),
                molecule: molecules[1].clone(),
            },
        ];
        let request = Request::new(
            ProfileSnapshot::manufacturer_context().unwrap(),
            reference::Reaction::BP,
        );
        let ordinary = simulate(&inputs, &request).unwrap();
        let selection = ordinary
            .documents
            .iter()
            .find(|document| document.molecule.sequence.len() == 2553)
            .unwrap()
            .selection
            .clone();
        let (captured, receipts) =
            simulate_with_projections(&inputs, &request, &selection).unwrap();
        assert_eq!(
            serde_json::to_value(&ordinary).unwrap(),
            serde_json::to_value(&captured).unwrap()
        );
        let receipts = receipts.unwrap();
        assert_eq!(receipts.len(), 2);
        let insert = receipts.iter().find(|receipt| receipt.input == 0).unwrap();
        let offset = insert
            .reference
            .native
            .projection
            .source_feature_indices
            .iter()
            .position(|index| *index == 2)
            .unwrap();
        let destination = insert.reference.native.destination_feature_indices[offset];
        let selected = captured
            .documents
            .iter()
            .find(|document| document.selection == selection)
            .unwrap();
        // Literal crossover oracle: the retained insert is b[8:51], and
        // donor backbone is d[2902:] + d[:651]. The native circular graph
        // starts at the first retained insert edge. The separately authored
        // fixture target starts with the donor prefix instead, a651-base
        // rotation; its annotation coordinates are not this product's origin.
        let b = inputs[0].molecule.sequence.as_str();
        let d = inputs[1].molecule.sequence.as_str();
        let expected = format!("{}{}{}", &b[8..51], &d[2902..], &d[..651]);
        assert_eq!(
            selected.molecule.topology,
            collagene_model::Topology::Circular
        );
        assert_eq!(selected.molecule.sequence.as_str(), expected);
        assert_eq!(&b[25..34], "ATGCCGTAA");
        assert_eq!(&selected.molecule.sequence.as_str()[17..26], "ATGCCGTAA");
        assert_eq!(
            inputs[0].molecule.features[2].location.spans(),
            [(25, 34, false)]
        );
        assert_eq!(
            selected.molecule.features[destination].location.spans(),
            vec![(17, 26, false)]
        );
        assert_eq!(
            insert.reference.native.projection.features[offset],
            selected.molecule.features[destination]
        );
        // Verify every literal insert-feature call through the captured native
        // source map, not a feature label or a recomputed content match.
        for position in 25..34 {
            let map = insert
                .reference
                .map
                .iter()
                .find(|map| map.source_start <= position && position < map.source_end)
                .unwrap();
            assert!(!map.reverse);
            let mapped = map.destination_start + position - map.source_start;
            assert_eq!(mapped, position - 8);
            assert_eq!(
                b.as_bytes()[position],
                selected.molecule.sequence.as_str().as_bytes()[mapped]
            );
        }
    }
}
