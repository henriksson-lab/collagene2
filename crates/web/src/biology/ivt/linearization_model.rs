//! Editable picker policy and exact-event/settings capture, no cleavage algorithm.
use super::model::Draft as IvtDraft;
use crate::restriction::reaction_conditions_model::Draft as ConditionsDraft;
use collagene_core::{
    chemistry::tracked_duplex::GenerationAssumptions,
    cloning::restriction::{catalog, ivt_linearization as native},
    material::{lineage::EventBinding, run::MaterialOrigin},
    rna_chemistry::ivt,
};
use collagene_model::Molecule;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Draft {
    pub enzymes: Vec<String>,
    pub allow_blunt: bool,
    pub allow_five_prime: bool,
    pub allow_three_prime: bool,
    pub require_conditions: bool,
    pub require_resolved_ivt: bool,
    pub require_complete_origin: bool,
    pub comparisons: String,
    pub sites: String,
}
impl Default for Draft {
    fn default() -> Self {
        Self {
            enzymes: vec![],
            allow_blunt: true,
            allow_five_prime: true,
            allow_three_prime: false,
            require_conditions: false,
            require_resolved_ivt: false,
            require_complete_origin: true,
            comparisons: "100000000".into(),
            sites: "64".into(),
        }
    }
}
fn integer(value: &str, name: &str) -> Result<usize, String> {
    if value.is_empty() || value.len() > 20 || !value.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("Enter a whole number for {name}."));
    }
    value
        .parse()
        .map_err(|_| format!("{name} exceeds the supported range."))
}
fn origin_preclone(origin: &MaterialOrigin) -> Result<(), String> {
    struct Remaining(usize);
    impl std::io::Write for Remaining {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0 = self.0.checked_sub(bytes.len()).ok_or_else(|| {
                std::io::Error::other("Saved origin exceeds IVT picker request bound")
            })?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(Remaining(native::MAX_REQUEST_BYTES), origin)
        .map_err(|error| error.to_string())
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Snapshot {
    pub binding: EventBinding,
    pub epoch: u64,
    pub ivt: IvtDraft,
    pub picker: Draft,
    pub generation: GenerationAssumptions,
    pub origin_revision: Option<String>,
    pub conditions: ConditionsDraft,
}
#[derive(Clone)]
pub(crate) struct Capture {
    pub snapshot: Snapshot,
    pub request: native::Request,
    pub revision: String,
}
impl Capture {
    pub fn prepare(
        snapshot: Snapshot,
        source: &Molecule,
        origin: Option<&MaterialOrigin>,
    ) -> Result<Self, String> {
        let ivt = snapshot.ivt.request(source)?;
        if ivt.linearization.is_some() {
            return Err(
                "Clear the supplied cut offsets above before deriving restriction alternatives."
                    .into(),
            );
        }
        if snapshot.binding.revision != source.revision() {
            return Err("The IVT source event changed; inspect its current revision.".into());
        }
        let catalog = catalog::builtin_catalog().map_err(|e| e.to_string())?;
        let selection =
            catalog::select(&catalog, &snapshot.picker.enzymes).map_err(|e| e.to_string())?;
        let conditions = snapshot.conditions.options(&catalog, &selection)?;
        let mut allowed = Vec::new();
        for (enabled, kind) in [
            (snapshot.picker.allow_blunt, ivt::EndKind::Blunt),
            (
                snapshot.picker.allow_five_prime,
                ivt::EndKind::FivePrimeOverhang,
            ),
            (
                snapshot.picker.allow_three_prime,
                ivt::EndKind::ThreePrimeOverhang,
            ),
        ] {
            if enabled {
                allowed.push(kind);
            }
        }
        let source_input = match (snapshot.origin_revision.as_ref(), origin) {
            (Some(_), Some(origin)) if origin.event() == &snapshot.binding => {
                origin_preclone(origin)?;
                native::SourceInput::Saved {
                    origin: origin.clone(),
                }
            }
            (None, None) => native::SourceInput::DeclaredOriginal {
                generation: snapshot.generation.clone(),
            },
            _ => return Err("The saved origin no longer matches the exact source event.".into()),
        };
        let request = native::Request {
            ivt,
            source: source_input,
            selection,
            conditions,
            policy: native::Policy {
                allowed_end_kinds: allowed,
                require_known_conditions: snapshot.picker.require_conditions,
                require_resolved_ivt_checks: snapshot.picker.require_resolved_ivt,
                require_complete_origin: snapshot.picker.require_complete_origin,
            },
            limits: native::Limits {
                max_comparisons: integer(&snapshot.picker.comparisons, "shared comparison limit")?,
                max_sites_per_enzyme: integer(&snapshot.picker.sites, "sites per enzyme")?,
                max_report_bytes: native::MAX_REPORT_BYTES,
            },
        };
        let revision = native::revision(source, &request).map_err(|e| e.to_string())?;
        Ok(Self {
            snapshot,
            request,
            revision,
        })
    }
    pub fn guard(&self, current: &Snapshot) -> Result<(), String> {
        if self.snapshot != *current {
            return Err(
                "The source event, history or IVT/restriction settings changed; inspect again."
                    .into(),
            );
        }
        Ok(())
    }
    pub fn validate_report(&self, report: &native::Report) -> Result<(), String> {
        report
            .validate_binding(&self.revision)
            .map_err(|e| e.to_string())?;
        if report.source_revision != self.snapshot.binding.revision {
            return Err("The restriction alternatives belong to another source revision.".into());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn snapshot() -> Snapshot {
        Snapshot {
            binding: EventBinding {
                document_id: "document".into(),
                event_id: format!("event-{}", "1".repeat(64)),
                revision: "2".repeat(64),
            },
            epoch: 1,
            ivt: IvtDraft::default(),
            picker: Draft::default(),
            generation: GenerationAssumptions::default(),
            origin_revision: None,
            conditions: ConditionsDraft::default(),
        }
    }
    #[test]
    fn pending_acceptance_binds_event_epoch_cap_origin_generation_and_conditions() {
        let source = Molecule::dna(
            "IVT guard",
            "TAATACGACTCACTATAAGTTCCCGGGTT",
            collagene_model::Topology::Linear,
        )
        .unwrap();
        let mut original = snapshot();
        original.binding.revision = source.revision();
        original.ivt.first_base = "18".into();
        original.ivt.length = "7".into();
        original.picker.enzymes = vec!["smai".into()];
        let capture = Capture::prepare(original.clone(), &source, None).unwrap();
        assert!(capture.guard(&original).is_ok());
        let mut changes = Vec::new();
        let mut changed = original.clone();
        changed.epoch += 1;
        changes.push(changed);
        let mut changed = original.clone();
        changed.binding.event_id = format!("event-{}", "3".repeat(64));
        changes.push(changed);
        let mut changed = original.clone();
        changed.ivt.cap = ivt::CapChoice::CleanCapAgN7113;
        changes.push(changed);
        let mut changed = original.clone();
        changed.origin_revision = Some("4".repeat(64));
        changes.push(changed);
        let mut changed = original.clone();
        changed.generation.dam =
            collagene_core::chemistry::methylation::MethylationState::Unmodified;
        changes.push(changed);
        let mut changed = original;
        changed.conditions.buffer = "rCutSmart".into();
        changes.push(changed);
        for changed in changes {
            assert!(capture.guard(&changed).is_err());
        }
    }
    #[test]
    fn limits_require_literal_bounded_integer_entry() {
        for value in ["", "-1", "1.5", " 1", "999999999999999999999"] {
            assert!(integer(value, "work").is_err());
        }
        assert_eq!(integer("100000000", "work").unwrap(), 100_000_000);
    }
}
