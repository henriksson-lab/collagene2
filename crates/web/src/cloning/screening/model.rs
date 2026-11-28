use crate::{history::run_origin_model::PinnedOrigin, sources::identity::SourceInfo};
use collagene_core::{
    chemistry::{Strand, tracked_duplex::GenerationAssumptions},
    cloning::{restriction::catalog, screening as native},
    documents::workspace::Workspace,
    material::{
        derivation::SourceOrigin,
        lineage::{self, EventBinding},
    },
    sequence::analysis::Primer,
};
use collagene_model::{Alphabet, Molecule};

pub(super) fn enzyme_choices() -> Result<Vec<(String, String)>, String> {
    let catalogue = catalog::builtin_catalog().map_err(|error| error.to_string())?;
    Ok(catalogue
        .edition()
        .enzymes
        .iter()
        .map(|enzyme| (enzyme.id.clone(), enzyme.name.clone()))
        .collect())
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ComponentDraft {
    pub id: String,
    pub label: String,
    pub binding: EventBinding,
    pub origin: Option<PinnedOrigin>,
    pub declare_generation: bool,
    pub generation: GenerationAssumptions,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct SubjectDraft {
    pub id: String,
    pub label: String,
    pub role: native::Role,
    pub components: Vec<String>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct PcrDraft {
    pub name: String,
    pub forward: String,
    pub forward_tail: String,
    pub reverse: String,
    pub reverse_tail: String,
    pub bottom: bool,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum Candidate {
    Digest { label: String, enzymes: Vec<String> },
    Pcr(PcrDraft),
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Draft {
    pub components: Vec<ComponentDraft>,
    pub subjects: Vec<SubjectDraft>,
    pub candidates: Vec<Candidate>,
    pub tolerance: String,
    pub cleavage_work: String,
    pub pcr_work: String,
    pub panel_work: String,
    pub next_id: u32,
}
impl Default for Draft {
    fn default() -> Self {
        Self {
            components: vec![],
            subjects: vec![],
            candidates: vec![],
            tolerance: "0".into(),
            cleavage_work: "100000000".into(),
            pcr_work: "200000000".into(),
            panel_work: "100000".into(),
            next_id: 0,
        }
    }
}
impl Draft {
    pub fn add_component(&mut self, info: &SourceInfo) -> Result<(), String> {
        if self.components.len() >= native::MAX_COMPONENTS {
            return Err("At most 16 retained components are supported.".into());
        }
        self.next_id = self
            .next_id
            .checked_add(1)
            .ok_or("Component identifier limit reached.")?;
        let id = format!("component-{}", self.next_id);
        self.components.push(ComponentDraft {
            id: id.clone(),
            label: info.name.chars().take(128).collect(),
            binding: info.binding.clone(),
            origin: None,
            declare_generation: false,
            generation: GenerationAssumptions::default(),
        });
        if self.subjects.is_empty() {
            self.subjects.push(SubjectDraft {
                id: "intended".into(),
                label: "Intended product".into(),
                role: native::Role::Intended,
                components: vec![id],
            });
        }
        Ok(())
    }
    pub fn add_decoy(&mut self) -> Result<(), String> {
        if self.subjects.len() >= native::MAX_SUBJECTS {
            return Err("At most 16 named outcomes are supported.".into());
        }
        self.next_id = self
            .next_id
            .checked_add(1)
            .ok_or("Outcome identifier limit reached.")?;
        self.subjects.push(SubjectDraft {
            id: format!("decoy-{}", self.next_id),
            label: "Named alternative".into(),
            role: native::Role::Decoy(native::DecoyKind::Other),
            components: vec![],
        });
        Ok(())
    }
    pub fn remove_component(&mut self, id: &str) {
        self.components.retain(|c| c.id != id);
        for subject in &mut self.subjects {
            subject.components.retain(|component| component != id);
        }
    }
    pub fn add_catalogue_single_digests(&mut self) -> Result<usize, String> {
        let additions: Vec<_> = enzyme_choices()?
            .into_iter()
            .filter(|(id, _)| {
                !self.candidates.iter().any(|candidate| {
                    matches!(candidate, Candidate::Digest { enzymes, .. }
                    if enzymes.len() == 1 && enzymes[0] == *id)
                })
            })
            .collect();
        if additions.len() > native::MAX_TESTS.saturating_sub(self.candidates.len()) {
            return Err("The remaining catalogue digests do not fit the 64-test limit. Remove candidates or select individual enzymes.".into());
        }
        let count = additions.len();
        self.candidates
            .extend(additions.into_iter().map(|(id, name)| Candidate::Digest {
                label: name,
                enzymes: vec![id],
            }));
        Ok(count)
    }
}
pub(super) fn integer(value: &str, label: &str) -> Result<usize, String> {
    if value.is_empty() || value.len() > 20 || !value.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("Enter a whole number for {label}."));
    }
    value
        .parse()
        .map_err(|_| format!("{label} exceeds the supported range."))
}
pub(super) fn molecule<'a>(
    workspace: &'a Workspace,
    binding: &EventBinding,
) -> Result<&'a Molecule, String> {
    let event = workspace
        .history
        .iter()
        .find(|event| event.id == binding.event_id)
        .ok_or("The retained source event is unavailable.")?;
    if event.document_id != binding.document_id || event.revision != binding.revision {
        return Err("The source event belongs to a different document or revision.".into());
    }
    workspace
        .revisions
        .get(&binding.revision)
        .ok_or_else(|| "The source snapshot is unavailable.".into())
}
pub(super) fn info(workspace: &Workspace, event_id: &str) -> Option<SourceInfo> {
    let event = workspace
        .history
        .iter()
        .find(|event| event.id == event_id)?;
    let molecule = workspace.revisions.get(&event.revision)?;
    if molecule.sequence.alphabet() != Alphabet::Dna {
        return None;
    }
    Some(SourceInfo {
        binding: EventBinding {
            document_id: event.document_id.clone(),
            event_id: event.id.clone(),
            revision: event.revision.clone(),
        },
        name: molecule.name.clone(),
        length: molecule.sequence.len(),
        alphabet: molecule.sequence.alphabet(),
        topology: molecule.topology,
    })
}
/// Bounded visible choices; exact selected events survive filters and duplicate names.
pub(super) fn event_choices(workspace: &Workspace, query: &str, selected: &str) -> Vec<SourceInfo> {
    let query = query.to_lowercase();
    let mut choices: Vec<_> = info(workspace, selected).into_iter().collect();
    let remaining = 100 - choices.len();
    choices.extend(
        workspace
            .history
            .iter()
            .rev()
            .filter(|event| event.id != selected)
            .filter_map(|event| {
                let molecule = workspace.revisions.get(&event.revision)?;
                if molecule.sequence.alphabet() != Alphabet::Dna {
                    return None;
                }
                let source = SourceInfo {
                    binding: EventBinding {
                        document_id: event.document_id.clone(),
                        event_id: event.id.clone(),
                        revision: event.revision.clone(),
                    },
                    name: molecule.name.clone(),
                    length: molecule.sequence.len(),
                    alphabet: molecule.sequence.alphabet(),
                    topology: molecule.topology,
                };
                let searchable = format!(
                    "{} {} {}",
                    source.name, source.binding.document_id, source.binding.event_id
                )
                .to_lowercase();
                (query.is_empty() || searchable.contains(&query)).then_some(source)
            })
            .take(remaining),
    );
    choices
}
fn preclone(value: &impl serde::Serialize, remaining: &mut usize) -> Result<(), String> {
    struct Remaining<'a>(&'a mut usize);
    impl std::io::Write for Remaining<'_> {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            *self.0 = self.0.checked_sub(bytes.len()).ok_or_else(|| {
                std::io::Error::other("Screening sources exceed the request bound")
            })?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(Remaining(remaining), value).map_err(|e| e.to_string())
}
#[derive(Clone)]
pub(super) struct Capture {
    pub draft: Draft,
    pub source: EventBinding,
    pub epoch: u64,
    pub fingerprint: String,
    pub request: native::Request,
    pub revision: String,
}
impl Capture {
    pub fn prepare(
        workspace: &Workspace,
        source: &SourceInfo,
        epoch: u64,
        draft: Draft,
    ) -> Result<Self, String> {
        if draft.components.len() > native::MAX_COMPONENTS
            || draft.subjects.len() > native::MAX_SUBJECTS
            || draft.candidates.len() > native::MAX_TESTS
        {
            return Err("Screening setup exceeds its supported row limits.".into());
        }
        let mut remaining = native::MAX_REQUEST_BYTES;
        let mut bases = 0usize;
        // Admit all borrowed original assets before cloning any of them.
        for component in &draft.components {
            let molecule = molecule(workspace, &component.binding)?;
            bases = bases
                .checked_add(molecule.sequence.len())
                .ok_or("Source length limit exceeded.")?;
            if bases > 1_000_000 {
                return Err("Screening sources exceed one million combined bases.".into());
            }
            preclone(molecule, &mut remaining)?;
            if let Some(origin) = &component.origin {
                if origin.event() != &component.binding {
                    return Err("Saved material belongs to another exact source event.".into());
                }
                preclone(origin.as_ref(), &mut remaining)?;
                origin
                    .as_ref()
                    .validate_for_workspace(workspace, &component.binding)
                    .map_err(|e| e.to_string())?;
                if component.declare_generation {
                    return Err(
                        "Saved material cannot be re-declared as original source generation."
                            .into(),
                    );
                }
            }
        }
        let catalog = catalog::builtin_catalog().map_err(|e| e.to_string())?;
        let tests = draft
            .candidates
            .iter()
            .enumerate()
            .map(|(index, candidate)| {
                let (label, method) = match candidate {
                    Candidate::Digest { label, enzymes } => (
                        label.clone(),
                        native::Method::Digest {
                            selection: catalog::select(&catalog, enzymes)
                                .map_err(|e| e.to_string())?,
                        },
                    ),
                    Candidate::Pcr(pcr) => (
                        pcr.name.clone(),
                        native::Method::ExactPcr {
                            forward: Primer::new("Forward", &pcr.forward, &pcr.forward_tail)
                                .map_err(|e| e.to_string())?,
                            reverse: Primer::new("Reverse", &pcr.reverse, &pcr.reverse_tail)
                                .map_err(|e| e.to_string())?,
                            template_strand: if pcr.bottom {
                                Strand::Bottom
                            } else {
                                Strand::Top
                            },
                        },
                    ),
                };
                Ok(native::Test {
                    id: format!("test-{index}"),
                    label,
                    method,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        let request = native::Request {
            components: draft
                .components
                .iter()
                .map(|component| {
                    Ok(native::Component {
                        id: component.id.clone(),
                        label: component.label.clone(),
                        revision: component.binding.revision.clone(),
                        molecule: molecule(workspace, &component.binding)?.clone(),
                        origin: match &component.origin {
                            Some(origin) => SourceOrigin::Saved(origin.as_ref().clone()),
                            None => SourceOrigin::OriginalReference,
                        },
                        original_generation: component
                            .declare_generation
                            .then(|| component.generation.clone()),
                    })
                })
                .collect::<Result<Vec<_>, String>>()?,
            subjects: draft
                .subjects
                .iter()
                .map(|s| native::Subject {
                    id: s.id.clone(),
                    label: s.label.clone(),
                    role: s.role.clone(),
                    components: s.components.clone(),
                })
                .collect(),
            tests,
            readout: native::IdealSizeReadout {
                counterpart_tolerance_twice_bp: integer(
                    &draft.tolerance,
                    "size tolerance (half-bp units)",
                )?,
            },
            limits: native::Limits {
                max_cleavage_work: integer(&draft.cleavage_work, "shared cleavage work")?,
                max_pcr_work: integer(&draft.pcr_work, "shared PCR work")?,
                max_panel_checks: integer(&draft.panel_work, "panel comparison work")?,
                ..native::Limits::default()
            },
        };
        let revision = native::revision(&request).map_err(|e| e.to_string())?;
        Ok(Self {
            fingerprint: lineage::workspace_fingerprint(workspace).map_err(|e| e.to_string())?,
            source: source.binding.clone(),
            epoch,
            draft,
            request,
            revision,
        })
    }
    pub fn guard(
        &self,
        workspace: &Workspace,
        source: &SourceInfo,
        epoch: u64,
        draft: &Draft,
    ) -> Result<(), String> {
        if self.source != source.binding
            || self.epoch != epoch
            || self.draft != *draft
            || self.fingerprint
                != lineage::workspace_fingerprint(workspace).map_err(|e| e.to_string())?
        {
            return Err(
                "Sources, saved origins or screening settings changed; run the analysis again."
                    .into(),
            );
        }
        for component in &draft.components {
            molecule(workspace, &component.binding)?;
        }
        Ok(())
    }
    pub fn validate_report(&self, report: &native::Report) -> Result<(), String> {
        report
            .validate_binding(&self.revision)
            .map_err(|e| e.to_string())
    }
    pub fn export_json(&self, report: &native::Report) -> Result<Vec<u8>, String> {
        self.validate_report(report)?;
        #[derive(serde::Serialize)]
        struct ExportEvent<'a> {
            component_id: &'a str,
            event: &'a EventBinding,
        }
        #[derive(serde::Serialize)]
        struct Export<'a> {
            format: &'static str,
            workspace_fingerprint: &'a str,
            active_source: &'a EventBinding,
            component_events: Vec<ExportEvent<'a>>,
            report: &'a native::Report,
        }
        let value = Export {
            format: "collagene-screening-ui-evidence-1",
            workspace_fingerprint: &self.fingerprint,
            active_source: &self.source,
            component_events: self
                .draft
                .components
                .iter()
                .map(|c| ExportEvent {
                    component_id: c.id.as_str(),
                    event: &c.binding,
                })
                .collect(),
            report,
        };
        let mut remaining = native::MAX_REPORT_BYTES;
        preclone(&value, &mut remaining)?;
        serde_json::to_vec(&value).map_err(|error| error.to_string())
    }
    /// Only genuine retained native recognition coordinates on the active reference
    /// are navigated. A different source must be opened by the document owner first.
    pub fn recognition_ranges(
        &self,
        report: &native::Report,
        cell_index: usize,
        site_index: usize,
        active: &SourceInfo,
    ) -> Result<Vec<(usize, usize)>, String> {
        self.validate_report(report)?;
        let cell = report
            .cells
            .get(cell_index)
            .ok_or("Screening result is unavailable.")?;
        let component = self
            .draft
            .components
            .get(cell.component_index)
            .ok_or("Screening source is unavailable.")?;
        if component.binding != active.binding {
            return Err(
                "Open this exact source reference before selecting its recognition site.".into(),
            );
        }
        if component.origin.is_some() {
            return Err("Saved-material recognition coordinates are retained as text; no native reference projection was supplied.".into());
        }
        let evidence = cell
            .cleavage
            .as_ref()
            .ok_or("No native recognition evidence is available.")?;
        let site = evidence
            .sites
            .get(site_index)
            .ok_or("Recognition site is unavailable.")?;
        let enzyme = evidence
            .options
            .selection
            .enzymes
            .iter()
            .find(|e| e.enzyme.name == site.enzyme)
            .ok_or("Recognition identity is unavailable.")?;
        let start = site.recognition_start;
        let end = start
            .checked_add(enzyme.enzyme.recognition.len())
            .ok_or("Recognition coordinates exceed the supported range.")?;
        if start >= active.length || end > active.length.saturating_mul(2) {
            return Err("Recognition coordinates are outside this reference.".into());
        }
        if end <= active.length {
            Ok(vec![(start, end)])
        } else if active.topology == collagene_model::Topology::Circular {
            Ok(vec![(start, active.length), (0, end - active.length)])
        } else {
            Err("Recognition site extends outside the linear reference.".into())
        }
    }
}
