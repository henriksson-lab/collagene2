//! Editable fusion-position search settings. Parsing declares a bounded request;
//! it neither searches DNA nor selects a scientifically preferred junction.

pub const MAX_JUNCTIONS: usize = 12;
const MAX_FIELD_BYTES: usize = 20;
#[derive(Clone, PartialEq, Eq)]
pub struct WindowForm {
    pub minimum: String,
    pub maximum: String,
    /// Empty means enumerate the signed interval; otherwise pin one start.
    pub fixed: String,
}
impl Default for WindowForm {
    fn default() -> Self {
        Self {
            minimum: "-12".into(),
            maximum: "12".into(),
            fixed: String::new(),
        }
    }
}
#[derive(Clone, PartialEq, Eq)]
pub struct SearchForm {
    pub windows: Vec<WindowForm>,
    pub max_candidates: String,
    pub max_states: String,
    pub reject_palindromes: bool,
    pub distinct_reverse_complements: bool,
    pub minimum_hamming: String,
}
#[derive(Clone, PartialEq, Eq)]
pub struct Window {
    pub minimum: i32,
    pub maximum: i32,
    pub fixed: Option<i32>,
}
pub struct Parsed {
    pub windows: Vec<Window>,
    pub max_candidates: usize,
    pub max_states: usize,
}
impl SearchForm {
    pub fn request(
        &self,
        base: collagene_core::golden_gate_junction_edit::Request,
    ) -> Result<collagene_core::golden_gate_fusion_design::Request, String> {
        use collagene_core::golden_gate_fusion_design as fusion;
        let count = base.base.fragments.len().saturating_sub(1) + usize::from(base.base.circular);
        let parsed = self.parse(count)?;
        let minimum_hamming = self
            .minimum_hamming
            .parse::<usize>()
            .map_err(|_| "Enter a minimum Hamming distance.")?;
        Ok(fusion::Request {
            base,
            windows: parsed
                .windows
                .into_iter()
                .enumerate()
                .map(|(junction, window)| fusion::Window {
                    junction,
                    min_shift: window.fixed.unwrap_or(window.minimum),
                    max_shift: window.fixed.unwrap_or(window.maximum),
                })
                .collect(),
            constraints: fusion::Constraints {
                reject_palindromes: self.reject_palindromes,
                distinct_reverse_complements: self.distinct_reverse_complements,
                minimum_hamming,
            },
            limits: fusion::Limits {
                max_states: parsed.max_states,
                max_candidates: parsed.max_candidates,
            },
        })
    }
    pub fn new(junctions: usize) -> Self {
        Self {
            windows: vec![WindowForm::default(); junctions.min(MAX_JUNCTIONS)],
            max_candidates: "64".into(),
            max_states: "10000".into(),
            reject_palindromes: true,
            distinct_reverse_complements: true,
            minimum_hamming: "1".into(),
        }
    }
    pub fn parse(&self, junctions: usize) -> Result<Parsed, String> {
        if junctions == 0 || junctions > MAX_JUNCTIONS || self.windows.len() != junctions {
            return Err("Choose 1–12 declared junctions for fusion-position search.".into());
        }
        if self
            .windows
            .iter()
            .flat_map(|window| [&window.minimum, &window.maximum, &window.fixed])
            .chain([
                &self.max_candidates,
                &self.max_states,
                &self.minimum_hamming,
            ])
            .any(|field| field.len() > MAX_FIELD_BYTES || field.chars().any(char::is_control))
        {
            return Err(
                "Fusion search fields allow at most20 bytes and no control characters.".into(),
            );
        }
        let windows = self
            .windows
            .iter()
            .enumerate()
            .map(|(index, window)| {
                let minimum = window.minimum.parse::<i32>().map_err(|_| {
                    format!("Enter a signed minimum offset for junction {}.", index + 1)
                })?;
                let maximum = window.maximum.parse::<i32>().map_err(|_| {
                    format!("Enter a signed maximum offset for junction {}.", index + 1)
                })?;
                if minimum > maximum {
                    return Err(format!(
                        "Junction {} minimum exceeds its maximum.",
                        index + 1
                    ));
                }
                let fixed = if window.fixed.is_empty() {
                    None
                } else {
                    Some(window.fixed.parse::<i32>().map_err(|_| {
                        format!("Enter a signed fixed offset for junction {}.", index + 1)
                    })?)
                };
                if fixed.is_some_and(|value| value < minimum || value > maximum) {
                    return Err(format!(
                        "Junction {} fixed offset lies outside its declared window.",
                        index + 1
                    ));
                }
                Ok(Window {
                    minimum,
                    maximum,
                    fixed,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        let max_candidates = self
            .max_candidates
            .parse::<usize>()
            .map_err(|_| "Enter a positive candidate limit.")?;
        let max_states = self
            .max_states
            .parse::<usize>()
            .map_err(|_| "Enter a positive search-state limit.")?;
        if max_candidates == 0 || max_states == 0 {
            return Err("Fusion search limits must be positive.".into());
        }
        // Native preflight defines operation-specific maximum spans and budgets.
        Ok(Parsed {
            windows,
            max_candidates,
            max_states,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        assembly_draft_model::Draft,
        assembly_role_ids::Role,
        cloning_methods_model::{Mode, OverlapPreparation, Treatment},
    };
    use collagene_core::{
        assembly_intent,
        assembly_junction_edit::Edit,
        assembly_recipe::Preparation,
        catalog_digest,
        chemistry::{EndRepairMethod, EndTreatment},
        golden_gate_junction_edit::Row,
    };

    fn base_fixture() -> (Draft, assembly_intent::Request, Vec<Edit>) {
        let mut draft = Draft::default();
        draft.circular = false;
        draft.methods.mode = Mode::GoldenGate;
        draft.slots[1].reverse = true;
        draft.sync_rows();
        draft.methods.source_roles = Some(vec![
            Role {
                input_id: "source-capture-a".into(),
                role: "fragment-bin-a".into(),
            },
            Role {
                input_id: "source-capture-b".into(),
                role: "fragment-bin-b".into(),
            },
        ]);
        draft.methods.golden_gate.forward_length = "20".into();
        draft.methods.golden_gate.reverse_length = "21".into();
        for (index, row) in draft.methods.rows.iter_mut().enumerate() {
            row.preparation = OverlapPreparation::ExactPcr;
            row.forward_name = format!("Manual forward {index}");
            row.reverse_name = format!("Manual reverse {index}");
            row.forward = "ACGTACGTACGT".into();
            row.reverse_primer = "TGCATGCATGCA".into();
            row.forward_tail = "GGTCTCAAAA".into();
            row.reverse_tail = "GGTCTCCCCC".into();
        }
        draft.methods.rows[1].right_fusion = "CGTA".into();
        let base = assembly_intent::Request {
            name: "Independent original intent".into(),
            circular: false,
            fragments: draft
                .methods
                .source_roles
                .as_ref()
                .unwrap()
                .iter()
                .enumerate()
                .map(
                    |(index, role)| collagene_core::overlap_primer_design::Fragment {
                        role: role.role.clone(),
                        input_id: role.input_id.clone(),
                        start: 0,
                        end: 80,
                        reverse: index == 1,
                    },
                )
                .collect(),
        };
        let edits = vec![Edit {
            junction: 0,
            remove_left: 1,
            remove_right: 2,
            inserted: "GCA".into(),
        }];
        (draft, base, edits)
    }
    #[test]
    fn signed_windows_and_manual_pins_remain_literal_and_do_not_choose_defaults() {
        let mut form = SearchForm::new(2);
        form.windows[0] = WindowForm {
            minimum: "-7".into(),
            maximum: "3".into(),
            fixed: "-2".into(),
        };
        let parsed = form.parse(2).unwrap();
        assert_eq!(parsed.windows[0].minimum, -7);
        assert_eq!(parsed.windows[0].fixed, Some(-2));
        assert!(parsed.windows[1].fixed.is_none());
        assert_eq!(form.windows[0].fixed, "-2");
        form.windows[0].fixed = "4".into();
        assert!(form.parse(2).is_err());
        form.windows[0].fixed = "".into();
        form.windows[0].minimum = "pending".into();
        assert!(form.parse(2).is_err());
        assert_eq!(form.windows[0].minimum, "pending");
    }
    #[test]
    fn malformed_dimensions_overflow_and_raw_field_budgets_reject_before_request() {
        let mut form = SearchForm::new(1);
        assert!(form.parse(0).is_err());
        assert!(form.parse(2).is_err());
        form.windows[0].minimum = "-99999999999999999999".into();
        assert!(form.parse(1).is_err());
        form.windows[0] = WindowForm::default();
        form.max_states = "1\n".into();
        assert!(form.parse(1).is_err());
        form.max_states = "0".into();
        assert!(form.parse(1).is_err());
        form.max_states = "50000".into();
        form.windows[0].fixed = "1".repeat(21);
        assert!(form.parse(1).is_err());
    }

    #[test]
    fn fixed_manual_pcr_digest_preserves_names_tails_pins_repairs_and_treatments() {
        let (mut draft, base, edits) = base_fixture();
        let pin = catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &["bsai".into()],
        )
        .unwrap();
        let row = &mut draft.methods.rows[0];
        row.preparation = OverlapPreparation::PcrCatalogDigest;
        row.catalog_selection = Some(pin.clone());
        row.fragment = "3".into();
        row.repair = Some(EndRepairMethod::Blunt);
        row.treatment = Treatment::Dephosphorylate;
        let original = draft.clone();
        let request = base_request(&draft, base, edits.clone(), &[false, true]).unwrap();
        let Row::Retain(retained) = &request.rows[0] else {
            panic!("manual row was regenerated");
        };
        let Preparation::PcrCatalogDigest {
            forward,
            reverse,
            selection,
            fragment_index,
            repair,
        } = &retained.source.preparation
        else {
            panic!("manual preparation changed");
        };
        assert_eq!(forward.name, "Manual forward 0");
        assert_eq!(reverse.name, "Manual reverse 0");
        assert_eq!(forward.tail, "GGTCTCAAAA");
        assert_eq!(reverse.tail, "GGTCTCCCCC");
        assert_eq!(*selection, pin);
        assert_eq!(*fragment_index, 2);
        assert_eq!(*repair, Some(EndRepairMethod::Blunt));
        assert_eq!(retained.treatment, Some(EndTreatment::Dephosphorylate));
        assert_eq!(request.edits, edits);
        assert!(draft == original);
        // Native materialization can reject incompatible repair/chemistry; the
        // conversion must preserve the request instead of silently repairing it.
    }

    #[test]
    fn only_explicitly_selected_rows_bypass_stale_invalid_primer_text() {
        let (mut draft, base, edits) = base_fixture();
        draft.methods.rows[1].forward = "NOT VALID DNA".into();
        draft.methods.rows[1].reverse_primer = "also unfinished".into();
        draft.methods.rows[1].fragment = "unparsed".into();
        let original = draft.clone();
        let request = base_request(&draft, base.clone(), edits.clone(), &[false, true]).unwrap();
        assert!(
            matches!(&request.rows[1],Row::DesignPcr{role,input_id,reverse:true,forward_length:20,reverse_length:21}
            if role=="fragment-bin-b"&&input_id=="source-capture-b")
        );
        assert!(base_request(&draft, base, edits, &[false, false]).is_err());
        assert!(draft == original);
    }

    #[test]
    fn linear_terminal_fusion_and_stable_role_bindings_are_never_filled_from_target() {
        let (mut draft, base, edits) = base_fixture();
        let original = draft.clone();
        let request = base_request(&draft, base.clone(), edits.clone(), &[true, true]).unwrap();
        assert_eq!(request.terminal_right_overhang.as_deref(), Some("CGTA"));
        let source_a = collagene_core::Molecule::dna(
            "A",
            &"ACGT".repeat(20),
            collagene_core::Topology::Linear,
        )
        .unwrap();
        let source_b = collagene_core::Molecule::dna(
            "B",
            &"TGCA".repeat(20),
            collagene_core::Topology::Linear,
        )
        .unwrap();
        let originals = [
            ("source-capture-a", &source_a),
            ("source-capture-b", &source_b),
        ];
        assert!(collagene_core::golden_gate_junction_edit::preflight(&originals, &request).is_ok());
        assert_eq!(request.base.fragments[0].input_id, "source-capture-a");
        assert_eq!(request.base.fragments[1].role, "fragment-bin-b");
        assert!(
            matches!(&request.rows[0],Row::DesignPcr{input_id,..} if input_id=="source-capture-a")
        );
        assert!(draft == original);
        draft.methods.rows[1].right_fusion.clear();
        let empty = base_request(&draft, base.clone(), edits.clone(), &[true, true]).unwrap();
        assert_eq!(empty.terminal_right_overhang.as_deref(), Some(""));
        assert!(collagene_core::golden_gate_junction_edit::preflight(&originals, &empty).is_err());
        let retained_last = base_request(&draft, base, edits, &[true, false]).unwrap();
        assert!(retained_last.terminal_right_overhang.is_none());
    }
}

/// Shared raw-draft conversion for ordinary junction design and fusion placement.
/// Retained manual preparations remain literal; selected rows request native PCR.
pub fn base_request(
    before: &crate::assembly_draft_model::Draft,
    base: collagene_core::assembly_intent::Request,
    edits: Vec<collagene_core::assembly_junction_edit::Edit>,
    regenerate: &[bool],
) -> Result<collagene_core::golden_gate_junction_edit::Request, String> {
    use crate::cloning_methods_model::{Mode, OverlapInput};
    use collagene_core::golden_gate_junction_edit::{Request, Row};
    if before.methods.mode != Mode::GoldenGate
        || regenerate.len() != before.slots.len()
        || before.methods.rows.len() != before.slots.len()
        || before.slots.is_empty()
        || before.slots.len() > MAX_JUNCTIONS
    {
        return Err(
            "Choose Golden Gate and an explicit preparation choice for every source.".into(),
        );
    }
    let mut retained = before.methods.clone();
    for (index, selected) in regenerate.iter().enumerate() {
        if *selected {
            retained.rows[index] = OverlapInput {
                document: index,
                reverse: before.slots[index].reverse,
                ..Default::default()
            };
        }
    }
    let (_, fragments) = crate::cloning_preparations::prepared_rows(&retained, before.slots.len())?;
    let gate = &before.methods.golden_gate;
    let (forward_length, reverse_length) = if regenerate.iter().any(|selected| *selected) {
        (
            gate.forward_length
                .parse()
                .map_err(|_| "Enter a Golden Gate forward annealing length.")?,
            gate.reverse_length
                .parse()
                .map_err(|_| "Enter a Golden Gate reverse annealing length.")?,
        )
    } else {
        (0, 0)
    };
    let rows = fragments
        .into_iter()
        .enumerate()
        .map(|(index, fragment)| {
            if regenerate[index] {
                Row::DesignPcr {
                    role: fragment.role,
                    input_id: fragment.input_id,
                    reverse: fragment.reverse,
                    forward_length,
                    reverse_length,
                }
            } else {
                Row::Retain(collagene_core::ligation_recipe::Fragment {
                    source: fragment,
                    treatment: before.methods.rows[index].treatment.operation(),
                })
            }
        })
        .collect();
    let terminal_right_overhang = if !before.circular && regenerate.last() == Some(&true) {
        Some(
            before
                .methods
                .rows
                .last()
                .ok_or("Missing final source preparation.")?
                .right_fusion
                .clone(),
        )
    } else {
        None
    };
    Ok(Request {
        base,
        edits,
        rows,
        enzyme: gate.enzyme,
        clamp: gate.clamp.clone(),
        spacer: gate.spacer.clone(),
        terminal_right_overhang,
        allow_nicked: before.methods.allow_nicked,
    })
}
