//! Parse editable declarations without changing source DNA or choosing a result.
use collagene_core::{
    Molecule, Topology, catalog_digest, codon_design, domestication as design,
    domestication_preparation as physical,
};
#[derive(Clone, PartialEq, Eq)]
pub struct Form {
    pub enzymes: String,
    pub scope_start: String,
    pub scope_end: String,
    pub protected: String,
    pub features: Vec<usize>,
    pub retained: Vec<design::SitePin>,
    pub beam: String,
    pub work: String,
    pub max_patches: String,
}
impl Form {
    pub fn new(source: &Molecule, feature: Option<usize>) -> Self {
        Self {
            enzymes: "bsai".into(),
            scope_start: "1".into(),
            scope_end: source.sequence.len().to_string(),
            protected: String::new(),
            features: feature
                .filter(|&i| {
                    source
                        .features
                        .get(i)
                        .is_some_and(|f| f.kind.eq_ignore_ascii_case("CDS"))
                })
                .into_iter()
                .collect(),
            retained: vec![],
            beam: "128".into(),
            work: "20000000".into(),
            max_patches: "64".into(),
        }
    }
    pub fn request(&self, source: &Molecule) -> Result<design::Request, String> {
        if self.enzymes.len() > 2048
            || self.protected.len() > 8192
            || self.features.len() > design::MAX_FEATURES
            || self.retained.len() > design::MAX_SITES
            || [
                &self.scope_start,
                &self.scope_end,
                &self.beam,
                &self.work,
                &self.max_patches,
            ]
            .iter()
            .any(|s| s.len() > 20 || s.chars().any(char::is_control))
        {
            return Err(
                "Domestication declarations exceed their field or selection bounds.".into(),
            );
        }
        let ids = self
            .enzymes
            .split(|c: char| c == ',' || c.is_ascii_whitespace())
            .filter(|s| !s.is_empty())
            .map(str::to_ascii_lowercase)
            .collect::<Vec<_>>();
        if ids.is_empty() || ids.len() > 16 {
            return Err("Choose 1–16 enzyme catalog IDs.".into());
        }
        let selection = catalog_digest::select(
            &catalog_digest::builtin_catalog().map_err(|e| e.to_string())?,
            &ids,
        )
        .map_err(|e| e.to_string())?;
        let scope = range(source, &self.scope_start, &self.scope_end)?;
        let mut protected = Vec::new();
        for text in self
            .protected
            .split(|c: char| c == ',' || c.is_ascii_whitespace())
            .filter(|s| !s.is_empty())
        {
            if protected.len() >= 256 {
                return Err("Use at most256 protected ranges.".into());
            }
            let (start, end) = text
                .split_once('-')
                .ok_or("Protected ranges use start-end, with 1-based inclusive coordinates.")?;
            protected.push(range(source, start, end)?);
        }
        let mut seen = std::collections::BTreeSet::new();
        let features = self
            .features
            .iter()
            .map(|&i| {
                if !seen.insert(i) {
                    return Err("A coding feature was selected twice.".to_string());
                }
                design::feature_pin(source, i).map_err(|e| e.to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let work = self
            .work
            .parse::<u64>()
            .map_err(|_| "Enter a positive work limit.")?;
        let request = design::Request {
            source_revision: source.revision(),
            selection,
            scope: vec![scope],
            retained_sites: self.retained.clone(),
            coding_features: features,
            protected,
            options: codon_design::DesignOptions {
                beam_width: self.beam.parse().map_err(|_| "Enter a beam width.")?,
                max_work: work,
                max_homopolymer: 30000,
                ..Default::default()
            },
            limits: design::Limits {
                max_patches: self
                    .max_patches
                    .parse()
                    .map_err(|_| "Enter a maximum substitution count.")?,
                max_work: work,
            },
        };
        design::preflight(source, &request).map_err(|e| e.to_string())?;
        Ok(request)
    }
}
pub fn range(source: &Molecule, start: &str, end: &str) -> Result<design::SourceRange, String> {
    if start.len() > 20 || end.len() > 20 {
        return Err("Coordinates exceed their field bounds.".into());
    }
    let first = start
        .parse::<usize>()
        .map_err(|_| "Enter a 1-based start coordinate.")?;
    let last = end
        .parse::<usize>()
        .map_err(|_| "Enter a 1-based end coordinate.")?;
    let n = source.sequence.len();
    if first == 0 || first > n || last == 0 || last > n {
        return Err("Coordinates must name existing source bases.".into());
    }
    let last = if last < first {
        if source.topology != Topology::Circular {
            return Err("Only a circular source can cross its origin.".into());
        }
        last.checked_add(n).ok_or("Coordinate overflow")?
    } else {
        last
    };
    Ok(design::SourceRange {
        start: first - 1,
        end: last,
    })
}
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PreparationKind {
    InversePcr,
    OverlapExtension,
}
/// Editable raw values survive route switches, including unfinished numbers.
#[derive(Clone, PartialEq, Eq)]
pub struct OverlapForm {
    pub values: [String; 16],
}
impl Default for OverlapForm {
    fn default() -> Self {
        let o = collagene_core::overlap_extension_primer_design::Options::default();
        let p = o.primer_options;
        Self {
            values: [
                p.min_length.to_string(),
                p.max_length.to_string(),
                p.target_tm_c.to_string(),
                p.max_tm_difference_c.to_string(),
                p.min_gc_fraction.to_string(),
                p.max_gc_fraction.to_string(),
                p.max_homopolymer.to_string(),
                p.max_three_prime_complementarity.to_string(),
                p.max_pairs.to_string(),
                o.homology_flank.to_string(),
                o.max_oligo_length.to_string(),
                o.min_paired_bases.to_string(),
                o.max_comparisons.to_string(),
                o.max_candidates.to_string(),
                p.conditions.sodium_molar.to_string(),
                p.conditions.total_strand_molar.to_string(),
            ],
        }
    }
}
impl OverlapForm {
    pub fn options(
        &self,
    ) -> Result<collagene_core::overlap_extension_primer_design::Options, String> {
        if self
            .values
            .iter()
            .any(|s| s.is_empty() || s.len() > 64 || s.chars().any(char::is_control))
        {
            return Err(
                "Complete bounded overlap-extension settings before preparing primers.".into(),
            );
        }
        let integer = |i: usize| {
            self.values[i]
                .parse::<usize>()
                .map_err(|_| "Enter whole-number overlap-extension bounds.".to_string())
        };
        let real = |i: usize| -> Result<f64, String> {
            let v = self.values[i]
                .parse::<f64>()
                .map_err(|_| "Enter finite primer model settings.".to_string())?;
            if !v.is_finite() {
                return Err("Primer model settings must be finite.".into());
            }
            Ok(v)
        };
        let mut o = collagene_core::overlap_extension_primer_design::Options::default();
        o.primer_options.min_length = integer(0)?;
        o.primer_options.max_length = integer(1)?;
        o.primer_options.target_tm_c = real(2)?;
        o.primer_options.max_tm_difference_c = real(3)?;
        o.primer_options.min_gc_fraction = real(4)?;
        o.primer_options.max_gc_fraction = real(5)?;
        o.primer_options.max_homopolymer = integer(6)?;
        o.primer_options.max_three_prime_complementarity = integer(7)?;
        o.primer_options.max_pairs = integer(8)?;
        o.homology_flank = integer(9)?;
        o.max_oligo_length = integer(10)?;
        o.min_paired_bases = integer(11)?;
        o.max_comparisons = integer(12)?;
        o.max_candidates = integer(13)?;
        o.primer_options.conditions.sodium_molar = real(14)?;
        o.primer_options.conditions.total_strand_molar = real(15)?;
        Ok(o)
    }
}
#[derive(Clone, PartialEq, Eq)]
pub struct RouteForm {
    pub kind: PreparationKind,
    pub overlap_extension: OverlapForm,
    pub start: String,
    pub end: String,
    pub forward: String,
    pub reverse: String,
    pub reverse_tail: String,
}
impl Default for RouteForm {
    fn default() -> Self {
        Self {
            kind: PreparationKind::InversePcr,
            overlap_extension: Default::default(),
            start: String::new(),
            end: String::new(),
            forward: "18".into(),
            reverse: "18".into(),
            reverse_tail: "0".into(),
        }
    }
}
impl RouteForm {
    pub fn request(
        &self,
        source: &Molecule,
        design: design::Request,
        selection: design::Selection,
    ) -> Result<physical::Request, String> {
        if [
            &self.start,
            &self.end,
            &self.forward,
            &self.reverse,
            &self.reverse_tail,
        ]
        .iter()
        .any(|s| s.len() > 20 || s.chars().any(char::is_control))
        {
            return Err("Preparation fields exceed their bounds.".into());
        }
        if self
            .overlap_extension
            .values
            .iter()
            .any(|s| s.len() > 64 || s.chars().any(char::is_control))
        {
            return Err("Overlap-extension raw settings exceed their bounds.".into());
        }
        if self.kind == PreparationKind::OverlapExtension {
            return Ok(physical::Request {
                design,
                selection,
                route: physical::Route::OverlapExtension {
                    options: self.overlap_extension.options()?,
                },
            });
        }
        let interval = range(source, &self.start, &self.end)?;
        Ok(physical::Request {
            design,
            selection,
            route: physical::Route::InversePcr {
                edit_start: interval.start,
                edit_end: interval.end,
                forward_length: self
                    .forward
                    .parse()
                    .map_err(|_| "Enter a forward annealing length.")?,
                reverse_length: self
                    .reverse
                    .parse()
                    .map_err(|_| "Enter a reverse annealing length.")?,
                reverse_tail_bases: self
                    .reverse_tail
                    .parse()
                    .map_err(|_| "Enter a reverse-tail split length.")?,
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{Feature, Location};
    fn source(topology: Topology) -> Molecule {
        let mut source = Molecule::dna(
            "annotated source",
            "ATGGGTCTCTAACCTGACGTAGCTACGATCGTAC",
            topology,
        )
        .unwrap();
        source.features.push(Feature {
            kind: "CDS".into(),
            location: Location::range(0, 12),
            qualifiers: vec![],
            needs_review: false,
        });
        source
    }
    #[test]
    fn circular_scopes_protected_ranges_and_feature_pins_use_original_coordinates() {
        let source = source(Topology::Circular);
        let before = source.clone();
        let mut form = Form::new(&source, Some(0));
        form.scope_start = "28".into();
        form.scope_end = "4".into();
        form.protected = "8-10, 31-2".into();
        form.retained.push(design::SitePin {
            enzyme_id: "bsai".into(),
            start: 3,
            reverse: false,
        });
        let request = form.request(&source).unwrap();
        assert_eq!(request.scope[0].start, 27);
        assert_eq!(request.scope[0].end, source.sequence.len() + 4);
        assert_eq!(request.protected[1].start, 30);
        assert_eq!(request.protected[1].end, source.sequence.len() + 2);
        assert_eq!(
            request.coding_features[0],
            design::feature_pin(&source, 0).unwrap()
        );
        assert_eq!(request.retained_sites, form.retained);
        assert_eq!(source, before);
        let mut changed = source.clone();
        changed.features[0]
            .qualifiers
            .push(("note".into(), Some("new qualifier".into())));
        assert_ne!(
            form.request(&changed).unwrap().coding_features[0],
            request.coding_features[0]
        );
    }
    #[test]
    fn literal_invalid_fields_and_out_of_source_or_linear_wrap_are_refused() {
        let source = source(Topology::Linear);
        let mut form = Form::new(&source, Some(0));
        form.scope_start = "28".into();
        form.scope_end = "4".into();
        assert!(form.request(&source).is_err());
        form.scope_start = "1".into();
        form.scope_end = source.sequence.len().to_string();
        form.work = "12\n".into();
        assert!(form.request(&source).is_err());
        assert_eq!(form.work, "12\n");
        form.work = "20000000".into();
        form.protected = "2:3".into();
        assert!(form.request(&source).is_err());
        form.protected = String::new();
        form.features = vec![0, 0];
        assert!(form.request(&source).is_err());
        assert!(range(&source, "0", "3").is_err());
        assert!(range(&source, "1", "99999999999999999999").is_err());
    }
    #[test]
    fn route_encodes_explicit_wrapped_interval_and_tail_split_without_mutating_intent() {
        let source = source(Topology::Circular);
        let request = Form::new(&source, Some(0)).request(&source).unwrap();
        let selection = design::Selection {
            context_revision: "a".repeat(64),
            candidate_id: "b".repeat(64),
            patches: vec![],
        };
        let route = RouteForm {
            kind: PreparationKind::InversePcr,
            overlap_extension: Default::default(),
            start: "31".into(),
            end: "2".into(),
            forward: "8".into(),
            reverse: "9".into(),
            reverse_tail: "2".into(),
        };
        let before = serde_json::to_vec(&request).unwrap();
        let built = route.request(&source, request.clone(), selection).unwrap();
        let physical::Route::InversePcr {
            edit_start,
            edit_end,
            forward_length,
            reverse_length,
            reverse_tail_bases,
        } = built.route
        else {
            panic!("inverse route was explicitly selected")
        };
        assert_eq!(
            (
                edit_start,
                edit_end,
                forward_length,
                reverse_length,
                reverse_tail_bases
            ),
            (30, source.sequence.len() + 2, 8, 9, 2)
        );
        assert_eq!(serde_json::to_vec(&built.design).unwrap(), before);
        // Geometry/patch coverage remains the native physical engine's job.
        let input = collagene_core::verification::PlanInput {
            id: physical::INPUT_ID.into(),
            revision: source.revision(),
            molecule: source,
        };
        assert!(physical::preflight(&input, &built).is_err());
    }
    #[test]
    fn route_switch_preserves_inactive_partial_text_and_explicit_oe_options() {
        let source = source(Topology::Linear);
        let request = Form::new(&source, Some(0)).request(&source).unwrap();
        let selection = design::Selection {
            context_revision: "a".repeat(64),
            candidate_id: "b".repeat(64),
            patches: vec![],
        };
        let mut form = RouteForm::default();
        form.start = "unfinished".into();
        form.kind = PreparationKind::OverlapExtension;
        form.overlap_extension.values[0] = "8".into();
        form.overlap_extension.values[1] = "12".into();
        let native = form
            .request(&source, request.clone(), selection.clone())
            .unwrap();
        assert!(
            matches!(native.route,physical::Route::OverlapExtension{options} if options.primer_options.min_length==8&&options.primer_options.max_length==12)
        );
        assert_eq!(form.start, "unfinished");
        form.kind = PreparationKind::InversePcr;
        assert!(
            form.request(&source, request.clone(), selection.clone())
                .is_err()
        );
        form.kind = PreparationKind::OverlapExtension;
        form.overlap_extension.values[2] = "6e".into();
        assert!(form.request(&source, request, selection).is_err());
        assert_eq!(form.overlap_extension.values[2], "6e");
    }
}
