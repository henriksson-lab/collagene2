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
#[derive(Clone, PartialEq, Eq)]
pub struct RouteForm {
    pub start: String,
    pub end: String,
    pub forward: String,
    pub reverse: String,
    pub reverse_tail: String,
}
impl Default for RouteForm {
    fn default() -> Self {
        Self {
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
        } = built.route;
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
}
