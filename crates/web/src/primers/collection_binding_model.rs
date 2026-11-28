//! Bounded collection input and current-result context, independent of DOM state.
use crate::sources::identity::SourceInfo;
use collagene_core::cloning::pcr::primers::ThermoConditions;
use collagene_model::Sequence;
use collagene_oligo_binding::{Budget, BudgetLimits, Collection, Hit, Oligo, Options, Status};

pub const MAX_TEXT_BYTES: usize = 768 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Field {
    Minimum,
    Extra,
    Mismatches,
    Prime,
    Match,
    Mismatch,
    Score,
    Work,
    Hits,
    Output,
    Sodium,
    Strands,
}
impl Field {
    pub const MAIN: [Self; 4] = [Self::Minimum, Self::Extra, Self::Mismatches, Self::Prime];
    pub const SCORES: [Self; 3] = [Self::Match, Self::Mismatch, Self::Score];
    pub const BUDGET: [Self; 3] = [Self::Work, Self::Hits, Self::Output];
    pub const THERMO: [Self; 2] = [Self::Sodium, Self::Strands];
    pub fn label(self) -> &'static str {
        match self {
            Self::Minimum => "Minimum annealed bases",
            Self::Extra => "Maximum additional unpaired 5′ bases",
            Self::Mismatches => "Mismatch limit per concrete realization",
            Self::Prime => "Required possible contiguous 3′ pairing",
            Self::Match => "Match score",
            Self::Mismatch => "Mismatch score",
            Self::Score => "Minimum alignment score",
            Self::Work => "Shared work allowance",
            Self::Hits => "Retained hit allowance",
            Self::Output => "Output allocation charge allowance (bytes)",
            Self::Sodium => "Sodium (mM)",
            Self::Strands => "Total strands (nM)",
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Draft {
    pub text: String,
    values: [String; 12],
}
impl Default for Draft {
    fn default() -> Self {
        Self {
            text: String::new(),
            values: [
                "12", "200", "2", "0", "2", "-3", "12", "20000000", "1000", "8388608", "50", "500",
            ]
            .map(String::from),
        }
    }
}
impl Draft {
    pub fn value(&self, field: Field) -> String {
        self.values[field as usize].clone()
    }
    pub fn set(&mut self, field: Field, value: String) {
        self.values[field as usize] = value;
    }
    fn number<T: std::str::FromStr>(&self, field: Field) -> Result<T, String> {
        self.values[field as usize]
            .trim()
            .parse()
            .map_err(|_| format!("{} is not a valid number.", field.label()))
    }
    pub fn prepare(&self) -> Result<(Collection, Options, BudgetLimits, ThermoConditions), String> {
        let collection = parse_collection(&self.text)?;
        let options = Options {
            minimum_annealing_length: self.number(Field::Minimum)?,
            maximum_extra_unpaired: self.number(Field::Extra)?,
            maximum_mismatches: self.number(Field::Mismatches)?,
            minimum_three_prime_pairing: self.number(Field::Prime)?,
            match_score: self.number(Field::Match)?,
            mismatch_score: self.number(Field::Mismatch)?,
            minimum_score: self.number(Field::Score)?,
        };
        options.validate().map_err(|e| e.to_string())?;
        let budget = BudgetLimits {
            work: self.number(Field::Work)?,
            hits: self.number(Field::Hits)?,
            output_charge: self.number(Field::Output)?,
        };
        Budget::new(budget).map_err(|e| e.to_string())?;
        let sodium: f64 = self.number(Field::Sodium)?;
        let strands: f64 = self.number(Field::Strands)?;
        let conditions = ThermoConditions {
            sodium_molar: sodium / 1000.,
            total_strand_molar: strands * 1e-9,
        };
        conditions.validate().map_err(|e| e.to_string())?;
        Ok((collection, options, budget, conditions))
    }
}
pub fn parse_collection(text: &str) -> Result<Collection, String> {
    if text.len() > MAX_TEXT_BYTES {
        return Err("Oligo input exceeds 768 KiB.".into());
    }
    let text = text.trim();
    if text.is_empty() {
        return Err("Enter an oligo collection first.".into());
    }
    let mut oligos = Vec::new();
    let mut add = |name: &str, sequence: &str, tail: usize| -> Result<(), String> {
        if oligos.len() >= collagene_oligo_binding::MAX_OLIGOS
            || sequence.len() > collagene_oligo_binding::MAX_OLIGO_BASES
            || name.len() > 256
        {
            return Err("Collection exceeds oligo count, name or sequence bounds.".into());
        }
        let oligo = Oligo {
            id: format!("oligo-{}", oligos.len() + 1),
            name: name.trim().into(),
            ordered: Sequence::dna(sequence).map_err(|e| e.to_string())?,
            declared_tail_length: tail,
        };
        oligo.validate().map_err(|e| e.to_string())?;
        oligos.push(oligo);
        Ok(())
    };
    if text.starts_with('>') {
        let mut header = None::<String>;
        let mut sequence = String::new();
        for line in text.lines().map(str::trim).filter(|line| !line.is_empty()) {
            if let Some(next) = line.strip_prefix('>') {
                if let Some(previous) = header.take() {
                    let (name, tail) = fasta_header(&previous)?;
                    add(name, &sequence, tail)?;
                    sequence.clear();
                }
                if next.len() > 280 {
                    return Err("FASTA header is too long.".into());
                }
                header = Some(next.to_owned());
            } else {
                if header.is_none()
                    || sequence.len() + line.len() > collagene_oligo_binding::MAX_OLIGO_BASES
                {
                    return Err("Missing FASTA header or oligo exceeds 200 bases.".into());
                }
                sequence.push_str(line);
            }
        }
        if let Some(previous) = header {
            let (name, tail) = fasta_header(&previous)?;
            add(name, &sequence, tail)?;
        }
    } else {
        for (index, line) in text
            .lines()
            .filter(|line| !line.trim().is_empty())
            .enumerate()
        {
            let fields: Vec<_> = line.splitn(4, '\t').collect();
            if !(2..=3).contains(&fields.len()) {
                return Err(format!(
                    "Row {}: use tab-separated name, sequence, optional declared tail length.",
                    index + 1
                ));
            }
            let tail = if fields.len() == 3 {
                fields[2]
                    .trim()
                    .parse()
                    .map_err(|_| format!("Row {} has invalid tail length.", index + 1))?
            } else {
                0
            };
            add(fields[0], fields[1].trim(), tail)?;
        }
    }
    let collection = Collection { oligos };
    collection.validate().map_err(|e| e.to_string())?;
    Ok(collection)
}
fn fasta_header(header: &str) -> Result<(&str, usize), String> {
    let header = header.trim();
    if let Some((name, last)) = header.rsplit_once(char::is_whitespace) {
        if let Some(tail) = last.strip_prefix("tail=") {
            return Ok((
                name.trim(),
                tail.parse().map_err(|_| "Invalid FASTA tail length.")?,
            ));
        }
    }
    Ok((header, 0))
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Capture {
    pub source: SourceInfo,
    pub epoch: u64,
    pub generation: u32,
    pub draft: Draft,
    pub context_revision: String,
    pub source_revision: String,
}
impl Capture {
    pub fn current(
        &self,
        source: &SourceInfo,
        epoch: u64,
        generation: u32,
        draft: &Draft,
        source_revision: &str,
    ) -> bool {
        self.source == *source
            && self.epoch == epoch
            && self.generation == generation
            && self.draft == *draft
            && self.source_revision == source_revision
    }
}
pub fn status_label(status: Status, hits: usize) -> String {
    match status {
        Status::Unsearched => "Unsearched; absence not assessed".into(),
        Status::Limited(limit) => format!("Incomplete ({limit:?}); absence not established"),
        Status::Complete if hits == 0 => "Complete; no sites within the chosen model/policy".into(),
        Status::Complete => format!("Complete; {hits} retained alignments"),
    }
}
/// Convert original plus-reference coordinates to ascending half-open spans.
/// Reverse orientation changes column order, never the viewer's coordinate system.
pub fn spans(hit: &Hit) -> Vec<(usize, usize)> {
    let mut positions: Vec<_> = hit.columns.iter().map(|c| c.source_position).collect();
    positions.sort_unstable();
    positions.dedup();
    let mut spans = Vec::<(usize, usize)>::new();
    for p in positions {
        if let Some(last) = spans.last_mut() {
            if last.1 == p {
                last.1 = p + 1;
                continue;
            }
        }
        spans.push((p, p + 1));
    }
    spans
}
/// Stream a complete bounded export without constructing a second JSON value tree.
pub fn export_json(
    capture: &Capture,
    report: &collagene_core::cloning::pcr::collection_binding::Report,
) -> Result<Vec<u8>, String> {
    #[derive(serde::Serialize)]
    struct Evidence<'a> {
        schema_version: u32,
        source: &'a SourceInfo,
        context_revision: &'a str,
        collection: &'a Collection,
        options: &'a Options,
        budget: BudgetLimits,
        conditions: &'a ThermoConditions,
        report: &'a collagene_core::cloning::pcr::collection_binding::Report,
    }
    struct Bytes(Vec<u8>);
    impl std::io::Write for Bytes {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            const MAX: usize = 32 * 1024 * 1024;
            let needed = self
                .0
                .len()
                .checked_add(bytes.len())
                .filter(|n| *n <= MAX)
                .ok_or_else(|| std::io::Error::other("Binding export exceeds 32 MiB"))?;
            if needed > self.0.capacity() {
                let capacity = needed
                    .max(self.0.capacity().saturating_mul(2))
                    .max(1024)
                    .min(MAX);
                self.0
                    .try_reserve_exact(capacity - self.0.len())
                    .map_err(|_| std::io::Error::other("Binding export allocation failed"))?;
            }
            self.0.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let (collection, options, budget, conditions) = capture.draft.prepare()?;
    let evidence = Evidence {
        schema_version: 1,
        source: &capture.source,
        context_revision: &capture.context_revision,
        collection: &collection,
        options: &options,
        budget,
        conditions: &conditions,
        report,
    };
    let mut bytes = Bytes(Vec::new());
    serde_json::to_writer(&mut bytes, &evidence).map_err(|e| e.to_string())?;
    Ok(bytes.0)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn input_formats_preserve_declared_tail_separately() {
        let c = parse_collection(">forward tail=2\nGGACGT\n>reverse\nACGT").unwrap();
        assert_eq!(c.oligos[0].declared_tail_length, 2);
        assert_eq!(c.oligos[0].ordered.as_str(), "GGACGT");
        assert_eq!(
            parse_collection("forward\tGGACGT\t2").unwrap().oligos[0],
            c.oligos[0]
        );
    }
    #[test]
    fn malformed_or_unbounded_input_rejected() {
        for input in [
            ">empty\n",
            ">x tail=no\nACG",
            "x\tACG\t3",
            "x\tACG\t-1",
            "x,ACG,0",
        ] {
            assert!(parse_collection(input).is_err(), "{input}");
        }
        assert!(parse_collection(&"A".repeat(MAX_TEXT_BYTES + 1)).is_err());
    }
    #[test]
    fn terminal_mismatch_admission_is_explicit_and_zero_budget_is_valid() {
        let mut d = Draft::default();
        d.text = "x\tACGT".into();
        assert_eq!(d.prepare().unwrap().1.minimum_three_prime_pairing, 0);
        d.set(Field::Hits, "0".into());
        assert_eq!(d.prepare().unwrap().2.hits, 0);
        d.set(Field::Prime, "201".into());
        assert!(d.prepare().is_err());
    }
    #[test]
    fn incomplete_zero_hits_never_claim_absence() {
        assert!(status_label(Status::Unsearched, 0).contains("not assessed"));
        assert!(
            status_label(Status::Limited(collagene_oligo_binding::Limit::Work), 0)
                .contains("not established")
        );
        assert!(status_label(Status::Complete, 0).contains("within the chosen model"));
    }
}

#[cfg(test)]
mod lifecycle_tests {
    use super::*;
    #[test]
    fn exact_event_epoch_form_generation_and_source_revision_all_bind() {
        let source = SourceInfo {
            binding: collagene_core::material::lineage::EventBinding {
                document_id: "d".into(),
                event_id: "e".into(),
                revision: "r".into(),
            },
            name: "plasmid".into(),
            length: 10,
            alphabet: collagene_model::Alphabet::Dna,
            topology: collagene_model::Topology::Circular,
        };
        let capture = Capture {
            source: source.clone(),
            epoch: 5,
            generation: 2,
            draft: Draft::default(),
            context_revision: "context".into(),
            source_revision: "r".into(),
        };
        assert!(capture.current(&source, 5, 2, &capture.draft, "r"));
        let mut other = source.clone();
        other.binding.event_id = "other".into();
        assert!(!capture.current(&other, 5, 2, &capture.draft, "r"));
        assert!(!capture.current(&source, 6, 2, &capture.draft, "r"));
        assert!(!capture.current(&source, 5, 3, &capture.draft, "r"));
        assert!(!capture.current(&source, 5, 2, &capture.draft, "changed"));
        let mut changed = capture.draft.clone();
        changed.set(Field::Sodium, "100".into());
        assert!(!capture.current(&source, 5, 2, &changed, "r"));
        changed = capture.draft.clone();
        changed.text = "x\tACGT".into();
        assert!(!capture.current(&source, 5, 2, &changed, "r"));
    }
    #[test]
    fn original_coordinates_sort_reverse_origin_crossing_into_two_spans() {
        let c = parse_collection("x\tTAC").unwrap();
        let mut options = Options::default();
        options.minimum_annealing_length = 3;
        options.maximum_mismatches = 0;
        options.minimum_score = 0;
        let report = collagene_oligo_binding::search(
            &Sequence::dna("TACCG").unwrap(),
            collagene_model::Topology::Circular,
            &c,
            &options,
            &mut Budget::new(BudgetLimits::default()).unwrap(),
        )
        .unwrap();
        let hit = report
            .hits
            .iter()
            .find(|h| {
                h.orientation == collagene_oligo_binding::Orientation::Reverse && h.crosses_origin
            })
            .unwrap();
        assert_eq!(
            hit.columns
                .iter()
                .map(|c| c.source_position)
                .collect::<Vec<_>>(),
            vec![1, 0, 4]
        );
        assert_eq!(spans(hit), vec![(0, 2), (4, 5)]);
    }
}
