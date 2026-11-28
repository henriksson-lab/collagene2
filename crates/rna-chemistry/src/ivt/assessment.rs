//! Deterministic checks on declared template intervals; no biochemical simulator.
use super::*;

data! {
    #[derive(Copy)]
    pub enum Verdict {
        Pass,
        Fail,
        Indeterminate,
        Advisory,
        Informational,
    }
}
data! {
    #[derive(Copy)]
    pub enum Rule {
        ConsensusPromoter,
        CapInitiator,
        DeclaredInitiator,
        ConflictingRequirements,
        Linearization,
        RunoffBoundary,
        PolyAEncoding,
        PolyATerminus,
    }
}
data! {
    pub struct Finding {
        pub rule: Rule,
        pub verdict: Verdict,
        pub message: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 8>")]
        pub source_ids: Vec<String>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, MappedSpan, 17>")]
        pub locations: Vec<MappedSpan>,
    }
}
data! {
    pub struct Observation {
        /// Coding-strand DNA, including retained IUPAC ambiguity.
        pub dna: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 17>")]
        pub original_positions: Vec<usize>,
        pub complete: bool,
    }
}
data! {
    #[derive(Copy)]
    pub enum EndKind {
        Blunt,
        FivePrimeOverhang,
        ThreePrimeOverhang,
    }
}
data! {
    pub struct EndGeometry {
        pub kind: EndKind,
        pub overhang_bases: usize,
        /// Original forward-source bond coordinates. Linear bond n is permitted;
        /// circular bond n is canonicalized to 0. Strand roles follow Request.
        pub coding_cut_bond: usize,
        pub template_cut_bond: usize,
        pub coding_cut_offset: usize,
        pub template_cut_offset: usize,
    }
}
data! {
    #[derive(Copy)]
    pub enum RunKind {
        Adenosine,
        OtherConcrete,
        Ambiguous,
    }
}
data! {
    pub struct PolyARun {
        pub kind: RunKind,
        pub start_offset: usize,
        pub end_offset: usize,
        /// Preview only. Offsets/counts and source retain the complete encoding.
        pub dna_preview: String,
        pub preview_truncated: bool,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, MappedSpan, 2>")]
        pub source_spans: Vec<MappedSpan>,
    }
}
data! {
    pub struct PolyAEncoding {
        pub region: PolyARegion,
        pub adenosines: usize,
        pub other_concrete_bases: usize,
        pub ambiguous_bases: usize,
        pub longest_exact_a_run: usize,
        pub terminal: bool,
        pub runs_total: usize,
        pub segmentation_complete: bool,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, PolyARun, 1024>")]
        pub runs: Vec<PolyARun>,
    }
}
data! {
    pub struct Report {
        pub engine_version: String,
        pub context_revision: String,
        pub source_revision: String,
        pub source_length: usize,
        pub source_topology: Topology,
        pub profile: Profile,
        pub request: Request,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, MappedSpan, 2>")]
        pub transcript_spans: Vec<MappedSpan>,
        pub promoter: Observation,
        pub initiator_context: Observation,
        #[serde(deserialize_with = "required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub linearization: Option<EndGeometry>,
        #[serde(deserialize_with = "required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub poly_a: Option<PolyAEncoding>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Finding, 16>")]
        pub findings: Vec<Finding>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 8>")]
        pub assumptions: Vec<String>,
        pub report_revision: String,
    }
}
impl Report {
    fn content_revision(&self) -> Result<String> {
        hash(&(
            &self.engine_version,
            &self.context_revision,
            &self.source_revision,
            self.source_length,
            self.source_topology,
            &self.profile,
            &self.request,
            &self.transcript_spans,
            &self.promoter,
            &self.initiator_context,
            &self.linearization,
            &self.poly_a,
            &self.findings,
            &self.assumptions,
        ))
    }
    /// Documentary integrity only. A hash cannot authenticate scientific claims
    /// supplied by an untrusted producer; use validate_for for fresh replay.
    pub fn validate_binding(&self, expected: &str) -> Result<()> {
        bounded(self, MAX_REPORT_BYTES)?;
        self.request.validate_shape()?;
        self.profile.validate()?;
        if self.engine_version != ENGINE
            || self.source_revision != self.request.source_revision
            || self.context_revision != expected
            || self.context_revision != hash(&(ENGINE, &self.source_revision, &self.request))?
            || self.source_length == 0
            || self.source_length > MAX_TEMPLATE_BASES
            || self.transcript_spans.is_empty()
            || self.transcript_spans.len() > 2
            || self.promoter.dna.len() > 17
            || self.promoter.original_positions.len() != self.promoter.dna.len()
            || self.initiator_context.dna.len() > 8
            || self.initiator_context.original_positions.len() != self.initiator_context.dna.len()
            || self.findings.len() > 16
            || self.assumptions.len() > 8
            || self.assumptions.iter().any(|text| text.len() > 2048)
            || self.report_revision != self.content_revision()?
        {
            return Err(Error::Binding);
        }
        let transcript = &self.request.transcript;
        if transcript.start >= self.source_length
            || transcript.length > self.source_length
            || (self.source_topology == Topology::Linear
                && match transcript.coding_strand {
                    CodingStrand::Forward => {
                        transcript.length > self.source_length - transcript.start
                    }
                    CodingStrand::Reverse => transcript.length > transcript.start + 1,
                })
            || self.transcript_spans.iter().try_fold(0usize, |sum, span| {
                sum.checked_add(span.end.saturating_sub(span.start))
            }) != Some(transcript.length)
            || (self.promoter.complete && self.promoter.dna.len() != 17)
            || self.linearization.is_some() != self.request.linearization.is_some()
            || self.poly_a.as_ref().map(|encoding| &encoding.region) != self.request.poly_a.as_ref()
        {
            return Err(Error::Binding);
        }
        let span_valid = |span: &MappedSpan| {
            span.start < span.end
                && span.end <= self.source_length
                && span.reverse == (self.request.transcript.coding_strand == CodingStrand::Reverse)
        };
        if self.transcript_spans.iter().any(|span| !span_valid(span))
            || self
                .promoter
                .original_positions
                .iter()
                .chain(&self.initiator_context.original_positions)
                .any(|position| *position >= self.source_length)
            || self.findings.iter().any(|finding| {
                finding.message.len() > 2048
                    || finding.source_ids.len() > 8
                    || finding.locations.len() > 17
                    || finding.locations.iter().any(|span| !span_valid(span))
                    || finding
                        .source_ids
                        .iter()
                        .any(|id| !self.profile.sources.iter().any(|source| &source.id == id))
            })
        {
            return Err(Error::Binding);
        }
        if let Some(poly_a) = &self.poly_a {
            if poly_a.runs.len() > MAX_RETAINED_RUNS
                || poly_a.runs.len() > poly_a.runs_total
                || poly_a.segmentation_complete != (poly_a.runs.len() == poly_a.runs_total)
                || poly_a.region.start_offset >= poly_a.region.end_offset
                || poly_a.region.end_offset > self.request.transcript.length
                || poly_a.runs.iter().any(|run| {
                    run.start_offset < poly_a.region.start_offset
                        || run.end_offset > poly_a.region.end_offset
                        || run.start_offset >= run.end_offset
                        || run.dna_preview.len() > 64
                        || run.source_spans.len() > 2
                        || run.source_spans.iter().any(|span| !span_valid(span))
                })
            {
                return Err(Error::Binding);
            }
        }
        Ok(())
    }
    /// Independent current computation, including complete source/options binding.
    pub fn validate_for(&self, source: &Molecule, request: &Request) -> Result<()> {
        let expected = revision(source, request)?;
        self.validate_binding(&expected)?;
        if self != &assess(source, request)? {
            return Err(Error::Binding);
        }
        Ok(())
    }
}

fn observe(source: &Molecule, transcript: &Transcript, start: isize, length: usize) -> Observation {
    let mut observation = Observation {
        dna: String::with_capacity(length),
        original_positions: Vec::with_capacity(length),
        complete: true,
    };
    for offset in start..start + length as isize {
        match (
            transcript.base(source, offset),
            transcript.position(source, offset),
        ) {
            (Some(base), Some(position)) => {
                observation.dna.push(base as char);
                observation.original_positions.push(position);
            }
            _ => observation.complete = false,
        }
    }
    observation
}
fn exact_requirement(observed: &Observation, expected: &str) -> Verdict {
    if !observed.complete {
        return Verdict::Indeterminate;
    }
    if observed
        .dna
        .bytes()
        .zip(expected.bytes())
        .any(|(actual, required)| mask(actual) & mask(required) == 0)
    {
        return Verdict::Fail;
    }
    if !observed.complete || observed.dna.len() != expected.len() || observed.dna != expected {
        return Verdict::Indeterminate;
    }
    Verdict::Pass
}
fn observation_locations(observation: &Observation, strand: CodingStrand) -> Vec<MappedSpan> {
    // These observations have at most17 bases. One-base intervals also describe
    // a clipped promoter honestly without inventing omitted positions.
    observation
        .original_positions
        .iter()
        .map(|position| MappedSpan {
            start: *position,
            end: position + 1,
            reverse: strand == CodingStrand::Reverse,
        })
        .collect()
}
fn finding(
    rule: Rule,
    verdict: Verdict,
    message: &str,
    source_ids: &[&str],
    locations: Vec<MappedSpan>,
) -> Finding {
    Finding {
        rule,
        verdict,
        message: message.into(),
        source_ids: source_ids.iter().map(|id| (*id).to_owned()).collect(),
        locations,
    }
}
fn poly_a(
    source: &Molecule,
    transcript: &Transcript,
    region: &PolyARegion,
) -> Result<PolyAEncoding> {
    let mut report = PolyAEncoding {
        region: region.clone(),
        adenosines: 0,
        other_concrete_bases: 0,
        ambiguous_bases: 0,
        longest_exact_a_run: 0,
        terminal: region.end_offset == transcript.length,
        runs_total: 0,
        segmentation_complete: true,
        runs: Vec::new(),
    };
    let mut at = region.start_offset;
    while at < region.end_offset {
        let first = transcript
            .base(source, at as isize)
            .ok_or(Error::Invalid("missing declared poly(A) base"))?;
        let classify = |base| {
            if base == b'A' {
                RunKind::Adenosine
            } else if b"CGT".contains(&base) {
                RunKind::OtherConcrete
            } else {
                RunKind::Ambiguous
            }
        };
        let kind = classify(first);
        let start = at;
        let mut preview = String::new();
        while at < region.end_offset {
            let base = transcript
                .base(source, at as isize)
                .ok_or(Error::Invalid("missing declared poly(A) base"))?;
            if classify(base) != kind {
                break;
            }
            if preview.len() < 64 && report.runs.len() < MAX_RETAINED_RUNS {
                preview.push(base as char);
            }
            at += 1;
        }
        let length = at - start;
        match kind {
            RunKind::Adenosine => {
                report.adenosines += length;
                report.longest_exact_a_run = report.longest_exact_a_run.max(length);
            }
            RunKind::OtherConcrete => report.other_concrete_bases += length,
            RunKind::Ambiguous => report.ambiguous_bases += length,
        }
        report.runs_total += 1;
        if report.runs.len() < MAX_RETAINED_RUNS {
            report.runs.push(PolyARun {
                kind,
                start_offset: start,
                end_offset: at,
                preview_truncated: preview.len() != length,
                dna_preview: preview,
                source_spans: transcript.map_interval(source, start, at)?,
            });
        } else {
            report.segmentation_complete = false;
        }
    }
    Ok(report)
}

pub fn assess(source: &Molecule, request: &Request) -> Result<Report> {
    let context_revision = revision(source, request)?;
    let profile = profiles::current()?;
    let transcript = &request.transcript;
    let transcript_spans = transcript.map_interval(source, 0, transcript.length)?;
    let promoter = observe(source, transcript, -17, 17);
    // Prevent a promoter wrapping into the same bases as the declared RNA.
    let promoter_separate = source.topology != Topology::Circular
        || transcript
            .length
            .checked_add(17)
            .is_some_and(|length| length <= source.sequence.len());
    let mut findings = vec![finding(
        Rule::ConsensusPromoter,
        if promoter_separate {
            exact_requirement(&promoter, &profile.promoter_dna)
        } else {
            Verdict::Fail
        },
        if promoter_separate {
            "Match to the explicitly selected T7 consensus at -17..-1; a failed match is a profile conflict, not proof of zero transcription."
        } else {
            "Declared transcript and upstream promoter overlap within a single source turn; this runoff-template geometry cannot retain a separate promoter."
        },
        &["neb-t7-ag-start"],
        observation_locations(&promoter, transcript.coding_strand),
    )];
    let initiator_context = observe(source, transcript, 0, transcript.length.min(8));
    let cap_observation = observe(source, transcript, 0, transcript.length.min(2));
    let (cap_verdict, cap_message) = match request.cap {
        CapChoice::CleanCapAgN7113 => (
            exact_requirement(&cap_observation, "AG"),
            "Original CleanCap AG N-7113 requires AG at the declared +1/+2; this does not establish cap incorporation or efficiency.",
        ),
        CapChoice::NotSpecified => (
            Verdict::Indeterminate,
            "Cap chemistry was not specified; no initiator compatibility can be established.",
        ),
        CapChoice::Uncapped | CapChoice::EnzymaticPostTranscriptional => (
            Verdict::Informational,
            "No co-transcriptional cap-analogue sequence rule is imposed by this declaration; initiation and actual capping remain separate.",
        ),
    };
    findings.push(finding(
        Rule::CapInitiator,
        cap_verdict,
        cap_message,
        &["trilink-n7113-ag", "neb-t7-ag-start"],
        observation_locations(&cap_observation, transcript.coding_strand),
    ));
    let required = request.initiator.required_dna();
    if !required.is_empty() {
        let observed = observe(source, transcript, 0, transcript.length.min(required.len()));
        findings.push(finding(Rule::DeclaredInitiator, exact_requirement(&observed, required), "Check the caller's additional explicit initiator constraint; this is not a universal T7 efficiency or necessity rule.", &["conrad-2020-context"], observation_locations(&observed, transcript.coding_strand)));
    }
    if request.cap == CapChoice::CleanCapAgN7113 && !required.is_empty() {
        findings.push(finding(Rule::ConflictingRequirements, Verdict::Fail, "The separately requested G-start constraint conflicts with the selected AG cap initiator; no source sequence can satisfy both at +1.", &["neb-t7-ag-start", "trilink-n7113-ag"], observation_locations(&cap_observation, transcript.coding_strand)));
    }
    let linearization = if let Some(geometry) = &request.linearization {
        let kind = match geometry
            .coding_cut_offset
            .cmp(&geometry.template_cut_offset)
        {
            std::cmp::Ordering::Equal => EndKind::Blunt,
            std::cmp::Ordering::Less => EndKind::FivePrimeOverhang,
            std::cmp::Ordering::Greater => EndKind::ThreePrimeOverhang,
        };
        findings.push(finding(Rule::Linearization, if kind == EndKind::ThreePrimeOverhang { Verdict::Advisory } else { Verdict::Pass }, "Supplied downstream-end geometry only: blunt/5-prime overhangs follow the supplier's runoff guidance; a 3-prime overhang is discouraged. Cleavage, intact promoter retention and physical template preparation are not established.", &["neb-runoff-ends"], vec![]));
        findings.push(finding(Rule::RunoffBoundary, if geometry.template_cut_offset == transcript.length { Verdict::Pass } else { Verdict::Fail }, "Compare declared transcript end with the supplied template-strand cut offset. Actual end homogeneity/non-templated additions are not predicted.", &["neb-runoff-ends", "conrad-2020-context"], vec![]));
        Some(EndGeometry {
            kind,
            overhang_bases: geometry
                .coding_cut_offset
                .abs_diff(geometry.template_cut_offset),
            coding_cut_bond: transcript.bond(source, geometry.coding_cut_offset),
            template_cut_bond: transcript.bond(source, geometry.template_cut_offset),
            coding_cut_offset: geometry.coding_cut_offset,
            template_cut_offset: geometry.template_cut_offset,
        })
    } else {
        findings.push(finding(Rule::Linearization, Verdict::Indeterminate, "No downstream strand-cut geometry was supplied. A circular source is not a linearized runoff template; a sequence record's linear topology alone does not establish the requested physical end.", &["neb-runoff-ends"], vec![]));
        None
    };
    let poly_a = request
        .poly_a
        .as_ref()
        .map(|region| poly_a(source, transcript, region))
        .transpose()?;
    if let Some(encoding) = &poly_a {
        let locations = transcript.map_interval(
            source,
            encoding.region.start_offset,
            encoding.region.end_offset,
        )?;
        findings.push(finding(Rule::PolyAEncoding, if encoding.ambiguous_bases > 0 { Verdict::Indeterminate } else { Verdict::Informational }, "Exact declared encoded-region composition and contiguous runs only. Segmentation does not predict host recombination, RNA half-life, expression or actual polyadenylation.", &[], locations.clone()));
        findings.push(finding(Rule::PolyATerminus, if encoding.terminal { Verdict::Pass } else { Verdict::Fail }, "A declared poly(A) tail must reach the declared transcript's 3-prime boundary; downstream template-encoded bases are not silently removed.", &[], locations));
    }
    let mut report = Report {
        engine_version: ENGINE.into(), context_revision, source_revision: source.revision(),
        source_length: source.sequence.len(), source_topology: source.topology,
        profile, request: request.clone(), transcript_spans, promoter, initiator_context,
        linearization, poly_a, findings,
        assumptions: vec![
            "Coordinates declare coding-strand 5-prime to 3-prime order on the exact original DNA; the complementary strand is the polymerase template. At most one source turn is used.".into(),
            "Only the selected consensus and original N-7113 AG rule are admitted. Other cap analogues, polymerases, promoter variants and chemical modifications require separate rules.".into(),
            "Supplied cut geometry is documentary and must be established independently by the restriction/material owners before use. This module performs no digestion, RNA reconstruction or library mutation.".into(),
            "No yield, cap occupancy, end homogeneity, dsRNA, immunogenicity, expression or propagation-stability prediction is supplied. GGG-start experimental observations are not transferred to AG capping.".into(),
        ], report_revision: String::new(),
    };
    report.report_revision = report.content_revision()?;
    bounded(&report, MAX_REPORT_BYTES)?;
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
