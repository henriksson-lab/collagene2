//! Portable presentation of a retained comparison. This does not rerun alignment
//! or turn a historical receipt or host freshness label into scientific authority.
use crate::{
    alignment::history::{self as pairwise_history, Record},
    formats::alignment::{
        self as alignment_formats, AlignmentFormat, AlignmentMetadata, AlignmentRow,
    },
    output::pdf::ReviewBlock,
};
use collagene_model::{Alphabet, Error};
use serde::{Deserialize, Serialize};
use std::io::{self, Write};

pub const MAX_EXPORT_BYTES: usize = 16 * 1024 * 1024;
const MAX_NAME_BYTES: usize = 4096;
const CONTRACT: &str = "pairwise-export-3";

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Format {
    AlignedFasta,
    Clustal,
    Html,
    Pdf,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Export {
    pub bytes: Vec<u8>,
    pub mime: String,
    pub extension: String,
    pub warnings: Vec<String>,
}
impl Export {
    /// Bounded response-shape guard, not validation of scientific correctness.
    pub fn validate(&self) -> Result<(), Error> {
        if self.bytes.is_empty()
            || self.bytes.len() > MAX_EXPORT_BYTES
            || self.warnings.len() > 64
            || self.warnings.iter().any(|s| s.len() > 65_536)
        {
            return Err(Error::TooLarge);
        }
        let valid = match (self.extension.as_str(), self.mime.as_str()) {
            ("fasta", "text/plain;charset=utf-8") => self.bytes.starts_with(b">reference "),
            ("aln", "text/plain;charset=utf-8") => self.bytes.starts_with(b"CLUSTAL "),
            ("html", "text/html;charset=utf-8") => self.bytes.starts_with(b"<!doctype html>"),
            ("pdf", "application/pdf") => self.bytes.starts_with(b"%PDF-"),
            _ => false,
        };
        if !valid || (self.extension != "pdf" && std::str::from_utf8(&self.bytes).is_err()) {
            return Err(bad("invalid export MIME, extension or content header"));
        }
        Ok(())
    }
}
fn bad(message: &str) -> Error {
    Error::Format(format!("pairwise export: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct Counter(usize);
impl Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|n| *n <= pairwise_history::MAX_SAVED_BYTES)
            .ok_or_else(|| io::Error::other("comparison input byte limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
/// Borrowed input guard before copying/encoding a comparison. Event ownership
/// and actual source letters remain the host's retained-workspace contract.
pub fn preflight(record: &Record, reference_name: &str, query_name: &str) -> Result<(), Error> {
    for name in [reference_name, query_name] {
        if name.trim().is_empty()
            || name.len() > MAX_NAME_BYTES
            || name.chars().any(char::is_control)
        {
            return Err(bad("invalid source display name"));
        }
    }
    if record.schema_version != pairwise_history::VERSION {
        return Err(Error::Unsupported("pairwise export record schema".into()));
    }
    for binding in [&record.reference, &record.query] {
        if binding.document_id.is_empty()
            || binding.document_id.len() > 128
            || !binding
                .document_id
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
            || !digest(&binding.revision)
            || !binding.event_id.strip_prefix("event-").is_some_and(digest)
        {
            return Err(bad("invalid exact source event"));
        }
    }
    let r = &record.report;
    let a = &r.comparison.alignment;
    if record.reference.revision != r.reference_revision
        || record.query.revision != r.query_revision
    {
        return Err(bad("record source revisions disagree with report"));
    }
    if r.comparison.warnings.len() > 4096
        || a.warnings.len() > 4096
        || r.comparison
            .warnings
            .iter()
            .chain(&a.warnings)
            .any(|s| s.len() > 65_536)
        || a.reference_aligned.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
        || a.query_aligned.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
        || a.operations.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
        || a.differences.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
    {
        return Err(Error::TooLarge);
    }
    serde_json::to_writer(&mut Counter(0), record).map_err(|_| Error::TooLarge)?;
    record.validate_documentary()
}
pub fn revision(
    record: &Record,
    reference_name: &str,
    query_name: &str,
    fresh: bool,
    format: Format,
) -> Result<String, Error> {
    preflight(record, reference_name, query_name)?;
    crate::encoding::bounded_json::hash(&(
        CONTRACT,
        record,
        reference_name,
        query_name,
        fresh,
        format,
    ))
}

struct Blocks {
    values: Vec<ReviewBlock>,
    bytes: usize,
    limit: usize,
}
impl Blocks {
    fn push(&mut self, block: ReviewBlock) -> Result<(), Error> {
        let size = match &block {
            ReviewBlock::Heading(s) | ReviewBlock::Paragraph(s) => s.len(),
            ReviewBlock::Fixed(lines) => lines.iter().map(String::len).sum(),
        };
        self.bytes = self
            .bytes
            .checked_add(size)
            .filter(|n| *n <= self.limit)
            .ok_or(Error::TooLarge)?;
        self.values.push(block);
        Ok(())
    }
    fn paragraph(&mut self, text: String) -> Result<(), Error> {
        self.push(ReviewBlock::Paragraph(text))
    }
    fn heading(&mut self, text: &str) -> Result<(), Error> {
        self.push(ReviewBlock::Heading(text.into()))
    }
}
fn semantic_blocks(
    out: &mut Blocks,
    evidence: &crate::alignment::design_consensus::Evidence,
) -> Result<(), Error> {
    crate::encoding::bounded_json::preflight(evidence, out.limit.saturating_sub(out.bytes))?;
    out.heading("Design versus supplied consensus: conditional evidence")?;
    out.paragraph(evidence.statement.clone())?;
    out.paragraph(format!(
        "Design event: {} / {} / {}. Supplied consensus event: {} / {} / {}.",
        evidence.request.design.document_id,
        evidence.request.design.event_id,
        evidence.request.design.revision,
        evidence.request.consensus.document_id,
        evidence.request.consensus.event_id,
        evidence.request.consensus.revision
    ))?;
    out.paragraph(format!("Axis assumption: {:?}. Semantic context: {}. Evidence: {}. Original workspace fingerprint: {}. Coding complete: {}. Selected path conditional: {}.",
        evidence.request.consensus_axis, evidence.context_revision, evidence.evidence_revision,
        evidence.workspace_fingerprint, evidence.coding_complete, evidence.selected_path_conditional))?;
    for source in &evidence.provenance {
        out.paragraph(format!("Source {:?}: topology {:?}. Original evidence {:?}, format {:?}, bytes {}, hash {:?}. {}",
            source.role, source.topology, source.original.state, source.original.format,
            source.original.byte_count, source.original.blake3, source.original.message))?;
        for field in &source.fields {
            out.paragraph(format!(
                "Recorded source field {} (value complete: {}): {}",
                field.key, field.value_complete, field.value
            ))?;
        }
    }
    out.paragraph(format!("Coding/mapping work: mapping {}, coding calls {}, register segments {}, consequences {}. Limits: {:?}.",
        evidence.work.mapping_work, evidence.work.coding_calls, evidence.work.segments,
        evidence.work.consequences, evidence.request.limits))?;
    for row in &evidence.cds {
        out.heading("Declared design CDS")?;
        out.paragraph(format!("ID: {}. Annotation {} / {} / {}. Rationale: {}. Genetic code: {}. Mapping complete: {}; register complete: {}; consequences complete: {}.",
            row.id, row.design_cds.feature_index, row.design_cds.source_revision,
            row.design_cds.annotation_revision, row.design_cds.rationale, row.code.descriptor.id,
            row.mapping_complete, row.register_complete, row.consequences_complete))?;
        out.paragraph(format!("Genetic-code pin: {}. Descriptor: {:?}. Ordinary amino acids in TCAG codon enumeration order: {}.",
            row.code.revision, row.code.descriptor, row.code.ordinary_amino_acids))?;
        for (role, translation) in [
            ("Design", &row.design_translation),
            ("Projected consensus", &row.projected_translation),
        ] {
            if let Some(translation) = translation {
                out.paragraph(format!("{role} translation: {}. Protein: {}. Codons {}, trailing {}, ambiguous {}, unresolved {}, possible stops {}.",
                    translation.code_name, translation.protein, translation.translated_codons,
                    translation.trailing_bases, translation.ambiguous_codons,
                    translation.unresolved_codons, translation.possible_stop_codons))?;
                for warning in &translation.warnings {
                    out.paragraph(format!("{role} translation warning: {warning}"))?;
                }
            } else {
                out.paragraph(format!("{role} translation unavailable."))?;
            }
        }
        for finding in &row.findings {
            out.paragraph(format!(
                "{:?} / {:?}: {}",
                finding.rule, finding.verdict, finding.message
            ))?;
        }
        for segment in &row.register_segments {
            out.paragraph(format!("Coding ordinals after codon-start prefix: design [{}, {}) / projected [{}, {}); offset modulo3 {}.",
                segment.design_start, segment.design_end, segment.projected_start, segment.projected_end, segment.offset_mod3))?;
        }
        for change in &row.consequences {
            out.paragraph(format!("Conditional codon consequence {:?}: projected {} {} ({}); design {:?} {:?} ({:?}). Original design spans {:?}; original consensus spans {:?}. {}",
                change.kind, change.projected_codon_index, change.projected_codon, change.projected_amino_acid,
                change.design_codon_index, change.design_codon, change.design_amino_acid,
                change.design_coordinates, change.consensus_coordinates, change.message))?;
        }
    }
    Ok(())
}
fn label(fresh: bool) -> &'static str {
    if fresh {
        "Fresh comparison (host context)"
    } else {
        "Historical comparison - not recomputed on opening"
    }
}
fn blocks(
    record: &Record,
    rname: &str,
    qname: &str,
    fresh: bool,
    pdf: bool,
) -> Result<Vec<ReviewBlock>, Error> {
    let report = &record.report;
    let c = &report.comparison;
    let a = &c.alignment;
    let limit = if pdf {
        crate::output::pdf::MAX_PDF_TEXT_BYTES / 2
    } else {
        MAX_EXPORT_BYTES / 2
    };
    // Refuse an oversized report before copying its large aligned/difference text.
    let mut minimum = a
        .reference_aligned
        .len()
        .checked_add(a.query_aligned.len())
        .ok_or(Error::TooLarge)?;
    for d in &a.differences {
        minimum = minimum
            .checked_add(d.reference_bases.len())
            .and_then(|n| n.checked_add(d.query_bases.len()))
            .ok_or(Error::TooLarge)?;
    }
    if minimum > limit {
        return Err(Error::TooLarge);
    }
    let mut b = Blocks {
        values: Vec::new(),
        bytes: 0,
        limit,
    };
    b.heading("Sequence comparison")?;
    b.paragraph(label(fresh).into())?;
    b.paragraph("This is recorded sequence-model evidence, not sample identity, read confidence or experimental success. Export does not rerun alignment. Displayed row coordinates use the declared rotation/orientation; native source mappings are stated separately. Missing local flanks are not verified.".into())?;
    for (role, name, event) in [
        ("Reference", rname, &record.reference),
        ("Query", qname, &record.query),
    ] {
        b.heading(role)?;
        b.paragraph(format!(
            "Name: {name}\nDocument: {}\nEvent: {}\nRevision: {}",
            event.document_id, event.event_id, event.revision
        ))?;
    }
    b.heading("Calculation and uncertainty")?;
    b.paragraph(format!(
        "Engine: {}\nAlgorithm: {}\nContext revision: {}\nOptions: {}",
        report.engine_version,
        a.algorithm,
        report.context_revision,
        serde_json::to_string(&report.options).map_err(|e| bad(&e.to_string()))?
    ))?;
    b.paragraph(format!("Source lengths: reference {}, query {}. Topologies: reference {:?}, query {:?}. Reference origin offset {}; query reverse complement {}. Evaluated candidates {}. Exhaustive requested reference origins {}; exhaustive requested query orientations {}; best score proven {}.",report.reference_length,report.query_length,c.reference_topology,c.query_topology,c.reference_origin,c.query_reversed,c.evaluated_candidates,c.exhaustive_reference_origins,c.exhaustive_query_orientations,c.optimality_proven))?;
    b.paragraph(format!("Actual cumulative DP cells across candidates: {}. Selected-path unrestricted optimum proven: {}. Selected-path storage {:?}; stop {:?}; attempts {}. No proof enumerates all tied correspondences.",c.evaluated_cells,a.search.unrestricted_optimum_proven,a.search.storage,a.search.stop_reason,a.search.attempts.len()))?;
    for (index, candidate) in c.candidate_evidence.iter().enumerate() {
        b.paragraph(format!("Candidate {}: origin {}; reverse {}; score {:?}; disposition {:?}; actual cells {}. Individual candidate optimum proven {}.",index+1,candidate.reference_origin,candidate.query_reversed,candidate.score,candidate.disposition,candidate.search.as_ref().map_or(0,|search|search.evaluated_cells),candidate.search.as_ref().is_some_and(|search|search.unrestricted_optimum_proven)))?;
        if let Some(search) = &candidate.search {
            for attempt in &search.attempts {
                b.paragraph(format!("  Domain {:?}; cells {}; score {}; omitted-path bound {:?}; own unrestricted optimum proven {}.",attempt.domain,attempt.evaluated_cells,attempt.score,attempt.off_band_bound,attempt.unrestricted_optimum_proven))?;
            }
        }
    }
    b.paragraph(format!("Score {}. Reported aligned intervals (zero-based half-open, oriented): reference [{}, {}), query [{}, {}). Paired coverage: reference {:.6}%, query {:.6}%. Certain matches {}; compatible ambiguous {}; disjoint mismatches {}; inserted query bases {}; deleted reference bases {}. Unaligned flanks: reference {} left / {} right; query {} left / {} right.",a.score,a.reference_start,a.reference_end,a.query_start,a.query_end,a.stats.reference_coverage*100.,a.stats.query_coverage*100.,a.stats.certain_matches,a.stats.ambiguous_compatible,a.stats.mismatches,a.stats.inserted_query_bases,a.stats.deleted_reference_bases,a.reference_start,report.reference_length-a.reference_end,a.query_start,report.query_length-a.query_end))?;
    for warning in c.warnings.iter().chain(&a.warnings) {
        b.paragraph(format!("Warning: {warning}"))?;
    }
    if let Some(evidence) = &record.design_consensus {
        semantic_blocks(&mut b, evidence)?;
    }
    b.heading("Complete aligned columns")?;
    b.paragraph("Rows are shown exactly as returned, including gaps. R/Q numbers are zero-based oriented boundaries. Marks: | certain match; ? compatible ambiguity; X disjoint mismatch; blank at gaps.".into())?;
    if a.reference_aligned.is_empty() {
        b.paragraph("No aligned columns were returned. Neither complete source nor its unaligned flanks are represented as matches.".into())?;
    }
    let (mut r, mut q) = (a.reference_start, a.query_start);
    for start in (0..a.reference_aligned.len()).step_by(60) {
        let end = (start + 60).min(a.reference_aligned.len());
        let rs = &a.reference_aligned[start..end];
        let qs = &a.query_aligned[start..end];
        let marks: String = rs
            .bytes()
            .zip(qs.bytes())
            .map(|(r, q)| {
                if r == b'-' || q == b'-' {
                    ' '
                } else {
                    let (r, q) = (
                        crate::sequence::analysis::base_mask(r),
                        crate::sequence::analysis::base_mask(q),
                    );
                    if r & q == 0 {
                        'X'
                    } else if r.count_ones() == 1 && q.count_ones() == 1 {
                        '|'
                    } else {
                        '?'
                    }
                }
            })
            .collect();
        b.push(ReviewBlock::Fixed(vec![
            format!("R {r:>9} {rs}"),
            format!("            {marks}"),
            format!("Q {q:>9} {qs}"),
        ]))?;
        r += rs.bytes().filter(|b| *b != b'-').count();
        q += qs.bytes().filter(|b| *b != b'-').count();
    }
    b.heading("Complete differences and native coordinates")?;
    b.paragraph("Coordinates below are zero-based half-open. An empty [x,x) interval is an insertion/deletion caret, not an invented base. Circular native reference spans retain ordered parts. Query native intervals retain their source-forward boundaries even when displayed in reverse-complement orientation.".into())?;
    for (index, d) in a.differences.iter().enumerate() {
        let rr = report.map_reference(d.reference_start, d.reference_end)?;
        let qq = report.map_query(d.query_start, d.query_end)?;
        b.paragraph(format!("Difference {}: {:?}. Oriented R [{}, {}), Q [{}, {}). Native R {:?}; native Q [{}, {}). Reference calls: {}. Query calls: {}.",index+1,d.kind,d.reference_start,d.reference_end,d.query_start,d.query_end,rr,qq.0,qq.1,if d.reference_bases.is_empty(){"(empty)"}else{&d.reference_bases},if d.query_bases.is_empty(){"(empty)"}else{&d.query_bases}))?;
    }
    b.heading("Complete CIGAR")?;
    b.paragraph(if a.cigar.is_empty() {
        "(empty)".into()
    } else {
        a.cigar.clone()
    })?;
    Ok(b.values)
}

struct Limited(Vec<u8>);
impl Write for Limited {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(bytes.len())
            .is_none_or(|n| n > MAX_EXPORT_BYTES)
        {
            return Err(io::Error::other("comparison output byte limit"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn escaped(out: &mut Limited, text: &str) -> Result<(), Error> {
    for c in text.chars() {
        let s = match c {
            '&' => "&amp;",
            '<' => "&lt;",
            '>' => "&gt;",
            '\"' => "&quot;",
            '\'' => "&#39;",
            _ => {
                let mut b = [0; 4];
                out.write_all(c.encode_utf8(&mut b).as_bytes())
                    .map_err(|_| Error::TooLarge)?;
                continue;
            }
        };
        out.write_all(s.as_bytes()).map_err(|_| Error::TooLarge)?;
    }
    Ok(())
}
fn html(blocks: &[ReviewBlock]) -> Result<Vec<u8>, Error> {
    let mut out = Limited(Vec::new());
    out.write_all(b"<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width\"><title>Collagene 2 sequence comparison</title><style>body{max-width:80ch;margin:2em auto;padding:0 1em;font:14px/1.5 sans-serif}pre{font:12px/1.5 monospace;white-space:pre;overflow:auto}p{white-space:pre-wrap;overflow-wrap:anywhere}h2{break-after:avoid}pre{break-inside:avoid}@media print{body{max-width:none;margin:0}pre{overflow:visible}}</style><body>").map_err(|_|Error::TooLarge)?;
    for block in blocks {
        match block {
            ReviewBlock::Heading(text) => {
                out.write_all(b"<h2>").map_err(|_| Error::TooLarge)?;
                escaped(&mut out, text)?;
                out.write_all(b"</h2>").map_err(|_| Error::TooLarge)?;
            }
            ReviewBlock::Paragraph(text) => {
                out.write_all(b"<p>").map_err(|_| Error::TooLarge)?;
                escaped(&mut out, text)?;
                out.write_all(b"</p>").map_err(|_| Error::TooLarge)?;
            }
            ReviewBlock::Fixed(lines) => {
                out.write_all(b"<pre>").map_err(|_| Error::TooLarge)?;
                for line in lines {
                    escaped(&mut out, line)?;
                    out.write_all(b"\n").map_err(|_| Error::TooLarge)?;
                }
                out.write_all(b"</pre>").map_err(|_| Error::TooLarge)?;
            }
        }
    }
    out.write_all(b"</body></html>")
        .map_err(|_| Error::TooLarge)?;
    Ok(out.0)
}

/// Export the actual reported aligned path. `fresh` is a guarded host presentation
/// label, not a newly computed verdict. The host resolves exact retained sources
/// and guards its form/workspace before and after asynchronous generation.
pub fn export(
    record: &Record,
    reference_name: &str,
    query_name: &str,
    fresh: bool,
    format: Format,
) -> Result<Export, Error> {
    preflight(record, reference_name, query_name)?;
    let mut warnings=vec!["Portable comparison output is sequence-model evidence, not experimental verification. Export does not recompute the alignment.".into()];
    let (bytes, mime, extension) = match format {
        Format::AlignedFasta | Format::Clustal => {
            let a = &record.report.comparison.alignment;
            if a.reference_aligned.is_empty() {
                return Err(Error::Unsupported("an empty local alignment has no exchangeable rows; use HTML/PDF or complete comparison JSON".into()));
            }
            let rows=[(&record.reference,"reference",&a.reference_aligned),(&record.query,"query",&a.query_aligned)].into_iter().enumerate().map(|(index,(source,id,aligned))|AlignmentRow{id:format!("computed-{id}"),source_row_index:index,identifier:id.into(),description:format!("document={} event={} revision={} oriented_start={} oriented_end={} reference_origin={} query_reverse={}",source.document_id,source.event_id,source.revision,if index==0{a.reference_start}else{a.query_start},if index==0{a.reference_end}else{a.query_end},record.report.comparison.reference_origin,record.report.comparison.query_reversed),aligned:aligned.clone()}).collect::<Vec<_>>();
            let alphabet = if a.reference_aligned.contains('U') || a.query_aligned.contains('U') {
                Alphabet::Rna
            } else {
                Alphabet::Dna
            };
            let result = alignment_formats::export_rows(
                &rows,
                &AlignmentMetadata::default(),
                alphabet,
                if format == Format::AlignedFasta {
                    AlignmentFormat::AlignedFasta
                } else {
                    AlignmentFormat::Clustal
                },
            )?;
            warnings.extend(result.losses);
            if record.design_consensus.is_some() {
                warnings.push("This aligned text format omits design/consensus semantic roles, provenance, axis assumptions, conditional CDS/register/consequence evidence and semantic pins. Retain the complete comparison JSON or HTML/PDF.".into());
            }
            warnings.push("Standard aligned rows retain gaps and the displayed orientation/origin, but omit native coordinate maps, source display names, complete options, statuses, warnings and other receipt evidence. Keep the HTML/PDF or full JSON with this file. Identifiers reference/query designate roles; they are not renamed native documents.".into());
            (
                result.text.into_bytes(),
                "text/plain;charset=utf-8",
                if format == Format::AlignedFasta {
                    "fasta"
                } else {
                    "aln"
                },
            )
        }
        Format::Html => (
            html(&blocks(record, reference_name, query_name, fresh, false)?)?,
            "text/html;charset=utf-8",
            "html",
        ),
        Format::Pdf => {
            let blocks = blocks(record, reference_name, query_name, fresh, true)?;
            let (bytes, _pages, extra) = crate::output::pdf::export_documentary_pdf(
                crate::output::pdf::DocumentaryKind::SequenceComparison,
                "Collagene 2 sequence comparison",
                label(fresh),
                &record.reference.revision,
                &record.query.revision,
                &blocks,
            )?;
            warnings.extend(extra);
            (bytes, "application/pdf", "pdf")
        }
    };
    let result = Export {
        bytes,
        mime: mime.into(),
        extension: extension.into(),
        warnings,
    };
    result.validate()?;
    Ok(result)
}
