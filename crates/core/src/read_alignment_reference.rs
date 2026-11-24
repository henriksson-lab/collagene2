//! Explicit source-bound checking of sealed imported SAM/BAM alignment assertions.
//! No reference-name/URI lookup, variant calling, consensus or uniqueness inference.
use crate::{Alphabet, Error, Molecule, Topology, analysis::base_mask};
use collagene_read_alignments::{
    CigarKind, Format, Header, QueryProjection, Record, SequenceSymbol,
};
use serde::{Deserialize, Serialize};
use std::io::Write;

pub const ENGINE_VERSION: &str = "read-alignment-reference-1";
const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
const MAX_MANAGED_BYTES: usize = 64 * 1024 * 1024;
macro_rules! dto {($($item:tt)*)=>{#[derive(Clone,Debug,PartialEq,Eq,Serialize,Deserialize)]#[serde(deny_unknown_fields)]$($item)*};}
macro_rules! enumeration {($($item:tt)*)=>{#[derive(Clone,Copy,Debug,PartialEq,Eq,Serialize,Deserialize)]$($item)*};}

dto! {pub struct Options {
    /// Explicit caller mapping of this dictionary name to the supplied molecule.
    pub reference_name:String,
    /// Must equal the supplied molecule's canonical revision; no name authorizes it.
    pub expected_reference_revision:String,
    pub max_header_bytes:usize,pub max_record_bytes:usize,
    pub max_query_bases:usize,pub max_runs:usize,pub max_columns:usize,
    pub max_work:u64,pub max_report_bytes:usize,
}}
impl Default for Options {
    fn default() -> Self {
        Self {
            reference_name: String::new(),
            expected_reference_revision: String::new(),
            max_header_bytes: 16 * 1024 * 1024,
            max_record_bytes: 4 * 1024 * 1024,
            max_query_bases: 1_000_000,
            max_runs: 100_000,
            max_columns: 100_000,
            max_work: 128_000_000,
            max_report_bytes: 16 * 1024 * 1024,
        }
    }
}
enumeration! {pub enum Incomplete {
    MissingExplicitReferenceMapping,ReferenceRevisionMismatch,UnsupportedReferenceAlphabet,
    MissingReferenceDictionary,AmbiguousReferenceName,RecordUsesDifferentReference,
    ReferenceLengthMismatch,ReferenceTopologyMismatch,Unmapped,UnavailableCoordinates,
    MissingSequence,UnresolvedCalls,ReferenceEqualityOutsideAlignment,ConflictingImportedAssertions,
}}
enumeration! {pub enum Operation {Match,Equal,Difference,Insertion,Deletion,Skip,SoftClip,HardClip,Pad}}
enumeration! {pub enum CallOrigin {ReportedIupac,ReferenceEqualResolved,UnresolvedReportedSymbol}}
enumeration! {pub enum Relation {
    CertainAgreement,CertainDisagreement,PossibleAgreement,
    /// Derived from this reference, not an independently observed match.
    ReferenceEqualResolved,InsertionCall,SoftClippedCall,Unresolved,
}}
dto! {pub struct Column {
    pub query_index:u64,pub implied_original_query:Option<u64>,
    /// Exact lexical stored symbol, including lowercase or '='; never replaced in Record.
    pub reported_symbol:char,pub query_symbol:Option<char>,pub call_origin:CallOrigin,
    pub reference_position:Option<usize>,pub unwrapped_reference_position:Option<u64>,
    pub reference_symbol:Option<char>,pub raw_phred:Option<u8>,pub relation:Relation,
    /// Only concrete observed calls can contradict CIGAR '=' / 'X'.
    pub contradicts_cigar:bool,
}}
dto! {pub struct Run {
    pub operation:Operation,pub query_start:u64,pub query_end:u64,
    pub reference_start:u64,pub reference_end:u64,
    pub implied_original_query_span:Option<(u64,u64)>,
    /// Half-open indexes into Report.columns. D/N/H/P produce no columns.
    pub column_start:usize,pub column_end:usize,
}}
dto! {#[derive(Default)] pub struct Summary {
    pub checked_aligned_calls:u64,pub certain_agreements:u64,pub certain_disagreements:u64,
    pub possible_agreements:u64,pub reference_equal_resolved:u64,pub unresolved_calls:u64,
    pub cigar_contradictions:u64,pub inserted_calls:u64,pub soft_clipped_calls:u64,
    pub deleted_reference_bases:u64,pub skipped_reference_bases:u64,
    pub hard_clipped_absent_bases:u64,pub padding_units:u64,
}}
dto! {pub struct Binding {
    pub reference_name:String,pub reference_revision:String,
    pub header_revision:String,pub record_original_blake3:String,pub options_revision:String,
    pub format:String,
}}
// Transport/historical DTO only: decoding cannot restore the private record or
// collection. Fresh checks still require the sealed source plus actual DNA.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine_version: String,
    pub revision: String,
    pub binding: Binding,
    pub reference_bound: bool,
    pub complete: bool,
    pub incomplete: Vec<Incomplete>,
    pub reference_length: usize,
    pub topology: Topology,
    pub reverse: bool,
    pub reported_flags: u16,
    pub reported_mapping_quality: Option<u8>,
    pub reported_query_length: u64,
    pub implied_original_length: Option<u64>,
    pub sequence_present: bool,
    pub qualities_present: bool,
    pub runs: Vec<Run>,
    pub columns: Vec<Column>,
    pub summary: Summary,
    pub warnings: Vec<String>,
}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("read alignment reference: {message}"))
}
fn add(a: usize, b: usize) -> Result<usize, Error> {
    a.checked_add(b).ok_or(Error::TooLarge)
}
fn mul(a: usize, b: usize) -> Result<usize, Error> {
    a.checked_mul(b).ok_or(Error::TooLarge)
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
impl Options {
    pub fn validate(&self) -> Result<(), Error> {
        if self.reference_name.len() > 4096
            || (!self.expected_reference_revision.is_empty()
                && !hex(&self.expected_reference_revision))
            || self.max_header_bytes == 0
            || self.max_header_bytes > 64 * 1024 * 1024
            || self.max_record_bytes == 0
            || self.max_record_bytes > 32 * 1024 * 1024
            || self.max_query_bases == 0
            || self.max_query_bases > 1_000_000
            || self.max_runs == 0
            || self.max_runs > 100_000
            || self.max_columns == 0
            || self.max_columns > 1_000_000
            || self.max_work == 0
            || self.max_work > 1_000_000_000
            || self.max_report_bytes < 16384
            || self.max_report_bytes > MAX_REPORT_BYTES
        {
            return Err(bad("invalid resource/mapping policy"));
        }
        Ok(())
    }
}
struct JsonHash {
    hash: blake3::Hasher,
    bytes: usize,
    max: usize,
}
impl Write for JsonHash {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|&n| n <= self.max)
            .ok_or_else(|| std::io::Error::other("bounded reference report"))?;
        self.hash.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn json_hash(value: &impl Serialize, max: usize) -> Result<(String, usize), Error> {
    let mut out = JsonHash {
        hash: blake3::Hasher::new(),
        bytes: 0,
        max,
    };
    serde_json::to_writer(&mut out, value).map_err(|_| Error::TooLarge)?;
    Ok((out.hash.finalize().to_hex().to_string(), out.bytes))
}
fn header_bytes(header: &Header) -> Result<usize, Error> {
    if header.references().len() > 100_000 {
        return Err(Error::TooLarge);
    }
    let mut bytes = add(
        header.sam_bytes().len(),
        header.original_bam_header().map_or(0, |b| b.len()),
    )?;
    for reference in header.references() {
        bytes = add(bytes, add(reference.name.len(), 16)?)?;
    }
    Ok(bytes)
}
fn header_revision(header: &Header) -> String {
    let mut hash = blake3::Hasher::new();
    hash.update(b"read-alignment-full-header-v1");
    hash.update(&(header.sam_bytes().len() as u64).to_le_bytes());
    hash.update(header.sam_bytes());
    hash.update(&[u8::from(header.original_bam_header().is_some())]);
    if let Some(bam) = header.original_bam_header() {
        hash.update(&(bam.len() as u64).to_le_bytes());
        hash.update(bam);
    }
    hash.update(&(header.references().len() as u64).to_le_bytes());
    for reference in header.references() {
        hash.update(&(reference.name.len() as u64).to_le_bytes());
        hash.update(reference.name.as_bytes());
        hash.update(&reference.length.to_le_bytes());
        hash.update(&[u8::from(reference.circular)]);
    }
    hash.finalize().to_hex().to_string()
}
fn operation(kind: CigarKind) -> Operation {
    match kind {
        CigarKind::Match => Operation::Match,
        CigarKind::Equal => Operation::Equal,
        CigarKind::Difference => Operation::Difference,
        CigarKind::Insertion => Operation::Insertion,
        CigarKind::Deletion => Operation::Deletion,
        CigarKind::Skip => Operation::Skip,
        CigarKind::SoftClip => Operation::SoftClip,
        CigarKind::HardClip => Operation::HardClip,
        CigarKind::Pad => Operation::Pad,
    }
}
fn incomplete(report: &mut Report, reason: Incomplete) {
    if !report.incomplete.contains(&reason) {
        report.incomplete.push(reason);
    }
    report.complete = false;
}
fn finish(report: Report, options: &Options) -> Result<Report, Error> {
    json_hash(&report, options.max_report_bytes)?;
    Ok(report)
}

/// Check one immutable imported record against one explicitly chosen actual DNA
/// reference. The host must retain its original Record/Header and source event;
/// this source-bound report neither changes them nor attests unique placement.
pub fn check_record(
    header: &Header,
    record: &Record,
    reference: &Molecule,
    options: &Options,
) -> Result<Report, Error> {
    options.validate()?;
    collagene_model::bounds::model_preflight(&[reference])?;
    let header_size = header_bytes(header)?;
    let fields = record.fields();
    let query_size = fields.sequence.as_ref().map_or(0, String::len);
    let run_count = fields.cigar.len();
    if header_size > options.max_header_bytes
        || record.original_bytes().len() > options.max_record_bytes
        || run_count > options.max_runs
        || query_size > options.max_query_bases
        || query_size > options.max_columns
    {
        return Err(Error::TooLarge);
    }
    let minimum_work = add(
        add(header_size, record.original_bytes().len())?,
        add(reference.sequence.len(), run_count)?,
    )?;
    if minimum_work as u64 > options.max_work {
        return Err(Error::TooLarge);
    }
    // Charge both the sealed CoordinateMap and new report storage before either
    // operation-map or query-column vectors can be allocated. No D/N expansion.
    let memory = add(
        mul(
            run_count,
            std::mem::size_of::<collagene_read_alignments::MapRun>() + std::mem::size_of::<Run>(),
        )?,
        mul(query_size, std::mem::size_of::<Column>())?,
    )?;
    let report_bound = add(16384, add(mul(run_count, 256)?, mul(query_size, 512)?)?)?;
    if memory > MAX_MANAGED_BYTES || report_bound > options.max_report_bytes {
        return Err(Error::TooLarge);
    }
    // Header/Record are sealed; the core rejects a foreign full header context.
    let map = record
        .coordinate_map(header)
        .map_err(|e| bad(&e.to_string()))?;
    let (reference_revision, source_bytes) = json_hash(reference, 32 * 1024 * 1024)?;
    let logs = if run_count < 2 {
        1
    } else {
        usize::BITS as usize - (run_count - 1).leading_zeros() as usize + 1
    };
    let work = add(
        add(minimum_work, source_bytes)?,
        add(
            add(mul(query_size, logs + 4)?, mul(run_count, 4)?)?,
            report_bound,
        )?,
    )?;
    if work as u64 > options.max_work {
        return Err(Error::TooLarge);
    }
    let options_revision = json_hash(options, 16384)?.0;
    let binding = Binding {
        reference_name: options.reference_name.clone(),
        reference_revision,
        header_revision: header_revision(header),
        record_original_blake3: record.original_blake3(),
        options_revision,
        format: match record.format() {
            Format::Sam => "SAM",
            Format::Bam => "BAM",
        }
        .into(),
    };
    let revision = json_hash(&(ENGINE_VERSION, &binding), 16384)?.0;
    let mut report=Report {engine_version:ENGINE_VERSION.into(),revision,binding,
        reference_bound:false,complete:true,incomplete:vec![],reference_length:reference.sequence.len(),topology:reference.topology,
        reverse:map.reverse(),reported_flags:fields.flags,reported_mapping_quality:fields.mapping_quality,
        reported_query_length:map.reported_query_length(),implied_original_length:map.implied_original_length(),
        sequence_present:map.sequence_present(),qualities_present:fields.qualities.is_some(),runs:vec![],columns:vec![],summary:Summary::default(),
        warnings:vec!["Only imported alignment geometry and actual bound reference symbols are checked. MAPQ, M5/UR, NM/MD, secondary/supplementary status and modification tags are not validated. No network/reference lookup, alignment uniqueness, consensus, coverage, variant truth or biological accuracy is inferred.".into(),
            "SEQ '=' is reference-dependent imported evidence: its new resolved symbol is not an independently observed matching nucleotide. FLAG16 uses stored alignment orientation; implied original indexes include absent hard clips, without generating their calls.".into()]};
    if reference.sequence.alphabet() != Alphabet::Dna {
        incomplete(&mut report, Incomplete::UnsupportedReferenceAlphabet);
    } else {
        reference.validate()?;
    }
    if options.reference_name.is_empty() || options.expected_reference_revision.is_empty() {
        incomplete(&mut report, Incomplete::MissingExplicitReferenceMapping);
    } else if options.expected_reference_revision != report.binding.reference_revision {
        incomplete(&mut report, Incomplete::ReferenceRevisionMismatch);
    }
    if fields.reference.as_deref() != Some(options.reference_name.as_str()) {
        incomplete(&mut report, Incomplete::RecordUsesDifferentReference);
    }
    let matches = header
        .references()
        .iter()
        .filter(|r| r.name == options.reference_name)
        .count();
    if matches == 0 || !map.dictionary_bound() {
        incomplete(&mut report, Incomplete::MissingReferenceDictionary);
    } else if matches != 1 {
        incomplete(&mut report, Incomplete::AmbiguousReferenceName);
    }
    if let Some(named) = header.reference(&options.reference_name) {
        if named.length as usize != reference.sequence.len() {
            incomplete(&mut report, Incomplete::ReferenceLengthMismatch);
        }
        if named.circular != (reference.topology == Topology::Circular) {
            incomplete(&mut report, Incomplete::ReferenceTopologyMismatch);
        }
    }
    if fields.flags & 4 != 0 {
        incomplete(&mut report, Incomplete::Unmapped);
    } else if !map.usable_reported_mapping() {
        incomplete(&mut report, Incomplete::UnavailableCoordinates);
    }
    if !report.complete {
        return finish(report, options);
    }
    report.reference_bound = true;
    if !map.sequence_present() {
        incomplete(&mut report, Incomplete::MissingSequence);
    }
    report.runs = Vec::with_capacity(run_count);
    report.columns = Vec::with_capacity(query_size);
    for run in map.runs() {
        let column_start = report.columns.len();
        match run.kind {
            CigarKind::Deletion => {
                report.summary.deleted_reference_bases += run.reference_end - run.reference_start
            }
            CigarKind::Skip => {
                report.summary.skipped_reference_bases += run.reference_end - run.reference_start
            }
            CigarKind::HardClip => {
                report.summary.hard_clipped_absent_bases +=
                    u64::from(fields.cigar[report.runs.len()].len)
            }
            CigarKind::Pad => {
                report.summary.padding_units += u64::from(fields.cigar[report.runs.len()].len)
            }
            _ => {}
        }
        if map.sequence_present() {
            for query in run.query_start..run.query_end {
                let query_index = usize::try_from(query).map_err(|_| Error::TooLarge)?;
                let symbol = record
                    .reported_symbol(query_index)
                    .ok_or_else(|| bad("sealed query length disagrees with stored calls"))?;
                let projected = map.query_position(query);
                let (position, unwrapped, original) = match projected {
                    QueryProjection::Aligned {
                        reference_position,
                        unwrapped_position,
                        implied_original_query,
                        ..
                    } => (
                        Some(reference_position as usize),
                        Some(unwrapped_position),
                        implied_original_query,
                    ),
                    QueryProjection::Inserted {
                        implied_original_query,
                        ..
                    }
                    | QueryProjection::SoftClipped {
                        implied_original_query,
                    } => (None, None, implied_original_query),
                    _ => {
                        return Err(bad(
                            "sealed query projection unavailable inside checked run",
                        ));
                    }
                };
                let reference_symbol =
                    position.map(|p| reference.sequence.as_str().as_bytes()[p] as char);
                let raw_phred = fields.qualities.as_ref().map(|q| q[query_index]);
                let (reported_symbol, query_symbol, origin) = match symbol {
                    SequenceSymbol::Iupac(base) => (
                        base as char,
                        Some(base.to_ascii_uppercase() as char),
                        CallOrigin::ReportedIupac,
                    ),
                    SequenceSymbol::ReferenceEqual => (
                        '=',
                        reference_symbol,
                        if position.is_some() {
                            CallOrigin::ReferenceEqualResolved
                        } else {
                            CallOrigin::UnresolvedReportedSymbol
                        },
                    ),
                    SequenceSymbol::Unresolved(base) => {
                        (base as char, None, CallOrigin::UnresolvedReportedSymbol)
                    }
                };
                let aligned = position.is_some();
                if aligned && reported_symbol == '=' && run.kind == CigarKind::Difference {
                    // Two imported claims conflict. Resolution is still clearly
                    // reference-derived, never an independently observed match.
                    incomplete(&mut report, Incomplete::ConflictingImportedAssertions);
                }
                let relation = if origin == CallOrigin::ReferenceEqualResolved {
                    report.summary.reference_equal_resolved += 1;
                    Relation::ReferenceEqualResolved
                } else if query_symbol.is_none() {
                    report.summary.unresolved_calls += 1;
                    incomplete(
                        &mut report,
                        if reported_symbol == '=' {
                            Incomplete::ReferenceEqualityOutsideAlignment
                        } else {
                            Incomplete::UnresolvedCalls
                        },
                    );
                    Relation::Unresolved
                } else if let Some(reference_symbol) = reference_symbol {
                    let a = base_mask(query_symbol.unwrap() as u8);
                    let b = base_mask(reference_symbol as u8);
                    if a & b == 0 {
                        report.summary.certain_disagreements += 1;
                        Relation::CertainDisagreement
                    } else if a.count_ones() == 1 && b.count_ones() == 1 {
                        report.summary.certain_agreements += 1;
                        Relation::CertainAgreement
                    } else {
                        report.summary.possible_agreements += 1;
                        Relation::PossibleAgreement
                    }
                } else if run.kind == CigarKind::Insertion {
                    Relation::InsertionCall
                } else {
                    Relation::SoftClippedCall
                };
                if aligned {
                    report.summary.checked_aligned_calls += 1;
                }
                if run.kind == CigarKind::Insertion {
                    report.summary.inserted_calls += 1;
                }
                if run.kind == CigarKind::SoftClip {
                    report.summary.soft_clipped_calls += 1;
                }
                let contradicts_cigar = run.kind == CigarKind::Equal
                    && relation == Relation::CertainDisagreement
                    || run.kind == CigarKind::Difference && relation == Relation::CertainAgreement;
                report.summary.cigar_contradictions += u64::from(contradicts_cigar);
                report.columns.push(Column {
                    query_index: query,
                    implied_original_query: original,
                    reported_symbol,
                    query_symbol,
                    call_origin: origin,
                    reference_position: position,
                    unwrapped_reference_position: unwrapped,
                    reference_symbol,
                    raw_phred,
                    relation,
                    contradicts_cigar,
                });
            }
        }
        report.runs.push(Run {
            operation: operation(run.kind),
            query_start: run.query_start,
            query_end: run.query_end,
            reference_start: run.reference_start,
            reference_end: run.reference_end,
            implied_original_query_span: run.original_query_span,
            column_start,
            column_end: report.columns.len(),
        });
    }
    finish(report, options)
}
