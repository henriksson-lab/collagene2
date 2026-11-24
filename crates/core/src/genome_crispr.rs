//! Streaming CRISPR geometry over immutable large-reference blocks.
//! No serialized progress or imported report grants coverage authority.
use crate::{
    Error, Molecule, Topology,
    crispr::{self, GuideOptions, GuideQuery, NucleaseProfile, PamSide},
    genome::{
        self, BlockVerifier, EncodedBlock, GenomeManifest, GenomicInterval, SemanticValidator,
    },
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "genome-crispr-1";
pub const MAX_INTERVAL_BASES: u64 = 200_000;
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
dto! {pub struct BlockRequest{pub contig_id:String,pub index:u64}}
dto! {pub struct GenomeGuideOptions{pub source_role:String,pub interval:GenomicInterval,pub guide:GuideOptions}}
dto! {pub struct GenomeGuideCandidate{pub id:String,pub contig_id:String,pub guide_sequence:String,pub pam_sequence:String,pub guide_positions:Vec<u64>,pub pam_positions:Vec<u64>,pub reverse:bool,pub pam_certain:bool,pub concrete:bool,pub gc_bp:Option<u16>,pub longest_homopolymer:Option<usize>,pub passes_filters:bool,pub nominal_top_cut:Option<u64>,pub nominal_bottom_cut:Option<u64>,pub rejections:Vec<String>}}
dto! {pub struct GenomeGuideReport{pub engine_version:String,pub manifest_hash:String,pub content_revision:String,pub options_revision:String,pub options:GenomeGuideOptions,pub candidates:Vec<GenomeGuideCandidate>,pub complete:bool,pub comparisons:u64,pub limits:Vec<String>,pub warnings:Vec<String>}}
dto! {pub struct GenomeIntended{pub reference_role:String,pub manifest_hash:String,pub content_revision:String,pub contig_id:String,pub guide_positions:Vec<u64>,pub pam_positions:Vec<u64>,pub reverse:bool}}
dto! {pub struct GenomeSearchOptions{pub reference_role:String,pub profile:NucleaseProfile,pub max_mismatches:usize,pub intended:Option<GenomeIntended>,pub max_total_comparisons:u64,pub max_index_entries_per_block:usize,pub max_hits:usize}}
impl Default for GenomeSearchOptions {
    fn default() -> Self {
        Self {
            reference_role: "reference".into(),
            profile: crispr::spcas9(),
            max_mismatches: 3,
            intended: None,
            max_total_comparisons: 1_000_000_000_000,
            max_index_entries_per_block: 1_000_000,
            max_hits: 4096,
        }
    }
}
dto! {pub struct GenomeOffTargetHit{pub id:String,pub contig_id:String,pub reference_guide:String,pub pam_sequence:String,pub guide_positions:Vec<u64>,pub pam_positions:Vec<u64>,pub reverse:bool,pub min_mismatches:usize,pub max_mismatches:usize,pub pam_certain:bool,pub definite:bool,pub is_intended:bool,pub nominal_top_cut:Option<u64>,pub nominal_bottom_cut:Option<u64>}}
dto! {pub struct GenomeSearchProgress{pub engine_version:String,pub manifest_hash:String,pub content_revision:String,pub query:GuideQuery,pub options:GenomeSearchOptions,pub blocks_processed:u64,pub blocks_total:u64,pub bases_verified:u64,pub comparisons:u64,pub hits_retained:usize,pub finished:bool,pub cancelled:bool,pub complete:bool,pub exhaustive_for_manifest:bool,pub expected_components_complete:bool,pub source_integrity_verified:bool,pub next_blocks:Vec<BlockRequest>,pub limits:Vec<String>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("genome CRISPR: {s}"))
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    Ok(
        blake3::hash(&postcard::to_allocvec(v).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
fn role(s: &str) -> Result<(), Error> {
    if s.is_empty() || s.len() > 256 || s.chars().any(char::is_control) {
        Err(bad("invalid explicit document/reference role"))
    } else {
        Ok(())
    }
}
fn contig<'a>(m: &'a GenomeManifest, id: &str) -> Result<&'a genome::ContigRecord, Error> {
    m.contigs
        .iter()
        .find(|c| c.id == id)
        .ok_or_else(|| bad("contig role absent from manifest"))
}
fn footprint(p: &NucleaseProfile) -> usize {
    p.guide_length + p.pam.len()
}
fn warnings() -> Vec<String> {
    vec!["Declared PAM and ungapped Hamming geometry only; chromatin, population variants, bulges and experimental activity are not assessed.".into(),"Complete manifest coverage does not attest that provider build labels name a complete real genome. Missing expected components prevent whole-build absence claims.".into()]
}
impl GenomeGuideOptions {
    pub fn validate(&self, m: &GenomeManifest) -> Result<(), Error> {
        m.validate()?;
        role(&self.source_role)?;
        self.guide.validate()?;
        let c = contig(m, &self.interval.contig_id)?;
        let width = self
            .interval
            .end
            .checked_sub(self.interval.start)
            .ok_or_else(|| bad("interval end precedes start"))?;
        if width == 0
            || width > MAX_INTERVAL_BASES
            || width > c.length
            || self.interval.start >= c.length
            || c.topology == Topology::Linear && self.interval.end > c.length
        {
            return Err(bad(
                "interval must have 1..200000 owner starts within one contig traversal; circular end may be unwrapped",
            ));
        }
        Ok(())
    }
}
fn requests(c: &genome::ContigRecord, start: u64, count: u64) -> Vec<BlockRequest> {
    let mut out = vec![];
    let mut at = start;
    let mut left = count;
    while left > 0 {
        let pos = at % c.length;
        let bi = pos / genome::BLOCK_BASES as u64;
        let b = &c.blocks[bi as usize];
        let n = left.min(b.start + u64::from(b.count) - pos);
        let r = BlockRequest {
            contig_id: c.id.clone(),
            index: bi,
        };
        if !out.contains(&r) {
            out.push(r)
        }
        at += n;
        left -= n;
    }
    out
}
fn window_count(c: &genome::ContigRecord, start: u64, owners: u64, p: &NucleaseProfile) -> u64 {
    let wanted = owners + footprint(p) as u64 - 1;
    if c.topology == Topology::Circular {
        wanted
    } else {
        wanted.min(c.length - start)
    }
}
pub fn interval_required_blocks(
    m: &GenomeManifest,
    o: &GenomeGuideOptions,
) -> Result<Vec<BlockRequest>, Error> {
    o.validate(m)?;
    let c = contig(m, &o.interval.contig_id)?;
    Ok(requests(
        c,
        o.interval.start,
        window_count(
            c,
            o.interval.start,
            o.interval.end - o.interval.start,
            &o.guide.profile,
        ),
    ))
}
fn verified_window(
    v: &BlockVerifier,
    c: &genome::ContigRecord,
    start: u64,
    count: u64,
    encoded: &[EncodedBlock],
) -> Result<Vec<u8>, Error> {
    let req = requests(c, start, count);
    if encoded.len() != req.len() {
        return Err(bad("provide exactly the requested owner/halo blocks"));
    }
    let mut blocks = vec![];
    for (r, e) in req.iter().zip(encoded) {
        if r.contig_id != e.contig_id || r.index != e.descriptor.index {
            return Err(bad("block input order or role mismatch"));
        }
        let block = v.verify(&r.contig_id, r.index, &e.encoded)?;
        if block.descriptor() != &e.descriptor {
            return Err(bad("supplied block descriptor differs from manifest"));
        }
        blocks.push(block)
    }
    let mut dna = Vec::with_capacity(count as usize);
    let mut at = start;
    let mut left = count;
    while left > 0 {
        let pos = at % c.length;
        let bi = pos / genome::BLOCK_BASES as u64;
        let block = blocks
            .iter()
            .find(|b| b.descriptor().index == bi)
            .ok_or_else(|| bad("missing required window block"))?;
        let offset = (pos - block.descriptor().start) as usize;
        let n = (left as usize).min(block.dna().len() - offset);
        dna.extend_from_slice(&block.dna()[offset..offset + n]);
        left -= n as u64;
        at += n as u64;
    }
    Ok(dna)
}
fn molecule(dna: &[u8]) -> Result<Molecule, Error> {
    Molecule::dna(
        "Verified reference window",
        std::str::from_utf8(dna).map_err(|_| bad("invalid decoded DNA"))?,
        Topology::Linear,
    )
}
fn global_positions(local: &[usize], start: u64, len: u64) -> Vec<u64> {
    local.iter().map(|p| (start + *p as u64) % len).collect()
}
fn cuts(
    p: &NucleaseProfile,
    forward_start: u64,
    reverse: bool,
    c: &genome::ContigRecord,
) -> (Option<u64>, Option<u64>) {
    let go = if p.pam_side == PamSide::ThreePrime {
        0
    } else {
        p.pam.len()
    } as i128;
    let f = forward_start as i128;
    let total = footprint(p) as i128;
    let t = go + i128::from(p.nominal_top_offset);
    let b = go + i128::from(p.nominal_bottom_offset);
    let (t, b) = if reverse {
        (f + total - b, f + total - t)
    } else {
        (f + t, f + b)
    };
    let project = |n: i128| {
        if c.topology == Topology::Circular {
            Some(n.rem_euclid(c.length as i128) as u64)
        } else if n >= 0 && n <= c.length as i128 {
            Some(n as u64)
        } else {
            None
        }
    };
    (project(t), project(b))
}
fn global_id(
    mh: &str,
    role: &str,
    c: &str,
    profile: &NucleaseProfile,
    g: &[u64],
    p: &[u64],
    reverse: bool,
) -> Result<String, Error> {
    hash(&(ENGINE_VERSION, mh, role, c, profile, g, p, reverse))
}
pub fn design_interval(
    m: &GenomeManifest,
    o: &GenomeGuideOptions,
    blocks: &[EncodedBlock],
) -> Result<GenomeGuideReport, Error> {
    o.validate(m)?;
    let v = BlockVerifier::begin(m)?;
    let mh = m.manifest_hash()?;
    let c = contig(m, &o.interval.contig_id)?;
    let owners = o.interval.end - o.interval.start;
    let count = window_count(c, o.interval.start, owners, &o.guide.profile);
    let dna = verified_window(&v, c, o.interval.start, count, blocks)?;
    let local = molecule(&dna)?;
    let local_revision = crispr::guide_revision(&local, &o.guide)?;
    let n = footprint(&o.guide.profile);
    let mut candidates = vec![];
    let mut comparisons = 0u64;
    let mut complete = true;
    let mut limits = vec![];
    'outer: for f in 0..owners as usize {
        if n as u64 > c.length || f + n > dna.len() {
            break;
        }
        for reverse in [false, true] {
            if comparisons + n as u64 > o.guide.max_comparisons as u64 {
                complete = false;
                limits.push("Interval guide comparison budget reached.".into());
                break 'outer;
            }
            comparisons += n as u64;
            let oriented = if reverse { dna.len() - f - n } else { f };
            if let Some(hit) =
                crispr::candidate_at(&local, &o.guide, oriented, reverse, &local_revision)?
            {
                if candidates.len() == o.guide.max_candidates {
                    complete = false;
                    limits.push("Interval candidate retention limit reached; returned prefix is incomplete.".into());
                    break 'outer;
                }
                let g = global_positions(&hit.guide_positions, o.interval.start, c.length);
                let p = global_positions(&hit.pam_positions, o.interval.start, c.length);
                let (top, bottom) = cuts(&o.guide.profile, o.interval.start + f as u64, reverse, c);
                candidates.push(GenomeGuideCandidate {
                    id: global_id(
                        &mh,
                        &o.source_role,
                        &c.id,
                        &o.guide.profile,
                        &g,
                        &p,
                        reverse,
                    )?,
                    contig_id: c.id.clone(),
                    guide_sequence: hit.guide_sequence,
                    pam_sequence: hit.pam_sequence,
                    guide_positions: g,
                    pam_positions: p,
                    reverse,
                    pam_certain: hit.pam_certain,
                    concrete: hit.concrete,
                    gc_bp: hit.gc_bp,
                    longest_homopolymer: hit.longest_homopolymer,
                    passes_filters: hit.passes_filters,
                    nominal_top_cut: top,
                    nominal_bottom_cut: bottom,
                    rejections: hit.rejections,
                });
            }
        }
    }
    let mut warnings = warnings();
    warnings.push("Interval blocks passed local layout and encoded-content checks; this operation does not verify the entire contig's declared canonical sequence hash.".into());
    Ok(GenomeGuideReport {
        engine_version: ENGINE_VERSION.into(),
        manifest_hash: mh,
        content_revision: m.content_revision.clone(),
        options_revision: hash(o)?,
        options: o.clone(),
        candidates,
        complete,
        comparisons,
        limits,
        warnings,
    })
}
/// Private live cursor; intentionally has no serialization/clone implementation.
pub struct GenomeSearchSession {
    manifest: GenomeManifest,
    verifier: BlockVerifier,
    semantic: Option<SemanticValidator>,
    query: GuideQuery,
    options: GenomeSearchOptions,
    manifest_hash: String,
    ci: usize,
    bi: usize,
    blocks_processed: u64,
    bases_verified: u64,
    comparisons: u64,
    hits: Vec<GenomeOffTargetHit>,
    finished: bool,
    cancelled: bool,
    source_integrity_verified: bool,
    limits: Vec<String>,
}
impl GenomeSearchOptions {
    pub fn validate(&self, m: &GenomeManifest, q: &GuideQuery) -> Result<(), Error> {
        m.validate()?;
        role(&self.reference_role)?;
        let o = crispr::OffTargetOptions {
            profile: self.profile.clone(),
            max_mismatches: self.max_mismatches,
            intended: None,
            reference_scope: "Verified owner block stream".into(),
            max_comparisons: 500_000_000,
            max_index_entries: self.max_index_entries_per_block,
            max_hits: self.max_hits,
        };
        o.validate()?;
        if self.max_total_comparisons == 0
            || self.max_total_comparisons > 10_000_000_000_000
            || q.id.is_empty()
            || q.id.len() > 256
            || q.id.chars().any(char::is_control)
            || q.sequence.len() != self.profile.guide_length
            || !q.sequence.bytes().all(|b| b"ACGT".contains(&b))
        {
            return Err(bad("invalid concrete guide query or total work bound"));
        }
        if let Some(i) = &self.intended {
            role(&i.reference_role)?;
            if i.reference_role != self.reference_role
                || i.manifest_hash != m.manifest_hash()?
                || i.content_revision != m.content_revision
                || i.guide_positions.len() != self.profile.guide_length
                || i.pam_positions.len() != self.profile.pam.len()
            {
                return Err(bad(
                    "intended tuple must bind exact reference role, manifest, content and profile lengths",
                ));
            }
            let c = contig(m, &i.contig_id)?;
            if i.guide_positions
                .iter()
                .chain(&i.pam_positions)
                .any(|p| *p >= c.length)
            {
                return Err(bad("intended position outside contig"));
            }
        } // A bounded ordinary window also verifies current built-in profile semantics.
        let dummy = Molecule::dna(
            "profile preflight",
            &"A".repeat(footprint(&self.profile)),
            Topology::Linear,
        )?;
        crispr::guide_preflight(
            &dummy,
            &GuideOptions {
                profile: self.profile.clone(),
                ..Default::default()
            },
        )?;
        Ok(())
    }
}
impl GenomeSearchSession {
    pub fn begin(
        m: &GenomeManifest,
        q: &GuideQuery,
        o: &GenomeSearchOptions,
    ) -> Result<Self, Error> {
        o.validate(m, q)?;
        Ok(Self {
            manifest: m.clone(),
            verifier: BlockVerifier::begin(m)?,
            semantic: Some(SemanticValidator::begin(m)?),
            query: q.clone(),
            options: o.clone(),
            manifest_hash: m.manifest_hash()?,
            ci: 0,
            bi: 0,
            blocks_processed: 0,
            bases_verified: 0,
            comparisons: 0,
            hits: vec![],
            finished: false,
            cancelled: false,
            source_integrity_verified: false,
            limits: vec![],
        })
    }
    pub fn progress(&self) -> GenomeSearchProgress {
        let exhaustive = self.finished
            && !self.cancelled
            && self.limits.is_empty()
            && self.source_integrity_verified;
        let next = if self.finished {
            vec![]
        } else {
            let c = &self.manifest.contigs[self.ci];
            let b = &c.blocks[self.bi];
            requests(
                c,
                b.start,
                window_count(c, b.start, u64::from(b.count), &self.options.profile),
            )
        };
        GenomeSearchProgress {
            engine_version: ENGINE_VERSION.into(),
            manifest_hash: self.manifest_hash.clone(),
            content_revision: self.manifest.content_revision.clone(),
            query: self.query.clone(),
            options: self.options.clone(),
            blocks_processed: self.blocks_processed,
            blocks_total: self
                .manifest
                .contigs
                .iter()
                .map(|c| c.blocks.len() as u64)
                .sum(),
            bases_verified: self.bases_verified,
            comparisons: self.comparisons,
            hits_retained: self.hits.len(),
            finished: self.finished,
            cancelled: self.cancelled,
            complete: exhaustive && self.manifest.expected_complete,
            exhaustive_for_manifest: exhaustive,
            expected_components_complete: self.manifest.expected_complete,
            source_integrity_verified: self.source_integrity_verified,
            next_blocks: next,
            limits: self.limits.clone(),
            warnings: warnings(),
        }
    }
    pub fn page(&self, offset: usize, limit: usize) -> Result<Vec<GenomeOffTargetHit>, Error> {
        if limit == 0 || limit > 256 || offset > self.hits.len() {
            return Err(bad("hit page requires a valid offset and 1..256 rows"));
        }
        Ok(self.hits[offset..self.hits.len().min(offset.saturating_add(limit))].to_vec())
    }
    pub fn cancel(&mut self) -> GenomeSearchProgress {
        if !self.finished {
            self.cancelled = true;
            self.finished = true;
            self.semantic = None;
            self.limits
                .push("Search cancelled; unsearched blocks remain.".into())
        }
        self.progress()
    }
    /// A failed block/semantic operation poisons the live cursor; it cannot later claim completeness.
    pub fn advance(&mut self, blocks: &[EncodedBlock]) -> Result<GenomeSearchProgress, Error> {
        if self.finished {
            return Err(bad("search is already finished"));
        }
        match self.advance_inner(blocks) {
            Ok(()) => Ok(self.progress()),
            Err(e) => {
                self.finished = true;
                self.semantic = None;
                self.limits.push(
                    "Block or canonical DNA verification failed; no complete coverage evidence."
                        .into(),
                );
                Err(e)
            }
        }
    }
    fn advance_inner(&mut self, blocks: &[EncodedBlock]) -> Result<(), Error> {
        let c = &self.manifest.contigs[self.ci];
        let b = &c.blocks[self.bi];
        let count = window_count(c, b.start, u64::from(b.count), &self.options.profile);
        let dna = verified_window(&self.verifier, c, b.start, count, blocks)?;
        let owner = self.verifier.verify(&c.id, b.index, &blocks[0].encoded)?;
        self.semantic
            .as_mut()
            .ok_or_else(|| bad("missing private semantic validator"))?
            .push(&owner)?;
        self.bases_verified += u64::from(b.count);
        let total = footprint(&self.options.profile);
        if c.length >= total as u64 && dna.len() >= total {
            let remaining = self
                .options
                .max_total_comparisons
                .saturating_sub(self.comparisons);
            if remaining == 0 {
                self.finished = true;
                self.limits
                    .push("Total search comparison budget reached.".into());
                return Ok(());
            }
            let local = molecule(&dna)?;
            let opts = crispr::OffTargetOptions {
                profile: self.options.profile.clone(),
                max_mismatches: self.options.max_mismatches,
                intended: None,
                reference_scope: "Bounded verified owner plus right halo".into(),
                max_comparisons: remaining.min(500_000_000) as usize,
                max_index_entries: self.options.max_index_entries_per_block,
                max_hits: 4096,
            };
            let refs = [crispr::ReferenceInput {
                id: c.id.clone(),
                molecule: local,
                build: None,
                accession: c.accession.clone(),
            }];
            let r = crispr::search_off_targets(&self.query, &refs, &opts)?;
            self.comparisons += r.comparisons as u64;
            for h in r.hits {
                let f = *h
                    .guide_positions
                    .iter()
                    .chain(&h.pam_positions)
                    .min()
                    .ok_or_else(|| bad("empty local hit"))?;
                if f >= b.count as usize {
                    continue;
                }
                if self.hits.len() == self.options.max_hits {
                    self.limits
                        .push("Hit retention limit reached; later loci remain unsearched.".into());
                    self.finished = true;
                    break;
                }
                let g = global_positions(&h.guide_positions, b.start, c.length);
                let p = global_positions(&h.pam_positions, b.start, c.length);
                let is_intended = self.options.intended.as_ref().is_some_and(|i| {
                    i.reference_role == self.options.reference_role
                        && i.manifest_hash == self.manifest_hash
                        && i.content_revision == self.manifest.content_revision
                        && i.contig_id == c.id
                        && i.guide_positions == g
                        && i.pam_positions == p
                        && i.reverse == h.reverse
                });
                let (top, bottom) = cuts(&self.options.profile, b.start + f as u64, h.reverse, c);
                self.hits.push(GenomeOffTargetHit {
                    id: hash(&(
                        global_id(
                            &self.manifest_hash,
                            &self.options.reference_role,
                            &c.id,
                            &self.options.profile,
                            &g,
                            &p,
                            h.reverse,
                        )?,
                        &self.query,
                    ))?,
                    contig_id: c.id.clone(),
                    reference_guide: h.reference_guide,
                    pam_sequence: h.pam_sequence,
                    guide_positions: g,
                    pam_positions: p,
                    reverse: h.reverse,
                    min_mismatches: h.min_mismatches,
                    max_mismatches: h.max_mismatches,
                    pam_certain: h.pam_certain,
                    definite: h.definite,
                    is_intended,
                    nominal_top_cut: top,
                    nominal_bottom_cut: bottom,
                });
            }
            if !r.complete {
                self.finished = true;
                self.limits.extend(
                    r.limits
                        .into_iter()
                        .map(|s| format!("Owner-block analysis incomplete: {s}")),
                );
            }
        }
        self.blocks_processed += 1;
        self.bi += 1;
        if self.bi == c.blocks.len() {
            self.ci += 1;
            self.bi = 0;
        }
        if self.ci == self.manifest.contigs.len() && !self.finished {
            let receipt = self
                .semantic
                .take()
                .ok_or_else(|| bad("missing semantic EOF validator"))?
                .finish()?;
            if receipt.manifest_hash != self.manifest_hash
                || receipt.content_revision != self.manifest.content_revision
            {
                return Err(bad("semantic receipt identity differs from live cursor"));
            }
            self.source_integrity_verified = true;
            self.finished = true;
        }
        Ok(())
    }
}
dto! {pub struct GenomeIntervalReadOptions{pub source_role:String,pub interval:GenomicInterval}}
dto! {pub struct GenomeIntervalRead{pub manifest_hash:String,pub content_revision:String,pub options:GenomeIntervalReadOptions,pub sequence:String,pub warnings:Vec<String>}}
impl GenomeIntervalReadOptions {
    pub fn validate(&self, m: &GenomeManifest) -> Result<(), Error> {
        GenomeGuideOptions {
            source_role: self.source_role.clone(),
            interval: self.interval.clone(),
            guide: GuideOptions::default(),
        }
        .validate(m)
    }
}
pub fn read_interval_required_blocks(
    m: &GenomeManifest,
    o: &GenomeIntervalReadOptions,
) -> Result<Vec<BlockRequest>, Error> {
    o.validate(m)?;
    Ok(requests(
        contig(m, &o.interval.contig_id)?,
        o.interval.start,
        o.interval.end - o.interval.start,
    ))
}
pub fn read_interval(
    m: &GenomeManifest,
    o: &GenomeIntervalReadOptions,
    blocks: &[EncodedBlock],
) -> Result<GenomeIntervalRead, Error> {
    o.validate(m)?;
    let c = contig(m, &o.interval.contig_id)?;
    let dna = verified_window(
        &BlockVerifier::begin(m)?,
        c,
        o.interval.start,
        o.interval.end - o.interval.start,
        blocks,
    )?;
    Ok(GenomeIntervalRead{manifest_hash:m.manifest_hash()?,content_revision:m.content_revision.clone(),options:o.clone(),sequence:String::from_utf8(dna).map_err(|_|bad("invalid canonical DNA"))?,warnings:vec!["Local encoded-block integrity only; whole-contig canonical sequence hash has not been verified by this interval read.".into()]})
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn strings(v: &[String]) -> bool {
    v.len() <= 128
        && v.iter()
            .all(|s| !s.is_empty() && s.len() <= 4096 && !s.chars().any(char::is_control))
}
impl GenomeGuideReport {
    /// Shape and request binding only. Never grants live scientific verification authority.
    pub fn validate(&self) -> Result<(), Error> {
        role(&self.options.source_role)?;
        self.options.guide.validate()?;
        if self.engine_version != ENGINE_VERSION
            || !digest(&self.manifest_hash)
            || !digest(&self.content_revision)
            || self.options_revision != hash(&self.options)?
            || self.candidates.len() > self.options.guide.max_candidates
            || self.comparisons > self.options.guide.max_comparisons as u64
            || self.complete != self.limits.is_empty()
            || !strings(&self.limits)
            || !strings(&self.warnings)
        {
            return Err(bad("invalid interval guide report shape or limit state"));
        }
        let p = &self.options.guide.profile;
        let mut ids = std::collections::BTreeSet::new();
        for h in &self.candidates {
            if !digest(&h.id)
                || !ids.insert(&h.id)
                || h.contig_id != self.options.interval.contig_id
                || h.guide_sequence.len() != p.guide_length
                || h.pam_sequence.len() != p.pam.len()
                || h.guide_positions.len() != p.guide_length
                || h.pam_positions.len() != p.pam.len()
                || h.gc_bp.is_some_and(|n| n > 10000)
                || h.longest_homopolymer.is_some_and(|n| n > p.guide_length)
                || h.passes_filters != h.rejections.is_empty()
                || !strings(&h.rejections)
                || !h
                    .guide_sequence
                    .bytes()
                    .chain(h.pam_sequence.bytes())
                    .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
            {
                return Err(bad("invalid global guide shape"));
            }
        }
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        m: &GenomeManifest,
        o: &GenomeGuideOptions,
    ) -> Result<(), Error> {
        self.validate()?;
        o.validate(m)?;
        if self.options != *o
            || self.manifest_hash != m.manifest_hash()?
            || self.content_revision != m.content_revision
        {
            return Err(bad(
                "interval report is bound to another request or manifest",
            ));
        }
        let c = contig(m, &o.interval.contig_id)?;
        for h in &self.candidates {
            if h.guide_positions
                .iter()
                .chain(&h.pam_positions)
                .any(|p| *p >= c.length)
                || h.nominal_top_cut
                    .into_iter()
                    .chain(h.nominal_bottom_cut)
                    .any(|p| {
                        if c.topology == Topology::Circular {
                            p >= c.length
                        } else {
                            p > c.length
                        }
                    })
                || h.id
                    != global_id(
                        &self.manifest_hash,
                        &o.source_role,
                        &c.id,
                        &o.guide.profile,
                        &h.guide_positions,
                        &h.pam_positions,
                        h.reverse,
                    )?
            {
                return Err(bad(
                    "guide coordinates or stable identity differ from bound contig",
                ));
            }
        }
        Ok(())
    }
}
impl GenomeIntervalRead {
    pub fn validate(&self) -> Result<(), Error> {
        role(&self.options.source_role)?;
        let n = self
            .options
            .interval
            .end
            .checked_sub(self.options.interval.start)
            .filter(|n| *n > 0 && *n <= MAX_INTERVAL_BASES)
            .ok_or_else(|| bad("invalid interval read length"))?;
        if !digest(&self.manifest_hash)
            || !digest(&self.content_revision)
            || self.sequence.len() as u64 != n
            || !self
                .sequence
                .bytes()
                .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
            || !strings(&self.warnings)
        {
            return Err(bad("invalid interval read shape"));
        }
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        m: &GenomeManifest,
        o: &GenomeIntervalReadOptions,
    ) -> Result<(), Error> {
        self.validate()?;
        o.validate(m)?;
        if self.options != *o
            || self.manifest_hash != m.manifest_hash()?
            || self.content_revision != m.content_revision
        {
            return Err(bad("interval read is bound to another request or manifest"));
        }
        Ok(())
    }
}
impl GenomeSearchProgress {
    pub fn validate(&self) -> Result<(), Error> {
        role(&self.options.reference_role)?;
        self.options.profile.validate()?;
        let exhaustive = self.finished
            && !self.cancelled
            && self.limits.is_empty()
            && self.source_integrity_verified;
        if self.engine_version != ENGINE_VERSION
            || !digest(&self.manifest_hash)
            || !digest(&self.content_revision)
            || self.query.id.is_empty()
            || self.query.id.len() > 256
            || self.query.sequence.len() != self.options.profile.guide_length
            || !self.query.sequence.bytes().all(|b| b"ACGT".contains(&b))
            || self.options.max_mismatches > 6
            || self.options.max_hits == 0
            || self.options.max_hits > 4096
            || self.options.max_index_entries_per_block == 0
            || self.options.max_index_entries_per_block > 2_000_000
            || self.options.max_total_comparisons == 0
            || self.options.max_total_comparisons > 10_000_000_000_000
            || self.comparisons > self.options.max_total_comparisons
            || self.hits_retained > self.options.max_hits
            || self.blocks_processed > self.blocks_total
            || self.blocks_total == 0
            || self.blocks_total > 250000
            || self.next_blocks.len() > 3
            || self.finished != self.next_blocks.is_empty()
            || self.cancelled && (!self.finished || self.limits.is_empty())
            || self.exhaustive_for_manifest != exhaustive
            || self.complete != (exhaustive && self.expected_components_complete)
            || self.source_integrity_verified
                && (!self.finished || self.blocks_processed != self.blocks_total)
            || !strings(&self.limits)
            || !strings(&self.warnings)
        {
            return Err(bad("invalid search progress counters or coverage state"));
        }
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        m: &GenomeManifest,
        q: &GuideQuery,
        o: &GenomeSearchOptions,
    ) -> Result<(), Error> {
        self.validate()?;
        o.validate(m, q)?;
        if self.manifest_hash != m.manifest_hash()?
            || self.content_revision != m.content_revision
            || self.query != *q
            || self.options != *o
            || self.blocks_total != m.contigs.iter().map(|c| c.blocks.len() as u64).sum::<u64>()
            || self.bases_verified > m.total_bases
            || self.expected_components_complete != m.expected_complete
        {
            return Err(bad("search progress bound to another manifest or request"));
        }
        for r in &self.next_blocks {
            let c = contig(m, &r.contig_id)?;
            if r.index >= c.blocks.len() as u64 {
                return Err(bad("next block outside manifest"));
            }
        }
        if self.source_integrity_verified && self.bases_verified != m.total_bases {
            return Err(bad("semantic coverage count differs from manifest"));
        }
        Ok(())
    }
}
