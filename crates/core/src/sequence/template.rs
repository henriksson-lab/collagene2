//! Literal homology-donor design with an independent source/target/candidate checker.
//! Sequence consistency is not a prediction of cellular HDR or experimental editing.
use crate::sequence::{
    annotations::{self, MappingSegment},
    crispr::{self, GuideCandidate, GuideOptions},
};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const ENGINE_VERSION: &str = "edit-template-1";
const MAX_SOURCE: usize = 1_000_000;
const MAX_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_DONOR: usize = 100_000;
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
dto! {pub struct GuideBinding{pub options:GuideOptions,pub candidate:GuideCandidate}}
dto! {pub struct EditTemplateOptions{pub edit_start:usize,pub delete_length:usize,pub replacement:String,pub left_arm:usize,pub right_arm:usize,pub reverse_donor:bool,pub max_cut_distance:usize}}
impl Default for EditTemplateOptions {
    fn default() -> Self {
        Self {
            edit_start: 0,
            delete_length: 0,
            replacement: "A".into(),
            left_arm: 30,
            right_arm: 30,
            reverse_donor: false,
            max_cut_distance: 100,
        }
    }
}
dto! {pub struct HomologyArm{pub sequence:String,pub source_positions:Vec<usize>}}
dto! {pub struct EditTemplateCheck{pub engine_version:String,pub source_revision:String,pub target_revision:String,pub options_revision:String,pub candidate_revision:String,pub batch_revision:String,pub donor_length:usize,pub left_matches:bool,pub replacement_matches:bool,pub right_matches:bool,pub target_matches:bool,pub changes_sequence:bool,pub guide_cut_distance:Option<usize>,pub guide_within_distance:Option<bool>,pub consistent:bool,pub warnings:Vec<String>}}
dto! {pub struct EditTemplateReport{pub engine_version:String,pub source_revision:String,pub target_revision:String,pub options_revision:String,pub batch_revision:String,pub options:EditTemplateOptions,pub guide:Option<GuideBinding>,pub left:HomologyArm,pub right:HomologyArm,pub deleted_positions:Vec<usize>,pub donor:Molecule,pub edited_reference:Molecule,pub check:EditTemplateCheck,pub warnings:Vec<String>}}
fn bad(text: &str) -> Error {
    Error::Unsupported(format!("edit template: {text}"))
}
struct Bound(usize);
impl Write for Bound {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > MAX_BYTES.saturating_sub(self.0) {
            return Err(std::io::Error::other("bounded edit-template evidence"));
        }
        self.0 += b.len();
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded(v: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(Bound(0), v).map_err(|_| Error::TooLarge)
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    bounded(v)?;
    Ok(
        blake3::hash(&serde_json::to_vec(v).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn location(l: &Location, depth: usize, nodes: &mut usize) -> Result<(), Error> {
    *nodes = nodes.saturating_add(1);
    if depth > 32 || *nodes > 100_000 {
        return Err(Error::TooLarge);
    }
    match l {
        Location::Complement(x) => location(x, depth + 1, nodes)?,
        Location::Join(xs) | Location::Order(xs) => {
            for x in xs {
                location(x, depth + 1, nodes)?
            }
        }
        _ => (),
    };
    Ok(())
}
fn molecule_preflight(m: &Molecule, max: usize) -> Result<(), Error> {
    if m.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("source and intended target must be DNA"));
    }
    if m.sequence.len() > max
        || m.features.len() > 4096
        || m.metadata.len() > 4096
        || m.name.len() > 4096
    {
        return Err(Error::TooLarge);
    }
    let mut nodes = 1 + m.features.len() + m.metadata.len();
    for f in &m.features {
        nodes = nodes.saturating_add(f.qualifiers.len());
        location(&f.location, 0, &mut nodes)?;
    }
    bounded(m)?;
    m.validate()?;
    if m.sequence.as_str().bytes().any(|b| !b"ACGT".contains(&b)) {
        return Err(bad(
            "literal homology design requires concrete source and intended-target DNA; unresolved bases are not synthesized by inference",
        ));
    }
    Ok(())
}
fn replacement(o: &EditTemplateOptions) -> Result<String, Error> {
    if o.replacement.is_empty() {
        Ok(String::new())
    } else {
        let s = Sequence::dna(&o.replacement)?;
        if s.as_str().bytes().any(|b| !b"ACGT".contains(&b)) {
            return Err(bad("replacement must be concrete DNA"));
        }
        Ok(s.as_str().into())
    }
}
impl EditTemplateOptions {
    pub fn validate(&self) -> Result<(), Error> {
        if self.edit_start > MAX_SOURCE
            || self.delete_length > MAX_SOURCE
            || self.replacement.len() > 2 * MAX_DONOR
            || self.left_arm == 0
            || self.right_arm == 0
            || self.left_arm > MAX_DONOR
            || self.right_arm > MAX_DONOR
            || self.max_cut_distance > MAX_SOURCE
        {
            return Err(bad("invalid edit/homology/distance bounds"));
        }
        let r = replacement(self)?;
        if self.delete_length == 0 && r.is_empty() {
            return Err(bad("empty edit has no intended change"));
        }
        if self
            .left_arm
            .checked_add(r.len())
            .and_then(|n| n.checked_add(self.right_arm))
            .is_none_or(|n| n > MAX_DONOR)
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
pub fn source_preflight(
    source: &Molecule,
    target: &Molecule,
    o: &EditTemplateOptions,
    guide: Option<&GuideBinding>,
) -> Result<(), Error> {
    o.validate()?;
    molecule_preflight(source, MAX_SOURCE)?;
    molecule_preflight(target, MAX_SOURCE + MAX_DONOR)?;
    let n = source.sequence.len();
    if o.delete_length >= n
        || o.left_arm
            .checked_add(o.delete_length)
            .and_then(|x| x.checked_add(o.right_arm))
            .is_none_or(|x| x > n)
    {
        return Err(bad(
            "homology arms must be distinct, outside the edit, and traverse the original source at most once",
        ));
    }
    match source.topology {
        Topology::Linear => {
            if o.edit_start > n
                || o.edit_start < o.left_arm
                || o.edit_start
                    .checked_add(o.delete_length)
                    .and_then(|x| x.checked_add(o.right_arm))
                    .is_none_or(|x| x > n)
            {
                return Err(bad("linear edit lacks the requested homology flanks"));
            }
        }
        Topology::Circular => {
            if o.edit_start >= n {
                return Err(bad(
                    "circular edit start must be an original source boundary",
                ));
            }
        }
    }
    if let Some(g) = guide {
        bounded(g)?;
        crispr::validate_candidate(source, &g.options, &g.candidate)?;
    }
    Ok(())
}
pub fn options_revision(o: &EditTemplateOptions) -> Result<String, Error> {
    o.validate()?;
    hash(o)
}
pub fn batch_revision(
    source: &Molecule,
    target: &Molecule,
    o: &EditTemplateOptions,
    guide: Option<&GuideBinding>,
) -> Result<String, Error> {
    source_preflight(source, target, o, guide)?;
    hash(&(
        ENGINE_VERSION,
        source.revision(),
        target.revision(),
        options_revision(o)?,
        guide,
    ))
}
fn candidate_sequence(s: &str) -> Result<Sequence, Error> {
    if s.len() > 2 * MAX_DONOR {
        return Err(Error::TooLarge);
    }
    let sequence = Sequence::dna(s)?;
    if sequence.len() > MAX_DONOR || sequence.as_str().bytes().any(|b| !b"ACGT".contains(&b)) {
        return Err(bad(
            "supplied donor must be concrete DNA within the donor bound",
        ));
    }
    Ok(sequence)
}
pub fn check_batch_revision(
    source: &Molecule,
    target: &Molecule,
    o: &EditTemplateOptions,
    guide: Option<&GuideBinding>,
    donor: &str,
) -> Result<String, Error> {
    let seq = candidate_sequence(donor)?;
    hash(&(batch_revision(source, target, o, guide)?, seq.as_str()))
}
fn same_reference(a: &str, b: &str, topology: Topology) -> bool {
    a.len() == b.len()
        && if topology == Topology::Circular {
            crate::sequence::analysis::rotation_offset(a, b).is_some()
        } else {
            a == b
        }
}
fn cut_distance(
    source: &Molecule,
    o: &EditTemplateOptions,
    guide: Option<&GuideBinding>,
) -> (Option<usize>, Option<bool>) {
    let Some(g) = guide else { return (None, None) };
    let Some(cut) = g.candidate.nominal_top_cut else {
        return (None, Some(false));
    };
    let n = source.sequence.len();
    let distance = if source.topology == Topology::Circular {
        let d = (cut + n - o.edit_start) % n;
        if d <= o.delete_length {
            0
        } else {
            let end = (o.edit_start + o.delete_length) % n;
            let near = |a: usize, b: usize| {
                let d = a.abs_diff(b);
                d.min(n - d)
            };
            near(cut, o.edit_start).min(near(cut, end))
        }
    } else {
        let end = o.edit_start + o.delete_length;
        if cut < o.edit_start {
            o.edit_start - cut
        } else {
            cut.saturating_sub(end)
        }
    };
    (Some(distance), Some(distance <= o.max_cut_distance))
}
/// Independent checker: inspect supplied donor bases and reconstruct edited source
/// by an original-coordinate walk, without invoking the designer or its arm builder.
pub fn check(
    source: &Molecule,
    target: &Molecule,
    o: &EditTemplateOptions,
    guide: Option<&GuideBinding>,
    donor: &str,
) -> Result<EditTemplateCheck, Error> {
    source_preflight(source, target, o, guide)?;
    let candidate = candidate_sequence(donor)?;
    let oriented = if o.reverse_donor {
        candidate.reverse_complement()?
    } else {
        candidate.clone()
    };
    let supplied = oriented.as_str().as_bytes();
    let old = source.sequence.as_str().as_bytes();
    let n = old.len();
    let insert = replacement(o)?;
    let expected_len = o.left_arm + insert.len() + o.right_arm;
    let length_ok = supplied.len() == expected_len;
    let left_start = if source.topology == Topology::Circular {
        (o.edit_start + n - o.left_arm) % n
    } else {
        o.edit_start - o.left_arm
    };
    let right_start = (o.edit_start + o.delete_length) % n;
    let left_matches =
        length_ok && (0..o.left_arm).all(|j| supplied[j] == old[(left_start + j) % n]);
    let replacement_matches =
        length_ok && (0..insert.len()).all(|j| supplied[o.left_arm + j] == insert.as_bytes()[j]);
    let right_matches = length_ok
        && (0..o.right_arm)
            .all(|j| supplied[o.left_arm + insert.len() + j] == old[(right_start + j) % n]);
    let wraps = source.topology == Topology::Circular && o.edit_start + o.delete_length > n;
    let mut edited = Vec::with_capacity(n - o.delete_length + insert.len());
    if wraps {
        edited.extend_from_slice(insert.as_bytes());
        let mut i = right_start;
        for _ in 0..n - o.delete_length {
            edited.push(old[i]);
            i = (i + 1) % n;
        }
        if i != o.edit_start {
            return Err(bad("circular independent coordinate walk failed"));
        }
    } else {
        for (i, b) in old.iter().enumerate() {
            if i == o.edit_start {
                edited.extend_from_slice(insert.as_bytes());
            }
            if i < o.edit_start || i >= o.edit_start + o.delete_length {
                edited.push(*b);
            }
        }
        if o.edit_start == n {
            edited.extend_from_slice(insert.as_bytes());
        }
    }
    let edited = String::from_utf8(edited).map_err(|_| bad("non-DNA coordinate walk"))?;
    let target_matches = target.topology == source.topology
        && same_reference(&edited, target.sequence.as_str(), source.topology);
    let changes_sequence = !same_reference(&edited, source.sequence.as_str(), source.topology);
    let (distance, within) = cut_distance(source, o, guide);
    let mut warnings=vec!["Sequence/homology consistency only: cellular HDR, guide activity, repair outcome, off-target risk, donor strandedness and laboratory conditions are not established.".into(),"The intended target is an independent requirement; its bases are never copied into the generated edit or homology donor.".into()];
    if guide.is_none() {
        warnings.push(
            "No guide is bound; this donor check establishes no cleavage-coordinate association."
                .into(),
        );
    }
    if within == Some(false) {
        warnings.push("The bound guide's nominal top-strand cut is unavailable or outside the declared edit-distance bound; this is not a validated editing design.".into());
    }
    if o.left_arm < 20 || o.right_arm < 20 {
        warnings.push("Short homology arms are retained as literal sequence parameters; no empirical repair-efficiency claim is made.".into());
    }
    let r = EditTemplateCheck {
        engine_version: ENGINE_VERSION.into(),
        source_revision: source.revision(),
        target_revision: target.revision(),
        options_revision: options_revision(o)?,
        candidate_revision: hash(&candidate)?,
        batch_revision: check_batch_revision(source, target, o, guide, donor)?,
        donor_length: candidate.len(),
        left_matches,
        replacement_matches,
        right_matches,
        target_matches,
        changes_sequence,
        guide_cut_distance: distance,
        guide_within_distance: within,
        consistent: left_matches
            && replacement_matches
            && right_matches
            && target_matches
            && changes_sequence
            && within != Some(false),
        warnings,
    };
    r.validate()?;
    Ok(r)
}
fn arm(source: &Molecule, start: usize, len: usize) -> HomologyArm {
    let n = source.sequence.len();
    let source_positions = (0..len).map(|i| (start + i) % n).collect::<Vec<_>>();
    HomologyArm {
        sequence: source_positions
            .iter()
            .map(|i| source.sequence.as_str().as_bytes()[*i] as char)
            .collect(),
        source_positions,
    }
}
fn maps(positions: &[usize], offset: usize) -> Vec<MappingSegment> {
    let mut result = vec![];
    let mut i = 0;
    while i < positions.len() {
        let start = i;
        i += 1;
        while i < positions.len() && positions[i] == positions[i - 1] + 1 {
            i += 1;
        }
        result.push(MappingSegment {
            source_start: positions[start],
            source_end: positions[i - 1] + 1,
            destination_start: offset + start,
            reverse: false,
        });
    }
    result
}
fn note_feature(label: &str, start: usize, end: usize, reverse: bool) -> Option<Feature> {
    if start >= end {
        return None;
    }
    let range = Location::range(start, end);
    Some(Feature {
        kind: "misc_feature".into(),
        location: if reverse {
            Location::Complement(Box::new(range))
        } else {
            range
        },
        qualifiers: vec![
            ("label".into(), Some(label.into())),
            (
                "note".into(),
                Some(
                    "Sequence-design component; experimental recombination is not established."
                        .into(),
                ),
            ),
        ],
        needs_review: true,
    })
}
fn project(
    source: &Molecule,
    destination: &mut Molecule,
    mapping: &[MappingSegment],
    warnings: &mut Vec<String>,
) -> Result<(), Error> {
    let p = annotations::project_features_with_topology(
        source,
        mapping,
        destination.sequence.len(),
        destination.topology,
    )?;
    destination.features = p.features;
    warnings.extend(p.warnings.into_iter().take(64));
    if !p.unmapped_features.is_empty() {
        warnings.push(format!("{} source annotations remain on the original source graph rather than being asserted on the changed sequence.",p.unmapped_features.len()));
    }
    Ok(())
}
fn metadata(
    m: &mut Molecule,
    role: &str,
    source: &Molecule,
    target: &Molecule,
    o: &EditTemplateOptions,
    guide: Option<&GuideBinding>,
    batch: &str,
) -> Result<(), Error> {
    for(key,value)in [("engine",ENGINE_VERSION.into()),("role",role.into()),("source_revision",source.revision()),("target_revision",target.revision()),("options_revision",options_revision(o)?),("design_revision",batch.into()),("options",serde_json::to_string(o).map_err(|_|Error::TooLarge)?),("guide",serde_json::to_string(&guide).map_err(|_|Error::TooLarge)?),("authority","Unverified producer sequence-design evidence; reconstruct from the exact source, target requirement and parameters. No physical or cellular outcome is established.".into())]{m.metadata.insert(format!("edit_template:{key}"),value);}
    Ok(())
}
pub fn design(
    source: &Molecule,
    target: &Molecule,
    o: &EditTemplateOptions,
    guide: Option<&GuideBinding>,
) -> Result<EditTemplateReport, Error> {
    source_preflight(source, target, o, guide)?;
    let n = source.sequence.len();
    let insert = replacement(o)?;
    let left_start = if source.topology == Topology::Circular {
        (o.edit_start + n - o.left_arm) % n
    } else {
        o.edit_start - o.left_arm
    };
    let right_start = (o.edit_start + o.delete_length) % n;
    let left = arm(source, left_start, o.left_arm);
    let right = arm(source, right_start, o.right_arm);
    let deleted_positions = (0..o.delete_length)
        .map(|i| (o.edit_start + i) % n)
        .collect();
    let forward = format!("{}{}{}", left.sequence, insert, right.sequence);
    let donor_seq = if o.reverse_donor {
        Sequence::dna(&forward)?
            .reverse_complement()?
            .as_str()
            .to_string()
    } else {
        forward
    };
    let mut donor = Molecule::dna("homology_donor", &donor_seq, Topology::Linear)?;
    let mut donor_maps = maps(&left.source_positions, 0);
    donor_maps.extend(maps(&right.source_positions, o.left_arm + insert.len()));
    if o.reverse_donor {
        for m in &mut donor_maps {
            m.destination_start =
                donor_seq.len() - m.destination_start - (m.source_end - m.source_start);
            m.reverse = true;
        }
    }
    let wraps = source.topology == Topology::Circular && o.edit_start + o.delete_length > n;
    let mut edited_maps = vec![];
    let edited = if wraps {
        let right = (o.edit_start + o.delete_length) % n;
        if right < o.edit_start {
            edited_maps.push(MappingSegment {
                source_start: right,
                source_end: o.edit_start,
                destination_start: insert.len(),
                reverse: false,
            });
        }
        format!(
            "{}{}",
            insert,
            &source.sequence.as_str()[right..o.edit_start]
        )
    } else {
        let end = o.edit_start + o.delete_length;
        if o.edit_start > 0 {
            edited_maps.push(MappingSegment {
                source_start: 0,
                source_end: o.edit_start,
                destination_start: 0,
                reverse: false,
            });
        }
        if end < n {
            edited_maps.push(MappingSegment {
                source_start: end,
                source_end: n,
                destination_start: o.edit_start + insert.len(),
                reverse: false,
            });
        }
        format!(
            "{}{}{}",
            &source.sequence.as_str()[..o.edit_start],
            insert,
            &source.sequence.as_str()[end..]
        )
    };
    let mut edited_reference =
        Molecule::dna("intended_edited_reference", &edited, source.topology)?;
    let mut warnings=vec!["Homology arms refer to distinct original source positions and are not predictions of HDR efficiency, genome specificity or experimental success.".into(),"Original raw quality/physical evidence stays with source revisions; changed donor/reference DNA does not inherit it.".into()];
    if wraps {
        warnings.push("The edit removes the original circular origin; the edited reference starts at the replacement, followed by retained source bases.".into());
    }
    project(source, &mut donor, &donor_maps, &mut warnings)?;
    project(source, &mut edited_reference, &edited_maps, &mut warnings)?;
    let orient = |start: usize, end: usize| {
        if o.reverse_donor {
            (donor_seq.len() - end, donor_seq.len() - start)
        } else {
            (start, end)
        }
    };
    for (label, a, b) in [
        ("Left homology arm", 0, o.left_arm),
        (
            "Intended replacement",
            o.left_arm,
            o.left_arm + insert.len(),
        ),
        (
            "Right homology arm",
            o.left_arm + insert.len(),
            donor_seq.len(),
        ),
    ] {
        let (a, b) = orient(a, b);
        if let Some(f) = note_feature(label, a, b, o.reverse_donor) {
            donor.features.push(f);
        }
    }
    let at = if wraps { 0 } else { o.edit_start };
    if let Some(f) = note_feature("Intended replacement", at, at + insert.len(), false) {
        edited_reference.features.push(f);
    }
    let batch = batch_revision(source, target, o, guide)?;
    metadata(&mut donor, "donor", source, target, o, guide, &batch)?;
    metadata(
        &mut edited_reference,
        "edited_reference",
        source,
        target,
        o,
        guide,
        &batch,
    )?;
    let result = EditTemplateReport {
        engine_version: ENGINE_VERSION.into(),
        source_revision: source.revision(),
        target_revision: target.revision(),
        options_revision: options_revision(o)?,
        batch_revision: batch,
        options: o.clone(),
        guide: guide.cloned(),
        left,
        right,
        deleted_positions,
        check: check(source, target, o, guide, &donor_seq)?,
        donor,
        edited_reference,
        warnings,
    };
    result.validate()?;
    Ok(result)
}
impl EditTemplateCheck {
    /// Historical shape only. Imported pass flags do not confer current authority.
    pub fn validate(&self) -> Result<(), Error> {
        if self.engine_version.is_empty()
            || self.engine_version.len() > 128
            || [
                &self.source_revision,
                &self.target_revision,
                &self.options_revision,
                &self.candidate_revision,
                &self.batch_revision,
            ]
            .iter()
            .any(|s| !hex(s))
            || self.donor_length == 0
            || self.donor_length > MAX_DONOR
            || self.guide_cut_distance.is_some_and(|d| d > MAX_SOURCE)
            || self.warnings.len() > 64
            || self.warnings.iter().any(|s| s.len() > 4096)
            || self.consistent
                != (self.left_matches
                    && self.replacement_matches
                    && self.right_matches
                    && self.target_matches
                    && self.changes_sequence
                    && self.guide_within_distance != Some(false))
        {
            return Err(bad("invalid check structure"));
        }
        bounded(self)
    }
    pub fn validate_bindings_for(
        &self,
        source: &Molecule,
        target: &Molecule,
        o: &EditTemplateOptions,
        guide: Option<&GuideBinding>,
        donor: &str,
    ) -> Result<(), Error> {
        source_preflight(source, target, o, guide)?;
        self.validate()?;
        let seq = candidate_sequence(donor)?;
        if self.engine_version != ENGINE_VERSION
            || self.source_revision != source.revision()
            || self.target_revision != target.revision()
            || self.options_revision != options_revision(o)?
            || self.candidate_revision != hash(&seq)?
            || self.donor_length != seq.len()
            || self.batch_revision != check_batch_revision(source, target, o, guide, donor)?
            || self.guide_within_distance.is_some() != guide.is_some()
        {
            return Err(bad(
                "check source/target/options/candidate binding mismatch",
            ));
        }
        Ok(())
    }
}
impl EditTemplateReport {
    /// Historical structure and internally consistent producer claims, not a replay.
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        self.check.validate()?;
        if self.engine_version.is_empty()
            || self.engine_version.len() > 128
            || [
                &self.source_revision,
                &self.target_revision,
                &self.options_revision,
                &self.batch_revision,
            ]
            .iter()
            .any(|s| !hex(s))
            || self.options_revision != options_revision(&self.options)?
            || self.left.source_positions.len() != self.options.left_arm
            || self.right.source_positions.len() != self.options.right_arm
            || self.left.sequence.len() != self.options.left_arm
            || self.right.sequence.len() != self.options.right_arm
            || self.deleted_positions.len() != self.options.delete_length
            || self
                .left
                .source_positions
                .iter()
                .chain(&self.right.source_positions)
                .chain(&self.deleted_positions)
                .any(|p| *p >= MAX_SOURCE)
            || self.donor.topology != Topology::Linear
            || self.donor.sequence.alphabet() != Alphabet::Dna
            || self.donor.sequence.len() > MAX_DONOR
            || self.edited_reference.sequence.len() > MAX_SOURCE + MAX_DONOR
            || self.edited_reference.sequence.alphabet() != Alphabet::Dna
            || self.check.source_revision != self.source_revision
            || self.check.target_revision != self.target_revision
            || self.check.options_revision != self.options_revision
            || self.check.candidate_revision != hash(&self.donor.sequence)?
            || self.check.batch_revision
                != hash(&(self.batch_revision.as_str(), self.donor.sequence.as_str()))?
            || self.warnings.len() > 256
            || self.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("invalid bounded report structure"));
        }
        self.donor.validate()?;
        self.edited_reference.validate()?;
        bounded(self)
    }
    pub fn validate_bindings_for(
        &self,
        source: &Molecule,
        target: &Molecule,
        o: &EditTemplateOptions,
        guide: Option<&GuideBinding>,
    ) -> Result<(), Error> {
        source_preflight(source, target, o, guide)?;
        self.validate()?;
        if self.engine_version != ENGINE_VERSION
            || self.source_revision != source.revision()
            || self.target_revision != target.revision()
            || self.options != *o
            || self.guide.as_ref() != guide
            || self.batch_revision != batch_revision(source, target, o, guide)?
            || self.edited_reference.topology != source.topology
            || self
                .left
                .source_positions
                .iter()
                .chain(&self.right.source_positions)
                .chain(&self.deleted_positions)
                .any(|p| *p >= source.sequence.len())
        {
            return Err(bad("report source/target/options/guide binding mismatch"));
        }
        self.check
            .validate_bindings_for(source, target, o, guide, self.donor.sequence.as_str())
    }
}
/// Full computational replay; loading a producer report does not invoke this implicitly.
pub fn replay(
    source: &Molecule,
    target: &Molecule,
    o: &EditTemplateOptions,
    guide: Option<&GuideBinding>,
    report: &EditTemplateReport,
) -> Result<(), Error> {
    report.validate_bindings_for(source, target, o, guide)?;
    if design(source, target, o, guide)? != *report {
        return Err(bad(
            "producer report differs from fresh source-bound design",
        ));
    }
    Ok(())
}
