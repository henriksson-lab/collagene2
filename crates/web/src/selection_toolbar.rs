//! Contextual actions over the existing reference selection. No edits or source
//! copies happen here; the app owns annotation/primer dialogs and history.
use collagene_core::{Alphabet, Location, Molecule, Topology};
#[cfg(target_arch = "wasm32")]
use leptos::prelude::*;

const MAX_PARTS: usize = 4096;
const MAX_LOCATION_NODES: usize = 8192;
const MAX_GC_BASES: usize = 10_000;

#[derive(Clone, PartialEq, Eq)]
struct Summary {
    heading: String,
    detail: String,
    annotate: bool,
    primer: bool,
    edit: bool,
}
impl Summary {
    fn limited(heading: String, detail: &str) -> Self {
        Self {
            heading,
            detail: detail.into(),
            annotate: false,
            primer: false,
            edit: false,
        }
    }
}
#[derive(Default)]
struct LocationFlags {
    compound: bool,
    unresolved: bool,
    fuzzy: bool,
    ordered: bool,
}
fn location_flags(location: &Location) -> Option<LocationFlags> {
    let mut stack = vec![(location, 0usize)];
    let mut nodes = 0usize;
    let mut flags = LocationFlags::default();
    while let Some((location, depth)) = stack.pop() {
        nodes = nodes.checked_add(1)?;
        if nodes > MAX_LOCATION_NODES || depth >= 32 {
            return None;
        }
        match location {
            Location::Range { before, after, .. } => flags.fuzzy |= *before || *after,
            Location::Complement(inner) => stack.push((inner, depth + 1)),
            Location::Join(parts) | Location::Order(parts) => {
                flags.compound = true;
                flags.ordered |= matches!(location, Location::Order(_));
                if parts.len() > MAX_LOCATION_NODES.saturating_sub(nodes + stack.len()) {
                    return None;
                }
                stack.extend(parts.iter().rev().map(|part| (part, depth + 1)));
            }
            Location::Unresolved(_) => flags.unresolved = true,
        }
    }
    Some(flags)
}

/// Only one exact source-forward interval, or its two ordered circular-origin
/// pieces, can seed an oligo. Return a normalized start and unrolled end.
fn primer_geometry(
    molecule: &Molecule,
    range: Option<(usize, usize)>,
    feature: Option<usize>,
    parts: &[(usize, usize)],
) -> Option<(usize, usize)> {
    if molecule.sequence.alphabet() != Alphabet::Dna || parts.len() > 2 {
        return None;
    }
    let n = molecule.sequence.len();
    let flags = if let Some(index) = feature {
        location_flags(&molecule.features.get(index)?.location)?
    } else {
        LocationFlags::default()
    };
    if flags.fuzzy || flags.unresolved || flags.ordered {
        return None;
    }
    let (start, end) = match parts {
        [] if feature.is_none() => range?,
        [one] if !flags.compound && range == Some(*one) => *one,
        [(start, first_end), (second_start, end)]
            if molecule.topology == Topology::Circular
                && *start > 0
                && *start < n
                && *first_end == n
                && *second_start == 0
                && *end > 0
                && *end <= *start
                && range
                    .is_none_or(|current| current == parts[0] || current == (*start, n + *end)) =>
        {
            (*start, n + *end)
        }
        _ => return None,
    };
    (start < n
        && start < end
        && end - start <= n
        && end - start <= collagene_core::primers::MAX_OLIGO
        && (molecule.topology == Topology::Circular || end <= n))
        .then_some((start, end))
}
fn summarize(
    molecule: &Molecule,
    range: Option<(usize, usize)>,
    feature: Option<usize>,
    parts: &[(usize, usize)],
) -> Option<Summary> {
    if feature.is_none() && range.is_none() && parts.is_empty() {
        return None;
    }
    let len = molecule.sequence.len();
    if feature.is_none() {
        if let Some((a, _)) = range.filter(|(a, b)| {
            a == b && *b <= len && (parts.is_empty() || (parts.len() == 1 && parts[0] == (*a, *b)))
        }) {
            return Some(Summary {
                heading: format!("Insertion caret · boundary {a} of {len}"),
                detail: "Type reference symbols to insert at this boundary.".into(),
                annotate: false,
                primer: false,
                edit: true,
            });
        }
    }
    let (heading, flags) = if let Some(index) = feature {
        let Some(feature) = molecule.features.get(index) else {
            return Some(Summary::limited(
                "Selection changed".into(),
                "Select the current reference again.",
            ));
        };
        let heading = format!(
            "Annotation · {}",
            feature.label().chars().take(96).collect::<String>()
        );
        let Some(flags) = location_flags(&feature.location) else {
            return Some(Summary::limited(
                heading,
                "This location exceeds the selection summary limit. Inspect the source annotation.",
            ));
        };
        (heading, flags)
    } else {
        ("Reference selection".into(), LocationFlags::default())
    };
    // Check borrowed counts first. Never extract joined DNA or clone a molecule.
    if parts.len() > MAX_PARTS {
        return Some(Summary::limited(
            heading,
            "Many reference parts selected. Inspect the annotation or select a smaller reference range.",
        ));
    }
    if parts.is_empty() {
        if feature.is_none() {
            if let Some((a, b)) = range.filter(|(a, b)| a < b && *b <= len) {
                return summarize(molecule, range, None, &[(a, b)]);
            }
        }
        return Some(Summary::limited(
            heading,
            "No known local reference coordinates; unresolved locations are retained in the source.",
        ));
    }
    if parts.iter().any(|&(a, b)| a >= b || b > len) {
        return Some(Summary::limited(
            heading,
            "The selection no longer matches this source. Select the reference again.",
        ));
    }
    // At most 4096 interval pairs. Merge only for unique reference-length/GC
    // statistics; original parts and their order stay with the parent callback.
    let mut intervals = parts.to_vec();
    intervals.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::with_capacity(intervals.len());
    for (a, b) in intervals {
        if let Some(last) = merged.last_mut().filter(|last| a <= last.1) {
            last.1 = last.1.max(b);
        } else {
            merged.push((a, b));
        }
    }
    let Some(size) = merged
        .iter()
        .try_fold(0usize, |sum, &(a, b)| sum.checked_add(b - a))
    else {
        return Some(Summary::limited(
            heading,
            "Selected reference length exceeds the summary limit.",
        ));
    };
    let coordinates = parts
        .iter()
        .take(3)
        .map(|&(a, b)| format!("{}–{}", a + 1, b))
        .collect::<Vec<_>>()
        .join(", ");
    let mut detail = if parts.len() == 1 {
        format!(
            "{coordinates} · {size} {}",
            molecule.sequence.alphabet().short_unit()
        )
    } else {
        format!(
            "{} reference parts · {size} {} · {coordinates}{}",
            parts.len(),
            molecule.sequence.alphabet().short_unit(),
            if parts.len() > 3 { ", …" } else { "" }
        )
    };
    if flags.unresolved {
        detail.push_str(" · known portions only; unresolved coordinates remain");
    }
    if flags.fuzzy {
        detail.push_str(" · annotation endpoints are fuzzy");
    }
    if molecule.sequence.alphabet() != Alphabet::Protein
        && size <= MAX_GC_BASES
        && !flags.unresolved
    {
        let bytes = molecule.sequence.as_str().as_bytes();
        let (mut known, mut gc) = (0usize, 0usize);
        for &(a, b) in &merged {
            for &base in &bytes[a..b] {
                known += usize::from(b"ACGTU".contains(&base));
                gc += usize::from(b"GC".contains(&base));
            }
        }
        if known > 0 {
            detail.push_str(&format!(
                " · {:.1}% GC among {known} called A/C/G/T/U bases",
                gc as f64 / known as f64 * 100.
            ));
        } else {
            detail.push_str(" · GC unavailable for ambiguous bases");
        }
    }
    let contiguous = parts.len() == 1
        && !flags.compound
        && !flags.unresolved
        && !flags.fuzzy
        && range.is_some_and(|span| span == parts[0]);
    if contiguous
        && molecule.sequence.alphabet() == Alphabet::Dna
        && size > collagene_core::primers::MAX_OLIGO
    {
        detail.push_str(&format!(
            " · primer limit {} nt",
            collagene_core::primers::MAX_OLIGO
        ));
    }
    Some(Summary {
        heading,
        detail,
        // The parent preserves the exact selected feature Location, or every
        // known part as a Join. A first-part-only annotation is never implied.
        annotate: !flags.unresolved,
        primer: primer_geometry(molecule, range, feature, parts).is_some(),
        edit: contiguous,
    })
}

/// Reuse exactly the toolbar's eligibility rules for editor selection capture.
pub fn primer_range(
    molecule: &Molecule,
    range: Option<(usize, usize)>,
    feature: Option<usize>,
    parts: &[(usize, usize)],
) -> Result<(usize, usize), String> {
    primer_geometry(molecule, range, feature, parts).ok_or_else(||
        "Select one exact contiguous DNA domain of at most 200 bases and one circular turn; gaps, repeated parts, fuzzy, ordered and unresolved locations cannot seed a primer.".into())
}

/// A mutation intent may use a caret or a range longer than an oligo domain.
pub fn mutation_range(
    molecule: &Molecule,
    range: Option<(usize, usize)>,
    feature: Option<usize>,
    parts: &[(usize, usize)],
) -> Result<(usize, usize), String> {
    if molecule.sequence.alphabet() == Alphabet::Dna
        && summarize(molecule, range, feature, parts).is_some_and(|s| s.edit)
    {
        range.ok_or_else(|| "Select a DNA range or insertion caret.".into())
    } else {
        Err("Choose one exact DNA range or insertion caret; compound, fuzzy and unresolved locations cannot define a mutation intent.".into())
    }
}

#[cfg(target_arch = "wasm32")]
#[component]
pub fn SelectionToolbar(
    molecule: Memo<Molecule>,
    selection: RwSignal<Option<(usize, usize)>>,
    selected_feature: RwSignal<Option<usize>>,
    highlighted_spans: Memo<Vec<(usize, usize)>>,
    locked: Memo<bool>,
    on_annotate: Callback<()>,
    on_primer: Callback<bool>,
    on_mutation: Callback<()>,
    on_edit: Callback<()>,
    on_clear: Callback<()>,
) -> impl IntoView {
    let summary = Memo::new(move |_| {
        molecule.with(|molecule| {
            highlighted_spans
                .with(|parts| summarize(molecule, selection.get(), selected_feature.get(), parts))
        })
    });
    view! {
        <Show when=move||summary.with(Option::is_some)>
            <div class="selection-toolbar" role="region" aria-label="Selection actions">
                <div class="selection-toolbar-summary" role="status">
                    <strong>{move||summary.with(|value|value.as_ref().map(|s|s.heading.clone()).unwrap_or_default())}</strong>
                    <span>{move||summary.with(|value|value.as_ref().map(|s|s.detail.clone()).unwrap_or_default())}</span>
                </div>
                <div class="selection-toolbar-actions">
                    <button class="text-tool" disabled=move||locked.get()||!summary.with(|s|s.as_ref().is_some_and(|s|s.annotate)) title="Annotate the selected known reference parts" on:click=move|_|if !locked.get_untracked()&&summary.with_untracked(|s|s.as_ref().is_some_and(|s|s.annotate)){on_annotate.run(());}>"Annotate selection"</button>
                    <button class="text-tool" disabled=move||locked.get()||!summary.with(|s|s.as_ref().is_some_and(|s|s.primer)) title="Select one exact contiguous DNA range to make a primer" on:click=move|_|if !locked.get_untracked()&&summary.with_untracked(|s|s.as_ref().is_some_and(|s|s.primer)){on_primer.run(false);}>"Make forward primer"</button>
                    <button class="text-tool" disabled=move||locked.get()||!summary.with(|s|s.as_ref().is_some_and(|s|s.primer)) title="Use the reverse complement of one exact DNA range, written 5′→3′" on:click=move|_|if !locked.get_untracked()&&summary.with_untracked(|s|s.as_ref().is_some_and(|s|s.primer)){on_primer.run(true);}>"Make reverse primer"</button>
                    <button class="text-tool" disabled=move||locked.get()||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Dna)||!summary.with(|s|s.as_ref().is_some_and(|s|s.edit)) title="Declare a replacement, deletion or insertion as a new intended reference" on:click=move|_|if !locked.get_untracked()&&molecule.with_untracked(|m|m.sequence.alphabet()==Alphabet::Dna)&&summary.with_untracked(|s|s.as_ref().is_some_and(|s|s.edit)){on_mutation.run(());}>"Mutate selection"</button>
                    <button class="text-tool" disabled=move||locked.get()||!summary.with(|s|s.as_ref().is_some_and(|s|s.edit)) title="Edit one contiguous reference range or insertion caret" on:click=move|_|if !locked.get_untracked()&&summary.with_untracked(|s|s.as_ref().is_some_and(|s|s.edit)){on_edit.run(());}>"Edit bases"</button>
                    <button class="text-tool" disabled=move||locked.get() on:click=move|_|if !locked.get_untracked(){on_clear.run(());}>"Clear selection"</button>
                </div>
            </div>
        </Show>
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn circular_origin_pieces_seed_one_domain_without_enabling_direct_base_edit() {
        let molecule = Molecule::dna("circle", "ACGTACGTAA", Topology::Circular).unwrap();
        let parts = [(8, 10), (0, 3)];
        assert_eq!(
            primer_range(&molecule, None, None, &parts).unwrap(),
            (8, 13)
        );
        let summary = summarize(&molecule, None, None, &parts).unwrap();
        assert!(summary.primer && !summary.edit);
        assert!(mutation_range(&molecule, None, None, &parts).is_err());
        assert!(mutation_range(&molecule, Some((8, 10)), None, &parts).is_err());
        assert_eq!(
            primer_range(&molecule, Some((8, 10)), None, &parts).unwrap(),
            (8, 13)
        );
        assert_eq!(
            primer_range(&molecule, None, None, &[(3, 10), (0, 3)]).unwrap(),
            (3, 13)
        );
    }
    #[test]
    fn linear_gaps_multiturn_unknown_and_fuzzy_locations_remain_ineligible() {
        let mut molecule = Molecule::dna("circle", "ACGTACGTAA", Topology::Circular).unwrap();
        for parts in [
            vec![(8, 9), (0, 3)],
            vec![(8, 10), (1, 3)],
            vec![(2, 10), (0, 3)],
            vec![(8, 10), (0, 3), (8, 10)],
            vec![(0, 3), (8, 10)],
        ] {
            assert!(primer_range(&molecule, None, None, &parts).is_err());
        }
        let parts = [(8, 10), (0, 3)];
        for location in [
            Location::Order(vec![Location::range(8, 10), Location::range(0, 3)]),
            Location::Join(vec![
                Location::Range {
                    start: 8,
                    end: 10,
                    before: true,
                    after: false,
                },
                Location::range(0, 3),
            ]),
            Location::Join(vec![
                Location::range(8, 10),
                Location::Unresolved("remote".into()),
                Location::range(0, 3),
            ]),
        ] {
            molecule.features = vec![collagene_core::Feature {
                kind: "misc_feature".into(),
                location,
                qualifiers: vec![],
                needs_review: false,
            }];
            assert!(primer_range(&molecule, Some((8, 10)), Some(0), &parts).is_err());
        }
        molecule.features.clear();
        molecule.topology = Topology::Linear;
        assert!(primer_range(&molecule, None, None, &parts).is_err());
        assert_eq!(
            primer_range(&molecule, Some((1, 4)), None, &[(1, 4)]).unwrap(),
            (1, 4)
        );
        assert_eq!(
            mutation_range(&molecule, Some((1, 4)), None, &[(1, 4)]).unwrap(),
            (1, 4)
        );
    }
}
