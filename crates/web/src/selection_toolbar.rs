//! Contextual actions over the existing reference selection. No edits or source
//! copies happen here; the app owns annotation/primer dialogs and history.
use collagene_core::{Alphabet, Location, Molecule};
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
        primer: contiguous
            && molecule.sequence.alphabet() == Alphabet::Dna
            && size <= collagene_core::primers::MAX_OLIGO,
        edit: contiguous,
    })
}

#[component]
pub fn SelectionToolbar(
    molecule: Memo<Molecule>,
    selection: RwSignal<Option<(usize, usize)>>,
    selected_feature: RwSignal<Option<usize>>,
    highlighted_spans: Memo<Vec<(usize, usize)>>,
    locked: Memo<bool>,
    on_annotate: Callback<()>,
    on_primer: Callback<()>,
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
                    <button class="text-tool" disabled=move||locked.get()||!summary.with(|s|s.as_ref().is_some_and(|s|s.primer)) title="Select one exact contiguous DNA range to make a primer" on:click=move|_|if !locked.get_untracked()&&summary.with_untracked(|s|s.as_ref().is_some_and(|s|s.primer)){on_primer.run(());}>"Make primer"</button>
                    <button class="text-tool" disabled=move||locked.get()||!summary.with(|s|s.as_ref().is_some_and(|s|s.edit)) title="Edit one contiguous reference range or insertion caret" on:click=move|_|if !locked.get_untracked()&&summary.with_untracked(|s|s.as_ref().is_some_and(|s|s.edit)){on_edit.run(());}>"Edit bases"</button>
                    <button class="text-tool" disabled=move||locked.get() on:click=move|_|if !locked.get_untracked(){on_clear.run(());}>"Clear selection"</button>
                </div>
            </div>
        </Show>
    }
}
