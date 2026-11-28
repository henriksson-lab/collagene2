//! Explicit bounded initiation-pattern and ORF declarations in the existing QC owner.
use super::{construct_qc_model::Draft, construct_qc_view::feature_choice};
use collagene_core::{construct_qc as qc, sequence::translation};
use collagene_model::Molecule;
use leptos::prelude::*;

#[component]
pub(crate) fn ScanDeclarations(
    draft: RwSignal<Draft>,
    rows: Memo<Vec<(usize, String, bool)>>,
    get_molecule: Callback<(), Molecule>,
    status: RwSignal<String>,
) -> impl IntoView {
    let feature = RwSignal::new(String::new());
    let minimum = RwSignal::new("100".to_owned());
    let maximum = RwSignal::new("64".to_owned());
    let code_starts = RwSignal::new(false);
    let partials = RwSignal::new("complete".to_owned());
    let both = RwSignal::new(true);
    let add_context = Callback::new(move |(): ()| {
        let result = (|| {
            let index = feature
                .get_untracked()
                .parse::<usize>()
                .map_err(|_| "Choose a current CDS annotation.".to_owned())?;
            let profile =
                qc::initiation_context::kozak_1986_pattern().map_err(|e| e.to_string())?;
            let mut next = draft.get_untracked();
            if next
                .options
                .initiation_contexts
                .iter()
                .any(|claim| claim.cds_feature == index)
            {
                return Err(
                    "Remove the existing initiation-context declaration for this CDS first."
                        .to_owned(),
                );
            }
            next.options.initiation_contexts.push(qc::initiation_context::ContextClaim {
                cds_feature: index,
                profile,
                rationale: "User selected a literal comparison to the documented mammalian ACCATGG context; this is not a universal expression requirement.".into(),
            });
            next.options_for(&get_molecule.run(()))?;
            draft.set(next);
            Ok::<_, String>(())
        })();
        if let Err(error) = result {
            status.set(error);
        }
    });
    let apply_scan = Callback::new(move |(): ()| {
        let result = (|| {
            let minimum_amino_acids = minimum
                .get_untracked()
                .parse()
                .map_err(|_| "Enter a whole positive ORF length.".to_owned())?;
            let max_candidates = maximum
                .get_untracked()
                .parse()
                .map_err(|_| "Enter a whole candidate limit from 1 to 128.".to_owned())?;
            let incomplete = match partials.get_untracked().as_str() {
                "complete" => translation::IncompletePolicy::CompleteOnly,
                "linear" => translation::IncompletePolicy::IncludeLinearPartials,
                "all" => translation::IncompletePolicy::IncludeAllPartials,
                _ => return Err("Choose an explicit incomplete-ORF policy.".to_owned()),
            };
            let mut next = draft.get_untracked();
            next.options.orf_scan = Some(qc::unannotated_orfs::ScanOptions {
                start_policy: if code_starts.get_untracked() {
                    translation::StartPolicy::CodeStarts
                } else {
                    translation::StartPolicy::AtgOnly
                },
                incomplete,
                minimum_amino_acids,
                max_candidates,
                both_strands: both.get_untracked(),
            });
            next.options_for(&get_molecule.run(()))?;
            draft.set(next);
            Ok::<_, String>(())
        })();
        if let Err(error) = result {
            status.set(error);
        }
    });
    view! {
        <details>
            <summary>"Initiation context"</summary>
            <p class="helper">"Compare a declared CDS start to Kozak's experimentally studied ACCATGG sequence. This literal comparison is scoped to a declared mammalian context; conformity is not an initiation-strength or expression prediction."</p>
            {feature_choice(rows, feature, "CDS for initiation-context comparison")}
            <button class="button secondary" on:click=move |_| add_context.run(())>"Declare ACCATGG comparison"</button>
            {move || draft.with(|d| d.options.initiation_contexts.iter().enumerate().map(|(index, claim)| {
                let label = format!("CDS annotation {} · {}", claim.cds_feature + 1, claim.profile.id);
                view! {<p>{label}<button class="text-tool" on:click=move |_| draft.update(|d| {
                    if index < d.options.initiation_contexts.len() { d.options.initiation_contexts.remove(index); }
                })>"Remove context"</button></p>}
            }).collect_view())}
        </details>
        <details>
            <summary>"Unannotated ORF scan"</summary>
            <p class="helper">"Enumerate definite-start ORFs under the selected genetic code and compare them to declared CDS coverage. An ORF is a sequence candidate, not proof of an expressed protein. Limits and ambiguous boundaries remain explicit."</p>
            <label>"Minimum amino acids"<input type="number" min="1" step="1" aria-label="QC minimum ORF amino acids"
                prop:value=move || minimum.get() on:input=move |e| minimum.set(event_target_value(&e))/></label>
            <label>"Maximum retained candidates"<input type="number" min="1" max="128" step="1" aria-label="QC ORF candidate limit"
                prop:value=move || maximum.get() on:input=move |e| maximum.set(event_target_value(&e))/></label>
            <label class="checkbox-label"><input type="checkbox" prop:checked=move || code_starts.get()
                on:change=move |e| code_starts.set(event_target_checked(&e))/>"Use all starts in the selected genetic code"</label>
            <label class="checkbox-label"><input type="checkbox" prop:checked=move || both.get()
                on:change=move |e| both.set(event_target_checked(&e))/>"Scan both strands"</label>
            <label>"Incomplete ORFs"<select prop:value=move || partials.get()
                on:change=move |e| partials.set(event_target_value(&e))>
                <option value="complete">"Complete only"</option>
                <option value="linear">"Include linear edge candidates"</option>
                <option value="all">"Include linear and circular partial candidates"</option>
            </select></label>
            <button class="button secondary" on:click=move |_| apply_scan.run(())>"Apply ORF scan settings"</button>
            <button class="text-tool" on:click=move |_| draft.update(|d| d.options.orf_scan = None)>"Disable ORF scan"</button>
            <p class="helper">{move || draft.with(|d| match &d.options.orf_scan {
                Some(options) => format!("Committed: minimum {} aa · at most {} candidates · {:?} · {:?} · both strands {}", options.minimum_amino_acids, options.max_candidates, options.start_policy, options.incomplete, options.both_strands),
                None => "ORF scan is disabled.".into(),
            })}</p>
        </details>
    }
}
