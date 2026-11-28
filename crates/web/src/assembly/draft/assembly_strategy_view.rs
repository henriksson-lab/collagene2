//! Presentation-only strategy choice. Starting a draft does not verify a method.
use crate::cloning::preparation::cloning_methods_model::Mode;
use leptos::prelude::*;

#[component]
pub fn AssemblyStrategyPicker(initial: Mode, on_start: Callback<Mode>) -> impl IntoView {
    let selected = RwSignal::new(initial);
    view! {
        <section class="assembly-strategy-picker" aria-label="Choose assembly method">
            <fieldset>
                <legend>"How do you want to assemble the construct?"</legend>
                {[
                    (Mode::Ligation, "Digest and ligate", "Join prepared fragments with compatible physical ends."),
                    (Mode::Overlap, "Gibson / HiFi / In-Fusion", "Prepare fragments with overlaps and review their junctions."),
                    (Mode::GoldenGate, "Type IIS / Golden Gate", "Choose ordered Type IIS restriction/ligation or an explicit finite Golden Gate reaction."),
                    (Mode::Tailing, "TA / GC cloning", "Prepare two source fragments, add explicit 3′ tails and verify ligation against your declared junctions."),
                    (Mode::MultiSiteReference, "BP/LR reference exchange", "Reconstruct parent exchange under a declared reference profile."),
                ].into_iter().map(move |(mode, title, description)| view! {
                    <label class="assembly-strategy-choice" class:chosen=move||selected.get()==mode>
                        <input type="radio" name="assembly-strategy" prop:checked=move||selected.get()==mode on:change=move|_|selected.set(mode)/>
                        <span><strong>{title}</strong><small>{description}</small></span>
                    </label>
                }).collect_view()}
            </fieldset>
            <button class="button primary" on:click=move|_|on_start.run(selected.get_untracked())>"Start design"</button>
        </section>
    }
}
