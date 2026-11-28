//! Editable source preparation controls and composition of independent review owners.
use crate::cloning::preparation::{
    cloning_methods_model::*,
    review::{MethodsReview, TargetContext},
};
use collagene_core::{
    chemistry::{EndRepairMethod, FivePrime, Terminus, ThreePrime},
    cloning::{
        ligation::{
            tailing::TailingBase,
            tailing_recipe::{Action, Profile},
        },
        overlap::profiles::OverlapMethod,
    },
    documents::transfer::Candidate,
    sequence::analysis::Primer,
};
use collagene_model::Molecule;
use leptos::{html, prelude::*};
use wasm_bindgen::JsCast;
fn preparation_label(mode: Mode) -> &'static str {
    match mode {
        Mode::Tailing => "Tailing",
        Mode::Ligation => "Ligation",
        Mode::GoldenGate => "Type IIS / Golden Gate",
        Mode::Overlap => "Overlap",
        Mode::MultiSiteReference => "Reference",
    }
}
#[component]
pub fn CloningMethodsPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    source_info: Callback<usize, Option<crate::sources::identity::SourceInfo>>,
    target_info: Memo<TargetContext>,
    selected_region: Callback<(), Option<(usize, usize, usize)>>,
    config: RwSignal<MethodsState>,
    epoch: RwSignal<u64>,
    on_candidate: Callback<Candidate>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
    primer_sources: Memo<Vec<crate::sources::identity::SourceInfo>>,
    list_primers: Callback<
        collagene_core::material::lineage::EventBinding,
        Result<collagene_core::cloning::pcr::reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primers::primer_picker::ResolveRequest,
        Result<collagene_core::cloning::pcr::reuse::Resolved, String>,
    >,
    #[prop(optional_no_strip)] get_workspace: Option<
        Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    >,
    #[prop(optional_no_strip)] on_finalize: Option<Callback<Candidate, Result<(), String>>>,
    #[prop(optional_no_strip)] on_source: Option<Callback<(usize, usize), Result<(), String>>>,
    // View state is separate from preparation/verification identity.
    #[prop(optional_no_strip)] active_row: Option<Memo<usize>>,
    #[prop(optional_no_strip)] on_edit_row: Option<Callback<usize>>,
    #[prop(default = false)] fixed_sources: bool,
    #[prop(default = false)] disable_plain_primer_design: bool,
) -> impl IntoView {
    let row_nodes = StoredValue::new_local(
        (0..12)
            .map(|_| NodeRef::<html::Div>::new())
            .collect::<Vec<_>>(),
    );
    let focus_row = Callback::new(move |index: usize| {
        if index >= config.with_untracked(|c| c.rows.len()).min(12) {
            return;
        }
        if let Some(select) = on_edit_row {
            select.run(index);
        }
        // The selected editor mounts in a reactive update. Focus after that
        // update; an old frame must not focus a different source/workspace.
        let expected_epoch = epoch.get_untracked();
        let expected_config = config.with_untracked(MethodsState::preparation_context);
        let frame = wasm_bindgen::closure::Closure::once_into_js(move || {
            if epoch.try_get_untracked() != Some(expected_epoch)
                || config.try_with_untracked(MethodsState::preparation_context)
                    != Some(expected_config)
                || active_row.is_some_and(|row| row.try_get_untracked() != Some(index))
            {
                return;
            }
            if let Some(node) = row_nodes
                .try_with_value(|nodes| nodes.get(index).and_then(|n| n.get()))
                .flatten()
            {
                node.scroll_into_view();
                let _ = node.focus();
            }
        });
        if let Some(window) = web_sys::window() {
            let _ = window.request_animation_frame(frame.unchecked_ref());
        }
    });
    let status = RwSignal::new(String::new());
    view! {<div class="inspector-content cloning-methods-panel" class:active-fragment-preparation=active_row.is_some()><Show when=move||!fixed_sources><h2>"Cloning methods"</h2></Show><p class="helper">{if fixed_sources {"Verify the prepared fragments against the separately reviewed intended construct."} else {"The active document is the independently chosen target. Method evidence never selects a target or establishes experimental success."}}</p><Show when=move||!fixed_sources>
        <label>"Method workflow"<select aria-label="Cloning workflow" disabled=disable_plain_primer_design prop:value=move||match config.get().mode{Mode::Tailing=>"tailing",Mode::Overlap=>"overlap",Mode::Ligation=>"ligation",Mode::GoldenGate=>"golden_gate",Mode::MultiSiteReference=>"multisite"} on:change=move|e|config.update(|c|c.mode=match event_target_value(&e).as_str(){"overlap"=>Mode::Overlap,"ligation"=>Mode::Ligation,"golden_gate"=>Mode::GoldenGate,_=>Mode::Tailing})><option value="tailing">"TA / GC tailing and ligation"</option><option value="overlap">"Gibson / HiFi / In-Fusion overlap"</option><option value="ligation">"Ordered multi-fragment ligation"</option><option value="golden_gate">"Type IIS / Golden Gate — ordered or finite reaction"</option></select></label></Show>
        <Show when=move||config.get().mode==Mode::Tailing><p class="helper">"Keep reconstructed existing ends, or add one declared nucleotide to both actual 3′ ends of an established blunt duplex. A pre-tailed vector is kept unchanged unless you explicitly choose an addition. Polymerase activity, mixtures and proprietary kit chemistry are not inferred."</p><label>"End-family profile"<select aria-label="Tailing end-family profile" prop:value=move||match config.get().tailing_profile{Profile::Ta=>"ta",Profile::Gc=>"gc"} on:change=move|e|config.update(|c|c.tailing_profile=if event_target_value(&e)=="gc"{Profile::Gc}else{Profile::Ta})><option value="ta">"TA · complementary 3′ T/A ends"</option><option value="gc">"GC · complementary 3′ C/G ends"</option></select></label><div class="method-shortcuts"><button class="text-tool" on:click=move|_|config.update(|c|{c.tailing_profile=Profile::Ta;c.vector.action=Action::AddTail{base:TailingBase::T};c.insert.action=Action::AddTail{base:TailingBase::A};})>"Declare T/A additions"</button><button class="text-tool" on:click=move|_|config.update(|c|{c.tailing_profile=Profile::Gc;c.vector.action=Action::AddTail{base:TailingBase::C};c.insert.action=Action::AddTail{base:TailingBase::G};})>"Declare C/G additions"</button></div><Show when=move||config.with(|c|c.rows.len()!=2)><p role="alert">"TA/GC cloning needs exactly two prepared fragments. Add the missing fragment or explicitly remove extra rows."</p></Show><p class="helper">"Choose two ordered prepared fragments below. The family assessment checks actual prepared ends separately from phosphate/OH eligibility and the independently declared target."</p><label class="checkbox-label"><input type="checkbox" aria-label="Allow retained tailing ligation nicks" prop:checked=move||config.get().allow_nicked on:change=move|e|config.update(|c|c.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label><p class="helper">"Exact PCR has 5′ OH termini. A sealed ligation may need explicit phosphorylation. Factory-vector chemistry and host nick repair are not inferred."</p></Show>
        <Show when=move||config.get().mode==Mode::Overlap><HomologyControls config/></Show>
        <Show when=move||config.get().mode==Mode::Overlap><label>"Published homology guidance"<select aria-label="Overlap profile method" prop:value=move||format!("{:?}",config.get().method) on:change=move|e|config.update(|c|c.method=match event_target_value(&e).as_str(){"Gibson"=>OverlapMethod::Gibson,"InFusion"=>OverlapMethod::InFusion,_=>OverlapMethod::HiFi})><option value="Gibson">"Gibson"</option><option value="HiFi">"HiFi"</option><option value="InFusion">"In-Fusion"</option></select></label><p class="helper">"Fragments are listed in product order with explicit orientation. Homology geometry is checked exactly; published length guidance is advisory. No reaction completion or unique mixture product is established."</p></Show>
            <Show when=move||matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate)><p class="helper">"List fragments in intended ligation order. Preparation, end repair and terminal treatment are explicit source-bound operations. Whole linear DNA retains unknown terminal chemistry; its sequence alone does not imply cleaved or PCR ends. Unknown ends can withhold ligation even after phosphorylation."</p><label class="checkbox-label"><input type="checkbox" aria-label="Allow retained ordered ligation nicks" prop:checked=move||config.get().allow_nicked on:change=move|e|config.update(|c|c.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label></Show>
            <h3>{if fixed_sources{"Prepare active fragment"}else{"1. Choose source fragments and regions"}}</h3>
            <Show when=move||config.with(|c|c.rows.len()>12)><p role="alert">"This preparation exceeds the supported 12 rows. No partial row list is shown."</p></Show><For each=move||config.with(|c|if c.rows.len()<=12{(0..c.rows.len()).collect::<Vec<_>>()}else{Vec::new()}) key=|i|*i children=move|i|view!{<Show when=move||active_row.is_none_or(|row|row.get()==i)><div node_ref=row_nodes.with_value(|nodes|nodes[i]) tabindex="-1" class="method-preparation-focus"><OverlapFields config index=i query_options primer_sources list_primers resolve_primer epoch selected_region fixed_sources edited_junctions=disable_plain_primer_design get_query source_info get_workspace/><Show when=move||config.get().mode==Mode::Tailing&&i<2><TailFields vector=i==0 config/></Show></div></Show>}/>
            <Show when=move||!fixed_sources&&config.with(|c|c.mode!=Mode::Tailing||c.rows.len()<2)><button class="button secondary full" disabled=move||{fixed_sources||config.with(|c|c.rows.len()>=12||(c.mode==Mode::Tailing&&c.rows.len()>=2))} on:click=move|_|{let mut result=Ok(());config.update(|c|result=c.push_row(OverlapInput::default()));if let Err(error)=result{status.set(error);}}>{move||if config.get().mode==Mode::Tailing{"Add TA/GC fragment"}else if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Add ligation fragment"}else{"Add overlap fragment"}}</button><p class="helper">"This panel supports up to 12 ordered fragments; the typed core has its separate documented work bounds."</p></Show><Show when=move||config.get().mode==Mode::Overlap><label>"Overlap lengths in junction order (bp)"<input aria-label="Profile overlap lengths" prop:value=move||config.get().overlaps on:input=move|e|config.update(|c|c.overlaps=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Scan repeated overlap words" prop:checked=move||config.get().scan_repeats on:change=move|e|config.update(|c|c.scan_repeats=event_target_checked(&e))/>"Search supplied fragments for repeated homology words"</label><label>"Repeat comparison budget"<input type="number" min="1" aria-label="Overlap repeat budget" prop:value=move||config.get().comparisons on:input=move|e|config.update(|c|c.comparisons=event_target_value(&e))/></label><label>"Maximum reported homology occurrences"<input type="number" min="1" aria-label="Overlap occurrence limit" prop:value=move||config.get().occurrences on:input=move|e|config.update(|c|c.occurrences=event_target_value(&e))/></label></Show>
        <Show when=move||config.get().mode==Mode::GoldenGate><p class="helper">"Ordered Type IIS restriction/ligation reconstructs the specified fragment order with positional cleavage checks. It does not model cycling or enumerate a reaction mixture. The separate finite Golden Gate reaction below digests explicit tube substrates and evaluates candidate recutting; its stated search limits remain explicit."</p></Show>
        <MethodsReview molecule get_query source_info target_info config epoch on_candidate on_products on_finalize on_source focus_row active_row fixed_sources status disable_plain_primer_design/>
    </div>}
}

#[component]
fn TailFields(vector: bool, config: RwSignal<MethodsState>) -> impl IntoView {
    let label = if vector { "Fragment 1" } else { "Fragment 2" };
    view! {<fieldset class="tail-source-fields"><legend>{format!("{label} · end processing")}</legend>
        <label>"End action"<select aria-label=format!("{label} tailing action") prop:value=move||config.with(|state|match state.input(vector).action{Action::KeepExisting=>"keep",Action::AddTail{base}=>base.as_str()}) on:change=move|event|config.update(|state|state.input_mut(vector).action=match event_target_value(&event).as_str(){"A"=>Action::AddTail{base:TailingBase::A},"T"=>Action::AddTail{base:TailingBase::T},"G"=>Action::AddTail{base:TailingBase::G},"C"=>Action::AddTail{base:TailingBase::C},_=>Action::KeepExisting})><option value="keep">"Keep existing prepared ends"</option><option value="A">"Add A to both 3′ ends"</option><option value="T">"Add T to both 3′ ends"</option><option value="G">"Add G to both 3′ ends"</option><option value="C">"Add C to both 3′ ends"</option></select></label>
        <label>"Terminal treatment after the selected action"<select aria-label=format!("{label} tailing terminal treatment") prop:value=move||config.with(|state|match state.input(vector).treatment{Treatment::Keep=>"keep",Treatment::Phosphorylate=>"phosphorylate",Treatment::Dephosphorylate=>"dephosphorylate"}) on:change=move|event|config.update(|state|state.input_mut(vector).treatment=match event_target_value(&event).as_str(){"phosphorylate"=>Treatment::Phosphorylate,"dephosphorylate"=>Treatment::Dephosphorylate,_=>Treatment::Keep})><option value="keep">"Keep reconstructed 5′ chemistry"</option><option value="phosphorylate">"Add 5′ phosphorylation"</option><option value="dephosphorylate">"Add 5′ dephosphorylation"</option></select></label>
        <p class="helper">"Keep existing emits no tailing step. Adding a tail requires freshly established blunt duplex ends; adding a second tail to an already tailed vector is not automatic."</p>
    </fieldset>}
}

#[component]
fn SourceGenerationFields(config: RwSignal<MethodsState>, index: usize) -> impl IntoView {
    let row = move || config.with(|c| c.rows.get(index).cloned().unwrap_or_default());
    view! {<details class="assembly-fragment-options"><summary>"Original source generation and modification assumptions"</summary>
        <Show when=move||row().run_origin.is_some()><p class="helper">"Saved Run/Derived material inherits its freshly reconstructed positional state. Original-source declarations are inactive here and cannot replace that state."</p></Show>
        <Show when=move||row().run_origin.is_none()><label class="checkbox-label"><input type="checkbox" aria-label=format!("Declare original source generation {}",index+1) prop:checked=move||row().source_generation.is_some() on:change=move|e|{let enabled=event_target_checked(&e);config.update(|c|if let Some(r)=c.rows.get_mut(index){if r.run_origin.is_none(){r.source_generation=if enabled{Some(Default::default())}else{None};}})}/>"Declare original source generation under the Dam/Dcm/CpG model"</label>
        <p class="helper">"Without a declaration, eligible source positions remain unknown. Declaring this closed model assigns states at original Dam/Dcm/CpG motifs and treats other eligible positions as unmodified; other modification systems are outside its scope. These are explicit source assumptions, not inferred host identity, vendor chemistry or imported metadata."</p>
        {[("Dam",0usize),("Dcm",1usize),("CpG",2usize)].into_iter().map(|(label,which)|view!{<label>{format!("{label} original-source state")}<select aria-label=format!("Fragment {} {label} source generation",index+1) disabled=move||row().run_origin.is_some()||row().source_generation.is_none() prop:value=move||{let form=row().source_generation.unwrap_or_default();match which {0=>form.dam,1=>form.dcm,_=>form.cpg}} on:change=move|e|{let value=event_target_value(&e);config.update(|c|if let Some(r)=c.rows.get_mut(index){if r.run_origin.is_none(){if let Some(form)=&mut r.source_generation{match which {0=>form.dam=value,1=>form.dcm=value,_=>form.cpg=value};}}})}><option value="unknown">"Unknown"</option><option value="unmodified">"Explicitly unmodified"</option><option value="fully">"Explicitly fully modified at eligible original motifs"</option><option value="mixed">"Hemi/mixed — unresolved under this model"</option></select></label>}).collect_view()}
        </Show>
        <Show when=move||row().has_pcr_amplicon()><p class="helper">"The declaration precedes this PCR. Standard newly synthesized PCR strands are unmodified under the model; original positional modifications are not copied onto them."</p></Show>
        <Show when=move||config.get().mode==Mode::GoldenGate><p class="helper">"Golden Gate preparations and recut checks use positional cleavage assessment. Unknown, impaired, unsupported or limited evidence can withhold a determined pool; no candidate or isolated clone is inferred from that uncertainty."</p></Show>
    </details>}
}

#[component]
fn OverlapFields(
    config: RwSignal<MethodsState>,
    index: usize,
    query_options: Memo<Vec<(usize, String)>>,
    primer_sources: Memo<Vec<crate::sources::identity::SourceInfo>>,
    list_primers: Callback<
        collagene_core::material::lineage::EventBinding,
        Result<collagene_core::cloning::pcr::reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primers::primer_picker::ResolveRequest,
        Result<collagene_core::cloning::pcr::reuse::Resolved, String>,
    >,
    epoch: RwSignal<u64>,
    selected_region: Callback<(), Option<(usize, usize, usize)>>,
    fixed_sources: bool,
    edited_junctions: bool,
    get_query: Callback<usize, Option<Molecule>>,
    source_info: Callback<usize, Option<crate::sources::identity::SourceInfo>>,
    get_workspace: Option<
        Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    >,
) -> impl IntoView {
    let row = move || config.with(|c| c.rows.get(index).cloned().unwrap_or_default());
    let region_status = RwSignal::new(String::new());

    view! {<fieldset class="overlap-source-row"><legend>{format!("Fragment {}",index+1)}</legend>
        <Show when=move||fixed_sources><p class="preparation-source-summary">{move||{
            let r=row();
            let name=query_options.with(|options|options.iter().find(|(i,_)|*i==r.document).map(|(_,name)|name.clone()).unwrap_or_else(||"Source unavailable".into()));
            if r.interval_start.is_empty() && r.interval_end.is_empty(){format!("{} · {} prepared fragment",name,if r.reverse{"reversed"}else{"forward"})}else{format!("{} · source {}–{} (1-based, unrolled) · {}",name,r.interval_start,if r.interval_end.is_empty(){"end"}else{&r.interval_end},if r.reverse{"reverse"}else{"forward"})}
        }}</p></Show>
        <Show when=move||!fixed_sources>
        <label>"Original source"<select disabled=fixed_sources aria-label=move||format!("{} source {}",preparation_label(config.get().mode),index+1) prop:value=move||row().document.to_string() on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){let next=event_target_value(&e).parse().unwrap_or(0);if r.document!=next{r.source_generation=None;}r.document=next})>{move||query_options.get().into_iter().map(|(n,label)|view!{<option value=n.to_string()>{label}</option>}).collect_view()}</select></label>
        <Show when=move||matches!(config.get().mode,Mode::Overlap|Mode::GoldenGate)><button class="button secondary" disabled=fixed_sources aria-label=format!("Set fragment {} from sequence selection",index+1) on:click=move|_|{match selected_region.run(()){Some((document,start,end)) if start<end&&start.checked_add(1).is_some()=>{config.update(|c|if let Some(row)=c.rows.get_mut(index){if row.document!=document{row.source_generation=None;}row.document=document;row.interval_start=(start+1).to_string();row.interval_end=end.to_string();});region_status.set("Selected source region copied. Review preparation and orientation; primers have not been changed.".into());},_=>region_status.set("Select one contiguous region in the sequence view first. Compound or origin-spanning feature selections need an explicit interval.".into())}} >"Set from sequence selection"</button><label>"Primer design interval start (1-based inclusive)"<input disabled=fixed_sources type="number" min="1" aria-label=format!("Overlap interval start {}",index+1) prop:value=move||row().interval_start on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.interval_start=event_target_value(&e)})/></label><label>"Primer design interval end (1-based inclusive; blank = source end)"<input disabled=fixed_sources type="number" min="1" aria-label=format!("Overlap interval end {}",index+1) prop:value=move||row().interval_end on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.interval_end=event_target_value(&e)})/></label><p class="helper">"Intervals declare intended regions for primer proposals. For a circular source, an end beyond its length describes one unrolled origin crossing. Existing whole-source, digest and manual PCR preparations are not silently sliced."</p></Show>
        </Show>
        <p class="helper preparation-region-status" role="status">{move||region_status.get()}</p>
        <Show when=move||config.get().mode==Mode::GoldenGate&&!edited_junctions><label>"Right fusion end (assembly direction, 5′→3′)"<input maxlength="4" aria-label=format!("Golden Gate right fusion {}",index+1) prop:value=move||row().right_fusion on:input=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.right_fusion=event_target_value(&e)})/></label><p class="helper">{move||format!("Enter {} concrete bases explicitly. The left fusion comes from this oriented source region; enter the next fragment’s first fusion bases explicitly for a matching junction. The right fusion is added design DNA. Existing-site compatibility is checked from its actual ends.",if config.with(|c|c.golden_gate.enzyme)==collagene_core::cloning::golden_gate::GoldenGateEnzyme::SapI{3}else{4})}</p></Show>
        <SourceGenerationFields config index/>
        <label>"Prepare fragment"<select aria-label=move||format!("{} preparation {}",preparation_label(config.get().mode),index+1) prop:value=move||match row().preparation{OverlapPreparation::WholeLinear=>"whole",OverlapPreparation::DeclaredDuplex=>"declared",OverlapPreparation::Digest=>"digest",OverlapPreparation::CatalogDigest=>"catalog",OverlapPreparation::ExactPcr=>"pcr",OverlapPreparation::MismatchPcr=>"mismatch_pcr",OverlapPreparation::MismatchPcrCatalogDigest=>"mismatch_pcr_catalog",OverlapPreparation::PcrCatalogDigest=>"pcr_catalog",OverlapPreparation::AnnealedOligos=>"anneal",OverlapPreparation::OverlapExtensionPcr=>"extension_pcr",OverlapPreparation::OverlapExtensionPcrCatalogDigest=>"extension_pcr_catalog"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.preparation=match event_target_value(&e).as_str(){"declared"=>OverlapPreparation::DeclaredDuplex,"digest"=>OverlapPreparation::Digest,"catalog"=>OverlapPreparation::CatalogDigest,"pcr"=>OverlapPreparation::ExactPcr,"mismatch_pcr"=>OverlapPreparation::MismatchPcr,"mismatch_pcr_catalog"=>OverlapPreparation::MismatchPcrCatalogDigest,"pcr_catalog"=>OverlapPreparation::PcrCatalogDigest,"anneal"=>OverlapPreparation::AnnealedOligos,"extension_pcr"=>{r.overlap_extension_pcr.get_or_insert_with(Default::default);OverlapPreparation::OverlapExtensionPcr},"extension_pcr_catalog"=>{r.overlap_extension_pcr.get_or_insert_with(Default::default);OverlapPreparation::OverlapExtensionPcrCatalogDigest},_=>OverlapPreparation::WholeLinear}})><option value="anneal" disabled=move||row().run_origin.is_some()||config.get().mode!=Mode::Ligation||(fixed_sources&&row().anneal_bottom_document.is_none())>"Annealed oligo pair"</option><option value="whole">"Whole linear DNA"</option><option value="declared" disabled=move||row().run_origin.is_some()>"Declare intact linear duplex"</option><option value="digest">"Manual starter-enzyme digest"</option><option value="catalog">"Pinned catalog digest fragment"</option><option value="pcr">"Exact PCR product"</option><option value="mismatch_pcr">"Intentional-mismatch PCR product"</option><option value="extension_pcr">"Overlap-extension PCR"</option><option value="extension_pcr_catalog">"Overlap-extension PCR → pinned catalog digest"</option><option value="pcr_catalog">"Exact PCR → pinned catalog digest"</option><option value="mismatch_pcr_catalog">"Intentional-mismatch PCR → pinned catalog digest"</option></select></label>
        {get_workspace.map(|snapshot|view!{<Show when=move||!fixed_sources&&(row().run_origin.is_some()||!matches!(row().preparation,OverlapPreparation::DeclaredDuplex|OverlapPreparation::AnnealedOligos))><crate::history::run_origin_view::RunOriginPicker config index source_info epoch get_workspace=snapshot/></Show>})}
        <Show when=move||row().preparation==OverlapPreparation::WholeLinear><p class="helper">"Use an already linear source. Circular DNA must first be prepared by restriction digest or exact PCR; its origin is not an assumed cut."</p></Show>
        <Show when=move||row().preparation==OverlapPreparation::AnnealedOligos><crate::cloning::preparation::anneal_preparation_view::AnnealPreparationFields config index query_options get_query source_info epoch fixed_sources/></Show>
        <Show when=move||row().preparation==OverlapPreparation::DeclaredDuplex><p class="helper">"Operator assumption for an original, intact linear duplex: declare its two terminal groups explicitly. This does not measure chemistry, cut or repair DNA, change its reference sequence, or reuse producer metadata. Circular, single-stranded and imported physical/reference-only intermediate sources are rejected by the core."</p><DeclaredEndFields config index left=true/><DeclaredEndFields config index left=false/><p class="helper">"All groups initially remain unknown. Orientation is applied after this declaration; any ligation terminal treatment remains a separate explicit step."</p></Show>
        <Show when=move||matches!(row().preparation,OverlapPreparation::Digest|OverlapPreparation::CatalogDigest|OverlapPreparation::PcrCatalogDigest|OverlapPreparation::MismatchPcrCatalogDigest|OverlapPreparation::OverlapExtensionPcrCatalogDigest)><Show when=move||row().preparation==OverlapPreparation::Digest>
            <label>"Restriction enzymes (comma-separated)"<input aria-label=move||format!("{} digest enzymes {}",preparation_label(config.get().mode),index+1) prop:value=move||row().enzymes on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.enzymes=event_target_value(&e)})/></label>
            <label>"Digest fragment number (1-based)"<input type="number" min="1" aria-label=move||format!("{} digest fragment {}",preparation_label(config.get().mode),index+1) prop:value=move||row().fragment on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.fragment=event_target_value(&e)})/></label>
            </Show><Show when=move||matches!(row().preparation,OverlapPreparation::PcrCatalogDigest|OverlapPreparation::MismatchPcrCatalogDigest|OverlapPreparation::OverlapExtensionPcrCatalogDigest)><crate::cloning::preparation::pcr_digest_view::PcrDigestPicker config index get_query source_info epoch get_workspace/></Show><Show when=move||matches!(row().preparation,OverlapPreparation::CatalogDigest|OverlapPreparation::PcrCatalogDigest|OverlapPreparation::MismatchPcrCatalogDigest|OverlapPreparation::OverlapExtensionPcrCatalogDigest)>{move||row().catalog_selection.map(|selection|view!{<p>{format!("Selected fragment {} · {}",row().fragment,selection.enzymes.iter().map(|e|e.enzyme.name.as_str()).collect::<Vec<_>>().join(", "))}</p><details><summary>"Pinned catalog identity"</summary><p>{format!("Edition {} · content {}",selection.catalog_edition,selection.catalog_content_hash)}</p>{selection.enzymes.iter().map(|entry|view!{<p>{format!("{} · {} · {} · cuts {}/{}",entry.id,entry.enzyme.name,entry.enzyme.recognition,entry.enzyme.top_cut,entry.enzyme.bottom_cut)}</p>}).collect_view()}</details>})}<p class="helper">"The exact catalog selection and fragment index are retained without starter-enzyme name substitution. For an original-source fragment, choose it in Virtual digest. PCR-product fragments use the explicit PCR/digest settings here. Saved pins describe a request; they do not restore current digest status or prove cleavage."</p></Show>
            <label>"Optional explicit end repair"<select aria-label=move||format!("{} end repair {}",preparation_label(config.get().mode),index+1) prop:value=move||match row().repair{None=>"none",Some(EndRepairMethod::FillFivePrimeOverhangs)=>"fill",Some(EndRepairMethod::RemoveThreePrimeOverhangs)=>"remove",Some(EndRepairMethod::Blunt)=>"blunt"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.repair=match event_target_value(&e).as_str(){"fill"=>Some(EndRepairMethod::FillFivePrimeOverhangs),"remove"=>Some(EndRepairMethod::RemoveThreePrimeOverhangs),"blunt"=>Some(EndRepairMethod::Blunt),_=>None}})><option value="none">"None — retain cohesive ends"</option><option value="fill">"Fill 5′ overhangs"</option><option value="remove">"Remove 3′ overhangs"</option><option value="blunt">"Blunt ends"</option></select></label>
            <p class="helper">"Fragments are reconstructed in core digest order using the stated complete-cleavage assumption. This is not evidence of enzyme activity, methylation permission or reagent compatibility."</p>
        </Show>
        <Show when=move||matches!(row().preparation,OverlapPreparation::ExactPcr|OverlapPreparation::PcrCatalogDigest)><label>"PCR template strand"<select aria-label=format!("Fragment {} PCR template strand",index+1) prop:value=move||match row().template_strand{collagene_core::chemistry::Strand::Top=>"top",collagene_core::chemistry::Strand::Bottom=>"bottom"} on:change=move|event|{let strand=match event_target_value(&event).as_str(){"bottom"=>collagene_core::chemistry::Strand::Bottom,_=>collagene_core::chemistry::Strand::Top};config.update(|c|if let Some(row)=c.rows.get_mut(index){row.template_strand=strand;});}><option value="top">"Top · original 5′→3′ coordinates"</option><option value="bottom" disabled=move||row().run_origin.is_none()>"Bottom · saved strand 5′→3′ coordinates"</option></select></label><p class="helper">"The whole selected strand is checked for competing primer sites. Selecting a contribution does not crop the template; nicked or unsupported source geometry remains a diagnostic."</p></Show>
        <Show when=move||row().is_overlap_extension_pcr()><crate::cloning::preparation::overlap_extension_preparation_view::OverlapExtensionPreparationFields config index/></Show>
        <Show when=move||row().is_pcr()>
            <crate::primers::primer_picker::PrimerPicker label=format!("Fragment {} forward",index+1) sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_name=primer.name;r.forward=primer.annealing.as_str().into();r.forward_tail=primer.tail;}))/>
            <label>"Forward primer name"<input maxlength="256" aria-label=move||format!("{} forward primer name {}",preparation_label(config.get().mode),index+1) prop:value=move||row().forward_name on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_name=event_target_value(&e)})/></label>
            <crate::primers::primer_picker::PrimerPicker label=format!("Fragment {} reverse",index+1) sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_name=primer.name;r.reverse_primer=primer.annealing.as_str().into();r.reverse_tail=primer.tail;}))/>
            <label>"Reverse primer name"<input maxlength="256" aria-label=move||format!("{} reverse primer name {}",preparation_label(config.get().mode),index+1) prop:value=move||row().reverse_name on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_name=event_target_value(&e)})/></label>

            <label>"Forward annealing sequence 5′→3′"<input aria-label=move||format!("{} forward primer {}",preparation_label(config.get().mode),index+1) prop:value=move||row().forward on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward=event_target_value(&e)})/></label>
            <label>"Forward primer 5′ extension"<input aria-label=move||format!("{} forward extension {}",preparation_label(config.get().mode),index+1) prop:value=move||row().forward_tail on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_tail=event_target_value(&e)})/></label>
            <label>"Reverse annealing sequence 5′→3′"<input aria-label=move||format!("{} reverse primer {}",preparation_label(config.get().mode),index+1) prop:value=move||row().reverse_primer on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_primer=event_target_value(&e)})/></label>
            <label>"Reverse primer 5′ extension"<input aria-label=move||format!("{} reverse extension {}",preparation_label(config.get().mode),index+1) prop:value=move||row().reverse_tail on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_tail=event_target_value(&e)})/></label>
            <Show when=move||matches!(row().preparation,OverlapPreparation::MismatchPcr|OverlapPreparation::MismatchPcrCatalogDigest)>
                <div class="primer-settings-grid">
                    <label>"Maximum mismatches per annealing domain"<input type="number" min="0" max="200" aria-label=format!("PCR maximum mismatches {}",index+1) prop:value=move||row().mismatch_max_mismatches on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.mismatch_max_mismatches=event_target_value(&e)})/></label>
                    <label>"Required exact 3′ bases"<input type="number" min="0" max="200" aria-label=format!("PCR exact three prime bases {}",index+1) prop:value=move||row().mismatch_three_prime_exact on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.mismatch_three_prime_exact=event_target_value(&e)})/></label>
                </div>
                <details><summary>"PCR search limits"</summary><div class="primer-grid"><label>"Maximum binding hits per primer"<input type="number" min="1" aria-label=format!("Mismatch PCR max_hits {}",index+1) prop:value=move||row().mismatch_max_hits on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.mismatch_max_hits=event_target_value(&e)})/></label><label>"Binding comparison limit per primer"<input type="number" min="1" aria-label=format!("Mismatch PCR max_work {}",index+1) prop:value=move||row().mismatch_max_work on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.mismatch_max_work=event_target_value(&e)})/></label><label>"Primer pair check limit"<input type="number" min="1" aria-label=format!("Mismatch PCR max_pairs {}",index+1) prop:value=move||row().mismatch_max_pairs on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.mismatch_max_pairs=event_target_value(&e)})/></label><label>"PCR product limit"<input type="number" min="1" aria-label=format!("Mismatch PCR max_products {}",index+1) prop:value=move||row().mismatch_max_products on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.mismatch_max_products=event_target_value(&e)})/></label><label>"Total PCR product base limit"<input type="number" min="1" aria-label=format!("Mismatch PCR max_product_bases {}",index+1) prop:value=move||row().mismatch_max_product_bases on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.mismatch_max_product_bases=event_target_value(&e)})/></label></div><p class="helper">"These limits are retained in inspection and preparation. A capped search can leave alternatives unresolved; any conditional geometry retains that evidence during verification."</p></details>
                <p class="helper">"Primer letters specify the introduced changes. Verification compares the reconstructed product with your independently declared target; primer settings never change that target. Internal changes must already be present in the intended reference. Binding compatibility does not establish amplification."</p>
            </Show>
            <details><summary>"Review current primer domains"</summary><pre>{move||{let r=row();format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}",r.forward_name,r.forward_tail,r.forward,r.forward,r.forward_tail,r.reverse_name,r.reverse_tail,r.reverse_primer,r.reverse_primer,r.reverse_tail)}}</pre></details>
            <p class="helper">"Enter annealing domains and added tails separately. Both oligos are written 5′→3′; no primer design or unique amplicon is assumed."</p>
        </Show>
        <Show when=move||!fixed_sources><label class="checkbox-label"><input disabled=fixed_sources type="checkbox" aria-label=format!("Reverse overlap source {}",index+1) prop:checked=move||row().reverse on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse=event_target_checked(&e)})/>"Reverse complement prepared fragment"</label></Show>
        <Show when=move||matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate)><label>"Explicit 5′ terminal treatment"<select aria-label=format!("Ligation terminal treatment {}",index+1) prop:value=move||match row().treatment{Treatment::Keep=>"none",Treatment::Phosphorylate=>"phosphorylate",Treatment::Dephosphorylate=>"dephosphorylate"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.treatment=match event_target_value(&e).as_str(){"phosphorylate"=>Treatment::Phosphorylate,"dephosphorylate"=>Treatment::Dephosphorylate,_=>Treatment::Keep}})><option value="none">"None — keep reconstructed chemistry"</option><option value="phosphorylate">"Add 5′ phosphorylation"</option><option value="dephosphorylate">"Add 5′ dephosphorylation"</option></select></label><p class="helper">"Treatment does not resolve unknown 3′ chemistry, repair ends or imply host nick repair."</p></Show>
        <p class="helper">{move||row().summary()}</p>
        <Show when=move||!fixed_sources><div class="method-row-actions"><button class="text-tool" disabled=move||fixed_sources||index==0 on:click=move|_|config.update(|c|{if index>0&&index<c.rows.len(){c.swap_rows(index,index-1);}})>"Move up"</button><button class="text-tool" disabled=move||fixed_sources||config.with(|c|index+1>=c.rows.len()) on:click=move|_|config.update(|c|{if index+1<c.rows.len(){c.swap_rows(index,index+1);}})>"Move down"</button><button class="text-tool" disabled=move||fixed_sources||config.with(|c|c.rows.len()<=1||(c.mode==Mode::Tailing&&c.rows.len()<=2)) on:click=move|_|config.update(|c|{if (c.mode!=Mode::Tailing||c.rows.len()>2)&&c.rows.len()>1&&index<c.rows.len(){c.remove_row(index);}})>"Remove fragment"</button></div></Show>
    </fieldset>}
}

#[component]
fn DeclaredEndFields(config: RwSignal<MethodsState>, index: usize, left: bool) -> impl IntoView {
    let end = move || {
        config.with(|c| {
            c.rows
                .get(index)
                .map(|row| {
                    if left {
                        row.left_terminus
                    } else {
                        row.right_terminus
                    }
                })
                .unwrap_or(Terminus::UNKNOWN)
        })
    };
    let side = if left { "Left" } else { "Right" };
    view! {<fieldset><legend>{format!("{side} source end")}</legend><p class="helper">{if left{"Top-strand 5′ terminus and bottom-strand 3′ terminus, before orientation."}else{"Bottom-strand 5′ terminus and top-strand 3′ terminus, before orientation."}}</p>
        <label>"Declared 5′ group"<select aria-label=format!("Fragment {} {side} declared 5 prime",index+1) prop:value=move||match end().five_prime{FivePrime::Phosphate=>"phosphate",FivePrime::Hydroxyl=>"hydroxyl",FivePrime::Unknown=>"unknown"} on:change=move|e|{let value=match event_target_value(&e).as_str(){"phosphate"=>FivePrime::Phosphate,"hydroxyl"=>FivePrime::Hydroxyl,_=>FivePrime::Unknown};config.update(|c|if let Some(row)=c.rows.get_mut(index){if left{row.left_terminus.five_prime=value}else{row.right_terminus.five_prime=value}})}><option value="unknown">"Unknown"</option><option value="phosphate">"Phosphate (P)"</option><option value="hydroxyl">"Hydroxyl (OH)"</option></select></label>
        <label>"Declared 3′ group"<select aria-label=format!("Fragment {} {side} declared 3 prime",index+1) prop:value=move||match end().three_prime{ThreePrime::Hydroxyl=>"hydroxyl",ThreePrime::Blocked=>"blocked",ThreePrime::Unknown=>"unknown"} on:change=move|e|{let value=match event_target_value(&e).as_str(){"hydroxyl"=>ThreePrime::Hydroxyl,"blocked"=>ThreePrime::Blocked,_=>ThreePrime::Unknown};config.update(|c|if let Some(row)=c.rows.get_mut(index){if left{row.left_terminus.three_prime=value}else{row.right_terminus.three_prime=value}})}><option value="unknown">"Unknown"</option><option value="hydroxyl">"Hydroxyl (OH)"</option><option value="blocked">"Blocked"</option></select></label>
    </fieldset>}
}

#[component]
fn HomologyControls(config: RwSignal<MethodsState>) -> impl IntoView {
    let form = move || config.with(|c| c.homology.clone().unwrap_or_default());
    view! {<details class="assembly-fragment-options"><summary>"Reconstruct product strands"</summary>
        <p class="helper">"Choose the processing assumptions separately from the published homology-length guidance. Reference checks do not process strands."</p>
        <label>"Strand processing"<select aria-label="Assembly strand processing route" prop:value=move||config.with(|c|match(c.homology.as_ref().is_some_and(|f|f.enabled),c.annealing.as_ref().is_some_and(|f|f.enabled)){(true,false)=>"five",(false,true)=>"three",(false,false)=>"reference",_=>"conflict"}) on:change=move|e|{let mode=event_target_value(&e);config.update(|c|{
            if let Some(form)=&mut c.homology{form.enabled=false;}
            if let Some(form)=&mut c.annealing{form.enabled=false;}
            match mode.as_str(){"five"=>c.homology.get_or_insert_with(Default::default).enabled=true,"three"=>c.annealing.get_or_insert_with(Default::default).enabled=true,_=>{}}
        });}><option value="reference">"Reference-only overlap check"</option><option value="five">"5′ resection → annealing, fill and eligible nick sealing"</option><option value="three">"3′ resection → 5′-tail annealing (unsealed)"</option><option value="conflict" disabled>"Conflicting saved modes — choose one"</option></select></label>
        <Show when=move||config.with(|c|c.method==OverlapMethod::InFusion)><p class="helper">"In-Fusion length guidance may be paired with the explicit 3′ resection route or the reference-only check. The generic 5′ fill/seal route does not establish In-Fusion chemistry; cellular rescue is not modeled."</p></Show>
        <Show when=move||config.with(|c|c.rows.len()==1&&!c.circular&&c.homology.as_ref().is_some_and(|f|f.enabled))><p class="helper">"A single linear fragment has no joining boundary. Leave junction homology empty; outer-end processing remains explicit below."</p></Show>

        <Show when=move||form().enabled&&config.with(|c|c.homology.is_some())>
            <label>"5′ resection at joining ends"<select prop:value=move||if form().exact_homology{"homology"}else{"explicit"} on:change=move|e|config.update(|c|c.homology.get_or_insert_with(Default::default).exact_homology=event_target_value(&e)=="homology")><option value="homology">"Exactly the declared junction homology"</option><option value="explicit">"Explicit cuts for every fragment"</option></select></label>
            <Show when=move||!form().exact_homology><label>"5′ cuts, left:right per oriented fragment"<input type="text" placeholder="20:20, 20:20" prop:value=move||form().cuts on:input=move|e|config.update(|c|c.homology.get_or_insert_with(Default::default).cuts=event_target_value(&e))/></label><p class="helper">"One comma-separated pair per prepared fragment, in product order. Left removes the top-strand 5′ end; right removes the bottom-strand 5′ end. Values refer to the oriented prepared fragment."</p></Show>
            <Show when=move||!config.get().circular><label>"Linear outer-end processing"<select prop:value=move||if form().all_accessible_ends{"all"}else{"joining"} on:change=move|e|config.update(|c|c.homology.get_or_insert_with(Default::default).all_accessible_ends=event_target_value(&e)=="all")><option value="joining">"Idealized scope: process joining ends only"</option><option value="all">"Also resect accessible outer ends"</option></select></label><p class="helper">"Joining-ends-only assumes outer termini are spared; it is an explicit idealization. Processing outer ends retains the resulting real 3′ overhangs rather than declaring a blunt product."</p><Show when=move||form().all_accessible_ends><div class="primer-grid"><label>"Left outer 5′ cut (bp)"<input type="text" prop:value=move||form().outer_left on:input=move|e|config.update(|c|c.homology.get_or_insert_with(Default::default).outer_left=event_target_value(&e))/></label><label>"Right outer 5′ cut (bp)"<input type="text" prop:value=move||form().outer_right on:input=move|e|config.update(|c|c.homology.get_or_insert_with(Default::default).outer_right=event_target_value(&e))/></label></div></Show></Show>
            <details><summary>"Terminal pairing and search limits"</summary><div class="primer-grid">
                <label>"Minimum paired terminal bases"<input type="text" prop:value=move||form().min_paired_bases on:input=move|e|config.update(|c|c.homology.get_or_insert_with(Default::default).min_paired_bases=event_target_value(&e))/></label>
                <label>"Terminal comparison limit"<input type="text" prop:value=move||form().max_comparisons on:input=move|e|config.update(|c|c.homology.get_or_insert_with(Default::default).max_comparisons=event_target_value(&e))/></label>
                <label>"Retained terminal candidates"<input type="text" prop:value=move||form().max_candidates on:input=move|e|config.update(|c|c.homology.get_or_insert_with(Default::default).max_candidates=event_target_value(&e))/></label>
            </div></details>
            <p class="helper">"Review every removed 5′ segment, joining nick and gap below. Unsupported gaps withhold a reusable product. Save a reconstructed product with its run history, then choose that exact saved material in a later preparation."</p>
        </Show>
        <Show when=move||config.with(|c|c.annealing.as_ref().is_some_and(|f|f.enabled))><AnnealingControls config/></Show>
    </details>}
}

#[component]
fn AnnealingControls(config: RwSignal<MethodsState>) -> impl IntoView {
    let form = move || config.with(|c| c.annealing.clone().unwrap_or_default());
    view! {<section class="three-prime-annealing-options">
        <p class="helper">"Shorten declared accessible 3′ ends, then pair the exposed 5′ homology tails. Surviving covalent strands remain separate. No gap fill, nick sealing or host rescue is performed."</p>
        <label><input type="checkbox" prop:checked=move||form().planned_ecoli_rescue on:change=move|e|config.update(|c|c.annealing.get_or_insert_with(Default::default).planned_ecoli_rescue=event_target_checked(&e))/>
            "Plan E. coli transformation and cellular rescue after annealing"
        </label>
        <p class="helper">"Required for the In-Fusion workflow. This records a planned next step; it does not simulate host repair, seal nicks or establish a reusable physical product."</p>
        <label>"3′ cuts, left:right per oriented prepared fragment"<input type="text" maxlength="1024" aria-label="Explicit 3 prime resection source cuts" placeholder="0:20, 20:0" prop:value=move||form().cuts on:input=move|e|config.update(|c|c.annealing.get_or_insert_with(Default::default).cuts=event_target_value(&e))/></label>
        <p class="helper">"One pair per fragment in product order. Left shortens the stored bottom-strand 3′ suffix; right shortens the stored top-strand 3′ suffix. Include any outer-end cuts explicitly. Original 5′ groups are retained; positive cuts require an accessible 3′ hydroxyl."</p>
        <details><summary>"Annealing and work limits"</summary><div class="primer-grid">
            <label>"Minimum paired terminal bases"<input type="text" maxlength="20" prop:value=move||form().min_paired_bases on:input=move|e|config.update(|c|c.annealing.get_or_insert_with(Default::default).min_paired_bases=event_target_value(&e))/></label>
            <label>"Terminal comparison limit"<input type="text" maxlength="20" prop:value=move||form().max_comparisons on:input=move|e|config.update(|c|c.annealing.get_or_insert_with(Default::default).max_comparisons=event_target_value(&e))/></label>
            <label>"Retained terminal candidates"<input type="text" maxlength="20" prop:value=move||form().max_candidates on:input=move|e|config.update(|c|c.annealing.get_or_insert_with(Default::default).max_candidates=event_target_value(&e))/></label>
        </div></details>
        <Show when=move||config.with(|c|c.rows.len()==1&&!c.circular)><p role="alert">"An annealing recipe needs two or more fragments, or one fragment for self-circularization. A standalone linear 3′ resection is a separate source-processing operation."</p></Show>
        <p class="helper">"This ordered model accepts intact linear sources with explicitly accounted blunt or staggered ends. Existing 5′ extensions contribute real homology; any 3′ extension outside the joined contour must be explicitly removed. It reports separate strands, gaps and nicks against a source-derived reference. The assembled reference stays incomplete and cannot be reused as physical material."</p>
    </section>}
}
