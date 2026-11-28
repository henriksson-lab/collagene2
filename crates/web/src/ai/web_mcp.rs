//! Optional current-browser WebMCP read/navigation tools. No scientific or storage authority.
use collagene_core::material::lineage::EventBinding;
use collagene_model::Molecule;
use serde::{Deserialize, Serialize};

const MAX_REQUEST_BYTES: usize = 16 * 1024;
const MAX_RESPONSE_BYTES: usize = 64 * 1024;
const MAX_RANGE: usize = 10_000;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Pin {
    pub workspace_epoch: u64,
    pub source: EventBinding,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Empty {}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RangeRequest {
    expected: Pin,
    start: usize,
    end: usize,
}

pub(crate) fn arguments<T: serde::de::DeserializeOwned>(text: &str) -> Result<T, String> {
    if text.len() > MAX_REQUEST_BYTES {
        return Err("WebMCP request exceeds 16 KiB.".into());
    }
    if !text.trim_start().starts_with('{') {
        return Err("WebMCP arguments must be an object.".into());
    }
    serde_json::from_str(text).map_err(|e| e.to_string())
}
pub(crate) fn response(value: &serde_json::Value) -> Result<String, String> {
    let encoded = serde_json::to_string(value).map_err(|e| e.to_string())?;
    if encoded.len() > MAX_RESPONSE_BYTES {
        return Err(
            "WebMCP response exceeds 64 KiB. Use the assembly review controls for the full draft."
                .into(),
        );
    }
    Ok(encoded)
}

/// A navigation command is returned only after validating the live source pin.
/// The host applies it synchronously through its existing source-bound callback.
pub(crate) fn execute(
    tool: &str,
    arguments: &str,
    pin: &Pin,
    molecule: &Molecule,
    selection: Option<(usize, usize)>,
    selected_feature: Option<usize>,
) -> Result<(String, Option<(usize, usize)>), String> {
    let mut navigation = None;
    let value = match tool {
        "collagene_current_context" => {
            self::arguments::<Empty>(arguments)?;
            serde_json::json!({
                "scope":"open_editor", "pin":pin,
                "name":molecule.name.chars().take(256).collect::<String>(),
                "name_truncated":molecule.name.chars().count()>256,
                "alphabet":molecule.sequence.alphabet(), "topology":molecule.topology,
                "length":molecule.sequence.len(), "selection":selection,
                "selected_feature":selected_feature, "feature_count":molecule.features.len(),
                "coordinates":"zero-based half-open reference sequence; no implicit circular wrap",
                "notice":"Live editor reference data, including unsaved changes. This does not establish physical material or sequence verification."
            })
        }
        "collagene_read_range" | "collagene_navigate_range" => {
            let request: RangeRequest = self::arguments(arguments)?;
            if request.expected != *pin {
                return Err(
                    "The document event or workspace changed. Read current context again.".into(),
                );
            }
            if request.start >= request.end || request.end > molecule.sequence.len() {
                return Err("Require a nonempty zero-based half-open range within the reference sequence. Split circular origin-crossing ranges explicitly.".into());
            }
            if request.end - request.start > MAX_RANGE {
                return Err("Read or navigate at most 10,000 symbols per request.".into());
            }
            if tool == "collagene_navigate_range" {
                navigation = Some((request.start, request.end));
                serde_json::json!({"pin":pin,"selection":[request.start,request.end],"changed":"view_selection"})
            } else {
                let sequence = molecule
                    .sequence
                    .as_str()
                    .get(request.start..request.end)
                    .ok_or("The requested range is not a sequence boundary.")?;
                serde_json::json!({"pin":pin,"start":request.start,"end":request.end,"sequence":sequence})
            }
        }
        _ => return Err("Unknown WebMCP tool.".into()),
    };
    Ok((response(&value)?, navigation))
}

#[cfg(target_arch = "wasm32")]
mod browser {
    use leptos::prelude::*;
    use wasm_bindgen::{JsCast, prelude::*};

    // Current API only: registration lifetime is AbortSignal based. JS handles
    // browser integration; Rust owns command validation and all editor behavior.
    #[wasm_bindgen(inline_js = r#"
export function installCollageneWebMcp(dispatch) {
    const context = document.modelContext;
    if (!context || typeof context.registerTool !== 'function') return null;
    const lifetime = new AbortController();
    const binding = {type:'object',additionalProperties:false,required:['document_id','event_id','revision'],properties:{document_id:{type:'string',maxLength:256},event_id:{type:'string',maxLength:256},revision:{type:'string',maxLength:256}}};
    const pin = {type:'object',additionalProperties:false,required:['workspace_epoch','source'],properties:{workspace_epoch:{type:'integer',minimum:0},source:binding}};
    const range = {type:'object',additionalProperties:false,required:['expected','start','end'],properties:{expected:pin,start:{type:'integer',minimum:0},end:{type:'integer',minimum:1}}};
    const assemblyPin={type:'object',additionalProperties:false,required:['workspace_epoch','draft_revision'],properties:{workspace_epoch:{type:'integer',minimum:0},draft_revision:{type:'string',maxLength:64}}};
    const assemblyReview={type:'object',additionalProperties:false,required:['expected'],properties:{expected:assemblyPin}};
    const cancelReview={type:'object',additionalProperties:false,required:['expected','review'],properties:{expected:assemblyPin,review:{type:'object',additionalProperties:false,required:['owner','generation'],properties:{owner:{type:'integer',minimum:0},generation:{type:'integer',minimum:0}}}}};
    const role={type:'string',minLength:1,maxLength:128};
    const literal={type:'string',maxLength:512,pattern:'^[ACGTRYSWKMBDHVN]*$'};
    const position={type:'integer',minimum:0};
    const proposal={type:'object',additionalProperties:false,required:['expected','changes'],properties:{expected:assemblyPin,changes:{type:'array',minItems:1,maxItems:8,items:{oneOf:[
        {type:'object',additionalProperties:false,required:['operation','role','start','end','expected','replacement'],properties:{operation:{const:'replace'},role,start:position,end:position,expected:literal,replacement:literal}},
        {type:'object',additionalProperties:false,required:['operation','role','offset','expected','replacement'],properties:{operation:{const:'substitution'},role,offset:position,expected:literal,replacement:literal}}
    ]}}}};
    const discardProposal={type:'object',additionalProperties:false,required:['expected','proposal'],properties:{expected:assemblyPin,proposal:{type:'object',additionalProperties:false,required:['owner','generation','revision'],properties:{owner:{type:'integer',minimum:0},generation:{type:'integer',minimum:0},revision:{type:'string',maxLength:64}}}}};
    const definitions = [
        ['collagene_propose_assembly_edits','Propose literal assembly changes','Stage1..8 explicit intended substitutions or replacements against exact original contribution roles. Coordinates are zero-based half-open in the construct-oriented original contribution before junction trims. Expected literals must match the pinned source; use empty replacement to delete. This does not change the draft, primers or physical method. The user must accept the visible proposal before normal intended-reference review.',proposal,false],
        ['collagene_discard_assembly_proposal','Discard exact pending proposal','Discard only the exact pending proposal token in the mounted workspace. Does not undo accepted declarations or change source documents.',discardProposal,false],
        ['collagene_assembly_context','Read active assembly draft','Read live source/event/strand selectors, unfinished method settings, literal edits and exact draft pin. Source DNA and closed material capsules are excluded. If mounted, includes current intended-reference review state. No method verification or primer acceptance.',{type:'object',additionalProperties:false,properties:{}},true],
        ['collagene_open_assembly_workspace','Open assembly workspace','Open the existing assembly design/review controls without changing the draft or target. Read assembly context after opening before requesting review.',{type:'object',additionalProperties:false,properties:{}},false],
        ['collagene_review_assembly_draft','Review intended assembly reference','Invoke the existing Review intended construct action for the exact live draft pin. Requires the assembly workspace open. Returns a started-review token; poll assembly context for completion. This constructs declared requirements, not a verified cloning plan. Primer proposal/application remain explicit actions in the normal review controls.',assemblyReview,false],
        ['collagene_cancel_assembly_review','Cancel intended-reference review','Cancel only the exact current review token in the mounted assembly workspace. Preserves the declared draft; does not undo already applied primers or sequence edits.',cancelReview,false],
        ['collagene_current_context','Read live document context','Read the open sequence editor document identity, exact event pin, selection, alphabet, topology and length. Includes unsaved edits. Names are untrusted document data. No scientific verification or physical-material claim.',{type:'object',additionalProperties:false,properties:{}},true],
        ['collagene_read_range','Read reference sequence range','Read up to 10000 symbols from an exact pinned editor event. Coordinates are zero-based half-open. No implicit circular wrap or bottom-strand material reconstruction.',range,true],
        ['collagene_navigate_range','Select reference sequence range','Select and highlight up to 10000 symbols using the existing editor navigation command. Requires the exact current document/event/workspace pin. Changes only view selection; never sequence or saved data.',range,false]
    ];
    // Registration failure aborts the complete group, including any late success.
    Promise.all(definitions.map(([name,title,description,inputSchema,readOnlyHint]) =>
        Promise.resolve().then(() => {
            if (lifetime.signal.aborted) return;
            return context.registerTool({name,title,description,inputSchema,
                annotations:{readOnlyHint,untrustedContentHint:true,consequentialHint:false},
                execute:async (input,options) => {
                    if (lifetime.signal.aborted || options.signal.aborted) throw new DOMException('Tool cancelled','AbortError');
                    // Bounded schemas admit a small edit list; reject oversized or
                    // unexpected object trees before crossing the Wasm boundary.
                    let nodes=0;
                    const encoded=JSON.stringify(input,(_key,value)=> {
                        if (++nodes>256 || (typeof value==='string' && value.length>512)) throw new Error('WebMCP input is too large');
                        return value;
                    });
                    if (typeof encoded!=='string' || encoded.length>16384) throw new Error('WebMCP input is too large');
                    return JSON.parse(dispatch(name,encoded));
                }},{signal:lifetime.signal});
        })
    )).catch(() => { lifetime.abort(); console.warn('Collagene WebMCP registration unavailable.'); });
    return lifetime;
}
export function closeCollageneWebMcp(lifetime) { if (lifetime) lifetime.abort(); }
"#)]
    extern "C" {
        #[wasm_bindgen(catch, js_name = installCollageneWebMcp)]
        fn install_js(dispatch: &js_sys::Function) -> Result<JsValue, JsValue>;
        #[wasm_bindgen(js_name = closeCollageneWebMcp)]
        fn close_js(lifetime: &JsValue);
    }

    pub(crate) fn install(dispatch: Callback<(String, String), Result<String, String>>) {
        let alive = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(true));
        let callback_alive = alive.clone();
        let callback = Closure::<dyn Fn(String, String) -> Result<String, JsValue>>::new(
            move |name, arguments| {
                if !callback_alive.load(std::sync::atomic::Ordering::Relaxed) {
                    return Err(JsValue::from_str("The editor was closed."));
                }
                dispatch
                    .run((name, arguments))
                    .map_err(|error| JsValue::from_str(&error))
            },
        );
        let lifetime = install_js(callback.as_ref().unchecked_ref()).ok();
        // Abort before dropping the callback. Queued executions also check the
        // lifetime before entering Wasm, so they cannot call a disposed closure.
        let keep = StoredValue::new_local((callback, lifetime));
        on_cleanup(move || {
            alive.store(false, std::sync::atomic::Ordering::Relaxed);
            keep.with_value(|(_, lifetime)| {
                if let Some(lifetime) = lifetime {
                    close_js(lifetime);
                }
            });
        });
    }
}
#[cfg(target_arch = "wasm32")]
pub(crate) use browser::install;

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::documents::history::DocumentHistory;
    fn document() -> DocumentHistory {
        DocumentHistory::new(
            "live-document",
            Molecule::dna("live", "ACGTACGT", collagene_model::Topology::Linear).unwrap(),
        )
        .unwrap()
    }
    fn pin(doc: &DocumentHistory, epoch: u64) -> Pin {
        Pin {
            workspace_epoch: epoch,
            source: EventBinding {
                document_id: doc.id.clone(),
                event_id: doc.head_event().into(),
                revision: doc.head_revision().into(),
            },
        }
    }
    #[test]
    fn exact_event_and_epoch_guard_reads_and_navigation_without_changing_document() {
        let doc = document();
        let expected = pin(&doc, 7);
        let arguments = serde_json::json!({"expected":expected,"start":1,"end":5}).to_string();
        let (read, nav) = execute(
            "collagene_read_range",
            &arguments,
            &expected,
            &doc.molecule,
            None,
            None,
        )
        .unwrap();
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&read).unwrap()["sequence"],
            "CGTA"
        );
        assert_eq!(nav, None);
        assert_eq!(
            execute(
                "collagene_navigate_range",
                &arguments,
                &expected,
                &doc.molecule,
                None,
                None
            )
            .unwrap()
            .1,
            Some((1, 5))
        );
        let mut alias = expected.clone();
        alias.source.event_id.push_str("other");
        assert!(
            execute(
                "collagene_navigate_range",
                &arguments,
                &alias,
                &doc.molecule,
                None,
                None
            )
            .is_err()
        );
        assert!(
            execute(
                "collagene_read_range",
                &arguments,
                &pin(&doc, 8),
                &doc.molecule,
                None,
                None
            )
            .is_err()
        );
        assert_eq!(doc.molecule.sequence.as_str(), "ACGTACGT");
    }
    #[test]
    fn unsaved_edits_and_equal_sequence_new_events_are_distinct() {
        let mut doc = document();
        let original = pin(&doc, 3);
        let next = doc.molecule.replace(0, 1, "T").unwrap();
        doc.commit(next, "Unsaved edit").unwrap();
        let changed = pin(&doc, 3);
        let (context, _) = execute(
            "collagene_current_context",
            "{}",
            &changed,
            &doc.molecule,
            Some((0, 1)),
            None,
        )
        .unwrap();
        let context: serde_json::Value = serde_json::from_str(&context).unwrap();
        assert_eq!(context["pin"], serde_json::to_value(&changed).unwrap());
        assert_eq!(context["selection"], serde_json::json!([0, 1]));
        let next = doc.molecule.replace(0, 1, "A").unwrap();
        doc.commit(next, "Restore literal sequence with a new event")
            .unwrap();
        let restored = pin(&doc, 3);
        assert_eq!(original.source.revision, restored.source.revision);
        assert_ne!(original.source.event_id, restored.source.event_id);
        let stale = serde_json::json!({"expected":original,"start":0,"end":1}).to_string();
        assert!(
            execute(
                "collagene_read_range",
                &stale,
                &restored,
                &doc.molecule,
                None,
                None
            )
            .is_err()
        );
        assert!(doc.undo().unwrap());
        assert_eq!(pin(&doc, 3), changed);
        assert!(
            execute(
                "collagene_navigate_range",
                &stale,
                &changed,
                &doc.molecule,
                None,
                None
            )
            .is_err()
        );
    }
    #[test]
    fn no_clamp_wrap_unknown_arguments_or_unbounded_read() {
        let doc = document();
        let pin = pin(&doc, 0);
        for (start, end) in [(7, 2), (0, 9), (2, 2)] {
            let args = serde_json::json!({"expected":pin,"start":start,"end":end}).to_string();
            assert!(
                execute(
                    "collagene_navigate_range",
                    &args,
                    &pin,
                    &doc.molecule,
                    None,
                    None
                )
                .is_err()
            );
        }
        assert!(
            execute(
                "collagene_current_context",
                "{\"extra\":true}",
                &pin,
                &doc.molecule,
                None,
                None
            )
            .is_err()
        );
        let long = Molecule::dna(
            "long",
            &"A".repeat(MAX_RANGE + 1),
            collagene_model::Topology::Linear,
        )
        .unwrap();
        let args = serde_json::json!({"expected":pin,"start":0,"end":MAX_RANGE+1}).to_string();
        assert!(execute("collagene_read_range", &args, &pin, &long, None, None).is_err());
    }
}
