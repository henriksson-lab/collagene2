//! Desktop navigation metadata: the viewer/digest workflow and extensible secondary tools.
use leptos::{html, prelude::*};
#[derive(Clone, Copy)]
pub enum AuxiliaryAction {
    ImportUrl,
    Alignments,
    Oligos,
}
#[derive(Clone, Copy)]
enum Target {
    Panel(&'static str),
    Auxiliary(AuxiliaryAction),
}
impl Target {
    fn id(self) -> &'static str {
        match self {
            Self::Panel(id) => id,
            Self::Auxiliary(AuxiliaryAction::ImportUrl) => "import_url",
            Self::Auxiliary(AuxiliaryAction::Alignments) => "alignments",
            Self::Auxiliary(AuxiliaryAction::Oligos) => "oligos",
        }
    }
}

#[derive(Clone, Copy)]
struct Tool {
    label: &'static str,
    group: &'static str,
    target: Target,
}
const TOOLS: &[Tool] = &[
    Tool {
        label: "PCR",
        group: "Cloning",
        target: Target::Panel("pcr"),
    },
    Tool {
        label: "Primers",
        group: "Cloning",
        target: Target::Panel("primers"),
    },
    Tool {
        label: "Golden Gate",
        group: "Cloning",
        target: Target::Panel("golden_gate"),
    },
    Tool {
        label: "TOPO",
        group: "Cloning",
        target: Target::Panel("topo"),
    },
    Tool {
        label: "Recombination",
        group: "Cloning",
        target: Target::Panel("recombination"),
    },
    Tool {
        label: "Mutagenesis",
        group: "Cloning",
        target: Target::Panel("mutagenesis"),
    },
    Tool {
        label: "Remove internal sites",
        group: "Cloning",
        target: Target::Panel("domestication"),
    },
    Tool {
        label: "Plan",
        group: "Cloning",
        target: Target::Panel("plan"),
    },
    Tool {
        label: "End repair / annealing",
        group: "Cloning",
        target: Target::Panel("ends"),
    },
    Tool {
        label: "Cloning methods",
        group: "Cloning",
        target: Target::Panel("methods"),
    },
    Tool {
        label: "Methylation",
        group: "Sequence analysis",
        target: Target::Panel("methylation"),
    },
    Tool {
        label: "CRISPR",
        group: "Sequence analysis",
        target: Target::Panel("crispr"),
    },
    Tool {
        label: "Biology",
        group: "Sequence analysis",
        target: Target::Panel("biology"),
    },
    Tool {
        label: "Codon design",
        group: "Sequence analysis",
        target: Target::Panel("codons"),
    },
    Tool {
        label: "Compare",
        group: "Sequence analysis",
        target: Target::Panel("compare"),
    },
    Tool {
        label: "Dot plot",
        group: "Sequence analysis",
        target: Target::Panel("dotplot"),
    },
    Tool {
        label: "Gel",
        group: "Sequence analysis",
        target: Target::Panel("gel"),
    },
    Tool {
        label: "Alignments",
        group: "Sequence analysis",
        target: Target::Auxiliary(AuxiliaryAction::Alignments),
    },
    Tool {
        label: "Read validation",
        group: "Reads and references",
        target: Target::Panel("reads"),
    },
    Tool {
        label: "Read assembly",
        group: "Reads and references",
        target: Target::Panel("contigs"),
    },
    Tool {
        label: "Genome references",
        group: "Reads and references",
        target: Target::Panel("genomes"),
    },
    Tool {
        label: "Feature libraries",
        group: "Reads and references",
        target: Target::Panel("library"),
    },
    Tool {
        label: "Annotation files",
        group: "Reads and references",
        target: Target::Panel("annotation_files"),
    },
    Tool {
        label: "Import URL",
        group: "Files and history",
        target: Target::Auxiliary(AuxiliaryAction::ImportUrl),
    },
    Tool {
        label: "Oligo table",
        group: "Files and history",
        target: Target::Auxiliary(AuxiliaryAction::Oligos),
    },
    Tool {
        label: "History",
        group: "Files and history",
        target: Target::Panel("history"),
    },
];
const GROUPS: &[&str] = &[
    "Cloning",
    "Sequence analysis",
    "Reads and references",
    "Files and history",
];
#[component]
pub fn InspectorNavigation(
    panel: RwSignal<&'static str>,
    on_auxiliary: Callback<AuxiliaryAction>,
    on_assemble: Callback<()>,
) -> impl IntoView {
    let menu = NodeRef::<html::Details>::new();
    view! {
        <nav class="inspector-tabs mvp-inspector-tabs" aria-label="Analysis panels">
            <button data-panel-id="features" class:active=move||panel.get()=="features" on:click=move|_|{panel.set("features");if let Some(menu)=menu.get(){menu.set_open(false);}}>"Features"</button>
            <button data-panel-id="restriction" class:active=move||panel.get()=="restriction" on:click=move|_|{panel.set("restriction");if let Some(menu)=menu.get(){menu.set_open(false);}}>"RE sites"</button>
            <button data-panel-id="enzymes" class:active=move||panel.get()=="enzymes" on:click=move|_|{panel.set("enzymes");if let Some(menu)=menu.get(){menu.set_open(false);}}>"Digest"</button>
            <button data-panel-id="pcr" class:active=move||panel.get()=="pcr" on:click=move|_|{panel.set("pcr");if let Some(menu)=menu.get(){menu.set_open(false);}}>"PCR"</button>
            <button data-panel-id="methods" class:active=move||panel.get()=="methods" on:click=move|_|{on_assemble.run(());if let Some(menu)=menu.get(){menu.set_open(false);}}>"Assemble"</button>
            <button data-panel-id="plan" class:active=move||panel.get()=="plan" on:click=move|_|{panel.set("plan");if let Some(menu)=menu.get(){menu.set_open(false);}}>"Plan"</button>
            <details class="more-tools" node_ref=menu class:active=move||!matches!(panel.get(),"features"|"enzymes"|"restriction"|"pcr"|"methods"|"plan")>
                <summary aria-label="More tools">"More tools"<span class="active-tool-name">{move||TOOLS.iter().find_map(|t|match t.target{Target::Panel(id)if id==panel.get()=>Some(t.label),_=>None}).unwrap_or("")}</span></summary>
                <div class="more-tools-content">{GROUPS.iter().map(|group|view!{<section><h3>{*group}</h3>{TOOLS.iter().filter(|tool|tool.group==*group).map(|tool|{let target=tool.target;view!{<button class="secondary-tool" data-tool-id=target.id() on:click=move|_|{match target{Target::Panel(id)=>panel.set(id),Target::Auxiliary(action)=>on_auxiliary.run(action)}if let Some(menu)=menu.get(){menu.set_open(false);}}>{tool.label}</button>}}).collect_view()}</section>}).collect_view()}</div>
            </details>
        </nav>
    }
}
