//! Desktop navigation metadata: the viewer/digest workflow and extensible secondary tools.
use leptos::{html, prelude::*};
#[derive(Clone, Copy)]
pub enum AuxiliaryAction {
    ImportUrl,
    Alignments,
    Oligos,
    GoldenGateWorkflow,
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
            Self::Auxiliary(AuxiliaryAction::GoldenGateWorkflow) => "golden_gate_workflow",
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
        target: Target::Auxiliary(AuxiliaryAction::GoldenGateWorkflow),
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
pub fn panel_title(id: &str) -> &'static str {
    match id {
        "features" => "Annotations",
        "restriction" => "Restriction sites",
        "enzymes" => "Virtual digest",
        _ => TOOLS
            .iter()
            .find_map(|t| match t.target {
                Target::Panel(key) if key == id => Some(t.label),
                _ => None,
            })
            .unwrap_or("Tools"),
    }
}
fn rail_path(id: &str) -> &'static str {
    match id {
        "features" => "M4 5h16M4 12h11M4 19h16M7 3v4M18 10v4M8 17v4",
        "restriction" => {
            "M9 8l11 12M9 16 20 4M4 4a3 3 0 1 0 0 6 3 3 0 0 0 0-6M4 14a3 3 0 1 0 0 6 3 3 0 0 0 0-6"
        }
        "enzymes" => "M8 3h8M10 3v7l-6 9v2h16v-2l-6-9V3M7 16h10",
        "primers" => "M3 7h15l-4-4M18 7l-4 4M21 17H6l4-4M6 17l4 4",
        "history" => "M3 11a9 9 0 1 1 2 7M3 4v7h7M12 7v6l4 2",
        _ => "M12 4v.01M12 12v.01M12 20v.01",
    }
}
#[component]
pub fn InspectorNavigation(
    panel: RwSignal<&'static str>,
    on_auxiliary: Callback<AuxiliaryAction>,
    on_assemble: Callback<()>,
) -> impl IntoView {
    let menu = NodeRef::<html::Details>::new();
    view! {
        <nav class="viewer-tool-rail" aria-label="Sequence tools">
            {[("features","Annotations"),("restriction","Restriction sites"),("enzymes","Virtual digest"),("primers","Primers"),("history","History")].into_iter().map(|(id,label)|view!{
                <button class="viewer-rail-button" class:active=move||panel.get()==id title=label aria-label=label aria-expanded=move||(panel.get()==id).to_string() data-panel-id=id on:click=move |_|{panel.update(|current|*current=if *current==id{""}else{id});if let Some(menu)=menu.get(){menu.set_open(false);}}><svg viewBox="0 0 24 24" aria-hidden="true"><path d=rail_path(id)/></svg></button>
            }).collect_view()}
            <details class="viewer-tools-menu" data-transient-menu="" node_ref=menu><summary class="viewer-rail-button" title="More sequence tools" aria-label="More sequence tools"><svg viewBox="0 0 24 24" aria-hidden="true"><circle cx="12" cy="5" r="1.5"/><circle cx="12" cy="12" r="1.5"/><circle cx="12" cy="19" r="1.5"/></svg></summary><div class="viewer-tool-menu-content">
                <button class="secondary-tool" on:click=move |_|{on_assemble.run(());if let Some(menu)=menu.get(){menu.set_open(false);}}>"Assembly workspace"</button>
                {GROUPS.iter().map(|group|view!{<section><h3>{*group}</h3>{TOOLS.iter().filter(|tool|tool.group==*group).map(|tool|{let target=tool.target;view!{<button class="secondary-tool" data-tool-id=target.id() on:click=move |_|{match target{Target::Panel(id)=>panel.set(id),Target::Auxiliary(action)=>on_auxiliary.run(action)}if let Some(menu)=menu.get(){menu.set_open(false);}}>{tool.label}</button>}}).collect_view()}</section>}).collect_view()}
            </div></details>
        </nav>
    }
}
