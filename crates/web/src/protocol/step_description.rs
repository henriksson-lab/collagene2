use collagene_core::documents::transfer;
use collagene_core::planning::plan::Operation;
use collagene_core::planning::plan::Source;
use collagene_core::{chemistry::tracked_duplex as modified_duplex, verification};

use leptos::prelude::*;

pub(crate) fn describe_step(step: &Operation) -> String {
    let source = |s: &Source| match s {
        Source::Input(id) => format!("input {id}"),
        Source::Step(i) => i
            .checked_add(1)
            .map(|number| format!("product of step {number}"))
            .unwrap_or_else(|| format!("invalid step index {i}")),
    };
    match step {
        Operation::OverlapExtension {
            left,
            right,
            options,
        } => format!(
            "Reciprocal overlap extension of {} ({}) / {} ({})\nSelected terminal homology: {} bp\nRetain actual source strands and extend paired 3′ ends; alternative terminal pairings remain diagnostic.",
            source(&left.source),
            if left.reverse {
                "reverse complement"
            } else {
                "forward"
            },
            source(&right.source),
            if right.reverse {
                "reverse complement"
            } else {
                "forward"
            },
            options.overlap
        ),
        Operation::ReconstructMaterial { template, origin } => format!(
            "Reconstruct reusable PCR/digest material from {} · exact event {}",
            source(template),
            origin.product.event_id
        ),

        Operation::ReconstructRunProduct { template, origin } => format!(
            "Reconstruct saved material from {}\nOrigin run: {} · product {}\nUse the unchanged originating recipe and target; retain nicks, termini, positional modifications and all incomplete findings.",
            source(template),
            origin.selection.run_id,
            origin.selection.step_index + 1
        ),
        Operation::CatalogDigest {
            template,
            selection,
            fragment_index,
        } => format!(
            "Catalog-pinned complete digest of {}\nEnzymes: {}\nCatalog edition: {}\nZero-based fragment: {}\nCleavage assumes compatible conditions and no blocking methylation; no experimental success is asserted.",
            source(template),
            selection
                .enzymes
                .iter()
                .map(|e| e.enzyme.name.as_str())
                .collect::<Vec<_>>()
                .join(", "),
            selection.catalog_edition,
            fragment_index
        ),
        Operation::DeclareSourceTermini {
            template,
            left,
            right,
        } => format!(
            "Original intact linear duplex: {}\nDeclared left terminus: {:?}\nDeclared right terminus: {:?}\nOperator assumptions only; no measurement, synthesis or physical graph repair is inferred.",
            source(template),
            left,
            right
        ),
        Operation::TopoActions { sources, request } => format!(
            "Accounted TOPO chemical actions\n{}\nLiteral declarations: {}\nIdeal hydrolytic profiles are explicit model assumptions; excised material remains evidence, not a claimed isolated product.",
            sources
                .iter()
                .map(|s| format!("{} = {}", s.role, source(&s.source)))
                .collect::<Vec<_>>()
                .join("\n"),
            serde_json::to_string_pretty(request).unwrap_or_default()
        ),
        Operation::CheckedTopo {
            sources,
            request,
            selection,
            policy,
        } => format!(
            "Final source-backed TOPO graph with explicit independent-target acceptance\n{}\nPreparation/transfers: {}\nSelected component/seed: {}\nAcceptance policy: {}\nAccepted allowances do not repair nicks, eliminate other components or establish experimental success. The sequence view remains chemically restricted.",
            sources
                .iter()
                .map(|entry| format!("{} = {}", entry.role, source(&entry.source)))
                .collect::<Vec<_>>()
                .join("\n"),
            serde_json::to_string_pretty(request).unwrap_or_default(),
            serde_json::to_string(selection).unwrap_or_default(),
            serde_json::to_string_pretty(policy).unwrap_or_default()
        ),
        Operation::Topo {
            sources,
            request,
            selection,
        } => format!(
            "Source-backed TOPO strand transfers\n{}\nPreparation and transfer recipe: {}\nSelected component: {}\nPhysical graph and paired-contour evidence are retained in the reconstruction report. A reference preview does not establish repaired, intact duplex chemistry or experimental success.",
            sources
                .iter()
                .map(|s| format!("{} = {}", s.role, source(&s.source)))
                .collect::<Vec<_>>()
                .join("\n"),
            serde_json::to_string_pretty(request).unwrap_or_default(),
            serde_json::to_string(selection).unwrap_or_default()
        ),
        Operation::MultiSiteRecombineReference {
            sources,
            request,
            selection,
        } => format!(
            "Multi-parent {:?} reference exchange\n{}\nProfile {} / {}\nSpecificities {:?}\nSelected candidate {}, component {} from batch {}. All reciprocal components retained; no chemical strand authority.",
            request.reaction.reaction,
            sources
                .iter()
                .map(|s| format!("{} = {}", s.role, source(&s.source)))
                .collect::<Vec<_>>()
                .join("\n"),
            request.reaction.profile.id,
            request.reaction.profile.edition,
            request.reaction.specificities,
            selection.candidate_id + 1,
            selection.component_id + 1,
            selection.batch_revision
        ),
        Operation::RecombineReference {
            sources,
            request,
            selection,
        } => format!(
            "Reference-only {:?} exchange of {} and {}\nExplicit profile {} / {} (revision {})\nSelected candidate {}, component {} from batch {}.\nAll reciprocal components are computed under the declared full-site profile and bounded limits. Selection does not establish physical strand state, biological repair, efficacy or an experimental winner.",
            request.reaction,
            source(&sources[0]),
            source(&sources[1]),
            request.profile.profile.id,
            request.profile.profile.edition,
            request.profile.revision,
            selection.candidate_id,
            selection.component_id,
            selection.batch_revision
        ),
        Operation::DeclareSourceModifications {
            template,
            assumptions,
        } => format!(
            "Declare source generation for {}\nDam: {}\nDcm: {}\nCpG: {}\nClosed Dam/Dcm/CpG model: modification is assigned only to eligible positions in the original source motifs. Other original positions and newly synthesized bases are unmodified under this scope. Descendant motifs are not rescanned to infer inheritance. This is an assumption, not imported physical evidence.",
            source(template),
            crate::cloning::methods::checked_cloning::state_label(assumptions.dam),
            crate::cloning::methods::checked_cloning::state_label(assumptions.dcm),
            crate::cloning::methods::checked_cloning::state_label(assumptions.cpg)
        ),
        Operation::CheckedDigest {
            template,
            selection,
            options,
            fragment_index,
        } => format!(
            "Checked restriction digest of {}\nEnzymes: {}\nPinned positional supplier profile revision: {}\nSelect fragment {} in ascending reference top-cut order. Unknown or impaired sites withhold undetermined products; blocked sites remain uncut under this model. Reagent conditions and experimental success remain outside the check.",
            source(template),
            selection
                .enzymes
                .iter()
                .map(|e| e.enzyme.name.as_str())
                .collect::<Vec<_>>()
                .join(", "),
            format!(
                "{} · catalog {} · {} comparisons / {} sites",
                modified_duplex::profile_revision(),
                selection.catalog_content_hash,
                options.max_comparisons,
                options.max_sites
            ),
            fragment_index
                .checked_add(1)
                .map(|n| n.to_string())
                .unwrap_or_else(|| "out of range".into())
        ),
        Operation::Tail3Prime { template, base } => format!(
            "Substrate: {}\nAdd one non-templated {} to both actual 3′ ends. Source end chemistry and nicks are retained; enzyme activity and completion are assumptions.",
            source(template),
            base.as_str()
        ),
        Operation::ResectThreePrime {
            source: substrate,
            cuts,
        } => format!(
            "Explicit 3′ resection of {}\nLeft bottom-strand cut: {} bp · right top-strand cut: {} bp\nAccessible 3′ suffixes are removed; original 5′ groups remain. No fill or sealing is inferred.",
            source(substrate),
            cuts.left,
            cuts.right
        ),
        Operation::HomologyAnnealing {
            fragments, policy, ..
        } => format!(
            "3′ resection and homology annealing\n{} oriented fragments · {:?}\nJoining homology: {:?} bp\nExplicit cuts: {:?}\nThe joined reference remains incomplete. Unfilled gaps and unsealed nicks are retained; no reusable assembled material or host processing is inferred.",
            fragments.len(),
            policy.topology,
            policy.overlaps,
            policy.cuts
        ),
        Operation::HomologyAssembly {
            fragments, policy, ..
        } => format!(
            "Physical homology assembly\n{} oriented fragments · {:?}\nJoining homology: {:?} bp\n5′ resection: {:?}\nOuter-end policy: {:?}\nGeneric resection/annealing/gap-fill/ligation assumptions; inspect processing evidence before saving reusable material.",
            fragments.len(),
            policy.topology,
            policy.overlaps,
            policy.resection,
            policy.outer_end_policy
        ),
        Operation::ProfiledOverlap {
            fragments,
            overlaps,
            circular,
            options,
        } => format!(
            "{:?} overlap profile\nOrdered fragments: {}\nJunction homology: {:?}bp\nTopology: {}\nPublished length guidance is advisory; unique reaction outcomes and reagent activity are not established.",
            options.method,
            fragments
                .iter()
                .map(|f| format!(
                    "{} ({})",
                    source(&f.source),
                    if f.reverse {
                        "reverse complement"
                    } else {
                        "forward"
                    }
                ))
                .collect::<Vec<_>>()
                .join(" → "),
            overlaps,
            if *circular { "circular" } else { "linear" }
        ),
        Operation::RepairEnds { template, method } => format!(
            "Substrate: {}\nEnd repair: {:?}\nOnly modeled 3′ fill-in/exonuclease changes occur; 5′ chemistry and internal nicks remain explicit.",
            source(template),
            method
        ),
        Operation::AnnealOligos {
            top,
            bottom,
            options,
            candidate_index,
        } => format!(
            "Ordered 5′→3′ oligos: {} / {}\nExplicit zero-based annealing candidate {}\nDeclared chemistry/search policy: {:?}\nAnnealing creates no new phosphodiester bonds.",
            source(top),
            source(bottom),
            candidate_index,
            options
        ),
        Operation::PcrWithMismatches {
            template,
            forward,
            reverse,
            options,
        } => format!(
            "Template: {}\nOligos 5′→3′: {} / {}\nBinding/search policy: {:?}",
            source(template),
            forward.ordered_sequence(),
            reverse.ordered_sequence(),
            options
        ),
        Operation::SelectedGoldenGate {
            templates,
            options,
            candidate_index,
            directed_sequence_key,
            physical_sequence_key,
        } => format!(
            "Golden Gate inputs: {}\nExplicit zero-based candidate {}\nDirected sequence key {}\nPhysical sequence key {}\nPolicy: {:?}",
            templates.iter().map(source).collect::<Vec<_>>().join(", "),
            candidate_index,
            directed_sequence_key,
            physical_sequence_key,
            options
        ),
        Operation::GoldenGate {
            templates,
            options,
            candidate_index,
        } => format!(
            "Golden Gate inputs: {}\nExplicit zero-based candidate {}\nPolicy: {:?}",
            templates.iter().map(source).collect::<Vec<_>>().join(", "),
            candidate_index,
            options
        ),
        Operation::Pcr {
            template_strand,
            template,
            forward,
            reverse,
        } => format!(
            "Template: {} · {:?} strand frame\nForward primer {} (5′→3′): {}\nReverse primer {} (5′→3′): {}",
            source(template),
            template_strand,
            forward.name,
            forward.ordered_sequence(),
            reverse.name,
            reverse.ordered_sequence()
        ),
        Operation::Digest {
            template,
            enzymes,
            fragment_index,
        } => format!(
            "Template: {}\nComplete digest with: {}\nSelected zero-based fragment index: {} (ascending top-cut order)",
            source(template),
            enzymes.join(", "),
            fragment_index
        ),
        Operation::TreatEnds {
            template,
            treatment,
        } => format!(
            "Substrate: {}\nTerminal treatment: {:?}\nInternal nick repair is not assumed.",
            source(template),
            treatment
        ),
        Operation::Ligate {
            fragments,
            circular,
            allow_nicked,
            options,
        } => format!(
            "Ordered substrates: {}\nProduct topology: {}\nRetained nicks explicitly allowed: {}\nAlternative-enumeration bounds: {:?}",
            fragments
                .iter()
                .map(|f| format!(
                    "{} ({})",
                    source(&f.source),
                    if f.reverse { "reversed" } else { "forward" }
                ))
                .collect::<Vec<_>>()
                .join(" → "),
            if *circular { "circular" } else { "linear" },
            allow_nicked,
            options
        ),
        Operation::Overlap {
            fragments,
            overlaps,
            circular,
        } => format!(
            "Fragment order: {}\nJunction overlaps (bp): {}\nProduct topology: {}",
            fragments.iter().map(source).collect::<Vec<_>>().join(" → "),
            overlaps
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(", "),
            if *circular { "circular" } else { "linear" }
        ),
    }
}
