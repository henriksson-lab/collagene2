//! Checked protocol output owns fresh verification and immutable source bindings.
use super::{
    engine::verify,
    evidence::{PhysicalProcessingEvidence, ReferenceStepEvidence, VerificationReport},
};
use crate::{
    chemistry::{self, tracked_duplex as modified_duplex},
    planning::plan::{DesignPlan, Operation, Source},
};
use collagene_model::{Alphabet, Error, Molecule};
pub fn escape_html(value: &str) -> String {
    value
        .chars()
        .map(|c| match c {
            '&' => "&amp;".into(),
            '<' => "&lt;".into(),
            '>' => "&gt;".into(),
            '"' => "&quot;".into(),
            '\'' => "&#39;".into(),
            _ => c.to_string(),
        })
        .collect()
}

/// Resolve only a validated bound input or an actually reconstructed preceding
/// step. Invalid/forward/missing references must not invent ordered oligo bases.
pub(crate) fn protocol_oligo_sequence<'a>(
    plan: &'a DesignPlan,
    report: &'a VerificationReport,
    source: &Source,
    step_index: usize,
) -> Option<&'a str> {
    let molecule = match source {
        Source::Input(id) => {
            let mut matches = plan.inputs.iter().filter(|input| &input.id == id);
            let input = matches.next()?;
            if matches.next().is_some() || input.revision != input.molecule.revision() {
                return None;
            }
            &input.molecule
        }
        Source::Step(index) if *index < step_index => report.product(*index)?,
        Source::Step(_) => return None,
    };
    if molecule.sequence.alphabet() != Alphabet::Dna
        || molecule.sequence.len() > chemistry::MAX_ANNEAL_OLIGO_BASES
    {
        return None;
    }
    Some(molecule.sequence.as_str())
}

/// Native-generated documentary appendix. It cannot supply molecular status.
pub(crate) struct ProtocolAppendix {
    pub paragraphs: Vec<String>,
    pub tables: Vec<ProtocolTable>,
}
pub(crate) struct ProtocolTable {
    pub heading: String,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}
impl ProtocolAppendix {
    pub(crate) fn validate(&self) -> Result<(), Error> {
        let mut bytes = 0usize;
        let mut rows = 0usize;
        if self.paragraphs.len() > 8192 || self.tables.len() > 128 {
            return Err(Error::TooLarge);
        }
        let mut charge = |text: &str| -> Result<(), Error> {
            bytes = bytes
                .checked_add(text.len())
                .filter(|n| *n <= 512 * 1024)
                .ok_or(Error::TooLarge)?;
            Ok(())
        };
        for p in &self.paragraphs {
            charge(p)?;
        }
        for table in &self.tables {
            rows = rows
                .checked_add(table.rows.len())
                .filter(|n| *n <= 8192)
                .ok_or(Error::TooLarge)?;
            if table.columns.is_empty()
                || table.columns.len() > 8
                || table.rows.iter().any(|r| r.len() != table.columns.len())
            {
                return Err(Error::Format("invalid protocol appendix table".into()));
            }
            charge(&table.heading)?;
            for cell in table.columns.iter().chain(table.rows.iter().flatten()) {
                charge(cell)?;
            }
        }
        Ok(())
    }
}
/// Fresh in-process protocol context. It cannot be deserialized or constructed
/// from an imported report; immutable borrows keep plan/target paired with this
/// exact native reconstruction while several output formats are rendered.
pub(crate) struct CheckedProtocol<'a> {
    plan: &'a DesignPlan,
    target: &'a Molecule,
    report: VerificationReport,
    appendix: Option<&'a ProtocolAppendix>,
}
impl<'a> CheckedProtocol<'a> {
    pub(crate) fn with_appendix(mut self, appendix: &'a ProtocolAppendix) -> Result<Self, Error> {
        appendix.validate()?;
        self.appendix = Some(appendix);
        Ok(self)
    }
    pub(crate) fn appendix(&self) -> Option<&ProtocolAppendix> {
        self.appendix
    }

    pub(crate) fn plan(&self) -> &'a DesignPlan {
        self.plan
    }
    pub(crate) fn target(&self) -> &'a Molecule {
        self.target
    }
    pub(crate) fn report(&self) -> &VerificationReport {
        &self.report
    }
    pub(crate) fn new(plan: &'a DesignPlan, target: &'a Molecule) -> Self {
        // Native verify owns the source/operation guards and independently rebuilds
        // the product. Invalid inputs produce fresh draft findings, never authority.
        CheckedProtocol {
            plan,
            target,
            report: verify(plan, target),
            appendix: None,
        }
    }
    pub(crate) fn html(&self) -> String {
        let plan = self.plan();
        let target = self.target();
        let report = self.report();
        let e = escape_html;
        let mut body = format!(
            "<h1>{}</h1><p><strong>{:?}</strong> · computed sequence checks; laboratory conditions require review.</p><p>Target: {} · {} bases · {:?}</p><h2>Starting materials</h2><ul>",
            e(&plan.name),
            report.status,
            e(&target.name),
            target.sequence.len(),
            target.topology
        );
        for input in &plan.inputs {
            body.push_str(&format!(
            "<li>Input <strong>{}</strong>: {} · {} bases · {:?}<br>Revision <code>{}</code></li>",
            e(&input.id),
            e(&input.molecule.name),
            input.molecule.sequence.len(),
            input.molecule.topology,
            e(&input.revision)
        ));
        }
        body.push_str("</ul>");
        match crate::output::standard_review::lines(plan, report) {
            Ok(lines) if !lines.is_empty() => {
                body.push_str("<h2>Declared source-part standard</h2>");
                for line in lines {
                    body.push_str(&format!("<p>{}</p>", e(&line)));
                }
            }
            Ok(_) => (),
            Err(error) => body.push_str(&format!(
                "<p>Source-format advice unavailable: {}</p>",
                e(&error.to_string())
            )),
        }
        body.push_str("<h2>Cloning steps</h2><ol>");
        let source = |s: &Source| {
            e(&match s {
                Source::Input(id) => format!("input {id}"),
                Source::Step(i) => i
                    .checked_add(1)
                    .map(|number| format!("product of step {number}"))
                    .unwrap_or_else(|| format!("invalid step index {i}")),
            })
        };
        for (step_index, step) in plan.steps.iter().enumerate() {
            match step {
            Operation::ResectThreePrime { source: operand, cuts } => {
                body.push_str(&format!("<li>Explicit3′→5′ resection of {}: remove {} nt from bottom3′ left end and {} nt from top3′ right end. Positive cuts require accessible3′OH; original5′ groups and surviving modifications stay unchanged. No fill, ligase, unblocking or host repair.</li>",source(operand),cuts.left,cuts.right));
            },
            Operation::HomologyAnnealing {fragments,options,policy} => {
                body.push_str(&format!("<li>Generic3′ processing and exposed5′tail annealing of ordered substrates {}. Explicit policy <code>{}</code>; advisory profile <code>{}</code>. Every nick/gap retains separate covalent strands. The reference contour does not establish a reusable assembled molecule or proprietary In-Fusion/cellular rescue.</li>",fragments.iter().map(|f|format!("{} ({})",source(&f.source),if f.reverse{"reverse complement"}else{"forward"})).collect::<Vec<_>>().join(", "),e(&serde_json::to_string(policy).unwrap_or_default()),e(&serde_json::to_string(options).unwrap_or_default())));
                if let Some(rescue) = &policy.planned_rescue {
                    body.push_str(&format!("<li>Planned next step: {}. This is a documentary intention; transformation, host repair, isolation and experimental success are not simulated. The joined reference remains incomplete and is not reusable physical material.</li>", e(rescue.label())));
                }
            },
            Operation::HomologyAssembly { fragments, options, policy } => {
                body.push_str(&format!("<li>Physical homology assembly of ordered substrates {}. Declared geometry and 5′ processing/fill/sealing policy: <code>{}</code>. Advisory sequence profile: <code>{}</code>. Surviving source strands retain modifications; only actual template-supported synthesis is unmodified. Eligible paired 3′OH/5′P nicks may seal; unresolved chemistry stays explicit. Linear outer-end processing is retained, not silently blunted. This model does not establish experimental success or a proprietary formulation.</li>",fragments.iter().map(|f|format!("{} ({})",source(&f.source),if f.reverse{"reverse complement"}else{"forward"})).collect::<Vec<_>>().join(", "),e(&serde_json::to_string(policy).unwrap_or_default()),e(&serde_json::to_string(options).unwrap_or_default())));
            },
            Operation::OverlapExtension { left, right, options } => {
                body.push_str(&format!("<li>Denature overlapping duplexes {} ({}) and {} ({}), reanneal left top/right bottom strands across {} bp, and extend both paired 3′ ends. Minimum terminal pairing {} bp. Complementary synthesis is assumed unmodified; selected geometry does not establish isolation or experimental success.</li>", source(&left.source), if left.reverse {"reverse complement"} else {"forward"}, source(&right.source), if right.reverse {"reverse complement"} else {"forward"}, options.overlap, options.min_paired_bases));
            },
            Operation::ReconstructMaterial { template, origin } => {
                body.push_str(&format!("<li>Retrieve saved material {} by replaying its exact source and native derivation for event <code>{}</code>. Preserve actual strands, termini, modifications and inherited uncertainty. This retrieves material; it does not rerun historical reactions, create a target assessment or prove isolation.</li>",source(template),e(&origin.product.event_id)));
            },
            Operation::ReconstructRunProduct { template, origin } => {
                body.push_str(&format!("<li>Use saved material {} by reconstructing run <code>{}</code>, product {} against its unchanged original target <code>{}</code>. Preserve actual nicks, termini and positional modifications. Complete origin evidence is retained; this is not a newly performed reaction, inferred repair or proof of isolation.</li>", source(template), e(&origin.selection.run_id), origin.selection.step_index+1, e(&origin.target_event.revision)));
            },
            Operation::TopoActions{sources,request}=>{
                body.push_str(&format!("<li>Accounted TOPO chemical actions from original roles {}. Every original base and excised fragment remains in the full graph. Explicit selected component {} and source seed strand {}, index {} (zero-based); component selection is not isolation. Ideal hydrolytic profiles are declared model assumptions, not a named nuclease, kit or host-repair guarantee.<ul>",sources.iter().map(|r|format!("{} = {}",e(&r.role),source(&r.source))).collect::<Vec<_>>().join(", "),request.selection.component,request.selection.reference_seed.strand,request.selection.reference_seed.index));
                body.push_str(&topo_recipe_html(&request.preparation));
                body.push_str("</ul></li>");
            },

            Operation::Topo { sources, request, selection } | Operation::CheckedTopo { sources, request, selection, .. } => {
                body.push_str(&format!("<li>Source-backed TOPO elementary-transfer graph with explicit ordered roles: {}. Current native model limits <code>{}</code>; selection <code>{}</code>. Full physical evidence is retained; a paired-contour reference view does not establish nick repair, flap removal, supplier state, target acceptance or experimental success.<ul>",sources.iter().map(|s|format!("{} = {}",e(&s.role),source(&s.source))).collect::<Vec<_>>().join(", "),e(&serde_json::to_string(&request.engine_limits).unwrap_or_default()),e(&serde_json::to_string(selection).unwrap_or_default())));
                if let Operation::CheckedTopo { policy, .. } = step {
                    body.push_str(&format!("<li>Final-only independent-target physical policy: <code>{}</code>. Acceptance is conditional on this declared graph and fully specified DNA. Rotation/reverse comparison is allowed only as explicitly listed; source-coordinate annotation correspondence does not establish chemical ancestry, graph repair or experimental success.</li>", e(&serde_json::to_string(policy).unwrap_or_default())));
                }
                body.push_str(&topo_recipe_html(request));
                body.push_str("</ul></li>");
            },
            Operation::RecombineReference { sources, request, selection } => body.push_str(&format!("<li>Reference-only {:?} exchange of {} and {}. Pinned full profile/data, work/annotation limits: <code>{}</code>. Explicit component selection: <code>{}</code>. All reciprocal components are computed; selected reference DNA does not establish chemical state, efficacy, repair or experimental selection.</li>",request.reaction,source(&sources[0]),source(&sources[1]),e(&serde_json::to_string(request).unwrap_or_default()),e(&serde_json::to_string(selection).unwrap_or_default()))),
            Operation::PcrWithMismatches{template,forward,reverse,options}=>body.push_str(&format!("<li>Intentional-mismatch PCR of {}. Oligos 5′→3′: <code>{}</code> and <code>{}</code>. Binding/search policy: <code>{}</code>. Primer-derived changes are reconstructed; amplification is not established.</li>",source(template),e(&forward.ordered_sequence()),e(&reverse.ordered_sequence()),e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::GoldenGate{templates,options,candidate_index}=>body.push_str(&format!("<li>Golden Gate with {:?}, inputs {}. Explicit zero-based candidate {} under policy <code>{}</code>. Cleavage, orientation, cohesive ends and recognition-site removal are checked. Finite one-copy fragment search does not establish exhaustive mixture uniqueness; conditions and quantities are not specified by the molecular model.</li>",options.enzyme,templates.iter().map(source).collect::<Vec<_>>().join(", "),candidate_index,e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::SelectedGoldenGate{templates,options,candidate_index,directed_sequence_key,physical_sequence_key}=>body.push_str(&format!("<li>Key-bound Golden Gate with {:?}, inputs {}. Explicit zero-based candidate {} under policy <code>{}</code>. Cleavage, orientation, cohesive ends and recognition-site removal are checked. Finite one-copy fragment search does not establish exhaustive mixture uniqueness; conditions and quantities are not specified by the molecular model. Fresh enumeration must match directed key <code>{}</code> and physical key <code>{}</code>.</li>",options.enzyme,templates.iter().map(source).collect::<Vec<_>>().join(", "),candidate_index,e(&serde_json::to_string(options).unwrap_or_default()),e(directed_sequence_key),e(physical_sequence_key))),
            Operation::MultiSiteRecombineReference{sources,request,selection}=>body.push_str(&format!("<li>Declared multi-site {:?} reference exchange using sources {}. Profile <code>{}</code>; specificity order {:?}; select candidate {}, component {} from bound batch <code>{}</code>. Retain all reciprocal components. This is a reference design, not established chemical material, experimental selection or kit success.</li>",request.reaction.reaction,e(&sources.iter().map(|s|format!("{}: {}",s.role,source(&s.source))).collect::<Vec<_>>().join(", ")),e(&request.reaction.profile.id),request.reaction.specificities,selection.candidate_id,selection.component_id,e(&selection.batch_revision))),
            Operation::Pcr { template, forward, reverse, template_strand } => body.push_str(&format!("<li>PCR of {} in {:?} strand coordinates with {}: <code>{}</code> and {}: <code>{}</code>. Primers shown 5′→3′. The molecular model does not specify reaction conditions or quantities.</li>",source(template),template_strand,e(&forward.name),e(&forward.ordered_sequence()),e(&reverse.name),e(&reverse.ordered_sequence()))),
            Operation::DeclareSourceTermini{template,left,right}=>body.push_str(&format!("<li>Declare original intact linear duplex source termini for {}: left <code>{}</code>, right <code>{}</code>. Explicit operator hypothesis only; not measured chemistry, synthesis, repair, recovered imported physical state or experimental success. Original sequence/features/citations are retained; methylation remains unknown.</li>",source(template),e(&serde_json::to_string(left).unwrap_or_default()),e(&serde_json::to_string(right).unwrap_or_default()))),
            Operation::DeclareSourceModifications{template,assumptions}=>body.push_str(&format!("<li>Declare closed Dam/Dcm/CpG source generation for {} under explicit assumptions <code>{}</code>. Bind the original input revision; these are not measured modifications. Descendant strand states are propagated without resynthesizing methylation.</li>",source(template),e(&serde_json::to_string(assumptions).unwrap_or_default()))),
            Operation::CheckedDigest{template,selection,options,fragment_index,..}=>body.push_str(&format!("<li>Checked positional-modification digest of {} with {}. Profile <code>{}</code>; zero-based fragment {}. Analysis limits: {} comparisons, {} sites, {} fragments, {} report bytes. Blocked sites remain intact; unknown/impaired or incomplete evidence withholds fully determined products. Reaction conditions and experimental success remain unspecified.</li>",source(template),e(&selection.enzymes.iter().map(|e|e.enzyme.name.clone()).collect::<Vec<_>>().join(", ")),e(&modified_duplex::profile_revision()),fragment_index,options.max_comparisons,options.max_sites,options.max_fragments,options.max_report_bytes)),
            Operation::CatalogDigest{template,selection,fragment_index,..}=>body.push_str(&format!("<li>Catalog-pinned complete digest of {} under exact selection <code>{}</code>; select zero-based fragment {} in top-reference cut ordering. Identity/name/geometry are checked against the embedded catalog. Complete cleavage and compatible conditions remain assumptions; supplier observations do not establish kinetics, methylation compatibility or experimental completion.</li>",source(template),e(&serde_json::to_string(selection).unwrap_or_default()),fragment_index)),
            Operation::Digest{template,enzymes,fragment_index}=>body.push_str(&format!("<li>Complete restriction digest of {} with {}. Select zero-based fragment index {} in ascending top-reference cut order. Cleavage assumes compatible conditions and no blocking methylation.</li>",source(template),e(&enzymes.join(", ")),fragment_index)),
            Operation::TreatEnds{template,treatment}=>body.push_str(&format!("<li>Terminal {:?} treatment of {}. Internal nick repair is not assumed; the molecular model does not specify reaction conditions or quantities.</li>",treatment,source(template))),
            Operation::Tail3Prime{template,base}=>body.push_str(&format!("<li>Append one {} nucleotide to both strand 3′ ends of {}. Requires reconstructed blunt linear DNA with accessible 3′-OH. Five-prime chemistry and internal nicks remain unchanged; no template ancestry or polymerase/kit activity is inferred.</li>",base.as_str(),source(template))),
            Operation::ProfiledOverlap{fragments,overlaps,circular,options}=>body.push_str(&format!("<li>Profiled {:?} ordered overlap assembly: {}. Overlaps {:?} bp; {} product. Explicit source orientation is retained. Method recommendations are advisory, not evidence of reagent compatibility, experimental completion or uniqueness. Policy <code>{}</code>.</li>",options.method,fragments.iter().map(|f|format!("{} ({})",source(&f.source),if f.reverse{"reverse complement"}else{"forward"})).collect::<Vec<_>>().join(" → "),overlaps,if *circular{"circular"}else{"linear"},e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::RepairEnds{template,method}=>body.push_str(&format!("<li>End repair {:?} of {}. Only explicitly modeled 3′ fill-in/exonuclease changes occur. Five-prime chemistry and internal nicks are retained; polymerase activity, reagents and experimental success are assumptions.</li>",method,source(template))),
            Operation::AnnealOligos{top,bottom,options,candidate_index}=>body.push_str(&format!("<li>Anneal declared 5′→3′ ssDNA inputs {} and {}.<br>Top ordered oligo 5′→3′: <code>{}</code> (5′ {:?}; 3′ {:?}).<br>Bottom ordered oligo 5′→3′: <code>{}</code> (5′ {:?}; 3′ {:?}).<br>Select zero-based candidate {} under policy <code>{}</code>. Exact antiparallel pairing and placement alternatives are reconstructed; no covalent ligation or implicit phosphorylation occurs.</li>",source(top),source(bottom),e(protocol_oligo_sequence(plan,&report,top,step_index).unwrap_or("Unresolved or unsupported bound oligo sequence")),options.top_chemistry.five_prime,options.top_chemistry.three_prime,e(protocol_oligo_sequence(plan,&report,bottom,step_index).unwrap_or("Unresolved or unsupported bound oligo sequence")),options.bottom_chemistry.five_prime,options.bottom_chemistry.three_prime,candidate_index,e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::Ligate{fragments,circular,allow_nicked,options}=>body.push_str(&format!("<li>Intended ordered ligation: {}. Product {}; retained nicks {}. Enumerate one-copy alternatives under explicit bounds <code>{}</code>; retain the intended path independently. No probabilities or experimental isolation are inferred. Full extractable outcome paths/sequences are retained below.</li>",fragments.iter().map(|f|format!("{} ({})",source(&f.source),if f.reverse{"reverse"}else{"forward"})).collect::<Vec<_>>().join(", "),if *circular{"circular"}else{"linear"},if *allow_nicked{"allowed"}else{"not allowed"},e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::Overlap { fragments, overlaps, circular } => body.push_str(&format!("<li>Ordered overlap assembly, in this order: {}. Junction overlaps {:?} bp; {} product. The molecular model does not specify reagent conditions or quantities.</li>",fragments.iter().map(source).collect::<Vec<_>>().join(", "),overlaps,if *circular {"circular"} else {"linear"})),
        }
            if matches!(
                step,
                Operation::CatalogDigest { .. } | Operation::CheckedDigest { .. }
            ) {
                if let Some(evidence) = report
                    .step_reports
                    .iter()
                    .find(|entry| entry.step_index == step_index)
                    .and_then(|entry| entry.reaction_conditions.as_ref())
                {
                    match evidence.protocol_lines() {
                        Ok(lines) => { body.push_str("<li><ul>"); for line in lines { body.push_str(&format!("<li>{}</li>", e(&line))); } body.push_str("</ul></li>"); }
                        Err(_) => body.push_str("<li>Documentary reaction-condition evidence failed its retained integrity check.</li>"),
                    }
                } else {
                    body.push_str("<li>Documentary reaction conditions unassessed: no explicit formulation/buffer/temperature/time choices were declared. Complete cleavage remains an explicit molecular assumption.</li>");
                }
            }
            if matches!(step, Operation::Digest { .. }) {
                body.push_str("<li>Documentary reaction conditions unsupported for uncatalogued Digest: use an exact catalog selection to retain formulation-specific choices.</li>");
            }
            if matches!(
                step,
                Operation::GoldenGate { .. } | Operation::SelectedGoldenGate { .. }
            ) {
                body.push_str("<li>Golden Gate mixture documentary conditions unassessed: supplier formulations, buffer, temperature schedule, ligase, duration and star-activity compatibility are not assessed.</li>");
            }
        }
        body.push_str("</ol><h2>Reconstructed products</h2><ul>");
        for outcome in &report.steps {
            let index = outcome.step_index;
            let Some(product) = &outcome.product else {
                body.push_str(&format!(
                    "<li>Step {}: {:?}; direct blockers <code>{}</code></li>",
                    index + 1,
                    outcome.state,
                    e(&serde_json::to_string(&outcome.blocked_by).unwrap_or_default())
                ));
                continue;
            };
            body.push_str(&format!(
                "<li>Step {}: {} · {} bases · {:?}</li>",
                index + 1,
                e(&product.name),
                product.sequence.len(),
                product.topology
            ));
            if let Some(summary) = product
                .metadata
                .get("duplex_model")
                .and_then(|s| serde_json::from_str::<chemistry::DuplexSummary>(s).ok())
            {
                body.push_str(&format!(
                    "<li><pre style=\"white-space:pre-wrap\">{}</pre></li>",
                    e(&summary.describe())
                ));
            }
        }
        body.push_str("</ul><h2>Intended target map</h2>");
        match crate::output::svg::render_map(target) {
            Ok(svg) => body.push_str(&format!("<div class=\"protocol-map\">{svg}</div>")),
            Err(error) => body.push_str(&format!(
                "<p>Map unavailable: {}</p>",
                e(&error.to_string())
            )),
        }
        body.push_str("<h2>Verification findings</h2><ul>");
        for finding in &report.findings {
            body.push_str(&format!("<li>{}<ul>", e(&finding.summary())));
            for citation in &finding.citations {
                body.push_str(&format!(
                    "<li>{} · {}</li>",
                    e(&citation.source),
                    e(&citation.locator)
                ));
            }
            body.push_str("</ul></li>");
        }
        body.push_str("</ul>");
        if report.step_reports.iter().any(|step| step.topo.is_some()) {
            body.push_str("<h2>TOPO physical step evidence</h2><ul>");
            for step in &report.step_reports {
                if let Some(topo) = &step.topo {
                    let physical = &topo.projection.physical;
                    let graph = &physical.product;
                    let unpaired = graph
                        .unpaired
                        .iter()
                        .map(|run| run.end - run.start)
                        .sum::<usize>();
                    let contour = match &topo.projection.outcome {
                        crate::cloning::topo::physical::projection::Outcome::Proven(view) => {
                            format!(
                                "Proven {:?} paired contour with {} retained nicks; known terminal chemistry {}; all requested transfers applied {}. This is no target-acceptance or repair claim.",
                                view.topology,
                                view.nicks.len(),
                                view.fully_known_terminal_chemistry,
                                view.all_requested_steps_applied
                            )
                        }
                        crate::cloning::topo::physical::projection::Outcome::Withheld(reason) => {
                            format!(
                                "Contour withheld: {reason:?}. No target-derived path, flap removal or nick repair was attempted."
                            )
                        }
                    };
                    body.push_str(&format!("<li>Step {}: {} components, {} actual chains, {} graph nicks and {} unpaired nucleotides retained. {} The full bounded graph and source-call maps remain in the downloadable verification JSON.<ul>",step.step_index+1,graph.components.len(),graph.chains.len(),graph.nicks.len(),unpaired,e(&contour)));
                    body.push_str(&format!("<li>{} retained source steps; {} newly requested steps. Native source origins: {}. A consumed activation is never made available again.</li>",topo.preparation.retained_steps.len(),physical.steps.len().saturating_sub(topo.preparation.retained_steps.len()),topo.preparation.graph_origins.len()));
                    for (i, chain) in graph.chains.iter().enumerate() {
                        body.push_str(&format!("<li>Actual chain {}5′→3′: <code>{}</code>; {}; source-call mappings {:?}.</li>",i,e(chain.sequence.as_str()),e(&format!("{:?}",chain.ends)),chain.sources));
                    }
                    for nick in &graph.nicks {
                        body.push_str(&format!(
                        "<li>Actual unsealed nick: source {:?}3′ {} to source {:?}5′ {:?}.</li>",
                        nick.three,
                        e(&format!("{:?}", nick.three_prime)),
                        nick.five,
                        nick.five_prime
                    ));
                    }
                    if let Some(annotations) = &topo.annotations {
                        body.push_str(&format!("<li>{} projected source annotations; correspondence and review evidence retained in this canonical product, without inherited-function or chemistry claims.</li>",annotations.features.len()));
                    }
                    for assumption in &topo.preparation.assumptions {
                        body.push_str(&format!("<li>{}</li>", e(assumption)));
                    }
                    body.push_str("</ul></li>");
                }
            }
            body.push_str("</ul>");
        }
        if !report.reference_step_reports.is_empty() {
            body.push_str("<h2>Whole reference-network evidence</h2><ul>");
            for step in &report.reference_step_reports {
                let ReferenceStepEvidence::MultiSiteRecombination { selection, batch } =
                    &step.evidence;
                body.push_str(&format!("<li>Step {}: {} original parents; {} candidate networks; {} reciprocal components retained. Complete within declared model: {}. Selected candidate {}, component {}. No target determined a pairing; reference agreement does not establish physical chemistry.<ul>",step.step_index+1,batch.sources.len(),batch.reference.candidates.len(),batch.documents.len(),batch.reference.complete,selection.candidate_id+1,selection.component_id+1));
                for note in batch
                    .reference
                    .warnings
                    .iter()
                    .chain(&batch.reference.limits)
                    .chain(&batch.reference.unsupported)
                    .chain(&batch.warnings)
                {
                    body.push_str(&format!("<li>{}</li>", e(note)));
                }
                body.push_str("</ul></li>");
            }
            body.push_str("</ul>");
        }
        if !report.physical_processing_reports.is_empty() {
            body.push_str("<h2>Accounted physical processing evidence</h2><ul>");
            for step in &report.physical_processing_reports {
                match &step.evidence {
                    PhysicalProcessingEvidence::LigationOutcomes(evidence) => {
                        body.push_str(&format!("<li>Step {} ligation outcomes: intended path retained; {} alternatives, {} checked paths, {} incompatible and {} unresolved. Enumeration complete {}, exhaustive within one-copy scope {}. Limits: {}. These are possible native paths, not probabilities, a complete mixture model or proof of isolation.<ol>",step.step_index+1,evidence.report.alternatives.len(),evidence.report.paths_checked,evidence.report.incompatible_paths,evidence.report.unresolved_paths,evidence.report.enumeration_complete,evidence.report.exhaustive_within_scope,e(&format!("{:?}",evidence.report.limits))));
                        for (i, product) in std::iter::once(&evidence.report.intended)
                            .chain(&evidence.report.alternatives)
                            .enumerate()
                        {
                            body.push_str(&format!("<li>{}: {} · {:?}, {} bases · material <code>{}</code>; {} retained nicks.</li>",if i==0{"Intended".into()}else{format!("Alternative {i}")},e(&product.path.iter().map(|p|format!("{} ({})",p.key,if p.reverse{"reverse"}else{"forward"})).collect::<Vec<_>>().join(" → ")),product.topology,product.sequence.len(),e(&product.modifications.state_revision),product.duplex.nicks.len()));
                        }
                        body.push_str(
                            "</ol><details><summary>Full extractable outcome FASTA</summary><pre>",
                        );
                        body.push_str(&e(&evidence.report.fasta().unwrap_or_else(|error| {
                            format!("Outcome extraction refused: {error}")
                        })));
                        body.push_str("</pre></details><details><summary>Full typed outcome paths, strands, modifications and bounds (JSON)</summary><pre>");
                        body.push_str(&e(&serde_json::to_string_pretty(evidence)
                            .unwrap_or_else(|error| format!("Outcome JSON refused: {error}"))));
                        body.push_str("</pre></details></li>");
                    }
                    PhysicalProcessingEvidence::ThreePrimeResection(geometry) => {
                        body.push_str(&format!("<li>Step {} actual3′ resection: bottom left {} nt / top right {} nt; retained top {} nt, bottom {} nt at offset {}. Exact copied/removed source spans and unchanged5′ groups remain in JSON evidence.</li>",step.step_index+1,geometry.cuts.left,geometry.cuts.right,geometry.product.top.len(),geometry.product.bottom.len(),geometry.product.bottom_start));
                    }
                    PhysicalProcessingEvidence::HomologyAnnealing(geometry) => {
                        body.push_str(&format!("<li>Step {} annealed diagnostic: {} actual surviving strands, {} junctions; paired gapless diagnostic {}; terminal search complete {}. No joining bond is sealed and no host repair is assumed.<ul>",step.step_index+1,geometry.segments.len(),geometry.junctions.len(),geometry.paired_diagnostic.is_some(),geometry.terminal_search_complete));
                        for j in &geometry.junctions {
                            body.push_str(&format!("<li>Sources {} → {}: {} bp homology; top discontinuity [{}, {}) 3′ {:?}/5′ {:?}; bottom [{}, {}) 5′ {:?}/3′ {:?}. Equal boundaries mean an unsealed nick; nonempty intervals remain unfilled gaps.</li>",j.left_source+1,j.right_source+1,j.overlap,j.top.start,j.top.end,j.top.three_prime,j.top.five_prime,j.bottom.start,j.bottom.end,j.bottom.five_prime,j.bottom.three_prime));
                        }
                        body.push_str("</ul></li>");
                    }
                    PhysicalProcessingEvidence::HomologyAssembly(geometry) => {
                        body.push_str(&format!("<li>Step {} physical homology assembly: {} junctions, {} resected strand spans; paired diagnostic product available {}; incomplete {}; terminal search complete {}; {} comparisons.<ul>",step.step_index+1,geometry.junctions.len(),geometry.resected.len(),geometry.product.is_some(),geometry.incomplete,geometry.terminal_search_complete,geometry.comparisons));
                        for junction in &geometry.junctions {
                            body.push_str(&format!("<li>Sources {} → {}: {} bp homology; gaps top/bottom {}/{}; seals top/bottom {}/{}; top nick at {} and bottom nick at {} (zero-based product boundaries).</li>",junction.left_source+1,junction.right_source+1,junction.overlap,junction.top_gap,junction.bottom_gap,junction.top_sealed,junction.bottom_sealed,junction.top_nick,junction.bottom_nick));
                        }
                        for note in geometry.warnings.iter().chain(&geometry.limits) {
                            body.push_str(&format!("<li>{}</li>", e(note)));
                        }
                        body.push_str("</ul></li>");
                    }
                    PhysicalProcessingEvidence::OverlapExtension(geometry) => {
                        body.push_str(&format!("<li>Step {} reciprocal extension: selected homology {} bp; retained terminal pairings {:?}; complete search {}; unique within terminal model {}; {} comparisons. Original source geometry revisions: {} / {}. Exact copied/synthesized strand mappings remain in the reconstruction report.</li>", step.step_index+1, geometry.options.overlap, geometry.terminal_pairings.overlaps, geometry.terminal_pairings.complete, geometry.terminal_pairings.unique, geometry.terminal_pairings.comparisons, e(&geometry.source_geometry_revisions[0]), e(&geometry.source_geometry_revisions[1])));
                    }
                    PhysicalProcessingEvidence::DerivedOrigin(evidence) => {
                        body.push_str("<li>Origin report documentary reaction-condition provenance: these retained choices do not declare current reaction conditions or establish compatibility.<ul>");
                        for conditions in &evidence.reaction_conditions {
                            match conditions.protocol_lines() {
                                Ok(lines) => for line in lines { body.push_str(&format!("<li>{}</li>", e(&line))); },
                                Err(_) => body.push_str("<li>Retained documentary conditions failed integrity validation.</li>"),
                            }
                        }
                        body.push_str("</ul></li>");
                        let material = &evidence.material;
                        body.push_str(&format!("<li>Step {}: replayed derivation event <code>{}</code>; material <code>{}</code>. {} Origin incomplete: {}; showing {} of {} findings{}.<ul>",step.step_index+1,e(&evidence.origin.product.event_id),e(&material.material_revision),e(&material.duplex.describe()),material.origin_incomplete,material.origin_findings.len(),material.origin_findings_total,if material.origin_findings_truncated {" (truncated)"} else {""}));
                        for finding in &material.origin_findings {
                            body.push_str(&format!("<li>{}</li>", e(finding)));
                        }
                        body.push_str("</ul></li>");
                    }
                    PhysicalProcessingEvidence::RunProduct(evidence) => {
                        body.push_str(&format!("<li>Step {}: saved run <code>{}</code>, product {}, actual material summary <code>{}</code>. Full originating report is retained with the reconstruction recipe.<ul>", step.step_index+1, e(&evidence.origin.selection.run_id), evidence.origin.selection.step_index+1, e(&evidence.material_revision)));
                        if let Ok(origin_report) = evidence.decode_report() {
                            body.push_str("<li>Origin report documentary reaction-condition provenance: these retained choices do not declare current reaction conditions or establish compatibility.<ul>");
                            match crate::material::derivation::reaction_conditions::from_origin_report(&evidence.origin, &origin_report) {
                                Ok(records) => for conditions in records {
                                    match conditions.protocol_lines() {
                                        Ok(lines) => for line in lines { body.push_str(&format!("<li>{}</li>", e(&line))); },
                                        Err(_) => body.push_str("<li>Retained documentary conditions failed integrity validation.</li>"),
                                    }
                                },
                                Err(_) => body.push_str("<li>Origin documentary reaction conditions could not be retained within current bounds.</li>"),
                            }
                            body.push_str("</ul></li>");
                            body.push_str(&format!("<li>Origin assessment: {:?}; engine {}. Reuse does not upgrade incomplete findings.</li>", origin_report.status,e(&origin_report.engine_version)));
                            for finding in &origin_report.findings {
                                body.push_str(&format!("<li>{}</li>", e(&finding.summary())));
                            }
                            for assumption in &origin_report.assumptions {
                                body.push_str(&format!(
                                    "<li>Origin assumption: {}</li>",
                                    e(assumption)
                                ));
                            }
                        }
                        body.push_str("</ul></li>");
                    }
                    PhysicalProcessingEvidence::TopoActions(report) => {
                        let graph = &report.projection.physical.product;
                        body.push_str(&format!("<li>Step {}: {} original source calls remain across {} components, {} chains and {} unpaired calls. Faithful selected duplex available: {}. Excised material is retained, not discarded from evidence.<ul>",step.step_index+1,report.projection.physical.source_bases,graph.components.len(),graph.chains.len(),graph.unpaired.iter().map(|r|r.end-r.start).sum::<usize>(),report.tracked.is_some()));
                        for action in &report.projection.physical.steps {
                            body.push_str(&format!("<li>{}</li>", e(&format!("{action:?}"))));
                        }
                        for note in &report.bridge_notes {
                            body.push_str(&format!("<li>{}</li>", e(note)));
                        }
                        body.push_str("</ul></li>");
                    }
                }
            }
            body.push_str("</ul>");
        }

        let annotations = report
            .step_reports
            .iter()
            .filter_map(|step| {
                step.topo
                    .as_ref()
                    .and_then(|product| product.annotations.as_ref())
                    .map(|evidence| (step.step_index, evidence))
            })
            .collect::<Vec<_>>();
        if !annotations.is_empty() {
            body.push_str("<h2>TOPO annotation correspondence</h2><p>Projected source coordinates establish neither inherited chemical ancestry nor biological feature validity. Unmapped annotations remain original source evidence.</p><ul>");
            for (step_index, annotation) in annotations {
                let mut offset = 0usize;
                for role in &annotation.sources {
                    let end = offset.saturating_add(role.projected);
                    // The projector emits each role's features in source-role order.
                    let reviewed = annotation
                        .features
                        .get(offset..end)
                        .unwrap_or(&[])
                        .iter()
                        .filter(|feature| feature.needs_review)
                        .count();
                    body.push_str(&format!("<li>Step {}, role {}: {} projected ({} require review), {} original annotations retained unmapped, {} outside selected correspondence. Projection supported: {}.</li>",step_index+1,e(&role.role),role.projected,reviewed,role.unmapped.len(),role.dropped,role.projection_supported));
                    offset = end;
                }
            }
            body.push_str("</ul>");
        }
        body.push_str("<h2>Assumptions and limits</h2><ul>");
        for assumption in &report.assumptions {
            body.push_str(&format!("<li>{}</li>", e(assumption)));
        }
        body.push_str(&format!("</ul><h2>Provenance</h2><p>Engine {}<br>Plan <code>{}</code><br>Target <code>{}</code></p>",e(&report.engine_version),e(&report.plan_revision),e(&report.target_revision)));
        if let Some(appendix) = self.appendix() {
            body.push_str("<section><h2>Declared reaction quantities</h2><p>Arithmetic for the declared subset is separate from the molecular verification above. Unlisted conditions and reagents remain unspecified.</p>");
            for line in &appendix.paragraphs {
                body.push_str(&format!("<p>{}</p>", e(line)));
            }
            for table in &appendix.tables {
                body.push_str(&format!("<h3>{}</h3><table><thead><tr>", e(&table.heading)));
                for column in &table.columns {
                    body.push_str(&format!("<th>{}</th>", e(column)));
                }
                body.push_str("</tr></thead><tbody>");
                for row in &table.rows {
                    body.push_str("<tr>");
                    for cell in row {
                        body.push_str(&format!("<td>{}</td>", e(cell)));
                    }
                    body.push_str("</tr>");
                }
                body.push_str("</tbody></table>");
            }
            body.push_str("</section>");
        }
        format!(
            "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width\"><title>{}</title><style>body{{font:16px/1.6 system-ui;max-width:850px;margin:40px auto;padding:0 24px;color:#18302d}}code{{overflow-wrap:anywhere}}.protocol-map>svg{{display:block;max-width:100%;height:auto}}h1,h2{{line-height:1.2}}li{{margin:12px 0}}@media print{{body{{margin:0;max-width:none}}h2{{break-after:avoid}}li{{break-inside:avoid}}}}</style>{body}</html>",
            e(&plan.name)
        )
    }
}
/// Rechecks the plan at export time; the caller cannot attach an unrelated report.
fn topo_recipe_html(request: &crate::cloning::topo::Request) -> String {
    let roots = request.recipe.materials.iter().filter_map(|m| {
        if let crate::cloning::topo::MaterialRecipe::SavedGraph { origin, .. } = m {
            Some(crate::material::capsule::NativeOriginRef::TopoGraph(origin))
        } else {
            None
        }
    });
    if crate::material::capsule::native_origins_preflight(roots).is_err() {
        return "<li>Native graph source exceeds its current bounds or is invalid; no opaque payload is printed.</li>".into();
    }
    fn describe(r: &crate::cloning::topo::Request, depth: usize, out: &mut String) {
        if depth > 32 {
            return;
        }
        for m in &r.recipe.materials {
            match m{
   crate::cloning::topo::MaterialRecipe::Prepared{declaration_id,assumption,strands,..}=>{out.push_str(&format!("<li>Original preparation hypothesis {}: {}<ul>",escape_html(declaration_id),escape_html(assumption)));for s in strands{out.push_str(&format!("<li>Source strand {}:5′ {:?},3′ {}; original intervals {}.</li>",escape_html(&s.id),s.five_prime,escape_html(&format!("{:?}",s.three_prime)),s.segments.iter().map(|p|format!("{}:{}+{} {}",escape_html(&p.source_role),p.start,p.length,if p.reverse{"reverse"}else{"forward"})).collect::<Vec<_>>().join("; ")));}out.push_str("</ul></li>");},
   crate::cloning::topo::MaterialRecipe::ExactPcr{material_id,template_role,forward,reverse,tail}=>out.push_str(&format!("<li>New exact-PCR synthesis {} from {}: ordered5′→3′ oligos <code>{}</code> / <code>{}</code>; explicit3′ tail {:?}.</li>",escape_html(material_id),escape_html(template_role),escape_html(&forward.ordered_sequence()),escape_html(&reverse.ordered_sequence()),tail)),
   crate::cloning::topo::MaterialRecipe::SavedGraph{material_id,origin}=>{out.push_str(&format!("<li>Retrieve saved whole native graph {} (record <code>{}</code>, state <code>{}</code>), preserving all components and consumed activations. The following retained transcript is source history, not instructions to repeat historical PCR or transfers.<ul>",escape_html(material_id),escape_html(&origin.selection.record_hash),escape_html(&origin.selection.state_revision)));match origin.validate_shallow(){Ok(record)=>{describe(&record.request.request,depth+1,out);},Err(error)=>out.push_str(&format!("<li>Native graph unavailable: {}</li>",escape_html(&error.to_string())))}out.push_str("</ul></li>");}
  }
        }
        for (index, step) in r.steps.iter().enumerate() {
            use crate::cloning::topo::physical::{Step, processing::Action};
            let call = |c: crate::cloning::topo::physical::Call| {
                format!("source strand {} call {}", c.strand, c.index)
            };
            let text = match step {
                Step::Transfer(t) => format!(
                    "Transfer activated3′ {} to5′ {} under {:?} geometry; no implicit ligase or end repair.",
                    call(t.donor),
                    call(t.acceptor),
                    t.geometry
                ),
                Step::Process(Action::CleaveThreePrimeFlap {
                    bond,
                    displacement_step,
                    profile,
                    assumption,
                }) => format!(
                    "Cleave the3′ flap bond from {} to {}, displaced at transcript step {}, under {:?}. Declared conditions: {}. Retain the excised fragment; hydrolysis creates3′OH and5′phosphate, not TOPO activation.",
                    call(bond.upstream),
                    call(bond.downstream),
                    displacement_step + 1,
                    profile,
                    assumption
                ),
                Step::Process(Action::SealNick {
                    selection,
                    ligase_assumption,
                }) => format!(
                    "Seal the explicitly selected nick from3′ {} to5′ {}. Declared ligase conditions: {}. Fresh reconstruction checks exposed3′OH/5′phosphate and continuous template pairing; no other nick is repaired.",
                    call(selection.three),
                    call(selection.five),
                    ligase_assumption
                ),
            };
            out.push_str(&format!("<li>Reaction step {}: {} All source calls and components remain accounted; this declaration does not establish experimental isolation.</li>",index+1,escape_html(&text)));
        }
    }
    let mut out = String::new();
    describe(request, 0, &mut out);
    out
}
pub fn export_protocol_html(plan: &DesignPlan, target: &Molecule) -> String {
    CheckedProtocol::new(plan, target).html()
}
