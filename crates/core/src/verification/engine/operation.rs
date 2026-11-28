//! One staged native operation. Only the execution owner commits its product
//! and private tracked state; failed branches may retain explicitly scoped evidence.
use super::*;
use std::cell::Cell;

pub(super) struct MaterialBudget {
    pub product_bases: usize,
    pub product_bytes: usize,
    pub step_evidence_bytes: usize,
}
impl MaterialBudget {
    pub fn new() -> Self {
        Self {
            product_bases: 0,
            product_bytes: 0,
            step_evidence_bytes: 0,
        }
    }
}
pub(super) struct StepAttempt {
    pub result: Result<(Molecule, Option<TrackedDuplex>, bool), Error>,
    pub rule: RuleId,
    pub notes: Vec<String>,
    pub incomplete: Vec<String>,
    pub step_reports: Vec<StepReport>,
    pub physical_processing_reports: Vec<PhysicalProcessingStepReport>,
    pub reference_step_reports: Vec<ReferenceStepReport>,
    pub golden_gate: Option<crate::cloning::golden_gate::GoldenGateReport>,
}
impl StepAttempt {
    #[allow(clippy::too_many_arguments)]
    pub fn run(
        plan: &DesignPlan,
        target: &Molecule,
        index: usize,
        plan_hash: &str,
        inputs: &BTreeMap<String, &Molecule>,
        steps: &[StepOutcome],
        physical: &[Option<TrackedDuplex>],
        reference_only: &[bool],
        budget: &mut MaterialBudget,
        retained_processing: usize,
        cleavage_budget: &mut modified_duplex::CleavageBudget,
        pcr_budget: &mut analysis::ExactSearchBudget,
    ) -> Self {
        let step = &plan.steps[index];
        let rule = Cell::new(operation_rule(step));
        let mut notes = Vec::new();
        let findings = &mut notes;
        let mut incomplete = Vec::new();
        let mut step_reports = Vec::new();
        let mut physical_processing_reports = Vec::new();
        let mut reference_step_reports = Vec::new();
        let mut golden_gate = None;
        let final_golden_gate = &mut golden_gate;
        let mut product_bases = budget.product_bases;
        let mut product_bytes = budget.product_bytes;
        let mut step_evidence_bytes = budget.step_evidence_bytes;
        let result = (|| -> Result<(Molecule, Option<TrackedDuplex>, bool), Error> {
            rule.set(RuleId::TopoTargetPolicy);
            if matches!(step, Operation::CheckedTopo { .. }) && index + 1 != plan.steps.len() {
                return Err(Error::Unsupported("CheckedTopo independent-target acceptance is supported only as the final plan step; intermediate target acceptance and downstream physical graph consumption are not modeled.".into()));
            }
            rule.set(RuleId::ReferenceOnlyAuthority);
            let inherited_reference_only = step.sources().iter().any(|source| match source {
                Source::Step(i) => reference_only.get(*i).copied().unwrap_or(false),
                Source::Input(_) => false,
            });
            if inherited_reference_only
                && matches!(
                    step,
                    Operation::Pcr { .. }
                        | Operation::PcrWithMismatches { .. }
                        | Operation::AnnealOligos { .. }
                        | Operation::Topo { .. }
                        | Operation::CheckedTopo { .. }
                        | Operation::TopoActions { .. }
                )
            {
                return Err(Error::Unsupported("Reference-only ancestry cannot establish a physical PCR/annealing/TOPO substrate. Its full originating graph must be consumed by an explicit supported physical operation; no materialization or graph repair is inferred.".into()));
            }
            if inherited_reference_only
                && matches!(
                    step,
                    Operation::HomologyAssembly { .. }
                        | Operation::HomologyAnnealing { .. }
                        | Operation::ResectThreePrime { .. }
                )
            {
                return Err(Error::Unsupported("Reference-only ancestry cannot provide a physical homology substrate; reconstruct actual original strands first.".into()));
            }
            let resolve = |source: &Source| -> Result<&Molecule, Error> {
                match source {
                    Source::Input(id) => inputs.get(id).copied().ok_or_else(|| {
                        rule.set(RuleId::SourceDependency);
                        Error::Format(format!("Unknown input {id}"))
                    }),
                    Source::Step(i) => steps
                        .get(*i)
                        .and_then(|step| step.product.as_ref())
                        .ok_or_else(|| {
                            rule.set(RuleId::SourceDependency);
                            Error::Format(format!(
                                "Step {index} has no reconstructed dependency {i}"
                            ))
                        }),
                }
            };
            let resolve_tracked = |source: &Source| -> Result<TrackedDuplex, Error> {
                let m = resolve(source)?;
                let result = match source {
                    Source::Input(_) => Duplex::from_molecule(m).and_then(|duplex| TrackedDuplex::from_unknown_duplex(&duplex)),
                    Source::Step(i) => physical.get(*i).and_then(Option::as_ref).cloned().ok_or_else(|| { rule.set(RuleId::SourceChemistry); Error::Unsupported(format!("Step {i} has no established physical strand/modification state; reconstruct a supported physical graph")) }),
                };
                if result.is_err() {
                    rule.set(RuleId::SourceChemistry);
                }
                result
            };
            let resolve_physical =
                |source: &Source| -> Result<Duplex, Error> { Ok(resolve_tracked(source)?.duplex) };
            rule.set(operation_rule(step));
            let mut next_physical = None;
            let mut diagnostic_reference_only = false;
            let mut product = match step {
                Operation::ResectThreePrime { source, cuts } => {
                    rule.set(RuleId::EndRepair);
                    rule.set(RuleId::EvidenceBudget);
                    if retained_processing + physical_processing_reports.len() >= 64 {
                        return Err(Error::TooLarge);
                    }
                    rule.set(operation_rule(step));
                    let source_model = resolve(source)?;
                    let tracked = resolve_tracked(source)?;
                    rule.set(RuleId::EndRepair);
                    let (processed, geometry) = tracked.resect_three_prime(cuts)?;
                    geometry.validate_documentary()?;
                    let mut product = processed
                        .duplex
                        .molecule(format!("3′-resected operand {}", index + 1))?;
                    let top_end = product.sequence.len();
                    if source_model.sequence.as_str().get(..top_end)
                        != Some(product.sequence.as_str())
                    {
                        return Err(Error::Format(
                            "Resected top reference does not match the actual original source domain"
                                .into(),
                        ));
                    }
                    rule.set(RuleId::AnnotationProjection);
                    project_into(
                        source_model,
                        &[MappingSegment {
                            source_start: 0,
                            source_end: top_end,
                            destination_start: 0,
                            reverse: false,
                        }],
                        &mut product,
                        findings,
                    )?;
                    carry_unmapped_evidence(source_model, &mut product);
                    findings.extend(geometry.warnings.iter().take(100).cloned());
                    findings.push(format!("Step {} removes bottom3′ suffix {} nt and top3′ suffix {} nt. Original5′ end groups and every surviving strand/modification remain; any bottom5′ extension beyond the displayed top reference stays in the actual duplex. No filling, nick sealing or host repair is inferred.",index+1,cuts.left,cuts.right));
                    rule.set(RuleId::EvidenceBudget);
                    step_evidence_bytes = step_evidence_bytes
                        .checked_add(evidence_bytes(&geometry, 16 * 1024 * 1024)?)
                        .filter(|n| *n <= 32 * 1024 * 1024)
                        .ok_or(Error::TooLarge)?;
                    rule.set(operation_rule(step));
                    physical_processing_reports.push(PhysicalProcessingStepReport {
                        step_index: index,
                        evidence: PhysicalProcessingEvidence::ThreePrimeResection(geometry),
                    });
                    next_physical = Some(processed);
                    product
                }
                Operation::HomologyAnnealing {
                    fragments,
                    options,
                    policy,
                } => {
                    rule.set(RuleId::HomologyReconstruction);
                    rule.set(RuleId::EvidenceBudget);
                    if retained_processing + physical_processing_reports.len() >= 64 {
                        return Err(Error::TooLarge);
                    }
                    rule.set(operation_rule(step));
                    policy.validate_shape(fragments.len())?;
                    policy.validate_profile(options.method)?;
                    options.validate()?;
                    if fragments.len() == 1 && policy.topology != Topology::Circular {
                        return Err(Error::Unsupported("Single linear3′ processing uses ResectThreePrime, not a joining assembly operation".into()));
                    }
                    let mut tracked = Vec::with_capacity(fragments.len());
                    let mut oriented = Vec::with_capacity(fragments.len());
                    for fragment in fragments {
                        let material = resolve_tracked(&fragment.source)?;
                        let bottom_end = material
                            .duplex
                            .bottom_start
                            .checked_add(material.duplex.bottom.len() as i64)
                            .ok_or(Error::TooLarge)?;
                        tracked.push(if fragment.reverse {
                            material.reverse()?
                        } else {
                            material
                        });
                        // Keep annotations in their real original top frame. Reversing
                        // unequal strands does not create a reference stock of bottom DNA.
                        oriented.push((
                            resolve(&fragment.source)?.clone(),
                            fragment.reverse,
                            bottom_end,
                        ));
                    }
                    let (_separate_processed_operands, geometry) =
                        modified_duplex::homology_anneal(
                            &tracked.iter().collect::<Vec<_>>(),
                            policy,
                        )?;
                    geometry.validate_documentary()?;
                    // The reference follows actual 5′ termini, including real bottom-only
                    // calls. A top-document concatenation would discard staggered ends.
                    let mut product = Molecule::dna(
                        format!("Annealed reference {}", index + 1),
                        geometry.reference.as_str(),
                        policy.topology,
                    )?;
                    for placement in geometry.source_placements()? {
                        let (source, reverse, bottom_end) =
                            oriented.get(placement.source).ok_or_else(|| {
                                Error::Format("missing annealed source placement".into())
                            })?;
                        let (source_start, source_end, frame_start) = if *reverse {
                            let low =
                                (*bottom_end - placement.contour_length as i64).max(0) as usize;
                            let high = (*bottom_end).max(0) as usize;
                            let high = high.min(source.sequence.len());
                            (
                                low.min(high),
                                high,
                                usize::try_from(*bottom_end - high as i64)
                                    .map_err(|_| Error::TooLarge)?,
                            )
                        } else {
                            (0, source.sequence.len().min(placement.contour_length), 0)
                        };
                        let mut map = Vec::new();
                        let mut cursor = 0;
                        let span_length = source_end - source_start;
                        while cursor < span_length {
                            let absolute = placement
                                .reference_start
                                .checked_add(frame_start)
                                .and_then(|n| n.checked_add(cursor))
                                .ok_or(Error::TooLarge)?;
                            let destination = if policy.topology == Topology::Circular {
                                absolute % product.sequence.len()
                            } else {
                                absolute
                            };
                            let span = (span_length - cursor).min(
                                product.sequence.len().checked_sub(destination).ok_or_else(
                                    || {
                                        Error::Format(
                                            "annealing placement exceeds linear contour".into(),
                                        )
                                    },
                                )?,
                            );
                            if span == 0 || map.len() >= 3 {
                                return Err(Error::TooLarge);
                            }
                            let (start, end) = if *reverse {
                                (source_end - cursor - span, source_end - cursor)
                            } else {
                                (source_start + cursor, source_start + cursor + span)
                            };
                            map.push(MappingSegment {
                                source_start: start,
                                source_end: end,
                                destination_start: destination,
                                reverse: *reverse,
                            });
                            cursor += span;
                        }
                        rule.set(RuleId::AnnotationProjection);
                        project_into_with_collapse(
                            source,
                            &map,
                            &mut product,
                            findings,
                            policy.topology == Topology::Circular && oriented.len() == 1,
                        )?;
                        carry_unmapped_evidence(source, &mut product);
                    }
                    product.name = format!("Annealed reference {}", index + 1);
                    product.metadata.insert(
                        "homology_annealing:policy".into(),
                        serde_json::to_string(policy).map_err(|e| Error::Format(e.to_string()))?,
                    );
                    product.metadata.insert(
                        "homology_annealing:profile".into(),
                        serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?,
                    );
                    findings.extend(
                        geometry
                            .warnings
                            .iter()
                            .chain(&geometry.limits)
                            .take(100)
                            .cloned(),
                    );
                    incomplete.push(format!("Step {} retains {} separate source strands and {} unsealed joining discontinuities. This is a diagnostic sequence contour only: no synthesis, ligation or cellular rescue establishes a reusable assembled product.",index+1,geometry.segments.len(),2*geometry.junctions.len()));
                    if !geometry.terminal_search_complete {
                        incomplete.push(format!(
                            "Step {} bounded exposed-terminal alternative search is incomplete.",
                            index + 1
                        ));
                    }
                    rule.set(RuleId::EvidenceBudget);
                    step_evidence_bytes = step_evidence_bytes
                        .checked_add(evidence_bytes(&geometry, 16 * 1024 * 1024)?)
                        .filter(|n| *n <= 32 * 1024 * 1024)
                        .ok_or(Error::TooLarge)?;
                    rule.set(operation_rule(step));
                    physical_processing_reports.push(PhysicalProcessingStepReport {
                        step_index: index,
                        evidence: PhysicalProcessingEvidence::HomologyAnnealing(geometry),
                    });
                    diagnostic_reference_only = true;
                    product
                }
                Operation::HomologyAssembly {
                    fragments,
                    options,
                    policy,
                } => {
                    rule.set(RuleId::HomologyReconstruction);
                    rule.set(RuleId::EvidenceBudget);
                    if retained_processing + physical_processing_reports.len() >= 64 {
                        return Err(Error::TooLarge);
                    }
                    rule.set(operation_rule(step));
                    policy.validate_shape(fragments.len())?;

                    options.validate()?;
                    // Resolve actual native states first. Reference-only steps and
                    // imported producer metadata cannot supply this operation.
                    let mut tracked = Vec::with_capacity(fragments.len());
                    for fragment in fragments {
                        let material = resolve_tracked(&fragment.source)?;
                        tracked.push(if fragment.reverse {
                            material.reverse()?
                        } else {
                            material
                        });
                    }
                    let refs = tracked.iter().collect::<Vec<_>>();
                    rule.set(RuleId::HomologyReconstruction);
                    let (material, geometry) = modified_duplex::homology_assemble(&refs, policy)?;
                    rule.set(RuleId::EvidenceBudget);
                    step_evidence_bytes = step_evidence_bytes
                        .checked_add(evidence_bytes(&geometry, 32 * 1024 * 1024)?)
                        .filter(|n| *n <= 32 * 1024 * 1024)
                        .ok_or(Error::TooLarge)?;
                    rule.set(operation_rule(step));
                    findings.extend(geometry.warnings.iter().take(100).cloned());
                    findings.extend(geometry.limits.iter().take(100).cloned());
                    if geometry.incomplete {
                        incomplete.push(format!("Step {} physical homology assembly retains unresolved nick chemistry, terminal alternatives or bounded-search limits.",index+1));
                    }
                    let reference_top_start = geometry.reference_top_start;
                    let reference_top_end = geometry.reference_top_end;
                    let reference_length = geometry.reference_length;
                    let diagnostic = if material.is_none() {
                        geometry
                            .product
                            .as_ref()
                            .map(|d| d.molecule(format!("Homology diagnostic {}", index + 1)))
                            .transpose()?
                    } else {
                        None
                    };
                    physical_processing_reports.push(PhysicalProcessingStepReport {
                        step_index: index,
                        evidence: PhysicalProcessingEvidence::HomologyAssembly(geometry),
                    });
                    let mut product = if let Some(material) = &material {
                        material
                            .duplex
                            .molecule(format!("Homology assembly {}", index + 1))?
                    } else {
                        diagnostic_reference_only = true;
                        incomplete.push(format!("Step {} retains paired diagnostic strands but withholds reusable material: internal nick chemistry cannot be faithfully represented by the current private duplex model.",index+1));
                        diagnostic.ok_or_else(|| Error::Unsupported("Physical homology assembly withheld its product: no unsupported gap fill or reference-derived material is invented.".into()))?
                    };
                    // Advisory sequence guidance is checked on the exact source
                    // molecules, after actual strand reconstruction. It cannot
                    // choose a physical path or create the private material above.
                    let sources = fragments
                        .iter()
                        .map(|f| resolve(&f.source).cloned())
                        .collect::<Result<Vec<_>, _>>()?;
                    let reverse = fragments.iter().map(|f| f.reverse).collect::<Vec<_>>();
                    let circular = policy.topology == Topology::Circular;
                    let profile = crate::cloning::overlap::profiles::check_oriented_fragments(
                        &sources,
                        &reverse,
                        &policy.overlaps,
                        circular,
                        options,
                    )?;
                    if !profile.geometry_valid {
                        return Err(Error::Format("Source homology guidance disagrees with freshly reconstructed material geometry".into()));
                    }
                    if !profile.complete {
                        incomplete.push(format!("Step {} advisory homology repeat scan reached its bound; alternatives remain unsearched.",index+1));
                    }
                    findings.extend(profile.findings.iter().take(100).cloned());
                    findings.extend(profile.warnings.iter().take(100).cloned());
                    let mut oriented = Vec::with_capacity(sources.len());
                    for (source, reverse) in sources.iter().zip(&reverse) {
                        if *reverse {
                            let mut model = Molecule::dna(
                                format!("Oriented {}", source.name),
                                source.sequence.reverse_complement()?.as_str(),
                                source.topology,
                            )?;
                            rule.set(RuleId::AnnotationProjection);
                            project_into(
                                source,
                                &[MappingSegment {
                                    source_start: 0,
                                    source_end: source.sequence.len(),
                                    destination_start: 0,
                                    reverse: true,
                                }],
                                &mut model,
                                findings,
                            )?;
                            carry_unmapped_evidence(source, &mut model);
                            oriented.push(model);
                        } else {
                            oriented.push(source.clone());
                        }
                    }
                    let refs = oriented.iter().collect::<Vec<_>>();
                    let reference =
                        assemble_overlap(&refs, &policy.overlaps, circular, index, findings)?;
                    if reference.sequence.len() != reference_length
                        || reference_top_start >= reference_top_end
                        || reference
                            .sequence
                            .as_str()
                            .get(reference_top_start..reference_top_end)
                            != Some(product.sequence.as_str())
                    {
                        return Err(Error::Format("Source-only annotation reference differs from independently constructed physical strands".into()));
                    }
                    rule.set(RuleId::AnnotationProjection);
                    project_into(
                        &reference,
                        &[MappingSegment {
                            source_start: reference_top_start,
                            source_end: reference_top_end,
                            destination_start: 0,
                            reverse: false,
                        }],
                        &mut product,
                        findings,
                    )?;
                    carry_unmapped_evidence(&reference, &mut product);
                    product.metadata.insert(
                        "homology_assembly:policy".into(),
                        serde_json::to_string(policy).map_err(|e| Error::Format(e.to_string()))?,
                    );
                    product.metadata.insert(
                        "homology_assembly:profile".into(),
                        serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?,
                    );
                    findings.push(format!("Step {} reconstructs retained source strands, declared 5′ cuts, template-supported synthesis and eligible nick seals before target comparison. This explicit ideal model is not experimental success, host rescue or proprietary kit chemistry.",index+1));
                    next_physical = material;
                    product
                }
                Operation::OverlapExtension {
                    left,
                    right,
                    options,
                } => {
                    rule.set(RuleId::OverlapGeometry);
                    rule.set(RuleId::EvidenceBudget);
                    if retained_processing + physical_processing_reports.len() >= 64 {
                        return Err(Error::TooLarge);
                    }
                    rule.set(operation_rule(step));
                    let mut a = resolve_tracked(&left.source)?;
                    let mut b = resolve_tracked(&right.source)?;
                    if left.reverse {
                        a = a.reverse()?;
                    }
                    if right.reverse {
                        b = b.reverse()?;
                    }
                    rule.set(RuleId::OverlapGeometry);
                    let (extended, geometry) = a.overlap_extend(&b, options)?;
                    let mut product = extended
                        .duplex
                        .molecule(format!("Overlap extension {}", index + 1))?;
                    let offset = a.duplex.top.len() - options.overlap;
                    for (oriented, destination_start) in [(left, 0), (right, offset)] {
                        let source = resolve(&oriented.source)?;
                        rule.set(RuleId::AnnotationProjection);
                        project_into(
                            source,
                            &[MappingSegment {
                                source_start: 0,
                                source_end: source.sequence.len(),
                                destination_start,
                                reverse: oriented.reverse,
                            }],
                            &mut product,
                            findings,
                        )?;
                    }
                    findings.extend(geometry.warnings.iter().take(100).cloned());
                    findings.push(format!("Step {} denatures two intact duplexes, pairs the left top strand with the right bottom strand over {} bp and extends their paired 3′ ends. Retained strands preserve modifications; synthesized complements are assumed unmodified. Selection does not prove isolation, amplification efficiency or reaction-mixture uniqueness.", index+1, options.overlap));
                    if !geometry.terminal_pairings.complete || !geometry.terminal_pairings.unique {
                        incomplete.push(format!("Step {} selected terminal pairing is not uniquely established within the bounded terminal-pair model; alternative products remain possible.", index+1));
                    }
                    rule.set(RuleId::EvidenceBudget);
                    step_evidence_bytes = step_evidence_bytes
                        .checked_add(evidence_bytes(&geometry, 32 * 1024 * 1024)?)
                        .filter(|n| *n <= 32 * 1024 * 1024)
                        .ok_or(Error::TooLarge)?;
                    rule.set(operation_rule(step));
                    physical_processing_reports.push(PhysicalProcessingStepReport {
                        step_index: index,
                        evidence: PhysicalProcessingEvidence::OverlapExtension(geometry),
                    });
                    next_physical = Some(extended);
                    product
                }
                Operation::ReconstructMaterial { template, origin } => {
                    rule.set(RuleId::MaterialOrigin);
                    if !matches!(template, Source::Input(_)) {
                        return Err(Error::Format(
                            "Derived material reconstruction requires an exact original input event"
                                .into(),
                        ));
                    }
                    let fresh = crate::material::derivation::replay(origin, resolve(template)?)?;
                    crate::material::derivation::reaction_conditions::validate(
                        fresh.reaction_conditions(),
                    )?;
                    let reaction_conditions = fresh.reaction_conditions().to_vec();
                    let condition_bytes = evidence_bytes(&reaction_conditions, 4 * 1024 * 1024)?;
                    let modification_summary = fresh.tracked().summary()?;
                    let duplex_summary = fresh.tracked().duplex.summary()?;
                    let bytes =
                        evidence_bytes(origin, crate::material::derivation::MAX_CAPSULE_BYTES)?
                            .checked_add(evidence_bytes(&modification_summary, 32 * 1024 * 1024)?)
                            .and_then(|bytes| bytes.checked_add(condition_bytes))
                            .ok_or(Error::TooLarge)?;
                    rule.set(RuleId::EvidenceBudget);
                    step_evidence_bytes = step_evidence_bytes
                        .checked_add(bytes)
                        .filter(|n| *n <= 32 * 1024 * 1024)
                        .ok_or(Error::TooLarge)?;
                    rule.set(operation_rule(step));
                    let (
                        reference,
                        tracked,
                        origin_incomplete,
                        origin_findings,
                        origin_findings_total,
                        origin_findings_truncated,
                    ) = fresh.into_parts();
                    let evidence = DerivedOriginEvidence {
                        reaction_conditions,
                        origin: origin.clone(),
                        material: MaterialStateEvidence {
                            material_revision: modification_summary.state_revision.clone(),
                            duplex: duplex_summary,
                            modifications: modification_summary,
                            origin_incomplete,
                            origin_findings,
                            origin_findings_total,
                            origin_findings_truncated,
                        },
                    };
                    evidence.validate_documentary()?;
                    findings.push(format!("Step {} freshly reconstructs saved material derivation event {} from its exact original source and typed recipe. Actual strand termini and positional modifications are retained; this is material retrieval, not a new target check or proof of experimental isolation.",index+1,origin.product.event_id));
                    if origin_incomplete {
                        incomplete.push(format!("Step {} retains incomplete origin findings; saved material reuse does not upgrade them.",index+1));
                    }
                    physical_processing_reports.push(PhysicalProcessingStepReport {
                        step_index: index,
                        evidence: PhysicalProcessingEvidence::DerivedOrigin(evidence),
                    });
                    next_physical = Some(tracked);
                    reference
                }
                Operation::ReconstructRunProduct { template, origin } => {
                    rule.set(RuleId::SavedMaterial);
                    if !matches!(template, Source::Input(_)) {
                        return Err(Error::Format("Saved material reconstruction requires an explicit saved input event, not an enclosing step projection".into()));
                    }
                    let fresh = crate::material::capsule::replay(origin, resolve(template)?)?;
                    let bytes = evidence_bytes(
                        &fresh.evidence,
                        crate::material::capsule::MAX_REPORT_BYTES,
                    )?;
                    rule.set(RuleId::EvidenceBudget);
                    step_evidence_bytes = step_evidence_bytes
                        .checked_add(bytes)
                        .filter(|n| *n <= 32 * 1024 * 1024)
                        .ok_or(Error::TooLarge)?;
                    rule.set(operation_rule(step));
                    findings.push(format!("Step {} freshly reconstructed saved run {} product {} against its original independent target; strands, nicks, termini, positional modifications and generation provenance remain actual tracked state. The selected event is a declared source, not proof of isolation or experimental success.", index+1, origin.selection.run_id, origin.selection.step_index+1));
                    if fresh.incomplete {
                        incomplete.push(format!("Step {} originating run retains incomplete findings; material reuse does not upgrade that assessment. Complete originating evidence is retained.", index+1));
                    }
                    physical_processing_reports.push(PhysicalProcessingStepReport {
                        step_index: index,
                        evidence: PhysicalProcessingEvidence::RunProduct(fresh.evidence),
                    });
                    next_physical = Some(fresh.tracked);
                    fresh.reference
                }
                Operation::MultiSiteRecombineReference {
                    sources,
                    request,
                    selection,
                } => {
                    rule.set(RuleId::RecombinationReference);
                    let resolved = sources
                        .iter()
                        .map(|source| resolve(&source.source))
                        .collect::<Result<Vec<_>, _>>()?;
                    crate::cloning::recombination::multisite::source_preflight(&resolved, request)?;
                    let inputs = sources
                        .iter()
                        .zip(&resolved)
                        .map(
                            |(source, molecule)| crate::cloning::recombination::reference::Input {
                                role: source.role.clone(),
                                molecule: (*molecule).clone(),
                            },
                        )
                        .collect::<Vec<_>>();
                    let fresh =
                        crate::cloning::recombination::multisite::simulate(&inputs, request)?;
                    let bytes = evidence_bytes(&fresh, 32 * 1024 * 1024)?;
                    rule.set(RuleId::EvidenceBudget);
                    step_evidence_bytes = step_evidence_bytes
                        .checked_add(bytes)
                        .filter(|n| *n <= 32 * 1024 * 1024)
                        .ok_or(Error::TooLarge)?;
                    rule.set(operation_rule(step));
                    let selected = fresh.selected(selection);
                    findings.push(format!("Step {} retains a whole {}-parent multi-site reference exchange: {} candidate networks and every reciprocal component. No target selects pairings; reference provenance is not physical strand/modification authority.",index+1,inputs.len(),fresh.reference.candidates.len()));
                    reference_step_reports.push(ReferenceStepReport {
                        step_index: index,
                        evidence: ReferenceStepEvidence::MultiSiteRecombination {
                            selection: selection.clone(),
                            batch: fresh,
                        },
                    });
                    let selected = selected?;
                    incomplete.push(format!("Step {} is a declared multi-site reference design. Recombination chemistry, source-strand inheritance, reaction efficacy and mixture selection remain unestablished.",index+1));
                    selected.molecule
                }
                Operation::TopoActions { sources, request } => {
                    rule.set(RuleId::TopoReconstruction);
                    let molecules = sources
                        .iter()
                        .map(|s| resolve(&s.source))
                        .collect::<Result<Vec<_>, _>>()?;
                    validate_topo_shape(sources, &request.preparation)?;
                    collagene_model::bounds::model_preflight(&molecules)?;
                    molecules.iter().try_fold(0usize, |sum, molecule| {
                        sum.checked_add(molecule.sequence.len())
                            .filter(|n| *n <= 100_000)
                            .ok_or(Error::TooLarge)
                    })?;
                    let mut templates = Vec::new();
                    for (source, molecule) in sources.iter().zip(&molecules) {
                        match &source.source {
                            Source::Step(_) => {
                                // Only actual fresh private material can enter the
                                // PCR-only synthesis route. No imported summary or
                                // matching sequence establishes this permission.
                                templates
                                    .push((source.role.clone(), resolve_tracked(&source.source)?));
                            }
                            Source::Input(_) => {
                                if chemistry::has_unmapped_physical_evidence(molecule)
                                    || molecule.metadata.contains_key("tracked_modifications")
                                    || molecule.metadata.contains_key(
                                        crate::cloning::recombination::REFERENCE_ONLY_KEY,
                                    )
                                    || molecule.metadata.keys().any(|key| key.starts_with("topo:"))
                                {
                                    return Err(Error::Unsupported("TopoActions requires original reference sources with explicit preparation declarations; imported physical/reference-only summaries cannot restore their lost graph or tracked modification state.".into()));
                                }
                            }
                        }
                    }
                    // Native input cloning remains bounded before preparation. The
                    // private adapter checks each Step role is PCR-only and agrees
                    // exactly with its tracked template; source metadata is retained.
                    let inputs = sources
                        .iter()
                        .zip(&molecules)
                        .map(|(source, molecule)| crate::cloning::topo::Input {
                            role: source.role.clone(),
                            molecule: (*molecule).clone(),
                        })
                        .collect::<Vec<_>>();
                    let tracked_templates = templates
                        .iter()
                        .map(|(role, material)| (role.clone(), material))
                        .collect::<Vec<_>>();
                    let fresh = crate::cloning::topo::actions::inspect_reconstructed(
                        &inputs,
                        request,
                        &tracked_templates,
                    )?;
                    retain_topo_source_incompleteness(&fresh.preparation, index, &mut incomplete);
                    if !templates.is_empty() {
                        findings.push(format!("Step {} uses {} fresh upstream material role(s) only as templates for explicitly declared new PCR synthesis. Template atoms, internal modifications and generation history are not inherited into the new strands; arbitrary Prepared correspondence from prior products remains unsupported.", index + 1, templates.len()));
                    }
                    let fresh_bytes =
                        evidence_bytes(&fresh, crate::cloning::topo::actions::MAX_REPORT_BYTES)?;
                    rule.set(RuleId::EvidenceBudget);
                    step_evidence_bytes = step_evidence_bytes
                        .checked_add(fresh_bytes)
                        .filter(|n| *n <= 32 * 1024 * 1024)
                        .ok_or(Error::TooLarge)?;
                    rule.set(operation_rule(step));
                    let selected_bases = fresh
                        .product
                        .as_ref()
                        .map_or(0, |m| m.sequence.len())
                        .checked_add(
                            fresh
                                .tracked
                                .as_ref()
                                .map_or(0, |t| t.top.len() + t.bottom.len()),
                        )
                        .ok_or(Error::TooLarge)?;
                    let tracked_bytes = fresh
                        .tracked
                        .as_ref()
                        .map(|t| evidence_bytes(t, 32 * 1024 * 1024))
                        .transpose()?
                        .unwrap_or(0);
                    let selected_document_bytes = fresh
                        .product
                        .as_ref()
                        .map(|m| evidence_bytes(m, 32 * 1024 * 1024))
                        .transpose()?
                        .unwrap_or(0);
                    if product_bases
                        .checked_add(selected_bases)
                        .is_none_or(|n| n > 32_000_000)
                        || product_bytes
                            .checked_add(step_evidence_bytes)
                            .and_then(|n| n.checked_add(tracked_bytes))
                            .and_then(|n| n.checked_add(selected_document_bytes))
                            .is_none_or(|n| n > 64 * 1024 * 1024)
                    {
                        return Err(Error::Unsupported("TopoActions selected material plus full retained graph exceeds the plan material/evidence budget before duplicating its selected product.".into()));
                    }
                    let graph = &fresh.projection.physical.product;
                    findings.push(format!("Step {} fresh TopoActions retains {} components, {} original source calls, {} remaining graph nicks and {} unpaired calls. Explicit component selection is not physical isolation or a mixture winner; all other components and released activations remain typed evidence.",index+1,graph.components.len(),fresh.projection.physical.source_bases,graph.nicks.len(),graph.unpaired.iter().map(|r|r.end-r.start).sum::<usize>()));
                    if !fresh.projection.physical.all_requested_steps_applied {
                        findings.push(format!("Step {} retains an unapplied/unknown requested action under its explicit projection options; target agreement cannot establish that it occurred.",index+1));
                    }
                    findings.extend(fresh.bridge_notes.iter().cloned());
                    let actual = fresh.product.as_ref().cloned();
                    let tracked = fresh.tracked.as_ref().cloned();
                    physical_processing_reports.push(PhysicalProcessingStepReport {
                        step_index: index,
                        evidence: PhysicalProcessingEvidence::TopoActions(fresh),
                    });
                    let (Some(actual), Some(tracked)) = (actual, tracked) else {
                        return Err(Error::Unsupported(format!(
                            "Step {} TopoActions material is withheld: no faithfully representable paired duplex. Full preparation, parent/sealed graphs, flaps, components and terminal evidence are retained; no target-derived reference product was substituted.",
                            index + 1
                        )));
                    };
                    rule.set(RuleId::ProductBudget);
                    product_bytes = product_bytes
                        .checked_add(tracked_bytes)
                        .ok_or(Error::TooLarge)?;
                    next_physical = Some(tracked);
                    actual
                }

                Operation::Topo {
                    sources,
                    request,
                    selection,
                }
                | Operation::CheckedTopo {
                    sources,
                    request,
                    selection,
                    ..
                } => {
                    rule.set(RuleId::TopoReconstruction);
                    let molecules = sources
                        .iter()
                        .map(|source| resolve(&source.source))
                        .collect::<Result<Vec<_>, _>>()?;
                    validate_topo_shape(sources, request)?;
                    collagene_model::bounds::model_preflight(&molecules)?;
                    let inputs = sources
                        .iter()
                        .zip(&molecules)
                        .map(|(s, m)| crate::cloning::topo::Input {
                            role: s.role.clone(),
                            molecule: (*m).clone(),
                        })
                        .collect::<Vec<_>>();
                    let mut tracked_templates = vec![];
                    for source in sources {
                        if matches!(source.source, Source::Step(_)) {
                            tracked_templates
                                .push((source.role.clone(), resolve_tracked(&source.source)?));
                        }
                    }
                    let templates = tracked_templates
                        .iter()
                        .map(|(r, t)| (r.clone(), t))
                        .collect::<Vec<_>>();
                    let task = if let Operation::CheckedTopo { policy, .. } = step {
                        crate::cloning::topo::run_preparation::Task::Check {
                            request: request.clone(),
                            selection: selection.clone(),
                            target: target.clone(),
                            policy: policy.clone(),
                        }
                    } else {
                        crate::cloning::topo::run_preparation::Task::Project {
                            request: request.clone(),
                            selection: selection.clone(),
                        }
                    };
                    let reply = crate::cloning::topo::inspect_reconstructed_task(
                        &inputs, &templates, &task,
                    )
                    .map_err(|error| {
                        Error::Format(format!(
                            "TOPO source/recipe reconstruction rejected: {error}"
                        ))
                    })?;
                    let (fresh, outcome) = match reply {
                        crate::cloning::topo::run_preparation::Reply::Projected(p) => (p, None),
                        crate::cloning::topo::run_preparation::Reply::Checked(a) => {
                            (a.product, Some(a.outcome))
                        }
                        _ => return Err(Error::Format("native TOPO reply kind changed".into())),
                    };
                    retain_topo_source_incompleteness(&fresh.preparation, index, &mut incomplete);
                    let fully_specified = fresh.projection.physical.fully_specified_dna;
                    let physical_bytes =
                        evidence_bytes(&fresh, crate::cloning::topo::MAX_REPORT_BYTES)?;
                    rule.set(RuleId::EvidenceBudget);
                    step_evidence_bytes = step_evidence_bytes
                        .checked_add(physical_bytes)
                        .filter(|n| *n <= 32 * 1024 * 1024)
                        .ok_or(Error::TooLarge)?;
                    rule.set(operation_rule(step));
                    let preview = fresh.reference.clone();
                    step_reports.push(StepReport {
                        step_index: index,
                        reaction_conditions: None,
                        topo: Some(fresh),
                    });
                    findings.push(format!("Step {} retains its complete freshly reconstructed TOPO preparation, transfer graph, all physical components and paired-contour proof. A selected reference view does not seal nicks, remove flaps or establish experimental success.",index+1));
                    rule.set(RuleId::TopoTargetPolicy);
                    match &outcome {
                        Some(crate::cloning::topo::TargetOutcome::Rejected { reasons }) => return Err(Error::Format(format!("CheckedTopo independent-target policy rejected: {}", reasons.join("; ")))),
                        Some(crate::cloning::topo::TargetOutcome::Withheld { reason }) => return Err(Error::Unsupported(format!("CheckedTopo target assessment withheld: {reason}"))),
                        Some(crate::cloning::topo::TargetOutcome::AcceptedUnderDeclaredGraph) if !fully_specified => return Err(Error::Unsupported("CheckedTopo graph retains ambiguous DNA; target agreement cannot establish fully determined sequence.".into())),
                        _ => {}
                    }
                    rule.set(RuleId::TopoSelection);
                    let Some(preview) = preview else {
                        return Err(Error::Unsupported(format!(
                            "Step {} has no uniquely proven paired TOPO contour. Its full physical graph is retained; no product DNA or target-derived path is invented.",
                            index + 1
                        )));
                    };
                    if outcome.is_none() {
                        incomplete.push(format!("Step {} TOPO contour matches only an explicitly selected reference path. No physical target-acceptance policy was supplied, and downstream graph consumption is not modeled; nicks, other components and transfer outcomes remain in full step evidence.",index+1));
                    }
                    preview
                }
                Operation::RecombineReference {
                    sources,
                    request,
                    selection,
                } => {
                    rule.set(RuleId::RecombinationReference);
                    let molecules = [resolve(&sources[0])?, resolve(&sources[1])?];
                    // Borrowed preflight happens before cloning the immutable source snapshots.
                    let inputs = recombination_inputs(sources, molecules, request)?;
                    let fresh = crate::cloning::recombination::simulate(&inputs, request)?;
                    rule.set(RuleId::RecombinationReference);
                    if fresh.batch_revision != selection.batch_revision {
                        return Err(Error::Format("Reference component selection belongs to a different source/profile/options request".into()));
                    }
                    findings.push(format!("Step {} reference search: {} candidate(s), {} site(s), complete under declared profile: {}. Profile snapshot {} retains all reciprocal components; selected candidate {} component {} only.", index+1, fresh.reference.candidates.len(), fresh.reference.sites.len(), fresh.reference.complete, request.profile.revision, selection.candidate_id, selection.component_id));
                    findings.extend(fresh.reference.diagnostics);
                    findings.extend(fresh.warnings);
                    let selected = fresh
                        .documents
                        .into_iter()
                        .find(|p| &p.selection == selection)
                        .ok_or_else(|| {
                            Error::Format(
                                "Selected reference candidate/component was not reconstructed"
                                    .into(),
                            )
                        })?;
                    findings.extend(selected.warnings.iter().cloned());
                    incomplete.push(format!("Step {} selects candidate {} component {} under pinned reference profile {} / {} ({}). All reciprocal components were reconstructed; this establishes reference sequence only, not recombinase efficacy, physical strand/modification states, cellular repair, selection or chemical uniqueness.", index + 1, selection.candidate_id, selection.component_id, request.profile.profile.id, request.profile.profile.edition, request.profile.revision));
                    selected.molecule
                }
                Operation::Pcr {
                    template,
                    forward,
                    reverse,
                    template_strand,
                } => {
                    rule.set(RuleId::PcrBinding);
                    let substrate = resolve(template)?;
                    rule.set(RuleId::PcrBinding);
                    pcr_budget.charge(substrate.sequence.len())?;
                    pcr_budget.charge(substrate.features.len())?;
                    rule.set(RuleId::SourceChemistry);
                    let frame = match template {
                        Source::Input(_) => {
                            if *template_strand != chemistry::Strand::Top {
                                return Err(Error::Unsupported(
                                    "Bottom PCR requires a freshly reconstructed native source step".into(),
                                ));
                            }
                            crate::cloning::pcr::source::Template::reference(substrate)?
                        }
                        Source::Step(producer) => {
                            let tracked = physical.get(*producer)
                                .and_then(Option::as_ref)
                                .ok_or_else(|| Error::Unsupported(
                                    "PCR requires actual physical source state; reference-only evidence is not a template".into(),
                                ))?;
                            rule.set(RuleId::PcrBinding);
                            pcr_budget.charge(tracked.top.len())?;
                            pcr_budget.charge(tracked.bottom.len())?;
                            rule.set(RuleId::SourceChemistry);
                            crate::cloning::pcr::source::Template::fresh(
                                substrate,
                                tracked,
                                *template_strand,
                            )?
                        }
                    };
                    rule.set(RuleId::PcrBinding);
                    let (matches, _) =
                        frame.geometries_with_budget(forward, reverse, pcr_budget)?;
                    if matches.is_empty() {
                        return Err(Error::Format(format!("step {}: no PCR product", index + 1)));
                    }
                    rule.set(RuleId::PcrSelection);
                    if matches.len() != 1 {
                        return Err(Error::Unsupported(format!(
                            "step {}: {} PCR products require disambiguation",
                            index + 1,
                            matches.len()
                        )));
                    }
                    let geometry = &matches[0];
                    rule.set(RuleId::ProductIntegrity);
                    let duplex = Duplex::exact_pcr(&geometry.product.sequence)?;
                    let mut product = duplex.molecule(format!("PCR {}", index + 1))?;
                    rule.set(RuleId::AnnotationProjection);
                    project_pcr_frame_annotations(frame, geometry, &mut product, findings)?;
                    next_physical = Some(TrackedDuplex::exact_pcr(&duplex.top)?);
                    product
                }
                Operation::PcrWithMismatches {
                    template,
                    forward,
                    reverse,
                    options,
                } => {
                    rule.set(RuleId::MismatchPcrPolicy);
                    let substrate = resolve(template)?;
                    rule.set(RuleId::SourceChemistry);
                    let physical_template = match template {
                        Source::Input(_) => {
                            if chemistry::has_unmapped_physical_evidence(substrate) {
                                return Err(Error::Unsupported("Mismatched PCR of imported physical evidence requires its originating graph".into()));
                            }
                            None
                        }
                        Source::Step(i) => {
                            let duplex = physical.get(*i).and_then(Option::as_ref).ok_or_else(|| {
                                Error::Unsupported("Mismatched PCR requires an established physical template from the originating step; a reference-only assembly cannot be materialized by PCR".into())
                            })?;
                            if duplex.bottom_start != 0
                                || duplex.bottom.len() != duplex.top.len()
                                || !duplex.nicks.is_empty()
                            {
                                return Err(Error::Unsupported(
                                    "Mismatched PCR of staggered/nicked templates is not modeled"
                                        .into(),
                                ));
                            }
                            Some(duplex)
                        }
                    };
                    rule.set(RuleId::MismatchPcrPolicy);
                    let report = match physical_template {
                        Some(duplex) => {
                            crate::cloning::pcr::mutagenesis::simulate_pcr_reconstructed(
                                substrate, duplex, forward, reverse, options,
                            )?
                        }
                        None => crate::cloning::pcr::mutagenesis::simulate_pcr(
                            substrate, forward, reverse, options,
                        )?,
                    };
                    findings.extend(
                        report
                            .warnings
                            .iter()
                            .take(20)
                            .map(|s| s.chars().take(2048).collect()),
                    );
                    if report.products.len() != 1 {
                        rule.set(
                            if report.enumeration_complete && report.products.is_empty() {
                                RuleId::PcrBinding
                            } else {
                                RuleId::PcrSelection
                            },
                        );
                        return Err(
                            if report.enumeration_complete && report.products.is_empty() {
                                Error::Format(
                                    "No mismatched PCR product under the recorded binding policy"
                                        .into(),
                                )
                            } else {
                                Error::Unsupported(format!(
                                    "{} mismatched PCR geometries require disambiguation; enumeration complete: {}",
                                    report.products.len(),
                                    report.enumeration_complete
                                ))
                            },
                        );
                    }
                    let geometry = &report.products[0];
                    rule.set(RuleId::PcrSelection);
                    if !geometry.certain || geometry.overlapping_domains {
                        return Err(Error::Unsupported("Ambiguous or overlapping-primer products remain previews, not established amplified duplexes".into()));
                    }
                    if !report.enumeration_complete {
                        incomplete.push(
                            "Mismatched PCR alternative-product search or geometry is incomplete"
                                .to_string(),
                        );
                    }
                    rule.set(RuleId::ProductIntegrity);
                    let duplex = Duplex::exact_pcr(&geometry.sequence)?;
                    let mut product = duplex.molecule(format!("Mismatched PCR {}", index + 1))?;
                    rule.set(RuleId::AnnotationProjection);
                    project_mutant_pcr(substrate, geometry, &mut product, findings)?;
                    findings.push(format!("Mismatched PCR introduced {} template-base changes under explicit 3′/mismatch policy; amplification efficiency is not established.",geometry.substitutions.len()));
                    next_physical = Some(TrackedDuplex::exact_pcr(&duplex.top)?);
                    product
                }
                Operation::GoldenGate {
                    templates,
                    options,
                    candidate_index,
                }
                | Operation::SelectedGoldenGate {
                    templates,
                    options,
                    candidate_index,
                    ..
                } => {
                    rule.set(RuleId::GoldenGateAssembly);
                    let sources = templates
                        .iter()
                        .map(|s| resolve(s).cloned())
                        .collect::<Result<Vec<_>, _>>()?;
                    let substrates = templates
                        .iter()
                        .map(&resolve_tracked)
                        .collect::<Result<Vec<_>, _>>()?;
                    let fresh_gate = crate::cloning::golden_gate::simulate_tracked(
                        &sources,
                        &substrates,
                        options,
                        Some(*candidate_index),
                        cleavage_budget,
                    )?;
                    let selected_tracked = fresh_gate.selected;
                    let mut local_report = Some(fresh_gate.report);
                    let report = if index + 1 == plan.steps.len() {
                        *final_golden_gate = local_report.take();
                        final_golden_gate.as_ref().expect("fresh final enumeration")
                    } else {
                        local_report
                            .as_ref()
                            .expect("fresh intermediate enumeration")
                    };
                    rule.set(RuleId::GoldenGateAssembly);
                    let selected=report.candidates.get(*candidate_index).ok_or_else(||if report.enumeration_complete{Error::Format("Golden Gate candidate index does not exist".into())}else{Error::Unsupported("Golden Gate search is incomplete and the selected candidate is unavailable".into())})?;
                    if let Operation::SelectedGoldenGate {
                        directed_sequence_key,
                        physical_sequence_key,
                        ..
                    } = step
                    {
                        if directed_sequence_key.len() != 64
                            || physical_sequence_key.len() != 64
                            || selected.directed_sequence_key != *directed_sequence_key
                            || selected.physical_sequence_key != *physical_sequence_key
                        {
                            return Err(Error::Format("Golden Gate explicit candidate keys disagree with the fresh native enumeration".into()));
                        }
                    }
                    if !selected.cleavage.complete {
                        return Err(Error::Unsupported("Golden Gate candidate site absence is unproven: bounded recut evidence is incomplete; no site-free cycling-product claim is established".into()));
                    }
                    if !selected.restriction_site_free {
                        return Err(Error::Format("Golden Gate candidate retains recognition sites; it is not a site-free cycling product".into()));
                    }
                    if !report.enumeration_complete {
                        incomplete.push(format!(
                            "Golden Gate enumeration incomplete: {}",
                            report.truncation_reasons.join(" ")
                        ));
                    }
                    findings.push(format!("Golden Gate checked candidate {} of {} enumerated paths; {} stable physical sequence groups, {} uncut carryover sources. This finite one-copy model does not establish exhaustive reaction-mixture uniqueness.",candidate_index+1,report.candidates.len(),report.stable_physical_product_count,report.carryover_sources.len()));
                    findings.extend(
                        report
                            .warnings
                            .iter()
                            .take(20)
                            .map(|s| s.chars().take(2048).collect()),
                    );
                    let duplex = selected.duplex.clone();
                    let mut product = duplex.molecule(format!("Golden Gate {}", index + 1))?;
                    let mut offset = 0usize;
                    for oriented in &selected.path {
                        let fragment = &report.fragments[oriented.fragment];
                        let source = &sources[fragment.source_index];
                        let mut intermediate =
                            fragment.duplex.molecule("Golden Gate source fragment")?;
                        rule.set(RuleId::AnnotationProjection);
                        let map = region_map(
                            source.sequence.len(),
                            fragment.top_start,
                            fragment.duplex.top.len(),
                            0,
                        )?;
                        rule.set(RuleId::AnnotationProjection);
                        project_into(source, &map, &mut intermediate, findings)?;
                        let oriented_duplex = if oriented.reverse {
                            fragment.duplex.reverse()?
                        } else {
                            fragment.duplex.clone()
                        };
                        let map = if oriented.reverse {
                            reverse_reference_map(
                                intermediate.sequence.len(),
                                fragment.duplex.bottom_start + fragment.duplex.bottom.len() as i64,
                                offset,
                                product.sequence.len(),
                                product.topology,
                            )?
                        } else {
                            vec![MappingSegment {
                                source_start: 0,
                                source_end: intermediate.sequence.len(),
                                destination_start: offset,
                                reverse: false,
                            }]
                        };
                        rule.set(RuleId::AnnotationProjection);
                        project_into(&intermediate, &map, &mut product, findings)?;
                        offset += oriented_duplex.top.len();
                    }
                    next_physical=Some(selected_tracked.ok_or_else(||Error::Unsupported("Selected Golden Gate material was not freshly reconstructed under checked cleavage".into()))?);
                    product
                }
                Operation::Overlap {
                    fragments,
                    overlaps,
                    circular,
                } => {
                    rule.set(RuleId::OverlapGeometry);
                    if fragments.is_empty()
                        || fragments.len() > 1000
                        || overlaps.len() != fragments.len() - 1 + usize::from(*circular)
                    {
                        return Err(Error::Format(
                            "ordered assembly requires one overlap per junction".into(),
                        ));
                    }
                    for source in fragments {
                        if let Source::Input(_) = source {
                            let input = resolve(source)?;
                            if chemistry::has_unmapped_physical_evidence(input) {
                                return Err(Error::Unsupported("overlap assembly of imported physical-intermediate evidence requires reconstruction of its originating graph".into()));
                            }
                        }
                        if let Source::Step(i) = source {
                            if let Some(Some(duplex)) = physical.get(*i) {
                                let ends = duplex.ends()?;
                                if duplex.left.three_prime == chemistry::ThreePrime::Blocked
                                    || duplex.right.three_prime == chemistry::ThreePrime::Blocked
                                {
                                    return Err(Error::Unsupported("overlap assembly of a known blocked3′ terminus requires an explicitly modeled unblocking/resection mechanism".into()));
                                }
                                if !ends.0.overhang.is_empty()
                                    || !ends.1.overhang.is_empty()
                                    || !duplex.nicks.is_empty()
                                {
                                    return Err(Error::Unsupported("overlap assembly of staggered or nicked physical substrates is not yet modeled".into()));
                                }
                            }
                        }
                    }
                    let fragments = fragments
                        .iter()
                        .map(resolve)
                        .collect::<Result<Vec<_>, _>>()?;
                    assemble_overlap(&fragments, overlaps, *circular, index, findings)?
                }
                Operation::DeclareSourceTermini {
                    template,
                    left,
                    right,
                } => {
                    rule.set(RuleId::SourceChemistry);
                    if !matches!(template, Source::Input(_)) {
                        return Err(Error::Unsupported("Source termini may only be declared on an original input; descendants require their actual physical graph, not a new intact-source hypothesis.".into()));
                    }
                    let substrate = resolve(template)?;
                    declared_source_termini_preflight(substrate)?;
                    // Start from the same guarded intact reference assumption as
                    // other input chemistry; summaries never supply strand state.
                    let mut duplex = Duplex::from_molecule(substrate)?;
                    duplex.left = *left;
                    duplex.right = *right;
                    duplex.validate()?;
                    let tracked = TrackedDuplex::from_unknown_duplex(&duplex)?;
                    // DNA and all source feature/citation evidence are unchanged.
                    let mut product = substrate.clone();
                    product.metadata.insert(
                        "duplex_model".into(),
                        serde_json::to_string(&duplex.summary()?)
                            .map_err(|e| Error::Format(e.to_string()))?,
                    );
                    product.metadata.insert("source:declared_termini".into(),serde_json::to_string(&(substrate.revision(),left,right,"Operator hypothesis for intact original linear DNA; not measured chemistry, synthesis, repair or proof of material state.")).map_err(|e|Error::Format(e.to_string()))?);
                    findings.push(format!("Step {} explicitly assumes original intact linear duplex terminal groups: left {:?}, right {:?}. This is an operator hypothesis, not measured chemistry, synthesis, nick repair or recovery of imported physical state. Original DNA/features/citations are retained; positional modification states remain unknown.",index+1,left,right));
                    next_physical = Some(tracked);
                    product
                }
                Operation::DeclareSourceModifications {
                    template,
                    assumptions,
                } => {
                    rule.set(RuleId::SourceChemistry);
                    if !matches!(template, Source::Input(_)) {
                        return Err(Error::Unsupported("Source-generation assumptions may only bind an original input, never rescan a reconstructed descendant".into()));
                    }
                    let substrate = resolve(template)?;
                    let tracked = TrackedDuplex::from_input(substrate, assumptions)?;
                    let mut product = tracked
                        .duplex
                        .molecule(format!("Declared source generation {}", index + 1))?;
                    rule.set(RuleId::AnnotationProjection);
                    project_into(
                        substrate,
                        &[MappingSegment {
                            source_start: 0,
                            source_end: substrate.sequence.len(),
                            destination_start: 0,
                            reverse: false,
                        }],
                        &mut product,
                        findings,
                    )?;
                    findings.push(format!("Step {} records a closed Dam/Dcm/CpG source-generation assumption on exact input revision {}: {}. This does not measure modifications or infer any other methylase activity.",index+1,substrate.revision(),serde_json::to_string(assumptions).map_err(|e|Error::Format(e.to_string()))?));
                    next_physical = Some(tracked);
                    product
                }
                Operation::CheckedDigest {
                    template,
                    selection,
                    options,
                    conditions,
                    fragment_index,
                } => {
                    rule.set(RuleId::RestrictionProfile);
                    options.validate()?;
                    let tracked = resolve_tracked(template)?;
                    if let Some(options) = conditions {
                        let substrate = resolve(template)?;
                        let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
                        rule.set(RuleId::ReactionConditions);
                        let evidence =
                            crate::cloning::restriction::conditions::source::SourceReport::assess(
                                &catalog, &substrate, &tracked, selection, options,
                            )?;
                        rule.set(RuleId::EvidenceBudget);
                        step_evidence_bytes = step_evidence_bytes
                            .checked_add(evidence_bytes(
                                &evidence,
                                crate::cloning::restriction::conditions::MAX_REPORT_BYTES + 1024,
                            )?)
                            .filter(|n| *n <= 32 * 1024 * 1024)
                            .ok_or(Error::TooLarge)?;
                        rule.set(operation_rule(step));
                        if !evidence.satisfies_request() {
                            incomplete.push(format!("Step {} has an unresolved declared documentary reaction-condition constraint: {:?}; reconstructed sequence products remain independent of this evidence.", index + 1, evidence.report.documentary_status));
                        }
                        step_reports.push(StepReport {
                            step_index: index,
                            topo: None,
                            reaction_conditions: Some(evidence),
                        });
                    }
                    rule.set(RuleId::RestrictionCleavage);
                    let report = modified_duplex::checked_digest_with_budget(
                        &tracked,
                        &CheckedDigestOptions {
                            selection: selection.clone(),
                            max_comparisons: usize::try_from(options.max_comparisons)
                                .map_err(|_| Error::TooLarge)?,
                            max_sites: options.max_sites,
                        },
                        cleavage_budget,
                    )?;
                    evidence_bytes(&report, options.max_report_bytes)?;
                    if report.fragments.len() > options.max_fragments {
                        return Err(Error::TooLarge);
                    }
                    for site in report.sites.iter().take(100) {
                        findings.push(format!(
                            "Step {} checked {} site at source base {}: {:?}. {}",
                            index + 1,
                            site.enzyme,
                            site.recognition_start + 1,
                            site.outcome,
                            site.reason
                        ));
                    }
                    findings.extend(report.warnings.iter().take(100).cloned());
                    if !report.complete || !report.products_determined {
                        return Err(Error::Unsupported(format!(
                            "Checked digest cannot establish complete products under the pinned positional profile; unknown/impaired or bounded evidence remains. {}",
                            report.limits.join(" ")
                        )));
                    }
                    rule.set(RuleId::DigestSelection);
                    let fragment=report.fragments.get(*fragment_index).ok_or_else(||Error::Format(format!("checked digest has {} fragments; selected zero-based index {fragment_index} is unavailable",report.fragments.len())))?;
                    let mut product = fragment.duplex.duplex.molecule(format!(
                        "Checked digest {} · fragment {}",
                        index + 1,
                        fragment_index + 1
                    ))?;
                    let substrate = resolve(template)?;
                    rule.set(RuleId::AnnotationProjection);
                    project_into(
                        substrate,
                        &region_map(
                            substrate.sequence.len(),
                            fragment.top_start,
                            fragment.duplex.duplex.top.len(),
                            0,
                        )?,
                        &mut product,
                        findings,
                    )?;
                    findings.push(format!("Step {} selected fragment {} of {} fully determined positional-profile digest products; supplier conditions and experimental completion remain assumptions.",index+1,fragment_index+1,report.fragments.len()));
                    next_physical = Some(fragment.duplex.clone());
                    product
                }
                Operation::Digest {
                    template,
                    fragment_index,
                    ..
                }
                | Operation::CatalogDigest {
                    template,
                    fragment_index,
                    ..
                } => {
                    rule.set(RuleId::RestrictionProfile);
                    let source = resolve_tracked(template)?;
                    let enzymes = match step {
                        Operation::CatalogDigest {
                            selection,
                            conditions,
                            ..
                        } => {
                            let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
                            if let Some(options) = conditions {
                                let substrate = resolve(template)?;
                                rule.set(RuleId::ReactionConditions);
                                let evidence = crate::cloning::restriction::conditions::source::SourceReport::assess(&catalog, &substrate, &source, selection, options)?;
                                rule.set(RuleId::EvidenceBudget);
                                step_evidence_bytes = step_evidence_bytes
                                    .checked_add(evidence_bytes(
                                        &evidence,
                                        crate::cloning::restriction::conditions::MAX_REPORT_BYTES
                                            + 1024,
                                    )?)
                                    .filter(|n| *n <= 32 * 1024 * 1024)
                                    .ok_or(Error::TooLarge)?;
                                rule.set(operation_rule(step));
                                if !evidence.satisfies_request() {
                                    incomplete.push(format!("Step {} has an unresolved declared documentary reaction-condition constraint: {:?}; reconstructed sequence products remain independent of this evidence.", index + 1, evidence.report.documentary_status));
                                }
                                step_reports.push(StepReport {
                                    step_index: index,
                                    topo: None,
                                    reaction_conditions: Some(evidence),
                                });
                            }
                            let resolved =
                                selection.resolved_enzymes(&catalog).map_err(|error| {
                                    Error::Format(format!("Catalog digest pin rejected: {error}"))
                                })?;
                            findings.push(format!("Step {} uses explicit embedded catalog edition {} / content {} and exact ordered enzyme IDs {}. Selected identity/geometry is validated against the catalog; documentary supplier conditions do not establish complete cleavage.",index+1,selection.catalog_edition,selection.catalog_content_hash,selection.enzymes.iter().map(|e|e.id.as_str()).collect::<Vec<_>>().join(", ")));
                            resolved
                        }
                        Operation::Digest { enzymes, .. } => {
                            let catalog = analysis::starter_enzymes();
                            enzymes
                                .iter()
                                .map(|name| {
                                    catalog
                                        .iter()
                                        .find(|e| e.name.eq_ignore_ascii_case(name))
                                        .cloned()
                                        .ok_or_else(|| {
                                            Error::Unsupported(format!(
                                                "enzyme {name} is not in the embedded catalog"
                                            ))
                                        })
                                })
                                .collect::<Result<Vec<_>, _>>()?
                        }
                        _ => unreachable!("matched digest operation"),
                    };
                    if !source.generation.is_empty() {
                        let state = source.summary()?;
                        if state.known_modified > 0 || state.unknown_eligible > 0 {
                            incomplete.push(format!("Step {} uses uncatalogued assumed-complete cleavage on explicitly tracked modified/unknown material. Use CheckedDigest to establish cleavage under the positional profile.",index+1));
                        }
                    }
                    rule.set(RuleId::RestrictionCleavage);
                    let fragments = modified_duplex::digest(&source, &enzymes)?;
                    rule.set(RuleId::DigestSelection);
                    let fragment=fragments.get(*fragment_index).ok_or_else(||Error::Format(format!("digest has {} fragments; selected zero-based index {fragment_index} does not exist",fragments.len())))?;
                    let mut product = fragment.duplex.molecule(format!(
                        "Digest {} · fragment {}",
                        index + 1,
                        fragment_index + 1
                    ))?;
                    let template_molecule = resolve(template)?;
                    rule.set(RuleId::AnnotationProjection);
                    let map = region_map(
                        template_molecule.sequence.len(),
                        fragment.top_start,
                        fragment.duplex.top.len(),
                        0,
                    )?;
                    rule.set(RuleId::AnnotationProjection);
                    project_into(template_molecule, &map, &mut product, findings)?;
                    findings.push(format!("Step {} selects fragment {} of {} from the complete digest (top reference boundary {}); top/bottom lengths {}/{}.",index+1,fragment_index+1,fragments.len(),fragment.top_start,fragment.duplex.top.len(),fragment.duplex.bottom.len()));
                    next_physical = Some(fragment.duplex.clone());
                    product
                }
                Operation::TreatEnds {
                    template,
                    treatment,
                } => {
                    rule.set(RuleId::EndTreatment);
                    let tracked = resolve_tracked(template)?.treat(*treatment)?;
                    let duplex = &tracked.duplex;
                    let mut product = duplex.molecule(format!("End treatment {}", index + 1))?;
                    let source = resolve(template)?;
                    rule.set(RuleId::AnnotationProjection);
                    project_into(
                        source,
                        &[MappingSegment {
                            source_start: 0,
                            source_end: source.sequence.len(),
                            destination_start: 0,
                            reverse: false,
                        }],
                        &mut product,
                        findings,
                    )?;
                    next_physical = Some(tracked);
                    product
                }
                Operation::Tail3Prime { template, base } => {
                    rule.set(RuleId::TailAddition);
                    let source = resolve_physical(template)?;
                    let report =
                        crate::cloning::ligation::tailing::tail_blunt_duplex(&source, *base)?;
                    let product = tailed_document(
                        resolve(template)?,
                        &report,
                        format!("3′ tailing {}", index + 1),
                        findings,
                    )?;
                    findings.push(format!("Step {} explicitly appends one {} nucleotide to both strand 3′ termini. Added nucleotides have no source-template ancestry; five-prime chemistry and internal nicks are retained.",index+1,base.as_str()));
                    findings.extend(report.warnings.iter().take(100).cloned());
                    next_physical = Some(resolve_tracked(template)?.tail(*base)?);
                    product
                }
                Operation::ProfiledOverlap {
                    fragments,
                    overlaps,
                    circular,
                    options,
                } => {
                    rule.set(RuleId::HomologyPolicy);
                    let mut sources = Vec::new();
                    let mut reverse = Vec::new();
                    for oriented in fragments {
                        let source = resolve(&oriented.source)?;
                        if let Source::Input(_) = &oriented.source {
                            if chemistry::has_unmapped_physical_evidence(source) {
                                return Err(Error::Unsupported("profiled overlap assembly of imported physical evidence requires reconstruction of its originating graph".into()));
                            }
                        }
                        if let Source::Step(i) = &oriented.source {
                            if let Some(Some(d)) = physical.get(*i) {
                                let ends = d.ends()?;
                                if d.left.three_prime == chemistry::ThreePrime::Blocked
                                    || d.right.three_prime == chemistry::ThreePrime::Blocked
                                {
                                    return Err(Error::Unsupported("profiled overlap assembly of a known blocked3′ terminus requires an explicitly modeled unblocking/resection mechanism".into()));
                                }
                                if !ends.0.overhang.is_empty()
                                    || !ends.1.overhang.is_empty()
                                    || !d.nicks.is_empty()
                                {
                                    return Err(Error::Unsupported("profiled overlap assembly of staggered or nicked physical substrates is not modeled".into()));
                                }
                            }
                        }
                        sources.push(source.clone());
                        reverse.push(oriented.reverse);
                    }
                    let report = crate::cloning::overlap::profiles::check_oriented_fragments(
                        &sources, &reverse, overlaps, *circular, options,
                    )?;
                    if !report.geometry_valid {
                        return Err(Error::Format(format!(
                            "profiled overlap geometry is invalid: {}",
                            report.findings.join(" ")
                        )));
                    }
                    if !report.complete {
                        incomplete.push(format!("Step {} profile repeat scan reached its configured bound; alternative placements remain unsearched.",index+1));
                    }
                    findings.extend(report.findings.iter().take(100).cloned());
                    findings.extend(report.warnings.iter().take(100).cloned());
                    findings.push(format!("Step {} uses {:?} sequence recommendations. Recommendation compliance {:?} is advisory and does not establish reagent activity, experimental outcome or uniqueness.",index+1,options.method,report.recommendation_compliant));
                    let mut oriented_sources = Vec::new();
                    for (source, reverse) in sources.iter().zip(&reverse) {
                        if *reverse {
                            let mut oriented = Molecule::dna(
                                format!("Oriented {}", source.name),
                                source.sequence.reverse_complement()?.as_str(),
                                source.topology,
                            )?;
                            rule.set(RuleId::AnnotationProjection);
                            project_into(
                                source,
                                &[MappingSegment {
                                    source_start: 0,
                                    source_end: source.sequence.len(),
                                    destination_start: 0,
                                    reverse: true,
                                }],
                                &mut oriented,
                                findings,
                            )?;
                            carry_unmapped_evidence(source, &mut oriented);
                            oriented_sources.push(oriented);
                        } else {
                            oriented_sources.push(source.clone());
                        }
                    }
                    let refs = oriented_sources.iter().collect::<Vec<_>>();
                    let mut product =
                        assemble_overlap(&refs, overlaps, *circular, index, findings)?;
                    product.metadata.insert(
                        "overlap_profile:options".into(),
                        serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?,
                    );
                    product.metadata.insert(
                        "overlap_profile:source_revisions".into(),
                        serde_json::to_string(&report.source_revisions)
                            .map_err(|e| Error::Format(e.to_string()))?,
                    );
                    product.metadata.insert(
                        "overlap_profile:orientations".into(),
                        serde_json::to_string(&reverse)
                            .map_err(|e| Error::Format(e.to_string()))?,
                    );
                    product
                }
                Operation::RepairEnds { template, method } => {
                    rule.set(RuleId::EndRepair);
                    let source = resolve_physical(template)?;
                    rule.set(RuleId::EndRepair);
                    let repair = chemistry::repair_ends(&source, *method)?;
                    let product = repaired_document(
                        resolve(template)?,
                        &repair,
                        format!("End repair {}", index + 1),
                        findings,
                    )?;
                    findings.push(format!("Step {} end repair {:?}: top added/removed {}/{} bases; bottom added/removed {}/{} bases. Five-prime chemistry and internal nicks are retained.",index+1,method,repair.top_added_3prime.len(),repair.top_removed_3prime.len(),repair.bottom_added_3prime.len(),repair.bottom_removed_3prime.len()));
                    findings.extend(repair.warnings.iter().take(100).cloned());
                    next_physical = Some(resolve_tracked(template)?.repair(*method)?);
                    product
                }
                Operation::AnnealOligos {
                    top,
                    bottom,
                    options,
                    candidate_index,
                } => {
                    rule.set(RuleId::OligoAnnealing);
                    if !matches!(top, Source::Input(_)) || !matches!(bottom, Source::Input(_)) {
                        return Err(Error::Unsupported("annealing requires explicitly declared input oligos; duplex-producing steps do not implicitly denature into ssDNA".into()));
                    }
                    let top = resolve(top)?;
                    let bottom = resolve(bottom)?;
                    rule.set(RuleId::OligoAnnealing);
                    let preview = simulate_oligo_annealing(top, bottom, options)?;
                    let candidate = preview
                        .annealing
                        .candidates
                        .get(*candidate_index)
                        .ok_or_else(|| {
                            Error::Format("selected annealing candidate is unavailable".into())
                        })?;
                    let product = annealed_document_inner(
                        top,
                        bottom,
                        &preview.annealing,
                        *candidate_index,
                        findings,
                    )?;
                    findings.push(format!("Step {} anneals declared 5′→3′ ssDNA oligos at bottom offset {}, with {} paired bases. Explicit terminal chemistry is a plan assumption; annealing creates no covalent bonds.",index+1,candidate.duplex.bottom_start,candidate.paired_bases));
                    if !preview.annealing.complete {
                        incomplete.push(format!("Step {} annealing placement search was incomplete; unsearched alignments cannot establish uniqueness.",index+1));
                    }
                    if !preview.annealing.unique {
                        incomplete.push(format!("Step {} annealing has competing or unproven placements; an explicit selected candidate does not establish a unique annealed product.",index+1));
                    }
                    findings.extend(preview.annealing.warnings.iter().take(100).cloned());
                    next_physical = Some(TrackedDuplex::annealed_unmodified(&candidate.duplex)?);
                    product
                }
                Operation::Ligate {
                    fragments,
                    circular,
                    allow_nicked,
                    options,
                } => {
                    rule.set(RuleId::LigationGeometry);
                    if fragments.is_empty() || fragments.len() > 1000 {
                        return Err(Error::Format(
                            "ligation requires1–1000 explicitly oriented fragments".into(),
                        ));
                    }
                    rule.set(RuleId::LigationAlternatives);
                    options.validate()?;
                    rule.set(RuleId::LigationGeometry);
                    if fragments.len() > crate::cloning::ligation::outcomes::MAX_FRAGMENTS {
                        return Err(Error::TooLarge);
                    }
                    let mut raw = Vec::with_capacity(fragments.len());
                    let mut oriented_lengths = Vec::with_capacity(fragments.len());
                    let mut bases = 0usize;
                    for source in fragments {
                        let duplex = resolve_tracked(&source.source)?;
                        bases = bases
                            .saturating_add(duplex.top.len())
                            .saturating_add(duplex.bottom.len());
                        if bases > 2 * MAX_SEQUENCE_LEN {
                            return Err(Error::TooLarge);
                        }
                        oriented_lengths.push(if source.reverse {
                            duplex.bottom.len()
                        } else {
                            duplex.top.len()
                        });
                        raw.push(duplex);
                    }
                    let (duplex, outcomes_receipt) = {
                        rule.set(RuleId::EvidenceBudget);
                        if retained_processing + physical_processing_reports.len() >= 64 {
                            return Err(Error::TooLarge);
                        }
                        rule.set(operation_rule(step));
                        let keys = (0..fragments.len())
                            .map(|i| format!("fragment-{}", i + 1))
                            .collect::<Vec<_>>();
                        let operands =
                            keys.iter()
                                .zip(&raw)
                                .map(|(key, material)| {
                                    crate::cloning::ligation::outcomes::Operand { key, material }
                                })
                                .collect::<Vec<_>>();
                        let request = crate::cloning::ligation::outcomes::Request {
                            intended: keys
                                .iter()
                                .zip(fragments)
                                .map(|(key, source)| {
                                    crate::cloning::ligation::outcomes::Placement {
                                        key: key.clone(),
                                        reverse: source.reverse,
                                    }
                                })
                                .collect(),
                            circular: *circular,
                            allow_nicked: *allow_nicked,
                            options: options.clone(),
                        };
                        rule.set(RuleId::LigationChemistry);
                        let (outcomes, actual) =
                            crate::cloning::ligation::outcomes::analyze_fresh(&operands, &request)?
                                .into_parts();
                        let bindings = keys
                            .iter()
                            .zip(fragments)
                            .map(|(key, source)| {
                                Ok(LigationOutcomeBinding {
                                    key: key.clone(),
                                    source: source.source.clone(),
                                    source_revision: resolve(&source.source)?.revision(),
                                    material_revision: outcomes
                                        .sources
                                        .iter()
                                        .find(|s| &s.key == key)
                                        .ok_or_else(|| {
                                            Error::Format(
                                                "Missing native ligation source identity".into(),
                                            )
                                        })?
                                        .state_revision
                                        .clone(),
                                })
                            })
                            .collect::<Result<Vec<_>, Error>>()?;
                        let evidence = LigationOutcomeEvidence {
                            bindings,
                            report: outcomes,
                        };
                        rule.set(RuleId::LigationAlternatives);
                        evidence.validate_for_operation(
                            fragments,
                            *circular,
                            *allow_nicked,
                            options,
                        )?;
                        rule.set(RuleId::EvidenceBudget);
                        step_evidence_bytes = step_evidence_bytes.checked_add(evidence_bytes(&evidence,32*1024*1024)?).filter(|n|*n<=32*1024*1024).ok_or_else(||Error::Unsupported("Complete ligation outcomes plus prior processing evidence exceed the32MiB plan evidence bound; reduce explicit outcome bounds before capturing this plan. The intended chemistry was not rejected and no alternatives were silently dropped.".into()))?;
                        rule.set(operation_rule(step));
                        let outcomes_receipt = (
                            crate::encoding::bounded_json::hash(&evidence)?,
                            evidence.report.context_revision.clone(),
                        );
                        findings.push(format!("Step {} retains the declared intended ligation and {} alternative one-copy paths; enumeration complete {}, exhaustive within that scope {}. Alternative existence does not invalidate the intended path or establish isolation; other copy numbers and incomplete assemblies remain outside this scope.",index+1,evidence.report.alternatives.len(),evidence.report.enumeration_complete,evidence.report.exhaustive_within_scope));
                        physical_processing_reports.push(PhysicalProcessingStepReport {
                            step_index: index,
                            evidence: PhysicalProcessingEvidence::LigationOutcomes(evidence),
                        });
                        (actual, outcomes_receipt)
                    };
                    findings.push(format!("Step {} ligation reconstructs both strands with {} retained nicks; covalently closed circular duplex in this model: {}. No nick repair is assumed.",index+1,duplex.nicks.len(),duplex.topology==Topology::Circular && duplex.nicks.is_empty()));
                    let mut product = duplex.molecule(format!("Ligation {}", index + 1))?;
                    {
                        let (hash, context) = outcomes_receipt;
                        product
                            .metadata
                            .insert("source:ligation_outcomes:evidence_hash".into(), hash);
                        product
                            .metadata
                            .insert("source:ligation_outcomes:context".into(), context);
                    }
                    let mut destination = 0usize;
                    for (oriented_source, oriented_top_length) in
                        fragments.iter().zip(&oriented_lengths)
                    {
                        let source = resolve(&oriented_source.source)?;
                        let map = if oriented_source.reverse {
                            let original = resolve_physical(&oriented_source.source)?;
                            let end = original.bottom_start + original.bottom.len() as i64;
                            reverse_reference_map(
                                source.sequence.len(),
                                end,
                                destination,
                                product.sequence.len(),
                                product.topology,
                            )?
                        } else {
                            vec![MappingSegment {
                                source_start: 0,
                                source_end: source.sequence.len(),
                                destination_start: destination,
                                reverse: false,
                            }]
                        };
                        rule.set(RuleId::AnnotationProjection);
                        project_into(source, &map, &mut product, findings)?;
                        destination += *oriented_top_length;
                    }
                    next_physical = Some(duplex);
                    product
                }
            };
            let is_reference_only = diagnostic_reference_only
                || inherited_reference_only
                || matches!(
                    step,
                    Operation::RecombineReference { .. }
                        | Operation::MultiSiteRecombineReference { .. }
                        | Operation::Topo { .. }
                        | Operation::CheckedTopo { .. }
                );
            rule.set(RuleId::ReferenceOnlyAuthority);
            if is_reference_only {
                if next_physical.is_some() {
                    return Err(Error::Unsupported("Reference-only ancestry cannot acquire chemical authority without a modeled materialization assumption".into()));
                }
                // Restrictive portable evidence only. Actual graph authority lives in the private flag.
                product.metadata.insert(
                    crate::cloning::recombination::REFERENCE_ONLY_KEY.into(),
                    "true".into(),
                );
                product.metadata.insert(
                    crate::cloning::recombination::CHEMISTRY_UNMAPPED_KEY.into(),
                    "true".into(),
                );
            }
            rule.set(RuleId::ProductIntegrity);
            if let Some(tracked) = &next_physical {
                let summary = tracked.summary()?;
                product.metadata.insert(
                    "tracked_modifications".into(),
                    serde_json::to_string(&summary).map_err(|e| Error::Format(e.to_string()))?,
                );
            }
            rule.set(RuleId::ProductIntegrity);
            product.validate()?;
            rule.set(RuleId::ProductBudget);
            product_bases = product_bases
                .saturating_add(product.sequence.len())
                .saturating_add(
                    next_physical
                        .as_ref()
                        .map_or(0, |d| d.top.len() + d.bottom.len()),
                );
            if product_bases > 32_000_000 {
                return Err(Error::Unsupported(
                    "plan products exceed the 32 million base materialization budget".into(),
                ));
            }
            product
                .metadata
                .insert("plan_revision".into(), plan_hash.into());
            product
                .metadata
                .insert("engine_version".into(), ENGINE_VERSION.into());
            rule.set(RuleId::ProductBudget);
            product_bytes =
                product_bytes.saturating_add(evidence_bytes(&product, 64 * 1024 * 1024)?);
            if product_bytes.saturating_add(step_evidence_bytes) > 64 * 1024 * 1024 {
                return Err(Error::Unsupported(
                    "plan products exceed the 64 MiB serialized annotation/evidence budget".into(),
                ));
            }
            Ok((product, next_physical, is_reference_only))
        })();
        budget.product_bases = product_bases;
        budget.product_bytes = product_bytes;
        budget.step_evidence_bytes = step_evidence_bytes;
        Self {
            result,
            rule: rule.get(),
            notes,
            incomplete,
            step_reports,
            physical_processing_reports,
            reference_step_reports,
            golden_gate,
        }
    }
}

pub(super) fn operation_rule(step: &Operation) -> RuleId {
    match step {
        Operation::ResectThreePrime { .. } => RuleId::EndRepair,
        Operation::HomologyAnnealing { .. } => RuleId::HomologyReconstruction,
        Operation::HomologyAssembly { .. } => RuleId::HomologyReconstruction,
        Operation::OverlapExtension { .. } => RuleId::OverlapGeometry,
        Operation::ReconstructMaterial { .. } => RuleId::MaterialOrigin,
        Operation::ReconstructRunProduct { .. } => RuleId::SavedMaterial,
        Operation::MultiSiteRecombineReference { .. } => RuleId::RecombinationReference,
        Operation::TopoActions { .. } => RuleId::TopoReconstruction,
        Operation::Topo { .. } => RuleId::TopoReconstruction,
        Operation::RecombineReference { .. } => RuleId::RecombinationReference,
        Operation::Pcr { .. } => RuleId::PcrBinding,
        Operation::PcrWithMismatches { .. } => RuleId::MismatchPcrPolicy,
        Operation::GoldenGate { .. } => RuleId::GoldenGateAssembly,
        Operation::Overlap { .. } => RuleId::OverlapGeometry,
        Operation::DeclareSourceTermini { .. } => RuleId::SourceChemistry,
        Operation::DeclareSourceModifications { .. } => RuleId::SourceChemistry,
        Operation::CheckedDigest { .. } => RuleId::RestrictionProfile,
        Operation::Digest { .. } => RuleId::RestrictionProfile,
        Operation::TreatEnds { .. } => RuleId::EndTreatment,
        Operation::Tail3Prime { .. } => RuleId::TailAddition,
        Operation::ProfiledOverlap { .. } => RuleId::HomologyPolicy,
        Operation::RepairEnds { .. } => RuleId::EndRepair,
        Operation::AnnealOligos { .. } => RuleId::OligoAnnealing,
        Operation::Ligate { .. } => RuleId::LigationGeometry,
        Operation::CheckedTopo { .. } => RuleId::TopoReconstruction,
        Operation::CatalogDigest { .. } => RuleId::RestrictionProfile,
        Operation::SelectedGoldenGate { .. } => RuleId::GoldenGateAssembly,
    }
}
