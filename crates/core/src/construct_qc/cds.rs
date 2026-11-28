//! CDS facts reuse the one current translation owner, including exceptions and
//! INSDC biological-order extraction. No alternative codon decoder is introduced.
use super::*;
use translation::{InitiationPolicy, TranslationOptions};
pub(super) struct CodingFrame {
    pub offset: usize,
    pub protein: String,
}

fn unresolved(checks: &mut Checks<'_>, i: usize, reason: &str) -> Result<(), Error> {
    for rule in [
        Rule::CdsFrame,
        Rule::CdsStart,
        Rule::CdsTerminalStop,
        Rule::CdsInternalStop,
        Rule::StoredTranslation,
    ] {
        checks.push(rule, Verdict::Indeterminate, &[i], reason)?;
    }
    Ok(())
}
pub(super) fn inspect(checks: &mut Checks<'_>) -> Result<(), Error> {
    let m = checks.m;
    let o = checks.o;
    let indices = m
        .features
        .iter()
        .enumerate()
        .filter(|(_, f)| f.kind.eq_ignore_ascii_case("CDS"))
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    if indices.is_empty() {
        checks.push(
            Rule::AnnotationBasis,
            Verdict::Indeterminate,
            &[],
            "No CDS role is annotated; this does not establish absence of coding content.",
        )?;
    }
    for i in indices {
        if !checks.assumed(i) {
            checks.push(
                Rule::AnnotationBasis,
                Verdict::Indeterminate,
                &[i],
                "CDS role/coordinates lack a current explicit annotation assumption or require review; imported feature labels alone cannot establish the role.",
            )?;
            unresolved(
                checks,
                i,
                "CDS checks depend on a current declared annotation assumption; no structural pass is issued for unreviewed or undeclared annotation.",
            )?;
            continue;
        }
        checks.push(
            Rule::AnnotationBasis,
            Verdict::Informational,
            &[i],
            "CDS role and coordinates are explicitly assumed for these checks. This declaration does not verify biological annotation truth or functional identity.",
        )?;
        let f = &m.features[i];
        let mut effective = f.clone();
        if !effective
            .qualifiers
            .iter()
            .any(|(k, _)| k == "transl_table")
        {
            effective
                .qualifiers
                .push(("transl_table".into(), Some(o.genetic_code.to_string())));
        }
        let translated =
            match translation::translate_cds_validated(m, &effective, &checks.source_revision) {
                Ok(t) => t,
                Err(error) => {
                    unresolved(
                        checks,
                        i,
                        &format!("Existing CDS translator cannot resolve this annotation: {error}"),
                    )?;
                    continue;
                }
            };
        let partial =
            f.qualifiers.iter().any(|(k, _)| k == "partial") || translated.codon_start != 1;
        if !partial && translated.translation.trailing_bases == 0 {
            checks.frames.insert(
                i,
                CodingFrame {
                    offset: usize::from(translated.codon_start - 1),
                    protein: translated.translation.protein.clone(),
                },
            );
        }
        if partial {
            checks.push(
                Rule::CdsFrame,
                Verdict::Indeterminate,
                &[i],
                "Partial CDS or codon_start offset does not assert a complete coding interval.",
            )?;
        } else {
            checks.push(
                Rule::CdsFrame,
                if translated.translation.trailing_bases == 0 {
                    Verdict::Pass
                } else {
                    Verdict::Fail
                },
                &[i],
                &format!(
                    "Declared CDS has {} trailing base(s) outside complete codons in genetic code {}.",
                    translated.translation.trailing_bases,
                    translated.translation.options.genetic_code,
                ),
            )?;
        }
        let oriented = f.location.extract(&m.sequence)?;
        let offset = usize::from(translated.codon_start - 1).min(oriented.len());
        let codons = (oriented.len() - offset) / 3;
        if codons == 0 || partial {
            for rule in [Rule::CdsStart, Rule::CdsTerminalStop, Rule::CdsInternalStop] {
                checks.push(
                    rule,
                    Verdict::Indeterminate,
                    &[i],
                    "Complete CDS boundaries are not declared or contain no complete codon; full-CDS integrity is unresolved.",
                )?;
            }
        } else {
            let first = Sequence::new(m.sequence.alphabet(), &oriented[offset..offset + 3])?;
            let start = translation::translate(
                &first,
                &TranslationOptions {
                    genetic_code: translated.translation.options.genetic_code,
                    frame: 0,
                    reverse: false,
                    initiation: InitiationPolicy::RequireRecognizedStart,
                },
            );
            let (verdict, reason) = match start {
                Ok(_) => (
                    Verdict::Pass,
                    "Every expansion of the first codon is a recognized initiator in the declared genetic code.",
                ),
                Err(Error::Unsupported(_)) => (
                    Verdict::Indeterminate,
                    "The first codon can be an initiator, but ambiguity prevents definite initiation.",
                ),
                Err(_) => (
                    Verdict::Fail,
                    "The first codon is not a recognized initiator in the declared genetic code.",
                ),
            };
            checks.push(Rule::CdsStart, verdict, &[i], reason)?;
            let terminal_offset = offset + 3 * (codons - 1);
            let terminal = Sequence::new(
                m.sequence.alphabet(),
                &oriented[terminal_offset..terminal_offset + 3],
            )?;
            let terminal = translation::translate(
                &terminal,
                &TranslationOptions {
                    genetic_code: translated.translation.options.genetic_code,
                    ..Default::default()
                },
            )?;
            let terminal_verdict = if translated.translation.trailing_bases > 0 {
                Verdict::Indeterminate
            } else if terminal.protein == "*" {
                Verdict::Pass
            } else if terminal.possible_stop_codons > 0 {
                Verdict::Indeterminate
            } else {
                Verdict::Fail
            };
            checks.push(
                Rule::CdsTerminalStop,
                terminal_verdict,
                &[i],
                "The final complete codon is tested by the existing translation owner; context-dependent or ambiguous termination is unresolved. Trailing bases are reported separately.",
            )?;
            let internal_possible = translated
                .translation
                .possible_stop_codons
                .saturating_sub(terminal.possible_stop_codons);
            let internal_verdict = if translated.internal_stops > 0 {
                Verdict::Fail
            } else if internal_possible > 0 {
                Verdict::Indeterminate
            } else {
                Verdict::Pass
            };
            checks.push(
                Rule::CdsInternalStop,
                internal_verdict,
                &[i],
                &format!(
                    "{} definite internal stop(s) and {internal_possible} potentially terminating internal codon(s) in the declared CDS.",
                    translated.internal_stops,
                ),
            )?;
        }
        let (verdict, reason) = match translated.stored_translation_matches {
            None => (
                Verdict::Indeterminate,
                "No stored /translation qualifier is available for a literal comparison.",
            ),
            Some(true) if translated.translation.unresolved_codons == 0 => (
                Verdict::Pass,
                "Stored /translation literally agrees with the independently computed CDS protein (one definite terminal stop omitted).",
            ),
            Some(true) => (
                Verdict::Indeterminate,
                "Stored /translation agrees with unresolved computed calls; agreement is not a resolved protein identity.",
            ),
            Some(false) => (
                Verdict::Informational,
                "Stored /translation differs literally from the independently computed CDS protein; this is a documentary disagreement, not a biological outcome prediction.",
            ),
        };
        checks.push(Rule::StoredTranslation, verdict, &[i], reason)?;
    }
    Ok(())
}
