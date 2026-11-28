//! Human coordinate entry and exact-event guards for IVT template inspection.
use collagene_core::{
    material::lineage::EventBinding,
    rna_chemistry::{Provenance, ivt},
};
use collagene_model::Molecule;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Draft {
    pub first_base: String,
    pub length: String,
    pub reverse: bool,
    pub cap: ivt::CapChoice,
    pub initiator: ivt::InitiatorRequirement,
    pub coding_cut: String,
    pub template_cut: String,
    pub poly_a_first: String,
    pub poly_a_last: String,
}
impl Default for Draft {
    fn default() -> Self {
        Self {
            first_base: String::new(),
            length: String::new(),
            reverse: false,
            cap: ivt::CapChoice::NotSpecified,
            initiator: ivt::InitiatorRequirement::NoAdditionalRequirement,
            coding_cut: String::new(),
            template_cut: String::new(),
            poly_a_first: String::new(),
            poly_a_last: String::new(),
        }
    }
}
fn number(text: &str, label: &str) -> Result<usize, String> {
    let text = text.trim();
    if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("Enter a whole number for {label}."));
    }
    text.parse()
        .map_err(|_| format!("{label} exceeds the supported range."))
}
fn one_based(text: &str, label: &str) -> Result<usize, String> {
    number(text, label)?
        .checked_sub(1)
        .ok_or_else(|| format!("{label} starts at 1."))
}
fn pair_present(first: &str, last: &str, label: &str) -> Result<bool, String> {
    match (first.trim().is_empty(), last.trim().is_empty()) {
        (true, true) => Ok(false),
        (false, false) => Ok(true),
        _ => Err(format!("Enter both {label} values, or leave both empty.")),
    }
}
impl Draft {
    pub fn request(&self, source: &Molecule) -> Result<ivt::Request, String> {
        ivt::source_preflight(source).map_err(|error| error.to_string())?;
        let first = one_based(&self.first_base, "first transcribed source base")?;
        let length = number(&self.length, "transcript length")?;
        let provenance = Provenance {
            statement: "User-declared IVT coding strand, transcript interval and template requirements; supplied cut offsets are geometry assumptions, not observed cleavage.".into(),
            reference: None,
        };
        let linearization = if pair_present(&self.coding_cut, &self.template_cut, "cut offset")? {
            Some(ivt::Linearization {
                coding_cut_offset: number(&self.coding_cut, "coding-strand cut offset")?,
                template_cut_offset: number(&self.template_cut, "template-strand cut offset")?,
                provenance: provenance.clone(),
            })
        } else {
            None
        };
        let poly_a = if pair_present(&self.poly_a_first, &self.poly_a_last, "poly(A) interval")? {
            Some(ivt::PolyARegion {
                start_offset: one_based(&self.poly_a_first, "first poly(A) transcript base")?,
                // An inclusive one-based end is the corresponding zero-based exclusive end.
                end_offset: number(&self.poly_a_last, "last poly(A) transcript base")?,
            })
        } else {
            None
        };
        let profile = ivt::profiles()
            .map_err(|e| e.to_string())?
            .into_iter()
            .next()
            .ok_or("No sourced IVT template profile is installed.")?
            .reference;
        let request = ivt::Request {
            source_revision: source.revision(),
            profile,
            transcript: ivt::Transcript {
                start: first,
                length,
                coding_strand: if self.reverse {
                    ivt::CodingStrand::Reverse
                } else {
                    ivt::CodingStrand::Forward
                },
            },
            cap: self.cap,
            initiator: self.initiator,
            linearization,
            poly_a,
            provenance,
        };
        request.validate_for(source).map_err(|e| e.to_string())?;
        Ok(request)
    }
}

#[derive(Clone)]
pub(crate) struct Capture {
    pub binding: EventBinding,
    pub epoch: u64,
    pub draft: Draft,
    pub request: ivt::Request,
    pub revision: String,
}
impl Capture {
    pub fn prepare(
        binding: EventBinding,
        epoch: u64,
        draft: Draft,
        source: &Molecule,
    ) -> Result<Self, String> {
        // Draft conversion performs the shared borrowed source admission before
        // constructing its revision pin; never hash an unadmitted source here.
        let request = draft.request(source)?;
        if binding.revision != request.source_revision {
            return Err("The IVT source event changed; inspect the current document again.".into());
        }
        let revision = ivt::revision(source, &request).map_err(|e| e.to_string())?;
        Ok(Self {
            binding,
            epoch,
            draft,
            request,
            revision,
        })
    }
    pub fn guard(&self, binding: &EventBinding, epoch: u64, draft: &Draft) -> Result<(), String> {
        if self.binding != *binding || self.epoch != epoch || self.draft != *draft {
            return Err(
                "The document, history event or IVT settings changed; run the checks again.".into(),
            );
        }
        Ok(())
    }
    pub fn validate_report(&self, report: &ivt::Report) -> Result<(), String> {
        report
            .validate_binding(&self.revision)
            .map_err(|e| e.to_string())?;
        if report.source_revision != self.binding.revision || report.request != self.request {
            return Err("The IVT response does not match its captured source and settings.".into());
        }
        Ok(())
    }
}
