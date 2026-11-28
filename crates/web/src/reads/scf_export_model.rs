//! Pure scope/acknowledgment guard for a prepared chromatogram export.
use collagene_core::{formats::scf::SamplePrecision, review::reads::traces::TraceOrientation};
/// Select a portable receipt basename without changing the SCF sample NAME.
pub fn safe_filename(name: &str, source_hash: &str, orientation: TraceOrientation) -> String {
    let candidate = format!("{name}.scf");
    let device = name
        .split('.')
        .next()
        .unwrap_or_default()
        .trim_end()
        .to_ascii_uppercase();
    let reserved = matches!(device.as_str(), "CON" | "PRN" | "AUX" | "NUL")
        || (device.len() == 4
            && (device.starts_with("COM") || device.starts_with("LPT"))
            && matches!(device.as_bytes()[3], b'1'..=b'9'));
    if !name.is_empty()
        && name.trim() == name
        && !name.ends_with('.')
        && candidate.len() <= 180
        && !reserved
        && name
            .as_bytes()
            .first()
            .is_some_and(u8::is_ascii_alphanumeric)
        && name
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"._- ".contains(&b))
    {
        return candidate;
    }
    let hash = if source_hash.len() == 64 && source_hash.bytes().all(|b| b.is_ascii_hexdigit()) {
        source_hash[..16].to_ascii_lowercase()
    } else {
        blake3::hash(source_hash.as_bytes()).to_hex()[..16].to_string()
    };
    let orientation = match orientation {
        TraceOrientation::AsStored => "forward",
        TraceOrientation::ReverseComplement => "reverse",
    };
    format!("trace-{hash}-{orientation}.scf")
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Scope {
    /// Opaque host workspace/document/event scope plus the molecule revision.
    pub document_revision: String,
    pub source_hash: String,
    pub orientation: TraceOrientation,
    pub precision: SamplePrecision,
}
#[derive(Clone, Debug)]
pub struct Review {
    scope: Scope,
    context: String,
    has_losses: bool,
}
impl Review {
    pub fn new(scope: Scope, context: String, has_losses: bool) -> Self {
        Self {
            scope,
            context,
            has_losses,
        }
    }
    pub fn context(&self) -> &str {
        &self.context
    }
    pub fn permit(&self, current: &Scope, acknowledged: bool) -> Result<(), String> {
        if current != &self.scope {
            return Err("Trace document, orientation or export precision changed. Prepare and review again.".into());
        }
        if self.has_losses && !acknowledged {
            return Err(
                "Review and acknowledge the listed documentary omissions before saving.".into(),
            );
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn portable_trace_filename_retains_extension_without_rewriting_sample_name() {
        let hash = "ab".repeat(32);
        assert_eq!(
            safe_filename("Sample 17", &hash, TraceOrientation::AsStored),
            "Sample 17.scf"
        );
        for name in [
            "x".repeat(1024),
            "sample/../other:trace".into(),
            "Séquencé".into(),
            "CON".into(),
        ] {
            let original = name.clone();
            let filename = safe_filename(&name, &hash, TraceOrientation::ReverseComplement);
            assert_eq!(filename, "trace-abababababababab-reverse.scf");
            assert!(filename.len() <= 180 && filename.ends_with(".scf"));
            assert_eq!(name, original);
        }
    }
    fn scope() -> Scope {
        Scope {
            document_revision: "11".repeat(32),
            source_hash: "22".repeat(32),
            orientation: TraceOrientation::AsStored,
            precision: SamplePrecision::Smallest,
        }
    }
    #[test]
    fn loss_review_needs_separate_explicit_acknowledgment() {
        let s = scope();
        let review = Review::new(s.clone(), "33".repeat(32), true);
        assert!(review.permit(&s, false).is_err());
        assert!(review.permit(&s, true).is_ok());
        let complete = Review::new(s.clone(), "44".repeat(32), false);
        assert!(complete.permit(&s, false).is_ok());
    }
    #[test]
    fn document_change_invalidates_even_when_attachment_is_identical() {
        let s = scope();
        let r = Review::new(s.clone(), "33".repeat(32), false);
        let mut changed = s;
        changed.document_revision = "55".repeat(32);
        assert!(r.permit(&changed, true).is_err());
    }
    #[test]
    fn identical_trace_in_another_document_or_event_does_not_reuse_review() {
        let mut first = scope();
        let revision = first.document_revision.clone();
        first.document_revision = format!("workspace7/documentA/event1:{revision}");
        let review = Review::new(first.clone(), "33".repeat(32), true);
        let mut another = first.clone();
        another.document_revision = format!("workspace7/documentB/event1:{revision}");
        assert!(review.permit(&another, true).is_err());
        another.document_revision = format!("workspace7/documentA/event2:{revision}");
        assert!(review.permit(&another, true).is_err());
        another.document_revision = format!("workspace8/documentA/event1:{revision}");
        assert!(review.permit(&another, true).is_err());
        assert!(review.permit(&first, true).is_ok());
    }
    #[test]
    fn orientation_and_precision_change_require_new_preparation() {
        let s = scope();
        let r = Review::new(s.clone(), "33".repeat(32), true);
        let mut changed = s.clone();
        changed.orientation = TraceOrientation::ReverseComplement;
        assert!(r.permit(&changed, true).is_err());
        changed = s;
        changed.precision = SamplePrecision::SixteenBit;
        assert!(r.permit(&changed, true).is_err());
    }
}
