//! Exact UI ownership of one asynchronous chromatogram preparation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Capture {
    pub scope: String,
    pub revision: String,
    pub generation: u32,
}
impl Capture {
    pub fn guard(&self, scope: &str, revision: &str, generation: u32) -> Result<(), String> {
        if self.scope != scope || self.revision != revision || self.generation != generation {
            Err("The trace document/event or preparation changed.".into())
        } else {
            Ok(())
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn equal_sequence_documents_foreign_events_and_late_generations_do_not_retarget() {
        let capture = Capture {
            scope: "epoch7:documentA:event2".into(),
            revision: "a".repeat(64),
            generation: 9,
        };
        capture.guard(&capture.scope, &capture.revision, 9).unwrap();
        for scope in [
            "epoch7:documentB:event2",
            "epoch8:documentA:event2",
            "epoch7:documentA:event3",
        ] {
            assert!(capture.guard(scope, &capture.revision, 9).is_err());
        }
        assert!(capture.guard(&capture.scope, &"b".repeat(64), 9).is_err());
        assert!(
            capture
                .guard(&capture.scope, &capture.revision, 10)
                .is_err()
        );
    }
}
