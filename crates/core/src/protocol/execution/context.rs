//! Borrowed current-request authority for domain admission. Never serialized or imported.
use collagene_model::{Error, Molecule};

pub(super) struct RequestContext<'a> {
    pub(super) revision: &'a String,
}
impl RequestContext<'_> {
    pub(super) fn admit(&self, molecule: Option<&Molecule>) -> Result<(), Error> {
        if let Some(molecule) = molecule {
            molecule.validate()?;
            if molecule.revision() != *self.revision {
                return Err(Error::Format("worker input revision mismatch".into()));
            }
        }
        Ok(())
    }
}
