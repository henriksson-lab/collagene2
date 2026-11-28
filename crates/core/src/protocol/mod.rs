//! Versioned bounded binary messages and validated document serialization.
pub mod bindings;
pub mod envelope;
pub mod execution;
pub mod request;
pub mod response;
pub mod schema;
pub mod session;
mod validation;

use self::envelope::{decode, encode};
use collagene_model::{Error, Molecule};
pub fn save_document(mol: &Molecule) -> Result<Vec<u8>, Error> {
    mol.validate()?;
    encode(b"CG2\0", mol)
}
pub fn load_document(bytes: &[u8]) -> Result<Molecule, Error> {
    let mol: Molecule = decode(b"CG2\0", bytes)?;
    mol.validate()?;
    Ok(mol)
}
