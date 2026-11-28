//! Experimental screening hypotheses against explicit named products/mixtures.
//! This owner predicts ideal size signatures, never gel migration or success.
mod assessment;
mod model;
pub use assessment::{assess, revision};
pub use model::*;
#[cfg(test)]
mod tests;
