pub mod analysis;
pub mod assembly;
pub mod catalog;
pub mod cloning;
pub mod genome;
pub mod protocol;
pub mod sequence;
pub mod workspace;

pub(crate) fn default_inspection_page() -> usize {
    20
}
