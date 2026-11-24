//! Compact immutable column-to-residue index. No full-alignment scan during painting.
use collagene_core::{Error, alignment_formats::AlignmentDocument};
const BLOCK: usize = 256;
#[derive(Debug)]
pub struct PreparedAlignment {
    document: AlignmentDocument,
    prefixes: Vec<Vec<u32>>,
}
impl PreparedAlignment {
    pub fn new(document: AlignmentDocument) -> Result<Self, Error> {
        document.validate()?;
        let prefixes = document
            .rows
            .iter()
            .map(|row| {
                let mut count = 0u32;
                let mut p = vec![0];
                for block in row.aligned.as_bytes().chunks(BLOCK) {
                    count += block.iter().filter(|b| **b != b'-').count() as u32;
                    p.push(count);
                }
                p
            })
            .collect();
        Ok(Self { document, prefixes })
    }
    pub fn document(&self) -> &AlignmentDocument {
        &self.document
    }
    /// Zero-based index in the ungapped row as stored; a gap has no residue index.
    pub fn position(&self, row: usize, column: usize) -> Result<Option<usize>, Error> {
        let row_data = self
            .document
            .rows
            .get(row)
            .ok_or_else(|| Error::Coordinates("unknown alignment row".into()))?;
        let bases = row_data.aligned.as_bytes();
        if column >= bases.len() {
            return Err(Error::Coordinates("alignment column outside row".into()));
        }
        if bases[column] == b'-' {
            return Ok(None);
        }
        let block = column / BLOCK;
        Ok(Some(
            self.prefixes[row][block] as usize
                + bases[block * BLOCK..column]
                    .iter()
                    .filter(|b| **b != b'-')
                    .count(),
        ))
    }
    pub fn residue_count(&self, row: usize) -> Option<usize> {
        self.prefixes.get(row)?.last().map(|n| *n as usize)
    }
}

/// Keep native CSS scroll dimensions below browser implementation limits while
/// retaining the complete logical column range (up to eight million columns).
pub const MAX_SCROLL_WIDTH: f64 = 8_000_000.;
pub fn spacer_width(columns: usize) -> f64 {
    (180. + columns as f64 * 10.).min(MAX_SCROLL_WIDTH)
}
pub fn logical_left(native: f64, columns: usize, viewport: f64) -> f64 {
    let native_range = (spacer_width(columns) - viewport).max(0.);
    let logical_range = (180. + columns as f64 * 10. - viewport).max(0.);
    if native_range == 0. {
        0.
    } else {
        native.clamp(0., native_range) * logical_range / native_range
    }
}
pub fn native_left(logical: f64, columns: usize, viewport: f64) -> f64 {
    let native_range = (spacer_width(columns) - viewport).max(0.);
    let logical_range = (180. + columns as f64 * 10. - viewport).max(0.);
    if logical_range == 0. {
        0.
    } else {
        logical.clamp(0., logical_range) * native_range / logical_range
    }
}
