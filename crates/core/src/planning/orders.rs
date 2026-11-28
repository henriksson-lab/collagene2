//! Literal named oligos converted into the existing bounded order/collection format.
//! No binding, synthesis, physical material or experimental success is inferred.
use crate::{
    formats::oligos::{
        self as oligo_formats, ColumnMapping, Delimiter, OligoImportOptions, OligoTable,
    },
    sequence::analysis::Primer,
};
use collagene_model::{Alphabet, Error};
use std::io::{self, Write};

const TABLE_CELL_BYTES: usize = 8192; // Existing oligo_formats cell profile.
struct Limited(Vec<u8>);
impl Write for Limited {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.0.len().saturating_add(bytes.len()) > oligo_formats::MAX_TABLE_BYTES {
            return Err(io::Error::other("oligo table byte limit"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn csv_error(error: csv::Error) -> Error {
    if error.is_io_error() {
        Error::TooLarge
    } else {
        Error::Format(error.to_string())
    }
}
/// Preserve row order, names, complete ordered DNA and explicit domain/tail split.
/// Single oligos are single rows, duplicate names remain distinct rows, and long
/// oligos follow the existing table's evidence-only reusable-primer policy.
pub fn from_primers(primers: &[Primer]) -> Result<OligoTable, Error> {
    if primers.is_empty() || primers.len() > oligo_formats::MAX_TABLE_ROWS {
        return Err(Error::Format(
            "Select 1..10000 named oligos for the order table".into(),
        ));
    }
    for primer in primers {
        if primer.name.trim().is_empty()
            || primer.name.len() > 256
            || primer.name.chars().any(char::is_control)
        {
            return Err(Error::Format(
                "Oligo names must be nonempty, at most 256 bytes and free of control characters"
                    .into(),
            ));
        }
        if primer.annealing.alphabet() != Alphabet::Dna || primer.annealing.is_empty() {
            return Err(Error::Format(
                "Oligo annealing domains must contain DNA".into(),
            ));
        }
        if primer.annealing.len().saturating_add(primer.tail.len()) > TABLE_CELL_BYTES {
            return Err(Error::TooLarge);
        }
        // Revalidate public serialized fields with the same DNA normalization as
        // the table adapter. IUPAC DNA remains explicit degeneracy, not chemistry.
        Primer::new(&primer.name, primer.annealing.as_str(), &primer.tail)?;
    }
    let mut bytes = Limited(Vec::new());
    {
        let mut writer = csv::WriterBuilder::new()
            .quote_style(csv::QuoteStyle::Always)
            .from_writer(&mut bytes);
        writer
            .write_record([
                "Name",
                "Ordered sequence",
                "Annealing domain",
                "5prime tail",
            ])
            .map_err(csv_error)?;
        for primer in primers {
            writer
                .write_record([
                    primer.name.as_str(),
                    primer.ordered_sequence().as_str(),
                    primer.annealing.as_str(),
                    primer.tail.as_str(),
                ])
                .map_err(csv_error)?;
        }
        writer.flush().map_err(|_| Error::TooLarge)?;
    }
    let source = String::from_utf8(bytes.0).map_err(|e| Error::Format(e.to_string()))?;
    oligo_formats::import_oligos(
        &source,
        &OligoImportOptions {
            delimiter: Delimiter::Comma,
            has_header: true,
            columns: ColumnMapping {
                name: 0,
                ordered_sequence: Some(1),
                annealing_sequence: Some(2),
                tail: Some(3),
                modification: None,
                amount: None,
                unit: None,
            },
        },
    )
}
