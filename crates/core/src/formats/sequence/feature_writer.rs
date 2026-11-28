//! Independently bounded feature-table output. No upstream feature/qualifier
//! wrapper receives native strings, field widths, Unicode or escaped quotes.
use super::{
    MAX_IMPORT_BYTES, MAX_LOCATION_WORK, MAX_QUALIFIER_LINES, feature_kind_valid,
    feature_qualifiers, location_text,
};
use collagene_model::{Error, Feature, Molecule};

const WIDTH: usize = 79;
const INDENT: &str = "                     ";
const NO_QUOTE: &[&str] = &[
    "anticodon",
    "codon_start",
    "estimated_length",
    "number",
    "rpt_type",
    "transl_except",
    "transl_table",
];
const REVIEW_KEY: &str = "collagene_review_required";
const REVIEW_NOTE: &str = "Collagene: annotation requires review";

fn unsupported(message: &str) -> Error {
    Error::Unsupported(format!(
        "GenBank feature output: {message}; retain the native workspace"
    ))
}

pub(super) fn write(molecule: &Molecule) -> Result<String, Error> {
    if molecule.features.is_empty() {
        return Ok(String::new());
    }
    let mut writer = Writer {
        output: String::new(),
        length: molecule.sequence.len(),
        feature_qualifiers: 0,
        qualifier_lines: 0,
        location_work: 0,
    };
    writer.line("FEATURES             Location/Qualifiers")?;
    for feature in &molecule.features {
        writer.feature(feature)?;
    }
    Ok(writer.output)
}

struct Writer {
    output: String,
    length: usize,
    feature_qualifiers: usize,
    qualifier_lines: usize,
    location_work: usize,
}

impl Writer {
    fn line(&mut self, line: &str) -> Result<(), Error> {
        if line.len() > WIDTH || line.contains(['\r', '\n']) {
            return Err(unsupported(
                "a physical output line exceeds its safe profile",
            ));
        }
        if self
            .output
            .len()
            .checked_add(line.len())
            .and_then(|n| n.checked_add(1))
            .filter(|n| *n <= MAX_IMPORT_BYTES)
            .is_none()
        {
            return Err(Error::TooLarge);
        }
        self.output.push_str(line);
        self.output.push('\n');
        Ok(())
    }

    fn qualifier_line(&mut self, line: &str) -> Result<(), Error> {
        self.qualifier_lines = self
            .qualifier_lines
            .checked_add(1)
            .filter(|n| *n <= MAX_QUALIFIER_LINES)
            .ok_or(Error::TooLarge)?;
        self.line(line)
    }

    fn feature(&mut self, feature: &Feature) -> Result<(), Error> {
        self.feature_qualifiers = 0;
        if !feature_kind_valid(&feature.kind) || feature.kind.len() > 15 {
            return Err(unsupported("feature key does not fit its 15-column field"));
        }
        let location = location_text(&feature.location)?;
        super::parse_location(&location, self.length)?;
        if !location.is_ascii() || location.bytes().any(|b| b.is_ascii_control()) {
            return Err(unsupported(
                "a location must be a single-line ASCII expression",
            ));
        }
        let work = 1 + location
            .bytes()
            .filter(|b| matches!(b, b',' | b'(' | b':'))
            .count();
        self.location_work = self
            .location_work
            .checked_add(work)
            .filter(|n| *n <= MAX_LOCATION_WORK)
            .ok_or(Error::TooLarge)?;
        self.ascii(&format!("     {:<16}", feature.kind), &location, false)?;
        for (key, value) in &feature.qualifiers {
            if key != REVIEW_KEY {
                self.qualifier(key, value.as_deref())?;
            }
        }
        if feature.needs_review {
            self.qualifier(REVIEW_KEY, Some("true"))?;
            if !feature
                .qualifiers
                .iter()
                .any(|(key, value)| key == "note" && value.as_deref() == Some(REVIEW_NOTE))
            {
                self.qualifier("note", Some(REVIEW_NOTE))?;
            }
        }
        Ok(())
    }

    /// Locations and position qualifiers are ASCII; prefer a comma boundary
    /// without ever subtracting an unchecked prefix length or dropping input.
    fn ascii(&mut self, prefix: &str, mut text: &str, qualifier: bool) -> Result<(), Error> {
        let mut line = prefix.to_owned();
        loop {
            let available = WIDTH
                .checked_sub(line.len())
                .filter(|n| *n > 0)
                .ok_or_else(|| unsupported("field prefix leaves no output capacity"))?;
            let count = if text.len() <= available {
                text.len()
            } else {
                text[..available].rfind(',').map_or(available, |i| i + 1)
            };
            line.push_str(&text[..count]);
            if qualifier {
                self.qualifier_line(&line)?;
            } else {
                self.line(&line)?;
            }
            text = &text[count..];
            if text.is_empty() {
                return Ok(());
            }
            if qualifier && text.starts_with('/') {
                return Err(unsupported(
                    "a bare continuation would introduce another qualifier",
                ));
            }
            line = INDENT.to_owned();
        }
    }

    fn qualifier(&mut self, key: &str, value: Option<&str>) -> Result<(), Error> {
        self.feature_qualifiers = self
            .feature_qualifiers
            .checked_add(1)
            .filter(|n| *n <= feature_qualifiers::MAX_QUALIFIERS)
            .ok_or(Error::TooLarge)?;
        feature_qualifiers::validate_key(key)?;
        if key.len() > 54 {
            return Err(unsupported("qualifier key exceeds the output line width"));
        }
        let Some(value) = value else {
            return self.qualifier_line(&format!("{INDENT}/{key}"));
        };
        feature_qualifiers::validate_value(value)?;
        let quoted = !NO_QUOTE.contains(&key);
        if !quoted && value.contains(['\n', '"']) {
            return Err(unsupported(
                "unquoted qualifiers cannot contain line breaks or quotes",
            ));
        }
        if matches!(key, "anticodon" | "transl_except") {
            if !value.is_ascii() || value.bytes().any(|b| b.is_ascii_control()) {
                return Err(unsupported(
                    "position qualifiers require single-line ASCII descriptors",
                ));
            }
            return self.ascii(&format!("{INDENT}/{key}="), value, true);
        }
        let mut line = format!("{INDENT}/{key}={}", if quoted { "\"" } else { "" });
        let mut chars = value.chars().peekable();
        while let Some(character) = chars.next() {
            if character == '\n' {
                self.qualifier_line(&line)?;
                line = INDENT.to_owned();
                continue;
            }
            let mut buffer = [0u8; 4];
            let token = if quoted && character == '"' {
                "\"\""
            } else {
                character.encode_utf8(&mut buffer)
            };
            let reserved = usize::from(quoted && chars.peek().is_none());
            if line
                .len()
                .checked_add(token.len())
                .and_then(|n| n.checked_add(reserved))
                .filter(|n| *n <= WIDTH)
                .is_none()
            {
                self.qualifier_line(&line)?;
                line = INDENT.to_owned();
                if !quoted && token == "/" {
                    return Err(unsupported(
                        "a bare continuation would introduce another qualifier",
                    ));
                }
            }
            line.push_str(token);
        }
        if quoted {
            line.push('"');
        }
        self.qualifier_line(&line)
    }
}

#[cfg(test)]
mod tests {
    use super::{INDENT, MAX_LOCATION_WORK, MAX_QUALIFIER_LINES, WIDTH, Writer, write};
    use crate::formats::sequence::feature_qualifiers;
    use collagene_model::{Feature, Location, Molecule, Topology};

    fn molecule(qualifiers: Vec<(String, Option<String>)>, review: bool) -> Molecule {
        let mut molecule = Molecule::dna("literal", "ACGTACGT", Topology::Linear).unwrap();
        molecule.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::Range {
                start: 0,
                end: 8,
                before: false,
                after: false,
            },
            qualifiers,
            needs_review: review,
        });
        molecule
    }

    fn qualifiers(table: &str) -> Vec<(String, Option<String>)> {
        let raw = table
            .lines()
            .skip(2)
            .map(|line| line[INDENT.len()..].to_owned())
            .collect::<Vec<_>>();
        feature_qualifiers::parse(&raw).unwrap()
    }

    #[test]
    fn unicode_and_quote_at_old_wrapper_boundary_keep_every_character_without_panicking() {
        let value = format!("{}\"{}", "a".repeat(50), "β".repeat(40));
        let table = write(&molecule(vec![("note".into(), Some(value.clone()))], false)).unwrap();
        assert!(table.lines().all(|line| line.len() <= WIDTH));
        let parsed = qualifiers(&table);
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].1.as_ref().unwrap().replace('\n', ""), value);
    }

    #[test]
    fn flags_repeats_explicit_prose_lines_and_review_marker_have_one_current_meaning() {
        let table = write(&molecule(
            vec![
                ("pseudo".into(), None),
                ("note".into(), Some("first\n/second is quoted prose".into())),
                ("note".into(), Some("second".into())),
                ("collagene_review_required".into(), Some("false".into())),
            ],
            true,
        ))
        .unwrap();
        let parsed = qualifiers(&table);
        assert_eq!(parsed[0], ("pseudo".into(), None));
        assert_eq!(
            parsed[1].1.as_deref(),
            Some("first\n/second is quoted prose")
        );
        assert_eq!(parsed.iter().filter(|(k, _)| k == "note").count(), 3);
        assert_eq!(
            parsed
                .iter()
                .filter(|(k, _)| k == "collagene_review_required")
                .count(),
            1
        );
        assert_eq!(
            parsed
                .iter()
                .find(|(k, _)| k == "collagene_review_required")
                .unwrap()
                .1
                .as_deref(),
            Some("true")
        );
    }

    #[test]
    fn long_prefix_and_bare_line_injection_are_refused_without_calling_upstream_writer() {
        assert!(
            write(&molecule(
                vec![("x".repeat(55), Some("value".into()))],
                false
            ))
            .is_err()
        );
        assert!(
            write(&molecule(
                vec![("codon_start".into(), Some("1\n/note=hidden".into()))],
                false
            ))
            .is_err()
        );
        let value = "9".repeat(200);
        let table = write(&molecule(
            vec![("number".into(), Some(value.clone()))],
            false,
        ))
        .unwrap();
        assert!(table.lines().all(|line| line.len() <= WIDTH));
        assert_eq!(qualifiers(&table)[0].1.as_deref(), Some(value.as_str()));
    }

    #[test]
    fn emitted_qualifier_bound_counts_review_additions_and_ignores_removed_markers() {
        let flags = vec![("pseudo".into(), None); feature_qualifiers::MAX_QUALIFIERS];
        assert!(write(&molecule(flags.clone(), false)).is_ok());
        assert!(write(&molecule(flags.clone(), true)).is_err());
        let mut with_removed_marker = flags;
        with_removed_marker.push(("collagene_review_required".into(), Some("false".into())));
        assert!(write(&molecule(with_removed_marker, false)).is_ok());
        assert!(
            write(&molecule(
                vec![("pseudo".into(), None); feature_qualifiers::MAX_QUALIFIERS + 1],
                false
            ))
            .is_err()
        );
    }

    #[test]
    fn continuation_and_location_work_bounds_refuse_before_the_overflowing_output_line() {
        let mut writer = Writer {
            output: String::new(),
            length: 8,
            feature_qualifiers: 0,
            qualifier_lines: MAX_QUALIFIER_LINES - 1,
            location_work: MAX_LOCATION_WORK - 2,
        };
        assert!(writer.qualifier("note", Some("first\nsecond")).is_err());
        assert_eq!(writer.output.lines().count(), 1);
        let before = writer.output.clone();
        let mut feature = molecule(vec![], false).features.remove(0);
        feature.location = Location::Join(vec![
            Location::Range {
                start: 0,
                end: 4,
                before: false,
                after: false,
            },
            Location::Range {
                start: 4,
                end: 8,
                before: false,
                after: false,
            },
        ]);
        assert!(writer.feature(&feature).is_err());
        assert_eq!(writer.output, before);
    }
}
