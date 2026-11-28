//! Bounded complete INSDC qualifier statements, independent of the header reader.
//! Lines contain only the content after feature-table indentation. Literal
//! backslashes remain literal; INSDC quote escaping uses doubled double quotes.
use collagene_model::Error;

pub(super) const MAX_QUALIFIERS: usize = 16_384;
const MAX_KEY_BYTES: usize = 4096;
const MAX_VALUE_BYTES: usize = 1024 * 1024;

fn invalid(message: &str) -> Error {
    Error::Format(format!("GenBank feature qualifiers: {message}"))
}

/// The same field grammar applies to raw imports and native feature export.
pub(super) fn validate_key(key: &str) -> Result<(), Error> {
    if key.len() > MAX_KEY_BYTES {
        return Err(Error::TooLarge);
    }
    // Digits retain positional names such as the saved-primer `5prime_tail`.
    // ':' '.' '-' retain the application's explicit custom namespaces;
    // whitespace, quotes, slashes and structural delimiters are refused.
    if !key
        .as_bytes()
        .first()
        .is_some_and(|b| b.is_ascii_alphanumeric() || *b == b'_')
        || !key
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b':' | b'.' | b'-'))
    {
        return Err(invalid("invalid qualifier field key"));
    }
    Ok(())
}

/// Native values may contain Unicode and ordinary prose wrapping, but cannot
/// insert carriage returns or other control characters into a physical record.
pub(super) fn validate_value(value: &str) -> Result<(), Error> {
    if value.len() > MAX_VALUE_BYTES {
        return Err(Error::TooLarge);
    }
    if value
        .chars()
        .any(|c| c.is_control() && !matches!(c, '\n' | '\t'))
    {
        return Err(invalid(
            "qualifier value contains an unsupported control character",
        ));
    }
    Ok(())
}

/// Parse every physical content line, retaining qualifier order and repeated
/// keys. No prefix parse or unconsumed nonblank line can produce success.
pub(super) fn parse(lines: &[String]) -> Result<Vec<(String, Option<String>)>, Error> {
    let mut bytes = 0usize;
    for line in lines {
        bytes = bytes
            .checked_add(line.len())
            .and_then(|n| n.checked_add(1))
            .filter(|n| *n <= super::MAX_IMPORT_BYTES)
            .ok_or(Error::TooLarge)?;
        if line.bytes().any(|b| b < b' ' && b != b'\t' || b == 0x7f) {
            return Err(invalid(
                "physical content contains a newline or control byte",
            ));
        }
    }
    let mut text = String::with_capacity(bytes);
    for line in lines {
        text.push_str(line);
        text.push('\n');
    }
    Parser {
        text: &text,
        cursor: 0,
    }
    .finish()
}

struct Parser<'a> {
    text: &'a str,
    cursor: usize,
}

impl Parser<'_> {
    fn line_end(&self) -> usize {
        self.text[self.cursor..]
            .find('\n')
            .map_or(self.text.len(), |offset| self.cursor + offset)
    }

    fn finish(mut self) -> Result<Vec<(String, Option<String>)>, Error> {
        let mut qualifiers = Vec::new();
        while self.cursor < self.text.len() {
            let end = self.line_end();
            if self.text[self.cursor..end]
                .trim_matches([' ', '\t'])
                .is_empty()
            {
                self.cursor = (end + 1).min(self.text.len());
                continue;
            }
            if qualifiers.len() == MAX_QUALIFIERS {
                return Err(Error::TooLarge);
            }
            if self.text.as_bytes()[self.cursor] != b'/' {
                return Err(invalid(
                    "a qualifier must start with '/' at the content column",
                ));
            }
            self.cursor += 1;
            let statement = &self.text[self.cursor..end];
            let (raw_key, has_value) = statement
                .split_once('=')
                .map_or((statement, false), |(key, _)| (key, true));
            let key = if has_value {
                raw_key
            } else {
                raw_key.trim_end_matches([' ', '\t'])
            };
            validate_key(key)?;
            let key = key.to_owned();
            let value = if has_value {
                self.cursor += raw_key.len() + 1;
                Some(if self.text.as_bytes().get(self.cursor) == Some(&b'"') {
                    self.quoted(&key)?
                } else {
                    self.bare()?
                })
            } else {
                self.cursor = (end + 1).min(self.text.len());
                None
            };
            qualifiers.push((key, value));
        }
        Ok(qualifiers)
    }

    fn append(value: &mut String, part: &str) -> Result<(), Error> {
        if value
            .len()
            .checked_add(part.len())
            .filter(|n| *n <= MAX_VALUE_BYTES)
            .is_none()
        {
            return Err(Error::TooLarge);
        }
        value.push_str(part);
        Ok(())
    }

    fn bare(&mut self) -> Result<String, Error> {
        let mut value = String::new();
        loop {
            let end = self.line_end();
            let part = &self.text[self.cursor..end];
            if part.contains('"') {
                return Err(invalid("an unquoted value contains a quote"));
            }
            Self::append(&mut value, part)?;
            self.cursor = (end + 1).min(self.text.len());
            if self.cursor == self.text.len() || self.text.as_bytes()[self.cursor] == b'/' {
                return Ok(value);
            }
            // Bare sequence/citation values concatenate physical continuations,
            // matching INSDC wrapping without inventing a separator.
        }
    }

    fn quoted(&mut self, key: &str) -> Result<String, Error> {
        self.cursor += 1;
        let mut value = String::new();
        let bytes = self.text.as_bytes();
        let mut start = self.cursor;
        while self.cursor < bytes.len() {
            match bytes[self.cursor] {
                b'"' => {
                    Self::append(&mut value, &self.text[start..self.cursor])?;
                    let escaped = if bytes.get(self.cursor + 1) == Some(&b'"') {
                        Some(2)
                    } else if bytes.get(self.cursor + 1) == Some(&b'\n')
                        && bytes.get(self.cursor + 2) == Some(&b'"')
                    {
                        Some(3)
                    } else {
                        None
                    };
                    if let Some(consumed) = escaped {
                        Self::append(&mut value, "\"")?;
                        self.cursor += consumed;
                        start = self.cursor;
                        continue;
                    }
                    self.cursor += 1;
                    let end = self.line_end();
                    if !self.text[self.cursor..end]
                        .trim_matches([' ', '\t'])
                        .is_empty()
                    {
                        return Err(invalid("trailing content after a closing quote"));
                    }
                    self.cursor = (end + 1).min(self.text.len());
                    return Ok(value);
                }
                b'\n' => {
                    Self::append(&mut value, &self.text[start..self.cursor])?;
                    // These machine fields contain no literal newline semantics.
                    // JSON string newlines are escaped; physical INSDC wrapping
                    // must not insert a control character into the typed ledger.
                    if key != "translation" && key != crate::sequence::annotations::merge::QUALIFIER
                    {
                        Self::append(&mut value, "\n")?;
                    }
                    self.cursor += 1;
                    start = self.cursor;
                }
                _ => self.cursor += 1,
            }
        }
        Err(invalid("unterminated quoted value"))
    }
}

#[cfg(test)]
mod tests {
    use super::{MAX_QUALIFIERS, MAX_VALUE_BYTES, parse, validate_key};

    fn lines(values: &[&str]) -> Vec<String> {
        values.iter().map(|line| (*line).to_owned()).collect()
    }

    #[test]
    fn repeated_flags_values_and_multiline_sequences_keep_their_declared_meaning() {
        let parsed = parse(&lines(&[
            "/pseudo",
            "/note=\"first line",
            "second line\"",
            "/translation=\"MKWV",
            "TFIS\"",
            "/citation=[1]",
            "/citation=[2]",
            "/empty=",
            "/label=\"\"",
        ]))
        .unwrap();
        assert_eq!(
            parsed,
            vec![
                ("pseudo".into(), None),
                ("note".into(), Some("first line\nsecond line".into())),
                ("translation".into(), Some("MKWVTFIS".into())),
                ("citation".into(), Some("[1]".into())),
                ("citation".into(), Some("[2]".into())),
                ("empty".into(), Some("".into())),
                ("label".into(), Some("".into())),
            ]
        );
    }

    #[test]
    fn positional_primer_keys_are_safe_custom_fields() {
        let parsed = parse(&lines(&["/5prime_tail=\"GG\"", "/3prime_tail=\"AA\""])).unwrap();
        assert_eq!(
            parsed,
            vec![
                ("5prime_tail".into(), Some("GG".into())),
                ("3prime_tail".into(), Some("AA".into()))
            ]
        );
        for key in ["5prime tail", "5prime/tail", "5prime=tail", "5prime\ntail"] {
            assert!(validate_key(key).is_err(), "{key:?}");
        }
    }

    #[test]
    fn insdc_escapes_unicode_and_literal_backslashes_are_retained_without_new_statements() {
        let parsed = parse(&lines(&[
            "/note=\"β says \"\"hello\"\" C:\\data",
            "/gene=literal quoted content\"",
            "/note=\"split \"",
            "\"quote\"",
            "/rpt_unit_seq=ACGT",
            "TGCA",
        ]))
        .unwrap();
        assert_eq!(
            parsed[0].1.as_deref(),
            Some("β says \"hello\" C:\\data\n/gene=literal quoted content")
        );
        assert_eq!(parsed[1].1.as_deref(), Some("split \"quote"));
        assert_eq!(parsed[2].1.as_deref(), Some("ACGTTGCA"));
    }

    #[test]
    fn malformed_qualifiers_never_succeed_as_prefixes_or_drop_following_statements() {
        for bad in [
            vec!["/note=\"unterminated", "/gene=x"],
            vec!["/note=\"valid\" junk", "/gene=x"],
            vec!["/note=bare\"quote"],
            vec!["/note=\"valid\" /gene=x"],
            vec!["/pseudo", "orphan text"],
            vec!["/bad key=x"],
            vec!["/=x"],
            vec!["//note=x"],
            vec!["/note=\"ok\"", "-9223372036854775808"],
        ] {
            assert!(parse(&lines(&bad)).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn qualifier_counts_values_and_physical_lines_are_bounded() {
        assert!(parse(&vec!["/pseudo".into(); MAX_QUALIFIERS + 1]).is_err());
        assert!(parse(&[format!("/note=\"{}\"", "a".repeat(MAX_VALUE_BYTES + 1))]).is_err());
        assert!(parse(&["/note=bad\n/gene=hidden".into()]).is_err());
        assert!(parse(&["/note=bad\r/gene=hidden".into()]).is_err());
    }
}
