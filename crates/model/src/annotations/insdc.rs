//! Bounded, full-consumption INSDC locations. Unsupported semantics retain the
//! complete validated expression instead of being converted into precise bases.
use crate::{Error, Location, MAX_SEQUENCE_LEN};

const MAX_BYTES: usize = 65_536;
const MAX_DEPTH: usize = 32;
const MAX_NODES: usize = 4_096;

fn invalid(reason: &str) -> Error {
    Error::Format(format!("invalid INSDC location: {reason}"))
}

pub fn check_text(text: &str) -> Result<(), Error> {
    if text.is_empty()
        || text.len() > MAX_BYTES
        || !text.is_ascii()
        || text.bytes().any(|byte| byte.is_ascii_control())
    {
        return Err(invalid("empty, oversized or non-ASCII/control text"));
    }
    let mut depth = 0usize;
    for byte in text.bytes() {
        match byte {
            b'(' => {
                depth += 1;
                if depth >= MAX_DEPTH {
                    return Err(invalid("nesting exceeds 32 levels"));
                }
            }
            b')' => {
                depth = depth
                    .checked_sub(1)
                    .ok_or_else(|| invalid("unbalanced parentheses"))?;
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err(invalid("unbalanced parentheses"));
    }
    Ok(())
}

pub fn parse_location(text: &str, length: usize) -> Result<Location, Error> {
    let text = text.trim();
    let location = Parser::new(text, Some(length))?.complete()?;
    let location = location.unwrap_or_else(|| Location::Unresolved(text.to_owned()));
    location.validate(length)?;
    Ok(location)
}

/// Validate raw feature-table syntax before any general GenBank parser executes.
/// Remote positions and local positions are not compared to a record length here.
/// Positive coordinates are limited to i64::MAX for safe GenBank interoperability.
pub fn validate_for_genbank(text: &str) -> Result<(), Error> {
    Parser::new(text.trim(), None)?.complete().map(|_| ())
}

/// Render the current native location directly, with the same syntax and bounds
/// as parsing. Opaque validated descriptors retain their complete original text.
pub fn location_text(location: &Location) -> Result<String, Error> {
    let mut formatter = Formatter {
        text: String::new(),
        nodes: 0,
    };
    formatter.location(location, 0)?;
    // A raw unresolved child can contain operators of its own: validate the
    // assembled expression too, so nesting/work limits apply across boundaries.
    Parser::new(&formatter.text, None)?.complete()?;
    Ok(formatter.text)
}

struct Formatter {
    text: String,
    nodes: usize,
}
impl Formatter {
    fn append(&mut self, text: &str) -> Result<(), Error> {
        if self
            .text
            .len()
            .checked_add(text.len())
            .is_none_or(|length| length > MAX_BYTES)
        {
            return Err(Error::TooLarge);
        }
        self.text.push_str(text);
        Ok(())
    }
    fn location(&mut self, location: &Location, depth: usize) -> Result<(), Error> {
        self.nodes += 1;
        if self.nodes > MAX_NODES {
            return Err(Error::TooLarge);
        }
        if depth >= MAX_DEPTH {
            return Err(invalid("nesting exceeds 32 levels"));
        }
        match location {
            Location::Range {
                start,
                end,
                before,
                after,
            } => {
                if start >= end || *end > MAX_SEQUENCE_LEN {
                    return Err(Error::Coordinates(
                        "native location interval exceeds model bounds".into(),
                    ));
                }
                let first = start
                    .checked_add(1)
                    .ok_or_else(|| invalid("native coordinate overflow"))?;
                if !*before && !*after && first == *end {
                    return self.append(&first.to_string());
                }
                if *before {
                    self.append("<")?;
                }
                self.append(&first.to_string())?;
                self.append("..")?;
                if *after {
                    self.append(">")?;
                }
                self.append(&end.to_string())
            }
            Location::Complement(inner) => {
                self.append("complement(")?;
                self.location(inner, depth + 1)?;
                self.append(")")
            }
            Location::Join(parts) | Location::Order(parts) => {
                if parts.is_empty() {
                    return Err(invalid("empty compound location"));
                }
                self.append(if matches!(location, Location::Join(_)) {
                    "join("
                } else {
                    "order("
                })?;
                for (index, part) in parts.iter().enumerate() {
                    if index != 0 {
                        self.append(",")?;
                    }
                    self.location(part, depth + 1)?;
                }
                self.append(")")
            }
            Location::Unresolved(raw) => {
                Parser::new(raw, None)?.complete()?;
                self.append(raw)
            }
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Boundary {
    Exact,
    Before,
    After,
    Uncertain,
}
#[derive(Clone, Copy)]
struct Position {
    low: u64,
    high: u64,
    boundary: Boundary,
}

struct Parser<'a> {
    text: &'a str,
    offset: usize,
    nodes: usize,
    length: Option<usize>,
}
impl<'a> Parser<'a> {
    fn new(text: &'a str, length: Option<usize>) -> Result<Self, Error> {
        check_text(text)?;
        Ok(Self {
            text,
            offset: 0,
            nodes: 0,
            length,
        })
    }
    fn complete(mut self) -> Result<Option<Location>, Error> {
        let location = self.location(0, false)?;
        self.spaces();
        if self.offset != self.text.len() {
            return Err(invalid("unconsumed suffix or multiple top-level locations"));
        }
        Ok(location)
    }
    fn charge(&mut self, depth: usize) -> Result<(), Error> {
        self.nodes += 1;
        if self.nodes > MAX_NODES {
            return Err(Error::TooLarge);
        }
        if depth >= MAX_DEPTH {
            return Err(invalid("nesting exceeds 32 levels"));
        }
        Ok(())
    }
    fn spaces(&mut self) {
        while self.text.as_bytes().get(self.offset) == Some(&b' ') {
            self.offset += 1;
        }
    }
    fn take(&mut self, token: &str) -> bool {
        self.spaces();
        if self.text[self.offset..].starts_with(token) {
            self.offset += token.len();
            true
        } else {
            false
        }
    }
    fn expect(&mut self, token: &str) -> Result<(), Error> {
        if self.take(token) {
            Ok(())
        } else {
            Err(invalid("missing token or malformed operator"))
        }
    }
    fn number(&mut self) -> Result<u64, Error> {
        self.spaces();
        let begin = self.offset;
        let mut value = 0u64;
        while let Some(byte @ b'0'..=b'9') = self.text.as_bytes().get(self.offset).copied() {
            value = value
                .checked_mul(10)
                .and_then(|n| n.checked_add(u64::from(byte - b'0')))
                .filter(|n| *n <= i64::MAX as u64)
                .ok_or_else(|| invalid("coordinate exceeds the positive signed-64-bit bound"))?;
            self.offset += 1;
        }
        if self.offset == begin || value == 0 {
            return Err(invalid(
                "coordinates and gap lengths must be positive integers",
            ));
        }
        Ok(value)
    }
    fn position(&mut self, depth: usize) -> Result<Position, Error> {
        self.charge(depth)?;
        if self.take("one-of(") {
            let first = self.number()?;
            let (mut low, mut high) = (first, first);
            while self.take(",") {
                self.charge(depth + 1)?;
                let value = self.number()?;
                low = low.min(value);
                high = high.max(value);
            }
            self.expect(")")?;
            return Ok(Position {
                low,
                high,
                boundary: Boundary::Uncertain,
            });
        }
        if self.take("(") {
            let low = self.number()?;
            self.expect(".")?;
            let high = self.number()?;
            self.expect(")")?;
            if low > high {
                return Err(invalid("reversed uncertain position"));
            }
            return Ok(Position {
                low,
                high,
                boundary: Boundary::Uncertain,
            });
        }
        let boundary = if self.take("<") {
            Boundary::Before
        } else if self.take(">") {
            Boundary::After
        } else {
            Boundary::Exact
        };
        let value = self.number()?;
        Ok(Position {
            low: value,
            high: value,
            boundary,
        })
    }
    fn check_position(&self, position: Position, remote: bool) -> Result<(), Error> {
        if !remote
            && self
                .length
                .is_some_and(|length| position.high > length as u64)
        {
            return Err(Error::Coordinates(
                "INSDC position exceeds the local sequence length".into(),
            ));
        }
        Ok(())
    }
    fn local_range(
        &self,
        start: Position,
        end: Position,
        remote: bool,
    ) -> Result<Option<Location>, Error> {
        self.check_position(start, remote)?;
        self.check_position(end, remote)?;
        if start.low > end.high {
            return Err(invalid(
                "range end precedes its start; use an explicit circular join",
            ));
        }
        if remote
            || self.length.is_none()
            || !matches!(start.boundary, Boundary::Exact | Boundary::Before)
            || !matches!(end.boundary, Boundary::Exact | Boundary::After)
        {
            return Ok(None);
        }
        let zero_based = start
            .low
            .checked_sub(1)
            .ok_or_else(|| invalid("zero coordinate"))?;
        let start_coordinate =
            usize::try_from(zero_based).map_err(|_| invalid("local coordinate overflow"))?;
        let end_coordinate =
            usize::try_from(end.high).map_err(|_| invalid("local coordinate overflow"))?;
        Ok(Some(Location::Range {
            start: start_coordinate,
            end: end_coordinate,
            before: start.boundary == Boundary::Before,
            after: end.boundary == Boundary::After,
        }))
    }
    /// Only a lexical accession token may precede ':'. Operators, commas and
    /// parentheses cannot be swallowed into an accession and escape a delimiter.
    fn remote_prefix(&mut self) -> Result<bool, Error> {
        self.spaces();
        let begin = self.offset;
        let mut end = begin;
        while self
            .text
            .as_bytes()
            .get(end)
            .is_some_and(|b| b.is_ascii_alphanumeric() || matches!(*b, b'_' | b'.'))
        {
            end += 1;
        }
        if self.text.as_bytes().get(end) != Some(&b':') {
            return Ok(false);
        }
        let token = &self.text[begin..end];
        if token.len() > 128 {
            return Err(invalid("oversized remote accession"));
        }
        let (accession, version) = token
            .split_once('.')
            .map_or((token, None), |(a, b)| (a, Some(b)));
        if !accession
            .as_bytes()
            .first()
            .is_some_and(u8::is_ascii_alphabetic)
            || !accession
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'_')
            || version.is_some_and(|v| v.is_empty() || !v.bytes().all(|b| b.is_ascii_digit()))
        {
            return Err(invalid("invalid remote accession/version token"));
        }
        self.offset = end + 1;
        Ok(true)
    }
    fn location(&mut self, depth: usize, remote: bool) -> Result<Option<Location>, Error> {
        self.charge(depth)?;
        if self.remote_prefix()? {
            self.location(depth + 1, true)?;
            return Ok(None);
        }
        if self.take("complement(") {
            let inner = self.location(depth + 1, remote)?;
            self.expect(")")?;
            return Ok(inner.map(|location| Location::Complement(Box::new(location))));
        }
        for (operator, kind) in [("join(", 0), ("order(", 1), ("bond(", 2)] {
            if self.take(operator) {
                let mut children = Vec::new();
                let mut exact = true;
                loop {
                    match self.location(depth + 1, remote)? {
                        Some(child) => children.push(child),
                        None => exact = false,
                    }
                    if !self.take(",") {
                        break;
                    }
                }
                self.expect(")")?;
                return Ok(if !exact || kind == 2 {
                    None
                } else if kind == 0 {
                    Some(Location::Join(children))
                } else {
                    Some(Location::Order(children))
                });
            }
        }
        if self.take("gap(") {
            if !self.take(")") {
                if !self.take("unk100") {
                    self.number()?;
                }
                self.expect(")")?;
            }
            return Ok(None);
        }
        let start = self.position(depth)?;
        self.check_position(start, remote)?;
        if self.take("..") {
            let end = self.position(depth)?;
            return self.local_range(start, end, remote);
        }
        if self.take("^") {
            if start.boundary != Boundary::Exact {
                return Err(invalid("uncertain between-base boundary"));
            }
            let end = self.number()?;
            self.check_position(
                Position {
                    low: end,
                    high: end,
                    boundary: Boundary::Exact,
                },
                remote,
            )?;
            let adjacent = start.low.checked_add(1) == Some(end);
            let origin = end == 1
                && start.low > 1
                && (remote || self.length.is_none_or(|length| start.low == length as u64));
            if !adjacent && !origin {
                return Err(invalid("between-base positions are not adjacent"));
            }
            return Ok(None);
        }
        // Older within-position notation is retained as uncertainty, never an interval.
        if self.take(".") {
            let high = self.number()?;
            if start.boundary != Boundary::Exact || high < start.low {
                return Err(invalid("invalid uncertain position"));
            }
            self.check_position(
                Position {
                    low: start.low,
                    high,
                    boundary: Boundary::Uncertain,
                },
                remote,
            )?;
            return Ok(None);
        }
        if start.boundary != Boundary::Exact {
            return Ok(None);
        }
        self.local_range(start, start, remote)
    }
}
