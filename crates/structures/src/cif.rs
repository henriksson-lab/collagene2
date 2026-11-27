//! Streaming CIF1 tokenization. Unknown categories are retained in the source,
//! not interpreted as coordinates, dictionaries, assemblies or symmetry copies.
use crate::*;
const MAX_TOKENS: usize = 4_000_000;
const MAX_COLUMNS: usize = 512;
const MAX_TOKEN_BYTES: usize = 1 << 20;

#[derive(Clone, Copy)]
struct Token<'a> {
    value: &'a str,
    quoted: bool,
    line: usize,
}
impl Token<'_> {
    fn missing(self) -> bool {
        !self.quoted && matches!(self.value, "." | "?")
    }
    fn control(self) -> bool {
        !self.quoted
            && (self.value.starts_with('_')
                || self.value.eq_ignore_ascii_case("loop_")
                || self.value.eq_ignore_ascii_case("stop_")
                || self.value.eq_ignore_ascii_case("global_")
                || self.value.get(..5).is_some_and(|s| {
                    s.eq_ignore_ascii_case("data_") || s.eq_ignore_ascii_case("save_")
                }))
    }
    fn text(self) -> Result<Option<String>> {
        if self.missing() {
            Ok(None)
        } else {
            field(self.value, self.line).map(Some)
        }
    }
    fn required(self) -> Result<String> {
        self.text()?
            .filter(|s| !s.is_empty())
            .ok_or_else(|| bad(self.line, "missing required atom field"))
    }
    fn numeric(self) -> Result<f64> {
        if self.quoted || self.missing() {
            return Err(bad(self.line, "missing or quoted numeric field"));
        }
        number(self.value, self.line)
    }
}
struct Lexer<'a> {
    text: &'a str,
    pos: usize,
    line: usize,
    count: usize,
}
impl<'a> Lexer<'a> {
    fn bump(&mut self) {
        let b = self.text.as_bytes();
        if b[self.pos] == b'\r'
            || (b[self.pos] == b'\n' && (self.pos == 0 || b[self.pos - 1] != b'\r'))
        {
            self.line += 1;
        }
        self.pos += 1;
    }
    fn next(&mut self) -> Result<Option<Token<'a>>> {
        let b = self.text.as_bytes();
        let n = b.len();
        loop {
            while self.pos < n && b[self.pos].is_ascii_whitespace() {
                self.bump();
            }
            if self.pos >= n {
                return Ok(None);
            }
            if b[self.pos] != b'#' {
                break;
            }
            while self.pos < n && !matches!(b[self.pos], b'\n' | b'\r') {
                self.bump();
            }
        }
        self.count += 1;
        if self.count > MAX_TOKENS {
            return Err(Error::Limit("CIF tokens"));
        }
        let line = self.line;
        let start = self.pos;
        let quote = b[start];
        let at_line_start = start == 0 || matches!(b[start - 1], b'\n' | b'\r');
        let (value, quoted) = if quote == b';' && at_line_start {
            self.bump();
            let begin = self.pos;
            loop {
                if self.pos >= n {
                    return Err(bad(line, "unterminated CIF text field"));
                }
                if self.pos - start > MAX_TOKEN_BYTES {
                    return Err(Error::Limit("CIF token bytes"));
                }
                if b[self.pos] == b';' && self.pos > 0 && matches!(b[self.pos - 1], b'\n' | b'\r') {
                    let mut end = self.pos;
                    if end > begin && b[end - 1] == b'\n' {
                        end -= 1;
                    }
                    if end > begin && b[end - 1] == b'\r' {
                        end -= 1;
                    }
                    let value = &self.text[begin..end];
                    self.bump();
                    if self.pos < n && !b[self.pos].is_ascii_whitespace() {
                        return Err(bad(line, "text delimiter must be separated"));
                    }
                    break (value, true);
                }
                self.bump();
            }
        } else if matches!(quote, b'\'' | b'"') {
            self.bump();
            let begin = self.pos;
            loop {
                if self.pos >= n || matches!(b[self.pos], b'\n' | b'\r') {
                    return Err(bad(line, "unterminated CIF quoted value"));
                }
                if self.pos - start > MAX_TOKEN_BYTES {
                    return Err(Error::Limit("CIF token bytes"));
                }
                if b[self.pos] == quote
                    && (self.pos + 1 == n || b[self.pos + 1].is_ascii_whitespace())
                {
                    let value = &self.text[begin..self.pos];
                    self.bump();
                    break (value, true);
                }
                self.bump();
            }
        } else {
            if matches!(quote, b'[' | b']' | b'$') {
                return Err(Error::Unsupported("reserved CIF1 value delimiter"));
            }
            while self.pos < n && !b[self.pos].is_ascii_whitespace() {
                self.bump();
                if self.pos - start > MAX_TOKEN_BYTES {
                    return Err(Error::Limit("CIF token bytes"));
                }
            }
            (&self.text[start..self.pos], false)
        };
        Ok(Some(Token {
            value,
            quoted,
            line,
        }))
    }
}
struct Input<'a> {
    lexer: Lexer<'a>,
    pending: Option<Token<'a>>,
}
impl<'a> Input<'a> {
    fn peek(&mut self) -> Result<Option<Token<'a>>> {
        if self.pending.is_none() {
            self.pending = self.lexer.next()?;
        }
        Ok(self.pending)
    }
    fn next(&mut self) -> Result<Option<Token<'a>>> {
        match self.pending.take() {
            Some(t) => Ok(Some(t)),
            None => self.lexer.next(),
        }
    }
}
type Row<'a> = BTreeMap<String, Token<'a>>;
fn atom(row: &Row<'_>, out: &mut Structure, line: usize) -> Result<()> {
    let get = |name: &str| {
        row.get(name)
            .copied()
            .ok_or_else(|| bad(line, "missing required atom_site column"))
    };
    let opt = |name: &str| {
        row.get(name)
            .copied()
            .map(Token::text)
            .transpose()
            .map(Option::flatten)
    };
    let id = get("_atom_site.id")?.required()?;
    let model = match row.get("_atom_site.pdbx_pdb_model_num") {
        Some(t) => {
            if t.quoted || t.missing() {
                return Err(bad(t.line, "invalid model number"));
            }
            positive(t.value, t.line)?
        }
        None => 1,
    };
    let record = match get("_atom_site.group_pdb")?.value {
        "ATOM" => Record::Atom,
        "HETATM" => Record::Hetatm,
        _ => return Err(bad(line, "unknown atom group_PDB")),
    };
    let label = Identity {
        atom: Some(get("_atom_site.label_atom_id")?.required()?),
        component: Some(get("_atom_site.label_comp_id")?.required()?),
        chain: Some(get("_atom_site.label_asym_id")?.required()?),
        sequence: opt("_atom_site.label_seq_id")?,
    };
    let author = Identity {
        atom: opt("_atom_site.auth_atom_id")?,
        component: opt("_atom_site.auth_comp_id")?,
        chain: opt("_atom_site.auth_asym_id")?,
        sequence: opt("_atom_site.auth_seq_id")?,
    };
    let maybe_numeric = |name: &str| -> Result<Option<f64>> {
        match row.get(name) {
            Some(t) if !t.missing() => t.numeric().map(Some),
            _ => Ok(None),
        }
    };
    let formal_charge = match row.get("_atom_site.pdbx_formal_charge") {
        Some(t) if !t.missing() => {
            if t.quoted {
                return Err(bad(t.line, "quoted formal charge"));
            }
            Some(
                t.value
                    .parse()
                    .map_err(|_| bad(t.line, "invalid formal charge"))?,
            )
        }
        _ => None,
    };
    let symbol = opt("_atom_site.type_symbol")?;
    push_atom(
        out,
        Atom {
            id,
            model,
            record,
            label,
            author,
            entity: opt("_atom_site.label_entity_id")?,
            alternate: opt("_atom_site.label_alt_id")?,
            insertion: opt("_atom_site.pdbx_pdb_ins_code")?,
            position: [
                get("_atom_site.cartn_x")?.numeric()?,
                get("_atom_site.cartn_y")?.numeric()?,
                get("_atom_site.cartn_z")?.numeric()?,
            ],
            element: element(symbol.as_deref().unwrap_or(""), line)?,
            occupancy: maybe_numeric("_atom_site.occupancy")?,
            b_factor: maybe_numeric("_atom_site.b_iso_or_equiv")?,
            formal_charge,
        },
    )
}
fn project(category: &str, row: &Row<'_>, out: &mut Structure, line: usize) -> Result<()> {
    match category {
        "_atom_site" => atom(row, out, line),
        "_struct_conn" => {
            let mut fields = BTreeMap::new();
            for (k, v) in row {
                let value = if !v.quoted && v.value == "." {
                    CifValue::NotApplicable
                } else if !v.quoted && v.value == "?" {
                    CifValue::Unknown
                } else {
                    CifValue::Text(field(v.value, v.line)?)
                };
                fields.insert(k.clone(), value);
            }
            push_connection(out, Connection::Mmcif { fields })
        }
        _ => {
            if out.notes.unprojected_cif_categories.len() >= 4096
                && !out.notes.unprojected_cif_categories.contains(category)
            {
                return Err(Error::Limit("CIF categories"));
            }
            out.notes.unprojected_cif_categories.insert(category.into());
            Ok(())
        }
    }
}
fn category(tag: &str) -> &str {
    tag.split_once('.').map(|(s, _)| s).unwrap_or(tag)
}
fn tag_name(t: Token<'_>) -> Result<String> {
    if t.value.len() < 2 {
        return Err(bad(t.line, "empty CIF tag name"));
    }
    field(t.value, t.line).map(|s| s.to_ascii_lowercase())
}
pub(crate) fn parse(text: &str, out: &mut Structure) -> Result<()> {
    if text.starts_with("#\\#CIF_2.") {
        return Err(Error::Unsupported("CIF2"));
    }
    let mut input = Input {
        lexer: Lexer {
            text,
            pos: 0,
            line: 1,
            count: 0,
        },
        pending: None,
    };
    let first = input
        .next()?
        .ok_or_else(|| bad(1, "missing CIF data block"))?;
    if first.quoted
        || !first
            .value
            .get(..5)
            .is_some_and(|s| s.eq_ignore_ascii_case("data_"))
        || first.value.len() == 5
    {
        return Err(bad(first.line, "missing CIF data block"));
    }
    out.block = Some(field(&first.value[5..], first.line)?);
    let mut seen = BTreeSet::new();
    let mut scalar_rows: BTreeMap<String, Row<'_>> = BTreeMap::new();
    let mut looped = BTreeSet::new();
    while let Some(t) = input.next()? {
        if !t.quoted && t.value.eq_ignore_ascii_case("loop_") {
            let mut names = vec![];
            while let Some(tag) = input.peek()? {
                if tag.quoted || !tag.value.starts_with('_') {
                    break;
                }
                input.next()?;
                if names.len() >= MAX_COLUMNS {
                    return Err(Error::Limit("CIF loop columns"));
                }
                let name = tag_name(tag)?;
                if !seen.insert(name.clone()) {
                    return Err(bad(tag.line, "duplicate CIF tag"));
                }
                names.push(name);
            }
            if names.is_empty() {
                return Err(bad(t.line, "loop has no tags"));
            }
            let cat = category(&names[0]).to_string();
            if names.iter().any(|n| category(n) != cat) {
                return Err(Error::Unsupported("mixed-category CIF loop"));
            }
            if scalar_rows.contains_key(&cat) || !looped.insert(cat.clone()) {
                return Err(Error::Unsupported("split-category CIF tables"));
            }
            let mut rows = 0usize;
            while let Some(value) = input.peek()? {
                if value.control() {
                    break;
                }
                let mut row = Row::new();
                for name in &names {
                    let v = input
                        .next()?
                        .ok_or_else(|| bad(t.line, "incomplete CIF loop row"))?;
                    if v.control() {
                        return Err(bad(v.line, "incomplete CIF loop row"));
                    }
                    row.insert(name.clone(), v);
                }
                project(&cat, &row, out, value.line)?;
                rows += 1;
            }
            if rows == 0 {
                return Err(bad(t.line, "empty CIF loop"));
            }
        } else if !t.quoted && t.value.starts_with('_') {
            let name = tag_name(t)?;
            if !seen.insert(name.clone()) {
                return Err(bad(t.line, "duplicate CIF tag"));
            }
            let v = input
                .next()?
                .ok_or_else(|| bad(t.line, "missing CIF item value"))?;
            if v.control() {
                return Err(bad(v.line, "missing CIF item value"));
            }
            let cat = category(&name).to_string();
            if looped.contains(&cat) {
                return Err(Error::Unsupported("split-category CIF tables"));
            }
            // Only projected scalar categories need retained token maps.
            if matches!(cat.as_str(), "_atom_site" | "_struct_conn") {
                scalar_rows.entry(cat).or_default().insert(name, v);
            } else {
                if out.notes.unprojected_cif_categories.len() >= 4096
                    && !out.notes.unprojected_cif_categories.contains(&cat)
                {
                    return Err(Error::Limit("CIF categories"));
                }
                out.notes.unprojected_cif_categories.insert(cat);
            }
        } else {
            return Err(Error::Unsupported(
                "multiple data blocks, dictionary frames or reserved CIF control",
            ));
        }
        if seen.len() > 65_536 {
            return Err(Error::Limit("CIF tags"));
        }
    }
    for (cat, row) in scalar_rows {
        project(&cat, &row, out, 1)?;
    }
    Ok(())
}
