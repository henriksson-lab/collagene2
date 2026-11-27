use crate::*;

fn part(s: &str, start: usize, end: usize) -> &str {
    s.get(start..s.len().min(end)).unwrap_or("").trim()
}
fn optional(s: &str, line: usize) -> Result<Option<String>> {
    if s.is_empty() {
        Ok(None)
    } else {
        field(s, line).map(Some)
    }
}
fn maybe_number(s: &str, line: usize) -> Result<Option<f64>> {
    if s.is_empty() {
        Ok(None)
    } else {
        number(s, line).map(Some)
    }
}
pub(crate) fn parse(text: &str, out: &mut Structure) -> Result<()> {
    let mut model = None;
    let mut explicit_models = false;
    let mut ended = false;
    let mut seen_models = BTreeSet::new();
    // Parsing normalizes line separators only; retained source stays byte-exact.
    let normalized = text.replace("\r\n", "\n").replace('\r', "\n");
    for (i, s) in normalized.lines().enumerate() {
        let line = i + 1;
        if s.trim().is_empty() {
            continue;
        }
        if s.len() > 1024 {
            return Err(Error::Limit("PDB record length"));
        }
        let record = part(s, 0, 6);
        if ended && !s.trim().is_empty() {
            return Err(bad(line, "records after END"));
        }
        match record {
            "MODEL" => {
                if model.is_some() || (!explicit_models && !out.atoms.is_empty()) {
                    return Err(bad(line, "nested MODEL or implicit/explicit model mixture"));
                }
                let n = positive(part(s, 10, 14), line)?;
                if !seen_models.insert(n) {
                    return Err(bad(line, "duplicate MODEL"));
                }
                if seen_models.len() > MAX_MODELS {
                    return Err(Error::Limit("models"));
                }
                model = Some(n);
                explicit_models = true;
            }
            "ENDMDL" => {
                if model.take().is_none() {
                    return Err(bad(line, "ENDMDL without MODEL"));
                }
            }
            "ATOM" | "HETATM" => {
                if explicit_models && model.is_none() {
                    return Err(bad(line, "atom outside MODEL"));
                }
                if s.len() < 54 {
                    return Err(bad(line, "truncated atom coordinates"));
                }
                let serial = positive(part(s, 6, 11), line)?.to_string();
                let atom = part(s, 12, 16);
                let component = part(s, 17, 20);
                if atom.is_empty() || component.is_empty() || part(s, 22, 26).is_empty() {
                    return Err(bad(line, "missing atom/residue identity"));
                }
                let charge = part(s, 78, 80);
                let formal_charge = if charge.is_empty() {
                    None
                } else {
                    let b = charge.as_bytes();
                    if b.len() != 2
                        || !(b'1'..=b'9').contains(&b[0])
                        || !matches!(b[1], b'+' | b'-')
                    {
                        return Err(bad(line, "invalid PDB formal charge"));
                    }
                    Some((b[0] - b'0') as i8 * if b[1] == b'+' { 1 } else { -1 })
                };
                push_atom(
                    out,
                    Atom {
                        id: serial,
                        model: model.unwrap_or(1),
                        record: if record == "ATOM" {
                            Record::Atom
                        } else {
                            Record::Hetatm
                        },
                        label: Identity {
                            atom: None,
                            component: None,
                            chain: None,
                            sequence: None,
                        },
                        author: Identity {
                            atom: optional(atom, line)?,
                            component: optional(component, line)?,
                            chain: optional(part(s, 21, 22), line)?,
                            sequence: optional(part(s, 22, 26), line)?,
                        },
                        entity: None,
                        alternate: optional(part(s, 16, 17), line)?,
                        insertion: optional(part(s, 26, 27), line)?,
                        position: [
                            number(part(s, 30, 38), line)?,
                            number(part(s, 38, 46), line)?,
                            number(part(s, 46, 54), line)?,
                        ],
                        element: element(part(s, 76, 78), line)?,
                        occupancy: maybe_number(part(s, 54, 60), line)?,
                        b_factor: maybe_number(part(s, 60, 66), line)?,
                        formal_charge,
                    },
                )?;
            }
            "CONECT" => {
                let from = positive(part(s, 6, 11), line)?.to_string();
                for col in [11, 16, 21, 26] {
                    let to = part(s, col, col + 5);
                    if !to.is_empty() {
                        let to = positive(to, line)?.to_string();
                        if from == to {
                            return Err(bad(line, "self CONECT"));
                        }
                        push_connection(
                            out,
                            Connection::Pdb {
                                model,
                                from: from.clone(),
                                to,
                            },
                        )?;
                    }
                }
            }
            "END" => {
                if model.is_some() {
                    return Err(bad(line, "END inside MODEL"));
                }
                ended = true;
            }
            _ => {
                if out.notes.ignored_pdb_records.len() >= 4096
                    && !out.notes.ignored_pdb_records.contains_key(record)
                {
                    return Err(Error::Limit("PDB record kinds"));
                }
                *out.notes
                    .ignored_pdb_records
                    .entry(record.into())
                    .or_default() += 1;
            }
        }
    }
    if model.is_some() {
        return Err(bad(1, "unterminated MODEL"));
    }
    Ok(())
}
