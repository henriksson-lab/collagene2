//! GPU-ready source-indexed points and explicitly deposited connections.
//! No distance-based bonds, folding, symmetry expansion or sequence mapping.
use crate::*;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Chains {
    All,
    Label(String),
    Author(Option<String>),
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Alternates {
    /// Show all deposited alternate locations without claiming coexistence.
    All,
    /// Unlabelled atoms plus this exact alternate ID. No highest-occupancy guess.
    Named(String),
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Selection {
    pub model: u32,
    pub chains: Chains,
    pub alternates: Alternates,
    pub hydrogens: bool,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Point {
    pub atom_index: u32,
    pub position: [f32; 3],
    pub color: [f32; 3],
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Edge {
    /// Original atom indices, never positions in the filtered point array.
    pub atoms: [u32; 2],
    pub declaration_index: u32,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Scene {
    pub source_hash: String,
    pub selection: Selection,
    pub points: Vec<Point>,
    pub edges: Vec<Edge>,
    /// Original coordinates / coordinate_magnitude are centered at this value,
    /// then divided by scaled_radius. No source f64 coordinate is overwritten.
    pub scaled_center: [f64; 3],
    pub coordinate_magnitude: f64,
    pub scaled_radius: f64,
    pub unresolved_connections: usize,
    pub withheld_symmetry_or_relation: usize,
}
pub fn includes(a: &Atom, s: &Selection) -> bool {
    a.model == s.model
        && (s.hydrogens || !matches!(a.element.as_deref(), Some("H" | "D")))
        && match &s.chains {
            Chains::All => true,
            Chains::Label(c) => a.label.chain.as_ref() == Some(c),
            Chains::Author(c) => &a.author.chain == c,
        }
        && match &s.alternates {
            Alternates::All => true,
            Alternates::Named(n) => a.alternate.as_ref().is_none_or(|a| a == n),
        }
}
fn color(element: Option<&str>) -> [f32; 3] {
    match element {
        Some("C") => [0.32, 0.38, 0.46],
        Some("N") => [0.20, 0.39, 0.91],
        Some("O") => [0.90, 0.20, 0.22],
        Some("S") => [0.92, 0.76, 0.12],
        Some("P") => [0.92, 0.45, 0.14],
        Some("H" | "D") => [0.78, 0.81, 0.86],
        Some("F" | "CL") => [0.24, 0.72, 0.32],
        _ => [0.65, 0.42, 0.77],
    }
}
fn present(value: Option<&CifValue>) -> Option<&str> {
    match value {
        Some(CifValue::Text(s)) if !s.is_empty() => Some(s),
        _ => None,
    }
}
fn qualifier_matches(value: Option<&CifValue>, actual: Option<&str>) -> bool {
    match value {
        None => true,
        Some(CifValue::NotApplicable) => actual.is_none(),
        Some(CifValue::Unknown) => false,
        Some(CifValue::Text(s)) => actual == Some(s.as_str()),
    }
}
fn partner(
    fields: &BTreeMap<String, CifValue>,
    n: u8,
    atoms: &[Atom],
    included: &[usize],
) -> Option<usize> {
    let key = |suffix: &str| format!("_struct_conn.ptnr{n}_{suffix}");
    let label_chain = present(fields.get(&key("label_asym_id")));
    let label_atom = present(fields.get(&key("label_atom_id")));
    let author_chain = present(fields.get(&key("auth_asym_id")));
    let label_mode = label_chain.is_some() && label_atom.is_some();
    if !label_mode && author_chain.is_none() {
        return None;
    }
    let mut found = None;
    for &i in included {
        let a = &atoms[i];
        let mut matches = true;
        for (namespace, id) in [("label", &a.label), ("auth", &a.author)] {
            for (suffix, value) in [
                ("atom_id", &id.atom),
                ("comp_id", &id.component),
                ("asym_id", &id.chain),
                ("seq_id", &id.sequence),
            ] {
                matches &= qualifier_matches(
                    fields.get(&key(&format!("{namespace}_{suffix}"))),
                    value.as_deref(),
                );
            }
        }
        matches &= qualifier_matches(
            fields.get(&format!("_struct_conn.pdbx_ptnr{n}_label_alt_id")),
            a.alternate.as_deref(),
        );
        matches &= qualifier_matches(
            fields.get(&format!("_struct_conn.pdbx_ptnr{n}_pdb_ins_code")),
            a.insertion.as_deref(),
        );
        let model = a.model.to_string();
        matches &= qualifier_matches(
            fields.get(&format!("_struct_conn.pdbx_ptnr{n}_pdb_model_num")),
            Some(&model),
        );
        // Author-only records must identify an atom, not just a whole chain.
        if !label_mode && present(fields.get(&key("auth_atom_id"))).is_none() {
            return None;
        }
        if matches {
            if found.is_some() {
                return None;
            }
            found = Some(i);
        }
    }
    found
}
pub fn build(structure: &Structure, selection: Selection) -> Result<Scene> {
    structure.validate()?;
    if !structure.models().contains(&selection.model) {
        return Err(bad(1, "selected model is absent"));
    }
    match &selection.chains {
        Chains::Label(n) | Chains::Author(Some(n)) => {
            field(n, 1)?;
            if n.is_empty() {
                return Err(bad(1, "empty chain selector"));
            }
        }
        _ => {}
    }
    if let Alternates::Named(n) = &selection.alternates {
        if n.is_empty() {
            return Err(bad(1, "empty alternate selector"));
        }
        field(n, 1)?;
        if !structure
            .atoms
            .iter()
            .any(|a| a.model == selection.model && a.alternate.as_ref() == Some(n))
        {
            return Err(bad(1, "selected alternate is absent"));
        }
    }
    let included: Vec<_> = structure
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| includes(a, &selection))
        .map(|(i, _)| i)
        .collect();
    if included.is_empty() {
        return Err(bad(1, "selection has no atoms"));
    }
    // Scale first to avoid finite-source extremes overflowing subtraction or
    // squared distances, on native and WASM alike.
    let magnitude = included
        .iter()
        .flat_map(|&i| structure.atoms[i].position)
        .map(f64::abs)
        .fold(0.0, f64::max)
        .max(f64::MIN_POSITIVE);
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &i in &included {
        for axis in 0..3 {
            let v = structure.atoms[i].position[axis] / magnitude;
            lo[axis] = lo[axis].min(v);
            hi[axis] = hi[axis].max(v);
        }
    }
    let center = std::array::from_fn(|axis| (lo[axis] + hi[axis]) * 0.5);
    let radius = included
        .iter()
        .map(|&i| {
            structure.atoms[i]
                .position
                .iter()
                .enumerate()
                .map(|(axis, v)| (v / magnitude - center[axis]).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .fold(0.0, f64::max)
        .max(1e-12);
    let points = included
        .iter()
        .map(|&i| Point {
            atom_index: i as u32,
            position: std::array::from_fn(|axis| {
                ((structure.atoms[i].position[axis] / magnitude - center[axis]) / radius) as f32
            }),
            color: color(structure.atoms[i].element.as_deref()),
        })
        .collect();
    let index: BTreeMap<_, _> = included
        .iter()
        .map(|&i| (structure.atoms[i].id.as_str(), i))
        .collect();
    let mut scene = Scene {
        source_hash: structure.source_hash.clone(),
        selection,
        points,
        edges: vec![],
        scaled_center: center,
        coordinate_magnitude: magnitude,
        scaled_radius: radius,
        unresolved_connections: 0,
        withheld_symmetry_or_relation: 0,
    };
    let mut pairs = BTreeSet::new();
    // The CIF resolver uses an index below; never scan all atoms for every bond.
    let mut label_index: BTreeMap<(Option<&str>, Option<&str>), Vec<usize>> = BTreeMap::new();
    let mut author_index: BTreeMap<(Option<&str>, Option<&str>), Vec<usize>> = BTreeMap::new();
    for &i in &included {
        let a = &structure.atoms[i];
        label_index
            .entry((a.label.chain.as_deref(), a.label.atom.as_deref()))
            .or_default()
            .push(i);
        author_index
            .entry((a.author.chain.as_deref(), a.author.atom.as_deref()))
            .or_default()
            .push(i);
    }
    let mut work = 0usize;
    for (di, declaration) in structure.connections.iter().enumerate() {
        let endpoints = match declaration {
            Connection::Pdb { model, from, to } => {
                if model.is_some_and(|m| m != scene.selection.model) {
                    continue;
                }
                index
                    .get(from.as_str())
                    .copied()
                    .zip(index.get(to.as_str()).copied())
            }
            Connection::Mmcif { fields } => {
                if !matches!(
                    present(fields.get("_struct_conn.conn_type_id")),
                    Some("covale" | "disulf")
                ) || [1, 2].into_iter().any(|n| {
                    present(fields.get(&format!("_struct_conn.ptnr{n}_symmetry"))) != Some("1_555")
                }) {
                    scene.withheld_symmetry_or_relation += 1;
                    continue;
                }
                let mut resolve = |n: u8| -> Result<Option<usize>> {
                    let get = |s: &str| present(fields.get(&format!("_struct_conn.ptnr{n}_{s}")));
                    let candidates =
                        if get("label_asym_id").is_some() && get("label_atom_id").is_some() {
                            label_index.get(&(get("label_asym_id"), get("label_atom_id")))
                        } else {
                            author_index.get(&(get("auth_asym_id"), get("auth_atom_id")))
                        };
                    let candidates = candidates.map(Vec::as_slice).unwrap_or(&[]);
                    work = work
                        .checked_add(candidates.len())
                        .ok_or(Error::Limit("connection resolution work"))?;
                    if work > 4_000_000 {
                        return Err(Error::Limit("connection resolution work"));
                    }
                    Ok(partner(fields, n, &structure.atoms, candidates))
                };
                let left = resolve(1)?;
                let right = resolve(2)?;
                left.zip(right)
            }
        };
        if let Some((a, b)) = endpoints {
            // Distinct alternate locations are never joined into a hybrid molecule.
            let aa = &structure.atoms[a].alternate;
            let bb = &structure.atoms[b].alternate;
            if a == b || (aa.is_some() && bb.is_some() && aa != bb) {
                scene.unresolved_connections += 1;
                continue;
            }
            let pair = if a < b {
                [a as u32, b as u32]
            } else {
                [b as u32, a as u32]
            };
            if pairs.insert(pair) {
                scene.edges.push(Edge {
                    atoms: pair,
                    declaration_index: di as u32,
                });
            }
        } else {
            scene.unresolved_connections += 1;
        }
    }
    Ok(scene)
}
