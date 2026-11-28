//! Explicit documentary assignments; labels and row order never identify a construct.
use collagene_core::{
    documents::workspace::Workspace,
    material::lineage as run_lineage,
    output::layout::{self as native, Assignment, Plate, Request, Source, Well},
};
#[derive(Clone, PartialEq, Eq)]
pub struct Scope {
    pub epoch: u64,
    pub fingerprint: String,
    pub batch_hash: String,
}
impl Scope {
    pub fn capture(workspace: &Workspace, epoch: u64, batch_hash: String) -> Result<Self, String> {
        collagene_core::planning::combinatorial::batch::load(workspace, &batch_hash)
            .map_err(|e| e.to_string())?;
        Ok(Self {
            epoch,
            fingerprint: run_lineage::workspace_fingerprint(workspace)
                .map_err(|e| e.to_string())?,
            batch_hash,
        })
    }
    pub fn permit(&self, workspace: &Workspace, epoch: u64) -> Result<(), String> {
        if epoch != self.epoch
            || run_lineage::workspace_fingerprint(workspace).map_err(|e| e.to_string())?
                != self.fingerprint
        {
            return Err("Workspace or saved batch changed; reopen the layout editor.".into());
        }
        Ok(())
    }
}
#[derive(Clone)]
pub struct Choice {
    pub label: String,
    pub source: Source,
}
fn recorded_status(status: collagene_core::verification::evidence::Status) -> &'static str {
    match status {
        collagene_core::verification::evidence::Status::Invalid => "inconsistent",
        collagene_core::verification::evidence::Status::Incomplete => "incomplete",
        collagene_core::verification::evidence::Status::ConsistentUnderAssumptions => {
            "consistent under stated assumptions"
        }
    }
}
pub fn choices(workspace: &Workspace, hash: &str) -> Result<Vec<Choice>, String> {
    let record = collagene_core::planning::combinatorial::batch::load(workspace, hash)
        .map_err(|e| e.to_string())?;
    let mut choices = record
        .rows
        .into_iter()
        .filter_map(|row| {
            row.declaration_revision.map(|declaration_revision| Choice {
                label: format!(
                    "Planned · {} · historical status {}",
                    row.name,
                    row.status
                        .map(|s| recorded_status(s).to_string())
                        .unwrap_or_else(|| "not recorded".into())
                ),
                source: Source::Planned {
                    batch_hash: hash.into(),
                    construct_id: row.construct_id,
                    declaration_revision,
                },
            })
        })
        .collect::<Vec<_>>();
    for summary in collagene_core::planning::combinatorial::release::list(workspace)
        .map_err(|e| e.to_string())?
    {
        if summary.batch_hash != hash {
            continue;
        }
        let manifest = collagene_core::planning::combinatorial::release::load(
            workspace,
            &summary.manifest_hash,
        )
        .map_err(|e| e.to_string())?;
        for row in manifest.rows {
            choices.push(Choice {
                label: format!(
                    "Released · {} · recorded {} · {}",
                    row.construct_id,
                    recorded_status(row.output.run.recorded_status),
                    row.output.run.recorded_engine
                ),
                source: Source::Released {
                    manifest_hash: summary.manifest_hash.clone(),
                    construct_id: row.construct_id,
                    run_id: row.output.run.run_id,
                    product: row.output.construct,
                },
            });
        }
    }
    Ok(choices)
}

/// Display only: selection retains the complete native source DTO.
pub fn compact_identity(source: &Source) -> String {
    let short = |value: &str| value.chars().take(10).collect::<String>();
    match source {
        Source::Planned {
            construct_id,
            declaration_revision,
            ..
        } => format!(
            "Planned · {construct} · declaration {}…",
            short(declaration_revision),
            construct = construct_id.chars().take(24).collect::<String>()
        ),
        Source::Released {
            construct_id,
            manifest_hash,
            product,
            ..
        } => format!(
            "Released · {construct} · manifest {}… · event {}…",
            short(manifest_hash),
            short(&product.event_id),
            construct = construct_id.chars().take(24).collect::<String>()
        ),
    }
}
pub fn source_identity(source: &Source) -> String {
    match source {
        Source::Planned {
            batch_hash,
            construct_id,
            declaration_revision,
        } => format!(
            "Construct {construct_id} · batch {batch_hash} · declaration {declaration_revision}"
        ),
        Source::Released {
            manifest_hash,
            construct_id,
            run_id,
            product,
        } => format!(
            "Released {construct_id} · manifest {manifest_hash} · run {run_id} · document {} event {} revision {}",
            product.document_id, product.event_id, product.revision
        ),
    }
}
pub fn saved_layouts(
    workspace: &Workspace,
    batch_hash: &str,
) -> Result<Vec<(String, String)>, String> {
    let mut values = vec![];
    for (hash, asset) in &workspace.attachments {
        if asset.kind == native::KIND {
            let layout = native::load(workspace, hash).map_err(|e| e.to_string())?;
            let mut belongs = false;
            for assignment in &layout.request.assignments {
                match &assignment.source {
                    Source::Planned {
                        batch_hash: owner, ..
                    } => belongs |= owner == batch_hash,
                    Source::Released { manifest_hash, .. } => {
                        let manifest = collagene_core::planning::combinatorial::release::load(
                            workspace,
                            manifest_hash,
                        )
                        .map_err(|e| e.to_string())?;
                        belongs |= manifest.batch_hash == batch_hash;
                    }
                }
            }
            if belongs {
                values.push((hash.clone(), layout.request.name));
            }
        }
    }
    Ok(values)
}
pub fn fill_choices(choices: &[Choice], kind: &str) -> Result<Vec<Choice>, String> {
    if !matches!(kind, "planned" | "released") {
        return Err(
            "Choose planned declarations or released product events for the fill preview.".into(),
        );
    }
    Ok(choices
        .iter()
        .filter(|c| {
            matches!(
                (&c.source, kind),
                (Source::Planned { .. }, "planned") | (Source::Released { .. }, "released")
            )
        })
        .cloned()
        .collect())
}
pub fn dimensions(plate: Plate) -> (usize, usize) {
    match plate {
        Plate::Wells96 => (8, 12),
        Plate::Wells384 => (16, 24),
    }
}
pub fn well_name(well: &Well) -> String {
    format!(
        "{}{}",
        char::from(b'A' + well.row),
        usize::from(well.column) + 1
    )
}
pub fn preview_fill(
    plate: Plate,
    choices: &[Choice],
    replicates: &str,
) -> Result<Vec<Assignment>, String> {
    let copies = replicates
        .parse::<usize>()
        .map_err(|_| "Enter an explicit replicate count.")?;
    if !(1..=384).contains(&copies) || choices.is_empty() {
        return Err("Choose at least one resolved construct and 1–384 replicates.".into());
    }
    let (rows, columns) = dimensions(plate);
    let count = choices
        .len()
        .checked_mul(copies)
        .filter(|n| *n <= rows * columns)
        .ok_or("Assignments exceed the chosen plate; no wells changed.")?;
    Ok(choices
        .iter()
        .flat_map(|choice| std::iter::repeat_n(choice, copies))
        .enumerate()
        .take(count)
        .map(|(i, choice)| Assignment {
            id: format!("assignment-{}", i + 1),
            well: Well {
                row: (i / columns) as u8,
                column: (i % columns) as u8,
            },
            label: String::new(),
            source: choice.source.clone(),
        })
        .collect())
}
pub fn request(
    name: String,
    plate: Option<Plate>,
    assignments: Vec<Assignment>,
) -> Result<Request, String> {
    let plate = plate.ok_or("Choose an explicit 96- or 384-well plate.")?;
    if name.trim().is_empty() {
        return Err("Name the documentary layout.".into());
    }
    Ok(Request {
        name,
        plate,
        assignments,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    fn choice(id: &str) -> Choice {
        Choice {
            label: "same display label".into(),
            source: Source::Planned {
                batch_hash: "a".repeat(64),
                construct_id: id.into(),
                declaration_revision: "b".repeat(64),
            },
        }
    }
    #[test]
    fn explicit_replicates_keep_source_ids_and_boundary_wells() {
        let selected = vec![choice("one"), choice("two")];
        let rows = preview_fill(Plate::Wells96, &selected, "2").unwrap();
        assert_eq!(rows.len(), 4);
        assert_ne!(rows[0].id, rows[1].id);
        assert!(matches!(&rows[2].source,Source::Planned{construct_id,..} if construct_id=="two"));
        assert_eq!(well_name(&Well { row: 7, column: 11 }), "H12");
        assert_eq!(
            well_name(&Well {
                row: 15,
                column: 23
            }),
            "P24"
        );
        assert!(preview_fill(Plate::Wells96, &selected, "49").is_err());
        assert!(request("layout".into(), None, vec![]).is_err());
    }
    #[test]
    fn fill_preview_does_not_modify_existing_assignments() {
        let existing: Vec<Assignment> = vec![];
        let proposed = preview_fill(Plate::Wells384, &[choice("unchanged")], "3").unwrap();
        assert!(existing.is_empty());
        assert_eq!(proposed.len(), 3);
        assert_eq!(well_name(&proposed[2].well), "A3");
    }
    #[test]
    fn fill_requires_explicit_identity_kind_and_keeps_exact_release_event() {
        let planned = choice("same");
        let released = Choice {
            label: "Released · same".into(),
            source: Source::Released {
                manifest_hash: "c".repeat(64),
                construct_id: "same".into(),
                run_id: "d".repeat(64),
                product: collagene_core::material::lineage::EventBinding {
                    document_id: "released".into(),
                    event_id: "e".repeat(64),
                    revision: "f".repeat(64),
                },
            },
        };
        let all = vec![planned, released];
        assert!(fill_choices(&all, "").is_err());
        assert_eq!(fill_choices(&all, "planned").unwrap().len(), 1);
        let selected = fill_choices(&all, "released").unwrap();
        assert!(
            matches!(&selected[0].source,Source::Released{product,..} if product.document_id=="released"&&product.event_id=="e".repeat(64))
        );
        let filled = preview_fill(Plate::Wells96, &selected, "2").unwrap();
        assert_eq!(filled.len(), 2);
        assert_ne!(filled[0].id, filled[1].id);
    }
}
