//! Shared conversion of editable preparation rows to existing typed recipes.
use crate::cloning_methods_model::{MethodsState, OverlapPreparation};
use collagene_core::analysis::Primer;

pub(crate) fn prepared_rows(
    state: &MethodsState,
    input_count: usize,
) -> Result<
    (
        Vec<String>,
        Vec<collagene_core::assembly_recipe::FragmentSpec>,
    ),
    String,
> {
    use collagene_core::assembly_recipe::{FragmentSpec, Preparation};
    if input_count != state.rows.len() || input_count == 0 || input_count > 12 {
        return Err("Choose 1–12 ordered original sources".into());
    }
    let roles = state.role_ids(input_count)?;
    let ids = roles
        .iter()
        .map(|role| role.input_id.clone())
        .collect::<Vec<_>>();
    let mut fragments = Vec::with_capacity(input_count);
    for (i, row) in state.rows.iter().enumerate() {
        if row.fragment.len() > 20 || row.forward_name.len() > 256 || row.reverse_name.len() > 256 {
            return Err("Preparation names or fragment number exceed their field bounds".into());
        }
        if [
            &row.enzymes,
            &row.forward,
            &row.reverse_primer,
            &row.forward_tail,
            &row.reverse_tail,
        ]
        .iter()
        .any(|s| s.len() > 65536)
        {
            return Err("Preparation fields exceed the 64 KiB field limit".into());
        }
        let preparation = match row.preparation {
            OverlapPreparation::WholeLinear => Preparation::WholeLinear,
            OverlapPreparation::DeclaredDuplex => Preparation::DeclaredDuplex {
                left: row.left_terminus,
                right: row.right_terminus,
            },
            OverlapPreparation::CatalogDigest | OverlapPreparation::PcrCatalogDigest => {
                let selection = row.catalog_selection.as_ref().ok_or(
                    "Pin a catalog selection using Virtual digest or the PCR-product enzyme selector; then choose its fragment.",
                )?;
                selection.validate().map_err(|e| e.to_string())?;
                let fragment_index = row
                    .fragment
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or("Digest fragment numbers start at 1")?;
                if row.preparation == OverlapPreparation::PcrCatalogDigest {
                    Preparation::PcrCatalogDigest {
                        forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail)
                            .map_err(|e| e.to_string())?,
                        reverse: Primer::new(
                            &row.reverse_name,
                            &row.reverse_primer,
                            &row.reverse_tail,
                        )
                        .map_err(|e| e.to_string())?,
                        selection: selection.clone(),
                        fragment_index,
                        repair: row.repair,
                    }
                } else {
                    Preparation::CatalogDigest {
                        selection: selection.clone(),
                        fragment_index,
                        repair: row.repair,
                    }
                }
            }
            OverlapPreparation::Digest => {
                let enzymes = row
                    .enzymes
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(str::to_owned)
                    .collect::<Vec<_>>();
                if enzymes.is_empty() || enzymes.len() > 64 || enzymes.iter().any(|e| e.len() > 256)
                {
                    return Err("Choose 1–64 explicit restriction enzyme names".into());
                }
                let fragment_index = row
                    .fragment
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or("Digest fragment numbers start at 1")?;
                Preparation::Digest {
                    enzymes,
                    fragment_index,
                    repair: row.repair,
                }
            }
            OverlapPreparation::ExactPcr => Preparation::ExactPcr {
                forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail)
                    .map_err(|e| e.to_string())?,
                reverse: Primer::new(&row.reverse_name, &row.reverse_primer, &row.reverse_tail)
                    .map_err(|e| e.to_string())?,
            },
        };
        fragments.push(FragmentSpec {
            role: roles[i].role.clone(),
            input_id: ids[i].clone(),
            preparation,
            reverse: row.reverse,
        });
    }
    Ok((ids, fragments))
}
