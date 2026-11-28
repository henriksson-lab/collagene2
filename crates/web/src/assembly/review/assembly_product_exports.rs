//! Exchange the guarded actual final reference without installing it as material.
use crate::assembly::review::assembly_product_export_model as model;
use crate::assembly::review::assembly_product_export_model::Format;
use crate::assembly::review::assembly_product_export_model::Prepared;
use crate::io::download::ExportArtifact;
use collagene_model::Molecule;
use leptos::prelude::*;
use std::sync::Arc;

struct Bound {
    identity: String,
    epoch: u64,
    export: Prepared,
}

#[component]
pub fn AssemblyProductExports(
    capture: Callback<(), Result<Molecule, String>>,
    identity: Memo<String>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let prepared = RwSignal::new(None::<Arc<Bound>>);
    let status = RwSignal::new(String::new());
    Effect::new(move |_| {
        let _ = (identity.get(), epoch.get());
        prepared.set(None);
        status.set(String::new());
    });
    let prepare = Callback::new(move |format: Format| {
        prepared.set(None);
        let bound = identity.get_untracked();
        let bound_epoch = epoch.get_untracked();
        let result = (|| -> Result<Prepared, String> {
            let actual = capture.run(())?;
            let export = model::prepare(&actual, format)?;
            if identity.get_untracked() != bound || epoch.get_untracked() != bound_epoch {
                return Err("The product review changed; prepare its export again.".into());
            }
            export.require_source(&capture.run(())?)?;
            Ok(export)
        })();
        match result {
            Ok(export) => {
                prepared.set(Some(Arc::new(Bound {
                    identity: bound,
                    epoch: bound_epoch,
                    export,
                })));
                status.set("Review the format losses below before downloading.".into());
            }
            Err(error) => status.set(error),
        }
    });
    let save = Callback::new(move |loss_report: bool| {
        let Some(bound) = prepared.get_untracked() else {
            return;
        };
        let result = (|| -> Result<crate::io::download::ExportArtifact, String> {
            if identity.get_untracked() != bound.identity || epoch.get_untracked() != bound.epoch {
                return Err("The product review changed; prepare its export again.".into());
            }
            bound.export.require_source(&capture.run(())?)?;
            if loss_report {
                let bytes = bound.export.loss_report(&bound.identity)?;
                ExportArtifact::prepare(
                    &format!("{}.losses.json", bound.export.filename),
                    "application/json",
                    &bytes,
                )
            } else {
                ExportArtifact::prepare(
                    &bound.export.filename,
                    "text/plain;charset=utf-8",
                    bound.export.text.as_bytes(),
                )
            }
        })();
        let dispatch_bound = bound.clone();
        crate::io::download::complete_download_guarded(
            result,
            move || {
                if identity.try_get_untracked().as_ref() != Some(&dispatch_bound.identity)
                    || epoch.try_get_untracked() != Some(dispatch_bound.epoch)
                {
                    return Err("The product review changed before export acceptance.".into());
                }
                dispatch_bound.export.require_source(&capture.run(())?)
            },
            move |result| {
                if status.try_get_untracked().is_none()
                    || identity.try_get_untracked().as_ref() != Some(&bound.identity)
                    || epoch.try_get_untracked() != Some(bound.epoch)
                {
                    return;
                }
                match result {
                    Ok(()) => status.set(
                        if loss_report {
                            "Export loss report downloaded."
                        } else {
                            "Actual product reference downloaded."
                        }
                        .into(),
                    ),
                    Err(error) => {
                        prepared.set(None);
                        status.set(error);
                    }
                }
            },
        );
    });
    view! {
        <div class="assembly-product-exports">
            <p>"Actual final product reference · may include a diagnostic product that differs from the intended target. Export does not create a saved construct or accept the target."</p>
            <div class="assembly-actions">
                <button type="button" data-menu-keep-open="" on:click=move |_| prepare.run(Format::GenBank)>"Prepare GenBank"</button>
                <button type="button" data-menu-keep-open="" on:click=move |_| prepare.run(Format::Fasta)>"Prepare FASTA"</button>
            </div>
            <p role="status">{move || status.get()}</p>
            {move || prepared.get().map(|bound| {
                let label = format!("Download actual product {}", bound.export.format.name());
                let details = format!("{} · {} bytes", bound.export.filename, bound.export.text.len());
                view! {
                    <div>
                        <p>{details}</p>
                        <ul>{bound.export.warnings.iter().cloned().map(|warning| view! { <li>{warning}</li> }).collect_view()}</ul>
                        <div class="assembly-actions">
                            <button type="button" on:click=move |_| save.run(false)>{label}</button>
                            <button type="button" on:click=move |_| save.run(true)>"Download export loss report"</button>
                        </div>
                    </div>
                }
            })}
        </div>
    }
}
