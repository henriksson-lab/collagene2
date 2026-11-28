use collagene_core::verification;
use collagene_core::{
    cloning::restriction::{
        catalog as catalog_digest, inspection as restriction_inspection, setups as digest_setups,
    },
    documents::workspace,
};

use leptos::prelude::*;

use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::current_event_binding;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct RestrictionController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl RestrictionController {
    pub(crate) fn navigate_restriction(
        &self,
        assembly_visible: RwSignal<bool>,
        digest_boundary_pick: RwSignal<Option<crate::restriction::digest::BoundaryPick>>,
        molecule_revision: Memo<String>,
        restrictions: crate::restriction::restriction_view::Restrictions,
    ) -> Callback<collagene_core::cloning::restriction::inspection::InspectionSite, ()> {
        let session = self.session;
        let WorkspaceSession {
            documents, active, ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            search_regions,
            ..
        } = self.viewer;
        let ShellState { panel, .. } = self.shell;
        Callback::new(
            move |site: collagene_core::cloning::restriction::inspection::InspectionSite| {
                let parts = site
                    .spans
                    .iter()
                    .map(|p| (p.start, p.end))
                    .collect::<Vec<_>>();
                batch(|| {
                    selected_feature.set(None);
                    restrictions.selected.set(Some(site.index));
                    selection.set(if parts.len() == 1 {
                        parts.first().copied()
                    } else {
                        None
                    });
                    search_regions.set(Some((molecule_revision.get_untracked(), parts)));
                    if assembly_visible.get_untracked() {
                        let enzyme = restrictions.catalog.with_value(|catalog| {
                            catalog
                                .iter()
                                .find(|enzyme| enzyme.name == site.enzyme)
                                .map(|enzyme| {
                                    collagene_core::cloning::restriction::inspection::EnzymeSpec {
                                        name: enzyme.name.clone(),
                                        recognition: enzyme.recognition.clone(),
                                        top_cut: enzyme.top_cut,
                                        bottom_cut: enzyme.bottom_cut,
                                    }
                                })
                        });
                        if let Some(enzyme) = enzyme {
                            let source = documents.with_untracked(|docs| {
                                current_event_binding(&docs[active.get_untracked()])
                            });
                            digest_boundary_pick.set(Some(
                                crate::restriction::digest::BoundaryPick {
                                    source,
                                    enzyme,
                                    recognition_start: site.recognition_start,
                                    reverse: site.reverse,
                                },
                            ));
                        }
                        panel.set("enzymes");
                    } else {
                        panel.set("restriction");
                    }
                });
            },
        )
    }
    pub(crate) fn save_digest_setup(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<
        (
            String,
            collagene_core::cloning::restriction::setups::Recipe,
            collagene_core::material::derivation::digest::DigestPreview,
        ),
        Result<(), String>,
    > {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            attachments,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        Callback::new(
            move |(name, recipe, report): (
                String,
                collagene_core::cloning::restriction::setups::Recipe,
                collagene_core::material::derivation::digest::DigestPreview,
            )|
                  -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                if documents
                    .with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]))
                    != recipe.source
                {
                    return Err(
                        "The active digest source event changed; inspect again before saving."
                            .into(),
                    );
                }
                let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog()
                    .map_err(|e| e.to_string())?;
                let (staged, _) = collagene_core::cloning::restriction::setups::capture(
                    &catalog,
                    &contig_snapshot.run(())?,
                    &recipe,
                    &name,
                    &report,
                )
                .map_err(|e| e.to_string())?;
                staged.validate().map_err(|e| e.to_string())?;
                attachments.set(staged.attachments);
                message.set(
                "Named digest setup retained with its exact source event and historical evidence."
                    .into(),
            );
                Ok(())
            },
        )
    }
}

/// Digest source navigation and pending material-admission dialogs.
#[derive(Clone, Copy)]
pub(crate) struct DigestWorkspaceState {
    pub(crate) pending_digest:
        RwSignal<Option<crate::restriction::digest_product::DigestOpenRequest>>,
    pub(crate) pending_material:
        RwSignal<Option<crate::restriction::digest_product::MaterialOpenRequest>>,
    pub(crate) digest_boundary_pick: RwSignal<Option<crate::restriction::digest::BoundaryPick>>,
}
impl DigestWorkspaceState {
    pub(crate) fn new() -> Self {
        Self {
            pending_digest: RwSignal::new(
                None::<crate::restriction::digest_product::DigestOpenRequest>,
            ),
            pending_material: RwSignal::new(
                None::<crate::restriction::digest_product::MaterialOpenRequest>,
            ),
            digest_boundary_pick: RwSignal::new(None::<crate::restriction::digest::BoundaryPick>),
        }
    }
    pub(crate) fn restore(&self) {
        self.pending_digest.set(None);
        self.pending_material.set(None);
        self.digest_boundary_pick.set(None);
    }
}
