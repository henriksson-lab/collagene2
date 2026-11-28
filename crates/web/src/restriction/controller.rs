use crate::workspace::session::{MAX_WORKSPACE_BYTES, document_bytes, recombination_source_info};
use collagene_core::documents::transfer;

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
    pub(crate) fn saved_enzyme_sets(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Memo<
        Vec<(
            String,
            collagene_core::cloning::restriction::enzyme_sets::Set,
        )>,
    > {
        let WorkspaceSession {
            documents,
            active,
            attachments,
            workspace_name,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        Memo::new(move |_| {
            let frame = WorkspaceSession::snapshot(
                &documents.get(),
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &attachments.get(),
                &workspace_name.get_untracked(),
            );
            match frame.and_then(|frame| {
                collagene_core::cloning::restriction::enzyme_sets::list(&frame)
                    .map_err(|e| e.to_string())
            }) {
                Ok(sets) => sets,
                Err(error) => {
                    message.set(format!("Saved enzyme sets need attention: {error}"));
                    Vec::new()
                }
            }
        })
    }
    pub(crate) fn save_enzyme_set(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<collagene_core::cloning::restriction::enzyme_sets::Set, Result<(), String>> {
        let WorkspaceSession {
            workspace_restoring,
            attachments,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        Callback::new(
        move |set: collagene_core::cloning::restriction::enzyme_sets::Set| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let (staged, _) = collagene_core::cloning::restriction::enzyme_sets::store(
                &contig_snapshot.run(())?,
                &set,
            )
            .map_err(|e| e.to_string())?;
            staged.validate().map_err(|e| e.to_string())?;
            attachments.set(staged.attachments);
            message.set(format!("Saved enzyme set: {}", set.name));
            Ok(())
        },
    )
    }
    pub(crate) fn delete_enzyme_set(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<String, Result<(), String>> {
        let WorkspaceSession {
            workspace_restoring,
            attachments,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        Callback::new(move |hash: String| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let staged = collagene_core::cloning::restriction::enzyme_sets::delete(
                &contig_snapshot.run(())?,
                &hash,
            )
            .map_err(|e| e.to_string())?;
            staged.validate().map_err(|e| e.to_string())?;
            attachments.set(staged.attachments);
            message.set("Deleted saved enzyme set.".into());
            Ok(())
        })
    }
    pub(crate) fn saved_digests(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Memo<Vec<collagene_core::cloning::restriction::setups::Summary>> {
        let WorkspaceSession {
            documents,
            active,
            attachments,
            workspace_name,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        Memo::new(move |_| {
            let assets = attachments.get();
            if !assets
                .values()
                .any(|asset| asset.kind == collagene_core::cloning::restriction::setups::KIND)
            {
                return Vec::new();
            }
            let frame = WorkspaceSession::snapshot(
                &documents.get(),
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &assets,
                &workspace_name.get_untracked(),
            );
            match frame.and_then(|frame| {
                collagene_core::cloning::restriction::setups::list(&frame)
                    .map_err(|e| e.to_string())
            }) {
                Ok(mut setups) => {
                    setups.sort_by(|a, b| {
                        a.name
                            .cmp(&b.name)
                            .then(a.attachment_hash.cmp(&b.attachment_hash))
                    });
                    setups
                }
                Err(error) => {
                    message.set(format!("Saved digest setups need attention: {error}"));
                    Vec::new()
                }
            }
        })
    }
    pub(crate) fn load_digest_setup(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<String, Result<collagene_core::cloning::restriction::setups::Record, String>>
    {
        let WorkspaceSession {
            workspace_restoring,
            ..
        } = self.session;
        Callback::new(
        move |hash: String| -> Result<collagene_core::cloning::restriction::setups::Record, String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            collagene_core::cloning::restriction::setups::load(&contig_snapshot.run(())?, &hash)
                .map_err(|e| e.to_string())
        },
    )
    }
    pub(crate) fn delete_digest_setup(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<String, Result<(), String>> {
        let WorkspaceSession {
            workspace_restoring,
            attachments,
            ..
        } = self.session;
        Callback::new(move |hash: String| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let staged = collagene_core::cloning::restriction::setups::delete(
                &contig_snapshot.run(())?,
                &hash,
            )
            .map_err(|e| e.to_string())?;
            attachments.set(staged.attachments);
            Ok(())
        })
    }
    pub(crate) fn queue_digest_gel(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
        resolve_retained_source: Callback<
            collagene_core::material::lineage::EventBinding,
            Result<collagene_model::Molecule, String>,
        >,
        gel_digest_seeds: RwSignal<Vec<crate::restriction::gel::DigestLaneSeed>>,
    ) -> Callback<(String, collagene_core::cloning::restriction::setups::Recipe), Result<(), String>>
    {
        let WorkspaceSession {
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ShellState { panel, .. } = self.shell;
        Callback::new(
            move |(name, recipe): (
                String,
                collagene_core::cloning::restriction::setups::Recipe,
            )|
                  -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                recipe
                    .validate_for_workspace(&contig_snapshot.run(())?)
                    .map_err(|e| e.to_string())?;
                resolve_retained_source.run(recipe.source.clone())?;
                if gel_digest_seeds.with_untracked(Vec::len) >= 32 {
                    return Err("The pending gel queue reached its 32-lane limit.".into());
                }
                gel_digest_seeds.update(|seeds| {
                    seeds.push(crate::restriction::gel::DigestLaneSeed {
                        name,
                        recipe,
                        epoch: workspace_epoch.get_untracked(),
                    })
                });
                panel.set("gel");
                Ok(())
            },
        )
    }
    pub(crate) fn fork_digest_source(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<
        collagene_core::material::lineage::EventBinding,
        Result<crate::sources::identity::SourceInfo, String>,
    > {
        let WorkspaceSession {
            documents,
            active,
            attachments,
            workspace_name,
            workspace_restoring,
            ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            search_regions,
            ..
        } = self.viewer;
        let ShellState { message, .. } = self.shell;
        Callback::new(move |source:collagene_core::material::lineage::EventBinding| ->Result<crate::sources::identity::SourceInfo,String> {
        if workspace_restoring.get_untracked(){return Err("Wait for the workspace update to finish.".into());}
        let mut docs=documents.get_untracked();
        let original=docs.iter().find(|doc|doc.id==source.document_id).ok_or("The retained digest source is unavailable.")?;
        let molecule=original.revision_molecule(&source.event_id).map_err(|e|e.to_string())?;
        if molecule.revision()!=source.revision{return Err("The saved digest source revision differs from its retained event.".into());}
        let prefix=format!("digest-source-{}",&source.revision[..16]);let mut suffix=1usize;
        let id=loop{let id=format!("{prefix}-{suffix}");if !docs.iter().any(|doc|doc.id==id){break id;}suffix+=1;};
        let fork=original.fork_event(&source.event_id,id).map_err(|e|e.to_string())?;
        let info=recombination_source_info(&fork);let index=docs.len();docs.push(fork);
        if docs.len()>1000||docs.iter().map(|doc|document_bytes(&doc.molecule)).sum::<usize>()>MAX_WORKSPACE_BYTES {
            return Err("The source fork exceeds this workspace's document or size limit.".into());
        }
        let frame=WorkspaceSession::snapshot(&docs,index,candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),
            &attachments.get_untracked(),&workspace_name.get_untracked())?;
        frame.validate().map_err(|e|e.to_string())?;
        batch(||{documents.set(docs);active.set(index);selection.set(None);selected_feature.set(None);search_regions.set(None);});
        message.set("Opened an editable copy of the exact digest source; original source heads and saved evidence are unchanged.".into());
        Ok(info)
    })
    }
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
