use collagene_core::documents::history::transition::{AnalysisOrigin, AnalysisProducer};
use collagene_core::documents::transfer;
use collagene_core::formats::sequence as formats;
use collagene_core::material::lineage::EventBinding;
use collagene_core::planning::plan::DesignPlan;
use collagene_core::protocol;
use collagene_core::verification;
use collagene_core::verification::evidence::VerificationReport;
use collagene_core::{
    alignment,
    documents::workspace,
    formats::{
        alignment as alignment_formats, annotations as annotation_formats, oligos as oligo_formats,
        original::{self as legacy_formats, collagene_lower_strand as collagene_xml_compat},
        sequence::gcg as gcg_sequence,
    },
    material::lineage as run_lineage,
    review::{
        reads::{
            archive as read_alignment_archive, session as read_alignment_session,
            store as read_alignment_store,
        },
        variants::session as variants_session,
    },
};
use collagene_model::Alphabet;
use collagene_model::Molecule;
use leptos::prelude::*;

use crate::workspace::session::WorkspaceSession;
use wasm_bindgen_futures::JsFuture;
use wasm_bindgen_futures::spawn_local;
use web_sys::HtmlInputElement;

use crate::workspace::session::OpenDocument;

use crate::workspace::session::MAX_WORKSPACE_BYTES;

use crate::workspace::session::document_bytes;
use crate::workspace::session::open_histories;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct ImportController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl ImportController {
    pub(crate) fn open_file(
        self,
        annotation_format: RwSignal<collagene_core::formats::annotations::AnnotationFormat>,
        annotation_source: RwSignal<String>,
        candidate: RwSignal<Option<transfer::Candidate>>,
        candidate_events: RwSignal<
            Option<(
                String,
                u64,
                collagene_core::material::lineage::EventBinding,
                Vec<collagene_core::material::lineage::InputEventBinding>,
            )>,
        >,
        existing_oligo: RwSignal<Option<collagene_core::formats::oligos::OligoTable>>,
        fastq_encoding: RwSignal<Option<oligo_formats::QualityEncoding>>,
        last_plan: RwSignal<Option<(DesignPlan, Molecule)>>,
        last_plan_target: RwSignal<
            Option<(String, u64, collagene_core::material::lineage::EventBinding)>,
        >,
        pending_alignment: RwSignal<
            Option<crate::viewer::alignment::alignment_view::PendingAlignment>,
        >,
        pending_compat_import: RwSignal<
            Option<(
                std::sync::Arc<Vec<u8>>,
                String,
                crate::io::sequence_compat_import::Profile,
            )>,
        >,
        pending_fastq: RwSignal<Option<std::sync::Arc<Vec<u8>>>>,
        pending_oligo: RwSignal<Option<(String, oligo_formats::Delimiter)>>,
        pending_read_alignment: RwSignal<
            Option<(
                std::sync::Arc<Vec<u8>>,
                collagene_core::review::reads::session::FileFormat,
                String,
            )>,
        >,
        pending_sbol: RwSignal<Option<std::sync::Arc<Vec<u8>>>>,
        pending_variants: RwSignal<
            Option<(
                std::sync::Arc<Vec<u8>>,
                collagene_core::review::variants::session::FileFormat,
                String,
            )>,
        >,
        plan_generation: RwSignal<u32>,
        plan_report: RwSignal<Option<VerificationReport>>,
    ) -> impl Fn(web_sys::Event) + Copy + 'static {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        let ShellState {
            panel,
            message,
            modal,
            modal_error,
            ..
        } = self.shell;
        move |event: web_sys::Event| {
            let input = event_target::<HtmlInputElement>(&event);
            let Some(files) = input.files() else {
                return;
            };
            let files: Vec<_> = (0..files.length()).filter_map(|i| files.get(i)).collect();
            if workspace_restoring.get_untracked() {
                input.set_value("");
                return;
            }
            let opened_epoch = workspace_epoch.get_untracked();
            spawn_local(async move {
                let mut imported = vec![];
                let mut warnings = vec![];
                let file_count = files.len();
                let mut requested_target = None::<String>;
                for file in files {
                    if file.size() > (formats::MAX_IMPORT_BYTES + 42) as f64 {
                        message.set("File exceeds the 32 MiB payload import limit.".into());
                        return;
                    }
                    let bytes = match JsFuture::from(file.array_buffer()).await {
                        Ok(buffer) => js_sys::Uint8Array::new(&buffer).to_vec(),
                        Err(_) => {
                            message.set("Could not read file.".into());
                            return;
                        }
                    };
                    if workspace_restoring.get_untracked()
                        || workspace_epoch.get_untracked() != opened_epoch
                    {
                        message.set("The workspace changed while reading files. Reopen the files to import them here.".into());
                        return;
                    }
                    let extension = file
                        .name()
                        .rsplit('.')
                        .next()
                        .unwrap_or("")
                        .to_ascii_lowercase();
                    let lower_name = file.name().to_ascii_lowercase();
                    if extension == "sbol"
                        || ((extension == "jsonld" || extension == "json")
                            && crate::io::sbol_exchange::looks_like_sbol(&bytes))
                    {
                        if file_count != 1 {
                            message.set("Open an SBOL graph separately to review its component and sequence choices; no files were imported.".into());
                            return;
                        }
                        if bytes.len() > collagene_design_exchange::MAX_BYTES {
                            message.set("SBOL source exceeds the 16 MiB exchange bound.".into());
                            return;
                        }
                        pending_sbol.set(Some(std::sync::Arc::new(bytes)));
                        modal_error.set(String::new());
                        modal.set("sbol");
                        return;
                    }

                    let vcf_source = lower_name.ends_with(".vcf")
                        || lower_name.ends_with(".vcf.gz")
                        || bytes.starts_with(b"##fileformat=VCF");
                    let bcf_source =
                        lower_name.ends_with(".bcf") || lower_name.ends_with(".bcf.gz");
                    if vcf_source || bcf_source {
                        if file_count != 1 {
                            message.set("Open a VCF/BCF source separately to inspect its declared alleles; no files were imported.".into());
                            return;
                        }
                        if bytes.len() > 16 * 1024 * 1024 {
                            message.set(
                                "This browser variant profile accepts source files up to 16 MiB."
                                    .into(),
                            );
                            return;
                        }
                        let format = if bcf_source {
                            collagene_core::review::variants::session::FileFormat::Bcf
                        } else {
                            collagene_core::review::variants::session::FileFormat::Vcf
                        };
                        pending_variants.set(Some((
                            std::sync::Arc::new(bytes),
                            format,
                            file.name(),
                        )));
                        modal_error.set(String::new());
                        modal.set("variants");
                        return;
                    }
                    let sam_header = [
                        b"@HD\t".as_slice(),
                        b"@SQ\t".as_slice(),
                        b"@RG\t".as_slice(),
                        b"@PG\t".as_slice(),
                        b"@CO\t".as_slice(),
                    ]
                    .iter()
                    .any(|prefix| bytes.starts_with(prefix));
                    if matches!(extension.as_str(), "sam" | "bam")
                        || (sam_header && !matches!(extension.as_str(), "fastq" | "fq"))
                    {
                        if file_count != 1 {
                            message.set("Open a SAM/BAM source separately to inspect its records and explicitly bind DNA references; no files were imported.".into());
                            return;
                        }
                        if bytes.len() > 16 * 1024 * 1024 {
                            message.set("This browser read-alignment profile accepts source files up to 16 MiB.".into());
                            return;
                        }
                        let format = if extension == "bam" {
                            collagene_core::review::reads::session::FileFormat::Bam
                        } else {
                            collagene_core::review::reads::session::FileFormat::Sam
                        };
                        pending_read_alignment.set(Some((
                            std::sync::Arc::new(bytes),
                            format,
                            file.name(),
                        )));
                        modal_error.set(String::new());
                        modal.set("read_alignments");
                        return;
                    }
                    if bytes.first() == Some(&b'@') {
                        if file_count != 1 {
                            message.set("Open FASTQ separately to choose its supplied quality encoding; no files were imported.".into());
                            return;
                        }
                        pending_fastq.set(Some(std::sync::Arc::new(bytes)));
                        fastq_encoding.set(None);
                        modal_error.set(String::new());
                        modal.set("fastq");
                        return;
                    }
                    if matches!(
                        extension.as_str(),
                        "aln"
                            | "clustal"
                            | "sto"
                            | "stockholm"
                            | "afa"
                            | "afasta"
                            | "cg2aln"
                            | "phy"
                            | "phylip"
                            | "msf"
                            | "nex"
                            | "nexus"
                            | "pir"
                            | "nbrf"
                            | "gde"
                            | "slx"
                            | "selex"
                            | "sel"
                    ) || bytes.starts_with(b"CLUSTAL")
                        || bytes.starts_with(b"# STOCKHOLM")
                        || bytes.starts_with(b"!!AA_MULTIPLE_ALIGNMENT")
                        || bytes.starts_with(b"!!NA_MULTIPLE_ALIGNMENT")
                        || [b">P1;", b">F1;", b">DL;", b">RL;", b">DC;", b">RC;"]
                            .iter()
                            .any(|header| bytes.starts_with(*header))
                        || bytes
                            .get(..6)
                            .is_some_and(|header| header.eq_ignore_ascii_case(b"#NEXUS"))
                    {
                        if file_count != 1 {
                            message.set("Open an alignment separately so its alphabet can be chosen explicitly; no files were imported.".into());
                            return;
                        }
                        let format = if bytes.starts_with(b"CLUSTAL")
                            || matches!(extension.as_str(), "aln" | "clustal")
                        {
                            Some(collagene_core::formats::alignment::AlignmentFormat::Clustal)
                        } else if bytes.starts_with(b"# STOCKHOLM")
                            || matches!(extension.as_str(), "sto" | "stockholm")
                        {
                            Some(collagene_core::formats::alignment::AlignmentFormat::Stockholm)
                        } else if matches!(extension.as_str(), "afa" | "afasta") {
                            Some(collagene_core::formats::alignment::AlignmentFormat::AlignedFasta)
                        } else if extension == "msf"
                            || bytes.starts_with(b"!!AA_MULTIPLE_ALIGNMENT")
                            || bytes.starts_with(b"!!NA_MULTIPLE_ALIGNMENT")
                        {
                            Some(collagene_core::formats::alignment::AlignmentFormat::Msf)
                        } else if matches!(extension.as_str(), "nex" | "nexus")
                            || bytes
                                .get(..6)
                                .is_some_and(|header| header.eq_ignore_ascii_case(b"#NEXUS"))
                        {
                            Some(collagene_core::formats::alignment::AlignmentFormat::Nexus)
                        } else if matches!(extension.as_str(), "pir" | "nbrf")
                            || [b">P1;", b">F1;", b">DL;", b">RL;", b">DC;", b">RC;"]
                                .iter()
                                .any(|header| bytes.starts_with(*header))
                        {
                            Some(collagene_core::formats::alignment::AlignmentFormat::Pir)
                        } else if extension == "gde" {
                            Some(collagene_core::formats::alignment::AlignmentFormat::Gde)
                        } else if matches!(extension.as_str(), "slx" | "selex" | "sel") {
                            Some(collagene_core::formats::alignment::AlignmentFormat::Selex)
                        } else {
                            None
                        };
                        pending_alignment.set(Some(
                            crate::viewer::alignment::alignment_view::PendingAlignment {
                                name: file.name(),
                                bytes: std::sync::Arc::new(bytes),
                                format,
                                native: extension == "cg2aln",
                            },
                        ));
                        modal_error.set(String::new());
                        modal.set("alignments");
                        return;
                    }
                    if matches!(extension.as_str(), "gff3" | "gff" | "gtf" | "bed") {
                        if bytes.len() > collagene_core::formats::annotations::MAX_ANNOTATION_BYTES
                        {
                            message.set("Annotation source exceeds the 8 MiB import limit.".into());
                            return;
                        }
                        if file_count != 1 {
                            message.set("Open an annotation file separately so its reference mapping can be reviewed; no files were imported.".into());
                            return;
                        }
                        match String::from_utf8(bytes) {
                            Ok(text) => {
                                annotation_source.set(text);
                                annotation_format.set(match extension.as_str() {
                                    "bed" => {
                                        collagene_core::formats::annotations::AnnotationFormat::Bed
                                    }
                                    "gtf" => {
                                        collagene_core::formats::annotations::AnnotationFormat::Gtf
                                    }
                                    _ => {
                                        collagene_core::formats::annotations::AnnotationFormat::Gff3
                                    }
                                });
                                panel.set("annotation_files");
                                message.set("Choose the exact annotation reference ID and review the proposed coordinates before applying.".into());
                            }
                            Err(_) => message.set("Annotation files must be UTF-8 text.".into()),
                        }
                        return;
                    }
                    if extension == "csv" || extension == "tsv" {
                        if file_count != 1 {
                            message.set("Open an oligo table separately to map its columns; no files were imported.".into());
                            return;
                        }
                        match String::from_utf8(bytes) {
                            Ok(text) => {
                                pending_oligo.set(Some((
                                    text,
                                    if extension == "tsv" {
                                        oligo_formats::Delimiter::Tab
                                    } else {
                                        oligo_formats::Delimiter::Comma
                                    },
                                )));
                                existing_oligo.set(None);
                                modal_error.set(String::new());
                                modal.set("oligos");
                            }
                            Err(_) => message.set("Oligo tables must be UTF-8 text.".into()),
                        }
                        return;
                    }
                    if bytes.starts_with(b"CGP\0") {
                        if file_count != 1 {
                            message.set("Open a candidate separately from sequence files.".into());
                            return;
                        }
                        match transfer::load_candidate(&bytes) {
                            Ok(imported_candidate) => {
                                candidate_events.set(None);
                                candidate.set(Some(imported_candidate));
                                plan_generation.update(|g| *g = g.wrapping_add(1));
                                plan_report.set(None);
                                last_plan.set(None);
                                last_plan_target.set(None);
                                panel.set("plan");
                                message.set("Candidate opened for review. Select its intended target and verify against your open documents.".into());
                            }
                            Err(e) => message.set(e.to_string()),
                        }
                        return;
                    } else if bytes.starts_with(b"CGD\0") {
                        match transfer::load_workspace(&bytes) {
                            Ok(snapshot) => {
                                requested_target = Some(snapshot.target_revision);
                                imported.extend(snapshot.documents);
                            }
                            Err(e) => {
                                message.set(e.to_string());
                                return;
                            }
                        }
                    } else if bytes.starts_with(b"CG2\0") {
                        match protocol::load_document(&bytes) {
                            Ok(m) => imported.push(m),
                            Err(e) => {
                                message.set(e.to_string());
                                return;
                            }
                        }
                    } else {
                        let compat_profile = if collagene_core::formats::sequence::gcg::is_gcg(&bytes) {
                        Some(crate::io::sequence_compat_import::Profile::Gcg)
                    } else if std::str::from_utf8(&bytes).ok().is_some_and(|text| {
                        text.trim_start_matches('\u{feff}')
                            .trim_start()
                            .starts_with('>')
                    }) {
                        let hint = match extension.as_str() {
                            "faa" => Some(Alphabet::Protein),
                            "frn" => Some(Alphabet::Rna),
                            _ => None,
                        };
                        Some(crate::io::sequence_compat_import::Profile::Fasta(hint))
                    } else if collagene_core::formats::original::is_collagene_xml(&bytes)
                        && collagene_core::formats::original::collagene_lower_strand::needs_lower_only_review(&bytes)
                            .unwrap_or(false)
                    {
                        Some(crate::io::sequence_compat_import::Profile::CollageneXml)
                    } else {
                        None
                    };
                        if let Some(profile) = compat_profile {
                            if file_count != 1 {
                                message.set("Open this sequence file separately to review its import profile; no files were imported.".into());
                                return;
                            }
                            pending_compat_import.set(Some((
                                std::sync::Arc::new(bytes),
                                file.name(),
                                profile,
                            )));
                            modal_error.set(String::new());
                            modal.set("compat_import");
                            return;
                        }
                        let report = if collagene_core::formats::original::is_xdna(&bytes) {
                            collagene_core::formats::original::import_xdna(&bytes, &file.name())
                        } else {
                            formats::import_bytes(&bytes)
                        };
                        match report {
                            Ok(report) => {
                                imported.extend(report.molecules);
                                warnings.extend(report.warnings);
                            }
                            Err(e) => {
                                message.set(format!("{}: {e}", file.name()));
                                return;
                            }
                        }
                    }
                    if imported.len() > 1000
                        || imported.iter().map(document_bytes).sum::<usize>() > MAX_WORKSPACE_BYTES
                    {
                        message.set("Import exceeds this prototype's 1,000-document / 64 MiB workspace limit. Open a smaller selection of files.".into());
                        return;
                    }
                }
                if let Err(error) = session.add_documents(imported) {
                    message.set(error);
                    return;
                }
                if let Some(revision) = requested_target {
                    if let Some(index) = documents.with_untracked(|docs| {
                        docs.iter().position(|d| d.molecule.revision() == revision)
                    }) {
                        active.set(index);
                    }
                }
                message.set(warnings.join(" "));
            });
            input.set_value("");
        }
    }
    pub(crate) fn open_variant_preview(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
        variant_asset: impl Fn() -> Result<
            (
                String,
                collagene_core::documents::workspace::WorkspaceAttachment,
            ),
            String,
        > + Copy
        + Send
        + Sync
        + 'static,
    ) -> Callback<(Molecule, String), Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        let ShellState { message, modal, .. } = self.shell;
        Callback::new(
            move |(mut molecule, identity): (Molecule, String)| -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for workspace restoration.".into());
                }
                collagene_model::bounds::model_preflight(&[&molecule])
                    .map_err(|e| e.to_string())?;
                if molecule
                    .metadata
                    .get("source:legacy:chemistry_unmapped")
                    .map(String::as_str)
                    != Some("true")
                    || molecule
                        .metadata
                        .get("variant:reference_only")
                        .map(String::as_str)
                        != Some("true")
                {
                    return Err(
                        "The allele preview lacks its reference-only chemistry restriction.".into(),
                    );
                }
                let (hash, asset) = variant_asset()?;
                if molecule.metadata.get("variant:unverified_source") != Some(&hash) {
                    return Err(
                        "The variant source changed. Check the selected allele again.".into(),
                    );
                }
                molecule.name = format!(
                    "{} · allele design",
                    molecule.name.chars().take(200).collect::<String>()
                );
                molecule.validate().map_err(|e| e.to_string())?;
                let mut assets = attachments.get_untracked();
                if assets.get(&hash).is_some_and(|a| a.kind != asset.kind) {
                    return Err("Original bytes are retained under another source format.".into());
                }
                assets.entry(hash).or_insert(asset);
                let mut staged=documents.with_untracked(|docs|->Result<Vec<OpenDocument>,String>{
            if docs.len()>=1000||docs.iter().map(|d|document_bytes(&d.molecule)).sum::<usize>().saturating_add(document_bytes(&molecule))>MAX_WORKSPACE_BYTES{return Err("This allele design exceeds the open-document budget.".into());}
            let source=docs.iter().find(|doc|format!("{}:{}:{}",workspace_epoch.get_untracked(),doc.id,doc.head_event())==identity).ok_or("The exact reference document event changed. Check the allele again.")?;
            let mut suffix=1usize;let id=loop{let id=format!("variant-design-{}-{suffix}",&source.head_revision()[..12]);if !docs.iter().any(|d|d.id==id){break id;}suffix+=1;};
            let derived={
                let analysis_origin = AnalysisOrigin::preview(AnalysisProducer::SelectedAllele, molecule.revision());
                source.derive_event(source.head_event(), id, molecule, "Selected-allele design preview; file evidence is unverified and confers no prepared chemistry", analysis_origin)
            }.map_err(|e|e.to_string())?;
            let mut staged=docs.clone();staged.push(derived);Ok(staged)
        })?;
                let index = staged.len() - 1;
                let frame = WorkspaceSession::snapshot(
                    &staged,
                    index,
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &assets,
                    &workspace_name.get_untracked(),
                )?;
                frame.validate().map_err(|e| e.to_string())?;
                batch(|| {
                    documents.set(std::mem::take(&mut staged));
                    attachments.set(assets);
                    active.set(index);
                    modal.set("");
                    message.set("Allele design opened with exact reference-event ancestry and original variant bytes. It remains a reference-only design target.".into());
                });
                Ok(())
            },
        )
    }
    pub(crate) fn keep_read_source(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        pending_read_alignment: RwSignal<
            Option<(
                std::sync::Arc<Vec<u8>>,
                collagene_core::review::reads::session::FileFormat,
                String,
            )>,
        >,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<(), Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        Callback::new(move |()| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for workspace restoration.".into());
            }
            let Some((source, format, name)) = pending_read_alignment.get_untracked() else {
                return Err("No read-alignment source is open.".into());
            };
            let format = match format {
                collagene_core::review::reads::session::FileFormat::Sam => {
                    collagene_core::review::reads::store::Format::Sam
                }
                collagene_core::review::reads::session::FileFormat::Bam => {
                    collagene_core::review::reads::store::Format::Bam
                }
            };
            let (hash, asset) =
                collagene_core::review::reads::archive::make_attachment(&source, format, &name)
                    .map_err(|e| e.to_string())?;
            if attachments.with_untracked(|assets| assets.contains_key(&hash)) {
                return Ok(());
            }
            let mut assets = attachments.get_untracked();
            assets.insert(hash, asset);
            let frame = documents.with_untracked(|docs| {
                WorkspaceSession::snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &assets,
                    &workspace_name.get_untracked(),
                )
            })?;
            frame.validate().map_err(|e| e.to_string())?;
            attachments.set(assets);
            message.set("Original read-alignment source retained in workspace. Reopening reparses the source; reference checks are not saved as authority.".into());
            Ok(())
        })
    }
    pub(crate) fn keep_alignment(
        &self,
        alignments: RwSignal<Vec<crate::viewer::alignment::alignment_view::AlignmentRevision>>,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::viewer::alignment::alignment_view::AlignmentRevision, Result<String, String>>
    {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_name,
            attachments,
            ..
        } = self.session;
        Callback::new(
            move |item: crate::viewer::alignment::alignment_view::AlignmentRevision| {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for workspace restoration to finish.".into());
                }
                if let Some(key) = alignments.with_untracked(|items| {
                    items
                        .iter()
                        .find(|a| a.revision == item.revision)
                        .map(|a| a.asset_hash.clone())
                }) {
                    return Ok(key);
                }
                if alignments.with_untracked(Vec::len) >= 128 {
                    return Err("This web workspace retains at most128 alignment collections. Export a workspace and open fewer collections.".into());
                }
                let mut assets = attachments.get_untracked();
                let key = item.asset_hash.clone();
                assets.insert(
                    key.clone(),
                    collagene_core::documents::workspace::WorkspaceAttachment {
                        kind: crate::viewer::alignment::alignment_view::ATTACHMENT_KIND.into(),
                        name: item.prepared.document().name.clone(),
                        media_type: "application/json".into(),
                        bytes: item.bytes.as_ref().clone(),
                    },
                );
                let frame = documents.with_untracked(|docs| {
                    WorkspaceSession::snapshot(
                        docs,
                        active.get_untracked(),
                        candidate.get_untracked().as_ref(),
                        &stored_candidates.get_untracked(),
                        &assets,
                        &workspace_name.get_untracked(),
                    )
                })?;
                frame.validate().map_err(|e| e.to_string())?;
                batch(|| {
                    attachments.set(assets);
                    alignments.update(|items| items.push(item));
                });
                Ok(key)
            },
        )
    }
    pub(crate) fn import_oligos(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        oligo_tables: RwSignal<Vec<std::sync::Arc<collagene_core::formats::oligos::OligoTable>>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<(oligo_formats::OligoTable, Vec<usize>), Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        let ShellState { message, modal, .. } = self.shell;
        Callback::new(
            move |(table, rows): (oligo_formats::OligoTable, Vec<usize>)| {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for workspace restoration to finish.".into());
                }
                table.validate().map_err(|e| e.to_string())?;
                let report = table.selected_molecules(&rows).map_err(|e| e.to_string())?;
                let mut staged_docs = documents.get_untracked();
                let first = staged_docs.len();
                let additions = open_histories(report.molecules, &staged_docs)?;
                staged_docs.extend(additions);
                if staged_docs.len() > 1000
                    || staged_docs
                        .iter()
                        .map(|d| document_bytes(&d.molecule))
                        .sum::<usize>()
                        > MAX_WORKSPACE_BYTES
                {
                    return Err(
                        "Oligo selection exceeds the workspace document limit. Select fewer rows."
                            .into(),
                    );
                }
                let mut staged_attachments = attachments.get_untracked();
                let bytes = serde_json::to_vec(&table).map_err(|e| e.to_string())?;
                let label = format!(
                    "Oligo collection {}",
                    table.source_blake3.chars().take(12).collect::<String>()
                );
                staged_attachments.insert(
                    blake3::hash(&bytes).to_hex().to_string(),
                    collagene_core::documents::workspace::WorkspaceAttachment {
                        kind: "oligo-table-v1".into(),
                        name: label,
                        media_type: "application/json".into(),
                        bytes,
                    },
                );
                let staged_active = if rows.is_empty() {
                    active.get_untracked()
                } else {
                    first
                };
                let frame = WorkspaceSession::snapshot(
                    &staged_docs,
                    staged_active,
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &staged_attachments,
                    &workspace_name.get_untracked(),
                )?;
                frame.validate().map_err(|e| e.to_string())?;
                let table = std::sync::Arc::new(table);
                batch(|| {
                    documents.set(staged_docs);
                    active.set(staged_active);
                    attachments.set(staged_attachments);
                    oligo_tables.update(|tables| {
                        if let Some(old) = tables.iter_mut().find(|old| {
                            old.source_blake3 == table.source_blake3 && old.options == table.options
                        }) {
                            *old = table;
                        } else {
                            tables.push(table);
                        }
                    });
                    modal.set("");
                    message.set(format!(
                        "Oligo collection retained; {} selected primers opened. {}",
                        rows.len(),
                        report.warnings.join(" ")
                    ));
                });
                Ok(())
            },
        )
    }
}

/// Imported oligo collections and pending file-admission dialogs.
#[derive(Clone, Copy)]
pub(crate) struct ImportState {
    pub(crate) oligo_tables:
        RwSignal<Vec<std::sync::Arc<collagene_core::formats::oligos::OligoTable>>>,
    pub(crate) existing_oligo: RwSignal<Option<collagene_core::formats::oligos::OligoTable>>,
    pub(crate) pending_oligo:
        RwSignal<Option<(String, collagene_core::formats::oligos::Delimiter)>>,
    pub(crate) pending_fastq: RwSignal<Option<std::sync::Arc<Vec<u8>>>>,
    pub(crate) pending_sbol: RwSignal<Option<std::sync::Arc<Vec<u8>>>>,
    pub(crate) pending_variants: RwSignal<
        Option<(
            std::sync::Arc<Vec<u8>>,
            collagene_core::review::variants::session::FileFormat,
            String,
        )>,
    >,
    pub(crate) pending_read_alignment: RwSignal<
        Option<(
            std::sync::Arc<Vec<u8>>,
            collagene_core::review::reads::session::FileFormat,
            String,
        )>,
    >,
}
impl ImportState {
    pub(crate) fn new() -> Self {
        Self {
            oligo_tables: RwSignal::new(Vec::<
                std::sync::Arc<collagene_core::formats::oligos::OligoTable>,
            >::new()),
            existing_oligo: RwSignal::new(None::<collagene_core::formats::oligos::OligoTable>),
            pending_oligo: RwSignal::new(None::<(String, oligo_formats::Delimiter)>),
            pending_fastq: RwSignal::new(None::<std::sync::Arc<Vec<u8>>>),
            pending_sbol: RwSignal::new(None::<std::sync::Arc<Vec<u8>>>),
            pending_variants: RwSignal::new(
                None::<(
                    std::sync::Arc<Vec<u8>>,
                    collagene_core::review::variants::session::FileFormat,
                    String,
                )>,
            ),
            pending_read_alignment: RwSignal::new(
                None::<(
                    std::sync::Arc<Vec<u8>>,
                    collagene_core::review::reads::session::FileFormat,
                    String,
                )>,
            ),
        }
    }
    pub(crate) fn restore(
        &self,
        oligos: Vec<std::sync::Arc<collagene_core::formats::oligos::OligoTable>>,
    ) {
        self.oligo_tables.set(oligos);
        self.existing_oligo.set(None);
        self.pending_oligo.set(None);
        self.pending_fastq.set(None);
        self.pending_sbol.set(None);
        self.pending_variants.set(None);
        self.pending_read_alignment.set(None);
    }
}
