use crate::server::CollageneServer;
use crate::server::response::bad;
use crate::server::response::bounded_json;
use crate::server::response::output;
use crate::session::genome::NativeSearch;
use crate::session::genome::native_summary;
use crate::session::genome_store;
use base64::Engine as _;
use base64::engine::general_purpose::STANDARD;
use rmcp::ErrorData;
use rmcp::handler::server::wrapper::Parameters;
use rmcp::model::CallToolResult;
use rmcp::tool;
use rmcp::tool_router;
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub(crate) struct GenomeBeginQuery {
    pub format: genome_store::ImportFormat,
    pub expected_bytes: u64,
    pub expected_hash: Option<String>,
    pub build: Option<collagene_core::tracks::genome::BuildInfo>,
    pub options: Option<collagene_core::tracks::genome::IngestOptions>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeAppendQuery {
    pub import_id: String,
    pub offset: u64,
    pub data_base64: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeFinishQuery {
    pub import_id: String,
    #[serde(default = "genome_finish_default")]
    pub max_blocks: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeInspectQuery {
    pub genome: Option<collagene_core::tracks::genome_session::GenomeRef>,
    #[serde(default)]
    pub contig_start: usize,
    #[serde(default = "genome_page_default")]
    pub max_contigs: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeExportQuery {
    pub genome: collagene_core::tracks::genome_session::GenomeRef,
    #[serde(default)]
    pub offset: u64,
    #[serde(default = "genome_chunk_default")]
    pub max_bytes: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeGuidesQuery {
    pub genome: collagene_core::tracks::genome_session::GenomeRef,
    pub options: collagene_core::tracks::crispr::GenomeGuideOptions,
    #[serde(default)]
    pub candidate_start: usize,
    #[serde(default = "genome_page_default")]
    pub max_candidates: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeReadQuery {
    pub genome: collagene_core::tracks::genome_session::GenomeRef,
    pub options: collagene_core::tracks::crispr::GenomeIntervalReadOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeSearchQuery {
    pub genome: collagene_core::tracks::genome_session::GenomeRef,
    pub query: collagene_core::sequence::crispr::GuideQuery,
    #[serde(default)]
    pub options: collagene_core::tracks::crispr::GenomeSearchOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeAdvanceQuery {
    pub job_id: String,
    #[serde(default = "genome_advance_default")]
    pub max_blocks: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomePageQuery {
    pub job_id: String,
    #[serde(default)]
    pub offset: usize,
    #[serde(default = "genome_page_default")]
    pub limit: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeCancelQuery {
    pub id: String,
    pub genome: Option<collagene_core::tracks::genome_session::GenomeRef>,
}
#[tool_router(router = genome_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Begin a session-owned large-reference import from bounded FASTA or native .cgref bytes. Pin the exact declared input byte length and optional BLAKE3; build/options apply only to FASTA. Returns an opaque import handle, never a filesystem path. Data are streamed into managed scratchfiles outside the small-document store; no giant Molecule or whole-file RAM buffer. Finish requires full private semantic readback before publication."
    )]
    pub(crate) async fn begin_genome_import(
        &self,
        Parameters(args): Parameters<GenomeBeginQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();let id=g.native()?.begin(args.format,args.expected_bytes,args.expected_hash,args.build,args.options).map_err(bad)?;output(json!({"import_id":id,"expected_bytes":args.expected_bytes,"received":0,"max_chunk_bytes":1048576,"ready":false,"authority":"Unpublished private streaming import; labels or supplied hashes do not establish genome completeness."}))}).await
    }
    #[tool(
        description = "Append one exact-offset base64 chunk (at most 1 MiB decoded) to a live owned large-reference import. Offset/size/base64 checks precede writes. Invalid stream/parser/hash evidence aborts the private import; no partial dataset is published. No caller-selected filesystem paths."
    )]
    pub(crate) async fn append_genome_import(
        &self,
        Parameters(args): Parameters<GenomeAppendQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > 1048576usize.div_ceil(3) * 4 {
            return Err(bad("genome chunk exceeds 1 MiB decoded"));
        }
        let bytes = STANDARD.decode(&args.data_base64).map_err(bad)?;
        if bytes.is_empty() || bytes.len() > 1048576 {
            return Err(bad("genome chunk requires 1..1048576 bytes"));
        }
        let server = self.clone();
        self.compute(move || {
            let mut g = server.genomes.lock().unwrap();
            let received = g
                .existing()?
                .append(&args.import_id, args.offset, &bytes)
                .map_err(bad)?;
            output(json!({"import_id":args.import_id,"received":received,"ready":false}))
        })
        .await
    }
    #[tool(
        description = "Finalize framing and cooperatively verify at most max_blocks (1..64) canonical packed blocks per call. Repeat while done=false. Every block and full canonical contig hash must pass a private readback cursor before the dataset handle is published. Serialized receipts/progress cannot restore authority. No one-shot whole-genome RAM allocation or unbounded verification job."
    )]
    pub(crate) async fn finish_genome_import(
        &self,
        Parameters(args): Parameters<GenomeFinishQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_blocks == 0 || args.max_blocks > 64 {
            return Err(bad("finish max_blocks requires 1..64"));
        }
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();let p=g.existing()?.finish(&args.import_id,args.max_blocks).map_err(bad)?;output(json!({"done":p.done,"import_id":p.import_id,"source_bytes":p.source_bytes,"blocks_verified":p.blocks_verified,"blocks_total":p.blocks_total,"dataset":p.dataset.as_ref().map(native_summary),"authority":"Private complete readback establishes stored DNA integrity only; declared build labels do not attest a complete public reference genome."}))}).await
    }
    #[tool(
        description = "List live large-reference dataset summaries or inspect a bounded contig metadata page for an exact opaque dataset/manifest/content binding. Build/accession/expected-component declarations remain visible and are not independently verified whole-genome coverage. Returns no source DNA or arbitrary paths."
    )]
    pub(crate) async fn inspect_genome(
        &self,
        Parameters(args): Parameters<GenomeInspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_contigs == 0 || args.max_contigs > 100 {
            return Err(bad("contig pages require 1..100 rows"));
        }
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();if let Some(r)=args.genome{r.validate().map_err(bad)?;let d=g.existing()?.inspect(&r.id,&r.manifest_hash,&r.content_revision).map_err(bad)?;let m=&d.reference.manifest;if args.contig_start>m.contigs.len(){return Err(bad("contig page offset exceeds count"));}let end=args.contig_start.saturating_add(args.max_contigs).min(m.contigs.len());let rows=m.contigs[args.contig_start..end].iter().map(|c|json!({"id":c.id,"source_id":c.source_id,"accession":c.accession,"description":c.description.chars().take(256).collect::<String>(),"length":c.length,"topology":c.topology,"sequence_hash":c.sequence_hash,"block_count":c.blocks.len()})).collect::<Vec<_>>();output(json!({"dataset":native_summary(&d),"contig_start":args.contig_start,"next_contig_start":end,"done":end==m.contigs.len(),"contigs":rows,"available_profiles":[collagene_core::sequence::crispr::spcas9(),collagene_core::sequence::crispr::cas12a()]}))}else{let rows=g.native.as_ref().map(|n|n.list().iter().map(native_summary).collect::<Vec<_>>()).unwrap_or_default();output(json!({"datasets":rows,"storage":"Managed session scratchfiles, separate from the 64 MiB small-document budget."}))}}).await
    }
    #[tool(
        description = "Read a bounded raw .cgref export byte chunk (at most 1 MiB) from an exact live owned dataset. Concatenate decoded chunks in offset order; no whole-file/base64 artifact is materialized. The native framing and hash manifest are retained. Reimport validates all encoded/canonical sequence hashes before publication; no caller-chosen path is accepted."
    )]
    pub(crate) async fn read_genome_export(
        &self,
        Parameters(args): Parameters<GenomeExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_bytes == 0 || args.max_bytes > 1048576 {
            return Err(bad("export chunks require 1..1048576 bytes"));
        }
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();let bytes=g.existing()?.read_export(&args.genome.id,&args.genome.manifest_hash,&args.genome.content_revision,args.offset,args.max_bytes).map_err(bad)?;let next=args.offset.checked_add(bytes.len()as u64).ok_or_else(||bad("export offset overflow"))?;output(json!({"genome":args.genome,"offset":args.offset,"next_offset":next,"done":bytes.len()<args.max_bytes,"data_base64":STANDARD.encode(&bytes),"encoding":"base64","mime_type":"application/vnd.collagene.cgref","bytes":bytes.len()}))}).await
    }
    #[tool(
        description = "Read at most 16384 verified local reference bases under exact global contig coordinates, role and dataset hashes, without converting a large genome into Molecule. Required packed blocks are read and verified locally. This interval does not establish fresh whole-contig or whole-manifest canonical hash coverage."
    )]
    pub(crate) async fn read_genome_interval(
        &self,
        Parameters(args): Parameters<GenomeReadQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args
            .options
            .interval
            .end
            .checked_sub(args.options.interval.start)
            .is_none_or(|n| n > 16384)
        {
            return Err(bad("MCP sequence windows are limited to 16384 bases"));
        }
        let server = self.clone();
        self.compute(move || {
            let mut g = server.genomes.lock().unwrap();
            let n = g.existing()?;
            let d = n
                .inspect(
                    &args.genome.id,
                    &args.genome.manifest_hash,
                    &args.genome.content_revision,
                )
                .map_err(bad)?;
            let req = collagene_core::tracks::crispr::read_interval_required_blocks(
                &d.reference.manifest,
                &args.options,
            )
            .map_err(bad)?;
            let blocks = n
                .read_blocks(
                    &args.genome.id,
                    &args.genome.manifest_hash,
                    &args.genome.content_revision,
                    &req,
                )
                .map_err(bad)?;
            let report = collagene_core::tracks::crispr::read_interval(
                &d.reference.manifest,
                &args.options,
                &blocks,
            )
            .map_err(bad)?;
            report
                .validate_bindings_for(&d.reference.manifest, &args.options)
                .map_err(bad)?;
            output(json!({"genome":args.genome,"read":report}))
        })
        .await
    }
    #[tool(
        description = "Design a bounded page of source-bound CRISPR guide/PAM candidates over global interval owner coordinates, including verified right halo/circular seam context. The exact source role, manifest/content hashes and nuclease/filter options remain bound. Returns full immutable JSON evidence plus at most100 candidates per page. Local integrity is not whole-contig or build completeness, editing activity or biological safety."
    )]
    pub(crate) async fn design_genome_interval_guides(
        &self,
        Parameters(args): Parameters<GenomeGuidesQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_candidates == 0 || args.max_candidates > 100 {
            return Err(bad("guide pages require1..100 rows"));
        }
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();let n=g.existing()?;let d=n.inspect(&args.genome.id,&args.genome.manifest_hash,&args.genome.content_revision).map_err(bad)?;let req=collagene_core::tracks::crispr::interval_required_blocks(&d.reference.manifest,&args.options).map_err(bad)?;let blocks=n.read_blocks(&args.genome.id,&args.genome.manifest_hash,&args.genome.content_revision,&req).map_err(bad)?;let report=collagene_core::tracks::crispr::design_interval(&d.reference.manifest,&args.options,&blocks).map_err(bad)?;report.validate_bindings_for(&d.reference.manifest,&args.options).map_err(bad)?;if args.candidate_start>report.candidates.len(){return Err(bad("candidate offset exceeds count"));}let end=args.candidate_start.saturating_add(args.max_candidates).min(report.candidates.len());drop(g);let artifact_id=server.session.text_artifact(bounded_json(&report)?,"application/json","genome-guides.json")?;output(json!({"genome":args.genome,"manifest_hash":report.manifest_hash,"content_revision":report.content_revision,"engine_version":report.engine_version,"options":report.options,"options_revision":report.options_revision,"complete":report.complete,"comparisons":report.comparisons,"limits":report.limits,"warnings":report.warnings,"candidate_count":report.candidates.len(),"candidate_start":args.candidate_start,"next_candidate_start":end,"done":end==report.candidates.len(),"candidates":&report.candidates[args.candidate_start..end],"artifact_id":artifact_id}))}).await
    }
    #[tool(
        description = "Start a private cooperative whole-manifest CRISPR Hamming/PAM search against one exact live file-backed dataset and explicit reference role. Each next owner block plus verified right halo is processed only by advance. All canonical contig hashes must be verified before exhaustive-for-manifest/absence evidence; missing declared components prevent complete-build flags. Labels/checkpoints cannot grant authority. No external aligner, physical activity or whole-public-build attestation."
    )]
    pub(crate) async fn start_genome_off_target_search(
        &self,
        Parameters(args): Parameters<GenomeSearchQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();if g.searches.len()>=4{return Err(bad("at most four live genome searches"));}let d=g.existing()?.inspect(&args.genome.id,&args.genome.manifest_hash,&args.genome.content_revision).map_err(bad)?;let state=collagene_core::tracks::crispr::GenomeSearchSession::begin(&d.reference.manifest,&args.query,&args.options).map_err(bad)?;g.serial=g.serial.checked_add(1).ok_or_else(||bad("search serial overflow"))?;let id=blake3::hash(format!("{}:{}:{}",args.genome.id,args.genome.manifest_hash,g.serial).as_bytes()).to_hex().to_string();let progress=state.progress();g.searches.insert(id.clone(),NativeSearch{genome:args.genome.clone(),state});output(json!({"job_id":id,"genome":args.genome,"progress":progress,"authority":"Private live cursor only; serialized progress cannot resume or prove fresh coverage."}))}).await
    }
    #[tool(
        description = "Advance a live whole-manifest search by at most max_blocks1..8 owner steps. The backend reads exact required owner/halo blocks from its managed file; clients cannot submit a forged continuation or skip components. Each call is bounded and cooperative. Coverage and source-integrity flags remain incomplete until private ordered canonical-hash replay reaches EOF. Errors invalidate the private search instead of returning false absence."
    )]
    pub(crate) async fn advance_genome_off_target_search(
        &self,
        Parameters(args): Parameters<GenomeAdvanceQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_blocks == 0 || args.max_blocks > 8 {
            return Err(bad("search advance requires1..8 owner blocks"));
        }
        let server = self.clone();
        self.compute(move || {
            let mut g = server.genomes.lock().unwrap();
            let mut search = g
                .searches
                .remove(&args.job_id)
                .ok_or_else(|| bad("unknown live genome search"))?;
            for _ in 0..args.max_blocks {
                let p = search.state.progress();
                if p.finished {
                    break;
                }
                let blocks = g
                    .existing()?
                    .read_blocks(
                        &search.genome.id,
                        &search.genome.manifest_hash,
                        &search.genome.content_revision,
                        &p.next_blocks,
                    )
                    .map_err(bad)?;
                search.state.advance(&blocks).map_err(bad)?;
            }
            let progress = search.state.progress();
            let genome = search.genome.clone();
            g.searches.insert(args.job_id.clone(), search);
            output(json!({"job_id":args.job_id,"genome":genome,"progress":progress}))
        })
        .await
    }
    #[tool(
        description = "Read a bounded exact hit page from a live privately computed whole-manifest search. Every hit retains global contig coordinates, orientation, ambiguity bounds and exact intended-role/manifest/coordinate classification. Partial/cancelled/missing-component evidence cannot imply whole-genome absence. Imported JSON or checkpoint fields cannot populate this private state."
    )]
    pub(crate) async fn read_genome_off_target_hits(
        &self,
        Parameters(args): Parameters<GenomePageQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > 100 {
            return Err(bad("hit pages require1..100 rows"));
        }
        let server = self.clone();
        self.compute(move||{let g=server.genomes.lock().unwrap();let s=g.searches.get(&args.job_id).ok_or_else(||bad("unknown live genome search"))?;let hits=s.state.page(args.offset,args.limit).map_err(bad)?;let progress=s.state.progress();let next=args.offset.checked_add(hits.len()).ok_or_else(||bad("hit offset overflow"))?;output(json!({"job_id":args.job_id,"genome":s.genome,"progress":progress,"offset":args.offset,"next_offset":next,"done":next==s.state.progress().hits_retained,"hits":hits}))}).await
    }
    #[tool(
        description = "Dispose a session-owned genome import/search handle, or release a dataset when its exact genome binding is supplied. Releasing a dataset cancels related searches and deletes its managed scratchfile; arbitrary paths are never accepted. Private parser/search cursors cannot be restored from serialized progress. No persistent application workspace changes."
    )]
    pub(crate) async fn cancel_genome_job(
        &self,
        Parameters(args): Parameters<GenomeCancelQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let mut g = server.genomes.lock().unwrap();
            if let Some(r) = args.genome {
                if r.id != args.id {
                    return Err(bad("dataset disposal binding mismatch"));
                }
                g.existing()?
                    .delete_dataset(&r.id, &r.manifest_hash, &r.content_revision)
                    .map_err(bad)?;
                g.searches.retain(|_, s| s.genome.id != r.id);
                return output(json!({"id":args.id,"removed":true,"kind":"dataset"}));
            }
            if let Some(mut search) = g.searches.remove(&args.id) {
                let progress = search.state.cancel();
                return output(
                    json!({"id":args.id,"removed":true,"kind":"search","progress":progress}),
                );
            }
            let removed = g
                .native
                .as_mut()
                .is_some_and(|n| n.cancel_import(&args.id).is_ok());
            output(json!({"id":args.id,"removed":removed,"kind":"import"}))
        })
        .await
    }
}

pub(crate) fn genome_chunk_default() -> usize {
    65536
}
pub(crate) fn genome_finish_default() -> usize {
    8
}
pub(crate) fn genome_page_default() -> usize {
    25
}
pub(crate) fn genome_advance_default() -> usize {
    1
}
