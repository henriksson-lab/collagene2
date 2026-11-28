mod resources;
pub(crate) mod response;

use crate::server::response::bad;
use crate::session::Session;
use crate::session::genome::GenomeResources;
use crate::session::genome_store;
use rmcp::ErrorData;
use rmcp::handler::server::router::tool::ToolRouter;
use rmcp::model::CallToolResult;
use schemars::JsonSchema;
use serde::Deserialize;
use std::sync::Arc;
use std::sync::Mutex;

#[derive(Debug, Clone, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EmptyRequest {}
#[derive(Clone)]
pub struct CollageneServer {
    pub(crate) session: Session,
    pub(crate) tool_router: ToolRouter<Self>,
    pub(crate) compute_gate: Arc<tokio::sync::Semaphore>,
    pub(crate) source_client: Arc<tokio::sync::OnceCell<collagene_sources::SourceClient>>,
    pub(crate) genomes: Arc<Mutex<GenomeResources>>,
}
impl CollageneServer {
    /// Retire this scientific session: queued and future compute admission fails,
    /// while already running native work retains its permit until completion.
    pub fn close_compute(&self) {
        self.compute_gate.close();
    }

    /// After close_compute(), a true result prevents a hosted owner from
    /// releasing its retired session quota or managed storage prematurely.
    pub fn has_pending_compute(&self) -> bool {
        self.compute_gate.available_permits() == 0
    }

    /// Trusted native embedding only: import already selected archive bytes into
    /// this memory-only MCP session. This grants no path or database authority.
    pub fn preload_archive_bytes(&self, bytes: &[u8]) -> Result<CallToolResult, ErrorData> {
        self.session.import_archive_bytes(bytes)
    }
    /// Load a literal native .cg2 document without inventing workspace events,
    /// run receipts, or reusable physical-material authority.
    pub fn preload_document_bytes(&self, bytes: &[u8]) -> Result<CallToolResult, ErrorData> {
        let molecule = collagene_core::protocol::load_document(bytes).map_err(bad)?;
        self.session.insert_documents(vec![molecule], vec![])
    }
    pub(crate) async fn compute<F>(&self, job: F) -> Result<CallToolResult, ErrorData>
    where
        F: FnOnce() -> Result<CallToolResult, ErrorData> + Send + 'static,
    {
        let permit = self
            .compute_gate
            .clone()
            .acquire_owned()
            .await
            .map_err(|e| bad(e.to_string()))?;
        tokio::task::spawn_blocking(move || {
            let _permit = permit;
            job()
        })
        .await
        .map_err(|e| ErrorData::internal_error(e.to_string(), None))?
    }
    /// Configure managed file-backed session storage. This grants no filesystem
    /// tool or arbitrary-path read access; ordinary documents remain memory-only.
    pub fn with_storage_root(path: impl AsRef<std::path::Path>) -> std::io::Result<Self> {
        Self::with_storage_policy(collagene_local_storage::StorageRoot::existing(path)?)
    }
    pub fn with_storage_policy(
        root: collagene_local_storage::StorageRoot,
    ) -> std::io::Result<Self> {
        let native = genome_store::NativeGenomes::in_root(&root)?;
        Ok(Self {
            genomes: Arc::new(Mutex::new(GenomeResources {
                native: Some(native),
                ..Default::default()
            })),
            ..Self::default()
        })
    }
}

impl Default for CollageneServer {
    fn default() -> Self {
        Self {
            session: Session::default(),
            tool_router: Self::analysis_router()
                + Self::assembly_intents_router()
                + Self::assembly_router()
                + Self::catalog_router()
                + Self::cloning_router()
                + Self::genome_router()
                + Self::protocol_router()
                + Self::sequence_router()
                + Self::workspace_router()
                + Self::image_evidence_router(),
            compute_gate: Arc::new(tokio::sync::Semaphore::new(1)),
            source_client: Arc::new(tokio::sync::OnceCell::new()),
            genomes: Arc::default(),
        }
    }
}

/// Serve the embedded server in this process. stdout belongs exclusively to MCP;
/// this helper neither opens a GUI nor starts an external executable.
pub async fn serve_stdio(
    server: CollageneServer,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    use rmcp::ServiceExt;
    server
        .serve(rmcp::transport::stdio())
        .await?
        .waiting()
        .await?;
    Ok(())
}

#[cfg(test)]
mod compute_lifecycle_tests {
    use super::*;

    #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
    async fn retired_compute_remains_pending_after_caller_abort_until_native_task_finishes() {
        let server = CollageneServer::default();
        let worker = server.clone();
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (release_tx, release_rx) = std::sync::mpsc::channel();
        let task = tokio::spawn(async move {
            worker
                .compute(move || {
                    started_tx.send(()).unwrap();
                    release_rx.recv().unwrap();
                    Ok(CallToolResult::success(vec![]))
                })
                .await
        });
        started_rx
            .recv_timeout(std::time::Duration::from_secs(5))
            .unwrap();
        assert!(server.has_pending_compute());
        server.close_compute();
        task.abort();
        assert!(
            server.has_pending_compute(),
            "RPC cancellation cannot release running native work"
        );
        assert!(
            server
                .compute(|| Ok(CallToolResult::success(vec![])))
                .await
                .is_err(),
            "retired engine cannot admit queued/replacement work"
        );
        release_tx.send(()).unwrap();
        tokio::time::timeout(std::time::Duration::from_secs(5), async {
            while server.has_pending_compute() {
                tokio::task::yield_now().await;
            }
        })
        .await
        .unwrap();
        assert!(!server.has_pending_compute());
    }
}
