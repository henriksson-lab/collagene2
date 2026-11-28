use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};
static SERIAL: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "collagene-mcp-storage-test-{}-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
async fn call(client: &RunningService<RoleClient, ()>, name: &str, args: Value) -> Value {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
#[test]
fn cli_help_version_and_configuration_errors_keep_stdout_clean() {
    let binary = env!("CARGO_BIN_EXE_collagene-mcp");
    for flag in ["--help", "--version"] {
        let result = Command::new(binary).arg(flag).output().unwrap();
        assert!(result.status.success());
        assert!(result.stdout.is_empty());
        assert!(!result.stderr.is_empty());
    }
    let help = Command::new(binary).arg("--help").output().unwrap();
    assert!(
        String::from_utf8(help.stderr)
            .unwrap()
            .contains("--storage-root PATH")
    );
    let f = Fixture::new();
    for args in [
        vec!["--unknown".into()],
        vec!["--storage-root".into()],
        vec!["--storage-root=".into()],
        vec![
            "--storage-root".into(),
            f.0.join("missing").into_os_string(),
        ],
        vec![
            "--storage-root".into(),
            f.0.clone().into_os_string(),
            "--storage-root".into(),
            f.0.clone().into_os_string(),
        ],
    ] {
        let result = Command::new(binary).args(args).output().unwrap();
        assert_eq!(result.status.code(), Some(2));
        assert!(result.stdout.is_empty());
        let error: Value = serde_json::from_slice(&result.stderr).unwrap();
        assert!(
            error["error"].as_str().unwrap().contains("arguments")
                || error["error"] == "storage_configuration"
        );
    }
    assert_eq!(fs::read_dir(&f.0).unwrap().count(), 0);
}
#[tokio::test]
async fn configured_sdk_session_keeps_opaque_reference_lifecycle_beneath_selected_root() {
    let f = Fixture::new();
    fs::write(f.0.join("user-sentinel"), b"untouched").unwrap();
    let server = CollageneServer::with_storage_root(&f.0).unwrap();
    let sessions = fs::read_dir(&f.0)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.is_dir())
        .collect::<Vec<_>>();
    assert_eq!(sessions.len(), 1);
    let session = sessions[0].clone();
    assert!(fs::read_dir(&session).unwrap().next().is_none());
    let (a, b) = tokio::io::duplex(65536);
    let handle = tokio::spawn(async move {
        server.serve(a).await.unwrap().waiting().await.unwrap();
    });
    let client = ().serve(b).await.unwrap();
    let tools = client.list_all_tools().await.unwrap();
    assert!(!tools.iter().any(|t| matches!(
        t.name.as_ref(),
        "read_file" | "write_file" | "list_files" | "open_path"
    )));
    let rejected=client.call_tool(CallToolRequestParams::new("begin_genome_import".to_owned()).with_arguments(json!({"format":"Fasta","expected_bytes":4,"path":f.0.join("user-sentinel").display().to_string()}).as_object().unwrap().clone())).await;
    assert!(rejected.is_err() || rejected.unwrap().is_error == Some(true));
    assert_eq!(fs::read_dir(&session).unwrap().count(), 0);
    let input = b">chr\nACGTACGT\n";
    let upload = call(
        &client,
        "begin_genome_import",
        json!({"format":"Fasta","expected_bytes":input.len()}),
    )
    .await;
    let id = upload["import_id"].as_str().unwrap();
    assert_eq!(fs::read_dir(&session).unwrap().count(), 1);
    call(
        &client,
        "append_genome_import",
        json!({"import_id":id,"offset":0,"data_base64":STANDARD.encode(input)}),
    )
    .await;
    let finished = call(
        &client,
        "finish_genome_import",
        json!({"import_id":id,"max_blocks":64}),
    )
    .await;
    assert!(finished["done"].as_bool().unwrap());
    let genome = finished["dataset"]["genome"].clone();
    let exported = call(
        &client,
        "read_genome_export",
        json!({"genome":genome,"offset":0,"max_bytes":1024}),
    )
    .await;
    assert!(
        !STANDARD
            .decode(exported["data_base64"].as_str().unwrap())
            .unwrap()
            .is_empty()
    );
    call(
        &client,
        "cancel_genome_job",
        json!({"id":genome["id"],"genome":genome}),
    )
    .await;
    assert_eq!(fs::read_dir(&session).unwrap().count(), 0);
    client.cancel().await.unwrap();
    handle.await.unwrap();
    assert!(!session.exists());
    assert_eq!(fs::read(f.0.join("user-sentinel")).unwrap(), b"untouched");
}
