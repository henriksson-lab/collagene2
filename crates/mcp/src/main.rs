use rmcp::ServiceExt;
use std::{ffi::OsString, path::PathBuf, process::ExitCode};
const HELP: &str = r#"Usage: collagene-mcp [--storage-root PATH]

Session-only MCP server over stdio. Standard output is reserved for MCP.

Options:
  --storage-root PATH  Existing folder for managed genome data; creates a private
                       per-session child. Default: private child of OS temp folder.
  -h, --help           Show this help on stderr.
  -V, --version        Show version on stderr.

Examples:
  collagene-mcp --storage-root /path/to/existing/CollageneData
  collagene-mcp --storage-root "C:\Users\you\AppData\Local\CollageneData"

No arbitrary file paths, read/write tools or filesystem HTTP API are enabled.
Ordinary documents remain memory-only. This is managed application storage,
not an OS sandbox against other processes running as the same user.
"#;
enum Action {
    Serve(Option<PathBuf>),
    Help,
    Version,
}
fn parse(args: impl IntoIterator<Item = OsString>) -> Result<Action, String> {
    let mut args = args.into_iter();
    let mut root = None;
    while let Some(arg) = args.next() {
        if arg == "--help" || arg == "-h" {
            return Ok(Action::Help);
        }
        if arg == "--version" || arg == "-V" {
            return Ok(Action::Version);
        }
        let value = if arg == "--storage-root" {
            args.next()
                .filter(|v| !v.is_empty() && !v.to_string_lossy().starts_with("--"))
                .ok_or_else(|| "--storage-root requires an existing directory PATH".to_string())?
        } else if let Some(value) = arg.to_str().and_then(|s| s.strip_prefix("--storage-root=")) {
            if value.is_empty() {
                return Err("--storage-root requires a nonempty PATH".into());
            }
            OsString::from(value)
        } else {
            return Err(format!(
                "unrecognized argument: {} (use --help)",
                arg.to_string_lossy()
            ));
        };
        if root.replace(PathBuf::from(value)).is_some() {
            return Err("--storage-root may be supplied only once".into());
        }
    }
    Ok(Action::Serve(root))
}
fn diagnostic(category: &str, error: impl std::fmt::Display) {
    eprintln!(
        "{}",
        serde_json::json!({"error":category,"message":error.to_string()})
    );
}
#[tokio::main]
async fn main() -> ExitCode {
    // stdout belongs exclusively to MCP. Never launch clients/scientific executables.
    let root = match parse(std::env::args_os().skip(1)) {
        Ok(Action::Help) => {
            eprint!("{HELP}");
            return ExitCode::SUCCESS;
        }
        Ok(Action::Version) => {
            eprintln!("collagene-mcp {}", env!("CARGO_PKG_VERSION"));
            return ExitCode::SUCCESS;
        }
        Ok(Action::Serve(root)) => root,
        Err(error) => {
            diagnostic("invalid_arguments", error);
            return ExitCode::from(2);
        }
    };
    let server = match root {
        Some(root) => match collagene_mcp::CollageneServer::with_storage_root(&root) {
            Ok(server) => server,
            Err(error) => {
                diagnostic("storage_configuration", error);
                return ExitCode::from(2);
            }
        },
        None => collagene_mcp::CollageneServer::default(),
    };
    match server.serve(rmcp::transport::stdio()).await {
        Ok(service) => match service.waiting().await {
            Ok(_) => ExitCode::SUCCESS,
            Err(error) => {
                diagnostic("mcp_service", error);
                ExitCode::FAILURE
            }
        },
        Err(error) => {
            diagnostic("mcp_service", error);
            ExitCode::FAILURE
        }
    }
}
