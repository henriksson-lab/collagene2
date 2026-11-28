use crate::session::MAX_STORE_BYTES;
use rmcp::ErrorData;
use rmcp::model::CallToolResult;
use serde::Serialize;
use serde_json::Value;

pub(crate) fn bad(message: impl ToString) -> ErrorData {
    ErrorData::invalid_params(message.to_string(), None)
}
// Bound temporary JSON materialization before session accounting, including binary arrays.
pub(crate) fn bounded_json<T: Serialize>(value: &T) -> Result<String, ErrorData> {
    bounded_json_limit(value, MAX_STORE_BYTES)
}
pub(crate) fn bounded_json_limit<T: Serialize>(value: &T, cap: usize) -> Result<String, ErrorData> {
    struct Sink {
        bytes: Vec<u8>,
        cap: usize,
    }
    impl std::io::Write for Sink {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if bytes.len() > self.cap.saturating_sub(self.bytes.len()) {
                return Err(std::io::Error::other(
                    "JSON artifact exceeds session budget",
                ));
            }
            self.bytes.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut sink = Sink {
        bytes: Vec::new(),
        cap,
    };
    serde_json::to_writer(&mut sink, value).map_err(bad)?;
    String::from_utf8(sink.bytes).map_err(bad)
}
pub(crate) fn output(value: Value) -> Result<CallToolResult, ErrorData> {
    Ok(CallToolResult::structured(value))
}
