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
/// Count admission bytes without first allocating a copy of the serialized input.
pub(crate) fn bounded_json_size<T: Serialize>(value: &T, cap: usize) -> Result<usize, ErrorData> {
    struct Count {
        bytes: usize,
        cap: usize,
    }
    impl std::io::Write for Count {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.bytes = self
                .bytes
                .checked_add(bytes.len())
                .filter(|count| *count <= self.cap)
                .ok_or_else(|| std::io::Error::other("JSON admission exceeds session budget"))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut count = Count { bytes: 0, cap };
    serde_json::to_writer(&mut count, value).map_err(bad)?;
    Ok(count.bytes)
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
