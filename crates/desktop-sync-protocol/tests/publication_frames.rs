//! Current review requests and acknowledgements reject noncurrent encodings.
use collagene_desktop_sync_protocol as wire;
fn request(command: wire::Command) -> wire::RequestEnvelope {
    wire::RequestEnvelope {
        request_id: 17,
        operation_id: wire::OperationId::new(format!("op-{}", "a".repeat(32))).unwrap(),
        command,
    }
}
#[test]
fn current_review_request_result_and_generic_ack_reject_downgraded_headers() {
    let id = wire::ReviewId::new(format!("review-{}", "b".repeat(32))).unwrap();
    let q = request(wire::Command::DiscardPushReview { review: id.clone() });
    let mut frame = wire::encode_request(&q).unwrap();
    assert_eq!(&frame[4..6], &wire::SCHEMA_VERSION.to_le_bytes());
    frame[4..6].copy_from_slice(&1u16.to_le_bytes());
    assert!(wire::decode_request(&frame).is_err());
    let result = wire::ResponseEnvelope {
        request_id: 17,
        operation_id: q.operation_id.clone(),
        context: None,
        outcome: wire::Outcome::Success(wire::Success::PushReviewDiscarded { review: id }),
    };
    let mut encoded = wire::encode_response(&result).unwrap();
    wire::decode_response_for(&encoded, &q).unwrap();
    encoded[4..6].copy_from_slice(&1u16.to_le_bytes());
    assert!(wire::decode_response_for(&encoded, &q).is_err());
    let error = wire::ResponseEnvelope {
        outcome: wire::Outcome::Failure(wire::Failure::Error {
            code: wire::ErrorCode::Invalid,
            message: wire::library::Text::new("Not applicable").unwrap(),
        }),
        ..result
    };
    let mut bytes = wire::encode_response(&error).unwrap();
    bytes[4..6].copy_from_slice(&1u16.to_le_bytes());
    assert!(wire::decode_response_for(&bytes, &q).is_err());
}
