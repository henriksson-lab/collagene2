//! New review capabilities cannot travel under a legacy envelope.
use collagene_desktop_sync_protocol as wire;
fn request(command: wire::Command) -> wire::RequestEnvelope {
    wire::RequestEnvelope {
        request_id: 17,
        operation_id: wire::OperationId::new(format!("op-{}", "a".repeat(32))).unwrap(),
        command,
    }
}
#[test]
fn legacy_prefix_requests_get_legacy_replies_without_enabling_new_review_variants() {
    let old = request(wire::Command::Remotes);
    let mut frame = wire::encode_request(&old).unwrap();
    frame[4..6].copy_from_slice(&1u16.to_le_bytes());
    let parsed = wire::decode_request(&frame).unwrap();
    let reply = wire::ResponseEnvelope {
        request_id: 17,
        operation_id: parsed.operation_id.clone(),
        context: None,
        outcome: wire::Outcome::Success(wire::Success::Remotes(
            wire::library::Rows::new(vec![]).unwrap(),
        )),
    };
    let encoded = wire::encode_response_for_frame(&reply, &frame).unwrap();
    assert_eq!(&encoded[4..6], &1u16.to_le_bytes());
    wire::decode_response_for(&encoded, &parsed).unwrap();
}
#[test]
fn review_request_result_and_generic_ack_require_version_two() {
    let id = wire::ReviewId::new(format!("review-{}", "b".repeat(32))).unwrap();
    let q = request(wire::Command::DiscardPushReview { review: id.clone() });
    let mut frame = wire::encode_request(&q).unwrap();
    assert_eq!(&frame[4..6], &2u16.to_le_bytes());
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
