use super::*;
use crate::{Alphabet, Error};
use serde::de::value::{Error as ValueError, U8Deserializer};
use std::{
    cell::Cell,
    panic::{AssertUnwindSafe, catch_unwind},
    sync::{Arc, Barrier},
};

struct ClaimedSequence<'a> {
    length: usize,
    index: usize,
    known: bool,
    polled: &'a Cell<usize>,
}
impl<'de> SeqAccess<'de> for ClaimedSequence<'_> {
    type Error = ValueError;
    fn size_hint(&self) -> Option<usize> {
        self.known.then_some(self.length - self.index)
    }
    fn next_element_seed<S: DeserializeSeed<'de>>(
        &mut self,
        seed: S,
    ) -> Result<Option<S::Value>, ValueError> {
        if self.index == self.length {
            return Ok(None);
        }
        self.index += 1;
        self.polled.set(self.polled.get() + 1);
        seed.deserialize(U8Deserializer::<ValueError>::new(7))
            .map(Some)
    }
}
impl<'de> Deserializer<'de> for ClaimedSequence<'_> {
    type Error = ValueError;
    fn deserialize_any<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, ValueError> {
        visitor.visit_seq(self)
    }
    serde::forward_to_deserialize_any! {
        bool i8 i16 i32 i64 u8 u16 u32 u64 f32 f64 char str string bytes byte_buf option
        unit unit_struct newtype_struct seq tuple tuple_struct map struct enum identifier ignored_any
    }
}
fn claimed<'a>(length: usize, known: bool, polled: &'a Cell<usize>) -> ClaimedSequence<'a> {
    ClaimedSequence {
        length,
        index: 0,
        known,
        polled,
    }
}
fn molecule(name: &str, entries: usize) -> Molecule {
    let mut m = Molecule::dna(name, "ACGT", Topology::Linear).unwrap();
    m.metadata = (0..entries)
        .map(|index| (format!("key-{index}"), "v".into()))
        .collect();
    m
}
fn feature() -> Feature {
    Feature {
        kind: "gene".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("example".into()))],
        needs_review: false,
    }
}

#[test]
fn declared_collection_limits_fail_before_reading_elements_or_reserving_the_claim() {
    let polled = Cell::new(0);
    assert!(bounded_vec::<_, u8, 1000>(claimed(usize::MAX, true, &polled)).is_err());
    assert_eq!(polled.get(), 0);
    let mut budget = Budget::new();
    assert!(
        Features(&mut budget)
            .deserialize(claimed(100_001, true, &polled))
            .is_err()
    );
    assert_eq!(polled.get(), 0);
    assert!(
        Qualifiers(&mut budget)
            .deserialize(claimed(QUALIFIERS + 1, true, &polled))
            .is_err()
    );
    assert_eq!(polled.get(), 0);
    assert!(
        Parts {
            budget: &mut budget,
            depth: 0
        }
        .deserialize(claimed(PARTS + 1, true, &polled))
        .is_err()
    );
    assert_eq!(polled.get(), 0);
    budget.0.qualifier_count(999_999).unwrap();
    assert!(
        Qualifiers(&mut budget)
            .deserialize(claimed(2, true, &polled))
            .is_err()
    );
    assert_eq!(
        polled.get(),
        0,
        "aggregate qualifier allowance must constrain a later feature before its first tuple is read"
    );
}

#[test]
fn unknown_collection_lengths_still_stop_before_constructing_the_extra_value() {
    let polled = Cell::new(0);
    let error = bounded_vec::<_, u8, 3>(claimed(100, false, &polled)).unwrap_err();
    assert!(error.to_string().contains("resource budget"));
    assert_eq!(
        polled.get(),
        4,
        "only three values may be constructed; the fourth seed refuses before deserializing its u8"
    );
    polled.set(0);
    assert_eq!(
        bounded_vec::<_, u8, 3>(claimed(3, false, &polled)).unwrap(),
        vec![7, 7, 7]
    );
    assert_eq!(polled.get(), 3);
}

#[test]
fn current_layout_roundtrips_full_native_feature_capacity_and_optional_source_text() {
    let mut m = molecule("native-reference", 0);
    m.features = vec![feature(); 100_000];
    let bytes = postcard::to_allocvec(&m).unwrap();
    let decoded: Molecule =
        with_native_frame_budget(128 * 1024 * 1024, || postcard::from_bytes(&bytes))
            .value
            .unwrap();
    assert_eq!(decoded, m);
    assert_eq!(
        postcard::to_allocvec(&decoded).unwrap(),
        bytes,
        "the current encoding is unchanged"
    );
    let mut source = molecule("original-record", 0);
    source.original_genbank = Some("x".repeat(9 * 1024 * 1024));
    let bytes = postcard::to_allocvec(&source).unwrap();
    assert_eq!(
        postcard::from_bytes::<Molecule>(&bytes).unwrap(),
        source,
        "original records keep their separate32MiB allowance"
    );
    let json = serde_json::to_string(&feature()).unwrap();
    assert_eq!(serde_json::from_str::<Feature>(&json).unwrap(), feature());
    let minimal = r#"{"name":"reference","sequence":{"alphabet":"Dna","bases":"ACGT"},"topology":"Linear","features":[],"metadata":{}}"#;
    assert_eq!(
        serde_json::from_str::<Molecule>(minimal)
            .unwrap()
            .original_genbank,
        None
    );
}

#[test]
fn locations_bound_fanout_depth_and_standalone_qualifier_collections() {
    let location = Location::Join(vec![Location::range(0, 4); PARTS]);
    let bytes = postcard::to_allocvec(&location).unwrap();
    assert_eq!(postcard::from_bytes::<Location>(&bytes).unwrap(), location);
    let too_many = Location::Order(vec![Location::range(0, 4); PARTS + 1]);
    assert!(postcard::from_bytes::<Location>(&postcard::to_allocvec(&too_many).unwrap()).is_err());
    let mut nested = Location::range(0, 4);
    for _ in 0..31 {
        nested = Location::Complement(Box::new(nested));
    }
    assert!(postcard::from_bytes::<Location>(&postcard::to_allocvec(&nested).unwrap()).is_ok());
    nested = Location::Complement(Box::new(nested));
    assert!(postcard::from_bytes::<Location>(&postcard::to_allocvec(&nested).unwrap()).is_err());
    let mut f = feature();
    f.qualifiers = vec![("key".into(), None); QUALIFIERS + 1];
    assert!(postcard::from_bytes::<Feature>(&postcard::to_allocvec(&f).unwrap()).is_err());
}

#[test]
fn individually_valid_molecules_share_frame_qualifier_budget_and_error_restores_it() {
    let mut row = molecule("valid-owner", 0);
    row.features = vec![feature()];
    let bytes = postcard::to_allocvec(&vec![row.clone(), row]).unwrap();
    let result = with_native_frame_budget(128 * 1024 * 1024, || {
        // Leave room for only one legal molecule's qualifier, without allocating
        // an attack-sized fixture. The second molecule must share this budget.
        frame(|b| b.qualifier_count(999_999)).unwrap();
        postcard::from_bytes::<Vec<Molecule>>(&bytes)
    })
    .value;
    assert!(result.is_err());
    assert!(
        with_native_frame_budget(128 * 1024 * 1024, || postcard::from_bytes::<Vec<Molecule>>(
            &bytes
        ))
        .value
        .is_ok()
    );
}

#[test]
fn workspace_scale_metadata_and_multiple_retained_sources_fit_frame_profile() {
    let rows = (0..1000)
        .map(|_| molecule("workspace-reference", 5))
        .collect::<Vec<_>>();
    let bytes = postcard::to_allocvec(&rows).unwrap();
    assert_eq!(
        with_native_frame_budget(128 * 1024 * 1024, || postcard::from_bytes::<Vec<Molecule>>(
            &bytes
        ))
        .value
        .unwrap(),
        rows
    );
    with_native_frame_budget(128 * 1024 * 1024, || {
        // Multiple owners can legitimately retain more than32MiB collectively.
        // Reuse a small borrowed slice when charging to avoid a huge fixture.
        let chunk = "x".repeat(1024 * 1024);
        for _ in 0..33 {
            frame(|b| b.original(&chunk)).unwrap();
        }
    })
    .value;
    with_native_frame_budget(3, || {
        frame(|b| b.text("ab", 1024)).unwrap();
        assert!(frame(|b| b.original("cd")).is_err());
    })
    .value;
}

#[test]
fn retained_evidence_keys_preserve_large_attachments_without_generic_escape() {
    let keys = [
        "trace_original_base64".into(),
        "source:snapgene:original_base64".into(),
        "legacy:original_base64".into(),
        "fastq:original_base64".into(),
        "fastq:quality_base64".into(),
        format!("annotation_import:{}:original_base64", "a".repeat(64)),
    ];
    let value = "x".repeat(1024 * 1024 + 1);
    for key in keys {
        let mut m = molecule("retained-source", 0);
        m.metadata.insert(key, value.clone());
        crate::bounds::annotation_edit_preflight(&[&m]).unwrap();
        let bytes = postcard::to_allocvec(&m).unwrap();
        assert_eq!(postcard::from_bytes::<Molecule>(&bytes).unwrap(), m);
    }
    for key in [
        "arbitrary_base64".into(),
        format!("annotation_import:{}:original_base64", "G".repeat(64)),
    ] {
        let mut m = molecule("ordinary-text", 0);
        m.metadata.insert(key, value.clone());
        assert!(postcard::from_bytes::<Molecule>(&postcard::to_allocvec(&m).unwrap()).is_err());
    }
}

#[test]
fn nested_frame_guards_share_counts_and_restore_after_panic() {
    with_native_frame_budget(128 * 1024 * 1024, || {
        frame(|b| b.metadata_count(1_000_000)).unwrap();
        let result =
            with_native_frame_budget(128 * 1024 * 1024, || frame(|b| b.metadata_count(1))).value;
        assert!(matches!(result, Err(Error::TooLarge)));
    })
    .value;
    let result = catch_unwind(AssertUnwindSafe(|| {
        with_native_frame_budget(128 * 1024 * 1024, || {
            frame(|b| b.metadata_count(1_000_000)).unwrap();
            panic!("deliberate decoder interruption");
        })
        .value
    }));
    assert!(result.is_err());
    with_native_frame_budget(128 * 1024 * 1024, || frame(|b| b.metadata_count(1_000_000)))
        .value
        .unwrap();
}

#[test]
fn simultaneous_native_frames_have_independent_thread_budgets() {
    let barrier = Arc::new(Barrier::new(2));
    let threads = (0..2)
        .map(|_| {
            let barrier = barrier.clone();
            std::thread::spawn(move || {
                with_native_frame_budget(128 * 1024 * 1024, || {
                    frame(|b| b.metadata_count(1_000_000)).unwrap();
                    barrier.wait();
                    assert!(frame(|b| b.metadata_count(1)).is_err());
                })
                .value
            })
        })
        .collect::<Vec<_>>();
    for thread in threads {
        thread.join().unwrap();
    }
    assert!(
        frame(|b| b.metadata_count(1_000_000)).is_ok(),
        "outside a frame the thread-local scope is absent"
    );
}

#[test]
fn native_text_and_metadata_fields_reject_excess_and_duplicates_without_changing_alphabet() {
    let mut m = molecule("source", 4097);
    assert!(postcard::from_bytes::<Molecule>(&postcard::to_allocvec(&m).unwrap()).is_err());
    m = molecule("source", 0);
    m.name = "x".repeat(1025);
    assert!(postcard::from_bytes::<Molecule>(&postcard::to_allocvec(&m).unwrap()).is_err());
    let duplicate = r#"{"name":"reference","sequence":{"alphabet":"Dna","bases":"ACGT"},"topology":"Linear","features":[],"metadata":{"k":"one","k":"two"}}"#;
    assert!(serde_json::from_str::<Molecule>(duplicate).is_err());
    let mut m = Molecule::dna("source", "ACGT", Topology::Linear).unwrap();
    m.sequence = Sequence::new(Alphabet::Rna, "ACGU").unwrap();
    assert_eq!(
        postcard::from_bytes::<Molecule>(&postcard::to_allocvec(&m).unwrap()).unwrap(),
        m
    );
}

#[test]
fn literal_current_postcard_records_keep_the_original_field_and_variant_order() {
    // Current Molecule: name="x", sequence=(Dna,"A"), Linear, no features,
    // no metadata, no original source. These bytes are specified independently
    // of the derived serializer and contain no decoder budget wrapper fields.
    let native = [1, b'x', 0, 1, b'A', 0, 0, 0, 0];
    let expected = Molecule::dna("x", "A", Topology::Linear).unwrap();
    assert_eq!(postcard::from_bytes::<Molecule>(&native).unwrap(), expected);
    assert_eq!(postcard::to_allocvec(&expected).unwrap(), native);
    let range = [0, 0, 1, 0, 0];
    assert_eq!(
        postcard::from_bytes::<Location>(&range).unwrap(),
        Location::range(0, 1)
    );
    let complement = [1, 0, 0, 1, 0, 0];
    assert_eq!(
        postcard::from_bytes::<Location>(&complement).unwrap(),
        Location::Complement(Box::new(Location::range(0, 1)))
    );
    let join = [2, 1, 0, 0, 1, 0, 0];
    assert_eq!(
        postcard::from_bytes::<Location>(&join).unwrap(),
        Location::Join(vec![Location::range(0, 1)])
    );
}

#[test]
fn bounded_binary_octets_use_literal_current_bytes_and_canonical_json_arrays() {
    #[derive(serde::Serialize, serde::Deserialize)]
    struct Packet {
        #[serde(
            serialize_with = "super::serialize_bytes",
            deserialize_with = "super::bounded_bytes::<_, 5>"
        )]
        bytes: Vec<u8>,
    }
    let packet = Packet {
        bytes: vec![0, 127, 128, 255],
    };
    let encoded = postcard::to_allocvec(&packet).unwrap();
    assert_eq!(encoded, [4, 0, 127, 128, 255]);
    assert_eq!(encoded, postcard::to_allocvec(&packet.bytes).unwrap());
    assert_eq!(
        postcard::from_bytes::<Packet>(&encoded).unwrap().bytes,
        packet.bytes
    );
    let json = serde_json::to_vec(&packet).unwrap();
    assert_eq!(json, br#"{"bytes":[0,127,128,255]}"#);
    assert_eq!(
        serde_json::from_slice::<Packet>(&json).unwrap().bytes,
        packet.bytes
    );
    assert!(serde_json::from_slice::<Packet>(br#"{"bytes":"abc"}"#).is_err());
    assert!(postcard::from_bytes::<Packet>(&[6, 1, 2, 3, 4, 5, 6]).is_err());
    assert!(serde_json::from_slice::<Packet>(br#"{"bytes":[1,2,3,4,5,{"#).is_err());
}

#[test]
fn frame_outcomes_distinguish_resource_rejection_from_malformed_data_and_reset() {
    let mut m = molecule("large-ordinary-source", 0);
    m.metadata
        .insert("source-evidence".into(), "x".repeat(1024 * 1024 + 1));
    let bytes = postcard::to_allocvec(&m).unwrap();
    let outcome = with_native_frame_budget(32 * 1024 * 1024, || {
        postcard::from_bytes::<Molecule>(&bytes)
    });
    assert!(outcome.value.is_err());
    assert!(outcome.resource_exhausted);
    let malformed =
        with_native_frame_budget(32 * 1024 * 1024, || postcard::from_bytes::<Molecule>(&[0]));
    assert!(malformed.value.is_err());
    assert!(!malformed.resource_exhausted);
    let oversized_sequence =
        postcard::to_allocvec(&(Alphabet::Dna, "A".repeat(crate::MAX_SEQUENCE_LEN + 1))).unwrap();
    let sequence = with_native_frame_budget(32 * 1024 * 1024, || {
        postcard::from_bytes::<Sequence>(&oversized_sequence)
    });
    assert!(sequence.value.is_err());
    assert!(sequence.resource_exhausted);
    let invalid_sequence = postcard::to_allocvec(&(Alphabet::Dna, "!invalid")).unwrap();
    let sequence = with_native_frame_budget(32 * 1024 * 1024, || {
        postcard::from_bytes::<Sequence>(&invalid_sequence)
    });
    assert!(sequence.value.is_err());
    assert!(!sequence.resource_exhausted);
    let nested = with_native_frame_budget(32 * 1024 * 1024, || {
        let polled = Cell::new(0);
        let inner = with_native_frame_budget(32 * 1024 * 1024, || {
            bounded_vec::<_, u8, 1>(claimed(2, true, &polled))
        });
        assert!(inner.resource_exhausted);
        assert!(inner.value.is_err());
        assert_eq!(polled.get(), 0);
    });
    assert!(nested.resource_exhausted);
    assert!(!with_native_frame_budget(32 * 1024 * 1024, || ()).resource_exhausted);
    let panic = catch_unwind(AssertUnwindSafe(|| {
        with_native_frame_budget(1, || {
            assert!(frame(|b| b.text("xx", 2)).is_err());
            panic!("resource rejection followed by unwind");
        })
    }));
    assert!(panic.is_err());
    assert!(!with_native_frame_budget(32 * 1024 * 1024, || ()).resource_exhausted);
}

#[test]
fn binary_and_unknown_length_collection_admission_mark_resource_outcomes() {
    #[derive(serde::Deserialize)]
    struct Bytes {
        #[serde(deserialize_with = "super::bounded_bytes::<_, 1>")]
        _bytes: Vec<u8>,
    }
    let bytes = with_native_frame_budget(32, || postcard::from_bytes::<Bytes>(&[2, 1, 2]));
    assert!(bytes.resource_exhausted);
    assert!(bytes.value.is_err());
    let polled = Cell::new(0);
    let values =
        with_native_frame_budget(32, || bounded_vec::<_, u8, 1>(claimed(2, false, &polled)));
    assert!(values.resource_exhausted);
    assert!(values.value.is_err());
    let barrier = Arc::new(Barrier::new(2));
    let other = barrier.clone();
    let flagged = std::thread::spawn(move || {
        with_native_frame_budget(1, || {
            assert!(frame(|b| b.text("xx", 2)).is_err());
            other.wait();
        })
        .resource_exhausted
    });
    let clean = with_native_frame_budget(1, || barrier.wait());
    assert!(!clean.resource_exhausted);
    assert!(flagged.join().unwrap());
}
