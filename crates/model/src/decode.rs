//! Current-format resource admission before constructing native model collections.
//! This is independent of molecular validity and preserves the existing Serde layout.
use crate::{Feature, Location, Molecule, Sequence, Topology, bounds::ModelBudget};
use serde::{
    Deserialize, Deserializer,
    de::{self, DeserializeSeed, EnumAccess, MapAccess, SeqAccess, VariantAccess, Visitor},
};
use std::{cell::RefCell, collections::BTreeMap, fmt, marker::PhantomData, rc::Rc};

const PARTS: usize = 16_384;
const QUALIFIERS: usize = 16_384;

thread_local! {
    static FRAME: RefCell<Option<ModelBudget>> = const { RefCell::new(None) };
}
struct FrameGuard {
    root: bool,
    _not_send: PhantomData<Rc<()>>,
}
impl FrameGuard {
    fn enter(maximum_bytes: usize) -> Self {
        let root = FRAME.with(|slot| {
            let mut slot = slot.borrow_mut();
            if slot.is_some() {
                false
            } else {
                *slot = Some(ModelBudget::native_frame(maximum_bytes));
                true
            }
        });
        Self {
            root,
            _not_send: PhantomData,
        }
    }
}
impl Drop for FrameGuard {
    fn drop(&mut self) {
        if self.root {
            FRAME.with(|slot| *slot.borrow_mut() = None);
        }
    }
}
/// Share aggregate native resource limits across a synchronous decoded frame.
/// Each collection kind is capped at one million; annotation and retained source
/// text share the admitted envelope byte maximum. Per-model limits remain lower.
/// Nested calls share the outer budget. The private non-Send guard never escapes
/// this closure; error or panic restores the thread-local state. No borrow is
/// held while invoking a deserializer, and no guard crosses an await point.
pub fn with_native_frame_budget<T>(maximum_bytes: usize, decode: impl FnOnce() -> T) -> T {
    let _guard = FrameGuard::enter(maximum_bytes);
    decode()
}
fn frame(
    charge: impl FnOnce(&mut ModelBudget) -> Result<(), crate::Error>,
) -> Result<(), crate::Error> {
    FRAME.with(|slot| {
        let mut slot = slot.borrow_mut();
        match slot.as_mut() {
            Some(budget) => charge(budget),
            None => Ok(()),
        }
    })
}
struct Budget(ModelBudget);
impl Budget {
    fn new() -> Self {
        Self(ModelBudget::annotation_edit())
    }
    fn feature(&mut self) -> Result<(), crate::Error> {
        self.0.feature_count(1)?;
        frame(|b| b.feature_count(1))
    }
    fn qualifier(&mut self) -> Result<(), crate::Error> {
        self.0.qualifier_count(1)?;
        frame(|b| b.qualifier_count(1))
    }
    fn node(&mut self) -> Result<(), crate::Error> {
        self.0.node_count(1)?;
        frame(|b| b.node_count(1))
    }
    fn metadata(&mut self) -> Result<(), crate::Error> {
        self.0.metadata_count(1)?;
        frame(|b| b.metadata_count(1))
    }
    fn text(&mut self, s: &str, cap: usize, original: bool) -> Result<(), crate::Error> {
        if original {
            self.0.original(s)?;
            frame(|b| b.original(s))
        } else {
            self.0.text(s, cap)?;
            frame(|b| b.text(s, cap))
        }
    }
}
fn hint<E: de::Error>(hint: Option<usize>, maximum: usize) -> Result<(), E> {
    if hint.is_some_and(|size| size > maximum) {
        return Err(E::custom(
            "native decoded collection exceeds its resource budget",
        ));
    }
    Ok(())
}
struct Element<T> {
    index: usize,
    maximum: usize,
    _type: PhantomData<T>,
}
impl<'de, T: Deserialize<'de>> DeserializeSeed<'de> for Element<T> {
    type Value = T;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<T, D::Error> {
        if self.index >= self.maximum {
            return Err(de::Error::custom(
                "native decoded collection exceeds its resource budget",
            ));
        }
        T::deserialize(d)
    }
}
/// Decode a current native collection without trusting its reserve hint.
/// A known oversized length is refused before reading any element; an unknown
/// hint is bounded before constructing the first element beyond MAX.
pub fn bounded_vec<'de, D, T, const MAX: usize>(d: D) -> Result<Vec<T>, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de>,
{
    struct Values<T, const MAX: usize>(PhantomData<T>);
    impl<'de, T: Deserialize<'de>, const MAX: usize> Visitor<'de> for Values<T, MAX> {
        type Value = Vec<T>;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "at most {MAX} native elements")
        }
        fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
            hint(seq.size_hint(), MAX)?;
            let mut values = Vec::new();
            while let Some(value) = seq.next_element_seed(Element {
                index: values.len(),
                maximum: MAX,
                _type: PhantomData,
            })? {
                values.push(value);
            }
            Ok(values)
        }
    }
    d.deserialize_seq(Values::<T, MAX>(PhantomData))
}
struct Text<'a> {
    budget: &'a mut Budget,
    cap: usize,
    original: bool,
}
impl<'de> DeserializeSeed<'de> for Text<'_> {
    type Value = String;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<String, D::Error> {
        d.deserialize_string(self)
    }
}
impl<'de> Visitor<'de> for Text<'_> {
    type Value = String;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("bounded native text")
    }
    fn visit_str<E: de::Error>(self, value: &str) -> Result<String, E> {
        self.budget
            .text(value, self.cap, self.original)
            .map_err(E::custom)?;
        Ok(value.to_owned())
    }
    fn visit_borrowed_str<E: de::Error>(self, value: &'de str) -> Result<String, E> {
        self.visit_str(value)
    }
    fn visit_string<E: de::Error>(self, value: String) -> Result<String, E> {
        self.budget
            .text(&value, self.cap, self.original)
            .map_err(E::custom)?;
        Ok(value)
    }
}
struct OptionalText<'a> {
    budget: &'a mut Budget,
    cap: usize,
    original: bool,
}
impl<'de> DeserializeSeed<'de> for OptionalText<'_> {
    type Value = Option<String>;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        d.deserialize_option(self)
    }
}
impl<'de> Visitor<'de> for OptionalText<'_> {
    type Value = Option<String>;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("optional bounded native text")
    }
    fn visit_none<E: de::Error>(self) -> Result<Self::Value, E> {
        Ok(None)
    }
    fn visit_unit<E: de::Error>(self) -> Result<Self::Value, E> {
        Ok(None)
    }
    fn visit_some<D: Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        Text {
            budget: &mut *self.budget,
            cap: self.cap,
            original: self.original,
        }
        .deserialize(d)
        .map(Some)
    }
}
fn missing<E: de::Error>(name: &'static str) -> E {
    E::missing_field(name)
}
fn required<T, E: de::Error>(value: Option<T>, name: &'static str) -> Result<T, E> {
    value.ok_or_else(|| missing(name))
}

#[derive(Deserialize)]
#[serde(field_identifier, rename_all = "snake_case")]
enum MoleculeField {
    Name,
    Sequence,
    Topology,
    Features,
    Metadata,
    OriginalGenbank,
    #[serde(other)]
    Ignore,
}
struct MoleculeSeed<'a>(&'a mut Budget);
impl<'de> DeserializeSeed<'de> for MoleculeSeed<'_> {
    type Value = Molecule;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Molecule, D::Error> {
        d.deserialize_struct(
            "Molecule",
            &[
                "name",
                "sequence",
                "topology",
                "features",
                "metadata",
                "original_genbank",
            ],
            self,
        )
    }
}
impl<'de> Visitor<'de> for MoleculeSeed<'_> {
    type Value = Molecule;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("a bounded native molecule")
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Molecule, A::Error> {
        let name = required(
            seq.next_element_seed(Text {
                budget: &mut *self.0,
                cap: 1024,
                original: false,
            })?,
            "name",
        )?;
        let sequence = required(seq.next_element::<Sequence>()?, "sequence")?;
        let topology = required(seq.next_element::<Topology>()?, "topology")?;
        let features = required(seq.next_element_seed(Features(&mut *self.0))?, "features")?;
        let metadata = required(seq.next_element_seed(Metadata(&mut *self.0))?, "metadata")?;
        let original_genbank = required(
            seq.next_element_seed(OptionalText {
                budget: &mut *self.0,
                cap: 0,
                original: true,
            })?,
            "original_genbank",
        )?;
        Ok(Molecule {
            name,
            sequence,
            topology,
            features,
            metadata,
            original_genbank,
        })
    }
    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Molecule, A::Error> {
        let (mut name, mut sequence, mut topology, mut features, mut metadata, mut original) =
            (None, None, None, None, None, None);
        while let Some(field) = map.next_key::<MoleculeField>()? {
            match field {
                MoleculeField::Name => {
                    if name.is_some() {
                        return Err(de::Error::duplicate_field("name"));
                    }
                    name = Some(map.next_value_seed(Text {
                        budget: &mut *self.0,
                        cap: 1024,
                        original: false,
                    })?);
                }
                MoleculeField::Sequence => {
                    if sequence.is_some() {
                        return Err(de::Error::duplicate_field("sequence"));
                    }
                    sequence = Some(map.next_value::<Sequence>()?);
                }
                MoleculeField::Topology => {
                    if topology.is_some() {
                        return Err(de::Error::duplicate_field("topology"));
                    }
                    topology = Some(map.next_value::<Topology>()?);
                }
                MoleculeField::Features => {
                    if features.is_some() {
                        return Err(de::Error::duplicate_field("features"));
                    }
                    features = Some(map.next_value_seed(Features(&mut *self.0))?);
                }
                MoleculeField::Metadata => {
                    if metadata.is_some() {
                        return Err(de::Error::duplicate_field("metadata"));
                    }
                    metadata = Some(map.next_value_seed(Metadata(&mut *self.0))?);
                }
                MoleculeField::OriginalGenbank => {
                    if original.is_some() {
                        return Err(de::Error::duplicate_field("original_genbank"));
                    }
                    original = Some(map.next_value_seed(OptionalText {
                        budget: &mut *self.0,
                        cap: 0,
                        original: true,
                    })?);
                }
                MoleculeField::Ignore => {
                    map.next_value::<de::IgnoredAny>()?;
                }
            }
        }
        Ok(Molecule {
            name: required(name, "name")?,
            sequence: required(sequence, "sequence")?,
            topology: required(topology, "topology")?,
            features: required(features, "features")?,
            metadata: required(metadata, "metadata")?,
            original_genbank: original.unwrap_or(None),
        })
    }
}
pub(crate) fn molecule<'de, D: Deserializer<'de>>(d: D) -> Result<Molecule, D::Error> {
    MoleculeSeed(&mut Budget::new()).deserialize(d)
}
struct Features<'a>(&'a mut Budget);
impl<'de> DeserializeSeed<'de> for Features<'_> {
    type Value = Vec<Feature>;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        d.deserialize_seq(self)
    }
}
impl<'de> Visitor<'de> for Features<'_> {
    type Value = Vec<Feature>;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("bounded native features")
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
        hint(seq.size_hint(), self.0.0.remaining_features())?;
        let mut values = Vec::new();
        while let Some(value) = seq.next_element_seed(FeatureSeed(&mut *self.0))? {
            values.push(value);
        }
        Ok(values)
    }
}
#[derive(Deserialize)]
#[serde(field_identifier, rename_all = "snake_case")]
enum FeatureField {
    Kind,
    Location,
    Qualifiers,
    NeedsReview,
    #[serde(other)]
    Ignore,
}
struct FeatureSeed<'a>(&'a mut Budget);
impl<'de> DeserializeSeed<'de> for FeatureSeed<'_> {
    type Value = Feature;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Feature, D::Error> {
        self.0.feature().map_err(de::Error::custom)?;
        d.deserialize_struct(
            "Feature",
            &["kind", "location", "qualifiers", "needs_review"],
            self,
        )
    }
}
impl<'de> Visitor<'de> for FeatureSeed<'_> {
    type Value = Feature;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("a bounded native feature")
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Feature, A::Error> {
        let kind = required(
            seq.next_element_seed(Text {
                budget: &mut *self.0,
                cap: 256,
                original: false,
            })?,
            "kind",
        )?;
        let location = required(
            seq.next_element_seed(Node {
                budget: &mut *self.0,
                depth: 0,
            })?,
            "location",
        )?;
        let qualifiers = required(
            seq.next_element_seed(Qualifiers(&mut *self.0))?,
            "qualifiers",
        )?;
        let needs_review = required(seq.next_element::<bool>()?, "needs_review")?;
        Ok(Feature {
            kind,
            location,
            qualifiers,
            needs_review,
        })
    }
    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Feature, A::Error> {
        let (mut kind, mut location, mut qualifiers, mut needs_review) = (None, None, None, None);
        while let Some(field) = map.next_key::<FeatureField>()? {
            match field {
                FeatureField::Kind => {
                    if kind.is_some() {
                        return Err(de::Error::duplicate_field("kind"));
                    }
                    kind = Some(map.next_value_seed(Text {
                        budget: &mut *self.0,
                        cap: 256,
                        original: false,
                    })?);
                }
                FeatureField::Location => {
                    if location.is_some() {
                        return Err(de::Error::duplicate_field("location"));
                    }
                    location = Some(map.next_value_seed(Node {
                        budget: &mut *self.0,
                        depth: 0,
                    })?);
                }
                FeatureField::Qualifiers => {
                    if qualifiers.is_some() {
                        return Err(de::Error::duplicate_field("qualifiers"));
                    }
                    qualifiers = Some(map.next_value_seed(Qualifiers(&mut *self.0))?);
                }
                FeatureField::NeedsReview => {
                    if needs_review.is_some() {
                        return Err(de::Error::duplicate_field("needs_review"));
                    }
                    needs_review = Some(map.next_value::<bool>()?);
                }
                FeatureField::Ignore => {
                    map.next_value::<de::IgnoredAny>()?;
                }
            }
        }
        Ok(Feature {
            kind: required(kind, "kind")?,
            location: required(location, "location")?,
            qualifiers: required(qualifiers, "qualifiers")?,
            needs_review: required(needs_review, "needs_review")?,
        })
    }
}
pub(crate) fn feature<'de, D: Deserializer<'de>>(d: D) -> Result<Feature, D::Error> {
    FeatureSeed(&mut Budget::new()).deserialize(d)
}
struct Qualifiers<'a>(&'a mut Budget);
impl<'de> DeserializeSeed<'de> for Qualifiers<'_> {
    type Value = Vec<(String, Option<String>)>;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        d.deserialize_seq(self)
    }
}
impl<'de> Visitor<'de> for Qualifiers<'_> {
    type Value = Vec<(String, Option<String>)>;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("bounded native qualifiers")
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
        hint(
            seq.size_hint(),
            QUALIFIERS.min(self.0.0.remaining_qualifiers()),
        )?;
        let mut values = Vec::new();
        while let Some(value) = seq.next_element_seed(Qualifier {
            budget: &mut *self.0,
            index: values.len(),
        })? {
            values.push(value);
        }
        Ok(values)
    }
}
struct Qualifier<'a> {
    budget: &'a mut Budget,
    index: usize,
}
impl<'de> DeserializeSeed<'de> for Qualifier<'_> {
    type Value = (String, Option<String>);
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        if self.index >= QUALIFIERS {
            return Err(de::Error::custom(
                "native feature qualifier count exceeds16384",
            ));
        }
        self.budget.qualifier().map_err(de::Error::custom)?;
        d.deserialize_tuple(2, self)
    }
}
impl<'de> Visitor<'de> for Qualifier<'_> {
    type Value = (String, Option<String>);
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("native qualifier key/value pair")
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
        let key = required(
            seq.next_element_seed(Text {
                budget: &mut *self.budget,
                cap: 4096,
                original: false,
            })?,
            "qualifier key",
        )?;
        let value = required(
            seq.next_element_seed(OptionalText {
                budget: &mut *self.budget,
                cap: 1024 * 1024,
                original: false,
            })?,
            "qualifier value",
        )?;
        if seq.next_element::<de::IgnoredAny>()?.is_some() {
            return Err(de::Error::custom("extra qualifier fields"));
        }
        Ok((key, value))
    }
}
struct Metadata<'a>(&'a mut Budget);
impl<'de> DeserializeSeed<'de> for Metadata<'_> {
    type Value = BTreeMap<String, String>;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        d.deserialize_map(self)
    }
}
impl<'de> Visitor<'de> for Metadata<'_> {
    type Value = BTreeMap<String, String>;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("bounded native metadata")
    }
    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        hint(map.size_hint(), self.0.0.remaining_metadata())?;
        let mut values = BTreeMap::new();
        while let Some(key) = map.next_key_seed(MetadataKey(&mut *self.0))? {
            let value = map.next_value_seed(Text {
                budget: &mut *self.0,
                cap: 1024 * 1024,
                original: ModelBudget::retained_metadata(&key),
            })?;
            if values.insert(key, value).is_some() {
                return Err(de::Error::custom("duplicate native metadata key"));
            }
        }
        Ok(values)
    }
}
struct MetadataKey<'a>(&'a mut Budget);
impl<'de> DeserializeSeed<'de> for MetadataKey<'_> {
    type Value = String;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<String, D::Error> {
        self.0.metadata().map_err(de::Error::custom)?;
        Text {
            budget: &mut *self.0,
            cap: 4096,
            original: false,
        }
        .deserialize(d)
    }
}
#[derive(Deserialize)]
#[serde(field_identifier)]
enum NodeVariant {
    Range,
    Complement,
    Join,
    Order,
    Unresolved,
}
struct Node<'a> {
    budget: &'a mut Budget,
    depth: usize,
}
impl<'de> DeserializeSeed<'de> for Node<'_> {
    type Value = Location;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Location, D::Error> {
        if self.depth >= crate::bounds::MAX_DEPTH {
            return Err(de::Error::custom("location nesting exceeds32levels"));
        }
        self.budget.node().map_err(de::Error::custom)?;
        d.deserialize_enum(
            "Location",
            &["Range", "Complement", "Join", "Order", "Unresolved"],
            self,
        )
    }
}
impl<'de> Visitor<'de> for Node<'_> {
    type Value = Location;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("a bounded native location")
    }
    fn visit_enum<A: EnumAccess<'de>>(self, data: A) -> Result<Location, A::Error> {
        let (variant, value) = data.variant::<NodeVariant>()?;
        Ok(match variant {
            NodeVariant::Range => {
                value.struct_variant(&["start", "end", "before", "after"], Range)?
            }
            NodeVariant::Complement => {
                Location::Complement(Box::new(value.newtype_variant_seed(Node {
                    budget: &mut *self.budget,
                    depth: self.depth + 1,
                })?))
            }
            NodeVariant::Join => Location::Join(value.newtype_variant_seed(Parts {
                budget: &mut *self.budget,
                depth: self.depth + 1,
            })?),
            NodeVariant::Order => Location::Order(value.newtype_variant_seed(Parts {
                budget: &mut *self.budget,
                depth: self.depth + 1,
            })?),
            NodeVariant::Unresolved => Location::Unresolved(value.newtype_variant_seed(Text {
                budget: &mut *self.budget,
                cap: 65536,
                original: false,
            })?),
        })
    }
}
pub(crate) fn location<'de, D: Deserializer<'de>>(d: D) -> Result<Location, D::Error> {
    Node {
        budget: &mut Budget::new(),
        depth: 0,
    }
    .deserialize(d)
}
struct Parts<'a> {
    budget: &'a mut Budget,
    depth: usize,
}
impl<'de> DeserializeSeed<'de> for Parts<'_> {
    type Value = Vec<Location>;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        d.deserialize_seq(self)
    }
}
impl<'de> Visitor<'de> for Parts<'_> {
    type Value = Vec<Location>;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("bounded location parts")
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
        hint(seq.size_hint(), PARTS.min(self.budget.0.remaining_nodes()))?;
        let mut values = Vec::new();
        while let Some(value) = seq.next_element_seed(Part {
            budget: &mut *self.budget,
            depth: self.depth,
            index: values.len(),
        })? {
            values.push(value);
        }
        Ok(values)
    }
}
struct Part<'a> {
    budget: &'a mut Budget,
    depth: usize,
    index: usize,
}
impl<'de> DeserializeSeed<'de> for Part<'_> {
    type Value = Location;
    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Location, D::Error> {
        if self.index >= PARTS {
            return Err(de::Error::custom("native location parts exceed16384"));
        }
        Node {
            budget: &mut *self.budget,
            depth: self.depth,
        }
        .deserialize(d)
    }
}
#[derive(Deserialize)]
#[serde(field_identifier, rename_all = "snake_case")]
enum RangeField {
    Start,
    End,
    Before,
    After,
    #[serde(other)]
    Ignore,
}
struct Range;
impl<'de> Visitor<'de> for Range {
    type Value = Location;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("native range fields")
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Location, A::Error> {
        let start = required(seq.next_element::<usize>()?, "start")?;
        let end = required(seq.next_element::<usize>()?, "end")?;
        let before = required(seq.next_element::<bool>()?, "before")?;
        let after = required(seq.next_element::<bool>()?, "after")?;
        Ok(Location::Range {
            start,
            end,
            before,
            after,
        })
    }
    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Location, A::Error> {
        let (mut start, mut end, mut before, mut after) = (None, None, None, None);
        while let Some(field) = map.next_key::<RangeField>()? {
            match field {
                RangeField::Start => {
                    if start.is_some() {
                        return Err(de::Error::duplicate_field("start"));
                    }
                    start = Some(map.next_value::<usize>()?);
                }
                RangeField::End => {
                    if end.is_some() {
                        return Err(de::Error::duplicate_field("end"));
                    }
                    end = Some(map.next_value::<usize>()?);
                }
                RangeField::Before => {
                    if before.is_some() {
                        return Err(de::Error::duplicate_field("before"));
                    }
                    before = Some(map.next_value::<bool>()?);
                }
                RangeField::After => {
                    if after.is_some() {
                        return Err(de::Error::duplicate_field("after"));
                    }
                    after = Some(map.next_value::<bool>()?);
                }
                RangeField::Ignore => {
                    map.next_value::<de::IgnoredAny>()?;
                }
            }
        }
        Ok(Location::Range {
            start: required(start, "start")?,
            end: required(end, "end")?,
            before: required(before, "before")?,
            after: required(after, "after")?,
        })
    }
}

#[cfg(test)]
mod tests;
