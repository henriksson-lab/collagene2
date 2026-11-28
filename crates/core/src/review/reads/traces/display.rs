//! Source-bound chromatogram display preparation. Original attachments stay in
//! the document and never travel back as display or export authority.
use super::{
    MAX_TRACE_BASES, MAX_TRACE_SAMPLES, QualityEvidence, Trace, TraceChannels, TraceFormat,
    TraceOrientation,
};
use collagene_model::{Alphabet, Error, Molecule, Sequence};
use serde::{Deserialize, Serialize};
pub const ENGINE: &str = "trace-display-1";
pub const BASE_BIN: usize = 64;
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SignalPyramid {
    levels: Vec<Vec<(i32, i32)>>,
}
impl SignalPyramid {
    fn prepare(values: &[i32], cancel: &mut impl FnMut() -> bool) -> Result<Self, Error> {
        let mut base = Vec::with_capacity(values.len().div_ceil(BASE_BIN));
        for chunk in values.chunks(BASE_BIN) {
            if cancel() {
                return Err(Error::Unsupported("trace display cancelled".into()));
            }
            let mut low = i32::MAX;
            let mut high = i32::MIN;
            for &v in chunk {
                low = low.min(v);
                high = high.max(v);
            }
            base.push((low, high));
        }
        let mut levels = vec![base];
        while levels.last().is_some_and(|v| v.len() > 1) {
            if cancel() {
                return Err(Error::Unsupported("trace display cancelled".into()));
            }
            let next = levels
                .last()
                .unwrap()
                .chunks(2)
                .map(|pair| {
                    pair.iter()
                        .copied()
                        .fold((i32::MAX, i32::MIN), |(lo, hi), (a, b)| {
                            (lo.min(a), hi.max(b))
                        })
                })
                .collect();
            levels.push(next);
        }
        Ok(Self { levels })
    }
    pub fn levels(&self) -> &[Vec<(i32, i32)>] {
        &self.levels
    }
    fn validate_shape(&self, samples: usize) -> Result<(), Error> {
        let mut count = samples.div_ceil(BASE_BIN);
        let mut level = 0;
        loop {
            let values = self
                .levels
                .get(level)
                .ok_or_else(|| bad("missing pyramid level"))?;
            if values.len() != count || values.iter().any(|&(low, high)| low > high) {
                return Err(bad("pyramid dimensions or extrema changed"));
            }
            level += 1;
            if count == 1 {
                break;
            }
            count = count.div_ceil(2);
        }
        if self.levels.len() != level {
            return Err(bad("extra pyramid levels"));
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DisplayTrace {
    pub name: String,
    pub format: TraceFormat,
    pub orientation: TraceOrientation,
    pub sequence: Sequence,
    pub called_symbols: Vec<u8>,
    pub peak_positions: Vec<u32>,
    pub channels: TraceChannels,
    pub qualities: Option<QualityEvidence>,
    pub warnings: Vec<String>,
    pub source_hash: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlotData {
    pub engine: String,
    pub source_revision: String,
    pub trace: DisplayTrace,
    pub pyramids: [SignalPyramid; 4],
    pub minimum: i32,
    pub maximum: i32,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("trace display: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
impl PlotData {
    /// All bulk source parsing, orientation and extrema work occurs here, on
    /// the worker. No caller-supplied display summaries are accepted as inputs.
    pub fn prepare(molecule: &Molecule, cancel: &mut impl FnMut() -> bool) -> Result<Self, Error> {
        if cancel() {
            return Err(Error::Unsupported("trace display cancelled".into()));
        }
        molecule.validate()?;
        let source_revision = molecule.revision();
        let trace = super::trace_from_molecule(molecule)?;
        Self::from_trace(source_revision, trace, cancel)
    }
    fn from_trace(
        source_revision: String,
        trace: Trace,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let values = trace.channels.all();
        let mut pyramids = Vec::with_capacity(4);
        let mut minimum = 0i32;
        let mut maximum = 0i32;
        for channel in values {
            let pyramid = SignalPyramid::prepare(channel, cancel)?;
            let &(low, high) = pyramid
                .levels
                .last()
                .and_then(|v| v.first())
                .ok_or_else(|| bad("empty signal"))?;
            minimum = minimum.min(low);
            maximum = maximum.max(high);
            pyramids.push(pyramid);
        }
        let Trace {
            name,
            format,
            orientation,
            sequence,
            called_symbols,
            peak_positions,
            channels,
            qualities,
            warnings,
            source_hash,
            ..
        } = trace;
        let result = Self {
            engine: ENGINE.into(),
            source_revision,
            trace: DisplayTrace {
                name,
                format,
                orientation,
                sequence,
                called_symbols,
                peak_positions,
                channels,
                qualities,
                warnings,
                source_hash,
            },
            pyramids: pyramids
                .try_into()
                .map_err(|_| bad("four pyramids required"))?,
            minimum,
            maximum,
        };
        result.validate_binding(&result.source_revision)?;
        Ok(result)
    }
    /// Structural validation of a trusted local worker reply. It does not
    /// reparse original bytes or rescan raw samples/min-max bins on the UI
    /// thread, and grants no authority for trace export or scientific edits.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let t = &self.trace;
        let samples = t.channels.len();
        let n = t.sequence.len();
        if self.engine != ENGINE
            || self.source_revision != expected
            || !hex(expected)
            || !hex(&t.source_hash)
            || samples == 0
            || samples > MAX_TRACE_SAMPLES
            || n > MAX_TRACE_BASES
            || t.name.trim().is_empty()
            || t.name.len() > 1024
            || t.sequence.alphabet() != Alphabet::Dna
            || t.channels.all().iter().any(|v| v.len() != samples)
            || t.peak_positions.len() != n
            || t.called_symbols.len() != n
            || t.peak_positions.iter().any(|&p| p as usize >= samples)
            || t.peak_positions.windows(2).any(|v| v[0] > v[1])
            || self.minimum > 0
            || self.maximum < 0
            || self.minimum > self.maximum
            || t.warnings.len() > 64
            || t.warnings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad("source or display dimensions changed"));
        }
        let dash = matches!(t.format, TraceFormat::Scf { code_set: 0, .. });
        if super::calls(&t.called_symbols, dash)? != t.sequence {
            return Err(bad("call symbols changed"));
        }
        match (&t.format, &t.qualities) {
            (TraceFormat::Abif { version }, _) if version / 100 != 1 => {
                return Err(bad("unsupported trace version"));
            }
            (TraceFormat::Abif { .. }, Some(QualityEvidence::ScfAccuracy(_)))
            | (TraceFormat::Scf { .. }, Some(QualityEvidence::AbifBasecallerScores(_))) => {
                return Err(bad("quality source semantics changed"));
            }
            _ => {}
        }
        match &t.qualities {
            Some(QualityEvidence::AbifBasecallerScores(v)) if v.len() != n => {
                return Err(bad("quality count changed"));
            }
            Some(QualityEvidence::ScfAccuracy(q))
                if [&q.a, &q.c, &q.g, &q.t].iter().any(|v| v.len() != n) =>
            {
                return Err(bad("quality count changed"));
            }
            _ => {}
        }
        let mut minimum = 0i32;
        let mut maximum = 0i32;
        for pyramid in &self.pyramids {
            pyramid.validate_shape(samples)?;
            let (low, high) = pyramid.levels.last().unwrap()[0];
            minimum = minimum.min(low);
            maximum = maximum.max(high);
        }
        if self.minimum != minimum || self.maximum != maximum {
            return Err(bad("amplitude range changed"));
        }
        Ok(())
    }
    pub fn max_zoom(&self) -> f64 {
        (16_000_000. / self.trace.channels.len() as f64)
            .min(128.)
            .max(0.125)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn literal_extrema_keep_signed_values_partial_bins_and_hierarchy() {
        let mut values = vec![7; 129];
        values[0] = -11;
        values[63] = 23;
        values[64] = -31;
        values[127] = 47;
        values[128] = -53;
        let p = SignalPyramid::prepare(&values, &mut || false).unwrap();
        assert_eq!(
            p.levels,
            vec![
                vec![(-11, 23), (-31, 47), (-53, -53)],
                vec![(-31, 47), (-53, -53)],
                vec![(-53, 47)]
            ]
        );
        p.validate_shape(129).unwrap();
        assert!(p.validate_shape(128).is_err());
    }
    #[test]
    fn every_overview_bin_matches_independent_literal_source_slices() {
        let values = (0..10003)
            .map(|i| ((i * 97) % 503) as i32 - 251)
            .collect::<Vec<_>>();
        let p = SignalPyramid::prepare(&values, &mut || false).unwrap();
        for (level, bins) in p.levels.iter().enumerate() {
            let width = 64usize << level;
            for (i, &(low, high)) in bins.iter().enumerate() {
                let slice = &values[i * width..((i + 1) * width).min(values.len())];
                assert!(slice.iter().all(|&v| v >= low && v <= high));
                assert!(slice.contains(&low) && slice.contains(&high));
            }
        }
        assert_eq!(p.levels[0].len(), 157);
        assert!(SignalPyramid::prepare(&values, &mut || true).is_err());
    }
    #[test]
    fn maximum_display_shape_uses_only_its_explicit_large_reply_selector() {
        use crate::protocol::{
            response::Response,
            schema::{AnalysisResult, WORKER_SCHEMA_VERSION},
        };
        // Codec/resource fixture, not evidence of an original chromatogram.
        let values = vec![i32::MAX; MAX_TRACE_SAMPLES];
        let channels = TraceChannels {
            a: values.clone(),
            c: values.clone(),
            g: values.clone(),
            t: values,
        };
        let pyramids =
            std::array::from_fn(|_| SignalPyramid::prepare(&channels.a, &mut || false).unwrap());
        let revision = "a".repeat(64);
        let data = PlotData {
            engine: ENGINE.into(),
            source_revision: revision.clone(),
            trace: DisplayTrace {
                name: "bounded display".into(),
                format: TraceFormat::Abif { version: 100 },
                orientation: TraceOrientation::AsStored,
                sequence: Sequence::dna("A").unwrap(),
                called_symbols: vec![b'A'],
                peak_positions: vec![0],
                channels,
                qualities: None,
                warnings: vec![],
                source_hash: "b".repeat(64),
            },
            pyramids,
            minimum: 0,
            maximum: i32::MAX,
        };
        data.validate_binding(&revision).unwrap();
        let response = Response {
            id: 1,
            revision,
            result: AnalysisResult::TraceDisplayPrepared(data),
        };
        let bytes = response.encode().unwrap();
        assert!(bytes.len() > crate::formats::sequence::MAX_IMPORT_BYTES);
        assert!(matches!(
            Response::decode(&bytes).unwrap().result,
            AnalysisResult::TraceDisplayPrepared(_)
        ));
        let ((id, revision, _), tail) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        let mut payload = postcard::to_allocvec(&(id, revision, 0u32)).unwrap();
        payload.extend_from_slice(tail);
        let mut ordinary = bytes[..42].to_vec();
        ordinary[6..10].copy_from_slice(&(payload.len() as u32).to_le_bytes());
        ordinary[10..42].copy_from_slice(blake3::hash(&payload).as_bytes());
        ordinary.extend(payload);
        assert!(matches!(Response::decode(&ordinary), Err(Error::TooLarge)));
        let mut previous = bytes;
        previous[4..6].copy_from_slice(&(WORKER_SCHEMA_VERSION - 1).to_le_bytes());
        assert!(matches!(
            Response::decode(&previous),
            Err(Error::Unsupported(_))
        ));
    }
}
