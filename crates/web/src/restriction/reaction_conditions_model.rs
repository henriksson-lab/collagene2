//! Exact formulation choices and text parsing for documentary reaction review.
use collagene_core::{
    cloning::restriction::{
        catalog::{CatalogSelection, catalog},
        conditions::{ChosenFormulation, ReactionConditions},
    },
    material::lineage::EventBinding,
};
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Draft {
    pub enabled: bool,
    pub formulations: BTreeMap<String, Option<String>>,
    pub buffer: String,
    pub temperature: String,
    pub minutes: String,
    pub minimum: String,
    pub activity_context: String,
    pub temperature_context: String,
}
impl Default for Draft {
    fn default() -> Self {
        Self {
            enabled: false,
            formulations: BTreeMap::new(),
            buffer: String::new(),
            temperature: "37".into(),
            minutes: "60".into(),
            minimum: "50".into(),
            activity_context: String::new(),
            temperature_context: String::new(),
        }
    }
}
impl Draft {
    pub fn from_saved(value: Option<ReactionConditions>) -> Self {
        let Some(value) = value else {
            return Self::default();
        };
        Self {
            enabled: true,
            formulations: value
                .formulations
                .into_iter()
                .map(|f| (f.enzyme_id, f.formulation_id))
                .collect(),
            buffer: value.buffer,
            temperature: celsius(value.temperature_milli_c),
            minutes: value.incubation_minutes.to_string(),
            minimum: value.minimum_activity_percent.to_string(),
            activity_context: value.activity_context.unwrap_or_default(),
            temperature_context: value.temperature_context.unwrap_or_default(),
        }
    }
    pub fn options(
        &self,
        c: &catalog::Catalog,
        s: &CatalogSelection,
    ) -> Result<Option<ReactionConditions>, String> {
        if !self.enabled {
            return Ok(None);
        }
        let context = |value: &str| {
            if value.is_empty() {
                None
            } else {
                Some(value.to_string())
            }
        };
        let options = ReactionConditions {
            schema_version: 1,
            formulations: s
                .enzymes
                .iter()
                .map(|e| ChosenFormulation {
                    enzyme_id: e.id.clone(),
                    formulation_id: self.formulations.get(&e.id).cloned().flatten(),
                })
                .collect(),
            buffer: self.buffer.trim().into(),
            temperature_milli_c: parse_celsius(&self.temperature)?,
            incubation_minutes: self
                .minutes
                .trim()
                .parse()
                .map_err(|_| "Enter whole incubation minutes.".to_string())?,
            minimum_activity_percent: self
                .minimum
                .trim()
                .parse()
                .map_err(|_| "Enter a whole activity threshold percent.".to_string())?,
            activity_context: context(&self.activity_context),
            temperature_context: context(&self.temperature_context),
        };
        options.validate_for(c, s).map_err(|e| e.to_string())?;
        Ok(Some(options))
    }
    pub fn keep_selected(&mut self, selection: &CatalogSelection) {
        self.formulations
            .retain(|enzyme, _| selection.enzymes.iter().any(|e| &e.id == enzyme));
    }
    pub fn selected_forms<'a>(
        &self,
        c: &'a catalog::Catalog,
        s: &CatalogSelection,
    ) -> Vec<&'a catalog::Formulation> {
        s.enzymes
            .iter()
            .filter_map(|e| {
                let id = self.formulations.get(&e.id)?.as_ref()?;
                c.edition()
                    .formulations
                    .iter()
                    .find(|f| &f.id == id && f.enzyme == e.id)
            })
            .collect()
    }
    pub fn buffers(&self, c: &catalog::Catalog, s: &CatalogSelection) -> Vec<String> {
        let mut values = self
            .selected_forms(c, s)
            .into_iter()
            .flat_map(|f| f.evidence.iter())
            .filter_map(|fact| match &fact.condition {
                catalog::Condition::BufferActivity { buffer, .. } => Some(buffer.clone()),
                _ => None,
            })
            .collect::<Vec<_>>();
        values.sort();
        values.dedup();
        values
    }
    pub fn contexts(
        &self,
        c: &catalog::Catalog,
        s: &CatalogSelection,
        temperature: bool,
    ) -> Vec<String> {
        let mut values = self
            .selected_forms(c, s)
            .into_iter()
            .flat_map(|f| f.evidence.iter())
            .filter(|f| {
                if temperature {
                    matches!(&f.condition, catalog::Condition::IncubationTemperature)
                } else {
                    matches!(&f.condition, catalog::Condition::BufferActivity { .. })
                }
            })
            .map(|f| f.context.clone())
            .collect::<Vec<_>>();
        values.sort();
        values.dedup();
        values
    }
}
pub(super) fn parse_celsius(value: &str) -> Result<i32, String> {
    let invalid =
        || "Enter a temperature from 0 to 100 °C, with at most three decimal places.".to_string();
    let value = value.trim();
    let (whole, fraction) = value.split_once('.').map_or((value, ""), |(a, b)| (a, b));
    if whole.is_empty()
        || !whole.bytes().all(|b| b.is_ascii_digit())
        || fraction.len() > 3
        || !fraction.bytes().all(|b| b.is_ascii_digit())
        || (value.contains('.') && fraction.is_empty())
    {
        return Err(invalid());
    }
    let whole: i32 = whole.parse().map_err(|_| invalid())?;
    let fraction: i32 = if fraction.is_empty() {
        0
    } else {
        fraction.parse::<i32>().map_err(|_| invalid())? * 10i32.pow(3 - fraction.len() as u32)
    };
    let milli = whole
        .checked_mul(1000)
        .and_then(|v| v.checked_add(fraction))
        .ok_or_else(invalid)?;
    if !(0..=100_000).contains(&milli) {
        return Err(invalid());
    }
    Ok(milli)
}
pub(super) fn celsius(value: i32) -> String {
    let sign = if value < 0 { "-" } else { "" };
    let absolute = i64::from(value).abs();
    if absolute % 1000 == 0 {
        format!("{sign}{}", absolute / 1000)
    } else {
        format!("{sign}{}.{:03}", absolute / 1000, absolute % 1000)
            .trim_end_matches('0')
            .into()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Capture {
    pub source: EventBinding,
    pub epoch: u64,
    pub selection: CatalogSelection,
    pub options: ReactionConditions,
}
impl Capture {
    pub fn current(
        &self,
        source: &EventBinding,
        epoch: u64,
        selection: &CatalogSelection,
        options: &ReactionConditions,
    ) -> bool {
        self.source == *source
            && self.epoch == epoch
            && self.selection == *selection
            && self.options == *options
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::cloning::restriction::catalog;
    fn selected() -> (catalog::catalog::Catalog, CatalogSelection) {
        let c = catalog::builtin_catalog().unwrap();
        let ids = c
            .edition()
            .enzymes
            .iter()
            .filter(|e| e.name == "EcoRI")
            .map(|e| e.id.clone())
            .collect::<Vec<_>>();
        let s = catalog::select(&c, &ids).unwrap();
        (c, s)
    }
    #[test]
    fn celsius_decimal_is_exact_and_refuses_nonfinite_fractional_noise() {
        for (input, expected) in [
            ("37", 37000),
            ("37.5", 37500),
            ("0.001", 1),
            ("100.000", 100000),
        ] {
            assert_eq!(parse_celsius(input).unwrap(), expected);
            assert_eq!(parse_celsius(&celsius(expected)).unwrap(), expected);
        }
        for input in [
            "NaN",
            "inf",
            "1e2",
            "-0.1",
            "37.",
            "37.0001",
            "100.001",
            "2147483647",
        ] {
            assert!(parse_celsius(input).is_err());
        }
    }
    #[test]
    fn opt_in_requires_valid_values_and_never_guesses_a_formulation() {
        let (c, s) = selected();
        let mut d = Draft::default();
        assert_eq!(d.options(&c, &s).unwrap(), None);
        d.enabled = true;
        assert!(d.options(&c, &s).is_err());
        d.buffer = "unrecorded buffer".into();
        let options = d.options(&c, &s).unwrap().unwrap();
        assert!(
            options
                .formulations
                .iter()
                .all(|f| f.formulation_id.is_none())
        );
        assert!(d.buffers(&c, &s).is_empty());
        d.formulations
            .insert(s.enzymes[0].id.clone(), Some("invented-HF".into()));
        assert!(d.options(&c, &s).is_err());
    }
    #[test]
    fn saved_none_clears_draft_and_capture_rejects_event_epoch_and_option_changes() {
        let (c, s) = selected();
        let mut d = Draft::default();
        d.enabled = true;
        d.buffer = "review buffer".into();
        let o = d.options(&c, &s).unwrap().unwrap();
        assert_eq!(
            Draft::from_saved(Some(o.clone())).options(&c, &s).unwrap(),
            Some(o.clone())
        );
        assert!(!Draft::from_saved(None).enabled);
        let source = EventBinding {
            document_id: "d".into(),
            event_id: "event-1".into(),
            revision: "r".into(),
        };
        let capture = Capture {
            source: source.clone(),
            epoch: 1,
            selection: s.clone(),
            options: o.clone(),
        };
        assert!(capture.current(&source, 1, &s, &o));
        assert!(!capture.current(&source, 2, &s, &o));
        let mut changed = source.clone();
        changed.event_id = "event-2".into();
        assert!(!capture.current(&changed, 1, &s, &o));
        let mut changed = o;
        changed.incubation_minutes += 1;
        assert!(!capture.current(&source, 1, &s, &changed));
    }
}
