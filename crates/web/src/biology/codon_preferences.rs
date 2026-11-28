//! One active codon-preference source; custom drafts never supplement host profiles.
use collagene_core::sequence::codon_design::{CodonPreferences, CodonWeight, hosts};

#[derive(Clone, Debug, PartialEq, Eq, Default)]
enum Choice {
    #[default]
    Uniform,
    Custom,
    Host(hosts::Reference),
}

#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Preferences {
    choice: Choice,
    pub custom: String,
}
impl Preferences {
    pub fn key(&self) -> String {
        match &self.choice {
            Choice::Uniform => "uniform".into(),
            Choice::Custom => "custom".into(),
            Choice::Host(reference) => format!("host:{}", reference.id),
        }
    }
    pub fn is_host(&self) -> bool {
        matches!(self.choice, Choice::Host(_))
    }
    pub fn host<'a>(&self, catalogue: &'a [hosts::Summary]) -> Option<&'a hosts::Summary> {
        let Choice::Host(reference) = &self.choice else {
            return None;
        };
        catalogue.iter().find(|host| &host.reference == reference)
    }
    /// Returns the host's required genetic code. Unknown selections leave state unchanged.
    pub fn select(
        &mut self,
        key: &str,
        catalogue: &[hosts::Summary],
    ) -> Result<Option<u8>, String> {
        let (choice, code) = match key {
            "uniform" => (Choice::Uniform, None),
            "custom" => (Choice::Custom, None),
            _ => {
                let host = catalogue
                    .iter()
                    .find(|host| key == format!("host:{}", host.reference.id))
                    .ok_or("The selected host is unavailable. Choose a current catalogue entry.")?;
                (
                    Choice::Host(host.reference.clone()),
                    Some(host.genetic_code),
                )
            }
        };
        self.choice = choice;
        Ok(code)
    }
    pub fn edit(&mut self, text: String) {
        self.custom = text;
        self.choice = Choice::Custom;
    }
    pub fn build(
        &self,
        genetic_code: u8,
        catalogue: &[hosts::Summary],
    ) -> Result<CodonPreferences, String> {
        match &self.choice {
            Choice::Uniform => Ok(CodonPreferences::Uniform),
            Choice::Host(reference) => {
                let host = self.host(catalogue).ok_or("The pinned host identity is unavailable. Choose its current catalogue entry again.")?;
                if genetic_code != host.genetic_code {
                    return Err("The selected host requires its declared genetic code. Choose the host again or switch to uniform/custom preferences.".into());
                }
                Ok(CodonPreferences::ReferenceHost {
                    reference: reference.clone(),
                })
            }
            Choice::Custom => {
                if self.custom.len() > 4096 {
                    return Err("Custom codon preferences exceed the 4096-byte input limit.".into());
                }
                let mut weights = Vec::new();
                for line in self.custom.lines().filter(|line| !line.trim().is_empty()) {
                    if weights.len() == 64 {
                        return Err("Custom preferences accept at most 64 codons.".into());
                    }
                    let fields = line
                        .split(|c: char| c.is_ascii_whitespace() || c == ',' || c == '=')
                        .filter(|field| !field.is_empty())
                        .collect::<Vec<_>>();
                    if fields.len() != 2 {
                        return Err("Use one DNA codon and nonnegative integer weight per line, for example AAG 100.".into());
                    }
                    let codon = fields[0].to_ascii_uppercase();
                    if codon.len() != 3 || !codon.bytes().all(|base| b"ACGT".contains(&base)) {
                        return Err(
                            "Each custom preference must name one concrete DNA codon.".into()
                        );
                    }
                    if weights
                        .iter()
                        .any(|weight: &CodonWeight| weight.codon == codon)
                    {
                        return Err("Each custom codon may appear only once.".into());
                    }
                    let weight = fields[1].parse::<u32>().map_err(|_| "Codon weights must be nonnegative whole numbers within the 32-bit range.")?;
                    weights.push(CodonWeight { codon, weight });
                }
                if weights.is_empty() {
                    return Err("Enter at least one custom codon weight, or choose Uniform.".into());
                }
                Ok(CodonPreferences::Explicit { weights })
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn preference_sources_are_exclusive_and_host_identity_and_code_are_pinned() {
        let catalogue = hosts::catalogue().unwrap();
        let host = catalogue
            .first()
            .expect("at least one curated reference host");
        let mut editor = Preferences::default();
        editor.edit("AAG 0".into());
        let key = format!("host:{}", host.reference.id);
        assert_eq!(
            editor.select(&key, &catalogue).unwrap(),
            Some(host.genetic_code)
        );
        assert_eq!(
            editor.build(host.genetic_code, &catalogue).unwrap(),
            CodonPreferences::ReferenceHost {
                reference: host.reference.clone()
            }
        );
        assert!(editor.build(255, &catalogue).is_err());
        assert!(editor.build(host.genetic_code, &[]).is_err());
        let mut replaced_catalogue = catalogue.clone();
        replaced_catalogue[0].reference.identity.push_str("changed");
        assert!(
            editor
                .build(host.genetic_code, &replaced_catalogue)
                .is_err()
        );
        let before = editor.clone();
        assert!(editor.select("host:missing", &catalogue).is_err());
        assert_eq!(editor, before);
        editor.select("uniform", &catalogue).unwrap();
        assert_eq!(
            editor.build(1, &catalogue).unwrap(),
            CodonPreferences::Uniform
        );
        editor.edit("aaa=0".into());
        assert!(!editor.is_host());
        assert_eq!(
            editor.build(1, &catalogue).unwrap(),
            CodonPreferences::Explicit {
                weights: vec![CodonWeight {
                    codon: "AAA".into(),
                    weight: 0
                }]
            }
        );
    }
    #[test]
    fn custom_preferences_reject_empty_ambiguous_duplicate_and_excessive_entries() {
        let mut editor = Preferences::default();
        for invalid in [
            "",
            "AAN 2",
            "AAG -1",
            "AAG 4294967296",
            "AAG 1\naag 2",
            "AAG 1 extra",
        ] {
            editor.edit(invalid.into());
            assert!(editor.build(1, &[]).is_err(), "{invalid}");
        }
        let all_codons = (0..64)
            .map(|index| {
                let alphabet = b"ACGT";
                format!(
                    "{}{}{} 1\n",
                    alphabet[index / 16] as char,
                    alphabet[index / 4 % 4] as char,
                    alphabet[index % 4] as char
                )
            })
            .collect::<String>();
        editor.edit(all_codons.clone());
        assert!(editor.build(1, &[]).is_ok());
        editor.edit(format!("{all_codons}AAA 1\n"));
        assert!(editor.build(1, &[]).is_err());
        editor.edit(" ".repeat(4097));
        assert!(editor.build(1, &[]).is_err());
    }
}
