//! Editable graph choices and exact source/settings capture, never biological proof.
use crate::{biology::construct_qc_model::Draft as QcDraft, sources::identity::SourceInfo};
use collagene_core::construct_qc::{self as native, vector_class as vector};
use collagene_model::Molecule;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct BindingDraft {
    pub key: u64,
    pub id: String,
    pub role: String,
    pub feature: String,
    pub reference: String,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SpacingDraft {
    pub key: u64,
    pub id: String,
    pub upstream: String,
    pub downstream: String,
    pub minimum: String,
    pub maximum: String,
    pub rationale: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct Form {
    pub profile: String,
    pub bindings: Vec<BindingDraft>,
    pub first: String,
    pub last: String,
    pub reverse: bool,
    pub spacing: Vec<SpacingDraft>,
    next_key: u64,
}
fn whole(value: &str, label: &str) -> Result<usize, String> {
    let value = value.trim();
    if value.is_empty() || value.len() > 20 || !value.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(format!("Enter a whole number for {label}."));
    }
    value
        .parse()
        .map_err(|_| format!("{label} exceeds the supported range."))
}
impl Form {
    pub fn for_profile(profile: &vector::Profile) -> Result<Self, String> {
        let mut form = Self {
            profile: profile.reference.id.clone(),
            ..Default::default()
        };
        for role in &profile.definition.roles {
            form.add_binding(&role.id)?;
        }
        Ok(form)
    }
    fn key(&mut self) -> Result<u64, String> {
        let key = self.next_key;
        self.next_key = self
            .next_key
            .checked_add(1)
            .ok_or("Graph control identity limit reached.")?;
        Ok(key)
    }
    pub fn add_binding(&mut self, role: &str) -> Result<(), String> {
        if self.bindings.len() >= 64 {
            return Err("At most 64 role bindings are supported.".into());
        }
        let key = self.key()?;
        self.bindings.push(BindingDraft {
            key,
            id: format!("{role}-{}", key + 1),
            role: role.into(),
            feature: String::new(),
            reference: String::new(),
        });
        Ok(())
    }
    pub fn add_spacing(&mut self) -> Result<(), String> {
        if self.spacing.len() >= 32 {
            return Err("At most 32 spacing constraints are supported.".into());
        }
        let key = self.key()?;
        self.spacing.push(SpacingDraft {
            key,
            id: format!("spacing-{}", key + 1),
            upstream: String::new(),
            downstream: String::new(),
            minimum: "0".into(),
            maximum: String::new(),
            rationale: String::new(),
        });
        Ok(())
    }
    pub fn active_bindings(&self) -> impl Iterator<Item = &BindingDraft> {
        self.bindings
            .iter()
            .filter(|binding| !binding.feature.is_empty() && !binding.reference.is_empty())
    }
    pub fn request(&self, source: &Molecule, parent: &QcDraft) -> Result<vector::Request, String> {
        if self.bindings.len() > 64
            || self.spacing.len() > 32
            || self.profile.len() > 128
            || self.first.len() > 128
            || self.last.len() > 128
        {
            return Err("Vector graph controls exceed their resource limits.".into());
        }
        // Reuse the scientific owner's borrowed admission before options_for
        // creates current source/annotation pins or copies shared reference DNA.
        native::Options::default()
            .validate(source)
            .map_err(|error| error.to_string())?;
        let profile = vector::catalogue()
            .map_err(|error| error.to_string())?
            .into_iter()
            .find(|profile| profile.reference.id == self.profile)
            .ok_or("Choose an installed declared vector profile.")?;
        // Admit borrowed shared declarations before options_for clones them.
        parent
            .options
            .validate(source)
            .map_err(|error| error.to_string())?;
        let options = parent.options_for(source)?;
        let mut bindings = Vec::new();
        for binding in &self.bindings {
            match (binding.feature.is_empty(), binding.reference.is_empty()) {
                (true, true) => continue,
                (true, false) | (false, true) => {
                    return Err(format!(
                        "Binding '{}': choose both an annotation and a pinned reference, or clear both.",
                        binding.id
                    ));
                }
                _ => {}
            }
            if binding.id.len() > 128 || binding.reference.len() > 128 || binding.role.len() > 128 {
                return Err("Binding identifiers are limited to 128 characters.".into());
            }
            bindings.push(vector::RoleBinding {
                id: binding.id.clone(),
                role_id: binding.role.clone(),
                feature_index: whole(&binding.feature, "annotation index")?,
                reference_id: binding.reference.clone(),
            });
        }
        let boundaries = match (self.first.is_empty(), self.last.is_empty()) {
            (true, true) => None,
            (false, false) => Some(vector::BoundarySelection {
                start_binding: self.first.clone(),
                end_binding: self.last.clone(),
                reverse: self.reverse,
            }),
            _ => return Err("Choose both boundary bindings, or leave both unset.".into()),
        };
        let spacing = self
            .spacing
            .iter()
            .map(|spacing| {
                if spacing.id.len() > 128
                    || spacing.upstream.len() > 128
                    || spacing.downstream.len() > 128
                    || spacing.rationale.len() > 1024
                {
                    return Err(
                        "Spacing identifiers or rationale exceed their supported length.".into(),
                    );
                }
                Ok(vector::DeclaredSpacing {
                    id: spacing.id.clone(),
                    upstream_binding: spacing.upstream.clone(),
                    downstream_binding: spacing.downstream.clone(),
                    minimum_gap: whole(&spacing.minimum, "minimum gap")?,
                    maximum_gap: (!spacing.maximum.trim().is_empty())
                        .then(|| whole(&spacing.maximum, "maximum gap"))
                        .transpose()?,
                    rationale: spacing.rationale.clone(),
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        let request = vector::Request {
            source_revision: source.revision(),
            profile: profile.reference,
            annotations: options.annotations,
            references: options.references,
            bindings,
            boundaries,
            spacing,
            literal_parts: options.required_parts,
            literal_part_limits: options.required_part_limits,
        };
        request
            .validate(source)
            .map_err(|error| error.to_string())?;
        Ok(request)
    }
}

#[derive(Clone)]
pub(super) struct Capture {
    pub source: SourceInfo,
    pub epoch: u64,
    pub parent: QcDraft,
    pub form: Form,
    pub request: vector::Request,
    pub revision: String,
}
impl Capture {
    pub fn prepare(
        source: SourceInfo,
        epoch: u64,
        parent: QcDraft,
        form: Form,
        molecule: &Molecule,
    ) -> Result<Self, String> {
        let request = form.request(molecule, &parent)?;
        if source.binding.revision != request.source_revision {
            return Err(
                "The vector-check source changed; select the current document again.".into(),
            );
        }
        let revision = vector::revision(molecule, &request).map_err(|error| error.to_string())?;
        Ok(Self {
            source,
            epoch,
            parent,
            form,
            request,
            revision,
        })
    }
    pub fn guard(
        &self,
        source: &SourceInfo,
        epoch: u64,
        parent: &QcDraft,
        form: &Form,
    ) -> Result<(), String> {
        if self.source != *source
            || self.epoch != epoch
            || self.parent != *parent
            || self.form != *form
        {
            return Err("Source event, workspace or construct-QC/vector declarations changed; run the checks again.".into());
        }
        Ok(())
    }
    pub fn validate_report(&self, report: &vector::Report) -> Result<(), String> {
        report
            .validate_binding(&self.revision)
            .map_err(|error| error.to_string())?;
        if report.request != self.request || report.source_revision != self.source.binding.revision
        {
            return Err(
                "The vector-check response differs from its captured source and declarations."
                    .into(),
            );
        }
        Ok(())
    }
}
