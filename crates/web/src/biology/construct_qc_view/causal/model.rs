//! Capture explicit CDS assumptions and exact event identity before native replay.
use super::selection::{self, EventChoice, Form};
use collagene_core::{
    construct_qc::{self, causal},
    documents::workspace::Workspace,
    material::lineage::{self, EventBinding},
    sequence::translation,
};

fn assumption(
    workspace: &Workspace,
    event: &EventBinding,
    index: &str,
    rationale: &str,
) -> Result<construct_qc::AnnotationAssumption, String> {
    let index = index
        .parse::<usize>()
        .map_err(|_| "Choose a CDS annotation.")?;
    let molecule = selection::molecule(workspace, event)?;
    let feature = molecule
        .features
        .get(index)
        .ok_or("The CDS annotation is unavailable.")?;
    if !feature.kind.eq_ignore_ascii_case("CDS") {
        return Err("Choose a CDS annotation, rather than another feature type.".into());
    }
    Ok(construct_qc::AnnotationAssumption {
        feature_index: index,
        source_revision: event.revision.clone(),
        annotation_revision: construct_qc::annotation_revision(feature)
            .map_err(|error| error.to_string())?,
        rationale: rationale.to_owned(),
    })
}
impl Form {
    pub fn request(
        &self,
        workspace: &Workspace,
        choices: &[EventChoice],
        product: &EventBinding,
    ) -> Result<causal::Request, String> {
        if self.rationale.trim().is_empty()
            || self.rationale.len() > 1024
            || self.rationale.chars().any(char::is_control)
        {
            return Err("Describe why these source and product annotations represent the CDS roles being compared (at most 1024 bytes).".into());
        }
        let source = self.source(choices)?.clone();
        let genetic_code = self
            .genetic_code
            .parse::<u8>()
            .map_err(|_| "Choose a genetic code.")?;
        translation::CodeDescriptor::for_id(genetic_code).map_err(|error| error.to_string())?;
        Ok(causal::Request {
            schema_version: 1,
            limits: causal::Limits::default(),
            pairs: vec![causal::CdsPair {
                id: "selected-cds".into(),
                source_cds: assumption(workspace, &source, &self.source_feature, &self.rationale)?,
                product_cds: assumption(
                    workspace,
                    product,
                    &self.product_feature,
                    &self.rationale,
                )?,
                source,
                product: product.clone(),
                genetic_code,
            }],
        })
    }
}
#[derive(Clone)]
pub(super) struct Capture {
    pub product: EventBinding,
    pub epoch: u64,
    pub form: Form,
    pub request: causal::Request,
    pub revision: String,
    pub fingerprint: String,
}
impl Capture {
    pub fn prepare(
        workspace: &Workspace,
        choices: &[EventChoice],
        product: EventBinding,
        epoch: u64,
        form: Form,
    ) -> Result<Self, String> {
        let request = form.request(workspace, choices, &product)?;
        let revision = causal::revision(workspace, &request).map_err(|error| error.to_string())?;
        let fingerprint =
            lineage::workspace_fingerprint(workspace).map_err(|error| error.to_string())?;
        Ok(Self {
            product,
            epoch,
            form,
            request,
            revision,
            fingerprint,
        })
    }
    pub fn guard(
        &self,
        workspace: &Workspace,
        product: &EventBinding,
        epoch: u64,
        form: &Form,
    ) -> Result<(), String> {
        if self.product != *product
            || self.epoch != epoch
            || self.form != *form
            || lineage::workspace_fingerprint(workspace).map_err(|error| error.to_string())?
                != self.fingerprint
        {
            return Err("The workspace, selected event or CDS declarations changed; run the comparison again.".into());
        }
        Ok(())
    }
    pub fn validate_report(
        &self,
        workspace: &Workspace,
        report: &causal::Report,
    ) -> Result<(), String> {
        // Documentary/request validation only: native scientific replay stays in
        // the worker. Use the canonical owner for limits and pair/code bindings.
        report
            .validate_request(workspace, &self.request)
            .map_err(|error| error.to_string())?;
        if report.context_revision != self.revision
            || report.workspace_fingerprint != self.fingerprint
        {
            return Err("The comparison response belongs to another captured workspace.".into());
        }
        Ok(())
    }
}
