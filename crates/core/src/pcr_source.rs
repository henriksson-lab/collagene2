//! Borrowed native PCR frames. No synthesized reference stock or restored report.
use crate::{
    analysis::{self, PcrGeometryProduct, Primer, TemplateView},
    annotations::MappingSegment,
    chemistry::{self, Strand},
    modified_duplex::TrackedDuplex,
};
use collagene_model::{Error, Molecule, Sequence, Topology};
#[derive(Clone, Copy)]
pub(crate) struct Template<'a> {
    pub reference: &'a Molecule,
    pub sequence: &'a Sequence,
    pub topology: Topology,
    pub strand: Strand,
    tracked: Option<&'a TrackedDuplex>,
}
impl<'a> Template<'a> {
    pub(crate) fn reference(reference: &'a Molecule) -> Result<Self, Error> {
        reference.validate()?;
        if chemistry::has_unmapped_physical_evidence(reference) {
            return Err(Error::Unsupported(
                "PCR of physical-intermediate evidence requires its fresh native origin".into(),
            ));
        }
        Ok(Self {
            reference,
            sequence: &reference.sequence,
            topology: reference.topology,
            strand: Strand::Top,
            tracked: None,
        })
    }
    pub(crate) fn fresh(
        reference: &'a Molecule,
        tracked: &'a TrackedDuplex,
        strand: Strand,
    ) -> Result<Self, Error> {
        reference.validate()?;
        tracked.validate()?;
        if tracked.duplex.top != reference.sequence
            || tracked.duplex.topology != reference.topology
            || !tracked.duplex.nicks.is_empty()
        {
            return Err(Error::Unsupported("PCR requires exact reference-matching intact native strands; nick propagation is not modeled".into()));
        }
        let sequence = match strand {
            Strand::Top => &tracked.duplex.top,
            Strand::Bottom => &tracked.duplex.bottom,
        };
        Ok(Self {
            reference,
            sequence,
            topology: tracked.duplex.topology,
            strand,
            tracked: Some(tracked),
        })
    }
    pub(crate) fn with_strand(self, strand: Strand) -> Result<Self, Error> {
        match self.tracked {
            Some(tracked) => Self::fresh(self.reference, tracked, strand),
            None if strand == Strand::Top => Ok(self),
            None => Err(Error::Unsupported(
                "Bottom preparation needs actual freshly reconstructed native strands".into(),
            )),
        }
    }
    pub(crate) fn view(self) -> TemplateView<'a> {
        TemplateView {
            sequence: self.sequence,
            topology: self.topology,
            template_strand: self.strand,
        }
    }
    pub(crate) fn context_revision(self) -> Result<String, Error> {
        crate::assembly_recipe::hash(&(
            "native-pcr-frame-1",
            self.reference.revision(),
            self.strand,
            self.tracked,
        ))
    }
    pub(crate) fn geometries_with_usage(
        self,
        forward: &Primer,
        reverse: &Primer,
    ) -> Result<(Vec<PcrGeometryProduct>, usize), Error> {
        let work = self
            .comparison_bases()?
            .checked_mul(
                forward
                    .annealing
                    .len()
                    .checked_add(reverse.annealing.len())
                    .ok_or(Error::TooLarge)?,
            )
            .and_then(|n| n.checked_mul(4))
            .ok_or(Error::TooLarge)?;
        if work > 100_000_000 {
            return Err(Error::TooLarge);
        }
        let mut materialized = 0usize;
        if let Some(tracked) = self.tracked {
            let d = &tracked.duplex;
            if d.bottom_start != 0 || d.bottom.len() != d.top.len() {
                if d.topology != Topology::Linear {
                    return Err(Error::Unsupported(
                        "Staggered circular PCR requires an explicit supported linearization"
                            .into(),
                    ));
                }
                let bottom_end = d
                    .bottom_start
                    .checked_add(i64::try_from(d.bottom.len()).map_err(|_| Error::TooLarge)?)
                    .ok_or(Error::TooLarge)?;
                let start = 0i64.min(d.bottom_start);
                let end = (d.top.len() as i64).max(bottom_end);
                let length = usize::try_from(end.checked_sub(start).ok_or(Error::TooLarge)?)
                    .map_err(|_| Error::TooLarge)?;
                if length > collagene_model::MAX_SEQUENCE_LEN {
                    return Err(Error::TooLarge);
                }
                // Charge diagnostic source copies too, not just emitted PCR
                // products: complement, coverage and selected-frame orientation.
                materialized = materialized
                    .checked_add(d.bottom.len())
                    .and_then(|n| n.checked_add(length))
                    .and_then(|n| {
                        n.checked_add(if self.strand == Strand::Bottom {
                            length
                        } else {
                            0
                        })
                    })
                    .filter(|n| *n <= 64 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                let bottom_plus = d.bottom.reverse_complement()?;
                let mut bases = Vec::with_capacity(length);
                for i in 0..length {
                    let coordinate = start.checked_add(i as i64).ok_or(Error::TooLarge)?;
                    let base = if coordinate >= 0 && (coordinate as usize) < d.top.len() {
                        d.top.as_str().as_bytes()[coordinate as usize]
                    } else {
                        let p = coordinate
                            .checked_sub(d.bottom_start)
                            .and_then(|n| usize::try_from(n).ok())
                            .filter(|p| *p < bottom_plus.len())
                            .ok_or_else(|| {
                                Error::Unsupported(
                                    "PCR source coverage has an unsupported covalent gap".into(),
                                )
                            })?;
                        bottom_plus.as_str().as_bytes()[p]
                    };
                    bases.push(base);
                }
                // Coverage is a diagnostic search only. It never seeds polymerase or
                // becomes an input document: actual products use self.sequence below.
                let coverage = Sequence::dna(
                    std::str::from_utf8(&bases)
                        .map_err(|_| Error::Format("PCR coverage is not ASCII".into()))?,
                )?;
                let oriented = if self.strand == Strand::Bottom {
                    coverage.reverse_complement()?
                } else {
                    coverage
                };
                let offset = usize::try_from(match self.strand {
                    Strand::Top => start.checked_neg().ok_or(Error::TooLarge)?,
                    Strand::Bottom => end.checked_sub(bottom_end).ok_or(Error::TooLarge)?,
                })
                .map_err(|_| Error::TooLarge)?;
                for geometry in analysis::simulate_pcr_geometry_view(
                    TemplateView {
                        sequence: &oriented,
                        topology: Topology::Linear,
                        template_strand: self.strand,
                    },
                    forward,
                    reverse,
                )? {
                    materialized = materialized
                        .checked_add(geometry.product.sequence.len())
                        .ok_or(Error::TooLarge)?;
                    let begin = geometry.product.forward_start;
                    let stop = begin
                        .checked_add(geometry.template_length)
                        .ok_or(Error::TooLarge)?;
                    if begin < offset
                        || stop
                            > offset
                                .checked_add(self.sequence.len())
                                .ok_or(Error::TooLarge)?
                    {
                        return Err(Error::Unsupported("Actual native strands admit an additional PCR geometry outside the selected intact strand frame; no unique product, absence or implicit repair is asserted".into()));
                    }
                }
            }
        }
        let products = analysis::simulate_pcr_geometry_view(self.view(), forward, reverse)?;
        for product in &products {
            materialized = materialized
                .checked_add(product.product.sequence.len())
                .ok_or(Error::TooLarge)?;
        }
        if materialized > 64 * 1024 * 1024 {
            return Err(Error::TooLarge);
        }
        Ok((products, materialized))
    }
    pub(crate) fn geometries(
        self,
        forward: &Primer,
        reverse: &Primer,
    ) -> Result<Vec<PcrGeometryProduct>, Error> {
        self.geometries_with_usage(forward, reverse)
            .map(|(products, _)| products)
    }
    pub(crate) fn comparison_bases(self) -> Result<usize, Error> {
        let mut n = self.sequence.len();
        if let Some(t) = self.tracked {
            let d = &t.duplex;
            if d.bottom_start != 0 || d.bottom.len() != d.top.len() {
                let end = d
                    .bottom_start
                    .checked_add(d.bottom.len() as i64)
                    .ok_or(Error::TooLarge)?;
                let union = (d.top.len() as i64)
                    .max(end)
                    .checked_sub(0i64.min(d.bottom_start))
                    .ok_or(Error::TooLarge)?;
                n = n
                    .checked_add(usize::try_from(union).map_err(|_| Error::TooLarge)?)
                    .ok_or(Error::TooLarge)?;
            }
        }
        Ok(n)
    }
    pub(crate) fn annotation_map(
        self,
        geometry: &PcrGeometryProduct,
    ) -> Result<Vec<MappingSegment>, Error> {
        if geometry.product.template_strand != self.strand {
            return Err(Error::Format(
                "PCR geometry strand differs from its native frame".into(),
            ));
        }
        let n = self.sequence.len();
        let start = geometry.product.forward_start;
        let length = geometry.template_length;
        if start >= n
            || length == 0
            || length > n
            || self.topology == Topology::Linear
                && start.checked_add(length).is_none_or(|end| end > n)
        {
            return Err(Error::Coordinates(
                "PCR geometry exceeds selected native strand".into(),
            ));
        }
        let mut result = vec![];
        let first = length.min(n - start);
        let pieces = [(start, start + first, 0usize), (0, length - first, first)];
        for (a, b, output) in pieces {
            if a == b {
                continue;
            }
            let destination = geometry
                .forward_tail_length
                .checked_add(output)
                .ok_or(Error::TooLarge)?;
            if self.strand == Strand::Top {
                result.push(MappingSegment {
                    source_start: a,
                    source_end: b,
                    destination_start: destination,
                    reverse: false,
                });
            } else {
                let d = &self
                    .tracked
                    .ok_or_else(|| {
                        Error::Unsupported(
                            "Bottom PCR annotation coordinates require native source state".into(),
                        )
                    })?
                    .duplex;
                let end = d
                    .bottom_start
                    .checked_add(d.bottom.len() as i64)
                    .ok_or(Error::TooLarge)?;
                let lo = (end - b as i64)
                    .max(0)
                    .min(self.reference.sequence.len() as i64);
                let hi = (end - a as i64)
                    .max(0)
                    .min(self.reference.sequence.len() as i64);
                if lo < hi {
                    result.push(MappingSegment {
                        source_start: lo as usize,
                        source_end: hi as usize,
                        destination_start: destination
                            .checked_add(
                                usize::try_from(end - hi - a as i64)
                                    .map_err(|_| Error::TooLarge)?,
                            )
                            .ok_or(Error::TooLarge)?,
                        reverse: true,
                    });
                }
            }
        }
        Ok(result)
    }
}
