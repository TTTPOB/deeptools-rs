use crate::config::ReferencePoint;
use crate::io::BedRecord;
use crate::io::Strand;

/// Represents a single bin span within the reference-point window.
#[derive(Debug, Clone, Copy)]
pub struct ReferenceBin {
    pub start: i64,
    pub end: i64,
    pub beyond_region: bool,
}

/// Precomputed layout for reference-point binning.
#[derive(Debug, Clone)]
pub struct ReferencePointPlan {
    pub reference: i64,
    pub window_start: i64,
    pub window_end: i64,
    pub bins: Vec<ReferenceBin>,
}

impl ReferencePointPlan {
    pub fn reference_point(
        record: &BedRecord,
        reference_point: ReferencePoint,
        bin_size: u32,
        upstream_bins: usize,
        downstream_bins: usize,
    ) -> Self {
        let reference = reference_coordinate(record, reference_point);
        let bins = build_bins(record, reference, bin_size, upstream_bins, downstream_bins);

        let window_start = bins.iter().map(|bin| bin.start).min().unwrap_or(reference);
        let window_end = bins.iter().map(|bin| bin.end).max().unwrap_or(reference);

        Self {
            reference,
            window_start,
            window_end,
            bins,
        }
    }
}

fn build_bins(
    record: &BedRecord,
    reference: i64,
    bin_size: u32,
    upstream_bins: usize,
    downstream_bins: usize,
) -> Vec<ReferenceBin> {
    let total_bins = upstream_bins + downstream_bins;
    let mut bins = Vec::with_capacity(total_bins);
    for bin_index in 0..total_bins {
        let (start, end) =
            bin_boundaries(bin_index, upstream_bins, bin_size, reference, record.strand);
        let beyond_region = bin_beyond_region(record, start, end);
        bins.push(ReferenceBin {
            start,
            end,
            beyond_region,
        });
    }
    bins
}

fn reference_coordinate(record: &BedRecord, reference_point: ReferencePoint) -> i64 {
    let start = record.start as i64;
    let end = record.end as i64;
    match reference_point {
        ReferencePoint::Tss => match record.strand {
            Strand::Negative => end,
            _ => start,
        },
        ReferencePoint::Tes => match record.strand {
            Strand::Negative => start,
            _ => end,
        },
        ReferencePoint::Center => (start + end) / 2,
    }
}

fn bin_boundaries(
    bin_index: usize,
    upstream_bins: usize,
    bin_size: u32,
    reference: i64,
    strand: Strand,
) -> (i64, i64) {
    let bin_size = bin_size as i64;
    let offset_start = (bin_index as i64 - upstream_bins as i64) * bin_size;
    let offset_end = offset_start + bin_size;

    match strand {
        Strand::Negative => {
            let start = reference - offset_end;
            let end = reference - offset_start;
            if start <= end {
                (start, end)
            } else {
                (end, start)
            }
        }
        _ => {
            let start = reference + offset_start;
            let end = reference + offset_end;
            if start <= end {
                (start, end)
            } else {
                (end, start)
            }
        }
    }
}

fn bin_beyond_region(record: &BedRecord, bin_start: i64, bin_end: i64) -> bool {
    let region_start = record.start as i64;
    let region_end = record.end as i64;

    match record.strand {
        Strand::Negative => bin_end <= region_start,
        _ => bin_start >= region_end,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record(strand: Strand, start: u32, end: u32) -> BedRecord {
        BedRecord {
            chrom: "chr1".to_string(),
            start,
            end,
            name: None,
            score: None,
            strand,
            extra_fields: Vec::new(),
        }
    }

    #[test]
    fn reference_plan_positive_strand() {
        let record = build_record(Strand::Positive, 100, 200);
        let plan = ReferencePointPlan::reference_point(&record, ReferencePoint::Tss, 10, 2, 2);
        assert_eq!(plan.reference, 100);
        assert_eq!(plan.bins.len(), 4);
        assert_eq!(plan.bins[0].start, 80);
        assert_eq!(plan.bins[0].end, 90);
        assert_eq!(plan.bins[3].start, 110);
        assert_eq!(plan.bins[3].end, 120);
    }

    #[test]
    fn reference_plan_negative_strand() {
        let record = build_record(Strand::Negative, 100, 200);
        let plan = ReferencePointPlan::reference_point(&record, ReferencePoint::Tss, 10, 2, 2);
        assert_eq!(plan.reference, 200);
        assert_eq!(plan.bins[0].start, 210);
        assert_eq!(plan.bins[0].end, 220);
        assert_eq!(plan.bins[2].start, 190);
        assert_eq!(plan.bins[2].end, 200);
    }
}
