pub mod traits;
pub mod samples;
pub mod regions;
mod executor;
mod collector;
mod coalesce;
mod worker;

pub use traits::{
    SignalBin, ModeTag, RegionPlan, PipelineMode,
    ensure_positive, ensure_multiple,
};
pub use samples::{Sample, WorkerSamples};
pub use regions::{Group, RegionTask, load_groups, derive_sample_labels, normalize_sort_sample_indices};
pub use collector::{RowCollector, InMemoryCollector, FileCollector, GroupBucketCollector};
pub use executor::execute_mode;
pub use worker::compute_row;

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::io::readers::bed::Strand;
    use crate::io::BedRecord;

    use super::executor::{WorkItem, input_order_is_compute_sorted};
    use super::traits::{RegionPlan, SignalBin};
    use super::worker::{aggregate_slice, index_from_coordinate};
    use crate::config::AverageTypeBins;

    fn work_item(idx: usize, chrom: &str, start: i64, end: i64) -> WorkItem {
        WorkItem {
            orig_idx: idx,
            group_index: 0,
            record: Arc::new(BedRecord {
                chrom: Arc::from(chrom),
                start: start as u32,
                end: end as u32,
                name: None,
                score: None,
                score_raw: None,
                strand: Strand::Unstranded,
                strand_raw: None,
                extra_fields: vec![],
            }),
            query_start: start,
            query_end: end,
        }
    }

    #[test]
    fn test_input_is_compute_sorted_true() {
        let items = vec![
            work_item(0, "chr1", 100, 200),
            work_item(1, "chr1", 300, 400),
            work_item(2, "chr2", 50, 150),
        ];
        assert!(input_order_is_compute_sorted(&items));
    }

    #[test]
    fn test_input_is_compute_sorted_false() {
        let items = vec![
            work_item(0, "chr2", 100, 200),
            work_item(1, "chr1", 300, 400),
        ];
        assert!(!input_order_is_compute_sorted(&items));
    }

    #[test]
    fn test_input_is_compute_sorted_empty() {
        let items: Vec<WorkItem> = vec![];
        assert!(input_order_is_compute_sorted(&items));
    }

    #[test]
    fn test_input_is_compute_sorted_single() {
        let items = vec![work_item(0, "chr1", 100, 200)];
        assert!(input_order_is_compute_sorted(&items));
    }

    #[derive(Clone)]
    struct TestBin {
        start: i64,
        end: i64,
        beyond_region: bool,
    }

    impl SignalBin for TestBin {
        fn start(&self) -> i64 {
            self.start
        }

        fn end(&self) -> i64 {
            self.end
        }

        fn beyond_region(&self) -> bool {
            self.beyond_region
        }
    }

    struct TestPlan {
        start: i64,
        end: i64,
        bins: Vec<TestBin>,
    }

    impl RegionPlan for TestPlan {
        type Bin = TestBin;

        fn window_start(&self) -> i64 {
            self.start
        }

        fn window_end(&self) -> i64 {
            self.end
        }

        fn bins(&self) -> &[Self::Bin] {
            &self.bins
        }
    }

    #[test]
    fn index_from_coordinate_bounds_checks() {
        let base = 100;
        let window_len = 50;
        assert_eq!(index_from_coordinate(90, base, window_len), 0);
        assert_eq!(index_from_coordinate(100, base, window_len), 0);
        assert_eq!(index_from_coordinate(125, base, window_len), 25);
        assert_eq!(index_from_coordinate(200, base, window_len), window_len);
    }

    #[test]
    fn aggregate_slice_ignores_nans() {
        let data = [1.0, f64::NAN, 3.0, 5.0];
        let mean = aggregate_slice(&data, AverageTypeBins::Mean).unwrap();
        assert!((mean - 3.0).abs() < 1e-6);

        let max = aggregate_slice(&data, AverageTypeBins::Max).unwrap();
        assert_eq!(max, 5.0);

        let median = aggregate_slice(&data, AverageTypeBins::Median).unwrap();
        assert_eq!(median, 3.0);
    }

    #[test]
    fn aggregate_mean_f64_precision() {
        let val = 35.92f64;
        let data = vec![val; 9];
        let mean = aggregate_slice(&data, AverageTypeBins::Mean).unwrap();
        assert!(
            (mean - val).abs() < 1e-10,
            "f64 mean drift too large: expected {val}, got {mean}, delta {}",
            (mean - val).abs()
        );
    }
}
