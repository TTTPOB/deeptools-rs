use anyhow::Result;

use crate::io::writers::StreamingMatrixWriter;
use crate::pipeline::matrix::{MatrixData, MatrixHeader, MatrixRow};

pub trait RowCollector: Send {
    type Output: Send;

    fn on_row(&mut self, row: MatrixRow) -> Result<()>;
    fn finalize(self, header: MatrixHeader) -> Result<Self::Output>;

    /// Discard this collector without finalising.  The default implementation
    /// simply drops `self`, which is appropriate for in-memory collectors.
    /// File-backed collectors should override this to release I/O resources.
    fn abort(self) where Self: Sized {}
}

pub struct InMemoryCollector {
    rows: Vec<MatrixRow>,
    sample_count: usize,
    bin_count: usize,
}

impl InMemoryCollector {
    pub fn with_capacity(capacity: usize, sample_count: usize, bin_count: usize) -> Self {
        Self {
            rows: Vec::with_capacity(capacity),
            sample_count,
            bin_count,
        }
    }
}

impl RowCollector for InMemoryCollector {
    type Output = MatrixData;

    fn on_row(&mut self, row: MatrixRow) -> Result<()> {
        self.rows.push(row);
        Ok(())
    }

    fn finalize(self, header: MatrixHeader) -> Result<Self::Output> {
        Ok(MatrixData {
            header,
            rows: self.rows,
            bin_count: self.bin_count,
            sample_count: self.sample_count,
        })
    }
}

pub struct FileCollector {
    writer: StreamingMatrixWriter,
}

impl FileCollector {
    pub fn new(writer: StreamingMatrixWriter) -> Self {
        Self { writer }
    }

    /// Discard the underlying writer without finalising the gzip stream.
    pub fn abort(self) {
        self.writer.abort();
    }
}

impl RowCollector for FileCollector {
    type Output = ();

    fn on_row(&mut self, row: MatrixRow) -> Result<()> {
        self.writer.write_row(&row)
    }

    fn finalize(self, header: MatrixHeader) -> Result<Self::Output> {
        self.writer.finish(&header)
    }

    fn abort(self) {
        self.writer.abort();
    }
}

pub struct GroupBucketCollector {
    buckets: Vec<Vec<MatrixRow>>,
    sample_count: usize,
    bin_count: usize,
}

impl GroupBucketCollector {
    pub fn new(group_count: usize, sample_count: usize, bin_count: usize) -> Self {
        Self {
            buckets: (0..group_count).map(|_| Vec::new()).collect(),
            sample_count,
            bin_count,
        }
    }

    pub fn on_row_with_group(&mut self, group_index: usize, row: MatrixRow) -> Result<()> {
        self.buckets[group_index].push(row);
        Ok(())
    }

    pub fn finalize_grouped<F>(self, header_builder: F) -> Result<MatrixData>
    where
        F: FnOnce(Vec<usize>) -> Result<MatrixHeader>,
    {
        let group_counts: Vec<usize> = self.buckets.iter().map(|b| b.len()).collect();
        let total_rows: usize = group_counts.iter().sum();
        let mut rows = Vec::with_capacity(total_rows);
        for bucket in self.buckets {
            rows.extend(bucket);
        }
        let header = header_builder(group_counts)?;
        Ok(MatrixData {
            header,
            rows,
            bin_count: self.bin_count,
            sample_count: self.sample_count,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::readers::bed::{BedRecord, Strand};
    use crate::pipeline::matrix::MatrixHeader;
    use std::sync::Arc;

    fn dummy_row(name: &str, values: Vec<f64>) -> MatrixRow {
        MatrixRow {
            record: BedRecord {
                chrom: Arc::from("chr1"),
                start: 100,
                end: 200,
                name: Some(name.to_string()),
                score: None,
                score_raw: None,
                strand: Strand::Unstranded,
                strand_raw: None,
                extra_fields: vec![],
            },
            values,
            sample_count: 1,
            bin_count: 3,
            exon_coords: None,
        }
    }

    #[test]
    fn group_bucket_collector_groups_rows_by_index() {
        let mut collector = GroupBucketCollector::new(3, 1, 3);
        collector.on_row_with_group(1, dummy_row("b", vec![1.0, 2.0, 3.0])).unwrap();
        collector.on_row_with_group(0, dummy_row("a", vec![4.0, 5.0, 6.0])).unwrap();
        collector.on_row_with_group(2, dummy_row("c", vec![7.0, 8.0, 9.0])).unwrap();
        collector.on_row_with_group(0, dummy_row("d", vec![10.0, 11.0, 12.0])).unwrap();

        let data = collector.finalize_grouped(|counts| {
            Ok(MatrixHeader::default_for_test(counts))
        }).unwrap();

        assert_eq!(data.rows.len(), 4);
        // Group 0: a, d. Group 1: b. Group 2: c.
        assert_eq!(data.rows[0].record.name.as_deref(), Some("a"));
        assert_eq!(data.rows[1].record.name.as_deref(), Some("d"));
        assert_eq!(data.rows[2].record.name.as_deref(), Some("b"));
        assert_eq!(data.rows[3].record.name.as_deref(), Some("c"));
        assert_eq!(data.header.group_boundaries, vec![0, 2, 3, 4]);
    }
}
