use anyhow::Result;

use crate::io::writers::StreamingMatrixWriter;
use crate::pipeline::matrix::{MatrixData, MatrixHeader, MatrixRow};

pub trait RowCollector: Send {
    type Output: Send;

    fn on_row(&mut self, row: MatrixRow) -> Result<()>;
    fn finalize(self, header: MatrixHeader) -> Result<Self::Output>;
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
}

impl RowCollector for FileCollector {
    type Output = ();

    fn on_row(&mut self, row: MatrixRow) -> Result<()> {
        self.writer.write_row(&row)
    }

    fn finalize(self, header: MatrixHeader) -> Result<Self::Output> {
        self.writer.finish(&header)
    }
}
