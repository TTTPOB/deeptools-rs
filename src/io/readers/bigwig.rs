use std::path::Path;

use bigtools::utils::reopen::ReopenableFile;
use bigtools::{BBIReadError, BigWigRead, BigWigReadOpenError, CachedBBIFileRead, Summary};
use thiserror::Error;

pub use bigtools::ChromInfo;

pub type CachedBigWig = BigWigRead<CachedBBIFileRead<ReopenableFile>>;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BigWigValue {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

impl From<bigtools::Value> for BigWigValue {
    fn from(value: bigtools::Value) -> Self {
        Self {
            start: value.start,
            end: value.end,
            value: value.value,
        }
    }
}

#[derive(Debug, Error)]
pub enum BigWigReadError {
    #[error("failed to open bigWig file: {0}")]
    Open(#[from] BigWigReadOpenError),
    #[error("failed to read bigWig interval: {0}")]
    Read(#[from] BBIReadError),
    #[error("failed to summarize bigWig: {0}")]
    Summary(std::io::Error),
}

pub struct BigWigReader {
    inner: CachedBigWig,
}

impl BigWigReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, BigWigReadError> {
        let reader = BigWigRead::open_file(path).map_err(BigWigReadError::Open)?;
        Ok(Self {
            inner: reader.cached(),
        })
    }

    pub fn chroms(&self) -> &[ChromInfo] {
        self.inner.chroms()
    }

    pub fn summary(&mut self) -> Result<Summary, BigWigReadError> {
        self.inner.get_summary().map_err(BigWigReadError::Summary)
    }

    pub fn values(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<BigWigValue>, BigWigReadError> {
        let iterator = self
            .inner
            .get_interval(chrom, start, end)
            .map_err(BigWigReadError::Read)?;
        let mut values = Vec::new();
        for value in iterator {
            let value = value.map_err(BigWigReadError::Read)?;
            values.push(BigWigValue::from(value));
        }
        Ok(values)
    }
}
