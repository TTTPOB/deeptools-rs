use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};

use crate::io::{BigWigReader, SharedBigWigReader};

pub struct Sample {
    path: PathBuf,
    reader: BigWigReader,
}

impl Sample {
    pub fn open(path: &Path) -> Result<Self> {
        let reader = BigWigReader::open(path)
            .with_context(|| format!("Failed to open bigWig file '{}'", path.display()))?;

        Ok(Self {
            path: path.to_path_buf(),
            reader,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn reader(&self) -> &BigWigReader {
        &self.reader
    }

    pub fn reader_mut(&mut self) -> &mut BigWigReader {
        &mut self.reader
    }

    /// Create a Sample from an already-opened shared reader.  The underlying
    /// mmap and metadata are shared via Arc; only the per-worker caches are
    /// fresh.
    pub fn from_shared(path: PathBuf, shared: Arc<SharedBigWigReader>) -> Self {
        Self {
            path,
            reader: BigWigReader::from_shared(shared),
        }
    }

    pub fn chrom_length(&self, chrom: &str) -> Option<u32> {
        self.reader.shared().find_chrom_length(chrom)
    }

    /// Read bigWig values using caller-provided decompression buffers.
    pub fn values_with_bufs(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        work_buf: &mut Vec<u8>,
        decode_buf: &mut Vec<u8>,
    ) -> Result<&[crate::io::BigWigValue], anyhow::Error> {
        self.reader
            .values_with_bufs(chrom, start, end, work_buf, decode_buf)
            .map_err(anyhow::Error::new)
    }
}

pub struct WorkerSamples {
    samples: Result<Vec<Sample>, String>,
    work_buf: Vec<u8>,
    decode_buf: Vec<u8>,
}

impl WorkerSamples {
    pub fn new(paths: Arc<Vec<PathBuf>>) -> Self {
        let samples = open_samples(paths.as_ref()).map_err(|err| err.to_string());
        // Compute max uncompress_buf_size across all successfully opened samples
        let max_buf_size = match &samples {
            Ok(s) => s
                .iter()
                .map(|sample| sample.reader().shared().uncompress_buf_size())
                .max()
                .unwrap_or(0),
            Err(_) => 0,
        };
        Self {
            samples,
            work_buf: Vec::with_capacity(max_buf_size),
            decode_buf: Vec::new(),
        }
    }

    /// Create per-worker Sample instances from pre-opened shared readers.
    /// Each worker gets its own caches but shares the mmap-backed immutable
    /// state, avoiding redundant mmap entries per thread.
    pub fn from_shared(
        paths: Arc<Vec<PathBuf>>,
        shared_readers: Arc<Vec<Arc<SharedBigWigReader>>>,
    ) -> Self {
        let max_buf_size = shared_readers
            .iter()
            .map(|r| r.uncompress_buf_size())
            .max()
            .unwrap_or(0);
        let samples = paths
            .iter()
            .zip(shared_readers.iter())
            .map(|(path, shared)| Sample::from_shared(path.clone(), Arc::clone(shared)))
            .collect();
        Self {
            samples: Ok(samples),
            work_buf: Vec::with_capacity(max_buf_size),
            decode_buf: Vec::new(),
        }
    }

    pub fn samples(&mut self) -> Result<&mut Vec<Sample>> {
        match &mut self.samples {
            Ok(samples) => Ok(samples),
            Err(message) => Err(anyhow!(message.clone())),
        }
    }

    /// Return mutable references to both the sample list and the shared
    /// decompression buffers.  This avoids borrow-checker issues that
    /// would arise from separate accessors for samples and buffers.
    pub fn samples_and_bufs(&mut self) -> Result<(&mut Vec<Sample>, &mut Vec<u8>, &mut Vec<u8>)> {
        match &mut self.samples {
            Ok(samples) => Ok((samples, &mut self.work_buf, &mut self.decode_buf)),
            Err(message) => Err(anyhow!(message.clone())),
        }
    }
}

fn open_samples(paths: &[PathBuf]) -> Result<Vec<Sample>> {
    let mut samples = Vec::with_capacity(paths.len());
    for path in paths {
        samples.push(Sample::open(path)?);
    }
    Ok(samples)
}
