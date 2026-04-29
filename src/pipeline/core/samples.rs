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
}

pub struct WorkerSamples {
    samples: Result<Vec<Sample>, String>,
}

impl WorkerSamples {
    pub fn new(paths: Arc<Vec<PathBuf>>) -> Self {
        let samples = open_samples(paths.as_ref()).map_err(|err| err.to_string());
        Self { samples }
    }

    /// Create per-worker Sample instances from pre-opened shared readers.
    /// Each worker gets its own caches but shares the mmap-backed immutable
    /// state, avoiding redundant mmap entries per thread.
    pub fn from_shared(
        paths: Arc<Vec<PathBuf>>,
        shared_readers: Arc<Vec<Arc<SharedBigWigReader>>>,
    ) -> Self {
        let samples = paths
            .iter()
            .zip(shared_readers.iter())
            .map(|(path, shared)| Sample::from_shared(path.clone(), Arc::clone(shared)))
            .collect();
        Self {
            samples: Ok(samples),
        }
    }

    pub fn samples(&mut self) -> Result<&mut Vec<Sample>> {
        match &mut self.samples {
            Ok(samples) => Ok(samples),
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
