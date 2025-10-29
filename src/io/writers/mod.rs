/// Placeholder module for future matrix writers.
///
/// The Rust reimplementation of `computeMatrix` will eventually expose writers
/// mirroring the Python project's serialization helpers. Until the matrix data
/// model stabilizes we keep a minimal placeholder so downstream code can depend
/// on an ergonomic module path without locking in an API prematurely.
#[derive(Debug, Default)]
pub struct WriterPlaceholder;

impl WriterPlaceholder {
    /// Creates a new placeholder instance. Replace this once writer support is implemented.
    pub fn new() -> Self {
        Self::default()
    }
}
