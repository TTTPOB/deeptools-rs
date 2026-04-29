pub mod auxiliary;
mod formatting;
mod matrix_gz;

pub use matrix_gz::{
    StreamingMatrixWriter, build_padded_header_payload, ensure_streaming_header_capacity,
};
