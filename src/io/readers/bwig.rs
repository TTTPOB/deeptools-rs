use std::collections::VecDeque;
use std::fs::File;
use std::io;
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::Arc;

use flate2::Decompress;
use thiserror::Error;

use quick_cache::sync::Cache;

use super::block_cache::SharedBlockCache;

const BIGWIG_MAGIC: u32 = 0x888F_FC26;

#[derive(Debug, Default)]
pub struct BigWigReaderStats {
    pub values_calls: u64,
    pub values_returned: u64,
    pub blocks_per_query_total: u64,
    pub block_cache_hits: u64,
    pub block_cache_misses: u64,
    pub decoded_bytes: u64,
}

#[derive(Debug, Clone)]
pub struct ChromInfo {
    pub name: String,
    pub length: u32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BigWigValue {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

#[derive(Debug, Error)]
pub enum BigWigReadError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("invalid bigWig magic: {0:#X}")]
    InvalidMagic(u32),
    #[error("decompression failed: {0}")]
    Decompression(String),
    #[error("chromosome '{0}' not found")]
    ChromNotFound(String),
}

#[derive(Debug, Clone, Copy)]
struct Block {
    offset: u64,
    size: u64,
}

const SHARED_CIR_CACHE_ENTRIES: usize = 1000;

#[derive(Debug, Clone)]
struct CirNodeItem {
    start_chrom_id: u32,
    start_base: u32,
    end_chrom_id: u32,
    end_base: u32,
    data_offset: u64,
    data_size: u64,
}

#[derive(Debug, Clone)]
struct CachedCirNode {
    is_leaf: bool,
    items: Vec<CirNodeItem>,
}

// Little-endian readers from &[u8]
macro_rules! read_le {
    ($buf:ident, $n:expr, $t:ty) => {{
        if $buf.len() < $n {
            return Err(io::ErrorKind::UnexpectedEof.into());
        }
        let v = <$t>::from_le_bytes($buf[..$n].try_into().unwrap());
        *$buf = &$buf[$n..];
        v
    }};
}

fn read_u16(s: &mut &[u8]) -> io::Result<u16> {
    Ok(read_le!(s, 2, u16))
}
fn read_u32(s: &mut &[u8]) -> io::Result<u32> {
    Ok(read_le!(s, 4, u32))
}
fn read_u64(s: &mut &[u8]) -> io::Result<u64> {
    Ok(read_le!(s, 8, u64))
}
fn read_f32(s: &mut &[u8]) -> io::Result<f32> {
    Ok(read_le!(s, 4, f32))
}

/// Helper: read exactly `buf.len()` bytes at `offset` via pread.
fn pread_exact(file: &File, offset: u64, buf: &mut [u8]) -> io::Result<()> {
    file.read_exact_at(buf, offset)
}

// ── Shared immutable state: File + parsed metadata ───────────────────────
// One instance per bigWig file; shared across threads via Arc.

/// Binary-search the sorted chroms slice for the given chromosome name and
/// return its length, or None if not present.
fn binary_search_chrom_length(chroms: &[ChromInfo], name: &str) -> Option<u32> {
    chroms
        .binary_search_by(|c| c.name.as_str().cmp(name))
        .ok()
        .map(|idx| chroms[idx].length)
}

pub struct BigWigFile {
    file: File,
    uncompress_buf_size: usize,
    chroms: Vec<ChromInfo>,
    chrom_id_by_name: Vec<(String, u32)>,
    cir_tree_root: u64,
    block_cache: Arc<SharedBlockCache>,
    cir_cache: Cache<u64, Arc<CachedCirNode>>,
}

impl BigWigFile {
    /// Open a bigWig file with a per-file block cache of the given capacity.
    pub fn open_with_block_cache_capacity(
        path: impl AsRef<Path>,
        block_cache_capacity: usize,
    ) -> Result<Self, BigWigReadError> {
        Self::open_with_cache(
            path,
            Arc::new(SharedBlockCache::with_capacity(block_cache_capacity)),
        )
    }

    /// Open a bigWig file using the provided shared block cache.
    pub fn open_with_cache(
        path: impl AsRef<Path>,
        block_cache: Arc<SharedBlockCache>,
    ) -> Result<Self, BigWigReadError> {
        let file = File::open(path)?;

        // BBI header (64 bytes)
        let mut header_buf = [0u8; 64];
        pread_exact(&file, 0, &mut header_buf)?;
        let mut s = &header_buf[..];
        let magic = read_u32(&mut s)?;
        if magic != BIGWIG_MAGIC {
            return Err(BigWigReadError::InvalidMagic(magic));
        }
        let _version = read_u16(&mut s)?;
        let _zoom_levels = read_u16(&mut s)?;
        let chrom_tree_offset = read_u64(&mut s)?;
        let _data_offset = read_u64(&mut s)?;
        let cir_tree_offset = read_u64(&mut s)?;
        let _field_count = read_u16(&mut s)?;
        let _defined_field_count = read_u16(&mut s)?;
        let _auto_sql_offset = read_u64(&mut s)?;
        let _total_summary_offset = read_u64(&mut s)?;
        let uncompress_buf_size = read_u32(&mut s)? as usize;

        // Parse chromosome B+ tree from file.
        let (chroms, chrom_id_by_name) = Self::parse_chrom_tree(&file, chrom_tree_offset)?;

        // The CIR tree root node is at offset + 48
        let cir_tree_root = cir_tree_offset + 48;

        Ok(Self {
            file,
            uncompress_buf_size,
            chroms,
            chrom_id_by_name,
            cir_tree_root,
            block_cache,
            cir_cache: Cache::new(SHARED_CIR_CACHE_ENTRIES),
        })
    }

    pub fn chroms(&self) -> &[ChromInfo] {
        &self.chroms
    }

    pub fn find_chrom_id(&self, name: &str) -> Option<u32> {
        self.chrom_id_by_name
            .binary_search_by(|(n, _)| n.as_str().cmp(name))
            .ok()
            .map(|idx| self.chrom_id_by_name[idx].1)
    }

    pub fn find_chrom_length(&self, name: &str) -> Option<u32> {
        binary_search_chrom_length(&self.chroms, name)
    }

    pub fn uncompress_buf_size(&self) -> usize {
        self.uncompress_buf_size
    }

    /// Look up a CIR tree node from the shared cache, or read and parse it from
    /// the file and insert it into the cache.
    fn get_or_read_cir_node(&self, offset: u64) -> io::Result<Arc<CachedCirNode>> {
        if let Some(node) = self.cir_cache.get(&offset) {
            return Ok(node);
        }
        let parsed = Self::read_cir_node_raw(&self.file, offset)?;
        let arc_node = Arc::new(parsed);
        self.cir_cache.insert(offset, Arc::clone(&arc_node));
        Ok(arc_node)
    }

    /// Read and parse a CIR tree node from the file at the given offset.
    fn read_cir_node_raw(file: &File, offset: u64) -> io::Result<CachedCirNode> {
        // Read node header (4 bytes) to determine size
        let mut header = [0u8; 4];
        file.read_exact_at(&mut header, offset)?;

        let is_leaf = header[0];
        let count = u16::from_le_bytes([header[2], header[3]]) as usize;
        let item_size = if is_leaf == 0 { 24 } else { 32 };
        let total = 4 + count * item_size;

        // Read the full node data
        let mut node_data = vec![0u8; total];
        node_data[..4].copy_from_slice(&header);
        if total > 4 {
            file.read_exact_at(&mut node_data[4..], offset + 4)?;
        }

        let mut items = Vec::with_capacity(count);
        for i in 0..count {
            let item_start = 4 + i * item_size;
            let key = &node_data[item_start..item_start + 16];

            let start_chrom_id = u32::from_le_bytes([key[0], key[1], key[2], key[3]]);
            let start_base = u32::from_le_bytes([key[4], key[5], key[6], key[7]]);
            let end_chrom_id = u32::from_le_bytes([key[8], key[9], key[10], key[11]]);
            let end_base = u32::from_le_bytes([key[12], key[13], key[14], key[15]]);

            let val_start = item_start + 16;
            let data_offset = u64::from_le_bytes([
                node_data[val_start],
                node_data[val_start + 1],
                node_data[val_start + 2],
                node_data[val_start + 3],
                node_data[val_start + 4],
                node_data[val_start + 5],
                node_data[val_start + 6],
                node_data[val_start + 7],
            ]);
            let data_size = if is_leaf == 0 {
                0
            } else {
                u64::from_le_bytes([
                    node_data[val_start + 8],
                    node_data[val_start + 9],
                    node_data[val_start + 10],
                    node_data[val_start + 11],
                    node_data[val_start + 12],
                    node_data[val_start + 13],
                    node_data[val_start + 14],
                    node_data[val_start + 15],
                ])
            };

            items.push(CirNodeItem {
                start_chrom_id,
                start_base,
                end_chrom_id,
                end_base,
                data_offset,
                data_size,
            });
        }

        Ok(CachedCirNode {
            is_leaf: is_leaf != 0,
            items,
        })
    }

    /// Parse the chromosome B+ tree from the file.
    ///
    /// Supports multi-level B+ trees by iteratively traversing internal
    /// (non-leaf) nodes via a stack, collecting chromosome entries only
    /// from leaf nodes.
    fn parse_chrom_tree(
        file: &File,
        offset: u64,
    ) -> Result<(Vec<ChromInfo>, Vec<(String, u32)>), BigWigReadError> {
        // Read tree header (32 bytes)
        let mut header = [0u8; 32];
        pread_exact(file, offset, &mut header).map_err(BigWigReadError::Io)?;
        let mut s = &header[..];
        let _magic = read_u32(&mut s)?;
        let _block_size = read_u32(&mut s)?;
        let key_size = read_u32(&mut s)?;
        let val_size = read_u32(&mut s)?;

        let root_offset = offset + 32;
        let mut chroms = Vec::new();
        let mut id_by_name = Vec::new();

        // Stack-based traversal to support multi-level B+ trees.
        // For single-level trees (the common case) this processes only
        // the root node, matching the previous behaviour.
        let mut stack = vec![root_offset];
        while let Some(node_offset) = stack.pop() {
            // Read node header (4 bytes): isLeaf(1) + reserved(1) + count(2 LE)
            let mut node_header = [0u8; 4];
            pread_exact(file, node_offset, &mut node_header).map_err(BigWigReadError::Io)?;
            let is_leaf = node_header[0] != 0;
            let count = u16::from_le_bytes([node_header[2], node_header[3]]) as usize;

            if is_leaf {
                // Leaf node: each entry is keySize + valSize bytes
                let item_size = key_size as usize + val_size as usize;
                let node_size = 4 + count * item_size;
                let mut node_data = vec![0u8; node_size];
                node_data[..4].copy_from_slice(&node_header);
                if node_size > 4 {
                    pread_exact(file, node_offset + 4, &mut node_data[4..])
                        .map_err(BigWigReadError::Io)?;
                }

                for i in 0..count {
                    let entry_start = 4 + i * item_size;
                    let key = &node_data[entry_start..entry_start + key_size as usize];
                    let name_end = key.iter().position(|&b| b == 0).unwrap_or(key.len());
                    let name = std::str::from_utf8(&key[..name_end])
                        .unwrap_or("")
                        .to_string();

                    let val_start = entry_start + key_size as usize;
                    let chrom_id = u32::from_le_bytes([
                        node_data[val_start],
                        node_data[val_start + 1],
                        node_data[val_start + 2],
                        node_data[val_start + 3],
                    ]);
                    let length = u32::from_le_bytes([
                        node_data[val_start + 4],
                        node_data[val_start + 5],
                        node_data[val_start + 6],
                        node_data[val_start + 7],
                    ]);

                    chroms.push(ChromInfo {
                        name: name.clone(),
                        length,
                    });
                    id_by_name.push((name, chrom_id));
                }
            } else {
                // Internal node: each entry is keySize + 8 bytes (child offset)
                let item_size = key_size as usize + 8;
                let node_size = 4 + count * item_size;
                let mut node_data = vec![0u8; node_size];
                node_data[..4].copy_from_slice(&node_header);
                if node_size > 4 {
                    pread_exact(file, node_offset + 4, &mut node_data[4..])
                        .map_err(BigWigReadError::Io)?;
                }

                for i in 0..count {
                    let entry_start = 4 + i * item_size;
                    let child_offset_start = entry_start + key_size as usize;
                    let child_offset = u64::from_le_bytes(
                        node_data[child_offset_start..child_offset_start + 8]
                            .try_into()
                            .unwrap(),
                    );
                    stack.push(child_offset);
                }
            }
        }

        chroms.sort_by(|a, b| a.name.cmp(&b.name));
        id_by_name.sort_by(|a, b| a.0.cmp(&b.0));

        Ok((chroms, id_by_name))
    }
}

// ── Per-worker BigWig reader ──────────────────────────────────────────────
// Each rayon worker gets its own instance with private caches while sharing
// the file handle and parsed metadata via Arc<BigWigFile>.
pub struct BigWigReader {
    shared: Arc<BigWigFile>,
    values_buf: Vec<BigWigValue>,
    blocks_buf: Vec<Block>,
    remaining_buf: VecDeque<u64>,
    pub stats: BigWigReaderStats,
}

impl BigWigReader {
    /// Create a per-worker reader sharing the file handle and metadata of an
    /// already-opened [`BigWigFile`]. Each worker gets its own
    /// CIR-node and block caches.
    pub fn from_shared(shared: Arc<BigWigFile>) -> Self {
        Self {
            shared,
            values_buf: Vec::new(),
            blocks_buf: Vec::new(),
            remaining_buf: VecDeque::new(),
            stats: BigWigReaderStats::default(),
        }
    }

    pub fn shared(&self) -> &Arc<BigWigFile> {
        &self.shared
    }

    /// Read bigWig values for the given region, allocating temporary
    /// decompression buffers internally.  Prefer `values_with_bufs` on the
    /// hot path to reuse externally-owned buffers across calls.
    pub fn values(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<&[BigWigValue], BigWigReadError> {
        let mut work_buf = Vec::with_capacity(self.shared.uncompress_buf_size);
        let mut decode_buf = Vec::new();
        self.values_with_bufs(chrom, start, end, &mut work_buf, &mut decode_buf)
    }

    /// Read bigWig values using caller-provided decompression buffers.
    /// This avoids per-call allocation when buffers are reused across
    /// multiple samples and batches (the common case in `process_batch`).
    pub fn values_with_bufs(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        work_buf: &mut Vec<u8>,
        decode_buf: &mut Vec<u8>,
    ) -> Result<&[BigWigValue], BigWigReadError> {
        self.stats.values_calls += 1;

        let chrom_id = self
            .shared
            .find_chrom_id(chrom)
            .ok_or_else(|| BigWigReadError::ChromNotFound(chrom.to_string()))?;

        self.values_buf.clear();
        self.search_cir_tree(chrom_id, start, end)?;

        for i in 0..self.blocks_buf.len() {
            let (offset, size) = (self.blocks_buf[i].offset, self.blocks_buf[i].size);
            let data = self.get_or_cache_block_with_bufs(offset, size, work_buf, decode_buf)?;
            if data.is_empty() {
                continue;
            }
            parse_block_values(&data, start, end, &mut self.values_buf);
        }

        self.stats.values_returned += self.values_buf.len() as u64;
        self.stats.blocks_per_query_total += self.blocks_buf.len() as u64;

        Ok(&self.values_buf)
    }

    fn search_cir_tree(&mut self, chrom_ix: u32, start: u32, end: u32) -> io::Result<()> {
        let cir_tree_root = self.shared.cir_tree_root;

        self.blocks_buf.clear();
        self.remaining_buf.clear();
        self.remaining_buf.push_front(cir_tree_root);

        while let Some(node_offset) = self.remaining_buf.pop_front() {
            let node = self.shared.get_or_read_cir_node(node_offset)?;

            for item in &node.items {
                if item.end_chrom_id < chrom_ix || item.start_chrom_id > chrom_ix {
                    continue;
                }
                if item.start_chrom_id == item.end_chrom_id {
                    if item.end_base <= start || item.start_base >= end {
                        if item.start_chrom_id == chrom_ix {
                            continue;
                        }
                    }
                }

                if node.is_leaf {
                    self.blocks_buf.push(Block {
                        offset: item.data_offset,
                        size: item.data_size,
                    });
                } else {
                    self.remaining_buf.push_front(item.data_offset);
                }
            }
        }

        Ok(())
    }

    fn get_or_cache_block_with_bufs(
        &mut self,
        offset: u64,
        size: u64,
        work_buf: &mut Vec<u8>,
        decode_buf: &mut Vec<u8>,
    ) -> io::Result<Arc<[u8]>> {
        let key = (offset, size);
        if let Some(data) = self.shared.block_cache.get(&key) {
            self.stats.block_cache_hits += 1;
            return Ok(data);
        }

        self.stats.block_cache_misses += 1;
        let raw = read_and_decompress(&self.shared.file, offset, size, work_buf, decode_buf)?;
        self.stats.decoded_bytes += raw.len() as u64;
        let data: Arc<[u8]> = Arc::from(raw);

        if !data.is_empty() {
            self.shared.block_cache.insert(key, Arc::clone(&data));
        }

        Ok(data)
    }
}

fn read_and_decompress<'a>(
    file: &File,
    offset: u64,
    size: u64,
    read_buf: &mut Vec<u8>,
    decode_buf: &'a mut Vec<u8>,
) -> io::Result<&'a [u8]> {
    if size == 0 {
        return Ok(&[]);
    }

    let buf_len = size as usize;
    if buf_len > read_buf.len() {
        read_buf.resize(buf_len, 0);
    }
    file.read_exact_at(&mut read_buf[..buf_len], offset)?;

    let block = &read_buf[..buf_len];
    if block.is_empty() {
        return Ok(&[]);
    }

    if block[0] == 0x78 {
        // zlib compressed — decompress into reusable decode_buf
        decode_buf.clear();
        let capacity_hint = buf_len * 4;
        if decode_buf.capacity() < capacity_hint {
            decode_buf.reserve(capacity_hint - decode_buf.capacity());
        }
        let mut decoder = Decompress::new(true);
        decoder
            .decompress_vec(block, decode_buf, flate2::FlushDecompress::Finish)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
        Ok(decode_buf.as_slice())
    } else {
        // Uncompressed — copy into decode_buf so we can return a consistent lifetime
        decode_buf.clear();
        decode_buf.extend_from_slice(&read_buf[..buf_len]);
        Ok(decode_buf.as_slice())
    }
}

fn parse_block_values(raw: &[u8], query_start: u32, query_end: u32, values: &mut Vec<BigWigValue>) {
    let mut s = raw;
    if s.len() < 24 {
        return;
    }

    let _chrom_id = read_u32(&mut s).unwrap_or(0);
    let chrom_start = read_u32(&mut s).unwrap_or(0);
    let chrom_end = read_u32(&mut s).unwrap_or(0);
    let item_step = read_u32(&mut s).unwrap_or(0);
    let item_span = read_u32(&mut s).unwrap_or(0);
    let block_type = s.first().copied().unwrap_or(0);
    s = &s[2..];
    let item_count = read_u16(&mut s).unwrap_or(0);

    if chrom_end <= query_start || chrom_start >= query_end {
        return;
    }

    match block_type {
        1 => {
            for _ in 0..item_count {
                if s.len() < 12 {
                    break;
                }
                let s_val = read_u32(&mut s).unwrap_or(0);
                let e_val = read_u32(&mut s).unwrap_or(0);
                let v = read_f32(&mut s).unwrap_or(f32::NAN);
                if e_val > query_start && s_val < query_end {
                    values.push(BigWigValue {
                        start: s_val,
                        end: e_val,
                        value: v,
                    });
                }
            }
        }
        2 => {
            for _ in 0..item_count {
                if s.len() < 8 {
                    break;
                }
                let s_val = read_u32(&mut s).unwrap_or(0);
                let v = read_f32(&mut s).unwrap_or(f32::NAN);
                let e_val = s_val + item_span;
                if e_val > query_start && s_val < query_end {
                    values.push(BigWigValue {
                        start: s_val,
                        end: e_val,
                        value: v,
                    });
                }
            }
        }
        3 => {
            for i in 0..item_count {
                if s.len() < 4 {
                    break;
                }
                let v = read_f32(&mut s).unwrap_or(f32::NAN);
                let s_val = chrom_start + i as u32 * item_step;
                let e_val = s_val + item_span;
                if e_val > query_start && s_val < query_end {
                    values.push(BigWigValue {
                        start: s_val,
                        end: e_val,
                        value: v,
                    });
                }
            }
        }
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── binary_search_chrom_length ────────────────────────────────────────────

    #[test]
    fn binary_search_chrom_length_found() {
        let chroms = vec![
            ChromInfo {
                name: "chr1".to_string(),
                length: 1000,
            },
            ChromInfo {
                name: "chr2".to_string(),
                length: 2000,
            },
            ChromInfo {
                name: "chrX".to_string(),
                length: 3000,
            },
        ];
        assert_eq!(binary_search_chrom_length(&chroms, "chr1"), Some(1000));
        assert_eq!(binary_search_chrom_length(&chroms, "chr2"), Some(2000));
        assert_eq!(binary_search_chrom_length(&chroms, "chrX"), Some(3000));
    }

    #[test]
    fn binary_search_chrom_length_not_found() {
        let chroms = vec![
            ChromInfo {
                name: "chr1".to_string(),
                length: 1000,
            },
            ChromInfo {
                name: "chr2".to_string(),
                length: 2000,
            },
        ];
        assert_eq!(binary_search_chrom_length(&chroms, "chr3"), None);
        assert_eq!(binary_search_chrom_length(&chroms, ""), None);
    }

    #[test]
    fn binary_search_chrom_length_empty_vec() {
        let chroms: Vec<ChromInfo> = vec![];
        assert_eq!(binary_search_chrom_length(&chroms, "chr1"), None);
    }

    // Build a 24-byte block header with the given parameters.
    fn build_header(
        chrom_start: u32,
        chrom_end: u32,
        item_step: u32,
        item_span: u32,
        block_type: u8,
        item_count: u16,
    ) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&0u32.to_le_bytes()); // chrom_id (ignored)
        buf.extend_from_slice(&chrom_start.to_le_bytes());
        buf.extend_from_slice(&chrom_end.to_le_bytes());
        buf.extend_from_slice(&item_step.to_le_bytes());
        buf.extend_from_slice(&item_span.to_le_bytes());
        buf.push(block_type);
        buf.push(0); // reserved
        buf.extend_from_slice(&item_count.to_le_bytes());
        buf
    }

    // ── block type 1 (bedGraph) ──────────────────────────────────────────────

    #[test]
    fn bedgraph_normal_case() {
        // Single item fully within query range [100, 300)
        let mut raw = build_header(0, 1000, 0, 0, 1, 1);
        raw.extend_from_slice(&100u32.to_le_bytes()); // start
        raw.extend_from_slice(&200u32.to_le_bytes()); // end
        raw.extend_from_slice(&1.5f32.to_le_bytes()); // value

        let mut values = Vec::new();
        parse_block_values(&raw, 100, 300, &mut values);

        assert_eq!(values.len(), 1);
        assert_eq!(values[0].start, 100);
        assert_eq!(values[0].end, 200);
        assert!((values[0].value - 1.5).abs() < 1e-6);
    }

    #[test]
    fn bedgraph_items_partially_overlapping_query() {
        // Three items: before, overlapping, after
        let mut raw = build_header(0, 1000, 0, 0, 1, 3);
        // Item 1: completely before query [500, 700) → should NOT appear
        raw.extend_from_slice(&100u32.to_le_bytes());
        raw.extend_from_slice(&200u32.to_le_bytes());
        raw.extend_from_slice(&1.0f32.to_le_bytes());
        // Item 2: overlaps query → should appear
        raw.extend_from_slice(&600u32.to_le_bytes());
        raw.extend_from_slice(&650u32.to_le_bytes());
        raw.extend_from_slice(&2.0f32.to_le_bytes());
        // Item 3: completely after query → should NOT appear
        raw.extend_from_slice(&800u32.to_le_bytes());
        raw.extend_from_slice(&900u32.to_le_bytes());
        raw.extend_from_slice(&3.0f32.to_le_bytes());

        let mut values = Vec::new();
        parse_block_values(&raw, 500, 700, &mut values);

        assert_eq!(values.len(), 1);
        assert_eq!(values[0].start, 600);
    }

    #[test]
    fn bedgraph_items_outside_query_range() {
        // Items exist but query range does not overlap any of them.
        let mut raw = build_header(0, 1000, 0, 0, 1, 2);
        raw.extend_from_slice(&100u32.to_le_bytes());
        raw.extend_from_slice(&200u32.to_le_bytes());
        raw.extend_from_slice(&1.0f32.to_le_bytes());
        raw.extend_from_slice(&300u32.to_le_bytes());
        raw.extend_from_slice(&400u32.to_le_bytes());
        raw.extend_from_slice(&2.0f32.to_le_bytes());

        let mut values = Vec::new();
        parse_block_values(&raw, 500, 600, &mut values);

        assert!(values.is_empty());
    }

    #[test]
    fn bedgraph_empty_block() {
        // item_count = 0, so no values should be pushed.
        let raw = build_header(0, 1000, 0, 0, 1, 0);

        let mut values = Vec::new();
        parse_block_values(&raw, 0, 1000, &mut values);

        assert!(values.is_empty());
    }

    // ── block type 2 (variableStep) ─────────────────────────────────────────

    #[test]
    fn variable_step_normal_case() {
        // item_span = 50; one item at start=200 → end=250, within query [100, 300)
        let mut raw = build_header(0, 1000, 0, 50, 2, 1);
        raw.extend_from_slice(&200u32.to_le_bytes()); // start
        raw.extend_from_slice(&4.0f32.to_le_bytes()); // value

        let mut values = Vec::new();
        parse_block_values(&raw, 100, 300, &mut values);

        assert_eq!(values.len(), 1);
        assert_eq!(values[0].start, 200);
        assert_eq!(values[0].end, 250);
        assert!((values[0].value - 4.0).abs() < 1e-6);
    }

    #[test]
    fn variable_step_query_filtering() {
        // item_span = 10; two items, only one overlaps query [500, 600)
        let mut raw = build_header(0, 1000, 0, 10, 2, 2);
        // Item 1: [200, 210) → outside query
        raw.extend_from_slice(&200u32.to_le_bytes());
        raw.extend_from_slice(&1.0f32.to_le_bytes());
        // Item 2: [550, 560) → inside query
        raw.extend_from_slice(&550u32.to_le_bytes());
        raw.extend_from_slice(&5.0f32.to_le_bytes());

        let mut values = Vec::new();
        parse_block_values(&raw, 500, 600, &mut values);

        assert_eq!(values.len(), 1);
        assert_eq!(values[0].start, 550);
        assert_eq!(values[0].end, 560);
    }

    // ── block type 3 (fixedStep) ─────────────────────────────────────────────

    #[test]
    fn fixed_step_normal_case() {
        // chrom_start=100, item_step=50, item_span=50, 3 items
        // Positions: [100,150), [150,200), [200,250) — all within query [0, 300)
        let mut raw = build_header(100, 250, 50, 50, 3, 3);
        raw.extend_from_slice(&1.0f32.to_le_bytes());
        raw.extend_from_slice(&2.0f32.to_le_bytes());
        raw.extend_from_slice(&3.0f32.to_le_bytes());

        let mut values = Vec::new();
        parse_block_values(&raw, 0, 300, &mut values);

        assert_eq!(values.len(), 3);
        assert_eq!(values[0].start, 100);
        assert_eq!(values[0].end, 150);
        assert!((values[0].value - 1.0).abs() < 1e-6);
        assert_eq!(values[1].start, 150);
        assert_eq!(values[2].start, 200);
    }

    #[test]
    fn fixed_step_query_filtering() {
        // chrom_start=0, item_step=100, item_span=100, 5 items
        // Positions: [0,100), [100,200), [200,300), [300,400), [400,500)
        // Query [150, 300) overlaps items at [100,200) and [200,300) only:
        //   [100,200): e_val=200 > 150 and s_val=100 < 300 → included
        //   [200,300): e_val=300 > 150 and s_val=200 < 300 → included
        //   [300,400): s_val=300 is NOT < query_end=300 → excluded
        let mut raw = build_header(0, 500, 100, 100, 3, 5);
        for v in 1..=5u32 {
            raw.extend_from_slice(&(v as f32).to_le_bytes());
        }

        let mut values = Vec::new();
        parse_block_values(&raw, 150, 300, &mut values);

        assert_eq!(values.len(), 2);
        assert_eq!(values[0].start, 100);
        assert_eq!(values[0].end, 200);
        assert_eq!(values[1].start, 200);
        assert_eq!(values[1].end, 300);
    }

    // ── edge cases ───────────────────────────────────────────────────────────

    #[test]
    fn raw_data_too_short() {
        // Less than 24 bytes → early return, no values pushed.
        let raw = vec![0u8; 10];
        let mut values = Vec::new();
        parse_block_values(&raw, 0, 1000, &mut values);
        assert!(values.is_empty());
    }

    #[test]
    fn block_entirely_before_query() {
        // chrom_end (200) <= query_start (500) → early return
        let raw = build_header(100, 200, 0, 0, 1, 1);
        // Append one valid bedGraph item that would otherwise match
        let mut raw = raw;
        raw.extend_from_slice(&100u32.to_le_bytes());
        raw.extend_from_slice(&150u32.to_le_bytes());
        raw.extend_from_slice(&9.0f32.to_le_bytes());

        let mut values = Vec::new();
        parse_block_values(&raw, 500, 1000, &mut values);
        assert!(values.is_empty());
    }

    #[test]
    fn block_entirely_after_query() {
        // chrom_start (800) >= query_end (500) → early return
        let raw = build_header(800, 900, 0, 0, 1, 1);
        let mut raw = raw;
        raw.extend_from_slice(&800u32.to_le_bytes());
        raw.extend_from_slice(&850u32.to_le_bytes());
        raw.extend_from_slice(&9.0f32.to_le_bytes());

        let mut values = Vec::new();
        parse_block_values(&raw, 0, 500, &mut values);
        assert!(values.is_empty());
    }

    #[test]
    fn unknown_block_type() {
        // block_type = 99 → falls through to `_ => {}`, no values pushed
        let raw = build_header(0, 1000, 0, 0, 99, 1);
        let mut values = Vec::new();
        parse_block_values(&raw, 0, 1000, &mut values);
        assert!(values.is_empty());
    }
}
