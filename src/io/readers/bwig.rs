use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io;
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::Arc;

use thiserror::Error;
use flate2::Decompress;

use super::block_cache::SharedBlockCache;

const BIGWIG_MAGIC: u32 = 0x888F_FC26;

#[derive(Debug, Default)]
pub struct BigWigReaderStats {
    pub values_calls: u64,
    pub values_returned: u64,
    pub blocks_per_query_total: u64,
    pub cir_cache_hits: u64,
    pub cir_cache_misses: u64,
    pub cir_cache_clears: u64,
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

const MAX_CIR_CACHE_ENTRIES: usize = 1000;

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
pub struct SharedBigWigReader {
    file: File,
    uncompress_buf_size: usize,
    chroms: Vec<ChromInfo>,
    chrom_id_by_name: Vec<(String, u32)>,
    cir_tree_root: u64,
    block_cache: Arc<SharedBlockCache>,
}

impl SharedBigWigReader {
    /// Open a bigWig file with a private block cache (not shared with other files).
    pub fn open(path: impl AsRef<Path>) -> Result<Self, BigWigReadError> {
        Self::open_with_cache(path, Arc::new(SharedBlockCache::new()))
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
        let (chroms, chrom_id_by_name) =
            Self::parse_chrom_tree(&file, chrom_tree_offset)?;

        // The CIR tree root node is at offset + 48
        let cir_tree_root = cir_tree_offset + 48;

        Ok(Self {
            file,
            uncompress_buf_size,
            chroms,
            chrom_id_by_name,
            cir_tree_root,
            block_cache,
        })
    }

    pub fn block_cache(&self) -> &Arc<SharedBlockCache> {
        &self.block_cache
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

        // Read root node header (4 bytes) to determine number of entries
        let root_offset = offset + 32;
        let mut root_header = [0u8; 4];
        pread_exact(file, root_offset, &mut root_header).map_err(BigWigReadError::Io)?;
        let count =
            u16::from_le_bytes([root_header[2], root_header[3]]) as usize;

        let item_size = key_size as usize + val_size as usize;
        let root_size = 4 + count * item_size;
        let mut root_data = vec![0u8; root_size];
        root_data[..4].copy_from_slice(&root_header);
        if root_size > 4 {
            pread_exact(file, root_offset + 4, &mut root_data[4..])
                .map_err(BigWigReadError::Io)?;
        }

        let mut chroms = Vec::new();
        let mut id_by_name = Vec::new();

        for i in 0..count {
            let entry_start = 4 + i * item_size;
            let key = &root_data[entry_start..entry_start + key_size as usize];
            let name_end = key.iter().position(|&b| b == 0).unwrap_or(key.len());
            let name = std::str::from_utf8(&key[..name_end])
                .unwrap_or("")
                .to_string();

            let val_start = entry_start + key_size as usize;
            let chrom_id = u32::from_le_bytes([
                root_data[val_start],
                root_data[val_start + 1],
                root_data[val_start + 2],
                root_data[val_start + 3],
            ]);
            let length = u32::from_le_bytes([
                root_data[val_start + 4],
                root_data[val_start + 5],
                root_data[val_start + 6],
                root_data[val_start + 7],
            ]);

            chroms.push(ChromInfo {
                name: name.clone(),
                length,
            });
            id_by_name.push((name, chrom_id));
        }

        chroms.sort_by(|a, b| a.name.cmp(&b.name));
        id_by_name.sort_by(|a, b| a.0.cmp(&b.0));

        Ok((chroms, id_by_name))
    }
}

// ── Per-worker BigWig reader ──────────────────────────────────────────────
// Each rayon worker gets its own instance with private caches while sharing
// the file handle and parsed metadata via Arc<SharedBigWigReader>.
pub struct BigWigReader {
    shared: Arc<SharedBigWigReader>,
    cir_node_cache: HashMap<u64, Arc<CachedCirNode>>,
    work_buf: Vec<u8>,
    decode_buf: Vec<u8>,
    values_buf: Vec<BigWigValue>,
    blocks_buf: Vec<Block>,
    remaining_buf: VecDeque<u64>,
    pub stats: BigWigReaderStats,
}

impl BigWigReader {
    /// Open a bigWig file, wrapping it in a shared + per-worker pair.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, BigWigReadError> {
        let shared = SharedBigWigReader::open(path)?;
        Ok(Self::from_shared(Arc::new(shared)))
    }

    /// Create a per-worker reader sharing the file handle and metadata of an
    /// already-opened [`SharedBigWigReader`]. Each worker gets its own
    /// CIR-node and block caches.
    pub fn from_shared(shared: Arc<SharedBigWigReader>) -> Self {
        let uncompress_buf_size = shared.uncompress_buf_size;
        Self {
            shared,
            cir_node_cache: HashMap::new(),
            work_buf: Vec::with_capacity(uncompress_buf_size),
            decode_buf: Vec::new(),
            values_buf: Vec::new(),
            blocks_buf: Vec::new(),
            remaining_buf: VecDeque::new(),
            stats: BigWigReaderStats::default(),
        }
    }

    pub fn chroms(&self) -> &[ChromInfo] {
        self.shared.chroms()
    }

    pub fn values(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<&[BigWigValue], BigWigReadError> {
        self.stats.values_calls += 1;

        let chrom_id = self
            .shared
            .find_chrom_id(chrom)
            .ok_or_else(|| BigWigReadError::ChromNotFound(chrom.to_string()))?;

        self.values_buf.clear();
        self.search_cir_tree(chrom_id, start, end)?;

        // Index-based iteration because `self.get_or_cache_block()` borrows
        // &mut self, conflicting with &self.blocks_buf.
        for i in 0..self.blocks_buf.len() {
            let (offset, size) = (self.blocks_buf[i].offset, self.blocks_buf[i].size);
            let data = self.get_or_cache_block(offset, size)?;
            if data.is_empty() {
                continue;
            }
            parse_block_values(&data, start, end, &mut self.values_buf);
        }

        self.stats.values_returned += self.values_buf.len() as u64;
        self.stats.blocks_per_query_total += self.blocks_buf.len() as u64;

        Ok(&self.values_buf)
    }

    fn search_cir_tree(
        &mut self,
        chrom_ix: u32,
        start: u32,
        end: u32,
    ) -> io::Result<()> {
        let file = &self.shared.file;
        let cir_tree_root = self.shared.cir_tree_root;

        self.blocks_buf.clear();
        self.remaining_buf.clear();
        self.remaining_buf.push_front(cir_tree_root);

        // Local stat accumulators to avoid borrow conflicts with split fields.
        let mut cir_hits: u64 = 0;
        let mut cir_misses: u64 = 0;
        let mut cir_clears: u64 = 0;

        while let Some(node_offset) = self.remaining_buf.pop_front() {
            let node = if let Some(cached) = self.cir_node_cache.get(&node_offset) {
                cir_hits += 1;
                Arc::clone(cached)
            } else {
                cir_misses += 1;
                let parsed = SharedBigWigReader::read_cir_node_raw(file, node_offset)?;
                let arc_parsed = Arc::new(parsed);
                if self.cir_node_cache.len() >= MAX_CIR_CACHE_ENTRIES {
                    self.cir_node_cache.clear();
                    cir_clears += 1;
                }
                self.cir_node_cache.insert(node_offset, Arc::clone(&arc_parsed));
                arc_parsed
            };

            for item in &node.items {
                // Overlap check
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

        self.stats.cir_cache_hits += cir_hits;
        self.stats.cir_cache_misses += cir_misses;
        self.stats.cir_cache_clears += cir_clears;

        Ok(())
    }

    fn get_or_cache_block(
        &mut self,
        offset: u64,
        size: u64,
    ) -> io::Result<Arc<[u8]>> {
        let key = (offset, size);
        if let Some(data) = self.shared.block_cache.get(&key) {
            self.stats.block_cache_hits += 1;
            return Ok(data);
        }

        self.stats.block_cache_misses += 1;
        let raw = read_and_decompress(
            &self.shared.file,
            offset,
            size,
            &mut self.work_buf,
            &mut self.decode_buf,
        )?;
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

fn parse_block_values(
    raw: &[u8],
    query_start: u32,
    query_end: u32,
    values: &mut Vec<BigWigValue>,
) {
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
