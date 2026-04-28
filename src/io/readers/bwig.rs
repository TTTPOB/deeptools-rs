use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};
use thiserror::Error;
use zune_inflate::DeflateDecoder;

const BIGWIG_MAGIC: u32 = 0x888F_FC26;
const CHROM_TREE_MAGIC: u32 = 0x78CA_8C91;

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

#[derive(Debug, Clone)]
struct Block {
    offset: u64,
    size: u64,
}

const MAX_BLOCK_CACHE_ENTRIES: usize = 5000;

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

pub struct BigWigReader {
    mmap: Mmap,
    uncompress_buf_size: usize,
    chroms: Vec<ChromInfo>,
    chrom_id_by_name: Vec<(String, u32)>,
    cir_tree_root: u64,
    cir_node_cache: HashMap<u64, CachedCirNode>,
    block_cache: HashMap<(u64, u64), Vec<u8>>,
}

// Little-endian readers from &[u8]
macro_rules! read_le {
    ($buf:ident, $n:expr, $t:ty) => {{
        if $buf.len() < $n { return Err(io::ErrorKind::UnexpectedEof.into()); }
        let v = <$t>::from_le_bytes($buf[..$n].try_into().unwrap());
        *$buf = &$buf[$n..];
        v
    }};
}

fn read_u16(s: &mut &[u8]) -> io::Result<u16> { Ok(read_le!(s, 2, u16)) }
fn read_u32(s: &mut &[u8]) -> io::Result<u32> { Ok(read_le!(s, 4, u32)) }
fn read_u64(s: &mut &[u8]) -> io::Result<u64> { Ok(read_le!(s, 8, u64)) }
fn read_f32(s: &mut &[u8]) -> io::Result<f32> { Ok(read_le!(s, 4, f32)) }

impl BigWigReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, BigWigReadError> {
        let file = File::open(path)?;
        // SAFETY: The file is opened read-only and is not modified while the
        // mmap exists. The mmap is kept alive by BigWigReader for the duration
        // of its lifetime.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        // File can be closed; the mmap keeps the pages alive on Linux.
        drop(file);

        // BBI header (64 bytes) — read from mmap
        let mut s = &mmap[..64];
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

        // Parse chromosome B+ tree from mmap (no seek needed for zoom levels,
        // since all subsequent reads use absolute offsets from the header).
        let (chroms, chrom_id_by_name) =
            Self::parse_chrom_tree(&mmap[..], chrom_tree_offset)?;

        // The CIR tree root node is at offset + 48
        let cir_tree_root = cir_tree_offset + 48;

        Ok(Self {
            mmap,
            uncompress_buf_size,
            chroms,
            chrom_id_by_name,
            cir_tree_root,
            cir_node_cache: HashMap::new(),
            block_cache: HashMap::new(),
        })
    }

    pub fn chroms(&self) -> &[ChromInfo] {
        &self.chroms
    }

    pub fn values(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<BigWigValue>, BigWigReadError> {
        let chrom_id = self.find_chrom_id(chrom)
            .ok_or_else(|| BigWigReadError::ChromNotFound(chrom.to_string()))?;

        let blocks = self.search_cir_tree(chrom_id, start, end)?;

        let mut values = Vec::new();
        let mut work_buf = vec![0u8; self.uncompress_buf_size];

        for block in &blocks {
            let data = self.get_or_cache_block(block.offset, block.size, &mut work_buf)?;
            if data.is_empty() { continue; }
            parse_block_values(&data, start, end, &mut values);
        }

        Ok(values)
    }

    fn find_chrom_id(&self, name: &str) -> Option<u32> {
        self.chrom_id_by_name
            .binary_search_by(|(n, _)| n.as_str().cmp(name))
            .ok()
            .map(|idx| self.chrom_id_by_name[idx].1)
    }

    fn search_cir_tree(&mut self, chrom_ix: u32, start: u32, end: u32) -> io::Result<Vec<Block>> {
        let mut blocks = Vec::new();
        let mut remaining: VecDeque<u64> = VecDeque::with_capacity(2048);
        remaining.push_front(self.cir_tree_root);

        while let Some(node_offset) = remaining.pop_front() {
            let node = if let Some(cached) = self.cir_node_cache.get(&node_offset) {
                cached.clone()
            } else {
                let parsed = Self::read_cir_node_raw(&self.mmap[..], node_offset)?;
                self.cir_node_cache.insert(node_offset, parsed.clone());
                parsed
            };

            for item in &node.items {
                // Overlap check (same as original)
                if item.end_chrom_id < chrom_ix || item.start_chrom_id > chrom_ix {
                    continue;
                }
                if item.start_chrom_id == item.end_chrom_id {
                    if item.end_base <= start || item.start_base >= end {
                        if item.start_chrom_id == chrom_ix { continue; }
                    }
                }

                if node.is_leaf {
                    blocks.push(Block { offset: item.data_offset, size: item.data_size });
                } else {
                    remaining.push_front(item.data_offset);
                }
            }
        }

        Ok(blocks)
    }

    fn read_cir_node_raw(mmap_data: &[u8], offset: u64) -> io::Result<CachedCirNode> {
        let start = offset as usize;
        // Need at least 4 bytes for the node header
        if start + 4 > mmap_data.len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated CIR node header"));
        }

        let is_leaf = mmap_data[start];
        let count = u16::from_le_bytes([mmap_data[start + 2], mmap_data[start + 3]]) as usize;

        let item_size = if is_leaf == 0 { 24 } else { 32 };
        let total = 4 + count * item_size;
        if start + total > mmap_data.len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated CIR node data"));
        }
        let node_data = &mmap_data[start..start + total];

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
                node_data[val_start], node_data[val_start+1],
                node_data[val_start+2], node_data[val_start+3],
                node_data[val_start+4], node_data[val_start+5],
                node_data[val_start+6], node_data[val_start+7],
            ]);
            let data_size = if is_leaf == 0 {
                0
            } else {
                u64::from_le_bytes([
                    node_data[val_start+8], node_data[val_start+9],
                    node_data[val_start+10], node_data[val_start+11],
                    node_data[val_start+12], node_data[val_start+13],
                    node_data[val_start+14], node_data[val_start+15],
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

        Ok(CachedCirNode { is_leaf: is_leaf != 0, items })
    }

    fn get_or_cache_block(
        &mut self,
        offset: u64,
        size: u64,
        work_buf: &mut Vec<u8>,
    ) -> io::Result<Vec<u8>> {
        let key = (offset, size);
        if let Some(data) = self.block_cache.get(&key) {
            return Ok(data.clone());
        }

        let raw = read_and_decompress(&self.mmap[..], offset, size, work_buf)?;
        let data = raw.to_vec();

        if !data.is_empty() {
            if self.block_cache.len() >= MAX_BLOCK_CACHE_ENTRIES {
                self.block_cache.clear();
            }
            self.block_cache.insert(key, data.clone());
        }

        Ok(data)
    }

    fn parse_chrom_tree(
        mmap_data: &[u8],
        offset: u64,
    ) -> Result<(Vec<ChromInfo>, Vec<(String, u32)>), BigWigReadError> {
        let start = offset as usize;
        if start + 32 > mmap_data.len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated chrom tree header").into());
        }
        // Read tree header (32 bytes)
        let mut s = &mmap_data[start..start + 32];
        let _magic = read_u32(&mut s)?;
        let _block_size = read_u32(&mut s)?;
        let key_size = read_u32(&mut s)?;
        let val_size = read_u32(&mut s)?;

        // Read root node
        let root_offset = (offset + 32) as usize;
        if root_offset + 4 > mmap_data.len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated chrom tree root node header").into());
        }
        let count = u16::from_le_bytes([mmap_data[root_offset + 2], mmap_data[root_offset + 3]]) as usize;

        let item_size = key_size as usize + val_size as usize;
        let root_size = 4 + count * item_size;
        if root_offset + root_size > mmap_data.len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated chrom tree root node data").into());
        }
        let root_data = &mmap_data[root_offset..root_offset + root_size];

        let mut chroms = Vec::new();
        let mut id_by_name = Vec::new();

        for i in 0..count {
            let entry_start = 4 + i * item_size;
            let key = &root_data[entry_start..entry_start + key_size as usize];
            let name_end = key.iter().position(|&b| b == 0).unwrap_or(key.len());
            let name = std::str::from_utf8(&key[..name_end]).unwrap_or("").to_string();

            let val_start = entry_start + key_size as usize;
            let chrom_id = u32::from_le_bytes([
                root_data[val_start], root_data[val_start+1],
                root_data[val_start+2], root_data[val_start+3],
            ]);
            let length = u32::from_le_bytes([
                root_data[val_start+4], root_data[val_start+5],
                root_data[val_start+6], root_data[val_start+7],
            ]);

            chroms.push(ChromInfo { name: name.clone(), length });
            id_by_name.push((name, chrom_id));
        }

        chroms.sort_by(|a, b| a.name.cmp(&b.name));
        id_by_name.sort_by(|a, b| a.0.cmp(&b.0));

        Ok((chroms, id_by_name))
    }
}

fn read_and_decompress<'a>(
    mmap_data: &'a [u8],
    offset: u64,
    size: u64,
    work_buf: &'a mut Vec<u8>,
) -> io::Result<&'a [u8]> {
    let start = offset as usize;
    let end = start + size as usize;
    if end > mmap_data.len() || start >= mmap_data.len() {
        return Ok(&[]);
    }
    let block = &mmap_data[start..end];

    if block.is_empty() {
        return Ok(&[]);
    }

    if block[0] == 0x78 {
        // zlib compressed
        let decoded = DeflateDecoder::new(block)
            .decode_zlib()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
        let len = decoded.len();
        if len > work_buf.len() {
            work_buf.resize(len, 0);
        }
        work_buf[..len].copy_from_slice(&decoded);
        Ok(&work_buf[..len])
    } else {
        Ok(block)
    }
}

fn parse_block_values(
    raw: &[u8],
    query_start: u32,
    query_end: u32,
    values: &mut Vec<BigWigValue>,
) {
    let mut s = raw;
    if s.len() < 24 { return; }

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
                if s.len() < 12 { break; }
                let s_val = read_u32(&mut s).unwrap_or(0);
                let e_val = read_u32(&mut s).unwrap_or(0);
                let v = read_f32(&mut s).unwrap_or(f32::NAN);
                if e_val > query_start && s_val < query_end {
                    values.push(BigWigValue { start: s_val, end: e_val, value: v });
                }
            }
        }
        2 => {
            for _ in 0..item_count {
                if s.len() < 8 { break; }
                let s_val = read_u32(&mut s).unwrap_or(0);
                let v = read_f32(&mut s).unwrap_or(f32::NAN);
                let e_val = s_val + item_span;
                if e_val > query_start && s_val < query_end {
                    values.push(BigWigValue { start: s_val, end: e_val, value: v });
                }
            }
        }
        3 => {
            for i in 0..item_count {
                if s.len() < 4 { break; }
                let v = read_f32(&mut s).unwrap_or(f32::NAN);
                let s_val = chrom_start + i as u32 * item_step;
                let e_val = s_val + item_span;
                if e_val > query_start && s_val < query_end {
                    values.push(BigWigValue { start: s_val, end: e_val, value: v });
                }
            }
        }
        _ => {}
    }
}
