use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;

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

pub struct BigWigReader {
    file: File,
    uncompress_buf_size: usize,
    chroms: Vec<ChromInfo>,
    chrom_id_by_name: Vec<(String, u32)>,
    cir_tree_root: u64,
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

fn read_exact_at(file: &mut File, offset: u64, buf: &mut [u8]) -> io::Result<()> {
    file.seek(SeekFrom::Start(offset))?;
    file.read_exact(buf)
}

impl BigWigReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, BigWigReadError> {
        let mut file = File::open(path)?;

        // BBI header (64 bytes)
        let mut hdr = [0u8; 64];
        file.read_exact(&mut hdr)?;
        let mut s = &hdr[..];
        let magic = read_u32(&mut s)?;
        if magic != BIGWIG_MAGIC {
            return Err(BigWigReadError::InvalidMagic(magic));
        }
        let _version = read_u16(&mut s)?;
        let zoom_levels = read_u16(&mut s)?;
        let chrom_tree_offset = read_u64(&mut s)?;
        let _data_offset = read_u64(&mut s)?;
        let cir_tree_offset = read_u64(&mut s)?;
        let _field_count = read_u16(&mut s)?;
        let _defined_field_count = read_u16(&mut s)?;
        let _auto_sql_offset = read_u64(&mut s)?;
        let _total_summary_offset = read_u64(&mut s)?;
        let uncompress_buf_size = read_u32(&mut s)? as usize;

        // Skip zoom headers
        if zoom_levels > 0 {
            file.seek(SeekFrom::Current(zoom_levels as i64 * 24))?;
        }

        // Parse chromosome B+ tree
        let (chroms, chrom_id_by_name) =
            Self::parse_chrom_tree(&mut file, chrom_tree_offset)?;

        // The CIR tree root node is at offset + 48
        let cir_tree_root = cir_tree_offset + 48;

        Ok(Self {
            file,
            uncompress_buf_size,
            chroms,
            chrom_id_by_name,
            cir_tree_root,
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

        // Lazy CIR tree traversal — same algorithm as bigtools' search_cir_tree_inner
        let blocks = self.search_cir_tree(chrom_id, start, end)?;

        let mut values = Vec::new();
        let mut comp_buf = vec![0u8; self.uncompress_buf_size * 2];
        let mut work_buf = vec![0u8; self.uncompress_buf_size];

        for block in &blocks {
            let raw = read_and_decompress(&mut self.file, block.offset, &mut comp_buf, &mut work_buf)?;
            if raw.is_empty() { continue; }
            parse_block_values(raw, start, end, &mut values);
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
            let (child_offsets, node_blocks) =
                Self::process_cir_node(&mut self.file, node_offset, chrom_ix, start, end)?;

            for child in child_offsets.into_iter().rev() {
                remaining.push_front(child);
            }
            blocks.extend(node_blocks);
        }

        Ok(blocks)
    }

    fn process_cir_node(
        file: &mut File,
        offset: u64,
        chrom_ix: u32,
        start: u32,
        end: u32,
    ) -> io::Result<(Vec<u64>, Vec<Block>)> {
        let mut hdr = [0u8; 4];
        read_exact_at(file, offset, &mut hdr)?;
        let is_leaf = hdr[0];
        let count = u16::from_le_bytes([hdr[2], hdr[3]]) as usize;

        // CIR tree: key=16 bytes (id,start,id,end), val=16 bytes (offset,size)
        let item_size = if is_leaf == 0 { 16 + 8 } else { 16 + 16 };
        let total = 4 + count * item_size;
        let mut node_data = vec![0u8; total];
        read_exact_at(file, offset, &mut node_data)?;

        let mut child_offsets = Vec::new();
        let mut blocks = Vec::new();

        for i in 0..count {
            let item_start = 4 + i * item_size;
            let key = &node_data[item_start..item_start + 16];

            let start_chrom_id = u32::from_le_bytes([key[0], key[1], key[2], key[3]]);
            let start_base = u32::from_le_bytes([key[4], key[5], key[6], key[7]]);
            let end_chrom_id = u32::from_le_bytes([key[8], key[9], key[10], key[11]]);
            let end_base = u32::from_le_bytes([key[12], key[13], key[14], key[15]]);

            // Overlap check
            if end_chrom_id < chrom_ix || start_chrom_id > chrom_ix {
                continue;
            }
            if start_chrom_id == end_chrom_id {
                if end_base <= start || start_base >= end {
                    if start_chrom_id == chrom_ix { continue; }
                }
            }

            if is_leaf == 0 {
                let child_off_start = item_start + item_size - 8;
                let child_off = u64::from_le_bytes([
                    node_data[child_off_start], node_data[child_off_start+1],
                    node_data[child_off_start+2], node_data[child_off_start+3],
                    node_data[child_off_start+4], node_data[child_off_start+5],
                    node_data[child_off_start+6], node_data[child_off_start+7],
                ]);
                child_offsets.push(child_off);
            } else {
                let val_start = item_start + 16;
                let data_offset = u64::from_le_bytes([
                    node_data[val_start], node_data[val_start+1],
                    node_data[val_start+2], node_data[val_start+3],
                    node_data[val_start+4], node_data[val_start+5],
                    node_data[val_start+6], node_data[val_start+7],
                ]);
                let data_size = u64::from_le_bytes([
                    node_data[val_start+8], node_data[val_start+9],
                    node_data[val_start+10], node_data[val_start+11],
                    node_data[val_start+12], node_data[val_start+13],
                    node_data[val_start+14], node_data[val_start+15],
                ]);
                blocks.push(Block { offset: data_offset, size: data_size });
            }
        }

        Ok((child_offsets, blocks))
    }

    fn parse_chrom_tree(
        file: &mut File,
        offset: u64,
    ) -> Result<(Vec<ChromInfo>, Vec<(String, u32)>), BigWigReadError> {
        // Read tree header (32 bytes)
        let mut tree_hdr = [0u8; 32];
        read_exact_at(file, offset, &mut tree_hdr)?;
        let mut s = &tree_hdr[..];
        let _magic = read_u32(&mut s)?;
        let _block_size = read_u32(&mut s)?;
        let key_size = read_u32(&mut s)?;
        let val_size = read_u32(&mut s)?;

        // Read root node
        let root_offset = offset + 32;
        let mut node_hdr = [0u8; 4];
        read_exact_at(file, root_offset, &mut node_hdr)?;
        let count = u16::from_le_bytes([node_hdr[2], node_hdr[3]]) as usize;

        let item_size = key_size as usize + val_size as usize;
        let root_size = 4 + count * item_size;
        let mut root_data = vec![0u8; root_size];
        read_exact_at(file, root_offset, &mut root_data)?;

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
    file: &mut File,
    offset: u64,
    comp_buf: &'a mut Vec<u8>,
    work_buf: &'a mut Vec<u8>,
) -> io::Result<&'a [u8]> {
    file.seek(SeekFrom::Start(offset))?;
    let n = file.read(comp_buf)?;
    if n < 2 {
        return Ok(&[]);
    }

    if comp_buf[0] == 0x78 {
        // zlib compressed
        let decoded = DeflateDecoder::new(&comp_buf[..n])
            .decode_zlib()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
        let len = decoded.len();
        if len > work_buf.len() {
            work_buf.resize(len, 0);
        }
        work_buf[..len].copy_from_slice(&decoded);
        Ok(&work_buf[..len])
    } else {
        Ok(&comp_buf[..n])
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
