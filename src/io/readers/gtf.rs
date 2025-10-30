use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result, anyhow};
use bio::io::gff::{GffType, Reader};

use crate::config::GtfOptions;
use crate::io::readers::bed::{BedRecord, Strand as BedStrand};

#[derive(Debug, Clone)]
struct Transcript {
    id: String,
    chrom: String,
    start: u32,
    end: u32,
    strand: BedStrand,
    exons: Vec<(u32, u32)>,
}

pub fn load_gtf_records(path: &Path, options: &GtfOptions) -> Result<Vec<BedRecord>> {
    let mut reader = Reader::from_file(path, GffType::GTF2)
        .with_context(|| format!("Failed to open GTF file '{}'", path.display()))?;
    collect_transcripts(&mut reader, options)
        .with_context(|| format!("Failed to parse GTF file '{}'", path.display()))?
        .into_iter()
        .map(transcript_to_bed)
        .collect()
}

fn collect_transcripts<R: std::io::Read>(
    reader: &mut Reader<R>,
    options: &GtfOptions,
) -> Result<Vec<Transcript>> {
    let mut transcripts = Vec::new();
    let mut index_by_id = HashMap::new();
    let mut pending_exons: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

    for record in reader.records() {
        let record = record.map_err(|err| anyhow!("Encountered invalid GTF record: {err}"))?;

        let feature_type = record.feature_type();
        if feature_type == options.transcript_id {
            let Some(transcript_id) = record
                .attributes()
                .get(&options.transcript_id_designator)
                .map(|value| value.to_string())
            else {
                continue;
            };

            let start = record.start().checked_sub(1).ok_or_else(|| {
                anyhow!("Transcript start was smaller than 1 for {transcript_id}")
            })?;
            let start = u32::try_from(start)
                .with_context(|| format!("Transcript start overflow for {transcript_id}"))?;
            let end = u32::try_from(*record.end())
                .with_context(|| format!("Transcript end overflow for {transcript_id}"))?;

            if end <= start {
                continue;
            }

            let strand_symbol = record.strand().map(|value| value.to_string());
            let strand = match strand_symbol.as_deref() {
                Some("+") => BedStrand::Positive,
                Some("-") => BedStrand::Negative,
                _ => BedStrand::Unstranded,
            };

            let mut transcript = Transcript {
                id: transcript_id.clone(),
                chrom: record.seqname().to_string(),
                start,
                end,
                strand,
                exons: pending_exons.remove(&transcript_id).unwrap_or_default(),
            };

            transcript.exons.sort_by_key(|exon| exon.0);

            index_by_id.insert(transcript_id, transcripts.len());
            transcripts.push(transcript);
            continue;
        }

        if feature_type != options.exon_id {
            continue;
        }

        let Some(parent_id) = record
            .attributes()
            .get(&options.transcript_id_designator)
            .map(|value| value.to_string())
        else {
            continue;
        };

        let start = record
            .start()
            .checked_sub(1)
            .ok_or_else(|| anyhow!("Exon start was smaller than 1 for {parent_id}"))?;
        let start =
            u32::try_from(start).with_context(|| format!("Exon start overflow for {parent_id}"))?;
        let end = u32::try_from(*record.end())
            .with_context(|| format!("Exon end overflow for {parent_id}"))?;

        if end <= start {
            continue;
        }

        let target = if let Some(&index) = index_by_id.get(&parent_id) {
            &mut transcripts[index].exons
        } else {
            pending_exons.entry(parent_id.clone()).or_default()
        };

        target.push((start, end));
    }

    for transcript in &mut transcripts {
        transcript.exons.sort_by_key(|exon| exon.0);
    }

    Ok(transcripts)
}

fn transcript_to_bed(transcript: Transcript) -> Result<BedRecord> {
    let mut extra_fields = Vec::new();

    if !transcript.exons.is_empty() {
        let mut block_sizes = Vec::new();
        let mut block_starts = Vec::new();

        for (start, end) in &transcript.exons {
            let clipped_start = (*start).max(transcript.start);
            let clipped_end = (*end).min(transcript.end);
            if clipped_end <= clipped_start {
                continue;
            }

            block_sizes.push((clipped_end - clipped_start).to_string());
            let offset = clipped_start.saturating_sub(transcript.start);
            block_starts.push(offset.to_string());
        }

        if !block_sizes.is_empty() {
            let thick_start = transcript.start.to_string();
            let thick_end = transcript.end.to_string();
            let item_rgb = "0".to_string();
            let block_count = block_sizes.len().to_string();
            let block_sizes = with_trailing_comma(block_sizes.into_iter());
            let block_starts = with_trailing_comma(block_starts.into_iter());

            extra_fields = vec![
                thick_start,
                thick_end,
                item_rgb,
                block_count,
                block_sizes,
                block_starts,
            ];
        }
    }

    Ok(BedRecord {
        chrom: transcript.chrom,
        start: transcript.start,
        end: transcript.end,
        name: Some(transcript.id),
        score: None,
        strand: transcript.strand,
        extra_fields,
    })
}

fn with_trailing_comma<I>(values: I) -> String
where
    I: Iterator<Item = String>,
{
    let collected: Vec<String> = values.collect();
    if collected.is_empty() {
        return String::new();
    }

    let mut joined = collected.join(",");
    joined.push(',');
    joined
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn options() -> GtfOptions {
        GtfOptions {
            keep_exons: false,
            transcript_id: "transcript".to_string(),
            exon_id: "exon".to_string(),
            transcript_id_designator: "transcript_id".to_string(),
        }
    }

    #[test]
    fn parses_basic_transcript() {
        let data =
            b"chr1\tsrc\ttranscript\t5\t10\t.\t+\t.\ttranscript_id \"tx1\"; gene_id \"g1\";\n";
        let mut reader = Reader::new(Cursor::new(&data[..]), GffType::GTF2);
        let transcripts = collect_transcripts(&mut reader, &options()).expect("collect");
        assert_eq!(transcripts.len(), 1);
        assert_eq!(transcripts[0].chrom, "chr1");
        assert_eq!(transcripts[0].start, 4);
        assert_eq!(transcripts[0].end, 10);
        assert_eq!(transcripts[0].strand, BedStrand::Positive);

        let record = transcript_to_bed(transcripts.into_iter().next().unwrap()).expect("bed");
        assert_eq!(record.name.as_deref(), Some("tx1"));
        assert!(record.extra_fields.is_empty());
    }

    #[test]
    fn attaches_exons_even_if_seen_first() {
        let data = b"chr1\tsrc\texon\t5\t7\t.\t-\t.\ttranscript_id \"tx1\"; gene_id \"g1\";\nchr1\tsrc\ttranscript\t5\t10\t.\t-\t.\ttranscript_id \"tx1\"; gene_id \"g1\";\n";
        let mut reader = Reader::new(Cursor::new(&data[..]), GffType::GTF2);
        let transcripts = collect_transcripts(&mut reader, &options()).expect("collect");
        assert_eq!(transcripts.len(), 1);
        assert_eq!(transcripts[0].strand, BedStrand::Negative);
        assert_eq!(transcripts[0].exons.len(), 1);
        assert_eq!(transcripts[0].exons[0], (4, 7));

        let record = transcript_to_bed(transcripts.into_iter().next().unwrap()).expect("bed");
        assert_eq!(record.extra_fields.len(), 6);
        assert_eq!(record.extra_fields[3], "1"); // block count
        assert_eq!(record.extra_fields[4], "3,");
        assert_eq!(record.extra_fields[5], "0,");
    }

    #[test]
    fn skips_features_without_ids() {
        let data = b"chr1\tsrc\ttranscript\t5\t10\t.\t+\t.\tgene_id \"g1\";\nchr1\tsrc\texon\t5\t7\t.\t+\t.\tgene_id \"g1\";\n";
        let mut reader = Reader::new(Cursor::new(&data[..]), GffType::GTF2);
        let transcripts = collect_transcripts(&mut reader, &options()).expect("collect");
        assert!(transcripts.is_empty());
    }
}
