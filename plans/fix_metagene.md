# Metagene Mode Fix Investigation

## Problem Summary

The metagene test case fails with signal value differences between Python DeepTools and the Rust implementation. The test uses `scale-regions` mode with `--metagene` flag on a GTF file containing two transcripts with multiple exons.

## Status Update (2025-12-01)

- Implemented explicit interval masking for metagene plans so coverage is only filled for the zones that contribute to bins (exons + upstream/downstream), leaving introns as missing data (NaN or 0 when `--missingDataAsZero`).
- Added optional `included_intervals` to `ScaleRegionsPlan` and taught `compute_sample_bins()` to respect it; metagene bins now ignore intronic signal.
- Re-ran python-compatibility `metagene` test: **PASS**, max delta ≈ 1e-6.
- Root cause: intronic bigWig signal was being averaged into bins that spanned exon gaps; Python collapses introns before binning.

## Next Steps

- Run the full python-compatibility suite (all scenarios) after upcoming changes to ensure continued parity.
- Keep an eye on performance; per-interval filling is still cached but now happens per zone when metagene is enabled.

## Test Configuration

```bash
scale-regions -R test.gtf -S test1.bw.bw -a 300 -b 500 --unscaled5prime 20 --unscaled3prime 50 --bs 10 -p 1 --metagene
```

**GTF Data:**
- Transcript 1 (+ strand): exons at (0,50), (400,510), (980,1000) = 180 bp total
- Transcript 2 (- strand): exons at (100,150), (500,610), (1080,1100) = 180 bp total

**Bin Layout (187 total bins):**
- Bins 0-49: upstream (50 bins, 500 bp)
- Bins 50-51: unscaled 5' (2 bins, 20 bp)
- Bins 52-151: body (100 bins, scaled to 1000 bp)
- Bins 152-156: unscaled 3' (5 bins, 50 bp)
- Bins 157-186: downstream (30 bins, 300 bp)

## Observed Differences

### Before Fix

**Row 1 (+ strand):**
- Bin 153: python=27.960000, rust=27.959999 (floating point, acceptable)
- Bin 154: python=26.540000, rust=31.231524 (significant!)
- Bin 159: python=15.450000, rust=15.450001 (floating point, acceptable)

**Row 2 (- strand):**
- Bins 50-52: significant differences (11+ magnitude)
- Bins 152-156: significant differences (9-11 magnitude)

### After First Fix (swap unscaled zone order for - strand)

**Row 1 (+ strand):**
- Bin 154: python=26.540000, rust=31.231524 (still different)

**Row 2 (- strand):**
- Bin 52: python=17.850000, rust=28.162420 (still different, but fewer bins affected)

## Root Cause Analysis

### Python's Metagene Binning Approach

Python's `heatmapper.py` uses a **two-phase** approach:

1. **Signal Collection Phase (`coverage_from_big_wig`):**
   - Fetches base-pair level signal for each exon interval in a zone
   - Concatenates all values into a single contiguous array
   - NaN-pads for out-of-bounds regions

2. **Binning Phase (`coverage_from_array`):**
   - Uses `np.linspace()` to partition the concatenated array into bins
   - Each bin averages values from array indices, NOT genomic coordinates
   - Key insight: bins are computed on the **concatenated exon signal**, ignoring intron gaps

Example for body zone with intervals `[(20, 50), (400, 480)]` (110 bp total, 100 bins):
```python
# Python creates a 110-element array of signal values
values_array = concatenate(signal[20:50], signal[400:480])  # length 110

# Then partitions using linspace
pos_array = np.linspace(0, 110, 100, endpoint=False, dtype=int)  # [0,1,2,...,109]
pos_array = np.append(pos_array, 110)

# Bin 27 covers values_array[29:30] (NOT genomic 29:30!)
# Bin 28 covers values_array[30:31] (first bp of second exon's signal)
```

### Rust's Current Binning Approach

Rust's `zones.rs` uses a **coordinate-based** approach:

1. **Bin Planning Phase (`intervals_to_bins`):**
   - Computes genomic coordinates for each bin
   - Uses `coordinate_from_offset()` to map offsets to genomic positions
   - Produces bins like `[(20, 21), (21, 22), ..., (49, 50), (400, 401), ...]`

2. **Signal Collection Phase (`compute_sample_bins`):**
   - Fetches signal for entire window (including introns)
   - For each bin, extracts `coverage[start_idx:end_idx]` using genomic coordinates

The problem: When bins span exon boundaries, the coordinate mapping is correct, BUT the signal fetching reads from the **genomic coordinate space** which includes intron signal (or NaN for introns without signal).

### Specific Issue: Body Zone Boundary Bins

For the body zone `[(20, 50), (400, 480)]`:
- Total exon length: 110 bp
- Bins needed: 100

At the exon boundary (offsets 29-31):
- Bin 27: offset [29:30] → genomic [49:50] ✓ (last bp of first exon)
- Bin 28: offset [30:31] → genomic [400:401] ✓ (first bp of second exon)
- Bin 29: offset [31:33] → genomic [401:403] ✓

The coordinate mapping is correct! But the issue is that bins like 27 have only 1 bp of signal, while bins like 29 have 2 bp. Python's `linspace` creates similar uneven partitions, but both should produce the same values...

### Deeper Issue: Non-Contiguous Interval Binning

For zones with **non-contiguous intervals** (like unscaled 3' zone `[(480, 510), (980, 1000)]`):

Python concatenates signal from both intervals into a 50-element array, then partitions evenly.

Rust maps bin coordinates across the gap:
- Bin 154: offset [20:30] → should be genomic [500:510] (within first interval)

But wait - let me verify this is actually handled correctly...

## Fix Applied

### Issue 1: Wrong Zone Order and Bin Counts for - Strand

In `build_scale_negative()`, the original code was:
```rust
let (unscaled5, body, unscaled3, _, _) =
    chop_regions(exons, options.unscaled_3_prime, options.unscaled_5_prime);

// WRONG: appending in wrong order with wrong bin counts
append_scale_bins(bins, &unscaled3, unscaled3_bins);  // 20bp with 5 bins
append_scale_bins(bins, &body, body_bins);
append_scale_bins(bins, &unscaled5, unscaled5_bins);  // 50bp with 2 bins
```

Python's zone order for - strand: `[(upstream, 30), (unscaled5prime, 5), (body, 100), (unscaled3prime, 2), (downstream, 50)]`

Where:
- `unscaled5prime` = leftBins from `chopRegions()` = 50 bp (from genomic left = bio 3')
- `unscaled3prime` = rightBins from `chopRegions()` = 20 bp (from genomic right = bio 5')
- `b = unscaled_3_prime // bin_size = 5` bins for 50 bp
- `d = unscaled_5_prime // bin_size = 2` bins for 20 bp

**Fixed to:**
```rust
let (unscaled5, body, unscaled3, _, _) =
    chop_regions(exons, options.unscaled_3_prime, options.unscaled_5_prime);

// CORRECT: match Python's zone order and bin count assignment
append_scale_bins(bins, &unscaled5, unscaled3_bins);  // 50bp with 5 bins
append_scale_bins(bins, &body, body_bins);
append_scale_bins(bins, &unscaled3, unscaled5_bins);  // 20bp with 2 bins
```

## Historical Remaining Issue (resolved)

After the fix, Row 2 (- strand) is mostly correct but:
- Bin 52: python=17.850000, rust=28.162420

And Row 1 (+ strand) still has:
- Bin 154: python=26.540000, rust=31.231524

Both problematic bins are at **exon boundary transitions** in the body/unscaled zones.

## Historical Next Steps / Investigation (for reference)

### Hypothesis (pre-fix)

The remaining differences were likely due to how Rust handled bins that spanned exon boundaries when:
1. The bin coordinate range was correct
2. Signal was fetched for the entire window including introns
3. The index calculation `bin_coord - window_start` picked up intronic signal

### Proposed Fix (pre-fix)

- Restrict signal collection to exon/zone intervals in metagene mode, mirroring Python's concatenated-array binning.

### Files Considered

- `src/pipeline/zones.rs`
- `src/pipeline/core/mod.rs`

## Test Command

```bash
cargo run --release -- scale-regions -R deeptools/deeptools/test/test_data/test.gtf \
  -S deeptools/deeptools/test/test_data/test1.bw.bw \
  -a 300 -b 500 --unscaled5prime 20 --unscaled3prime 50 --bs 10 -p 1 --metagene \
  --outFileName /tmp/test_metagene.mat.gz
```

Compare with reference: `deeptools/deeptools/test/test_heatmapper/master_metagene.mat`
