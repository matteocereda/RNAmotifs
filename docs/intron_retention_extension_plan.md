# Extending RNAmotifs/RNAmaRs for Intron Retention Events

**Status: Implemented.** All steps below have been completed and tested with PTBP1 HepG2 RI data.

## Motivation

RNAmotifs currently handles only cassette exon (skipped exon, SE) events. Intron retention (IR) is an increasingly recognized mode of splicing regulation, particularly in neuronal differentiation, cancer, and stress response. Extending RNAmaRs to score IR events would broaden its utility for identifying RBPs driving splicing changes in phenotypic comparisons where IR is the dominant splicing change.

## Biological model

For cassette exons, the four analysis regions are:
- **R1**: upstream intron (regulatory elements promoting/repressing exon inclusion)
- **R2**: exon body (exonic splicing enhancers/silencers)
- **R3**: downstream intron (regulatory elements downstream of the exon)
- **R4**: distal flanking region

For intron retention, the regulatory landscape is fundamentally different. The key regulatory elements determining whether an intron is retained or spliced out reside **within the intron itself** — branch point sequences, polypyrimidine tracts, and intronic splicing regulatory elements. The proposed region mapping for IR is:

- **R1**: upstream flanking intron (upstream of the 5' splice site)
- **R2**: first half of the retained intron (5' splice site to midpoint)
- **R3**: second half of the retained intron (midpoint to 3' splice site)
- **R4**: downstream flanking intron (downstream of the 3' splice site)

This places the core regulatory regions (R2, R3) within the retained intron, where the branch point (typically 18-40 nt upstream of the 3' splice site, in R3) and polypyrimidine tract are located. R1 and R4 capture exonic splicing regulatory elements in the flanking exons near each splice site.

## Parameter reuse: `in_exon` and `in_intron`

The existing `in_exon` (default 30 bp) and `in_intron` (default 300 bp) parameters are reused with their natural semantic meaning:

- `in_exon` → extent into **exonic** sequence = R1 and R4 for IR (flanking exons near splice sites)
- `in_intron` → extent into **intronic** sequence = R2 and R3 for IR (into the retained intron from each splice site)

This means for IR events:
- **R1**: last `in_exon` bp of the upstream exon (approaching the 5' splice site)
- **R2**: first `min(in_intron, half_intron)` bp of the retained intron from the 5' splice site
- **R3**: last `min(in_intron, half_intron)` bp of the retained intron before the 3' splice site (branch point, polypyrimidine tract)
- **R4**: first `in_exon` bp of the downstream exon (after the 3' splice site)

For short retained introns (< 2 × `in_intron`), R2 and R3 meet at the midpoint, providing full intron coverage. No new parameters are needed — the same discovery mode parameter grid (`in_exon`, `in_intron`, `hw`, `ew`) tests the same biological question for both SE and IR events: how far from splice sites do regulatory elements extend?

## rMATS RI coordinate format

The rMATS `RI.MATS.JunctionCountOnly.txt` output defines IR events with:
- `upstreamES`, `upstreamEE` — upstream exon (5' exon flanking the retained intron)
- `downstreamES`, `downstreamEE` — downstream exon (3' exon flanking the retained intron)
- The retained intron spans from `upstreamEE` to `downstreamES`

## Implementation plan

### Step 1: Extend `data/prepare_mars_exons.py` for RI events

Add `parse_rmats_ri()` function alongside existing `parse_rmats_se()`.

**Input coordinate mapping** — encode IR events in the RNAmotifs format as:
```
row_id;second_id;chr;strand;upstreamES;upstreamEE;downstreamES;downstreamEE;dIRank;RI
```

The 10th field (`;RI`) marks the event type for the C++ code.

**dIRank classification for IR** — following the ENCODE convention where `IncLevelDifference = KD - Control`:
- `delta_PSI > 0` → retention increases upon knockdown → RBP normally suppresses retention → **silenced** (dIRank < -1)
- `delta_PSI < 0` → retention decreases upon knockdown → RBP promotes retention → **enhanced** (dIRank > 1)

This is consistent with the SE convention where positive dIRank means the event (inclusion for SE, retention for IR) is promoted by the RBP.

Add `--event-type {SE,RI,both}` CLI argument. Skip UCSC cassetteExon filter for RI events.

### Step 2: Add event_type to C++ Config and Exon structs

**File**: `src/cpp/rnamotifs_core.h`

Add `string event_type = "SE"` to the `Config` struct. The event type is auto-detected from the optional 10th field in each input line, or set globally via a CLI flag.

**File**: `src/cpp/rnamotifs_mars_score.cpp`

Add `string event_type = "SE"` to the `Exon` struct (line 41). Parse from the 10th field in `read_input_exons()`.

### Step 3: IR-specific region computation in `rnamotifs_core.cpp`

**File**: `src/cpp/rnamotifs_core.cpp`, function `process_tetramer_file()` (lines 305-410)

Currently, region boundaries are computed as:
```cpp
unsigned int left_intron_len  = in_start - skip_start;
unsigned int exon_len         = in_stop  - in_start;
unsigned int right_intron_len = skip_stop - in_stop;

unsigned int ie_actual = min(exon_len/2, in_exon);    // extent into "exon"
unsigned int ii_left   = min(left_intron_len/2, in_intron);
unsigned int ii_right  = min(right_intron_len/2, in_intron);
```

For IR events, the "exon body" between `v6` and `v7` is actually the retained intron. The `in_exon` and `in_intron` parameters swap their directional roles relative to SE: `in_exon` controls extent into flanking exons (R1, R4) and `in_intron` controls extent into the retained intron (R2, R3). Add an IR-specific code path:

```cpp
if (event_type == "RI") {
    // v5=upstreamES, v6=upstreamEE (5'SS), v7=downstreamES (3'SS), v8=downstreamEE
    unsigned int retained_intron_len = in_stop - in_start;  // v7 - v6
    unsigned int half_intron = retained_intron_len / 2;
    unsigned int upstream_exon_len = in_start - skip_start;   // v6 - v5
    unsigned int downstream_exon_len = skip_stop - in_stop;   // v8 - v7

    // Extent into retained intron from each splice site (capped at half-intron)
    unsigned int ii_actual = min(half_intron, ii);  // ii = in_intron parameter

    // Extent into flanking exons (capped at half-exon length)
    unsigned int ie_left  = min(upstream_exon_len / 2, ie);   // ie = in_exon parameter
    unsigned int ie_right = min(downstream_exon_len / 2, ie);

    // R1: upstream exon near 5'SS [v6 - ie_left, v6]
    r1s = in_start - ie_left;
    r1e = in_start;

    // R2: retained intron from 5'SS [v6, v6 + ii_actual]
    r2s = in_start;
    r2e = in_start + ii_actual;

    // R3: retained intron toward 3'SS [v7 - ii_actual, v7]
    r3s = in_stop - ii_actual;
    r3e = in_stop;

    // R4: downstream exon near 3'SS [v7, v7 + ie_right]
    r4s = in_stop;
    r4e = in_stop + ie_right;
}
```

For short retained introns (< 2 × `in_intron`), `ii_actual = half_intron` and R2/R3 meet exactly at the midpoint, covering the full intron. For long introns, R2 and R3 each scan `in_intron` bp from their respective splice sites, focusing on the regulatory hotspots (5'SS motifs in R2, branch point and polypyrimidine tract in R3).

The `map_to_splicing_map()` call and BED overlap logic remain unchanged — they just operate on different region boundaries.

### Step 4: IR-aware `count_per_regions()` in `rnamotifs_core.cpp`

Same principle: for IR events, the three counting regions become:
- Region 1: upstream/downstream flanking (R1 + R4)
- Region 2: first half of retained intron (R2)
- Region 3: second half of retained intron (R3)

### Step 5: IR-aware eCLIP overlap in `rnamotifs_mars_score.cpp`

**File**: `src/cpp/rnamotifs_mars_score.cpp`, function `compute_eclip_overlap()` (lines 294-346)

For IR events, the eCLIP overlap uses the same `in_exon`/`in_intron` semantics — `in_exon` scans into flanking exons (R1, R4), `in_intron` scans into the retained intron (R2, R3):

```cpp
if (e.event_type == "RI") {
    // For IR: swap the exon/intron direction at each splice site
    // R1: upstream exon (in_exon bp before 5'SS)
    // R2: retained intron (in_intron bp after 5'SS)
    // R3: retained intron (in_intron bp before 3'SS)
    // R4: downstream exon (in_exon bp after 3'SS)
    int positions[4] = {
        e.v6 + offset_ei,   // R1: from 5'SS into upstream exon (in_exon extent)
        e.v6 + offset_ie,   // R2: from 5'SS into retained intron (in_intron extent)
        e.v7 + offset_ei,   // R3: from 3'SS into retained intron (in_intron extent)
        e.v7 + offset_ie    // R4: from 3'SS into downstream exon (in_exon extent)
    };
}
```

The fixed window sizes (`in_exon` + `in_intron` positions per region) are appropriate — regulatory elements are concentrated near splice sites. For IR, `in_intron=300` scans 300 bp into the intron from each splice site, which captures the branch point (18-40 nt from 3'SS), polypyrimidine tract, and proximal ISE/ISS elements.

### Step 6: Update orchestrator

**File**: `rnamotifs-mars`

Add `--event-type` argument. When `--event-type RI` is specified:
- Use RI-specific input files from `prepare_mars_exons.py`
- Pass `--event-type RI` to the rnamotifs binary
- Score output and heatmaps are labeled as IR events

### Step 7: Update `read_input_exons()` in C++ to parse event_type

Parse optional 10th semicolon-delimited field. Default to `"SE"` if absent (backward compatible).

## Key design decisions

1. **Backward compatibility**: The 10th field is optional. Existing SE input files work without changes. The C++ code defaults to SE behavior.

2. **Parameter reuse**: `in_exon` controls extent into exonic regions (R1, R4 for IR), `in_intron` controls extent into intronic regions (R2, R3 for IR). No new parameters needed. The same discovery parameter grid applies to both SE and IR events.

3. **dIRank polarity**: Positive dIRank = RBP promotes the event (retention). This is consistent with the SE convention and the ENCODE rMATS sign convention.

4. **Separate runs for SE and RI**: SE and RI events should be analyzed in separate runs rather than mixed, because the region semantics differ. The orchestrator can chain both automatically.

## Visualization: all plots must support IR

All existing plots must work for IR events with appropriate axis labeling. The key difference is that R2/R3 represent intronic (not exonic) sequence, so labels and color semantics must reflect this.

### RNA splicing maps

The per-tetramer positional enrichment bar plots (splicing maps) must relabel the x-axis for IR events:
- R1: "upstream exon" (instead of "upstream intron")
- R2: "5' retained intron" (instead of "exon body")
- R3: "3' retained intron" (instead of "downstream intron")
- R4: "downstream exon" (instead of "distal flank")

**File**: `src/R/mars/config_RNAmars.R`, function `add_splicing_maps_to_final_ht()`

Add an `event_type` parameter. When `event_type == "RI"`, use IR-specific region labels and adjust the region boundary markers (vertical lines at splice sites rather than exon-intron junctions).

### Association heatmap

The heatmap (`generate_heatmap.R`) needs:
- Column annotations (R1/R2/R3 enrichment bars) relabeled for IR
- PEAK binding profile matrix column headers: `R1_enh` → "upstream exon enh", etc.
- Title/subtitle indicating "Intron Retention" instead of "Cassette Exon"

**File**: `src/R/mars/generate_heatmap.R`

Add `event_type` parameter passed via `--event_type` CLI argument. Propagate to `plot_association_heatmap()`.

### PEAK binding profile heatmap (left annotation)

The `plotROI_Heatmap()` function renders the normalized PEAK matrix. For IR, the column labels should read "R1 (exon)", "R2 (5'intron)", "R3 (3'intron)" instead of "R1", "R2", "R3".

**File**: `src/R/mars/config_RNAmars.R`, function `plotROI_Heatmap()`

### Per-tetramer enrichment PDFs

The individual tetramer enrichment plots (generated during Phase 1 by `selection_of_tetramers.R` and `selection_of_tetramers_old_subset.R`) show positional enrichment across the four regions. These need:
- IR-specific region boundary annotations
- Adjusted x-axis labels

**Files**: `src/R/mars/selection_of_tetramers.R`, `src/R/mars/selection_of_tetramers_old_subset.R`

### Implementation approach

Rather than duplicating plot code, add a single `get_region_labels(event_type)` helper function in `config_RNAmars.R`:

```r
get_region_labels <- function(event_type = "SE") {
  if (event_type == "RI") {
    list(R1 = "upstream exon", R2 = "5' retained intron",
         R3 = "3' retained intron", R4 = "downstream exon",
         center = "retained intron", flanks = "flanking exons")
  } else {
    list(R1 = "upstream intron", R2 = "exon body",
         R3 = "downstream intron", R4 = "distal flank",
         center = "cassette exon", flanks = "flanking introns")
  }
}
```

All plotting functions call this to get labels, keeping IR support centralized.

## Testing

1. **PTBP1 validation**: Use PTBP1 HepG2 RI.MATS data (~5,500 events). PTBP1 is a known regulator of intron retention in neuronal genes. Expect motif enrichment of CU-rich tetramers in R3 (near 3' splice site, polypyrimidine tract).

2. **U2AF2 validation**: U2AF2 recognizes the polypyrimidine tract. Expect strong eCLIP signal in R3 for silenced (spliced-out) introns.

3. **Cross-event comparison**: For RBPs that regulate both SE and IR, compare optimal parameters and binding profiles between event types.

4. **Visual validation**: Verify that all plots (splicing maps, heatmaps, PEAK profiles, per-tetramer enrichment) render correctly with IR labels and that region boundaries are drawn at the correct positions.

## Files to modify

| File | Change |
|------|--------|
| `data/prepare_mars_exons.py` | Add `parse_rmats_ri()`, `--event-type` argument |
| `src/cpp/rnamotifs_core.h` | Add `event_type` to Config |
| `src/cpp/rnamotifs_core.cpp` | IR-specific region computation in `process_tetramer_file()` and `count_per_regions()` |
| `src/cpp/rnamotifs_mars_score.cpp` | Add `event_type` to Exon struct, IR-aware `compute_eclip_overlap()` |
| `rnamotifs-mars` | Add `--event-type` CLI argument, pass to all sub-commands |
| `src/R/mars/config_RNAmars.R` | Add `get_region_labels()`, update `plot_association_heatmap()`, `plotROI_Heatmap()`, `add_splicing_maps_to_final_ht()` |
| `src/R/mars/generate_heatmap.R` | Add `--event_type` argument, propagate to plotting functions |
| `src/R/mars/selection_of_tetramers.R` | IR-aware region labels in per-tetramer plots |
| `src/R/mars/selection_of_tetramers_old_subset.R` | IR-aware region labels in group splicing map plots |

## Estimated effort

- Step 1 (prepare_mars_exons.py): ~2 hours
- Steps 2-4 (C++ core changes): ~4 hours
- Step 5 (mars_score eCLIP): ~2 hours
- Steps 6-7 (orchestrator + parsing): ~1 hour
- Visualization (all R plots): ~3 hours
- Testing + validation: ~3 hours
- **Total: ~15 hours**
