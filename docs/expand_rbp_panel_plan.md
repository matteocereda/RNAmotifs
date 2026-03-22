# Expand RBP Panel with Additional eCLIP Binding Data

## Context

RNAmaRs currently uses 15 HepG2 + 13 K562 = 23 unique RBPs (28 total cell-line-specific). The bottleneck is **eCLIP data availability**: rMATS knockdown data exists for 21 RBPs in EACH cell line, but eCLIP binding data is missing for several. ENCODE has eCLIP for ~150 RBPs across HepG2/K562. The goal is to expand the training panel and improve RBP identification coverage.

**Current gap analysis:**

| Cell line | rMATS KD | eCLIP | Missing eCLIP (have KD) |
|-----------|----------|-------|-------------------------|
| HepG2 | 21 RBPs | 15 RBPs | AGGF1, EFTUD2, FXR1, PUS1, RBM15, TARDBP |
| K562 | 21 RBPs | 13 RBPs | HNRNPC, HNRNPK, NCBP2, QKI, RBFOX2, SF3A3, RBM22, UCHL5 |

**Existing infrastructure:**
- `data/download_encode.sh` — queries ENCODE REST API for eCLIP BAMs
- `data/process_eclip.sh` — full iCount pipeline: BAM → crosslinks → peaks → 3nt merge → liftOver
- `data/eclip_qc.py` — quality control reports
- Both scripts have hardcoded RBP lists that need updating

---

## Phase 1: Fill immediate gaps (6 HepG2 + 8 K562 RBPs)

These RBPs already have knockdown data in both cell lines. We just need their eCLIP peaks.

### Step 1: Query ENCODE for available eCLIP experiments

**Script**: `data/query_encode_eclip.py` (new)

Query ENCODE REST API to check which of the 14 missing RBPs have eCLIP data available:

```
HepG2 missing: AGGF1, EFTUD2, FXR1, PUS1, RBM15, TARDBP
K562 missing:  HNRNPC, HNRNPK, NCBP2, QKI, RBFOX2, SF3A3, RBM22, UCHL5
```

API query: `https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP&biosample_ontology.term_name={cell_line}&target.label={RBP}&status=released&format=json`

For each hit, record: experiment accession, RBP, cell line, assembly, file accessions, IDR status, audit flags.

### Step 2: Download and process available eCLIP data

Update `data/download_encode.sh` and `data/process_eclip.sh` RBP lists to include newly found RBPs. Run the existing pipeline:

```bash
bash data/process_eclip.sh --rbp HNRNPC --cell-line K562 --assembly both --output data/eCLIP_processed
# ... repeat for each missing RBP with available data
```

Quality filters to apply:
- IDR ≤ 0.05 (ENCODE standard for reproducibility)
- ≥2 biological replicates
- Audit: no ERROR-level flags
- Peak saturation: keep only if ≥80% peaks retained at 50% read downsampling

### Step 3: Update hardcoded RBP lists

Files to update (4 locations):
- `rnamotifs-mars` → `get_rbps_for_cell_line()`
- `src/R/mars/generate_heatmap.R` → rbps vectors
- `src/cpp/rnamotifs_mars_score.cpp` → rbps vectors in main()
- `data/process_eclip.sh` → HEPG2_RBPS / K562_RBPS

**Better approach**: Refactor to read RBP lists from a single config file (`data/mars_reference/rbp_panels.json`) instead of hardcoding in 4 places.

### Step 4: Generate mars_exons and re-run discovery

```bash
# Already done for all 21 RBPs in both cell lines
# Re-run discovery with expanded RBP panel
./rnamotifs-mars dummy.txt --mode discovery --cell-line HepG2 ...
./rnamotifs-mars dummy.txt --mode discovery --cell-line K562 ...
```

---

## Phase 2: Expand beyond current 21 RBPs with new ENCODE eCLIP

ENCODE has eCLIP for ~103 RBPs in HepG2 and ~120 in K562. To add new RBPs beyond the current 21, we need BOTH eCLIP and knockdown data.

### Step 5: Query ENCODE for RBPs with both eCLIP + shRNA-seq

Query ENCODE for all RBPs that have:
1. eCLIP experiment (released, HepG2 or K562)
2. shRNA-seq knockdown RNA-seq experiment (released, same cell line)
3. Matching control RNA-seq

Cross-reference to find RBPs with complete data (eCLIP + KD + control). The ~70-80 additional RBPs from ENCODE that have both types should be candidates.

### Step 6: Download shRNA-seq and run rMATS

For new RBPs not in the current 21:
1. Download shRNA-seq + control BAM files from ENCODE
2. Run rMATS-turbo for skipped exon quantification
3. Generate mars_exons input files

### Step 7: Download and process new eCLIP

Same pipeline as Step 2, for the additional RBPs.

### Step 8: Run discovery mode for expanded panel

Run full discovery on the expanded panel. With ~50+ RBPs this will take ~30+ hours per cell line. Use crash-safe resumption.

---

## Phase 3: External CLIP-seq databases (optional, medium risk)

For RBPs without ENCODE eCLIP, alternative databases can provide binding site data:

| Database | RBPs | Data types | Risk |
|----------|------|-----------|------|
| POSTAR3 | 351 | PAR-CLIP, iCLIP, eCLIP | Medium — different peak calling |
| CLIPdb | 111 | Multiple CLIP types | Medium — unified processing |
| starBase/ENCORI | 200+ | CLIP-seq aggregated | Low-medium — curated |

**Key consideration**: Non-ENCODE CLIP data uses different peak calling and quality thresholds. Would need:
- Standardized re-processing through our iCount pipeline (if BAMs available)
- Or careful peak format conversion + additional QC
- Benchmark against ENCODE RBPs to calibrate

---

## Refactoring: Dynamic RBP panel

### Step 9: Replace hardcoded RBP lists with config file

Create `data/mars_reference/rbp_panels.json`:
```json
{
  "HepG2": {
    "training_rbps": ["HNRNPC", "HNRNPK", ...],
    "has_eclip": true,
    "has_knockdown": true
  },
  ...
}
```

Update all 4 code locations to read from this file. This makes adding new RBPs a config change rather than code change.

### Step 10: Auto-detect RBPs from available data

Instead of hardcoding, scan `data/eCLIP_processed/{cell_line}/hg19/` for available .bed files and `data/mars_exons/{cell_line}/` for available exon files. The intersection defines the usable training panel. This eliminates manual list maintenance entirely.

---

## Verification

1. After Phase 1: Verify new eCLIP peaks pass QC (`data/eclip_qc.py`)
2. After Step 3: Run application mode with PTBP1 test to ensure expanded panel works
3. After discovery: Verify optimal_parameters.csv has entries for all new RBPs
4. Cross-validate: For RBPs available in both cell lines, check parameter consistency
5. Benchmark: Compare RBP ranking accuracy before/after expansion using known splicing regulators

## Priority

**Phase 1 (immediate, ~1 day)**: Fill gaps for the 14 missing RBPs — highest value per effort since knockdown data already exists.

**Step 9-10 (immediate, ~2 hours)**: Refactor hardcoded lists — prevents ongoing maintenance pain.

**Phase 2 (~1 week)**: Expand to ~50+ RBPs — requires shRNA-seq download + rMATS processing.

**Phase 3 (future)**: External databases — only if specific RBPs are needed that ENCODE doesn't cover.
