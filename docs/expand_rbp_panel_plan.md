# Expand RBP Panel with Additional eCLIP Binding Data

> **Refined June 2026** in light of the with-binding (wb) re-run. The boxed
> **⚠️ Lessons** below are operational facts verified during that work — read them
> before executing any step; several would otherwise silently corrupt the panel.

## Context

RNAmaRs currently uses 15 HepG2 + 13 K562 = 23 unique RBPs (28 cell-line-specific).
The bottleneck is **eCLIP data availability**: rMATS knockdown data exists for 21 RBPs
in EACH cell line, but eCLIP binding data is missing for several. ENCODE has eCLIP for
~150 RBPs across HepG2/K562. The goal is to expand the training panel and improve RBP
identification coverage.

**Current gap analysis:**

| Cell line | rMATS KD | eCLIP | Missing eCLIP (have KD) |
|-----------|----------|-------|-------------------------|
| HepG2 | 21 RBPs | 15 RBPs | AGGF1, EFTUD2, FXR1, PUS1, RBM15, TARDBP |
| K562 | 21 RBPs | 13 RBPs | HNRNPC, HNRNPK, NCBP2, QKI, RBFOX2, SF3A3, RBM22, UCHL5 |

**Existing infrastructure:**
- `data/download_encode.sh` — queries ENCODE REST API for eCLIP BAMs
- `data/process_eclip.sh` — BAM → crosslinks → peaks → liftOver. **Target: a single
  iCount pipeline for both eCLIP and iCLIP** (with SMInput/IDR control); currently uses
  a bedtools fallback on the IP only (see Lesson 2).
- `data/eclip_qc.py` — quality control reports
- `data/build_wb_mars_exons.py` — builds the Methods-faithful (wb) RNAmotifs input
- Hardcoded RBP lists in 3 code locations (see Step 3) still need refactoring

---

## ⚠️ Lessons from the wb re-run (must-read before expanding)

1. **eCLIP peaks lose strand — `process_eclip.sh` Step 7 is buggy.** `final_merge`
   (`bedtools merge -d 0` with **no `-s`**, then `awk '{print $1,$2,$3}'`) outputs
   **BED3 and merges +/− strands together**. Steps 4–6 keep strand; Step 7 drops it.
   Consequence: any **with-binding exon selection** done against these peaks is
   **unstranded → ~2× over-selection** (measured: PTBP1/K562 116 vs the published 56).
   **Fix before processing new eCLIP:** make `final_merge` strand-aware
   (`bedtools merge -s -d 0 -c 5,6 -o sum,distinct`, emit BED6). There is **no stranded
   peak file on disk** to recover from, so new RBPs must be re-processed with the fix.
   *(For the current panel we sidestepped this by reusing the published wb exon lists.)*

2. **Process ALL CLIP (eCLIP + iCLIP) with iCount, and include the input/IDR control
   (REQUIRED).** Going forward both eCLIP and iCLIP are processed through a **single,
   unified iCount pipeline**: BAM → crosslink sites (cDNA 5′ ends) → iCount
   **significant crosslinks / peaks** (FDR-based segmentation). The current
   `process_eclip.sh` instead uses a **bedtools fallback** (`bedtools merge -s -d 3`,
   ≥2 crosslinks) "in the absence of a pre-computed iCount segmentation"
   (`docs/data_preprocessing_methods.md`) and processes the **IP only** — **no
   SMInput normalisation, no IDR** (`merge_replicates` is a peak *union*). Replace the
   fallback with proper iCount calling, and add the control:
   (a) **eCLIP**: input-normalise against the **size-matched input (SMInput)**
   (fold-enrichment + significance) before peak calling — the input BAM is already
   downloaded ("IP + size-matched input") but currently unused;
   (b) **iCLIP**: iCount significant-crosslink FDR (iCLIP has no SMInput; significance
   comes from iCount's per-segment randomisation);
   (c) **both**: **IDR ≤ 0.05** across biological replicates → keep reproducible peaks
   only. Without this, background shared by IP and input inflates the "bound" set and
   the with-binding exon selection.
   *(The current panel reuses the published wb exon lists, so this gap does not affect
   the in-flight wb run, but it must be fixed before processing any new RBP.)*

3. **`dIRank` is categorical, not ΔΨ.** RNAmotifs/mars_score classify exons by
   `|dIRank| ≥ 1` (regulated) / `|dIRank| ≤ 0.1` (control). Writing the real ΔΨ
   (∈ (-1,1)) into that column makes **every regulated exon read as control → empty
   run / all-zero SRR**. New `mars_exons` files **must use `{-1, 0, +1}`**.

4. **`mars_score.cpp` threshold fix is required and applied.** Previously line 130 used
   `> 1.0` (strict), which rejected `dIRank == ±1` → control. Now `>= 1.0` / `<= -1.0`
   (aligned with `core.cpp`'s `>= dIRO`). **Rebuild `build/rnamotifs_mars_score`** in any
   fresh checkout before running discovery on categorical input.

5. **SRR/PEAK reference must be recomputed per panel.** Discovery scoring reads per-RBP
   SRR from `<mars-dir>/Rdata/<cell>_AUC.tsv`; if absent it falls back to **dummy
   SRR = 1 (CS-only, unfaithful)**. Recompute with `compute-peak` (Step 4.5).

6. **Compute is memory-bandwidth-bound, not core-bound.** Measured grid pace ≈
   **107 min/RBP (HepG2) / 122 min/RBP (K562)** at `-c 11`, B=1000. One combo uses
   ~7.4 effective cores; running 2 combos in parallel at `-c 5` was *slower per combo*
   (8.0 vs 4.6 min/combo throughput). **Adding cores or parallel workers does not help.**
   Budget ~110 min/RBP/cell for the grid (see Step 8).

7. **Raw sweep dirs are not auto-pruned and balloon.** rnamotifs-mars prunes only the
   non-optimal *score* dirs; the raw RNAmotifs output `results/<rbp>_<combo>[suffix]/`
   (~160 MB each, ~3.8 GB/RBP) accumulates and is **redundant** (the sims read the
   tiny organized `MaRs_discovery/.../sweep`). Run `data/prune_wb_raw.sh`-style cleanup
   (adapt its glob) during long runs, or you will exhaust disk at 50+ RBPs.

---

## Phase 1: Fill immediate gaps (6 HepG2 + 8 K562 RBPs)

These RBPs already have knockdown data in both cell lines. We need their eCLIP peaks.

### Step 1: Query ENCODE for available eCLIP experiments

**Script**: `data/query_encode_eclip.py` (new). Check which of the 14 missing RBPs have
eCLIP:

```
HepG2 missing: AGGF1, EFTUD2, FXR1, PUS1, RBM15, TARDBP
K562 missing:  HNRNPC, HNRNPK, NCBP2, QKI, RBFOX2, SF3A3, RBM22, UCHL5
```

API: `https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP&biosample_ontology.term_name={cell_line}&target.label={RBP}&status=released&format=json`

Record per hit: experiment accession, RBP, cell line, assembly, file accessions, IDR
status, audit flags.

### Step 2: Download and process available eCLIP data

> **⚠️ Apply Lesson 1 first**: patch `process_eclip.sh` `final_merge` to preserve
> strand (BED6) before running, otherwise the wb exon selection in Step 4 will
> over-select ~2×.

Update the RBP lists in `data/download_encode.sh` / `data/process_eclip.sh`, then:

```bash
bash data/process_eclip.sh --rbp HNRNPC --cell-line K562 --assembly both --output data/eCLIP_processed
# ... repeat for each missing RBP with available data
```

Quality filters (target pipeline, **with the input/IDR control of Lesson 2**):
- Peaks **input-normalised against the SMInput** (fold-enrichment + significance)
- **IDR ≤ 0.05** across the two biological replicates (reproducible peaks only)
- ≥2 biological replicates
- Audit: no ERROR-level flags
- Peak saturation: keep only if ≥80% peaks retained at 50% read downsampling
- QC report via `data/eclip_qc.py`

### Step 3: Update hardcoded RBP lists

Three code locations (note: `process_eclip.sh` lists are data-prep, not panel):
- `rnamotifs-mars` → `get_rbps_for_cell_line()` (candidate list; discovery then
  intersects it with the exon files actually present in `--mars-exons-dir`)
- `src/R/mars/generate_heatmap.R` → `rbps` vectors (scoring/heatmap panel)
- `src/cpp/rnamotifs_mars_score.cpp` → `rbps` vectors in `main()` (scoring panel)

**Better approach**: see [Refactoring](#refactoring-dynamic-rbp-panel) (Steps 9–10).

### Step 4: Generate `mars_exons` (Methods-faithful, with-binding)

> This step replaces the old one-liner. The RNAmotifs input must be the **with-binding
> (wb) subset**, not all regulated exons (Methods §"Selection of RBP-regulated CEs";
> see `reference_exon_selection_methods` and `data/build_wb_mars_exons.py`).

For each new RBP:
1. **Alternatively spliced CEs**: rMATS SE with `|Δμ(Ψ)| > 0.1` & `FDR < 0.1` & UCSC
   hg19 `knownAlt` cassetteExon. **Controls**: `|Δμ(Ψ)| < 0.01` & `FDR > 0.1` &
   RefSeq-but-not-knownAlt.
2. **With-binding filter**: keep regulated CEs with **≥1 (stranded!) iCount peak within
   300 nt into introns / 30 nt into exons** of the splice sites (whole feature if
   intron/exon < 600/60 nt). This is the `02_exons_with_one_peak.R` logic; reuse it
   only after the strand fix (Lesson 1).
3. Require **≥ 50 regulated wb CEs** per RBP (the panel-membership criterion).
4. Write `data/mars_exons/<cell>/<rbp>_input_rnamotifs.txt` with **categorical
   `dIRank` `{-1,0,+1}`** (Lesson 3), `;`-delimited, columns
   `idx;ID;chr;strand;upstreamEE;exonStart;exonEnd;downstreamES;dIRank`.

### Step 4.5: Recompute the SRR/PEAK reference (NEW — required, Lesson 5)

```bash
for cell in HepG2 K562; do
  build/rnamotifs_mars_score --compute-peak \
    --mars-exons-dir data/mars_exons/$cell -e data/eCLIP_processed/$cell/hg19 \
    -c $cell --in-exon 30 --in-intron 300 -b 1000 -o data/mars_reference/Rdata/
  cp data/mars_reference/Rdata/${cell}_SRR.tsv data/mars_reference/Rdata/${cell}_AUC.tsv
done
```
Verify the SRR values are **non-zero** before discovery (zeros ⇒ the dIRank/threshold
issues of Lessons 3–4).

### Step 5: Re-run discovery

```bash
./rnamotifs-mars dummy.txt --mode discovery --cell-line HepG2 \
    --mars-exons-dir data/mars_exons/HepG2 --eclip-dir data/eCLIP_processed/HepG2/hg19 \
    --mars-dir data/mars_reference -n discovery_HepG2 -g hg19 \
    -c 11 -b 1000 --p-empirical 0.01 --p-fisher 0.1 --in-exon 30 --in-intron 300
# ... same for K562. Resumable via discovery_manifest.json.
```
Then `data/recalculate_auroc.py --cell-line <cell> --results-dir results/MaRs_discovery/<cell>_hg19`.

---

## Phase 2: Expand beyond current 21 RBPs with new ENCODE eCLIP

ENCODE has eCLIP for ~103 RBPs (HepG2) / ~120 (K562). New RBPs need BOTH eCLIP and KD.

### Step 6: Query ENCODE for RBPs with both eCLIP + shRNA-seq
RBPs with (1) released eCLIP, (2) released shRNA-seq KD, (3) matching control, same
cell line. The ~70–80 ENCODE RBPs with both types are candidates.

### Step 7: Download shRNA-seq, run rMATS, build wb mars_exons
1. Download shRNA-seq + control BAMs. 2. Run rMATS-turbo (SE). 3. Build `mars_exons`
via **Step 4** (wb filter + categorical dIRank).

### Step 8: Process new eCLIP and run discovery for the expanded panel

eCLIP per Step 2 (with the strand fix). Then discovery + Step 4.5 SRR precompute.

> **⚠️ Revised compute/disk budget (Lessons 6–7).** Measured grid ≈ **110 min/RBP**
> at `-c 11`, B=1000 — and it is **memory-bandwidth-bound**, so more cores/parallel
> workers do **not** speed it up. For a ~50-RBP/cell panel:
> - **Grid**: ~50 × 110 min ≈ **~92 h/cell** (~4 days/cell). The old "~30 h" estimate
>   was ~3× optimistic.
> - **Bayes** (if used): ~(8+`bo-n-iter`) evals/RBP × 110/20 min/eval; budget
>   accordingly (and consider a reduced `--bo-n-iter`, as done in the wb run).
> - **Disk**: raw sweeps ~3.8 GB/RBP unpruned → ~190 GB/cell. Run the
>   `data/prune_wb_raw.sh`-style janitor (it deletes already-converted raw dirs;
>   the sims read the small organized sweep).
>
> Use crash-safe resumption (the manifest skips completed RBPs/combos), and keep
> `-c` **fixed** across resumes (bootstrap RNG is per-thread-seeded → reproducible
> only at constant core count).

---

## Phase 3: External CLIP-seq databases (optional, medium risk)

| Database | RBPs | Data types | Risk |
|----------|------|-----------|------|
| POSTAR3 | 351 | PAR-CLIP, iCLIP, eCLIP | Medium — different peak calling |
| CLIPdb | 111 | Multiple CLIP types | Medium — unified processing |
| starBase/ENCORI | 200+ | CLIP-seq aggregated | Low-medium — curated |

Non-ENCODE CLIP uses different peak calling/thresholds. Would need standardized
re-processing through the (strand-fixed) crosslink pipeline, or careful peak conversion
+ QC, benchmarked against ENCODE RBPs. **Strand must be preserved** for the wb filter.

---

## Refactoring: Dynamic RBP panel

### Step 9: Replace hardcoded RBP lists with a config file
Create `data/mars_reference/rbp_panels.json`:
```json
{ "HepG2": { "training_rbps": ["HNRNPC", "HNRNPK", "..."],
             "has_eclip": true, "has_knockdown": true } }
```
Update the 3 code locations (Step 3) to read it. Adding an RBP becomes a config change.

### Step 10: Auto-detect RBPs from available data
Discovery **already** filters the candidate list to RBPs that have an exon file in
`--mars-exons-dir` (`available_rbps`). Extend this to the **scoring panel** too: define
the usable panel as the intersection of `data/eCLIP_processed/<cell>/hg19/*.bed`,
`data/mars_exons/<cell>/*` and an SRR entry in `<cell>_AUC.tsv` — eliminating the
hardcoded scoring vectors in `generate_heatmap.R` / `rnamotifs_mars_score.cpp`.

---

## Verification

1. **eCLIP**: peaks pass `data/eclip_qc.py`; spot-check a `.bed` is **stranded/BED6**
   (Lesson 1) before exon selection.
2. **mars_exons**: `dIRank` ∈ `{-1,0,1}` only; ≥50 regulated wb CEs/RBP.
3. **SRR**: `<cell>_AUC.tsv` non-zero for every new RBP (Lesson 5).
4. **Binary**: `mars_score.cpp` uses `>= 1.0` and is rebuilt (Lesson 4).
5. **Discovery**: `optimal_parameters` / manifest has entries for all new RBPs; run a
   PTBP1 application-mode sanity check.
6. **Cross-validate**: for RBPs in both cells, check parameter consistency.
7. **Benchmark**: RBP-ranking accuracy before/after expansion on known regulators.

## Priority

- **eCLIP-pipeline prerequisites (do *first*; without them new data is silently
  wrong):** add **input/IDR control** normalisation (Lesson 2) + **strand fix**
  (Lesson 1) to `process_eclip.sh`; recompute the **SRR reference** (Step 4.5); rebuild
  `mars_score` (Lesson 4).
- **Phase 1 (~1 day of work + ~1–2 days compute/cell)**: fill the 14 gaps.
- **Steps 9–10 (~half day)**: refactor hardcoded lists.
- **Phase 2 (multi-week, compute-dominated)**: ~50+ RBPs ⇒ ~4 days grid/cell; plan disk.
- **Phase 3 (future)**: external databases only for RBPs ENCODE lacks.
