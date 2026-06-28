# Reproducibility — script deposit map

Every script needed to reproduce the manuscript, grouped by layer and in execution
order. Paths are relative to the repo root. The **canonical pipeline** for the BMC/MSB
revision is the with-binding (`_wb`) track: build engine → preprocess → `run_wb_all.sh`
→ `run_wb_sims.sh` → figures.

Notation: discovery is driven in **paper notation** (`--param-grid-n`, `n = 2·hw`;
`--param-grid-e`); on-disk result keys remain `hw`-based (`hw = n/2`). See README
"Notation: paper ↔ code".

---

## Layer 0 — Build system & compute engine (C++)
Compiled to `build/` via CMake. These are the algorithmic core.

| File | Role |
|------|------|
| `CMakeLists.txt` | Build definition (8 executables). |
| `src/cpp/rnamotifs_core.cpp` | Shared core (I/O, regions, genome). |
| `src/cpp/tetramer.cpp` | k-mer (MRM) enumeration. |
| `src/cpp/counting.cpp` | Per-position cluster counting. |
| `src/cpp/rnamotifs_search.cpp` | Enrichment search. |
| `src/cpp/bootstrap_fdr.cpp` | Empirical-FDR null (per-iteration RNG seed 30580; core-count-independent). |
| `src/cpp/rnamotifs_selection.cpp` | Motif selection / Fisher + empirical filter. |
| `src/cpp/rnamotifs_mars_score.cpp` | MaRs association score (CS×SRR; `--score-mode full\|cs-only`; `--compute-peak`). |
| `src/cpp/rnamotifs_conservation.cpp` | 🚧 PhyloP conservation (optional). |
| `src/cpp/rnamotifs_structure.cpp` | 🚧 ViennaRNA structure (optional). |

## Layer 1 — CLI entry points (Python)
| File | Role |
|------|------|
| `rnamotifs` | Single-RBP motif discovery + splicing map (`--n` paper notation). |
| `rnamotifs-mars` | MaRs orchestrator: discovery (grid/bayes) + application (`--param-grid-n/e`, `--bo-n-range/e-range`, `--score-mode`). |
| `rnamotifs-extract` | Extract enriched-motif coordinates (TSV/BED) from a results folder. |

## Layer 2 — R support library (called by the engine: plotting, selection, scoring)
| File | Role |
|------|------|
| `src/R/config.R`, `src/R/selection.R` | Core selection / config. |
| `src/R/cluster_profile.R`, `src/R/ir_splicing_map.R` | Splicing-map profiles. |
| `src/R/conservation_profile.R`, `src/R/structure_profile.R` | 🚧 optional profiles. |
| `src/R/mars/config_RNAmotifs.R`, `src/R/mars/conf/...`, `src/R/mars/config_RNAmars.R` | MaRs config. |
| `src/R/mars/compute_association_scores.R` | CS/SRR/AS computation. |
| `src/R/mars/selection_of_tetramers.R` (+ `_old_subset`) | MRM selection. |
| `src/R/mars/generate_heatmap.R`, `src/R/mars/sign_reg_plot.R`, `src/R/mars/figure_parameter_optimization.R` | MaRs plots. |

---

## Layer 3 — Data acquisition & preprocessing
Run once to build inputs. Outputs feed Layer 4.

| File | Role |
|------|------|
| `genomes/hg19.download.sh`, `hg38.download.sh`, `mm9.download.sh`, `mm10.download.sh`, `mm39.download.sh` | Reference genomes. |
| `genomes/download_phylop.sh` | 🚧 PhyloP tracks (conservation). |
| `data/download_encode.sh` | ENCODE eCLIP + shRNA-seq. |
| `data/download_table_s4.sh` | ENCODE DESeq2 + rMATS (Table S4). |
| `data/process_eclip.sh` | eCLIP → iCount/IDR peaks (`data/eCLIP_processed`). |
| `data/prepare_mars_exons.py` | ENCODE → alt/constitutive exon lists. |
| `data/build_wb_mars_exons.py` | Published **with-binding** selection → RNAmotifs input (`data/mars_exons_wb`). |
| `data/eclip_qc.py` | eCLIP QC report. |

## Layer 4 — Discovery / simulation pipeline (canonical `_wb` run)
The reproducible heart of the revision. **Resource cap: 11 cores / ≤58 GB; seed 30580.**

| File | Role / order |
|------|------|
| `data/run_wb_all.sh` | **Master**: compute-peak (SRR/PEAK) → GRID → BAYES → sims. |
| `data/run_wb_grid.sh` | Grid discovery only (both cell lines). |
| `data/run_wb_sims.sh` | Stage 3: all rebuttal simulations (calls Layer 5). |
| `data/run_wb_recompute.sh` | Deterministic recompute of rebuttal analyses on the `_wb` grid. |
| `data/recalculate_auroc.py` | Panel-filtered pooled true-vs-other AUROC from a discovery folder. |
| `data/build_wb_recalc_table.py` | Full-grid AUROC table from the `_wb` manifest. |
| `data/merge_discovery_results.py` | Pick best `in_intron` per RBP across runs. |

## Layer 5 — Rebuttal / revision analyses (deterministic; invoked by `run_wb_sims.sh`)
| File | Reviewer comment |
|------|------|
| `data/ablation_analysis.py` | 1.4 — CS×SRR ablation (CS-only vs full). |
| `data/loro_cross_validation.py` + `data/rescore_all_combos.py` | 1.6 — Leave-One-RBP-Out CV. |
| `data/cross_cell_validation.py` | 1.5 / 2.14 — cross-cell-line transfer. |
| `data/rescore_bayes_idr.py` | 2.12 — grid vs Bayesian optimisation. |
| `data/collect_rebuttal_numbers.py` | Consolidates every number quoted in PbP.md / Supp §9.3. |

## Layer 6 — Benchmarks (Methods / Supplementary)
| File | Role |
|------|------|
| `bench/run_benchmark.py`, `bench/run_bo_benchmark.sh`, `bench/postprocess_benchmark.py` | BO vs grid benchmark. |
| `bench/run_kmer_bench.py`/`.sh`, `bench/run_kmer_bench_wb.py`, `bench/analyze_kmer_bench.py`, `bench/report_kmer_bench.py` | k-mer-size sweep. |
| `bench/run_seed_test.sh`, `bench/test_option1.py`, `bench/check_corecount.py` | Bootstrap determinism / core-count independence. |
| `bench/chain_crosscell_bayes.sh` | Sequencing helper (resource-safe chaining). |

## Layer 7 — Figure generation (R)
| Directory | Figures |
|------|------|
| `Paper/Figure_2/*.R` … `Paper/Figure_5/*.R` | Original main figures 2–5. |
| `Paper/revision/revised_figure{2,3,4}_*.R` | **Revised** figures for the rebuttal (canonical for this revision). |
| `Paper/ENCODE_validation/Scripts/*.R` | ENCODE validation panels. |
| `Paper/HNRNPK_silencing_72h/Scripts/*.R` | HNRNPK case study. |
| `Paper/conf/*.R` | Figure-pipeline configs. |

---

## Separate track — RBP-panel expansion (GEO; NOT part of this paper's deposit)
These support the future panel-expansion effort (`docs/expand_rbp_panel_plan.md`), not the
current manuscript's results. **Excluded from the reproducibility deposit.** Deposit only
if/when the expansion is in scope.

**GEO survey:** `data/geo_clip_survey.py`, `data/geo_clip_knockdown_survey.py`,
`data/geo_clip_splicing_hif.py`, `data/geo_clip_verify_target.py`,
`data/geo_clip_verify_rnaseq_design.py`, `data/geo_download_rnaseq.py`

**GEO RNA-seq processing (STAR/rMATS):** `data/setup_star_references.sh`,
`data/run_star_alignment.sh`, `data/run_all_star_alignments.sh`, `data/run_rmats_geo.sh`

## Not analysis — exclude from the reproducibility deposit
| File | Why |
|------|------|
| `data/apply_docx_revision_edits.py` | Edits the manuscript `.docx`; not part of the analysis. |
| `data/prune_wb_raw.sh`, `data/prune_wb_raw` housekeeping | Disk cleanup utilities. |
