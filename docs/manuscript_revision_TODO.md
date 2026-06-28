# Manuscript .docx revision checklist

Edits to apply to `Paper/Grieco_Becchi_Boscagli_et_al__MSB_20260217.docx` (and the Supplementary
Note .docx) **once `Paper/PbP.md` is definitive** — i.e. when porting the agreed
rebuttal changes into the manuscript. Each item must keep the manuscript text
consistent with the corresponding PbP reply (AI-assisted reviewers cross-check both).

---

## [x] Comment 2.1 + naming — APPLIED to docx 2026-06-19 (red font)

**Applied:** (a) "cognate" → "candidate" throughout incl. title (8 spans); (b) framework rename
"RNAMaRs" → "RNAmotifs v2.0" (65 spans; GitHub URL `theRNAMaRs` preserved), with the new **MaRs
module** named at the Abstract (¶15) and Overview (¶24); title now "RNAmotifs v2.0: … and their
candidate Regulators of Splicing". Backup `*.bak_20260619_085317_pre_rename.docx`.

**PbP.md naming harmonization — DONE 2026-06-20:** swept RNAMaRs → RNAmotifs v2.0 (0 residual),
rewrote Comment 2.2 to the "RNAmotifs v2.0 + MaRs module" framing, introduced the rename in the
cover letter and Comment 1.1 (so Reviewer 1 sees it explained), title fixed to "their candidate".
Also moved intron-retention / RNA-structure / conservation out of the 2.2 contributions list into a
work-in-progress note (GitHub), consistent with Discussion ¶61.

---

## [x] PbP NUMERIC REFRESH — DONE (all 1.4/1.5/1.6/2.12/2.14 refreshed on _wb) on the with-binding (_wb) results — Comments 1.4, 1.5, 1.6 (+ cross-refs 2.12, 2.14)

**Why:** every quantitative result in these replies was computed on the OLD **whole-exon
(all-regulated)** run. The Methods-faithful pipeline was re-run on the **with-binding (_wb)**
exon set after the suffix/SRR bug fix (2026-06-17), giving materially different AUROCs (e.g.
HepG2 wb mean 0.873 vs whole-exon 0.803; PTBP1 optimal moved hw_15_ew_30→hw_35_ew_50). All
numbers below MUST be regenerated on `_wb` before the PbP is final, and the PbP and manuscript
kept consistent. Inline `<!-- NUMBERS PENDING wb RECOMPUTE -->` markers flag each reply in PbP.md.

**Dependency / ordering:** needs the wb GRID complete for both cells (HepG2 done 2026-06-19;
K562 in progress) AND the per-combo SCORE1/SCORE2 diagnostics retained — the raw-sweep pruner is
currently OFF, so do the recompute BEFORE any cleanup re-prunes `results/MaRs_discovery/*_wb`.

- **[x] Comment 1.4 (ablation) — DONE 2026-06-20.** Fixed `data/ablation_analysis.py` (AUROC was a
  per-tetramer match; now POOLED true-vs-other, so Full-AS reproduces the Figure 4 grid AUROC
  exactly — verified all 28 OK, 0 blanks). **Finding INVERTED the old conclusion:** SRR weighting is
  EMPIRICALLY NEUTRAL on _wb (means Full 0.844 vs CS-only 0.855; helps 5/hurts 6/tie 17; binomial
  p=0.73). Reply rewritten (decision: keep SRR as default-on but OPTIONAL `--score-mode`, framed as
  no-cost robustness prior, neutral here, no over-claim on noisy data). Outputs: results/wb/ablation_{HepG2,K562}.tsv.
- **[x] Comment 2.12 SRR-motivation paragraph — DONE 2026-06-20** (refreshed to the neutral-safeguard
  framing, p=0.73; matches 1.4).
- **[x] Comment 2.12 grid-vs-Bayes BENCHMARK — DONE 2026-06-22.** Fixed `rescore_bayes_idr.py`
  dir-suffix bug (`<cell>_hg19_wb_bayes` → `..._bayes_bayes`, DOUBLE _bayes). wb result INVERTED the
  whole-exon finding: grid is now marginally BETTER than BO (mean grid 0.838 vs BO 0.787, grid≥BO
  21/27), ρ=0.885 (n=27; RBM15/K562 no valid BO score). PbP 2.12 rewritten (no over-claim — noted
  BO's reduced 8+15 budget; "Option A"), benchmark para fixed to 23 combos (was wrongly 38=8+30).
  Cross-refs ρ updated at ~L71 (1.3) and 2.14. Output: results/wb/bayes_vs_grid_idr.csv.
- **[x] Comment 1.5 / 2.14 (cross-cell) — DONE 2026-06-22.** Fixed `cross_cell_validation.py` AUROC
  (same per-tetramer→POOLED fix as ablation; also removed per-direction <5 skip). wb numbers are
  STRONGER than whole-exon: matched 0.885 vs cross 0.868 (Δ −0.017), min 0.69 (was 0.52), ρ 0.67
  (was 0.43), 11/14 within 0.10. PbP 1.5 table + 2.14 updated; pending-marker removed.
- **[x] Comment 1.3 (LORO) — FIXED + UPDATED 2026-06-22.** Bug: `loro_cross_validation.py` used
  `exclude_zeros=True` (rewards sparse combos → n_training=1 → all-zero OOS) AND sims ran the broken
  `--use-improved` path → every OOS AUROC was 0. Fix: exclude_zeros=False (default) + use manifest
  grid (dropped --use-improved in run_wb_sims.sh). New: ρ 0.76 HepG2 / 0.70 K562 (was 0.91/0.59),
  mean LORO 0.732/0.662, UCHL5 & PUS1 → 0 OOS (honest weak-signal losses). PbP 1.3 + 2.14 cross-ref
  updated; Fig 4E regenerated (revised_figure4_E.R repointed to results/wb/_misc/loro_cross_validation.tsv).
- **[x] Comment 1.6 (k=4 vs 5/6) — DONE 2026-06-22.** Ran `bench/run_kmer_bench_wb.py` on the _wb
  exon sets at the wb-optimal (hw,ew) for HNRNPC/QKI/PTBP1/U2AF2, **all k at 11 cores** (Option B) so
  k=4 reproduces the grid/Figure 4 exactly (verified: HNRNPC 1.00, QKI 1.00, PTBP1 0.94, U2AF2 1.00).
  CAUGHT a stale-cache bug first (the k=4 no-suffix dirs reused 8-core results → spurious "k5 improves";
  cleared + re-ran clean). Result: k=4 optimal, no improvement at k>4 (k5: HNRNPC 0.98, QKI 0.67,
  PTBP1 0.91, U2AF2 0.80; k6: HNRNPC 0.99, PTBP1 0.76); mean Δ −0.14 (k5), −0.10 (k6). PbP 1.6 +
  Supp §9.3.5 updated, markers removed. Output: bench/kmer_bench_results_wb.json.
  NOTE on reproducibility: bootstrap is core-count dependent (per-thread RNG); Option-1 permanent fix
  (per-iteration seeding in bootstrap_fdr.cpp) is DEFERRED post-revision (~2h code, ~6d to re-adopt).

---

## [x] Comment 1.2 — Sens/Spec equations (Figure 3 + Methods "ROC-guided selection of RNAmotifs parameters", ¶97) — APPLIED 2026-06-17

**RESOLUTION (done):** switched to the **rank form** that the code actually implements
(`rnamotifs-mars:842–853` grid, `:1237–1265` Bayes) — there is **no L threshold anywhere in the
code** (verified: the only `0.01` hits are unrelated p-value quantiles). Docx ¶97 rewritten in red:
the L-threshold/binary-classification/Sens/Spec narrative was removed and replaced with "ground
truth = knocked-down RBP; positives = target-RBP MRM scores, negatives = other-RBP scores (pooled);
AUROC = P(AS_target > AS_other) = (wins + ½·ties)/(N₊·N₋); per enh/sil, larger retained". Figure 3A
SVG (`Paper/Figure_03_panelA_revised.svg`): Sens/Spec box replaced by the rank-form block, ROC curve
kept as visualization, header typo n∈{…,75}→{…,70} fixed. PbP Comment 1.2 reply rewritten to match.
Backups: `*.bak_20260617_111250_pre_rankform.docx`. NOTE: the figure SVG still needs to be rendered
and dropped into the actual Figure 3 composite (panel A only changed); the docx text is final.

---

## [x] Comment 2.4 — 15-nt smoothing — APPLIED to docx 2026-06-24 (Methods p.21, red)
_(was:)_ ## Comment 2.4 — 15-nt smoothing window vs *n*/*e* (Methods ~line 474 / ¶75)

**Why:** the PbP Comment 2.4 reply states that the 15-nt scrolling-window smoothing is
applied **only to the eCLIP CES profile** (for the **BS** metric and the **Figure 2D–E**
per-position visualisation) and is **independent of *n* and *e***, and that the
association score **AS = CS × SRR is computed from the *raw, unsmoothed* −2·log(Fisher)
profiles** (SRR on the raw profile; CS as a cosine over unsmoothed profiles). This is
code-verified (`src/cpp/rnamotifs_mars_score.cpp`: `smooth15` used only in the
BS/Fig-2D-E path; line 794 "SRR cosine uses the raw … profile (NOT smoothed)"; SCORE2
cosine never smoothed).

**To do in the docx:**
1. Check the Methods passage that introduces the 15-nt smoothing (line ~474 / ¶75). If
   it states or implies the smoothing feeds **SRR** or **CS**, correct it — the
   smoothing applies only to the **BS** metric and the **Figure 2D–E** profile.
2. Add the explicit sentence that the 15-nt window is fixed, applies only to the eCLIP
   CES profile, and is **independent of *n* and *e*** (which act on the sequence-based
   enrichment step).
3. Confirm AS = CS × SRR is described as using the raw/unsmoothed profiles, so the
   manuscript and the PbP Comment 2.4 reply do not contradict each other.
