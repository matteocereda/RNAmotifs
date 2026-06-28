# Post-simulation TODO (code & docs follow-ups)

Apply **after the with-binding (wb) re-run pipeline is finished** (grid + Bayes +
sims). Do **not** start while the pipeline is running — both items touch
`rnamotifs_mars_score` / shared code, and a rebuild mid-run would disrupt in-flight
scoring. (Manuscript-text edits are tracked separately in
`manuscript_revision_TODO.md`.)

---

## [ ] 0. DIAGNOSTIC — Bayes AUROC = 0.0000 for some wb RBPs (SRSF1, SF3A3, …)

**Observed (2026-06-17, HepG2 wb Bayes):** SRSF1 and SF3A3 returned `AUROC = 0.0000`
across **all** BO init + iters. A flat 0.0000 (not 0.5) means `n_pos = 0` — i.e. **no
enriched MRMs were found** for that RBP under the explored parameters (in the code,
`auroc` is only computed when `true_scores` is non-empty; otherwise `score` stays at its
initialized 0.0). This is an enrichment/data signal, **not** a ranking outcome.

**NOT caused by the 30→15 BO-iteration reduction:** the iteration budget only controls
how many (n, e) points are sampled; it cannot create MRMs the enrichment test does not
call. If all 8 init points + iters yield zero MRMs, more iterations would not change it.
The likely cause is upstream: the **wb (with-binding) exon set is ~10× smaller** than the
all-regulated set (see [[project_wb_vs_allregulated_input]]), plus the enrichment
thresholds (`--p-empirical 0.01`, `--p-fisher 0.1`, ≥5-MRM rule, `--in-intron 300`).

**Checks to run once Bayes completes (decisive → cause):**
1. **Grid vs Bayes on wb:** does the exhaustive **grid** (20 combos) also give SRSF1/SF3A3
   = 0 / SKIP on the wb set? If **yes** → data/thresholds, not BO at all (and not the
   iteration count). If grid found MRMs but BO did not → BO *sampling* is implicated, and
   the lever is **init coverage**, not 15-vs-30 iters.
2. **Whole-exon vs wb:** were these RBPs non-zero in the original all-regulated run? If
   non-zero on whole-exon but zero on wb → confirms the wb-subset-size explanation
   (expected Methods-faithful tradeoff; disclosable in one sentence).
3. **wb exon counts** for the zero RBPs (enh+sil) and the **per-combo n_tet** (number of
   enriched MRMs) — confirm n_tet = 0 throughout.
4. Enumerate the **full list** of wb RBPs with AUROC = 0 in both cells before drafting any
   rebuttal sentence; if it's a systematic small set, note it as a limitation rather than
   a bug.

---

## [x] 1. Make SRR optional — DONE 2026-06-26 (`--score-mode full|cs-only`)

_Implemented: C++ `rnamotifs_mars_score --score-mode` (flatten_scores uses CS alone when cs-only); `rnamotifs-mars --score-mode` threaded to both grid+bayes scoring invocations; README + tutorial updated. Verified: PTBP1 full=0.940, cs-only=0.928 (== ablation). Default `full` unchanged._

### (original spec below)
## 1. Make SRR optional in the combined association score

**Why (reviewer-driven):** the Comment 1.4 ablation showed CS is the dominant term and
SRR is an *orthogonal robustness safeguard* (and Comment 2.12 motivates SRR as
inverse-variance-style weighting). Follow-up: expose SRR as **optional** so users can
score on **CS alone** (AS = CS) when they prefer it, with the full AS = CS × SRR as the
default.

**To do:**
- Add a toggle (e.g. `--score-mode {full,cs-only}` or `--no-srr` / `--srr-weight`) in
  `rnamotifs-mars` and the `rnamotifs_mars_score` binary so AS = CS × SRR collapses to
  AS = CS (SRR ≡ 1) when disabled. (The ablation already computes the CS-only condition
  as a diagnostic — formalise it as a first-class scoring option.)
- Thread it through: `src/cpp/rnamotifs_mars_score.cpp` (the SCORE1×SCORE2 combination /
  `flatten_scores`), the `rnamotifs-mars` CLI, and any R/config path that consumes AS.
- Default = full (CS × SRR), unchanged; document the option in `docs/tutorial/` +
  `docs/tutorial/rnamotifs-mars.md` + README.

## [x] 2. Align code & docs nomenclature with the paper — DONE 2026-06-26

_Scope decision: keep internal on-disk keys (`hw_X_ew_Y` manifest keys, `SCORE1_`/`SCORE2_`
diagnostic filenames) for stability — 12 scripts + all canonical results depend on them.
Applied the safe, user-facing alignment instead:_
- _`rnamotifs` CLI: added `--n` (paper notation, `n = 2·hw`, even-number validated; `--n 30` ⇔ `-w 15`); `-w`/`-e` help text annotated with paper `n`/`e`._
- _README.md: added a "Notation: paper ↔ code" mapping table (n↔hw, e↔ew, SRR↔SCORE1, CS↔SCORE2, AS↔combined, MRM↔tetramer) after the AS description._
- _Tutorials: `parameters.md` (`--n`/`n = 2·hw` on `-w`, `e` on `-e`, cheat-sheet), `rnamotifs-mars.md` (SRR/CS/MRM phase labels, AS = CS × SRR, grid notes), tutorial `README.md` (MRM introduced at the tetramer definition)._
- _Verified: `rnamotifs` parses; `--n 31` rejected (must be even), `--n 30` → hw=15._

### (original spec below)
## 2. Align code & docs nomenclature with the paper

**Why:** the code and the paper use different names for the same quantities, which is a
readability/consistency liability (and was flagged in spirit by the reviewers, e.g.
Comment 1.1 on abbreviations). Harmonise the code, CLI flags, log messages, and docs to
the **paper's nomenclature**.

**Key mappings (code → paper):**
| Code | Paper | Note |
|------|-------|------|
| `hw` / `--half-window` | **n** (clustering window) | paper reports **n = 2·hw**; decide a single convention (expose `n` and/or document `n = 2·hw` clearly) |
| `ew` / `--enrichment-window` | **e** (enrichment window) | |
| `SCORE1` | **SRR** (Signal Recovery Rate) | |
| `SCORE2` | **CS** (Cosine Similarity) | |
| combined score | **AS** (Association Score) = CS × SRR | |
| tetramer / motif cluster | **MRM** (multivalent RNA motif) | |
| `dIRank` | regulation class | keep code value, align prose |
| CES / BS | CES / BS | already aligned |

**To do:** rename or alias in `src/cpp/*`, `rnamotifs` / `rnamotifs-mars` CLI help and
output, `src/R/mars/*`, then sweep `README.md` and `docs/tutorial/*` to match. Be
careful with the **n = 2·hw** convention so the CLI, docs, and paper agree on whether a
user passes the half-window or the full window.
