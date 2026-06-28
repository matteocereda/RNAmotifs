# RNAmotifs-MaRs

`rnamotifs-mars` couples RNAmotifs with **eCLIP** RBP-binding data to score how well a
discovered motif set matches an RBP's binding — the **MRM–RBP association score**. It
extends the [RNAMaRs](https://github.com/ceredamatteo-lab/theRNAmars) framework by
running RNAmotifs motif discovery directly inside the scoring loop.

## How it works (three phases)

| Phase | Name | What it does |
|-------|------|--------------|
| **1** | Multi-parameter sweep | Runs RNAmotifs across `hw`/`ew` combinations (paper notation `n = 2·hw` / `e`) to find enriched MRMs (multivalent RNA motifs; `tetramer` in code) under different clustering stringencies. |
| **2** | Signal Recovery Rate — **SRR** (`SCORE1` in code) | eCLIP-downsampling robustness of the RBP's binding profile. |
| **3** | Cosine Similarity — **CS** (`SCORE2` in code) | Cosine similarity between the motif-enrichment profile and the eCLIP binding profile. |

The **Association Score (AS) = CS × SRR** (`SCORE2 × SCORE1` in code), aggregated into per-RBP heatmaps.
The SRR weight is optional: pass `--score-mode cs-only` to score on the cosine similarity alone
(AS = CS), or keep the default `--score-mode full` for AS = CS × SRR.

## Two modes

- **Discovery** (`--mode discovery`) — sweeps `hw`/`ew` per RBP (grid or Bayesian
  optimisation, `--optim grid|bayes`) and picks the parameters that best identify the
  true RBP, writing a crash-safe, **resumable** manifest.
- **Application** (`--mode application`) — scores a new splicing dataset against a
  trained reference panel.

## Discovery example

```bash
./rnamotifs-mars dummy.txt --mode discovery \
    --cell-line HepG2 --mars-exons-dir data/mars_exons/HepG2 \
    --eclip-dir data/eCLIP_processed/HepG2/hg19 \
    --mars-dir data/mars_reference -n discovery_HepG2 -g hg19 \
    -c 8 -b 1000 --p-empirical 0.01
```

Intermediate results: `results/MaRs_discovery/<cell>_<genome>/<rbp>/` with `sweep/`
(RNAmotifs runs) and `scores/` (association scores + `diagnostics/`).

## Key extra options

| Flag | Default | Meaning |
|------|---------|---------|
| `--mode` | `application` | `discovery` (train params) or `application` (score new data). |
| `--cell-line` | — | `HepG2` or `K562` (eCLIP panel). |
| `--eclip-dir` | — | Path to eCLIP peak files. |
| `--mars-dir` | — | RNAmars data dir (`Tables/`, `Rdata/`). |
| `--mars-exons-dir` | — | Per-RBP RNAmotifs input files (discovery). |
| `--param-grid-n` | `10 30 50 70` | Clustering-window grid `n` (paper notation, `n = 2·hw`; each value even). Back-compat alias `--param-grid-hw` (= `n/2`) still accepted. |
| `--param-grid-e` | `30 50 100 200 300` | Enrichment-window grid `e` (alias `--param-grid-ew`). |
| `--optim` | `grid` | `grid` or `bayes` (Bayesian optimisation). |
| `--make-heatmaps` | off | After discovery, run `generate_heatmap.R` for every trained RBP → the Figure-4A final heatmaps (AS heatmap + per-RBP RNA splicing maps), using the grid optima. Reuses existing sweeps (no re-search); writes `…/<cell>_<genome>/heatmaps/<cell>_<RBP>_<enh\|sil>.pdf`. Resumable. |
| `--heatmap-top-mrms` | `0` (all) | Crop the Phase-3 heatmaps to the top-N enriched MRMs (columns). `5` reproduces the Figure-4A top-5 crop. |
| `--heatmap-top-rbps` | `0` (all) | Crop the Phase-3 heatmaps to the top-N RBPs (rows). |
| `--bo-n-init` / `--bo-n-iter` | 8 / 30 | Bayesian-opt initial points / iterations. |
| `--bo-seed` | 42 | Bayesian-opt random seed. |
| `--in-intron` | `300` | One or more intron-extent values (grid-searchable). |

All shared RNAmotifs parameters (`-b`, `--p-fisher`, `--p-empirical`, `--in-exon`,
`-k`, …) work as documented in [Parameters](parameters.md). Discovery mode is
**resumable**: re-running the same command skips completed RBPs/combinations via the
`discovery_manifest.json`.

See the project README's **RNAmotifs-MaRs** section for the full end-to-end workflow.


---
[← Back to tutorial index](README.md)
