# RNAmotifs-MaRs

`rnamotifs-mars` couples RNAmotifs with **eCLIP** RBP-binding data to score how well a
discovered motif set matches an RBP's binding — the **MRM–RBP association score**. It
extends the [RNAMaRs](https://github.com/ceredamatteo-lab/theRNAmars) framework by
running RNAmotifs motif discovery directly inside the scoring loop.

## How it works (three phases)

| Phase | Name | What it does |
|-------|------|--------------|
| **1** | Multi-parameter sweep | Runs RNAmotifs across `hw`/`ew` combinations to find enriched tetramers under different clustering stringencies. |
| **2** | Signal Recovery Rate (SCORE1) | eCLIP-downsampling robustness of the RBP's binding profile. |
| **3** | Profile similarity (SCORE2) | Cosine similarity between the motif-enrichment profile and the eCLIP binding profile. |

The **association score (AS) = SCORE2 × SCORE1**, aggregated into per-RBP heatmaps.

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
| `--param-grid-hw` | `5 15 25 35` | Half-window grid. |
| `--param-grid-ew` | `30 50 100 200 300` | Enrichment-window grid. |
| `--optim` | `grid` | `grid` or `bayes` (Bayesian optimisation). |
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
