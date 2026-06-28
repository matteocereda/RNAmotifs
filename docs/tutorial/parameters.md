# Parameters

Full reference for every `rnamotifs` flag, grouped by role. See
[How the algorithm works](README.md#how-the-algorithm-works-in-five-steps) for the
concepts behind them.

## Required

| Flag | Meaning |
|------|---------|
| `input_file` | Positional: the splicing-change file ([Input Format](input-format.md)), or an rMATS file with `--from-rmats`. |
| `-n, --name` | Run name; becomes the results sub-folder and plot title. |
| `-g, --genome` | Reference assembly: `hg19`, `hg38`, `mm9`, `mm10`, `mm39`. The matching genome must be downloaded (`genomes/`). |

## Clustering — *what counts as a motif cluster*

| Flag | Default | What it controls | When to change |
|------|---------|------------------|----------------|
| `-w, --half-window` | 15 | Half-width `hw` (bp) for clustering motif copies; copies within `±hw` join one cluster. In paper notation this is the clustering window **n = 2·hw** (`--n` is an alias: `--n 30` ⇔ `-w 15`). **Larger → looser, longer-range multivalency.** | Increase to capture dispersed motif arrays; decrease for tight local clusters. The MaRs sweep tries `hw ∈ {5, 15, 25, 35}` (n ∈ {10, 30, 50, 70}). |
| `-m, --min-height` | 4 | Minimum number of motif copies stacked within a window to call a cluster. **Higher → stricter.** | Lower (2–3) for sparse motifs / small exon sets; raise for very abundant motifs. |
| `-p, --pth` | 0.5 | Percentage threshold used when binarising cluster occupancy along the map. | Rarely changed; leave at default unless reproducing a specific protocol. |

## Region geometry — *where you look*

| Flag | Default | What it controls |
|------|---------|------------------|
| `-e, --enrichment-window` | 30 | Enrichment window `ew` (bp; paper notation **e**) over which positional enrichment is aggregated. **Larger → smoother, broader signal.** The MaRs sweep tries `{30, 50, 100, 200, 300}`. |
| `--in-exon` | 30 | How far (bp) into the **exon** the analysed region extends from each splice site. |
| `--in-intron` | 300 | How far (bp) into the **intron** the analysed region extends from each splice site. |
| `--event-type` | `SE` | `SE` (skipped exon) or 🚧 `RI` (intron retention — *work in progress*, see below). |

> `hw` (paper **n = 2·hw**) and `ew` (paper **e**) are the two parameters most worth
> tuning — they trade sensitivity vs specificity. There is no universal optimum; the
> **[MaRs discovery mode](rnamotifs-mars.md)** grid- or Bayes-searches them per RBP.

## Statistics — *how strict the filter is*

| Flag | Default | What it controls |
|------|---------|------------------|
| `-b, --bootstraps` | 10000 | Bootstrap iterations for the empirical-FDR null. **1000 is the practical setting** (fast; floor ≈ 10⁻³). 10000 is ~10× slower for a lower floor that rarely changes the called set. |
| `--p-fisher` | 0.1 | Fisher's-exact p-value threshold for regulated-vs-control enrichment. |
| `--p-empirical` | 0.00005 | Empirical p-value threshold from the bootstrap null. *Only reachable with a large `-b`*; with `-b 1000` use a looser value such as `0.01`. |
| `-k, --kmer-size` | 4 | Motif length. `4` (tetramers) is the validated default. `5`/`6` enumerate exponentially more, sparser motifs and generally **reduce** discriminability at higher cost. |

## Performance & plotting

| Flag | Default | Meaning |
|------|---------|---------|
| `-c, --cores` | 1 | CPU threads (k-mer scan + bootstrap are parallelised). **The workload is largely memory-bandwidth-bound**, so throughput plateaus well below very high core counts. Bootstrap reproducibility holds only at a **fixed** `--cores` (per-thread RNG seeding). |
| `--top-n` | all | Plot only the top-N ranked tetramers (cosmetic; does not change the called set). |

## 🚧 Optional analyses — work in progress

These run but are **not yet validated** in this release; treat output as experimental.

| Flag | Default | Meaning |
|------|---------|---------|
| 🚧 `--event-type RI` | — | Intron-retention mode (different region geometry). |
| 🚧 `--structure` | off | RNA secondary-structure accessibility profiles for enriched motifs (**requires ViennaRNA**). |
| 🚧 `--structure-window` | 31 | Folding window size (bp) for structure profiling. |
| 🚧 `--conservation` | off | PhyloP conservation profiles (**requires** pre-downloaded `.phylop.bin`, see `genomes/download_phylop.sh`). |

## rMATS input (with `--from-rmats`)

See [Input Format → Feeding rMATS directly](input-format.md#feeding-rmats-directly-recommended)
for `--rmats-incl`, `--rmats-fdr`, `--rmats-constit`, `--rmats-max-constit`.

---

## Tuning cheat-sheet

| Symptom | Try |
|---------|-----|
| Nothing called (`no_enriched.txt`) | Lower `--min-height` (2–3); loosen `--p-fisher`/`--p-empirical`; widen `-e`; ensure ≥ ~50 regulated exons; **verify `dIRank` is `{-1,0,1}`, not raw ΔΨ**. |
| Too many noisy motifs | Raise `--min-height`; tighten `--p-fisher`/`--p-empirical`; increase `-b`. |
| Clusters too local / fragmented | Increase `-w` and/or `-e`. |
| Want longer footprints | Increase `-w`; only try `-k 5` with a strong prior (`-k 4` usually wins). |
| Runs too slow | Use `-b 1000`; raise `-c` (diminishing returns — memory-bandwidth bound). |
| Don't know best `hw`/`ew` (paper `n`/`e`) | Don't guess — use [`rnamotifs-mars --mode discovery`](rnamotifs-mars.md) to grid/Bayes-search. |


---
[← Back to tutorial index](README.md)
