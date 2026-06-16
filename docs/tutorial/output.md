# Output & Interpretation

Results land in
`results/<name>/<timestamp>_<genome>_w<hw>_..._ew<ew>_..._b<boot>_.../`.

## Key files

| File | What it is |
|------|------------|
| `*-tets-*.pdf` | **The RNA splicing map** — per-position motif-cluster frequency, enhanced (red) vs silenced (blue) vs control, for the enriched tetramers. The headline result. |
| `enriched_tetramers.txt` | Plain list of tetramers that passed both significance filters. |
| `*-tets-*.csv` | Per-tetramer enrichment table (region-wise Fisher and empirical p-values). |
| `All_group_{enh,sil,both}-*.csv` | Per-position occupancy curves for **all** candidate tetramers, split by group (feed the splicing map and the MaRs cosine-similarity step). |
| `Enh_*`, `Sil_*`, `ES_*`, `e-*`, `s-*` | Group-specific enrichment matrices (enhanced-only, silenced-only, both). |
| `regions.csv` | The R1/R2/R3 region classification per tetramer. |
| `cluster_ids.txt` | Cluster assignments for the enriched motifs. |
| `bootstrap_<b>.tsv` / `.Rdata` | The empirical-FDR null draws (reused on resume — do not delete if re-scoring). |
| `no_enriched.txt` | Written instead of plots when **nothing** passes the filters (usually too few exons or thresholds too strict). |
| `rnamotifs.log` | Full run log. |

## Reading the RNA splicing map

The map plots, for each enriched motif, the fraction of exons carrying a motif
cluster at each position around the splice sites, as three curves: **enhanced**,
**silenced**, **control**.

- A **peak of the enhanced (or silenced) curve that rises above the control curve**
  at a given region marks where that motif's clusters concentrate to drive
  regulation.
- Position matters: e.g. pyrimidine motifs peaking in the **upstream intron** near
  **silenced** exons is the classic RNAmotifs readout of intronic silencing.
- The region labels **R1 / R2 / R3** correspond to upstream intron / exon body /
  downstream intron (skipped exon).

A motif is reported only if it beats **both** `--p-fisher` (Fisher's exact,
regulated vs control) and `--p-empirical` (bootstrap null). Loosen these (and/or
`--min-height`) if you expect signal but see `no_enriched.txt`.

## Extracting motif coordinates

```bash
./rnamotifs-extract results/MYRBP/<run>             # all enriched, TSV
./rnamotifs-extract results/MYRBP/<run> -t YCAY --bed
./rnamotifs-extract results/MYRBP/<run> -c silenced -o sil_motifs.tsv
```

| Option | Description | Default |
|--------|-------------|---------|
| `results_dir` | Path to a completed results folder (positional) | required |
| `-t, --tetramer` | Extract for a specific tetramer | all enriched |
| `-c, --category` | `enhanced`, `silenced`, or `both` | both |
| `-o, --output` | Output file | stdout |
| `--bed` | Output in BED format | TSV |


---
[← Back to tutorial index](README.md)
