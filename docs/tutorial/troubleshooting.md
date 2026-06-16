# Troubleshooting / FAQ

**"It finds nothing" / `no_enriched.txt`.**
~90% of the time this is the **`dIRank` pitfall**: real ΔΨ values in (-1, 1) are read
as control because the regulated threshold is `|dIRank| ≥ 1`. Use categorical
`{-1, 0, +1}`, or `--from-rmats`. Otherwise: lower `--min-height` (2–3), loosen
`--p-fisher`/`--p-empirical`, widen `-e`, and check you have enough exons.
See [Input Format](input-format.md#dirank--the-regulation-class-read-this-carefully).

**How many exons do I need?**
Roughly **≥ 50 regulated exons** per direction; fewer gives unstable maps.

**Why 512 motifs, not 256?**
The tetramer set includes degenerate IUPAC motifs (e.g. `Y`, `W`, `S`), not only the
256 exact 4-mers.

**Empirical p never gets small.**
The floor is `1/bootstraps`-limited. Raise `-b` (e.g. 10000) for a lower floor, but
expect the *called set* to change little. With `-b 1000`, set `--p-empirical` to a
reachable value such as `0.01`.

**"Genome not found."**
Run the matching `genomes/<assembly>.download.sh` first (`hg19`, `hg38`, `mm9`,
`mm10`, `mm39`).

**Results differ slightly between runs.**
The bootstrap RNG is seeded **per thread**, so results are reproducible **only at a
fixed `--cores`**. Change the core count and borderline calls can shift. Keep
`--cores` constant for reproducible comparisons.

**It's slow.**
Use `--bootstraps 1000` for development. Raising `--cores` helps only up to a point —
the workload is largely **memory-bandwidth-bound**, so throughput plateaus well below
high core counts.

**Too many noisy motifs.**
Raise `--min-height`, tighten `--p-fisher`/`--p-empirical`, and/or increase `-b`.

**ViennaRNA / PhyloP errors with `--structure` / `--conservation`.**
These are 🚧 **work-in-progress** features. `--structure` needs ViennaRNA installed;
`--conservation` needs pre-downloaded `.phylop.bin` files
(`genomes/download_phylop.sh`). Treat their output as experimental for now.


---
[← Back to tutorial index](README.md)
