# RNAmotifs

Discover clusters of short RNA motifs (tetramers) that are **positionally enriched**
around alternatively spliced exons regulated by an RNA-binding protein (RBP), and
read the resulting **RNA splicing maps**.

> Method reference: Cereda M. *et al.* **RNAmotifs: prediction of multivalent RNA
> motifs that control alternative splicing.** *Genome Biology* 2014;15(1):R20.

## Tutorial pages
- **[Input Format](input-format.md)** — file layout and the critical `dIRank` convention
- **[Parameters](parameters.md)** — every flag explained, with tuning guidance
- **[Output & Interpretation](output.md)** — what each file is, reading the splicing map
- **[Examples](examples.md)** — copy-paste command recipes
- **[Troubleshooting](troubleshooting.md)** — FAQ and common pitfalls
- **[RNAmotifs-MaRs](rnamotifs-mars.md)** — coupling with eCLIP for MRM–RBP association scores

> 🚧 Intron retention (`--event-type RI`), secondary structure (`--structure`) and
> conservation (`--conservation`) are **work-in-progress** features — see the notes on
> the relevant pages.

---

## What RNAmotifs does

When you knock down an RBP and run RNA-seq, some cassette exons are **included more**
(enhanced) and some **less** (silenced). RNAmotifs asks:

> *Which short sequence motifs, clustered at which positions around those exons,
> distinguish the regulated exons from unaffected (control) exons?*

It takes three exon sets — **enhanced**, **silenced**, **control** — scans the
sequence around each exon's splice sites for every k-mer (tetramers by default),
finds positions where a motif forms a **multivalent cluster** (several copies close
together), and tests whether each motif's clusters are **enriched** in the regulated
sets vs control. The output is an **RNA splicing map** plus a statistically filtered
list of enriched motifs.

---

## 2-minute quick start

```bash
# Build binaries (see README "Installation"):  cmake -B build && cmake --build build
# Download your genome:                         bash genomes/hg19.download.sh

./rnamotifs my_exons.txt \
    --name MYRBP \
    --genome hg19 \
    --bootstraps 1000 \
    --cores 8
```

Writes a results folder under `results/MYRBP/<timestamped-run>/` with the RNA
splicing-map PDF, the enriched-tetramer table, and per-region/per-group CSVs (see
[Output & Interpretation](output.md)).

> **Tip:** use `--bootstraps 1000` for everyday runs (fast). `--bootstraps 10000` is
> ~10× slower for a lower empirical-p floor that rarely changes the called set.

---

## How the algorithm works (in five steps)

Understanding these makes every [parameter](parameters.md) intuitive.

1. **Define regions** around each exon's two splice sites — into the exon
   (`--in-exon` bp) and into the flanking introns (`--in-intron` bp). For a skipped
   exon: **R1** = upstream intron, **R2** = exon body, **R3** = downstream intron.
2. **Find motif clusters.** For each k-mer, slide along each region and count copies
   within a **half-window** (`--half-window`). Where the local count (cluster
   *height*) reaches `--min-height`, the position is "occupied" — this is what
   *multivalent* means.
3. **Positional occupancy per group** — fraction of enhanced / silenced / control
   exons with a cluster at each position (`--pth` binarises occupancy). Plotting this
   *is* the splicing map.
4. **Enrichment test** — Fisher's exact test, regulated vs control, per region;
   must beat `--p-fisher`.
5. **Empirical FDR** — resample control exons `--bootstraps` times for an empirical
   p-value; must also beat `--p-empirical`.

So `hw`/`min-height`/`pth` define *what counts as a cluster*; `in-exon`/`in-intron`
define *where you look*; `p-fisher`/`p-empirical`/`bootstraps` define *how strict the
filter is*.

---

## Supported genomes

| Genome | Species | Assembly |
|--------|---------|----------|
| `hg19` | Human | GRCh37 |
| `hg38` | Human | GRCh38 |
| `mm9`  | Mouse | NCBI37 |
| `mm10` | Mouse | GRCm38 |
| `mm39` | Mouse | GRCm39 |

The matching genome must be downloaded first (`genomes/<assembly>.download.sh`).
