# Input Format

RNAmotifs takes a **semicolon-delimited** text file, one row per exon, **no header**:

```
row_id;second_id;chrom;strand;upstream_exon_end;exon_start;exon_end;downstream_exon_start;dIRank
```

Worked example (skipped-exon / SE events):

```
1;1391;chr12;-;131286033;131288917;131289088;131291526;-1
2;19375;chr3;+;47899002;47908735;47908828;47912302;1
3;7781;chr1;+;200100;200250;200320;201000;0
```

| Field | Meaning |
|------|---------|
| `row_id` | Sequential index (1..N) |
| `second_id` | Any stable event ID (e.g. the rMATS event ID) |
| `chrom` | Chromosome (e.g. `chr12`) |
| `strand` | `+` or `-` |
| `upstream_exon_end` | Genomic end of the upstream constitutive exon |
| `exon_start` | Start of the **alternative** (cassette) exon (0-based) |
| `exon_end` | End of the alternative exon |
| `downstream_exon_start` | Start of the downstream constitutive exon |
| `dIRank` | Regulation class (see below) |

Coordinates are **ascending genomic** regardless of strand
(`upstream_exon_end < exon_start < exon_end < downstream_exon_start`).

---

## `dIRank` — the regulation class (read this carefully)

RNAmotifs classifies each exon **by the value in the last column**:

| dIRank | Class | Meaning |
|--------|-------|---------|
| `≥ 1`  | **enhanced** | inclusion promoted by the RBP |
| `≤ -1` | **silenced** | inclusion repressed by the RBP |
| `-0.1 … 0.1` | **control** | constitutive / not regulated |
| `0.1 … 1` or `-1 … -0.1` | *ignored* | between the regulated and control bands |

> ⚠️ **The #1 pitfall.** `dIRank` is a **categorical regulation label**, *not* a raw
> ΔΨ. The thresholds are `|dIRank| ≥ 1` (regulated) and `|dIRank| ≤ 0.1` (control).
> If you write the real percent-spliced-in difference (ΔΨ ∈ (-1, 1)) here, **every
> regulated exon falls below the threshold and is treated as control** — the run
> finds nothing. Use the **categorical convention** `{-1, 0, +1}` (silenced /
> control / enhanced), or any values with `|·| ≥ 1` for regulated and `|·| ≤ 0.1` for
> control. If you have rMATS output, use `--from-rmats` (below) and let RNAmotifs
> build the labels for you.

How many exons? Enrichment statistics need power — roughly **≥ 50 regulated exons**
per direction is a sensible floor.

---

## Feeding rMATS directly (recommended)

Instead of building the file by hand, point RNAmotifs at an rMATS output and let it
apply the ΔΨ/FDR thresholds, sample controls, and write the `dIRank` labels:

```bash
./rnamotifs SE.MATS.JC.txt --from-rmats \
    --rmats-incl 0.1 --rmats-fdr 0.05 \
    --name MYRBP --genome hg19 --bootstraps 1000 --cores 8
```

| Flag | Default | Meaning |
|------|---------|---------|
| `--from-rmats` | off | Treat `input_file` as an rMATS file (SE). |
| `--rmats-incl` | 0.1 | Min `|IncLevelDifference|` (ΔΨ) for a regulated exon. |
| `--rmats-fdr` | 0.05 | Max FDR for regulated exons. |
| `--rmats-constit` | 0.01 | Max `|ΔΨ|` for an exon to count as a constitutive control. |
| `--rmats-max-constit` | 5000 | Cap on control exons (random sample if exceeded). |

---

## 🚧 Intron retention (work in progress)

Intron-retention input (`--event-type RI`, including the RI column mapping and the
`;RI` row suffix) is **implemented but not yet validated** in this release. The
coordinate semantics differ from SE (the regions map onto the retained intron). Treat
RI results as experimental until this feature is finalised. For production use, stick
to skipped exons (SE).


---
[← Back to tutorial index](README.md)
