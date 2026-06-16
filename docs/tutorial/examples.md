# Examples

Copy-paste recipes. See [Parameters](parameters.md) for what each flag does.

## Basic SE run (everyday settings)
```bash
./rnamotifs my_exons.txt -n MYRBP -g hg19 -b 1000 -c 8
```

## From rMATS (let RNAmotifs build the labels)
```bash
./rnamotifs SE.MATS.JC.txt --from-rmats \
    --rmats-incl 0.1 --rmats-fdr 0.05 \
    -n MYRBP -g hg19 -b 1000 -c 8
```

Stricter inclusion cutoff:
```bash
./rnamotifs SE.MATS.JC.txt --from-rmats --rmats-incl 0.15 --rmats-fdr 0.05 \
    -n MYRBP -g hg19 -b 1000 -c 8
```

## Looser clustering, broader enrichment window
For a sparse motif with dispersed sites:
```bash
./rnamotifs my_exons.txt -n MYRBP -g hg19 -w 25 -m 3 -e 100 -b 1000 -c 8
```

## High-precision final run
```bash
./rnamotifs my_exons.txt -n MYRBP -g hg19 -b 10000 --p-empirical 0.00005 -c 8
```

## Mouse (different genome)
```bash
bash genomes/mm9.download.sh           # once
./rnamotifs nova_exons.txt -n NOVA -g mm9 -b 1000 -c 8
```

## Parameter search (don't guess `hw`/`ew`)
Use the MaRs discovery mode to search parameters per RBP — see
[RNAmotifs-MaRs](rnamotifs-mars.md).

---

## 🚧 Work-in-progress examples

The following features are implemented but **not yet validated**; commands are shown
for reference only and their output should be treated as experimental:

```bash
# Intron retention (WIP)
./rnamotifs RI.MATS.JunctionCountOnly.txt --from-rmats --event-type RI \
    --rmats-incl 0.1 --rmats-fdr 0.1 -n MYRBP_RI -g hg19 -b 1000 -c 4

# Structure + conservation (WIP; need ViennaRNA / PhyloP data)
./rnamotifs nova_exons.txt -n NOVA -g mm9 \
    --structure --structure-window 31 --conservation -c 12
```


---
[← Back to tutorial index](README.md)
