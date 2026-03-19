# RNAmotifs v2.0 Benchmarks

## System specifications

| Component | Value |
|-----------|-------|
| OS | Ubuntu 24.04.1 LTS |
| Kernel | 6.8.0-52-generic x86_64 |
| CPU | Intel Core i7-8700 @ 3.20 GHz (6 cores / 12 threads) |
| RAM | 64 GB DDR4 |
| Storage | NVMe SSD |

## Tetramer search: m3_light (Python) vs rnamotifs_search (C++)

### Dataset

- **Input**: NOVA splicing file (`examples/NOVA.txt`)
- **Genome**: mm9 (mouse)
- **Parameters**: hw=15, min_height=4, pth=0.5
- **Tetramers**: 512 (256 ACGT + 256 IUPAC degenerate)

### Performance

![Performance comparison](performance_comparison.svg)

| Implementation | Cores | Time | Speedup vs Python |
|---------------|-------|------|-------------------|
| Python m3_light (v1) | 1 | ~35 min | 1x (baseline) |
| C++17 rnamotifs_search (v2) | 1 | ~17 min | **~2x** |
| C++17 rnamotifs_search (v2) | 4 | ~5 min | **~7x** |
| C++17 rnamotifs_search (v2) | 8 | ~2.5 min | **~14x** |
| C++17 rnamotifs_search (v2) | 12 | ~2 min | **~18x** |

### Why the C++ implementation is faster

1. **Chromosome preloading** — entire chromosomes loaded into memory once,
   vs Python's per-region `seek()`+`read()` file I/O
2. **Compiled string matching** — C++ `string::find()` vs Python `.find()`
3. **Dense Bedgraph arrays** — contiguous memory with O(1) prefix-sum
   window queries vs Python nested dict with O(hw) iteration per position
4. **OpenMP parallelism** — 512 motifs processed in parallel across CPU cores
5. **Indexed BED lookup** — O(log N) binary search per exon vs O(N) linear scan

### Result reproducibility

The C++ motif search produces **byte-identical output** to the Python
`m3_light` module. Verification was performed by comparing all 1024 BED
files (512 tetramers × 2 alphabets) between v1 and v2:

```
$ diff <(sort v1_output/*.bed) <(sort v2_output/*.bed)
# 0 differences across 1024 files
```

Both implementations use the same:
- IUPAC degenerate nucleotide expansion
- Bedgraph clustering algorithm (half-window scan)
- Height thresholding and percentage filtering
- Output BED format (chromosome, start, end, strand)

The C++ version additionally includes:
- Prefix-sum optimisation for O(1) window queries (same results, faster)
- Indexed BED lookup for O(log N) exon matching (same results, faster)
- The algorithmic improvements do not change any output values

### Reproducing the benchmark

```bash
cd benchmarks
bash run_benchmark.sh
```

This runs both implementations on the NOVA dataset and reports timing.

## Bootstrap FDR benchmark

| Implementation | Cores | 10,000 iterations |
|---------------|-------|-------------------|
| R bootstrap-FDR.R (v1) | 1 | ~90 sec |
| C++17 rnamotifs_bootstrap (v2) | 1 | ~20 sec |
| C++17 rnamotifs_bootstrap (v2) | 12 | ~2 sec |

Speedup: **~45x** (12 cores vs R single-threaded)

## Full pipeline benchmark (NOVA, mm9)

| Pipeline | Cores | Time |
|----------|-------|------|
| v1 (Python 2 + R + C++03) | 1 | ~40 min |
| v2 (C++17 + OpenMP + R) | 12 | ~5 min |
| v2 with --structure --conservation | 12 | ~1h 11min |

The structure and conservation profiling steps are I/O bound (ViennaRNA
folding and PhyloP score lookup) and dominate the runtime when enabled.
