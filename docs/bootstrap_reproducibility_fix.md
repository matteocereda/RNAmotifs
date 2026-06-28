## Bootstrap FDR reproducibility fix

### Bug

The bootstrap FDR step in the RNAmotifs pipeline (`src/cpp/bootstrap_fdr.cpp`) produced non-deterministic results across runs, even with the same input data and parameters. This made AUROC values non-reproducible and prevented fair comparison between grid search and Bayesian optimisation discovery results.

### Root cause

The bootstrap loop was parallelised with OpenMP using `schedule(dynamic, 1)`:

```cpp
mt19937 rng(42 + tid * 10007);
// ...
#pragma omp for schedule(dynamic, 1)
for (int b = 0; b < n_boot; ++b) { ... }
```

Each thread maintained its own Mersenne Twister PRNG seeded deterministically by thread ID (`42 + tid * 10007`). However, `schedule(dynamic)` assigns loop iterations to threads on a first-come-first-served basis at runtime. This means iteration *b* may be executed by thread 0 in one run and thread 3 in another, drawing from a different PRNG stream each time. The resulting bootstrap resamples — and therefore the set of tetramers surviving FDR — changed between runs.

Because the downstream AUROC is computed over the surviving tetramer set, which is typically small (5–30 tetramers), even a single tetramer difference can shift the AUROC substantially. For example, HNRNPC achieved AUROC = 1.0 at all 20 grid points in the original grid run, but the same parameters evaluated in a subsequent Bayesian optimisation run yielded AUROC values between 0.93 and 0.99 — a discrepancy caused entirely by bootstrap non-reproducibility.

### Fix

Two changes in `src/cpp/bootstrap_fdr.cpp`:

1. **Deterministic scheduling.** Changed `schedule(dynamic, 1)` to `schedule(static)`, which assigns iterations to threads at compile time in fixed-size contiguous blocks. This guarantees that iteration *b* is always executed by the same thread regardless of runtime conditions.

2. **Explicit seed.** Changed the base seed from the ad-hoc value `42` to `30580` to clearly distinguish reproducible runs from earlier non-reproducible ones.

```cpp
// Before (non-reproducible)
mt19937 rng(42 + tid * 10007);
#pragma omp for schedule(dynamic, 1)

// After (reproducible)
mt19937 rng(30580 + tid * 10007);
#pragma omp for schedule(static)
```

### Impact

- **Grid search results** must be regenerated with the fixed binary, as the original results used the non-reproducible bootstrap.
- **Bayesian optimisation benchmark** (grid vs bayes comparison) is now valid: both methods see the same deterministic objective function for any given (hw, ew) pair.
- **No change to the statistical procedure itself** — the number of bootstraps, the Fisher exact test, and the BH correction are all unchanged. Only the assignment of PRNG streams to iterations is made deterministic.

### Validation

Reproducibility was verified by running the full HepG2 discovery panel (15 RBPs) through both grid search and Bayesian optimisation with the fixed bootstrap, and confirming that identical (hw, ew) parameters produce identical AUROC values. Results for HNRNPC and PTBP1 are reported in `bench/bayes_vs_grid.md`.

---

## Final step (2026-06-26): core-count-independent seeding (release)

The per-thread seed above (`30580 + tid*10007`, `schedule(static)`) is deterministic at a
*fixed* core count but the raw bootstrap output still depends on the thread count. The
released version seeds the RNG **per iteration** instead — `mt19937(30580 + b*2654435761)`
created inside the parallel loop — so iteration *b* draws the same stream at any core count.
Verified: full pipeline re-run under the new seeding reproduces every grid/Bayes/ablation/
cross-cell/LORO/k-mer result of the previous (per-thread) run exactly; the change only makes
the *intermediate* bootstrap bit-identical across hardware. This is the mainline as of the release.
