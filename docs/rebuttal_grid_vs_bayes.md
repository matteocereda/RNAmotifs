# Rebuttal — Grid vs Bayesian Optimisation (for Reviewer 2)

## Suggested addition to the existing Response to Reviewer 2

> **Reviewer comment (implicit/anticipated):** Why use a discrete grid search rather than a more sophisticated continuous optimisation strategy such as Bayesian optimisation?

**Response:**

To evaluate whether the discrete grid search used for RNAmotifs parameter selection could be improved by a continuous optimisation strategy, we benchmarked grid search against Bayesian optimisation (BO) across the full reference panel (15 RBPs in HepG2, 13 in K562; 28 RBP–cell line combinations in total). BO employed a Gaussian process surrogate with a Matérn-5/2 kernel and Expected Improvement acquisition, evaluating 38 parameter combinations per RBP (8 Latin hypercube initial points + 30 sequential iterations) over a continuous search space (n ∈ [5, 35], e ∈ [30, 300]), compared to the grid's 20 discrete combinations.

Across all 28 comparisons, grid search achieved a mean AUROC of 0.745 versus 0.705 for BO (mean Δ = −0.039 ± 0.108). Grid search outperformed or matched BO in 19/28 cases (68%), while BO outperformed grid in 9/28 cases (32%). Notably, in all 28 cases the BO optimum fell at parameter values outside the grid's discrete set, yet this finer-grained exploration did not translate into improved performance.

The relative performance of the two strategies was not random, but instead tracked the strength of the underlying eCLIP regulatory signal as independently characterised by Van Nostrand et al. (Nature, 2020; Extended Data Fig. 6). RBPs at the top of the splicing map hierarchy — those with strong, position-specific eCLIP enrichment at regulated exons (e.g., HNRNPC, PTBP1, HNRNPK, HNRNPU) — achieved high AUROC with both methods (mean AUROC > 0.95), with BO providing marginal gains in some cases (e.g., PTBP1 in HepG2: +0.11; HNRNPK in HepG2: +0.13). By contrast, RBPs at the bottom of the hierarchy — core spliceosome components and factors with diffuse or absent position-specific binding (e.g., RBM22, SF3A3, U2AF1, U2AF2) — had low AUROC with both methods, and BO consistently underperformed grid (e.g., RBM22 in HepG2: −0.32; SF3B4 in K562: −0.24). For these RBPs, the AUROC landscape is essentially flat, and the GP surrogate overfits noise, converging to spurious optima.

This pattern was consistent across both cell lines: HNRNPU achieved AUROC = 1.0 in both HepG2 and K562 with both methods, while RBPs with low AUROC in one cell line (e.g., FXR1, RBM15) similarly showed near-zero AUROC in the other, regardless of optimisation strategy. This cross-cell-type consistency is in line with the high same-RBP correlations reported in Van Nostrand et al. (Extended Data Fig. 6d).

Furthermore, BO required approximately 10 hours of computation per RBP (38 evaluations × ~15 min each) compared to ~1 hour for grid search (20 evaluations), representing a ~10× increase in computational cost with no systematic improvement in performance.

We therefore retain grid search as the default optimisation strategy in RNAMaRs. Grid search is robust, interpretable, and computationally efficient. The cases where BO outperforms grid correspond to RBPs that already achieve high AUROC with grid search, making the marginal gains practically irrelevant for downstream RBP prioritisation. Conversely, for RBPs where performance is poor, the limitation lies in the biology (absence of position-specific binding signal) rather than in the optimisation strategy.

Per-RBP convergence plots and the complete benchmark results are provided in Supplementary Figure [X] and Supplementary Table [X].


---

# Results section — Suggested addition/modification

## Insert after paragraph [58] (parameter optimization paragraph), or as a new paragraph within the Discussion

To assess whether a continuous optimisation strategy could improve upon the grid search, we compared grid search with Bayesian optimisation (BO) across the full reference panel (28 RBP–cell line combinations; Methods). Grid search outperformed or matched BO in 19/28 cases (68%), achieving a higher mean AUROC (0.745 vs 0.705; mean Δ = −0.039 ± 0.108). The relative performance of the two strategies was predictable from the strength of the underlying eCLIP signal: RBPs with strong position-specific binding (HNRNPC, PTBP1, HNRNPK, HNRNPU) achieved high AUROC (>0.95) with both methods, while core spliceosome components with diffuse binding (RBM22, SF3A3, U2AF1/2) had low AUROC regardless of strategy (Supplementary Fig. [X]). This pattern was consistent across HepG2 and K562, in agreement with the cross-cell-type conservation of splicing regulatory architectures reported in Van Nostrand et al. [7]. BO required ~10× more computation per RBP without systematic improvement, supporting grid search as the robust and efficient default for parameter selection in RNAMaRs.


---

# Methods section — Fill-in for docs/methods_bayes_opt.md placeholders

Use these values in the final paragraph of methods_bayes_opt.md:

- Number of RBPs: 28 (15 HepG2 + 13 K562)
- Mean bayes AUROC: 0.705
- Mean grid AUROC: 0.745
- Mean Δ: −0.039 ± 0.108
- RBPs finding optima outside grid: 28/28
- Mean wall-time ratio: ~10× (38 vs 20 evaluations; ~560 min vs ~50 min per RBP)
