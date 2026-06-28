# Parameter Optimization for the MRM-RBP Association Score Framework

## Overview

The RNAmaRs framework quantifies the association between RNA-binding proteins (RBPs) and significant RNA tetramers identified by RNAmotifs through a composite Association Score (AS). For each tetramer $t$ and RBP $r$, the score is defined as:

$$AS(t, r) = \text{SCORE1}(t, r) \times \text{SCORE2}(t, r)$$

where SCORE1 is the Signal Recovery Rate derived from eCLIP downsampling analysis, and SCORE2 is the cosine similarity between the tetramer's positional enrichment profile and the RBP's eCLIP crosslinking density profile across the splicing regulatory region. Two parameters govern the upstream motif discovery step and must be optimized per RBP: the clustering half-window $h_w$, which controls the genomic distance within which co-occurring tetramers are grouped into regulatory clusters, and the enrichment window $e_w$, which defines the region size over which positional motif frequency is tallied. Because different RBPs bind motifs of varying spatial extent and at different distances from regulated exons, no single parameter combination is universally optimal.

## Grid search design

We performed an exhaustive grid search over $h_w \in \{5, 15, 25, 35\}$ nucleotides and $e_w \in \{30, 50, 100, 200, 300\}$ nucleotides, yielding 20 parameter combinations. For each combination, the full RNAmotifs pipeline was executed: positional tetramer frequencies were computed across regulated and control exon sets, bootstrap-based FDR correction was applied to identify significantly enriched tetramers, and association scores were calculated for all RBPs in the evaluation panel. The optimal $(h_w, e_w)$ pair was then selected per RBP based on its capacity to maximize self-identification, as quantified by the area under the receiver operating characteristic curve (AUROC).

The grid search was parallelized per RBP but executed sequentially across RBPs to maintain memory safety under the constraint of large positional score matrices. A crash-safe JSON manifest tracked completed runs, enabling seamless resumption after interruption. Critically, all intermediate score matrices $S$ were persisted to disk, allowing post-hoc recalculation of evaluation metrics under different panel compositions or scoring criteria without re-running the computationally expensive motif discovery and eCLIP intersection steps.

## AUROC computation

For each parameter combination and regulatory direction $d \in \{\text{enhanced}, \text{silenced}\}$, a score matrix $S^{(d)} \in \mathbb{R}^{N \times K}$ was constructed, where $N$ is the number of panel RBPs and $K$ is the number of significant tetramers identified at that parameter setting. Each entry $S^{(d)}_{r,t}$ contains the association score $AS(t, r)$ for tetramer $t$ and RBP $r$.

To evaluate whether RBP $r^*$ is preferentially associated with its own significant tetramers, we treated the $K$ association scores of the true RBP, $\{S^{(d)}_{r^*,t}\}_{t=1}^{K}$, as the positive set and the scores of all other panel RBPs for the same tetramers as the negative set. The AUROC was computed as the Wilcoxon-Mann-Whitney statistic:

$$\text{AUROC}(r^*) = \frac{1}{n_{\text{pos}} \cdot n_{\text{neg}}} \sum_{i=1}^{n_{\text{pos}}} \sum_{j=1}^{n_{\text{neg}}} \left[ \mathbb{I}(s^+_i > s^-_j) + \frac{1}{2}\mathbb{I}(s^+_i = s^-_j) \right]$$

where $s^+_i$ and $s^-_j$ denote scores from the positive and negative sets, respectively, and $\mathbb{I}(\cdot)$ is the indicator function. The final AUROC for each parameter combination was taken as $\max(\text{AUROC}_{\text{enhanced}}, \text{AUROC}_{\text{silenced}})$, reflecting the dominant regulatory mode of the RBP. For example, a predominantly silencing factor such as PTBP1 typically achieves higher AUROC on the silenced direction, whereas an enhancing factor such as SRSF1 performs best on the enhanced direction (Figure S3a).

## Statistical considerations and panel design

### Non-independence of tetramer observations

The $K$ significant tetramers that constitute the positive set are not statistically independent. Overlapping tetramers — for instance, CUCU and UCUC, which share a three-nucleotide overlap — capture highly correlated binding signals from the same underlying sequence context. Consequently, the effective sample size $K_{\text{eff}} < K$, and the variance of the AUROC estimator is underestimated under an independence assumption. We mitigated this in two ways. First, we imposed a minimum threshold of $K \geq 5$ significant tetramers for a parameter combination to be considered valid (see below). Second, we note that the biological constraint of sequence-specific binding inherently limits the motif space: functionally relevant tetramers cluster in sequence similarity space, and this redundancy, while reducing statistical power, simultaneously reinforces the biological signal by requiring consistent enrichment across related $k$-mers.

### Panel composition effects

AUROC is a relative metric whose value depends on the composition of the negative set. Including RBPs for which no eCLIP data are available in the target cell line produces structurally zero association scores ($\text{SCORE1} = 0$ due to absence of crosslinking signal), creating degenerate negatives. These zero-score entries inflate AUROC for RBPs with strong binding signals (whose positive scores easily exceed zero) but deflate AUROC for RBPs with weak or diffuse binding through an excess of ties at zero. The recommended practice, and the one adopted throughout our analyses, is to restrict the evaluation panel to RBPs with matched eCLIP data in the same cell line, ensuring that all negative-set scores reflect genuine binding measurements rather than missing data (Figure S3b).

### Heterogeneity of the negative set

A fundamental limitation of the one-vs-rest AUROC design is that not all "other" RBPs are equally irrelevant to a given tetramer. RBPs with overlapping binding preferences — such as U2AF1 and U2AF2 at polypyrimidine tracts, or PTBP1 and HNRNPC at CU-rich elements — produce biologically meaningful cross-reactivity that the AUROC framework penalizes as false positives. This cross-reactivity places an effective ceiling on achievable AUROC that varies across RBPs depending on the uniqueness of their binding specificity within the panel. We regard this not as a deficiency but as an informative property: RBPs achieving high AUROC despite panel overlap possess genuinely distinctive regulatory signatures, while moderate AUROC for members of binding-preference families reflects the inherent ambiguity in attributing shared motifs to individual factors (Figure S3c).

### Complementary evaluation metrics

To address the limitations of AUROC in this setting, we computed two complementary metrics. Mean Reciprocal Rank (MRR) captures the average quality of the true RBP's ranking across tetramers:

$$\text{MRR}(r^*) = \frac{1}{K} \sum_{t=1}^{K} \frac{1}{\text{rank}(r^*, t)}$$

where $\text{rank}(r^*, t)$ is the rank of $r^*$ among all panel RBPs when sorted by $AS(t, \cdot)$ in descending order. The Top-1 fraction directly measures the probability of correct identification:

$$\text{Top1}(r^*) = \frac{1}{K} \sum_{t=1}^{K} \mathbb{I}\left[\text{rank}(r^*, t) = 1\right]$$

MRR is less sensitive to panel composition than AUROC because it depends on ordinal rank rather than pairwise comparisons, while Top-1 fraction provides the most interpretable measure of practical utility — the fraction of tetramers for which the correct RBP is the top-scoring candidate. Together, these three metrics (AUROC, MRR, Top-1) provide a comprehensive assessment of self-identification capacity under different assumptions about what constitutes successful identification (Figure S3d).

### Cross-cell-line validation

For RBPs with eCLIP data available in both HepG2 and K562 cell lines, we assessed the consistency of optimal parameter selections across cell lines. Concordant optimal $(h_w, e_w)$ pairs indicate that the identified parameters reflect intrinsic properties of the RBP's binding mode — spatial extent and positional preference — rather than artifacts of a particular dataset. Conversely, parameter divergence between cell lines may arise from cell-type-specific regulatory mechanisms, differences in the splicing event landscape, or variation in eCLIP library quality and depth. We observed broadly consistent parameter preferences for the majority of shared RBPs, with divergent cases enriched for RBPs exhibiting cell-type-specific expression of isoforms or cofactors that modulate binding specificity (Figure S3e).

## Minimum tetramer threshold

Parameter combinations yielding fewer than 5 significant tetramers ($K < 5$) were excluded from the optimization. With small $K$, the AUROC estimator exhibits high variance: for $K = 1$, the statistic degenerates to a single rank comparison, and for $K = 2\text{--}4$, confidence intervals span a substantial fraction of the $[0, 1]$ range. The threshold of 5 was chosen as a pragmatic lower bound that ensures a minimum of $5 \times (N-1)$ pairwise comparisons in the Wilcoxon-Mann-Whitney statistic, providing sufficient resolution to distinguish self-identification from chance performance (AUROC = 0.5). In practice, most informative parameter combinations yielded $K \gg 5$, and the threshold primarily filtered out extreme parameter settings (e.g., very small $e_w$ with large $h_w$) that produced overly stringent or overly permissive motif discovery.
