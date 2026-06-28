# Data preprocessing methods

## eCLIP data processing

Enhanced crosslinking and immunoprecipitation (eCLIP) data for 28 RNA-binding
proteins (15 in HepG2, 13 in K562) were obtained from the ENCODE portal
(https://www.encodeproject.org). For each RBP, aligned BAM files from
replicate eCLIP experiments were downloaded using the ENCODE REST API.

### Crosslink site identification

Crosslink sites were extracted from aligned BAM files as follows. For
paired-end eCLIP libraries, only the second-in-pair reads (SAM flag 128) were
retained, as these carry the crosslink information at their 5' end. Reads with
mapping quality below 10 or unmapped reads (SAM flag 4) were excluded. The
crosslink position was defined as the nucleotide immediately upstream of the
read 5' end: for forward-strand reads, position = alignment start − 1; for
reverse-strand reads, position = alignment start + read length − 1. Crosslink
sites were output in BED6 format (chromosome, start, end, name, score, strand)
and sorted by genomic coordinate.

### Peak calling and merging

In the absence of a pre-computed iCount segmentation file, crosslink sites
were merged into peaks using BEDTools v2.31.1 (Quinlan and Hall, 2010).
Strand-specific merging was performed with a maximum distance of 3 nt between
adjacent crosslinks (`bedtools merge -s -d 3`), producing peaks at 3-nucleotide
resolution. Crosslink counts within each merged interval were summed. Peaks
supported by fewer than 2 crosslinks were discarded.

For each RBP, peaks from biological replicates were merged by concatenating all
replicate peak files, sorting by genomic coordinate, and performing a final
strand-specific merge with `bedtools merge -s -d 0`, summing scores across
overlapping intervals.

### Genome assembly liftOver

Processed peak files were generated for both GRCh37/hg19 and GRCh38/hg38
genome assemblies. When the source BAM files were aligned to hg19, peaks were
lifted to hg38 using the UCSC liftOver tool (Hinrichs et al., 2006) with the
`hg19ToHg38.over.chain.gz` chain file. Conversely, GRCh38-aligned data were
lifted to hg19 using `hg38ToHg19.over.chain.gz`. Chain files were obtained from
the UCSC Genome Browser download server
(https://hgdownload.cse.ucsc.edu/goldenpath/). After liftOver, peaks were
re-sorted and merged to remove any overlapping intervals introduced by
coordinate conversion. Unmapped positions were discarded.

### Output format

The final processed files follow the naming convention
`{RBP}.merged.3nt_peaks_liftOver_ordered_merged.bed` and contain three-column
BED records (chromosome, start, end), consistent with the format used by the
RNAMaRs framework (Becchi, Boscagli, Grieco et al.).

### Data quality note

HNRNPU eCLIP data for HepG2 and K562 were processed from distinct ENCODE
experiments (ENCSR240MVJ and ENCSR520BZQ, respectively). The K562 HNRNPU data
were reprocessed after an initial download artefact that produced identical
files for both cell lines. The corrected files have distinct peak sets
(HepG2: 21,242,688 peaks; K562: 19,421,573 peaks).

## Differential gene expression and splicing quantification

Differential gene expression (DESeq2) and differential splicing (rMATS)
results were obtained from the ENCODE portal for 262 RBP knockdown experiments
across HepG2 and K562 cell lines. File accession identifiers and S3 download
URLs were extracted from Supplementary Table S4 of the RNAMaRs manuscript.
A total of 514 files were downloaded: 472 DESeq2 differential expression
quantification files (.tsv) and 42 rMATS differential splicing quantification
archives (.tar.gz).

DESeq2 files contain per-gene results including base mean expression,
log2 fold change, standard error, Wald statistic, p-value and adjusted p-value
(Benjamini–Hochberg). rMATS archives contain skipped exon (SE), mutually
exclusive exon (MXE), alternative 3' splice site (A3SS), alternative 5' splice
site (A5SS) and retained intron (RI) quantifications with junction count-based
significance testing.

## Exon classification for RNAmotifs input

Alternative and constitutive cassette exons were classified from rMATS skipped
exon (SE) output following the criteria described in the RNAMaRs manuscript:

1. **Exon inclusion quantification.** Exon inclusion was quantified using the
   percent spliced-in (Ψ) metric from rMATS junction count output
   (SE.MATS.JunctionCountOnly.txt). The inclusion level difference
   (IncLevelDifference = Ψ_KD − Ψ_control) represents the change in exon
   inclusion upon RBP knockdown relative to control.

2. **Alternative exon identification.** Cassette exons with significant changes
   in average Ψ (|ΔΨ| > 0.1 and FDR < 0.1) were classified as alternatively
   spliced. The direction of regulation was assigned based on the ENCODE shRNA
   convention: IncLevelDifference > 0 (increased inclusion upon knockdown)
   indicates that the RBP normally silences the exon; IncLevelDifference < 0
   (decreased inclusion upon knockdown) indicates that the RBP normally enhances
   exon inclusion.

3. **Cassette exon annotation.** Alternatively spliced exons were cross-referenced
   with the UCSC Genome Browser hg19 knownAlt table (Navarro Gonzalez et al.,
   2021) to retain only exons annotated as "cassetteExon" (164,453 annotations).

4. **RBP binding evidence.** Exons with at least one eCLIP iCount peak within
   300 nucleotides into introns and 30 nucleotides into exons from the splice
   sites were considered as having RBP binding evidence. For introns shorter
   than 600 nt or exons shorter than 60 nt, the full length was evaluated.

5. **Constitutive exon selection.** Exons with |ΔΨ| ≤ 0.01 were classified as
   constitutive. When more than 5,000 constitutive exons were available, a
   random sample of 5,000 was selected (seed = 42) for computational
   efficiency.

6. **RNAmotifs input format.** Classified exons were written in the RNAmotifs
   semicolon-delimited format: row_id; event_id; chromosome; strand;
   upstream_exon_end; exon_start; exon_end; downstream_exon_start; dIRank.
   Enhanced exons received positive dIRank values (1 + |ΔΨ| × 10), silenced
   exons received negative values (−(1 + |ΔΨ| × 10)), and constitutive exons
   received dIRank = 0.

## Genome reference data

Chromosome sequences for hg19 (GRCh37) and mm9 (NCBI37) were obtained from
the UCSC Genome Browser in 2bit format and converted to per-chromosome plain
text files using the UCSC twoBitToFa utility. Each chromosome file contains
the complete nucleotide sequence without headers or line breaks, enabling
direct random-access indexing by genomic position during the tetramer search
step.

PhyloP evolutionary conservation scores were obtained from the UCSC Genome
Browser as BigWig files (100-way vertebrate alignment for hg19/hg38; 30-way
for mm9; 60-way for mm10) and converted to per-chromosome binary float32
arrays for efficient position-level access during conservation profiling.

## Software versions

| Tool | Version | Purpose |
|------|---------|---------|
| samtools | 1.19.2 | BAM processing and read filtering |
| BEDTools | 2.31.1 | Peak merging and interval operations |
| UCSC liftOver | latest | Genome assembly coordinate conversion |
| iCount-Mini | 3.0.1 | Installed for crosslink analysis (fallback to BEDTools used) |
| R | 4.x | Statistical analysis and visualization |
| rMATS | (ENCODE) | Differential splicing quantification (pre-computed) |
| DESeq2 | (ENCODE) | Differential gene expression (pre-computed) |

## References

Hinrichs, A. S. et al. The UCSC Genome Browser Database: update 2006.
*Nucleic Acids Res.* **34**, D590–D598 (2006).

Navarro Gonzalez, J. et al. The UCSC Genome Browser database: 2021 update.
*Nucleic Acids Res.* **49**, D1046–D1057 (2021).

Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for
comparing genomic features. *Bioinformatics* **26**, 841–842 (2010).

Cereda, M. et al. RNAmotifs: prediction of multivalent RNA motifs that
control alternative splicing. *Genome Biol.* **15**, R20 (2014).
