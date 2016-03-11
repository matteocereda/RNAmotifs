## RNAmotifs v.2.0 ##


This is the version 2.0 of [RNAmotifs](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-1-r20]).

RNAmotifs is an integrated *R, python, C++ software* that evaluates the sequence around differentially regulated alternative exons to identify clusters of short sequences, referred to as multivalent RNA motifs, bound by RNA-binding proteins. From a list of alternatively spliced exons, RNAmotifs identifies clusters of short non-degenerate or degenerate tetramers that are enriched at specific positions around the enhanced and silenced exons. Moreover, RNAmotifs generates the RNA splicing map of entiched motifs. RNAmotifs has been successfully used for the identification of multivalent RNA motifs bound by RNA-binding proteins such as NOVA, PTBP1, hnRNP C, TARDBP, and TIA1 and TIAL1.

RNAmotifs is based on [GeCo++](http://bioinformatics.oxfordjournals.org/content/27/9/1313.long), a C++ class library that provides a class hierarchy for the development of bioinformatic algorithm when annotations of genomic elements (e.g. binding sites, mutations) to sequences are taken into account. GeCo++ is available @ http://bioinformatics.emedea.it/geco/


Please cite:

**Cereda M**, Pozzoli U, Rot G, Juvan P, Schweitzer A, Clark T, Ule J. *RNAmotifs: prediction of multivalent RNA motifs that control alternative splicing.* Genome Biol. 2014 Jan 31;15(1):R20. doi: 10.1186/gb-2014-15-1-r20. PMID: [24485098](http://www.ncbi.nlm.nih.gov/pubmed/24485098)


## Installation

- build GeCo++ software
```
cd gMotifs
mkdir build
cd build
cmake ../
make
```
- add root folder to Python path:
```
export PYTHONPATH=$PYTHONPATH:/path_to_RNAmotifs_root_folder
```

## Contributors

RNAmotifs has been designed by Dr *Matteo Cereda* and Prof Jernej Ule. 
Main developer: Matteo Cereda. Contributing developers Dr Gregor Rot and Dr Uberto Pozzoli.

Contributions are always welcome!

## License

Please read the [Licence](LICENSE) first