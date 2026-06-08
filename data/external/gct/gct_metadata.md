# gcT External Validation Dataset

Source: `gcT` dataset from the R package `seqinr` version 4.2-36.

Primary citation:
Galtier, N. & Lobry, J.R. (1997). Relationships between genomic G+C content,
RNA secondary structures, and optimal growth temperature in prokaryotes.
Journal of Molecular Evolution 44:632-636.

Dataset notes:
- `gct_species.csv` contains 772 bacterial species with optimal growth
  temperature (`Topt`) and genomic G+C content (`GCw`).
- RNA tables contain stem G+C fractions from the original `gcT` object and are
  not equivalent to the locally extracted `GC_tRNA`, `GC_16S`, or `GC_23S`
  features used in the TEMPURA/NCBI pipeline.
- This dataset is used only as external validation and is not merged into
  `data/processed/real_gc_data.csv`.

Download source:
https://cran.r-project.org/src/contrib/Archive/seqinr/seqinr_4.2-36.tar.gz
