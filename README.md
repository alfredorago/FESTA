# FESTA
Flexible Exon-based Splicing and Transcription Annotation

The FESTA script allows the detection of new splicing events starting only from exon-based expression data.
The only input files required are a count table and an annotation table that matches each exon to a gene.
The output contains expression values for both the main gene and its isoform.
The latter can be further normalized to proportions of the total gene expression.
For more information, refer to the manuscript in bioRxiv.

# Installation
The package can be installed from R via devtools:
```R

install.packages("devtools")
devtools::install_github("alfredorago/FESTA")
```

# Usage
Typical usage includes running the function `ClusterExons` followed by `AverageExons`. 

`ClusterExons` sorts exons within a gene in either constitutively expressed or alternatively expressed groups.

`AverageExons` calculate the expression of each exon group for downstream analysis.
