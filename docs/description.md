# Overview

Performs Leiden clustering on single-cell RNA-seq data using principal component embeddings from Dimensionality Reduction block. The algorithm identifies groups of cells with similar gene expression patterns, enabling the discovery of distinct cell populations and biological states.

The block constructs a neighborhood graph from PCA embeddings using k-nearest neighbors, then applies the Leiden algorithm to partition cells into clusters. Leiden improves upon the Louvain method by guaranteeing well-connected communities. The resolution parameter controls clustering granularity: lower values produce fewer, larger clusters, while higher values yield more numerous, smaller clusters, enabling exploration from broad cell types to fine-grained subpopulations.

The resulting cluster assignments can be visualized on UMAP or t-SNE embeddings to assess cluster quality and spatial relationships, and are used by downstream blocks for cluster marker identification, pseudotime inference, and differential expression analysis.

The block uses scanpy v1.10.1 and leidenalg v0.10.2 for clustering. When using this block in your research, cite both the scanpy publication (Wolf et al. 2018) and the Leiden algorithm paper (Traag et al. 2019) listed below.

The following publications describe the methodologies used:

> Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. _Genome Biology_ **19**, 15 (2018). [https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
> Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. _Scientific Reports_ **9**, 5233 (2019). [https://doi.org/10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)
