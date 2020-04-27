\name{DegNorm-package}
\alias{DegNorm-package}
\docType{package}
\title{DegNorm: degradation normalization for RNA-seq data}
\description{\code{DegNorm} is an R package for degradation normalication 
for bulk RNA-seq data.DegNorm, short for degradation normalization, is a 
bioinformatics pipeline  designed to correct for bias due to the heterogeneous 
patterns of transcript degradation in RNA-seq data. 
}
\details{
DegNorm is a data-driven approach for RNA-Seq normalization resulting in 
the adjusted read count matrix. This adjustment applies to each gene within 
each sample, accounting for sample- and gene-specific degradation bias while 
simultaneously controlling for the sequencing depth. The algorithm at the 
center of DegNorm is the rank-one over-approximation of a gene's coverage 
score matrix, which is comprised of the different samples' coverage score 
curves along the transcript for each gene. For each gene, DegNorm estimates 
(1) an envelope function representing the ideal shape of the gene's coverage
curve when no degradation is present, and (2) scale factors for each sample
(for said gene) that indicates the relative abundance of the gene within the 
sample.

functions:
\code{read_coverage_batch,
degnorm,
plot_coverage,
plot_heatmap,
plot_corr,
plot_boxplot}
}
\author{
Bin Xiong, Ji-Ping Wang

Maintainer: Ji-Ping Wang <jzwang@northwestern.edu>
}
\references{
DegNorm reference: 

Xiong, B., Yang, Y., Fineis, F. Wang, J.-P., DegNorm: normalization of 
generalized transcript degradation improves accuracy in RNA-seq analysis,
Genome Biology, 2019,20:75
}
\keyword{ RNA-seq, degradation, normalization }
