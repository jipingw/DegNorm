\name{degnorm}
\alias{degnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to perform degradation normalization.
}
\description{
\code{degnorm} calcualtes the degradation index score for each 
gene within each sample and return the degradation-normalized read counts.
}
\usage{
degnorm(read_coverage,counts,iteration,loop,down_sampling=1,grid_size=10,
cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{read_coverage}{a list of converage matrices, one per gene}
    \item{counts}{dataframe of read counts, each row for one gene, and column
    for sample. The order and number of genes must match the order in
    read_coverage matrices.}
    \item{iteration}{iteration number for degnorm algorithm. 5 is sufficient.}
    \item{loop}{iteration number inside of nonnegative matrix 
    factorization-over approximation. Default is 100.}
    \item{down_sampling}{1 for yes (default) and 0 for no. If yes, average
    coverage score is calcualted on a grid of size specified by grid_size 
    argument. The new coverage matrix formed by the grid average score will 
    be used for baseline selection. This increases the efficiency of algorithm 
    while maintaining comparable accuracy.}
    \item{grid_size}{default size is 10 bp.}
    \item{cores}{number of cores. Default number if 1. Users should input 
    the maximum possible number of cores for efficiency.}
}

\value{
\code{degnorm} outputs a list of following objects:
\item{counts}{a data.drame of read counts for each gene within each sample.}
\item{counts_normed}{a data.drame of degradation-normalized read counts for
each gene within each sample.}
\item{DI}{a matrix of degradation index scores for each gene within each 
sample.}
\item{K}{normalizing scale factor for each gene within each sample after
accounting for degradation normalization.}
\item{convergence}{convergence tag; 0 = degnorm was not done on this gene 
because smaller counts or too short length.1 = degnorm was done with baseline 
selection. 2 = degnorm done without baseline selection because gene length 
(after filtering out low count regions)<200 bp. 3= baseline was found, but 
DI score is too large. 4 = baseline selection didn't coverge. }
\item{envelop}{list of the envelop curves for all genes.}
}


\examples{
##coverage_res_chr21 is a \code{CoverageClass} object from DegNorm Package.
data(coverage_res_chr21)
res_DegNorm = degnorm(read_coverage = coverage_res_chr21[[1]],
                    counts = coverage_res_chr21[[2]],
                    iteration = 2,
                    down_sampling = 1,
                    grid_size=10,
                    loop = 20,
                    cores=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
