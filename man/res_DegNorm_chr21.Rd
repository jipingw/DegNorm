\name{res_DegNorm_chr21}
\alias{res_DegNorm_chr21}
\docType{data}
\title{
Example DegNormClass data 
}
\description{
Example of \code{DegNormClass} data  from DegNorm package. It is the 
output from \code{degnorm} function for human chromosome 21. 
}
\usage{data("res_DegNorm_chr21")}
\format{
A \code{DegNormClass} list of the following items:
\describe{
\item{\code{counts}}{a data.drame of read counts for each gene within 
each sample.}
\item{\code{counts_normed}}{a data.drame of degradation-normalized read 
counts for each gene within each sample.}
\item{\code{DI}}{a matrix of degradation index scores for each gene within 
each sample.}
\item{\code{K}}{normalizing scale factor for each gene within each sample 
after accounting for degradation normalization.}
\item{\code{convergence}}{convergence tag; 0 = degnorm was not done on this 
gene because smaller counts or too short length.1 = degnorm was done with 
baseline selection. 2 = degnorm done without baseline selection because 
gene length (after filtering out low count regions)<200 bp. 3= baseline 
was found, but DI score is too large. 4 = baseline selection didn't 
coverge.}
\item{\code{envelop}}{a list of the envelop curves for all genes.}
}}

\examples{
data(res_DegNorm_chr21)
summary_DegNormClass(res_DegNorm_chr21)
}
\keyword{datasets}
