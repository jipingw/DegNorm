\name{coverage_res_chr21}
\alias{coverage_res_chr21}
\docType{data}
\title{
Example CoverageClass data 
}
\description{
Example of \code{CoverageClass} data  from DegNorm package. 
It is the output from \code{read_coverage_batch} function for human 
chromosome 21. 
}
\usage{data(coverage_res_chr21)}

\format{
A \code{coverageClass} list of the following
\describe{
    \item{\code{coverage}}{a list of converage matrices for all 
    genes within each sample}
    \item{\code{counts}}{a data.frame of read counts for all genes 
    within each sample.}
}
}

\examples{
data(coverage_res_chr21)
summary_CoverageClass(coverage_res_chr21)
}
\keyword{datasets}
