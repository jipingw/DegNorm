\name{summary_DegNormClass}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{summary_DegNormClass}
\title{
Summary method for DegNormClass.
}
\description{
It prints a summary of the data objects contained in the list from degnorm 
function.
}
\usage{
summary_DegNormClass(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{DegNormClass from \code{degnorm} function.
}
}
\value{On-screen summary of \code{DegNormClass} object.}

\examples{
## Summary of degnorm output (DegNormlass)
data(res_DegNorm_chr21)
summary_DegNormClass(res_DegNorm_chr21)
}
