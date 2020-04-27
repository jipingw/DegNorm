\name{summary_CoverageClass}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{summary_CoverageClass}
\title{
Summary method for CoverageClass.
}
\description{
It prints a summary of the data objects contained in the list from 
\code{read_coverage_batch}.
}
\usage{
summary_CoverageClass(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{CoverageClass from code{read_coverage_batch.}
}}
\value{On-screen plot of summary of \code{CoverageClass} object.}
\examples{
## Summary of coverage_cal_batch output (CoverageClass)
data(coverage_res_chr21)
summary_CoverageClass(coverage_res_chr21)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
