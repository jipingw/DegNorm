\name{plot_coverage}
\alias{plot_coverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coverage plot functions for DegNorm
}
\description{
\code{plot_coverage} plots the before- and after-degradation coverage curves}
\usage{
plot_coverage(gene_name, coverage_output, degnorm_output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{gene_name}{the name of the gene whose coverage coverage to be plotted.}
\item{coverage_output}{CoverageClass object, the output from function 
\code{coverage_cal_batch}.}
\item{degnorm_output}{DegNormClass object, the output from function 
\code{DegNorm}.}
}
\value{
The coverage curve before and after degradation normalization.
}
\details{
\code{plot_coverage}{ outputs the coverage curves before- and after-degradation 
normalization.}
}


\examples{
## gene named "SOD1",plot coverage curves
data(coverage_res_chr21)
data(res_DegNorm_chr21)
plot_coverage(gene_name="SOD1", coverage_output=coverage_res_chr21, 
degnorm_output=res_DegNorm_chr21)

}
