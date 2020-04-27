\name{DegNorm-plot-functions}
\alias{plot_corr}
\alias{plot_heatmap}
\alias{plot_boxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Degradation index (DI) score plot functions
}
\description{
DegNorm provides three functions for visualization gene-/sample-wise 
degradation.
}
\usage{
plot_corr(DI)
plot_heatmap(DI)
plot_boxplot(DI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{DI}{a matrix or data.frame of degradation index (DI) scores with 
each row corresponding to one gene and each column for a sample.}
}
\value{
These functions return a boxplot of DI scores by sample, a heatmap 
of DIS scores of all genes in all samples and a correlation plot of 
DI scores between samples respectively.
}

\details{
\code{plot_corr}{ plots the correlation matrix of DI scores between samples.}
\code{plot_heatmap}{ plots the heatmap of DI  scores. Left is ploted in 
descending order of average DI scores of genes where each row corresponds 
to one gene. In the right plot, DI scores were sorted within each sample and 
plotted in descending order.}
\code{plot_boxplot}{plots the boxplot of DI scores by samples.}
}

\examples{
## res_DegNorm_chr21 is degnorm otuput stored in sysdata.Rda
data(res_DegNorm_chr21)
plot_boxplot(res_DegNorm_chr21$DI)
plot_heatmap(res_DegNorm_chr21$DI)
plot_corr(res_DegNorm_chr21$DI)
}
