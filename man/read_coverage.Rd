
\name{read_coverage}
\alias{read_coverage}
\title{
Function to calculate read coverage score for one bam file
}
\description{
This function judges whether bam file is single-end and paired-end,
and generate bam file index if needed. It calls function 
\code{paired_end_cov_by_ch} or \code{single_end_by_ch}.      
It takes multiple-core structure
for parallel computing for efficiency.}
\usage{
read_coverage(bam_file,all_genes,cores)
}
\arguments{
\item{bam_file}{The name of the bam file.}
\item{all_genes}{An GRangesList object. It's the parsed genes annotation 
file from GTF file.}
\item{cores}{number of cores to use.}
}
\details{
This function judges whether bam file is single-end and paired-end,
and generate bam file index if needed. It takes multiple-core structure
for parallel computing for efficiency.}
\value{
This function returns a coverageClass object.  It contains a list of: 
(1) a list of coverage score for each gene in RLE format and 
(2) a dataframe for read counts   
}

\seealso{
\code{\link[DegNorm]{read_coverage_batch}}
}
\keyword{internal}
