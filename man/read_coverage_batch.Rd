\name{read_coverage_batch}
\alias{read_coverage_batch}
\title{
Compute the read coverage score and read counts for all genes in batch mode
}
\description{
This function calls \code{read_coverage} to compute read coverage score and 
read counts for all genes and samples.

Notes: 
\enumerate{
  \item Coverage score is calculated per gene, i.e., concatenation of all
  exons from the same gene.
  \item We follow the HTseq protocol for counting valid reads or read pairs for each 
  gene.
  \item When reading the alignment file, \code{isSecondaryAlignment} is set to 
  \code{FALSE} to avoid redundant counting.
  \item For paired-end data, \code{isPaired} is set to \code{TRUE}. We do not recommend
  setting \code{isProperPair} to \code{TRUE} as some fragment lengths may exceed 200 bp.
  \item Users can modify \code{scanBamParam} in the R code as needed.
}
}
\usage{
read_coverage_batch(bam_file_list, gtf_file, cores = 1)
}
\arguments{
\item{bam_file_list}{A character vector of BAM file paths.}
\item{gtf_file}{The GTF file that RNA-seq reads were aligned to.}
\item{cores}{Number of processor cores to use. Default = 1.}
}
\value{
A list with the following components:
\item{coverage}{A list of coverage matrices for all genes within each sample.}
\item{counts}{A \code{data.frame} of read counts for all genes within each sample.}
}
\examples{
## Read BAM file and GTF file from the package
bam_file_list <- list.files(
    path = system.file("extdata", package = "DegNorm"),
    pattern = ".bam$", full.names = TRUE
)
gtf_file <- list.files(
    path = system.file("extdata", package = "DegNorm"),
    pattern = ".gtf$", full.names = TRUE
)

## Run read_coverage_batch to calculate read coverage curves and read counts
coverage_res <- read_coverage_batch(bam_file_list, gtf_file, cores = 2)
}
\seealso{
\code{\link[DegNorm]{read_coverage}}
}
