\name{read_coverage_batch}
\alias{read_coverage_batch}
\title{
Compute the read coverage score and read counts for all genes in batch mode.
}
\description{
This function calls \code{read_coverage} to compute read coverage socre and 
read counts for all genes and samples. 

Notes: 
1. Coverage score is calcualted per gene, i.e. concatenation of all
exons from the same gene.

2. We follow HTseq protocol for counting valid read or read pairs for each 
gene.

3. When reading alignment file, \code{isSecondaryAlignment} flag is set as 
FALSE to avoid possible redundant counting.

4. For paired-end data, isPaired is set as TRUE. We don't recommend
setting isProperPair as TRUE as some fragments length may exceed 200bp.

5. User can modify scanBamParam in the R codes below as needed.
}
\usage{
read_coverage_batch(bam_file_list,gtf_file,cores=1) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{bam_file_list}{a character vector of bam file names.}
\item{gtf_file}{the gtf file that RNA-seq reads were aligned with
reference to.}
\item{cores}{number of cores to be used. Default=1.}
}

\value{
A list of the following:
\item{coverage}{a list of converage matrices for all genes within each sample.}
\item{counts}{data.frame of read counts for all genes within each sample.}
}

\examples{
## read bam file and gtf file from the package
bam_file_list <- list.files(path=system.file("extdata",package="DegNorm")
    ,pattern=".bam$",full.names=TRUE)
gtf_file <- list.files(path=system.file("extdata",package="DegNorm"),
    pattern=".gtf$",full.names=TRUE)

# run read_coverage_batch to calculate read coverage curves and read counts
coverage_res=read_coverage_batch(bam_file_list, gtf_file,cores=2)
}
\seealso{
\code{\link[DegNorm]{read_coverage}}
}