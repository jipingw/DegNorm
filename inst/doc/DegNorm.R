## ----setup, include = FALSE---------------------------------------------------
library(DegNorm)
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("DegNorm")

## -----------------------------------------------------------------------------
## specify bam_files from RNA-seq, you should replace it by your own bam files
bam_file_list=list.files(path=system.file("extdata",package="DegNorm"),
                        pattern=".bam$",full.names=TRUE)

## -----------------------------------------------------------------------------
## gtf_file you used for RNA-seq alignment, replace it by your own gtf file
gtf_file=list.files(path=system.file("extdata",package="DegNorm"),
                    pattern=".gtf$",full.names=TRUE)

## ----eval=TRUE----------------------------------------------------------------
## calculate the read coverage score for all genes of all samples
coverage_res_chr21_sub=read_coverage_batch(bam_file_list, gtf_file,cores=2)

## ----eval=FALSE---------------------------------------------------------------
#  ## save the coverage results
#  save(coverage_res_chr21_sub,file="coverage_res_chr21_sub.Rda")

## -----------------------------------------------------------------------------
data("coverage_res_chr21")
## summarize the coverage results
summary_CoverageClass(coverage_res_chr21)

## -----------------------------------------------------------------------------
## extract coverage scores and counts from coverage_res
coverage_matrix=coverage_res_chr21$coverage
counts=coverage_res_chr21$counts

## ----eval=TRUE----------------------------------------------------------------
res_DegNorm_chr21 = degnorm(read_coverage = coverage_res_chr21[[1]],
                    counts = coverage_res_chr21[[2]],
                    iteration = 5,
                    down_sampling = 1,
                    grid_size=10,
                    loop = 100,
                    cores=2)

## ----eval=FALSE---------------------------------------------------------------
#  ## save the DegNorm results
#  save(res_DegNorm_chr21,file="res_DegNorm_chr21.Rda")

## -----------------------------------------------------------------------------
data("res_DegNorm_chr21")

## -----------------------------------------------------------------------------
## summary of the DegNorm output
summary_DegNormClass(res_DegNorm_chr21)

## -----------------------------------------------------------------------------
## extrac normalized read counts
counts_normed=res_DegNorm_chr21$counts_normed

## ----fig.width=6,fig.height=5,message=FALSE,eval=TRUE-------------------------
##gene named "SOD1"
plot_coverage(gene_name="SOD1", coverage_output=coverage_res_chr21, 
            degnorm_output=res_DegNorm_chr21,group=c(0,1,1))

## ----fig.width=6,fig.height=5,message=FALSE,warning=FALSE,eval=TRUE-----------
plot_boxplot(DI=res_DegNorm_chr21$DI)

## ----fig.width=6,fig.height=5,message=FALSE,eval=TRUE-------------------------
plot_heatmap(DI=res_DegNorm_chr21$DI)

## ----fig.width=6,fig.height=5,message=FALSE,warning=FALSE,eval=TRUE-----------
plot_corr(DI=res_DegNorm_chr21$DI)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

