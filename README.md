DegNorm R package
================

**Maintainer**: Ji-Ping Wang, \<<jzwang@northwestern.edu>\>

([**Python package
download**](https://nustatbioinfo.github.io/DegNorm/))\*\*

**Reference**: Xiong, B., Yang, Y., Fineis, F. Wang, J.-P., DegNorm:
normalization of generalized transcript degradation improves accuracy in
RNA-seq analysis, Genome Biology, 2019,20:75

## What is DegNorm?

DegNorm, short for **Deg**radation **Norm**alization, is a
bioinformatics pipeline designed to correct for bias due to the
heterogeneous patterns of transcript degradation in RNA-seq data.
DegNorm helps improve the accuracy of the differential expression
analysis by accounting for this degradation.

In practice, RNA samples are often more-or-less degraded, and the
degradation severity is not only sample-specific, but gene-specific as
well. It is known that longer genes tend to degrade faster than shorter
ones. As such, commonplace global degradation normalization approaches
that impose a single normalization factor on all genes within a sample
can be ineffective in correcting for RNA degradation bias.

## DegNorm pipline available formats

We’ve developed an R package and an indepedent Python package
([download](https://nustatbioinfo.github.io/DegNorm/)), both of which
allow to run the entire pipeline from the RNA-seq alignment (.bam)
files.

## DegNorm main features

DegNorm R package contains two major functions: (1) processing the
RNA-seq alignment file (.bam) to calculate the coverage; and (2) using a
core algorithm written in RcppArmadillo to perform rank-one
over-approximation on converage matrices for each gene to estimate the
degramation index (DI) score for each gene within each sample.

DegNorm outputs DI scores together with degradation-normalized read
counts (based on DI scores). It also provides supplementary functions
for visualization of degradation at both gene and sample level. The
following diagram illustrates the flow of DegNorm pipeline.

 

![A diagram of
DegNorm.](http://bioinfo.stats.northwestern.edu/jzwang/DegNorm/degnorm_logo.png)

 

The following vignette is intended to provide example codes for running
DegNorm R package. It presumes that you have successfully installed
DegNorm package. We illustrate below how to: 1) calculate the read
coverage curves for all genes within all samples, and 2) perform
degradation normalization on coverage curves. Either step is computing
intensive. Dependent upon the number of samples and the sequencing
depth, the total computing time may last a few hours. DegNorm utilizes
the parallel computing functionality of R and automatically detects the
number of cores on your computer to run jobs in parallel. Due to the
large size of bam file and limited computing power of personal computer,
we recommend users to run it in servers or computing
clusters.

### 1\. Compute coverage score based on alignment .bam files

##### Set up input file: .bam and .gtf files.

``` r
## specify bam_files from RNA-seq, you should replace it by your own bam files
bam_file_list=list.files(path=system.file("extdata",package="DegNorm"),
                        pattern=".bam$",full.names=TRUE)
```

The three bam files were subsetted from a specific region of chorosome
21 from the origianl bam for package size limitation. Original files can
be found from the included reference
above.

``` r
## gtf_file you used for RNA-seq alignment, replace it by your own gtf file
gtf_file=list.files(path=system.file("extdata",package="DegNorm"),
                    pattern=".gtf$",full.names=TRUE)
```

##### Run main function to create read coverage matrix and read counts

``` r
## calculate the read coverage score for all genes of all samples
coverage_res_chr21_sub=read_coverage_batch(bam_file_list, gtf_file,cores=2)
```

    ## Parse gtf file...done 
    ## Index SRR873822_chr21_9900000-10000000.bam 
    ## Computing coverage for SRR873822_chr21_9900000-10000000.bam ... 
    ##    SRR873822_chr21_9900000-10000000.bam is a paired-end bam file 
    ## SRR873822_chr21_9900000-10000000.bam is done! 
    ## 
    ## Index SRR873834_chr21_9900000-10000000.bam 
    ## Computing coverage for SRR873834_chr21_9900000-10000000.bam ... 
    ##    SRR873834_chr21_9900000-10000000.bam is a paired-end bam file 
    ## SRR873834_chr21_9900000-10000000.bam is done! 
    ## 
    ## Index SRR873838_chr21_9900000-10000000.bam 
    ## Computing coverage for SRR873838_chr21_9900000-10000000.bam ... 
    ##    SRR873838_chr21_9900000-10000000.bam is a paired-end bam file 
    ## SRR873838_chr21_9900000-10000000.bam is done!

`cores` argument specifies the number of cores to use. Users should try
to use as many as possible cores to maximize the computing efficiency.

``` r
## save the coverage results
save(coverage_res_chr21_sub,file="coverage_res_chr21_sub.Rda")
```

Function `read_coverage_batch` returns the coverage matrices as a list,
one per gene, and a dataframe for read counts, each row for one gene and
each column for one sample.

``` r
data("coverage_res_chr21")
## summarize the coverage results
summary_CoverageClass(coverage_res_chr21)
```

    ## CoverageClass from read_coverage_batch function 
    ## $coverage    : list of length 339 
    ## $counts  : data.frame of dimension 339 by 3 
    ## 
    ## Samples:      SRR873822_chr21.bam SRR873834_chr21.bam SRR873838_chr21.bam 
    ## Total number genes:   339

``` r
## extract coverage scores and counts from coverage_res
coverage_matrix=coverage_res_chr21$coverage
counts=coverage_res_chr21$counts
```

### 2\. DegNorm core algorithm

Run degnorm core algorithm for degradation normalization. DegNorm
purpose is for differential expression analysis. Thus genes with
extremely low read counts from all samples are filtered out. The current
filtering criterion is that if more than half of the samples have less
than 5 read count, that gene will not be considered in the degnorm
algorithm. In the following example, I am using downsamling to save time
below (default). Alternatively you can set down\_sampling = 0, which
takes longer time. If `down_samplin= 1`, read coverage scores are binned
with size by `grid_size` for baseline selection to achieve better
efficiency. The default `grid_size` is 10 bp. We recommend to use a
`grid_size` less than 50 bp. `iteration` specifies the big loop in
DegNorm algorithm and 5 is usually sufficient. `loop` specifies the
iteration number in the matrix factorization over-approximation.

``` r
res_DegNorm_chr21 = degnorm(read_coverage = coverage_res_chr21[[1]],
                    counts = coverage_res_chr21[[2]],
                    iteration = 5,
                    down_sampling = 1,
                    grid_size=10,
                    loop = 100,
                    cores=2)
```

If `down_sampling= 0`, then the argument `grid_size` is ignored.

``` r
## save the DegNorm results
save(res_DegNorm_chr21,file="res_DegNorm_chr21.Rda")
```

Function `degnorm` returns a list of multiple objects. counts\_normed is
the one with degradation normalized read counts for you to input DeSeq
or EdgeR for DE analysis.

``` r
data("res_DegNorm_chr21")
```

``` r
## summary of the DegNorm output
summary_DegNormClass(res_DegNorm_chr21)
```

    ## DegNormClass from DegNorm function: 
    ## $counts  : data.frame of dimension 132 by 3 
    ## $counts_normed   : data.frame of dimension 132 by 3 
    ## $DI      : matrix array of dimension 132 by 3 
    ## $K       : matrix array of dimension 132 by 3 
    ## $convergence     : integer of length 132 
    ## $envelop     : list of length 132 
    ## 
    ## Samples:      SRR873822_chr21.bam SRR873834_chr21.bam SRR873838_chr21.bam 
    ## Total number genes:   132

The difference of number of genes between `res_DegNorm` and
`coverage_res` is 207 (339-132). The 207 genes were filtered out from
`degnorm` degradation normalization because less than half of the
samples (3) have more than 5 read count.

``` r
## extrac normalized read counts
counts_normed=res_DegNorm_chr21$counts_normed
```

### 3\. Plot functions in DegNorm

DegNorm provides four plot functions for visualization of degradation
and sample quality diagnosis.

  - plot\_coverage
  - plot\_corr
  - plot\_heatmap
  - plot\_boxplot

##### – Plot the before-/after-degradation coverage curves

``` r
##gene named "SOD1"
plot_coverage(gene_name="SOD1", coverage_output=coverage_res_chr21, 
            degnorm_output=res_DegNorm_chr21,group=c(0,0,1))
```

##### – Boxplot of the degradation index(DI) scores

``` r
plot_boxplot(DI=res_DegNorm_chr21$DI)
```

##### – Heatmap plot of the degradation index(DI) scores

``` r
plot_heatmap(DI=res_DegNorm_chr21$DI)
```

#### – Correlation matrix plot of degradation index(DI) scores

``` r
plot_corr(DI=res_DegNorm_chr21$DI)
```
