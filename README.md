
#  DegNorm R package
Authors: Bin Xiong, Ji-Ping Wang
   
Maintainer: [Ji-Ping Wang](http://bioinfo.stats.northwestern.edu/~jzwang/) <<jzwang@northwestern.edu>>
   
**Python package** ([download](https://nustatbioinfo.github.io/DegNorm/))

**Reference**:
Xiong, B., Yang, Y., Fineis, F. Wang, J.-P., DegNorm: normalization of generalized transcript
degradation improves accuracy in RNA-seq analysis, Genome Biology, 2019,20:75

## What is DegNorm?

DegNorm, short for **Deg**radation **Norm**alization, is a bioinformatics pipeline designed to correct 
for bias due to the heterogeneous patterns of transcript degradation in RNA-seq data. DegNorm 
helps improve the accuracy of the differential expression analysis by accounting for this degradation. 

In practice, RNA samples are often more-or-less degraded, and the degradation severity is not only 
sample-specific, but gene-specific as well. It is known that longer genes tend to degrade faster 
than shorter ones. As such, commonplace global degradation normalization approaches that impose a 
single normalization factor on all genes within a sample can be ineffective in correcting for RNA degradation bias.


## DegNorm pipline available formats

We've developed an R package and an indepedent Python package ([download](https://nustatbioinfo.github.io/DegNorm/)), 
both of which allow to run the entire pipeline from the RNA-seq alignment (.bam) files.


## DegNorm diagram

![A diagram of DegNorm.](
http://bioinfo.stats.northwestern.edu/jzwang/DegNorm/degnorm_logo.png)


## DegNorm main features

DegNorm R package contains two major functions: (1) processing the RNA-seq alignment file (.bam) to calculate the coverage; 
and (2) using a core algorithm written in RcppArmadillo to perform rank-one over-approximation on converage matrices for each 
gene to estimate the degramation index (DI) score for each gene within each sample.

DegNorm outputs DI scores together with degradation-normalized read counts (based on DI scores). 
It also provides supplementary functions for visualization of degradation at both gene and sample level.

The following vignette is intended to provide example codes for running  DegNorm R package. 
It presumes that you have successfully installed DegNorm package. We illustrate below how to: 1) calculate the read coverage curves for all genes within  all samples, and 2) perform degradation normalization on coverage curves. 
Either step is computing intensive. Dependent upon the number of samples and the 
sequencing depth, the total computing time may last a few hours. DegNorm utilizes 
the parallel computing functionality of R and automatically detects the number of 
cores on your computer to run jobs in parallel. Due to the large size of bam file and 
limited computing power of personal computer, we recommend users to run it in servers or computing clusters.


## Examples

### 1. Compute coverage score based on alignment .bam files

##### Set up input file: .bam and .gtf files.
```{r,eval=FALSE}
## set up working directory
library(DegNorm)
setwd("~/degnorm") #working directory
```

```{r}
## specify bam_files from RNA-seq, you should replace it by your own bam files
bam_file_list=list.files(path=system.file("extdata",package="DegNorm"),pattern=".bam$",
                          full.names=TRUE)
```

```{r}
## gtf_file you used for RNA-seq alignment, you should replace it by your own gtf file
gtf_file=list.files(path=system.file("extdata",package="DegNorm"),pattern=".gtf$",
                    full.names=TRUE)
```

##### Run main function to create read coverage matrix and read counts

```{r,eval=FALSE}
## calculate the read coverage score for all genes of all samples
coverage_res=read_coverage_batch(bam_file_list, gtf_file,cores=2)
```
`cores` argument specifies the number of cores to use. Users should try 
to use as many as possible cores to maximize the computing efficiency.


```{r,eval=FALSE}
## save the coverage results
save(coverage_res,file="coverage_res_chr21.Rda")
```

Function `read_coverage_batch` returns the coverage matrices as a list, one per gene, and a dataframe for read counts, each row for one gene and each column for one sample.

```{r}
load("coverage_res_chr21.Rda")
## summarize the coverage results
summary(coverage_res)
```

```{r}
## extract coverage scores and counts from coverage_res
coverage_matrix=coverage_res$coverage
counts=coverage_res$counts
```
### 2. DegNorm core algorithm

Run degnorm core algorithm for degradation normalization. DegNorm purpose is for differential expression analysis. Thus genes with extremely low read counts from all samples are filtered out. The current filtering criterion is that if more than half of the samples have less than 5 read count, that gene will not be considered in the degnorm algorithm. In the following example, I am using downsamling to save time below (default). Alternatively you can set down_sampling = 0, which takes longer time.
If `down_samplin= 1`, read coverage scores are binned with size by `grid_size` for baseline selection to achieve better efficiency. The default `grid_size` is 10 bp. We recommend to use a `grid_size` less than 50 bp. `iteration` specifies the big loop in DegNorm algorithm and 5 is usually sufficient. `loop` specifies the iteration number in the matrix factorization over-approximation.

```{r,eval=FALSE}
res_DegNorm = degnorm(read_coverage = coverage_res[[1]],
                      counts = coverage_res[[2]],
                      iteration = 5,
                      down_sampling = 1,
                      grid_size=10,
                      loop = 100)
```
If `down_samplin= 0`, then the argument `grid_size` is ignored.


```{r,eval=FALSE}
## save the DegNorm results
save(res_DegNorm,file="res_DegNorm_chr21.Rda")
```

Function `degnorm` returns a list of multiple objects. counts_normed is the one with degradation
normalized read counts for you to input DeSeq or EdgeR for DE analysis.

```{r}
load("res_DegNorm_chr21.Rda")
```

```{r}
## summary of the DegNorm output
summary(res_DegNorm)
```

The difference of number of genes between `res_DegNorm` and `coverage_res` is 207 (339-132). The 207 genes were filtered out from `degnorm` degradation normalization because less than half of the samples (3) have more than 5 read count.


```{r}
## extrac normalized read counts
counts_normed=res_DegNorm$counts_normed
```

### 3. Plot functions in DegNorm

DegNorm provides four plot functions for visualization of degradation and sample quality diagnosis.

 * plot_coverage
 * plot_corr
 * plot_heatmap
 * plot_boxplot

##### -- Plot the before-/after-degradation coverage curves 

```{r,fig.width=6,fig.height=5,message=FALSE}
## gene named "SOD1"
plot_coverage(gene_name="SOD1", coverage_output=coverage_res, degnorm_output=res_DegNorm)
```

##### -- Boxplot of the degradation index(DI) scores

```{r,fig.width=6,fig.height=5,message=FALSE,warning=FALSE}
plot_boxplot(DI=res_DegNorm$DI)
```

##### -- Heatmap plot of the degradation index(DI) scores

```{r,fig.width=6,fig.height=5,message=FALSE}
plot_heatmap(DI=res_DegNorm$DI)
```

##### -- Correlation matrix plot of degradation index(DI) scores

```{r,fig.width=6,fig.height=5,message=FALSE,warning=FALSE}
plot_corr(DI=res_DegNorm$DI)
```
