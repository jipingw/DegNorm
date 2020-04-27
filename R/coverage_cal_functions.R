
###############################################################################
############# Coverage calculation for batch files ############################
## read_coverage_batch: function to calculate read coverage score in batches
##    input:
##       - bam_file_list: list of bam file names
##       - gtf_file: name of gtf file where the RNA-seq data were aligned with 
##                 reference to.
##    output:
##       - coverageClass object, which includes:
##         (1) list of coverage matrices 
##         (2) a dataframe for read counts          
##
## Notes: 
##   1. Coverage score is calcualted per gene, i.e. concatenation of all
##      exons from the same gene.
##   2. We follow HTseq for counting valid read/read pairs for each gene.
##   3. When reading alignment file, isSecondaryAlignment flag is set as 
##      FALSE to avoid possible redundant counting.
##   4. For paired-end data, isPaired is set as TRUE. We don't recommend
##      setting isProperPair as TRUE as some fragments length may exceed 200bp.
##   5. User can modify scanBamParam in the R codes below as needed.
###############################################################################
read_coverage_batch=function(bam_file_list,gtf_file){
    i=1; options(warn=-1)
    all_genes=.gtf_parse(gtf_file)
    ##simplify bam file names if path information present in file names
    sample_names=as.vector(vapply(bam_file_list, function(x) 
                        rev(strsplit(x,"[\\\\]|[/]")[[1]])[1],character(1)))
    for (k in seq_along(bam_file_list)){
        #if bam index not exist, generate index
        if(!file.exists(paste(sample_names[k],".bai",sep=""))){
            cat("Index", sample_names[k],"\n")
            indexBam(bam_file_list[k])
        }
        cat ("Computing coverage for", sample_names[k],"...", "\n")
        results=read_coverage(bam_file_list[k], all_genes)
        ##separate coverage and read_counts
        coverage_RLE=RleList()    ## coverage in Rle format
        read_counts=data.frame()  ## data frame for read-counts
        read_counts=do.call(rbind,lapply(results,`[[`,2))
        coverage_RLE=lapply(results,`[[`,1)
        ##convert Rle object into vector coverage
        cores <- detectCores(logical = FALSE)
        cl <- makeCluster(cores-1)
        registerDoParallel(cl)
        coverage_RLE=coverage_RLE[lengths(coverage_RLE)>0]
        coverage_vector=foreach(i = seq_along(coverage_RLE),
                .packages=c("Rsamtools","GenomicAlignments","GenomicFeatures"),
                            .inorder=TRUE,.multicombine=TRUE) %dopar%
                    {  res=lapply(coverage_RLE[[i]],as.vector)
                        return(res)}
        stopImplicitCluster()
        stopCluster(cl)
        coverage_vector=do.call(c,coverage_vector)
        cat (sample_names[k], "is done!", "\n\n")
        if (k==1){
                coverage_matrix=coverage_vector
                read_counts_matrix=read_counts
        }else {
                coverage_matrix=mapply(rbind,coverage_matrix,coverage_vector)
                read_counts_matrix=cbind(read_counts_matrix,read_counts)
        }   
        coverage_results=list(coverage_matrix,read_counts_matrix)
    }
    colnames(read_counts_matrix)=sample_names
    coverage_matrix=lapply(coverage_matrix,function(x) 
                            {row.names(x)=sample_names; x})
    output=list(coverage=coverage_matrix,counts=read_counts_matrix)
    class(output)="CoverageClass"
    return(output)
}

################################gtf file parsing ###############################
# .gtf_parse: function to parse gtf file by compiling all exons from the same
#            gene to generate a total transcript. Genes which mapped to multiple
#             chromosomes or multiple strands are excluded to avoid confusion.
##    input:
##       - gtf_file: name of gtf file where the RNA-seq data were aligned with 
##                  reference to.
##    ouput:
##       - GRangesList object
################################################################################
.gtf_parse<-function(gtf_file){
    ## parse gtf file to create GRangesList object
    cat("Parse gtf file...") ## parse gtf file to create GRangesList object
    suppressMessages(txdb <- makeTxDbFromGFF(gtf_file))
    all_genes <- exonsBy(txdb, by="gene")
    ##exclude genes on two strands or genes exist on multiple chromosomes
    all_genes=all_genes[vapply(unique(runValue(strand(all_genes))),
                            length,numeric(1))==1]
    all_genes=all_genes[vapply(unique(runValue(seqnames(all_genes))),
                            length,numeric(1))==1]
    ##sorting exons on - strand into ascending rank
    all_genes=endoapply(all_genes, function(x) if(runValue(strand(x))=="-") 
                {return(rev(x))}else{return(x)})
    cat("done", "\n")
    return(all_genes)
}



################################################################################
######### Converage calculation for one bam file ###############################
## read_coverage: function to calculate read coverage score for one bam file
##                This function determines whether it single-end and paired-end,
##                and generates bam file index if needed.It takes multiple-core 
##                structure for parallel computing for efficiency.
################################################################################
##    input:
##       - bam_file: bam file names
##       - all_genes: parsed genes.gtf file
##
##    output:
##       - coverageClass object, which includes:
##         (1) list of coverage score for each gene in RLE format and 
##         (2) a dataframe for read counts          
##
################################################################################


read_coverage=function(bam_file,all_genes){
    cores <- detectCores(logical = FALSE)
    cl <- makeCluster(cores-1)
    registerDoParallel(cl)

    ##Filter out empty chromosomes (i.e. no genes  on chromosomes)
    ##chrom=seqnames(seqinfo(all_genes))
    chrom=as.vector(unique(unlist(runValue(seqnames(all_genes)))))

    ##for memory issue, coverage is calculated by chromosomes

    sample_name=rev(strsplit(bam_file,"[\\\\]|[/]")[[1]])[1]
    i=1
    if(suppressMessages(testPairedEndBam(bam_file))==TRUE){
        cat("  ",sample_name, "is a paired-end bam file","\n")
        results <- foreach(i = seq_along(chrom),
                        .export=c(".paired_end_cov_by_ch",
                        ".single_end_cov_by_ch", ".IntersectionStrict2"),
                        .packages = c("Rsamtools", "GenomicAlignments",
                                    "GenomicFeatures","data.table"),
                        .inorder=TRUE,.multicombine=TRUE) %dopar%
                        { 
                        res=.paired_end_cov_by_ch(bam_file,all_genes,chrom[i])
                        return(res)
                        }
    }else{
        cat("  ",sample_name, "is a single-end bam file","\n")
        results <- foreach(i = seq_along(chrom),
                        .export=c(".paired_end_cov_by_ch",
                                ".single_end_cov_by_ch",".IntersectionStrict2"),
                        .packages = c("Rsamtools", "GenomicAlignments",
                                "GenomicFeatures"),
                        .inorder=TRUE,.multicombine=TRUE) %dopar%
                        { 
                        res=.single_end_cov_by_ch(bam_file,all_genes,chrom[i])
                        return(res)}
    }
    stopImplicitCluster()
    stopCluster(cl)
    return(results)
}

################################################################################
###### Coverage calculation for paired end for one chromosome###################
## .paired_end_cov_by_ch: function to calculate read coverage score for one
##  chromosome for paired-end data.
################################################################################
##    input:
##       - bam: bam file names
##       - all_genes: parsed genes.gtf file
##       - ch: chromosome name
##
##    output:
##       - coverageClass object, which includes:
##         (1) list of coverage score for each gene in RLE format and 
##         (2) a dataframe for read counts          
##
## Notes: 
##   1. Coverage score is calcualted per gene, i.e. concatenation of all
##      exons from the same gene.
##   2. We follow HTseq protocol to count valid reads/read pairs for each gene.
##   3. When reading alignment file, isSecondaryAlignment flag is set as 
##      FALSE to avoid possible redundant counting.
##   4. isPaired is set as TRUE. We don't recommend setting isProperPair 
##      as TRUE as some fragments length may exceed 200bp.
##   5. User can modify scanBamParam in the R codes below as needed.
##############################################################################

.paired_end_cov_by_ch=function(bam,all_genes,ch){
    ##function to calculate the coverage score and return in Rle objects
    bf = BamFile(bam, asMates = TRUE, qnameSuffixStart = ".")
    tx_compat_cvg=RleList()
    genes=all_genes[unlist(runValue(seqnames(all_genes)))==ch]
    gr = as(seqinfo(bf), "GRanges")
    if(ch %in% runValue(seqnames(gr))){
        param = ScanBamParam(which = gr[ch],flag=scanBamFlag(isPaired=TRUE, 
                            isSecondaryAlignment=FALSE))
        gal2=readGAlignmentPairs(bf, param = param)
        gal2=gal2[strand(gal2)!="*"]
        gal2=gal2[is.na(seqnames(gal2))==FALSE]
        grl <- as(gal2, "GRangesList")
        gal3=reduce(grl)
        #only keep genes on the current chromosome
        results=.IntersectionStrict2(genes,gal2)
        tx2reads <- setNames(as(t(results), "List"), names(genes))
        ## calculate read pair
        read_counts=data.frame(count=lengths(tx2reads))
        rownames(read_counts)=names(tx2reads)

        #relist gal3 based on compatibility
        temp0=data.table(group=as(gal3,"data.frame")$group,
                        index=seq_len(sum(lengths(gal3))))
        #group2indexMap=unique(setDT(temp0))[, .(values = index), by =group ] 

        idmatch=function(tx2reads,group2indexMap){
            return(group2indexMap[group2indexMap$group %in% tx2reads]$index)
        }
        #tx2reads1=lapply(tx2reads,idmatch,group2indexMap=group2indexMap)
        tx2reads1=lapply(tx2reads,idmatch,group2indexMap=temp0)
        gal4=unlist(gal3)
        rm(gal3)
        compat_reads_by_tx <- extractList(gal4, tx2reads1)
        tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx,
                                            genes,
                                            ignore.strand=TRUE)
    }else{
        tx2reads1=IntegerList(vector("list", length(genes)))
        left_gal=GAlignments() #empty alignment file
        names(tx2reads1)=names(genes)
        read_counts=data.frame(count=lengths(tx2reads1))
        rownames(read_counts)=names(tx2reads1)
        compat_reads_by_tx <- extractList(left_gal, tx2reads1)
        tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx,
                                            genes,
                                            ignore.strand=FALSE)
    }
    return(list(tx_compat_cvg,read_counts))
}

################################################################################
###### Coverage calculation for single end for one chromosome###################
## .single_end_cov_by_ch: function to calculate read coverage score for
##                       one chromosome for sigle-end data
################################################################################
##    input:
##       - bam: bam file names
##       - all_genes: parsed genes.gtf file
##       - ch: chromosome name
##
##    output:
##       - coverageClass object, which includes:
##         (1) list of coverage score for each gene in RLE format and 
##         (2) a dataframe for read counts          
##
## Notes: 
##   1. Coverage score is calcualted per gene, i.e. concatenation of all
##      exons from the same gene.
##   2. We follow HTseq protocol to count valid reads/read pairs for each gene.
##   3. When reading alignment file, isSecondaryAlignment flag is set as 
##      FALSE to avoid possible redundant counting.
##   4. User can modify scanBamParam in the R codes below as needed.
################################################################################
.single_end_cov_by_ch=function(bam,all_genes,ch){
    ##function to calculate the coverage score and return in Rle objects
    bf = BamFile(bam)
    tx_compat_cvg=RleList()
    read_counts=data.frame()

    genes=all_genes[unlist(runValue(seqnames(all_genes)))==ch]
    gr = as(seqinfo(bf), "GRanges")
    if(ch %in% runValue(seqnames(gr))){
        param = ScanBamParam(which = gr[ch],
                            flag=scanBamFlag(isSecondaryAlignment=FALSE))
        gal2=readGAlignments(bf, param = param)
        gal2=gal2[strand(gal2)!="*"]
        gal2=gal2[is.na(seqnames(gal2))==FALSE]

        #only keep genes on the current chromosome
        results=.IntersectionStrict2(genes,gal2)
        tx2reads <- setNames(as(t(results), "List"), names(genes))
        compat_reads_by_tx <- extractList(gal2, tx2reads)

        read_counts=data.frame(count=lengths(tx2reads))
        rownames(read_counts)=names(tx2reads)
        tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx,
                                            genes,
                                            ignore.strand=TRUE)
    }else{
        tx2reads1=IntegerList(vector("list", length(genes)))
        names(tx2reads1)=names(genes)
        left_gal=GAlignments()
        compat_reads_by_tx <- extractList(left_gal, tx2reads1)
        read_counts=data.frame(count=lengths(tx2reads1))
        colnames(read_counts)=names(tx2reads1)  

        tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx,
                                            genes,
                                            ignore.strand=FALSE)
    }
    return(list(tx_compat_cvg,read_counts))
}

################################################################################
###########Function to compile compatible reads as HtSeq for each gene####
## .IntersectionStrict2: function modified from GenomicAlignment package for
##                      selection of compatible read/read pairs for each genes 
##                      following HTseq.
################################################################################
## input: please refer to IntersectionStrict function from GenomicAlignments
##        In particular, ignore.strand should be set as TRUE for paired-end.
##        For single-end data, depending on the sequencing, users should
##        set it accordingly.
################################################################################
## output: 
##        reads IDs compatible with features.
################################################################################
.IntersectionStrict2= function (features, reads,
                                ignore.strand = TRUE, inter.feature = TRUE) {
    ##function adapted to compile compatible reads as HTseq

    ov <- findOverlaps(reads, features, type = "within", 
                    ignore.strand = ignore.strand)
    if (inter.feature) {
        reads_to_keep <- which(countQueryHits(ov) == 1L)
        ov <- ov[queryHits(ov) %in% reads_to_keep]
    }
    ov
}



