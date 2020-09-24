################################################################################
# summary function for CoverageClass and DegNormClass
#   - summary_CoverageClass
#   - summary_DegNormClass
# Plot functions:
#   - plot_corr(DI)
#   - plot_boxplot(DI)
#   - plot_heatmap
#   - plot_coverage(coverage)
################################################################################
summary_CoverageClass=function(object){
    cat("CoverageClass from read_coverage_batch function","\n")
    cat(paste("$",names(object)[1],sep=""),"\t:", class(object$coverage), 
        "of length", length(object$coverage),"\n")
    cat(paste("$",names(object)[2],sep=""),"\t:", class(object$counts), 
        "of dimension", dim(object$counts)[1],"by",dim(object$counts)[2],"\n")
    cat("\nSamples:\t\t",colnames(object[[2]]),"\n")
    cat("Total number genes:\t", nrow(object[[2]]),"\n")
}

summary_DegNormClass=function(object){
    cat("DegNormClass from DegNorm function:","\n")
    cat(paste("$",names(object)[1],sep=""),"\t:", class(object$counts), 
        "of dimension", dim(object$counts)[1],"by",dim(object$counts)[2],"\n")
    cat(paste("$",names(object)[2],sep=""),"\t:", class(object$counts_normed), 
        "of dimension", dim(object$counts_normed)[1],"by",
        dim(object$counts_normed)[2],"\n")
    cat(paste("$",names(object)[3],sep=""),"\t\t:", class(object$DI), 
        "of dimension", dim(object$DI)[1],"by",dim(object$DI)[2],"\n")
    cat(paste("$",names(object)[4],sep=""),"\t\t:", class(object$K), 
        "of dimension", dim(object$K)[1],"by",dim(object$K)[2],"\n")
    cat(paste("$",names(object)[5],sep=""),"\t:", class(object$convergence), 
        "of length", length(object$convergence),"\n")
    cat(paste("$",names(object)[6],sep=""),"\t:", class(object$envelop), 
        "of length", length(object$envelop),"\n")
    cat("\nSamples:\t\t",colnames(object[[2]]),"\n")
    cat("Total number genes:\t", nrow(object[[2]]),"\n")
}

################################################################################
# plot functions for heatmap, boxplot and coverage curves
################################################################################
################################################################################
# plot.coverage: function to plot before/after degradation normalization 
#               coverage curve
################################################################################
plot_coverage <- function(gene_name, coverage_output, degnorm_output, 
                        group=NULL){
    # coverage plot before and after normalization for one single gene
    # input:
    #    - gene_name: the name of the gene
    #    - coverage_res: the output from read_coverage_batch function
    #    - degnorm_output: the output from DegNorm function
    # output:
    #    - a ggplot2 object of coverage plot before/after DegNorm normalization
    #    - group: a vector of factors representing sample conditions in the 
    #      order of the coverage matrix columns, either integer or char.
    # retrieve raw coverage from read_coverage_batch output coverage_res
    raw.coverage.matrix = data.table(t(coverage_output$coverage[[gene_name]]))
    n.sample = dim(raw.coverage.matrix)[2]
    # retrieve abundance scale K and envelope function from DegNorm output 
    K=as.matrix(degnorm_output$K[which(row.names(degnorm_output$K)
                                        ==gene_name),])
    envelop=degnorm_output$envelop[[gene_name]]
    degnorm.coverage.matrix = data.table(t(K%*%envelop))
    n=nrow(raw.coverage.matrix)
    raw.coverage.matrix=cbind(position=rep(seq_len(n),n.sample),label=
                            rep("raw",n.sample*n),stack(raw.coverage.matrix))
    degnorm.coverage.matrix=cbind(position=rep(seq_len(n),n.sample),
                                label=rep("DegNorm",n.sample*n),
                                stack(degnorm.coverage.matrix))
    dat.curve = data.table(rbind(raw.coverage.matrix, degnorm.coverage.matrix))
    colnames(dat.curve)[3:4] = c("coverage", "sample")
    default.palette=c("#0072B2","#D55E00","#009E73","#999999", "#E69F00", 
                    "#56B4E9", "#F0E442",  "#CC79A7")
    if(is.null(group)){ #customize color palette
        if(n.sample <= length(default.palette)){
            custom_color = default.palette[seq_len(n.sample)]
        }else{
            custom_color = plasma(n.sample)
        }
    }else{
        n.group = length(unique(group))
        if(length(n.group) <= length(default.palette)){
            custom_color = rep(default.palette[seq_len(n.group)], n.group)
        }else{
            custom_color=rep(plasma(length(n.group), end=0.9),n.group)
        }
    }
    p = ggplot(data = dat.curve, aes(x = dat.curve$position,
                            y = dat.curve$coverage,col = dat.curve$sample))+
        xlab("Transcript position") + ylab("Coverage score") +
        scale_color_manual(values = custom_color) +
        geom_line(size = 0.8) + theme_light() + facet_grid(dat.curve$label~.)
    return(p)
}

################################################################################
# generate boxplots from DI scores
################################################################################

plot_boxplot<-function(DI){
    # input:
    #     DI: a data table with each row  gene, each column for sample
    # output file:
    #     a ggplot2 object for boxplot
    
    dat_DI = data.frame(stack(DI))
    colnames(dat_DI)=c("gene_names","sample","DI_score")
    
    p.boxplot = 
        ggplot(data = dat_DI, aes(x = dat_DI$sample, y = dat_DI$DI_score)) + 
        xlab("Sample") + ylab("DI score") +
        stat_boxplot(geom ='errorbar') +
        geom_boxplot(aes(x = dat_DI$sample, y = dat_DI$DI_score, 
                        col = dat_DI$sample), outlier.shape = NA) +
        theme_classic() + ylim(c(0,max(DI))) +
        theme(legend.position="none")
    return(p.boxplot)
}

################################################################################
# generate heatmaps for the correlation matrix from DI scores
plot_corr<-function(DI){
    # input:
    #     DI: a data table with each row  gene, each column for sample
    # output file:
    #     a ggplot2 object for correlation matrix heatmap 

    p = heatmaply(apply(cor(DI),2,rev), 
        dendrogram = 'none',  # without dendrogram
        margins = c(100,100),
        width = 1200, height = 1000)
return(p)
}

################################################################################
# generate heatmaps of the DI matrix
plot_heatmap<-function(DI){
    # input:
    #     DI: a data table with each row  gene, each column for sample
    # output file:
    #     a ggplot2 object for correlation matrix heatmap 

    p1 = plot_ly(x=colnames(DI),z = as.matrix(DI[order(rowSums(DI)),]), 
            #colorscale = colorRamp(c("blue", "red")), 
            type = "heatmap", showscale = TRUE)
    # sort each column
    DI = apply(DI,2,sort)
    p2 = plot_ly(x=colnames(DI),z = as.matrix(DI[order(rowSums(DI)),]), 
            #colorscale = colorRamp(c("blue", "red")), 
            type = "heatmap", showscale = FALSE)
    # join two figures
    p = subplot(p1, p2, nrows = 1, shareY = TRUE)
return(p)
}

