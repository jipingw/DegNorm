
################################################################################
# ratioSVD: function to calulate SVD for the raw coverage
#           used to identify genes with similar coverage pattern
#           across all samples (minimal degradation)
# input:    coverage matrix (f) with rows representing samples
# output:   DI scores (rho)
################################################################################

.ratioSVD <- function(f){
    model.svd = svd(f, nv = 1, nu = 1)
    fitted = model.svd$d[1] * model.svd$u %*% t(model.svd$v)
    fitted[fitted < f] = f[fitted < f]
    return(1-rowSums(f)/(rowSums(fitted)+1))
}

################################################################################
# degnorm: function to perform degradation normalization
# input:
#    - coverage: list of gene coverages, each row stands for one sample;
#    - counts: read count matrix with gene ordered same as the coverage list;
#    - iterations: number of iterations for DegNorm (default to be 5);
#    - loop: number of iterations for NMF-OA (default to be 100).
#    - down_sampling: 1 for yes and 0 or no. If yes, coverage score is averaged
#                     on an interval of size specificed by grid_size. The
#                   size-reduced coverage matrix is fed for baseline selection
#                     for efficiency consideration
#    - grid_size: the bin size for average coverage score
#    - cores: number of cores to use. Default is  maximum number cores -1
# output: a list containing the following
#    - DI: Degradation Index matrix;
#    - counts_normed: DegNorm normalized counts;
#    - convergence: convergence status indicating whether baseline is found.
#    - envelop: envelop curve which stands for the estimated commons
#               shape before degradation.
################################################################################
degnorm <- function(read_coverage, counts, iteration=5, loop=100,
                    down_sampling=1, grid_size=10,cores=1){
    # filtering out genes with low expressions
    if(length(read_coverage)!=nrow(counts)) stop("Error: number of genes from
                    coverage list is inconsistent with the counts matrix!")
    if((down_sampling %in% c(0,1))==FALSE) stop("Error: down_sampling should be
                    1 (down-sampling) or 0 (not down-sampling)!")
    if(grid_size %% 1!=0||grid_size<0) {
        stop("Error: grid_size must be a positive integer.")}
    # filtering out genes with low expressions
    message("Filtering out genes with low read counts (x<5)..")
    n = dim(counts)[2]
    filter = apply(counts, 1, function(x) length(x[x>5])>=floor(n/2))
    counts = counts[filter,];  m = dim(counts)[1]
    read_coverage = read_coverage[which(filter=="TRUE")]
    # Estimating library sizes / scaling factors from SVD
    message("Initial estimation of count normalization factors...")
    norm.factor = colSums(counts)/median(colSums(counts))
    # Register multiple nodes
    max.cores=detectCores(logical = FALSE)
    if(max.cores-1 < cores) cores <- max.cores-1
    cl <- makeCluster(cores);  registerDoParallel(cl)
    message("Initial degradation normalization without base-line selection...")
    rho = NULL
    i=1
    rho = foreach(i = seq_len(m),.combine = "rbind",
                    .export = ".ratioSVD") %dopar% {
                    .ratioSVD(f = read_coverage[[i]])}
    stopImplicitCluster()
    stopCluster(cl)
    norm.factor = colSums(counts[apply(rho,1,max) < 0.10,])
    norm.factor = norm.factor/stats::median(norm.factor)
    scale = norm.factor
    counts.weighted = sweep(counts, 2, norm.factor, FUN = "/")
    message("DegNorm core algorithm starts...")
    # Calculating Degradation Index
    output=.degnorm_baseline(read_coverage, counts, iteration, loop,
        down_sampling, grid_size,scale,counts.weighted,cores)
    message("DegNorm core algorithm done!")
    class(output)="DegNormClass"
    return(output)
}

##degnorm baseline selection algorithm
.degnorm_baseline=function(read_coverage, counts, iteration, loop,
                        down_sampling, grid_size,scale,counts.weighted,cores){
    for(i in seq_len(iteration)){
        cat(paste0("   DegNorm iteration ",i), "\n")
        cores <- detectCores(logical = FALSE) #utilize multiple cores
        #cl <- makeCluster(cores)
        cl <- makeCluster(cores-1)
        registerDoParallel(cl)
        res = NULL;j=1;m=dim(counts)[1];n=dim(counts)[2]
        if(down_sampling==0){
            res = foreach(j = seq_len(m), .multicombine = TRUE,
                    .export = c(".optiNMFCPP",".NMFCPP",".bin_drop")) %dopar%{
                    .optiNMFCPP(read_coverage[[j]], scale, loop)}
        }else if (down_sampling==1){
            res = foreach(j = seq_len(m), .multicombine = TRUE,
            .export = c(".optiNMFCPP_grid",".NMFCPP",".bin_drop")) %dopar%{
            .optiNMFCPP_grid(read_coverage[[j]], scale, loop, grid_size)}
        }else{
            stop("Error:down_sampling argument should either 0 or 1!")
        }
        stopImplicitCluster()
        stopCluster(cl)
        rho = t(do.call(cbind, lapply(res, function(x) x$rho)))
        # in case rho contains NA, replace entire row by 0
        rho[apply(rho,1,function(x) sum(is.na(x))>0),]=0
        K = t(do.call(cbind, lapply(res, function(x) x$K)))
        convergence = vapply(res, function(x) x$convergence,numeric(1))
        envelop= lapply(res, function(x) x$envelop)
        dimnames(rho)=dimnames(counts);dimnames(K)=dimnames(counts);
        names(convergence)=row.names(counts); names(envelop)=row.names(counts)
        # Capping the DI scores at 0.9 to adjust for over-inflation
        rho[rho > 0.9] = 0.9; rho[rho < 0] = 0
        # Getting iterative DegNorm normalized counts
        counts_normed = counts.weighted / (1 - rho)
        # Applying sample-wise overall DI scores to low counts genes
        # which did not go through NMF
        rho[apply(rho,1,max)==0,] =
        matrix(rep(1-colSums(counts.weighted)/colSums(counts_normed),
                        sum(apply(rho,1,max)==0)),byrow = TRUE, ncol = n)
        # Updating sequencing depth scale
        norm.factor = colSums(counts_normed)
        norm.factor = norm.factor / median(norm.factor)
        counts.weighted = sweep(counts.weighted, 2, norm.factor, FUN = "/")
        scale = scale * norm.factor
    }
    output = list("counts"=counts,"counts_normed"=round(counts_normed),
                    "DI"=rho,"K"=K,"convergence"=convergence,
                    "envelop"=envelop)
    return(output)
}
