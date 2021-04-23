#Function: ghap.blup
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: calculate GBLUP solution for each haplotype allele

ghap.blup<-function(
  gebvs, 
  haplo,
  invcov,
  gebvsweights = NULL,
  haploweights = NULL,
  only.active.alleles = TRUE,
  batchsize = NULL,
  ncores = 1,
  verbose = TRUE
){
  
  #General data check
  if (class(haplo) != "GHap.haplo") {
    stop("Argument haplo must be a GHap.haplo object.")
  }
  if (only.active.alleles == FALSE) {
    haplo$allele.in <- rep(TRUE, times = haplo$nalleles)
    haplo$nalleles.in <- length(which(haplo$allele.in))
  }
  activealleles <- which(haplo$allele.in)
  if (is.null(haploweights) == TRUE) {
    haploweights <- rep(1,times=haplo$nalleles.in)
  }
  if (length(haploweights) != haplo$nalleles.in) {
    stop("Vector of haplotype weights must have the same length as the number of haplotype alleles.")
  }
  names(haploweights) <- activealleles
  ids <- rep(NA, times = length(gebvs))
  for (i in 1:length(ids)) {
    ids[i] <- which(haplo$id == names(gebvs)[i])
  }
  if (length(which(names(gebvs) %in% haplo$id)) != length(gebvs)) {
    stop("All ids in the vector of GEBVs must be present in the GHap.haplo object.")
  }
  invcov <- invcov[names(gebvs),names(gebvs)]
  if(identical(colnames(invcov),names(gebvs)) == FALSE){
    stop("Names declared in random effects and covariance matrix do not match!")
  }
  if (is.null(gebvsweights) == TRUE) {
    gebvsweights <- rep(1,times=length(gebvs))
  }
  if (length(gebvsweights) != length(gebvsweights)) {
    stop("Vector of GEBV weights must have the same length as the number of observations.")
  }
  names(gebvsweights) <- names(gebvs)
  k <- invcov%*%Diagonal(x=gebvsweights/mean(gebvsweights))%*%gebvs
  
  # Get number of cores
  if(Sys.info()["sysname"] == "Windows"){
    if(ncores > 1 & verbose == TRUE){
      cat("\nParallelization not supported yet under Windows (using a single core).")
    }
    ncores <- 1
  }else{
    ncores <- min(c(detectCores(), ncores))
  }
  
  # Generate lookup table
  lookup <- rep(NA,times=256)
  lookup[1:2] <- c(0,1)
  d <- 10
  i <- 3
  while(i <= 256){
    b <- d + lookup[1:(i-1)]
    lookup[i:(length(b)+i-1)] <- b
    i <- i + length(b)
    d <- d*10
  }
  lookup <- sprintf(fmt="%08d", lookup)
  lookup <- sapply(lookup, function(i){intToUtf8(rev(utf8ToInt(i)))})
  
  #Generate batch index
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(haplo$nalleles.in/10)
  }
  if(batchsize > haplo$nalleles.in){
    batchsize <- haplo$nalleles.in
  }
  activealleles <- which(haplo$allele.in)
  nbatches <- round(haplo$nalleles.in/(batchsize),digits=0) + 1
  mybatch <- paste("B",1:nbatches,sep="")
  batch <- rep(mybatch,each=batchsize)
  batch <- batch[1:haplo$nalleles.in]
  mybatch <- unique(batch)
  nbatches <- length(mybatch)
  
  #Log message
  if(verbose == TRUE){
    cat("Processing ", haplo$nalleles.in, " HapAlleles in ", nbatches, " batches.\n", sep="")
    cat("Inactive alleles will be ignored.\n")
  }
  
  #Auxiliary functions
  varfun <- function(j){
    x <- hap.geno[j,]
    xname <- rownames(hap.geno)[j]
    out <- haploweights[xname]*var(x)
    return(out)
  }
  freqfun <- function(j){
    return(sum(hap.geno[j,])/(2*haplo$nsamples.in))
  }
  gblup.FUN <- function(j) {
    x <- hap.geno[j,]
    xname <- rownames(hap.geno)[j]
    cent <- mean(x)
    x <- x - cent
    b <- sum(haploweights[xname]*x*k)
    varxb <- var(x*b)
    return(c(b,varxb,cent))
  }
  
  #Iterate batches
  hapreg <- NULL
  hapreg$BLOCK <- haplo$block[haplo$allele.in]
  hapreg$CHR <- haplo$chr[haplo$allele.in]
  hapreg$BP1 <- haplo$bp1[haplo$allele.in]
  hapreg$BP2 <- haplo$bp2[haplo$allele.in]
  hapreg$ALLELE <- haplo$allele[haplo$allele.in]
  hapreg$SCORE <- rep(NA, times=haplo$nalleles.in)
  hapreg$FREQ <- rep(NA, times=haplo$nalleles.in)
  hapreg$VAR <- rep(NA, times=haplo$nalleles.in)
  hapreg$pVAR <- rep(NA, times=haplo$nalleles.in)
  hapreg$CENTER <- rep(NA, times=haplo$nalleles.in)
  hapreg$SCALE <- rep(NA, times=haplo$nalleles.in)
  sumalleles <- 0
  sumvar <- 0
  for(i in 1:nbatches){
    idx <- which(batch == mybatch[i])
    slice <- activealleles[idx]
    hap.geno <- ghap.hslice(haplo = haplo, ids = ids, alleles = slice,
                            index = TRUE, lookup = lookup, ncores = ncores)
    if(Sys.info()["sysname"] == "Windows"){
      sumvar <- sumvar + sum(unlist(lapply(X = 1:nrow(hap.geno), FUN = varfun)))
      hapreg$FREQ[idx] <- unlist(lapply(X = 1:nrow(hap.geno), FUN = freqfun))
      a <- unlist(lapply(X = 1:nrow(hap.geno), FUN = gblup.FUN))
      a <- data.frame(matrix(a, nrow=nrow(hap.geno), byrow=TRUE))
      hapreg$SCORE[idx] <- a[,1]
      hapreg$VAR[idx] <- a[,2]
      hapreg$CENTER[idx] <- a[,3]
      hapreg$SCALE[idx] <- 1
    }else{
      sumvar <- sumvar + sum(unlist(mclapply(X = 1:nrow(hap.geno), FUN = varfun, mc.cores = ncores)))
      hapreg$FREQ[idx] <- unlist(mclapply(X = 1:nrow(hap.geno), FUN = freqfun, mc.cores = ncores))
      a <- unlist(mclapply(X = 1:nrow(hap.geno), FUN = gblup.FUN, mc.cores = ncores))
      a <- data.frame(matrix(a, nrow=nrow(hap.geno), byrow=TRUE))
      hapreg$SCORE[idx] <- a[,1]
      hapreg$VAR[idx] <- a[,2]
      hapreg$CENTER[idx] <- a[,3]
      hapreg$SCALE[idx] <- 1
    }
    if(verbose == TRUE){
      sumalleles <- sumalleles + length(idx)
      cat(sumalleles, "HapAlleles processed.\r")
    }
  }
  hapreg$SCORE <- hapreg$SCORE/sumvar
  hapreg$VAR <- hapreg$VAR*(1/sumvar)^2
  hapreg$pVAR <- hapreg$VAR/sum(hapreg$VAR)
  
  #Return results
  hapreg <- data.frame(hapreg, stringsAsFactors = FALSE)
  return(hapreg)
}