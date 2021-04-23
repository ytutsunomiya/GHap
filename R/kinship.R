#Function: ghap.kinship
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute haplotype covariance matrix

ghap.kinship<-function(
  haplo,
  weights=NULL,
  batchsize=NULL,
  only.active.samples=TRUE,
  only.active.alleles=TRUE,
  ncores=1,
  verbose=TRUE
){
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if inactive alleles and samples should be reactived
  if(only.active.alleles == FALSE){
    haplo$allele.in <- rep(TRUE,times=haplo$nalleles)
    haplo$nalleles.in <- length(which(haplo$allele.in))
  }
  if(only.active.samples == FALSE){
    haplo$id.in <- rep(TRUE,times=haplo$nsamples)
    haplo$nsamples.in<-length(which(haplo$id.in))
  }
  
  #Check weights
  if (is.null(weights) == TRUE) {
    weights <- rep(1,times=haplo$nalleles.in)
  }
  if (length(weights) != haplo$nalleles.in) {
    stop("Vector of weights must have the same length as the number of haplotype alleles.")
  }
  
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
  
  #Initialize kinship matrix
  if(verbose == TRUE){
    cat("Preparing",haplo$nsamples.in,"x",haplo$nsamples.in,"kinship matrix.\n")
  }
  K <- Matrix(data = 0, nrow = haplo$nsamples.in, ncol = haplo$nsamples.in)
  K <- as(as(K,"dsyMatrix"),"dspMatrix")
  
  #Kinship iterate function
  kinship.FUN<-function(i){
    slice <- activealleles[batch == mybatch[i]]
    hap.geno <- ghap.hslice(haplo = haplo, ids = which(haplo$id.in), alleles = slice,
                            index = TRUE, lookup = lookup, ncores = ncores)
    Ztmp <- as(hap.geno, "dgeMatrix")
    Ztmp.mean <- apply(X = Ztmp,MARGIN = 1,FUN = mean)
    Ztmp <- (Ztmp - Ztmp.mean)
    K <- K + crossprod(Ztmp*sqrt(weights[batch == mybatch[i]]))
    return(K)
  }
  
  #Loop by batch
  sumalleles <- 0
  for(i in 1:nbatches){
    K <- kinship.FUN(i)
    if(verbose == TRUE){
      sumalleles <- sumalleles + length(which(batch == mybatch[i]))
      cat(sumalleles, "HapAlleles processed.\r")
    }
  }
  
  
  #Scale kinship matrix
  #K <- K/haplo$nalleles.in
  #varfun <- function(j) return(var(haplo$genotypes[j,haplo$id.in]))
  #q <- unlist(mclapply(X=activealleles,FUN = varfun))
  #q <- sum(q*weights)
  q <- mean(diag(K))
  K <- K/q
  colnames(K) <- haplo$id[haplo$id.in]
  rownames(K) <- colnames(K)
  
  #Return output
  return(K)
  
}