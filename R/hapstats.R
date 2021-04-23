#Function: ghap.hapstats
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Summary statistics for haplotype alleles

ghap.hapstats<-function(
  haplo,
  alpha=c(1,1),
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
    haplo$nsamples.in <- length(which(haplo$id.in))
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
  
  #Hapstats iterate function
  hapstats.FUN <- function(j){
    x <- hap.geno[j,]
    N <- sum(x)
    FREQ <- N/(2*length(x))
    O.HOM <- length(which(x == 2))
    O.HET <- length(which(x == 1))
    return(c(N,FREQ,O.HOM,O.HET))
  }
  
  #Iterate batches
  hapstats <- NULL
  hapstats$BLOCK <- haplo$block[haplo$allele.in]
  hapstats$CHR <- haplo$chr[haplo$allele.in]
  hapstats$BP1 <- haplo$bp1[haplo$allele.in]
  hapstats$BP2 <- haplo$bp2[haplo$allele.in]
  hapstats$ALLELE<- haplo$allele[haplo$allele.in]
  hapstats$N <- rep(NA, times=haplo$nalleles.in)
  hapstats$FREQ <- rep(NA, times=haplo$nalleles.in)
  hapstats$O.HOM <- rep(NA, times=haplo$nalleles.in)
  hapstats$O.HET <- rep(NA, times=haplo$nalleles.in)
  sumalleles <- 0
  for(i in 1:nbatches){
    idx <- which(batch == mybatch[i])
    slice <- activealleles[idx]
    hap.geno <- ghap.hslice(haplo = haplo, ids = which(haplo$id.in), alleles = slice,
                            index = TRUE, lookup = lookup, ncores = ncores)
    if(Sys.info()["sysname"] == "Windows"){
      a <- lapply(X = which(haplo$allele.in), FUN = hapstats.FUN)
    }else{
      a <- mclapply(FUN=hapstats.FUN, X=1:nrow(hap.geno), mc.cores = ncores)
    }
    a <- data.frame(matrix(unlist(a), nrow=nrow(hap.geno), byrow=TRUE))
    hapstats$N[idx] <- a[,1]
    hapstats$FREQ[idx] <- a[,2]
    hapstats$O.HOM[idx] <- a[,3]
    hapstats$O.HET[idx] <- a[,4]
    if(verbose == TRUE){
      sumalleles <- sumalleles + length(idx)
      cat(sumalleles, "HapAlleles processed.\r")
    }
  }
  hapstats$E.HOM <- (hapstats$FREQ^2)*haplo$nsamples.in
  hapstats$RATIO <- (hapstats$E.HOM+alpha[1])/(hapstats$O.HOM+alpha[2])
  hapstats$BIN.logP <- -1*pbinom(q = hapstats$O.HOM,size = haplo$nsamples.in,prob = hapstats$FREQ^2,lower.tail=TRUE,log.p = TRUE)/log(10)
  hapstats$POI.logP <- -1*ppois(q = hapstats$O.HOM,lambda = hapstats$E.HOM,lower.tail=TRUE,log.p = TRUE)/log(10)
  hapstats <- data.frame(hapstats,stringsAsFactors = FALSE)
  hapstats$TYPE <- NA
  for(i in unique(hapstats$BLOCK)){
    slice <- which(hapstats$BLOCK == i)
    freq <- hapstats$FREQ[slice]
    sumfreq <- sum(freq)
    nalleles <- length(slice)
    type <- rep("REGULAR",times=nalleles)
    type[freq == 0] <- "ABSENT"
    minfreq <- min(freq[freq != 0])
    maxfreq <- max(freq)
    if(sumfreq == 1 & nalleles > 2){
      type[which(freq == minfreq)[1]] <- "MINOR"
      type[which(freq == maxfreq)[1]] <- "MAJOR"
    }else if(sumfreq == 1 & nalleles == 2){
      if(freq[1] == freq[2]){
        type[1] <- "MINOR"
        type[2] <- "MAJOR"
      }else if(freq[1] != freq[2] & pmin(freq[1],freq[2]) != 0){
        type[which(freq == minfreq)] <- "MINOR"
        type[which(freq == maxfreq)] <- "MAJOR"
      }else if(maxfreq == 1){
        type[which(freq == 1)] <- "SINGLETON"
      }
    }else if(sumfreq == 1 & nalleles == 1){
      type <- "SINGLETON"
    # }else if(sumfreq != 1 & nalleles == 1){
    #   type <- "REGULAR"
    # }else{
    #   type[which(freq == maxfreq)[1]] <- "MAJOR"
    # }
    }else{
      type <- "REGULAR"
    }
    hapstats$TYPE[slice] <- type
  }

  #Return object
  return(hapstats)
  
}