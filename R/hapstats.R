#Function: ghap.hapstats
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Summary statistics for haplotype alleles

ghap.hapstats <- function(
  object,
  alpha=c(1,1),
  batchsize=NULL,
  only.active.samples=TRUE,
  only.active.alleles=TRUE,
  ncores=1,
  verbose=TRUE
){
  
  
  #Check if object is a GHap.haplo object
  if(inherits(object, "GHap.haplo") == FALSE){
    stop("Argument object must be a GHap.haplo object.")
  }
  
  #Check if inactive alleles and samples should be reactived
  if(only.active.alleles == FALSE){
    object$allele.in <- rep(TRUE,times=object$nalleles)
    object$nalleles.in <- length(which(object$allele.in))
  }
  if(only.active.samples == FALSE){
    object$id.in <- rep(TRUE,times=object$nsamples)
    object$nsamples.in <- length(which(object$id.in))
  }
  
  #Generate batch index
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(object$nalleles.in/10)
  }
  if(batchsize > object$nalleles.in){
    batchsize <- object$nalleles.in
  }
  activealleles <- which(object$allele.in)
  nbatches <- round(object$nalleles.in/(batchsize),digits=0) + 1
  mybatch <- paste("B",1:nbatches,sep="")
  batch <- rep(mybatch,each=batchsize)
  batch <- batch[1:object$nalleles.in]
  mybatch <- unique(batch)
  nbatches <- length(mybatch)
  
  #Log message
  if(verbose == TRUE){
    cat("Processing ", object$nalleles.in, " HapAlleles in ", nbatches, " batches.\n", sep="")
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
  ncores <- min(c(detectCores(), ncores))
  hapstats <- NULL
  hapstats$BLOCK <- object$block[object$allele.in]
  hapstats$CHR <- object$chr[object$allele.in]
  hapstats$BP1 <- object$bp1[object$allele.in]
  hapstats$BP2 <- object$bp2[object$allele.in]
  hapstats$ALLELE<- object$allele[object$allele.in]
  hapstats$N <- rep(NA, times=object$nalleles.in)
  hapstats$FREQ <- rep(NA, times=object$nalleles.in)
  hapstats$O.HOM <- rep(NA, times=object$nalleles.in)
  hapstats$O.HET <- rep(NA, times=object$nalleles.in)
  sumalleles <- 0
  for(i in 1:nbatches){
    idx <- which(batch == mybatch[i])
    slice <- activealleles[idx]
    hap.geno <- ghap.slice(object = object, ids = which(object$id.in),
                           variants = slice,
                           index = TRUE, ncores = ncores)
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      a <- unlist(parLapply(cl = cl, fun = hapstats.FUN, X = 1:nrow(hap.geno)))
      stopCluster(cl)
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
  hapstats$E.HOM <- (hapstats$FREQ^2)*object$nsamples.in
  hapstats$RATIO <- (hapstats$E.HOM+alpha[1])/(hapstats$O.HOM+alpha[2])
  hapstats$BIN.logP <- -1*pbinom(q = hapstats$O.HOM,size = object$nsamples.in,prob = hapstats$FREQ^2,lower.tail=TRUE,log.p = TRUE)/log(10)
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
