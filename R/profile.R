#Function: ghap.profile
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute individual profiles based on HapAllele or marker scores

ghap.profile <- function(
  object,
  score,
  only.active.samples=TRUE,
  batchsize=NULL,
  ncores=1,
  verbose=TRUE
){
  
  # Check if object is a GHap.phase object -------------------------------------
  obtypes <- c("GHap.phase","GHap.haplo","GHap.plink")
  if(inherits(object, obtypes) == FALSE){
    stop("\nInput data must be a valid GHap object (phase, haplo or plink).")
  }
  fac <- c(2,1,1)
  names(fac) <- obtypes
  ncores <- min(c(detectCores(), ncores))
  
  # Scoring for phase ----------------------------------------------------------
  if(inherits(object, c("GHap.phase","GHap.plink"))){
    
    #Check if inactive samples should be reactived
    if(only.active.samples == FALSE){
      object$id.in <- rep(TRUE,times=fac[class(object)]*object$nsamples)
      object$nsamples.in <- length(which(object$id.in))/fac[class(object)]
    }
    
    #Check if all markers exist
    idx <- which(score$MARKER %in% object$marker)
    if(length(idx) != nrow(score)){
      emsg <- paste("\nFrom", nrow(score), "markers declared in the score data frame only",
                    length(idx), "were present in the GHap object.")
      stop(emsg)
    }
    
    #Check if missing scores are present
    score <- score[which(is.na(score$SCORE) == FALSE),]
    if(length(idx) != nrow(score)){
      emsg <- paste("\nFrom", length(idx), "markers declared in the score data frame",
                    length(idx)-nrow(score), "had missing values.")
      stop(emsg)
    }
    
    #Check if alleles exist
    midx <- which(object$marker %in% score$MARKER)
    markers <- data.frame(MARKER = object$marker[midx], A0 = object$A0[midx], A1 = object$A1[midx],
                          stringsAsFactors = FALSE)
    score <- merge(x = score, y = markers, by = "MARKER", sort=FALSE)
    score$isA0 <- as.integer(score$A0 == score$ALLELE)
    score$isA1 <- as.integer(score$A1 == score$ALLELE)
    score$isSUM <- score$isA0 + score$isA1
    idx <- which(score$isSUM != 1)
    if(length(idx) > 0){
      emsg <- paste("\nFrom", nrow(score), "markers declared in the score data frame",
                    length(idx), "had alien alleles.")
      stop(emsg)
    }
    markers <- 1:object$nmarkers
    names(markers) <- object$marker
    score$IDX <- markers[score$MARKER]
    score <- score[order(score$IDX),]
    
    #Generate batch index
    nmarkers <- nrow(score)
    if(is.null(batchsize) == TRUE){
      batchsize <- ceiling(nmarkers/10)
    }
    if(batchsize > nmarkers){
      batchsize <- nmarkers
    }
    id1 <- seq(1,nmarkers,by=batchsize)
    id2 <- (id1+batchsize)-1
    id1 <- id1[id2<=nmarkers]
    id2 <- id2[id2<=nmarkers]
    id1 <- c(id1,id2[length(id2)]+1)
    id2 <- c(id2,nmarkers)
    if(id1[length(id1)] > nmarkers){
      id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
    }
    nbatches <- length(id1)
    
    #Log message
    if(verbose == TRUE){
      cat("Processing ", nmarkers, " markers in ", nbatches, " batches.\n", sep="")
    }
    
    #auxiliary function
    flip.FUN <- function(x){
      tmp <- x
      tmp[which(x == 0)] <- 2
      tmp[which(x == 2)] <- 0
      return(tmp)
    }
    
    #Iterate batches
    out <- unique(cbind(object$pop[object$id.in], object$id[object$id.in]))
    out <- data.frame(out, stringsAsFactors = FALSE)
    colnames(out) <- c("POP","ID")
    out$SCORE <- rep(0, times = object$nsamples.in)
    summarkers <- 0
    for(i in 1:nbatches){
      idx <- id1[i]:id2[i]
      scoretmpA0 <- score[idx,]
      scoretmpA1 <- scoretmpA0[which(scoretmpA0$isA1 == 1),]
      scoretmpA0 <- scoretmpA0[which(scoretmpA0$isA0 == 1),]
      if(nrow(scoretmpA0) > 0){
        Ztmp <- ghap.slice(object = object, ids = out$ID,
                           variants = scoretmpA0$MARKER, transposed = FALSE,
                           impute = TRUE, ncores = ncores, unphase = TRUE)
        Ztmp <- apply(X = Ztmp, MARGIN = 2, FUN = flip.FUN)
        Ztmp <- (Ztmp - scoretmpA0$CENTER)/scoretmpA0$SCALE
        out$SCORE <- out$SCORE + as.numeric(scoretmpA0$SCORE%*%Ztmp)
      }
      if(nrow(scoretmpA1) > 0){
        Ztmp <- ghap.slice(object = object, ids = out$ID,
                           variants = scoretmpA1$MARKER, transposed = FALSE,
                           impute = TRUE, ncores = ncores, unphase = TRUE)
        Ztmp <- (Ztmp - scoretmpA1$CENTER)/scoretmpA1$SCALE
        out$SCORE <- out$SCORE + as.numeric(scoretmpA1$SCORE%*%Ztmp)
      }
      if(verbose == TRUE){
        summarkers <- summarkers + nrow(scoretmpA0) + nrow(scoretmpA1)
        cat(summarkers, "markers processed.\r")
      }
    }
  }
  
  #Return object
  return(out)
  
}
