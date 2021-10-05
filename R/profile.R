#Function: ghap.profile
#License: GPLv3 or later
#Modification date: 14 May 2021
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
  if(class(object) %in% obtypes == FALSE){
    stop("\nInput data must be a valid GHap object (phase, haplo or plink).")
  }
  fac <- c(2,1,1)
  names(fac) <- obtypes
  ncores <- min(c(detectCores(), ncores))
  
  # Scoring for haplo ----------------------------------------------------------
  if(class(object) == "GHap.haplo"){
    
    #Check if inactive alleles and samples should be reactived
    if(only.active.samples == FALSE){
      object$id.in <- rep(TRUE,times=object$nsamples)
      object$nsamples.in<-length(which(object$id.in))
    }
    
    #Prepare score vectors
    score$HAPLOtmp <- paste(score$BLOCK,score$CHR,
                            score$BP1,score$BP2,score$ALLELE,sep="_")
    haps <- NULL
    haps$HAPLOtmp <- paste(object$block,object$chr,
                           object$bp1,object$bp2,object$allele,sep="_")
    haps$IDX <- 1:length(haps$HAPLOtmp)
    haps <- data.frame(haps,stringsAsFactors = FALSE)
    haps <- merge(x = haps, y = score, by.x = "HAPLOtmp", by.y = "HAPLOtmp",
                  all.x=TRUE,sort=FALSE)
    haps <- haps[is.na(haps$SCORE) == FALSE,]
    
    #Log message
    if(nrow(haps) != nrow(score)){
      emsg <- paste0("From ",nrow(score)," HapAlleles declared but ",
                     nrow(haps)," were found.\n")
      stop(emsg)
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
      cat("Processing ", object$nalleles.in, " HapAlleles in ",
          nbatches, " batches.\n", sep="")
      cat("Inactive alleles will be ignored.\n")
    }
    
    #auxiliary function
    haplo_score.FUN <- function(j){
      x <- hap.geno[j,]
      xname <- as.numeric(rownames(hap.geno)[j])
      row <- which(haps$IDX == xname)
      x <- (x - haps$CENTER[row])/haps$SCALE[row]
      return(x*haps$SCORE[row])
    }
    
    #Iterate batches
    out <- NULL
    out$POP <- object$pop[object$id.in]
    out$ID <- object$id[object$id.in]
    out$SCORE <- rep(0, times = object$nsamples.in)
    sumalleles <- 0
    for(i in 1:nbatches){
      idx <- which(batch == mybatch[i])
      slice <- activealleles[idx]
      hap.geno <- ghap.slice(object = object, ids = which(object$id.in),
                             variants = slice, index = TRUE, ncores = ncores)
      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(ncores)
        clusterEvalQ(cl, library(Matrix))
        varlist <- list("hap.geno","haps")
        clusterExport(cl = cl, varlist = varlist, envir=environment())
        a <- unlist(parLapply(cl = cl, fun = haplo_score.FUN, X = 1:nrow(hap.geno)))
        stopCluster(cl)
        a <- data.frame(matrix(a, nrow=nrow(hap.geno), byrow=TRUE))
        a <- colSums(a)
        out$SCORE <- out$SCORE + a
      }else{
        a <- unlist(mclapply(X = 1:nrow(hap.geno), FUN = haplo_score.FUN, mc.cores = ncores))
        a <- data.frame(matrix(a, nrow=nrow(hap.geno), byrow=TRUE))
        a <- colSums(a)
        out$SCORE <- out$SCORE + a
      }
      if(verbose == TRUE){
        sumalleles <- sumalleles + length(idx)
        cat(sumalleles, "HapAlleles processed.\r")
      }
    }
    out <- data.frame(out,stringsAsFactors = FALSE)
    
  }
  
  # Scoring for phase ----------------------------------------------------------
  if(class(object) %in% c("GHap.phase","GHap.plink")){
    
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
    marker_score.FUN <- function(j){
      x <- marker.geno[j,]
      xname <- rownames(marker.geno)[j]
      row <- which(score$MARKER == xname)
      if(score$isA0[row] == 1){
        tmp <- x
        tmp[which(x == 0)] <- 2
        tmp[which(x == 2)] <- 0
        x <- tmp
      }
      x <- (x - score$CENTER[row])/score$SCALE[row]
      return(x*score$SCORE[row])
    }
    
    #Iterate batches
    out <- unique(cbind(object$pop[object$id.in], object$id[object$id.in]))
    out <- data.frame(out, stringsAsFactors = FALSE)
    colnames(out) <- c("POP","ID")
    out$SCORE <- rep(0, times = object$nsamples.in)
    summarkers <- 0
    for(i in 1:nbatches){
      idx <- id1[i]:id2[i]
      marker.geno <- ghap.slice(object = object, ids = out$ID, variants = score$MARKER[idx],
                                impute = TRUE, ncores = ncores, unphase = TRUE)
      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(ncores)
        clusterEvalQ(cl, library(Matrix))
        varlist <- list("marker.geno","score")
        clusterExport(cl = cl, varlist = varlist, envir=environment())
        a <- unlist(parLapply(cl = cl, fun = marker_score.FUN, X = 1:nrow(marker.geno)))
        stopCluster(cl)
        a <- data.frame(matrix(a, nrow=nrow(marker.geno), byrow=TRUE))
        a <- colSums(a)
        out$SCORE <- out$SCORE + a
      }else{
        a <- unlist(mclapply(X = 1:nrow(marker.geno), FUN = marker_score.FUN, mc.cores = ncores))
        a <- data.frame(matrix(a, nrow=nrow(marker.geno), byrow=TRUE))
        a <- colSums(a)
        out$SCORE <- out$SCORE + a
      }
      if(verbose == TRUE){
        summarkers <- summarkers + length(idx)
        cat(summarkers, "markers processed.\r")
      }
    }
  }
  
  #Return object
  return(out)
  
}
