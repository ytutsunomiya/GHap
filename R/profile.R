#Function: ghap.profile
#License: GPLv3 or later
#Modification date: 28 Apr 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute individual profiles based on HapAllele or marker scores

ghap.profile2 <- function(
  phase = NULL,
  haplo = NULL,
  score,
  only.active.samples=TRUE,
  batchsize=NULL,
  ncores=1,
  verbose=TRUE
){
  
  #Check if phase or haplo was provided
  if(is.null(phase) & is.null(haplo)){
    stop("A GHap.phase or GHap.haplo must be provided.")
  }
  
  #Check if both phase and haplo were provided
  if(is.null(phase) == FALSE & is.null(haplo) == FALSE){
    stop("Provide only one GHap.phase or GHap.haplo object.")
  }
  
  #Scoring for haplo
  if(is.null(haplo) == FALSE){
    
    #Check if haplo is a GHap.haplo object
    if(class(haplo) != "GHap.haplo"){
      stop("Argument haplo must be a GHap.haplo object.")
    }
    
    #Check if inactive alleles and samples should be reactived
    if(only.active.samples == FALSE){
      haplo$id.in <- rep(TRUE,times=haplo$nsamples)
      haplo$nsamples.in<-length(which(haplo$id.in))
    }
    
    #Prepare score vectors
    score$HAPLOtmp <- paste(score$BLOCK,score$CHR,score$BP1,score$BP2,score$ALLELE,sep="_")
    haps <- NULL
    haps$HAPLOtmp <- paste(haplo$block,haplo$chr,haplo$bp1,haplo$bp2,haplo$allele,sep="_")
    haps$IDX <- 1:length(haps$HAPLOtmp)
    haps <- data.frame(haps,stringsAsFactors = FALSE)
    haps <- merge(x = haps, y = score, by.x = "HAPLOtmp", by.y = "HAPLOtmp", all.x=TRUE,sort=FALSE)
    haps <- haps[is.na(haps$SCORE) == FALSE,]
    
    #Log message
    if(nrow(haps) != nrow(score)){
      stop("From ",nrow(score)," HapAlleles declared but ",nrow(haps)," were found.\n")
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
    
    #auxiliary function
    score.FUN <- function(j){
      x <- hap.geno[j,]
      xname <- as.numeric(rownames(hap.geno)[j])
      row <- which(haps$IDX == xname)
      x <- (x - haps$CENTER[row])/haps$SCALE[row]
      return(x*haps$SCORE[row])
    }
    
    #Iterate batches
    ncores <- min(c(detectCores(), ncores))
    out <- NULL
    out$POP <- haplo$pop[haplo$id.in]
    out$ID <- haplo$id[haplo$id.in]
    out$SCORE <- rep(0, times = haplo$nsamples.in)
    sumalleles <- 0
    for(i in 1:nbatches){
      idx <- which(batch == mybatch[i])
      slice <- activealleles[idx]
      hap.geno <- ghap.hslice(haplo = haplo, ids = which(haplo$id.in), alleles = slice,
                              index = TRUE, lookup = lookup, ncores = ncores)
      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(ncores)
        a <- unlist(parLapply(cl = cl, fun = score.FUN, X = 1:nrow(hap.geno)))
        stopCluster(cl)
        a <- data.frame(matrix(a, nrow=nrow(hap.geno), byrow=TRUE))
        a <- colSums(a)
        out$SCORE <- out$SCORE + a
      }else{
        a <- unlist(mclapply(X = 1:nrow(hap.geno), FUN = score.FUN, mc.cores = ncores))
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
  
  #Scoring for phase
  if(is.null(phase) == FALSE){
    
    #Check if phase is a GHap.phase object
    if(is.null(phase) == FALSE & class(phase) != "GHap.phase"){
      stop("Argument phase must be a GHap.phase object.")
    }
    
    #Check if inactive samples should be reactived
    if(only.active.samples == FALSE){
      phase$id.in <- rep(TRUE,times=2*phase$nsamples)
      phase$nsamples.in <- length(which(phase$id.in))/2
    }
    
    #Check if all markers exist
    idx <- which(score$MARKER %in% phase$marker)
    if(length(idx) != nrow(score)){
      emsg <- paste("\nFrom", nrow(score), "markers declared in the score data frame only",
                    length(idx), "were present in the GHap.phase object.")
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
    midx <- which(phase$marker %in% score$MARKER)
    markers <- data.frame(MARKER = phase$marker[midx], A0 = phase$A0[midx], A1 = phase$A1[midx],
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
    
    # Initialize lookup table
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
    ncores <- min(c(detectCores(), ncores))
    
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
    score.FUN <- function(j, phase.geno, score){
      x <- phase.geno[j,]
      xname <- rownames(phase.geno)[j]
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
    out <- unique(cbind(phase$pop[phase$id.in], phase$id[phase$id.in]))
    out <- data.frame(out, stringsAsFactors = FALSE)
    colnames(out) <- c("POP","ID")
    out$SCORE <- rep(0, times = phase$nsamples.in)
    summarkers <- 0
    for(i in 1:nbatches){
      idx <- id1[i]:id2[i]
      phase.geno <- ghap.pslice(phase = phase, ids = out$ID, markers = score$MARKER[idx],
                                lookup = lookup, ncores = ncores, unphase = TRUE)
      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(ncores)
        a <- unlist(parLapply(cl = cl, fun = score.FUN, X = 1:nrow(phase.geno),
                              phase.geno = phase.geno, score = score[idx,]))
        stopCluster(cl)
        a <- data.frame(matrix(a, nrow=nrow(phase.geno), byrow=TRUE))
        a <- colSums(a)
        out$SCORE <- out$SCORE + a
      }else{
        a <- unlist(mclapply(X = 1:nrow(phase.geno), FUN = score.FUN, mc.cores = ncores,
                             phase.geno = phase.geno, score = score[idx,]))
        a <- data.frame(matrix(a, nrow=nrow(phase.geno), byrow=TRUE))
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
