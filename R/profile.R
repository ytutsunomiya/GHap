#Function: ghap.profile
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute individual profiles based on HapAllele scores

ghap.profile <- function(
  score,
  haplo,
  only.active.samples=TRUE,
  batchsize=NULL,
  ncores=1,
  verbose=TRUE
){
  
  
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
  
  #auxiliary function
  score.FUN <- function(j){
    x <- hap.geno[j,]
    xname <- as.numeric(rownames(hap.geno)[j])
    row <- which(haps$IDX == xname)
    x <- (x - haps$CENTER[row])/haps$SCALE[row]
    return(x*haps$SCORE[row])
  }
  
  #Iterate batches
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
      a <- unlist(lapply(X = 1:nrow(hap.geno), FUN = score.FUN))
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
  
  #Return object
  return(out)
  
}