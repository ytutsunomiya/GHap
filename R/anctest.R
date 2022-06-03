#Function: ghap.anctest
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Predict ancestry of haplotypes

ghap.anctest <- function(
  object,
  blocks = NULL,
  prototypes,
  test = NULL,
  only.active.samples = TRUE,
  only.active.markers = TRUE,
  ncores = 1,
  verbose = TRUE
){
  
  # Check if object is of class GHap.phase -------------------------------------
  if(inherits(object, "GHap.phase") == FALSE){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Check if inactive markers and samples should be reactived ------------------
  if(only.active.markers == FALSE){
    object$marker.in <- rep(TRUE,times=object$nmarkers)
    object$nmarkers.in <- length(which(object$marker.in))
  }
  if(only.active.samples == FALSE){
    object$id.in <- rep(TRUE,times=2*object$nsamples)
    object$nsamples.in <- length(which(object$id.in))/2
  }
  
  # Organize prototype dataframe -----------------------------------------------
  nprotmrk <- length(which(prototypes$MARKER %in% object$marker))
  if(nprotmrk != nrow(prototypes)){
    emsg <- "\nMarkers in the prototypes should also be present in the object."
    stop(emsg)
  }
  if(identical(prototypes$MARKER, object$marker) == FALSE){
    protmrk <- prototypes$MARKER
    tmp <- data.frame(IDX = 1:object$nmarkers, MARKER = object$marker,
                      stringsAsFactors = FALSE)
    prototypes <- merge(x = tmp, y = prototypes, by = "MARKER", all.x=TRUE)
    prototypes <- prototypes[order(prototypes$IDX),-2]
    object$marker.in <- object$marker %in% protmrk & object$marker.in == TRUE
  }
  
  # Map test samples -----------------------------------------------------------
  test.idx <- which(object$id %in% test & object$id.in == TRUE)
  
  # Check if blocks exist ------------------------------------------------------
  if(is.null(blocks) == TRUE){
    
    # Calculate marker density
    mrkdist <- diff(object$bp)
    mrkdist <- mrkdist[which(mrkdist > 0)]
    density <- mean(mrkdist)
    
    # Generate blocks for admixture events
    g <- 10
    window <- (100e+6)/(2*g)
    window <- ceiling(window/density)
    step <- ceiling(window/4)
    blocks <- ghap.blockgen(object, windowsize = window,
                            slide = step, unit = "marker")
  }
  
  # Initialize block iteration function ----------------------------------------
  blockfun <- function(b){
    
    #Get block info
    block.info <- blocks[b, c("BLOCK","CHR","BP1","BP2")]
    
    #SNPs in the block
    snps <- which(object$chr == block.info$CHR &
                    object$bp >= block.info$BP1 &
                    object$bp <= block.info$BP2 &
                    object$marker.in == TRUE)
    blocksize <- length(snps)
    
    #Get test haplotypes
    Mtst <- ghap.slice(object = object, ids = test.idx, variants = snps,
                       index = TRUE, verbose = FALSE)
    Mref <- prototypes[snps,-1]
    
    #Prediction
    pred <- rep(NA, times = ncol(Mtst))
    sq <- rep(NA, times = ncol(Mref))
    for(h in 1:ncol(Mtst)){
      for(m in 1:ncol(Mref)){
        sq[m] <- sum((Mtst[,h] - Mref[,m])^2)
      }
      pred[h] <- which(sq == min(sq))
    }
    pred <- colnames(Mref)[pred]
    ids <- object$id[test.idx]
    pops <- object$pop[test.idx]
    names(pred) <- ids
    
    #Make output
    ids <- ids[1:length(ids) %% 2 == 0]
    pops <- pops[1:length(pops) %% 2 == 0]
    out <- rep(NA, times=length(ids)*8)
    for(i in 1:length(ids)){
      haps <- pred[which(names(pred) == ids[i])]
      haps <- unlist(c(block.info,pops[i],ids[i],haps[1],haps[2]))
      haps <- as.vector(haps)
      out[(i*8 - 7):(i*8)] <- haps
    }
    
    #Return output
    return(out)
  }
  
  
  
  # Compute ancestry -----------------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  if(verbose == TRUE){
    cat("\nPredicting ancestry of haplotypes... ")
  }
  if(ncores == 1){
    results <- lapply(FUN = blockfun, X = 1:nrow(blocks))
  }else{
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      clusterEvalQ(cl, library(Matrix))
      varlist <- list("blocks","object","test.idx")
      clusterExport(cl = cl, varlist = varlist, envir=environment())
      results <- unlist(parLapply(cl = cl, fun = blockfun, X = 1:nrow(blocks)))
      stopCluster(cl)
    }else{
      results <- mclapply(FUN = blockfun, X = 1:nrow(blocks), mc.cores = ncores)
    }
  }
  if(verbose == TRUE){
    cat("Done.\n")
    cat("Assembling results... ")
  }
  results <- data.frame(matrix(unlist(results), ncol=8, byrow=TRUE),
                        stringsAsFactors = F)
  colnames(results) <- c("BLOCK","CHR","BP1","BP2","POP","ID","HAP1","HAP2")
  results$BP1 <- as.numeric(results$BP1)
  results$BP2 <- as.numeric(results$BP2)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  # Return results -------------------------------------------------------------
  return(results)
  
  
}
