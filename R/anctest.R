#Function: ghap.anctest
#License: GPLv3 or later
#Modification date: 26 Apr 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Predict ancestry of haplotypes

ghap.anctest <- function(
  phase,
  blocks,
  prototypes,
  test = NULL,
  only.active.samples = TRUE,
  only.active.markers = TRUE,
  ncores = 1,
  verbose = TRUE
){
  
  # Check if phase is a GHap.phase object-------------------------------------------------------------
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Check if inactive markers and samples should be reactived-----------------------------------------
  if(only.active.markers == FALSE){
    phase$marker.in <- rep(TRUE,times=phase$nmarkers)
    phase$nmarkers.in <- length(which(phase$marker.in))
  }
  if(only.active.samples == FALSE){
    phase$id.in <- rep(TRUE,times=2*phase$nsamples)
    phase$nsamples.in <- length(which(phase$id.in))/2
  }
  
  # Organize prototype dataframe ---------------------------------------------------------------------
  nprotmrk <- length(which(prototypes$MARKER %in% phase$marker))
  if(nprotmrk != nrow(prototypes)){
    stop("Markers listed in the prototypes dataframe should be present in the GHap.phase object.")
  }
  if(identical(prototypes$MARKER, phase$marker) == FALSE){
    protmrk <- prototypes$MARKER
    tmp <- data.frame(IDX = 1:phase$nmarkers, MARKER = phase$marker, stringsAsFactors = FALSE)
    prototypes <- merge(x = tmp, y = prototypes, by = "MARKER", all.x=TRUE)
    prototypes <- prototypes[order(prototypes$IDX),-2]
    phase$marker.in <- phase$marker %in% protmrk & phase$marker.in == TRUE
  }
  
  # Map test samples----------------------------------------------------------------------------------
  test.idx <- which(phase$id %in% test & phase$id.in == TRUE)
  
  # Initialize lookup table----------------------------------------------------------------------------
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
  
  # Initialize block iteration function---------------------------------------------------------------
  blockfun <- function(b){
    
    #Get block info
    block.info <- blocks[b, c("BLOCK","CHR","BP1","BP2")]
    
    #SNPs in the block
    snps <- which(phase$chr == block.info$CHR &
                    phase$bp >= block.info$BP1 &
                    phase$bp <= block.info$BP2 &
                    phase$marker.in == TRUE)
    blocksize <- length(snps)
    
    #Get test haplotypes
    Mtst <- ghap.pslice(phase = phase, ids = test.idx, markers = snps,
                        index = TRUE, lookup = lookup, verbose = FALSE)
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
    ids <- phase$id[test.idx]
    pops <- phase$pop[test.idx]
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
  
  # Compute ancestry----------------------------------------------------------------------------------
  if(verbose == TRUE){
    cat("\nPredicting ancestry of haplotypes... ")
  }
  ncores <- min(c(detectCores(), ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cl <- makeCluster(ncores)
    results <- unlist(parLapply(cl = cl, fun = blockfun, X = 1:nrow(blocks)))
    stopCluster(cl)
  }else{
    results <- mclapply(FUN = blockfun, X = 1:nrow(blocks), mc.cores = ncores)
  }
  if(verbose == TRUE){
    cat("Done.\n")
    cat("Assembling results... ")
  }
  results <- data.frame(matrix(unlist(results), ncol=8, byrow=TRUE), stringsAsFactors = F)
  colnames(results) <- c("BLOCK","CHR","BP1","BP2","POP","ID","HAP1","HAP2")
  results$BP1 <- as.numeric(results$BP1)
  results$BP2 <- as.numeric(results$BP2)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  # Return results------------------------------------------------------------------------------------
  return(results)
  
  
}
