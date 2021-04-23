#Function: ghap.pslice
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Get a slice of the phase object

ghap.pslice <- function(
  phase,
  ids,
  markers,
  index=FALSE,
  unphase=FALSE,
  lookup=NULL,
  ncores=1,
  verbose=TRUE
){
  
  #Check if phase is a GHap.phase object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Generate lookup table
  if(is.null(lookup)){
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
  }
  
  # Calculate offset and bitloss
  offset <- ceiling((2*phase$nsamples)/8)
  bitloss <- 8 - ((2*phase$nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  
  # Get indices
  if(index == TRUE){
    
    # Retrieve indices
    iidx <- ids
    midx <- markers
    
    # Check if ids and markers are within range
    if(max(iidx) > 2*phase$nsamples){
      stop("Some of the provided ids are out of range")
    }
    if(max(midx) > phase$nmarkers){
      stop("Some of the provided markers are out of range")
    }
    
    # Organize indices
    names(iidx) <- phase$id[iidx]
    names(midx) <- phase$marker[midx]
    
  }else{
    
    # Retrieve indices
    iidx <- which(phase$id %in% ids)
    midx <- which(phase$marker %in% markers)
    
    # Check if ids and markers exist
    if(length(iidx) != 2*length(ids)){
      stop("Some of the provided ids were not found")
    }
    if(length(midx) != length(markers)){
      stop("Some of the provided marker names were not found")
    }
    
    # Organize indices
    iidx1 <- iidx[1:length(iidx) %% 2 == 1]
    iidx2 <- iidx[1:length(iidx) %% 2 == 0]
    names(iidx1) <- phase$id[iidx1]
    names(iidx2) <- phase$id[iidx2]
    iidx1 <- iidx1[ids]
    iidx2 <- iidx2[ids]
    iidx <- rep(NA, times = length(iidx))
    iidx[1:length(iidx) %% 2 == 1] <- iidx1
    iidx[1:length(iidx) %% 2 == 0] <- iidx2
    names(iidx)[1:length(iidx) %% 2 == 1] <- names(iidx1)
    names(iidx)[1:length(iidx) %% 2 == 0] <- names(iidx2)
    names(midx) <- phase$marker[midx]
    midx <- midx[markers]
    
  }
  
  # Get bits
  X <- matrix(data = NA, nrow = length(midx), ncol = length(iidx))
  getBitFun <- function(i){
    phase.con <- file(phase$phase, "rb")
    a <- seek(con = phase.con, where =offset*(midx[i]-1), origin = 'start',rw = 'r')
    geno <- readBin(phase.con, what=integer(), size = 1, n = offset, signed = FALSE)
    geno <- paste(lookup[geno+1], collapse = "")
    geno <- unlist(strsplit(x = geno, split = ""))
    geno <- as.integer(geno[1:(length(geno)-bitloss)])
    close.connection(phase.con)
    return(geno[iidx])
  }
  ncores <- min(c(detectCores(), ncores))
  if (Sys.info()["sysname"] == "Windows") {
    if(ncores > 1 & verbose == TRUE){
      cat("\nParallelization not supported yet under Windows (using a single core).\n")
    }
    X <- unlist(lapply(X = 1:length(midx), FUN = getBitFun))
    X <- matrix(data = X, nrow = length(midx), ncol = length(iidx), byrow = TRUE)
  }else{
    X <- unlist(mclapply(X = 1:length(midx), FUN = getBitFun, mc.cores = ncores))
    X <- matrix(data = X, nrow = length(midx), ncol = length(iidx), byrow = TRUE)
  }
  colnames(X) <- names(iidx)
  rownames(X) <- names(midx)
  
  # Unphase genotypes
  if(unphase == TRUE){
    cols <- 1:ncol(X) %% 2
    X <- X[,which(cols == 1)] + X[,which(cols == 0)]
  }
  
  # Return X
  return(X)
  
}