#Function: ghap.bslice
#License: GPLv3 or later
#Modification date: 02 May 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Get a slice of a GHap.plink object

ghap.bslice <- function(
  plink,
  ids,
  markers,
  index=FALSE,
  transposed=FALSE,
  sparse=TRUE,
  impute=FALSE,
  lookup=NULL,
  ncores=1,
  verbose=TRUE
){
  
  # Check if plink is a GHap.plink object --------------------------------------
  if(class(plink) != "GHap.plink"){
    stop("Argument plink must be a GHap.plink object.")
  }
  
  # Generate lookup table ------------------------------------------------------
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
    lookup <- sapply(lookup, function(i){intToUtf8(rev(utf8ToInt(i)))})
  }
  
  # Calculate offset and bitloss -----------------------------------------------
  offset <- ceiling((2*plink$nsamples)/8)
  bitloss <- 8 - ((2*plink$nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  
  # Get indices ----------------------------------------------------------------
  if(index == TRUE){
    
    # Retrieve indices
    iidx <- ids
    midx <- markers
    
    # Check if ids and markers are within range
    if(max(iidx) > plink$nsamples){
      stop("Some of the provided sample indices are out of range")
    }
    if(max(midx) > plink$nmarkers){
      stop("Some of the provided marker indices are out of range")
    }
    
    # Organize indices
    names(iidx) <- plink$id[iidx]
    names(midx) <- plink$marker[midx]
    
  }else{
    
    # Retrieve indices
    iidx <- which(plink$id %in% ids)
    midx <- which(plink$marker %in% markers)
    
    # Check if ids and markers exist
    if(length(iidx) != length(ids)){
      stop("Some of the provided ids were not found")
    }
    if(length(midx) != length(markers)){
      stop("Some of the provided marker names were not found")
    }
    
    # Organize indices
    names(iidx) <- plink$id[iidx]
    iidx <- iidx[ids]
    names(midx) <- plink$marker[midx]
    midx <- midx[markers]
    
  }
  
  # Get bits
  getBitFun <- function(i){
    plink.con <- file(plink$plink, "rb")
    a <- seek(con = plink.con, where = 3 + offset*(midx[i]-1), origin = 'start',rw = 'r')
    geno <- readBin(plink.con, what=integer(), size = 1, n = offset, signed = FALSE)
    geno <- paste(lookup[geno+1], collapse = "")
    nc <- nchar(geno)
    n <- seq(1, nc, by = 2)
    geno <- substring(geno, n, c(n[-1]-1, nc))
    geno[which(geno == "00")] <- "2"
    geno[which(geno == "01")] <- "1"
    geno[which(geno == "11")] <- "0"
    if(impute == TRUE){
      geno[which(geno == "10")] <- "0"
    }else{
      geno[which(geno == "10")] <- NA
    }
    geno <- as.integer(geno[1:plink$nsamples])
    close.connection(plink.con)
    return(geno[iidx])
  }
  ncores <- min(c(detectCores(), ncores))
  if (Sys.info()["sysname"] == "Windows") {
    cl <- makeCluster(ncores) 
    X <- unlist(parLapply(cl = cl, fun = getBitFun, X = 1:length(midx)))
    stopCluster(cl)
  }else{
    X <- unlist(mclapply(X = 1:length(midx), FUN = getBitFun, mc.cores = ncores))
  }
  if(transposed == FALSE){
    X <- Matrix(data = X, nrow = length(midx), ncol = length(iidx), byrow = TRUE, sparse = sparse)
    colnames(X) <- names(iidx)
    rownames(X) <- names(midx)
  }else{
    X <- Matrix(data = X, ncol = length(midx), nrow = length(iidx), byrow = FALSE, sparse = sparse)
    rownames(X) <- names(iidx)
    colnames(X) <- names(midx)
  }
  return(X)
  
}
