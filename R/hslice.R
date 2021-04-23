#Function: ghap.hslice
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Get a slice of the haplo object

ghap.hslice <- function(
  haplo,
  ids,
  alleles,
  index=FALSE,
  lookup=NULL,
  ncores=1,
  verbose=TRUE
){
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if alleles are indices
  if(is.numeric(alleles) == FALSE & is.integer(alleles) == FALSE){
    stop("Vector of alleles must be of class integer or numeric")
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
    lookup <- sapply(lookup, function(i){intToUtf8(rev(utf8ToInt(i)))})
  }
  
  # Calculate offset and bitloss
  offset <- ceiling((2*haplo$nsamples)/8)
  bitloss <- 8 - ((2*haplo$nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  
  # Get indices
  if(index == TRUE){
    
    # Retrieve indices
    iidx <- ids
    aidx <- alleles
    
    # Check if ids and markers are within range
    if(max(iidx) > haplo$nsamples){
      stop("Some of the provided ids are out of range")
    }
    if(max(aidx) > haplo$nalleles){
      stop("Some of the provided alleles are out of range")
    }
    
    # Organize indices
    names(iidx) <- haplo$id[iidx]
    names(aidx) <- as.character(aidx)
    
  }else{
    
    # Retrieve indices
    iidx <- which(haplo$id %in% ids)
    aidx <- which(1:haplo$nalleles %in% alleles)
    
    # Check if ids and markers exist
    if(length(iidx) != length(ids)){
      stop("Some of the provided ids were not found")
    }
    if(length(aidx) != length(alleles)){
      stop("Some of the provided marker names were not found")
    }
    
    # Organize indices
    names(iidx) <- haplo$id[iidx]
    iidx <- iidx[ids]
    aidx <- alleles
    names(aidx) <- as.character(aidx)
    
  }
  
  # Get bits
  X <- matrix(data = NA, nrow = length(aidx), ncol = length(iidx))
  getBitFun <- function(i){
    haplo.con <- file(haplo$genotypes, "rb")
    a <- seek(con = haplo.con, where = 3 + offset*(aidx[i]-1), origin = 'start',rw = 'r')
    geno <- readBin(haplo.con, what=integer(), size = 1, n = offset, signed = FALSE)
    geno <- paste(lookup[geno+1], collapse = "")
    nc <- nchar(geno)
    n <- seq(1, nc, by = 2)
    geno <- substring(geno, n, c(n[-1]-1, nc))
    geno[which(geno == "00")] <- "0"
    geno[which(geno == "01")] <- "1"
    geno[which(geno == "11")] <- "2"
    geno <- as.integer(geno[1:haplo$nsamples])
    close.connection(haplo.con)
    return(geno[iidx])
  }
  ncores <- min(c(detectCores(), ncores))
  if (Sys.info()["sysname"] == "Windows") {
    if(ncores > 1 & verbose == TRUE){
      cat("\nParallelization not supported yet under Windows (using a single core).\n")
    }
    X <- unlist(lapply(X = 1:length(aidx), FUN = getBitFun))
    X <- matrix(data = X, nrow = length(aidx), ncol = length(iidx), byrow = TRUE)
  }else{
    X <- unlist(mclapply(X = 1:length(aidx), FUN = getBitFun, mc.cores = ncores))
    X <- matrix(data = X, nrow = length(aidx), ncol = length(iidx), byrow = TRUE)
  }
  colnames(X) <- names(iidx)
  rownames(X) <- names(aidx)
  return(X)
  
}