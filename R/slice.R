#Function: ghap.slice
#License: GPLv3 or later
#Modification date: 12 May 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Get a slice of a GHap object

ghap.slice <- function(
  object,
  ids,
  variants,
  index=FALSE,
  transposed=FALSE,
  sparse=TRUE,
  unphase=FALSE,
  impute=FALSE,
  lookup=NULL,
  ncores=1,
  verbose=TRUE
){
  
  # Check if object is a GHap.phase object -------------------------------------
  obtypes <- c("GHap.phase","GHap.haplo","GHap.plink")
  if(class(object) %in% obtypes == FALSE){
    stop("\nInput data must be a valid GHap object (phase, haplo or plink).")
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
    if(class(object) != "GHap.phase"){
      lookup <- sapply(lookup, function(i){intToUtf8(rev(utf8ToInt(i)))})
    }
  }
  
  # Calculate offset and bitloss -----------------------------------------------
  offset <- ceiling((2*object$nsamples)/8)
  bitloss <- 8 - ((2*object$nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  
  # Get indices ----------------------------------------------------------------
  if(index == TRUE){
    
    # Retrieve indices
    iidx <- ids
    vidx <- variants
    
    # Organize indices for ids
    if(max(iidx) > 2*object$nsamples & class(object) == "GHap.phase"){
      stop("Some of the provided ids are out of range")
    }else if(max(iidx) > object$nsamples & class(object) != "GHap.phase"){
      stop("Some of the provided ids are out of range")
    }
    names(iidx) <- object$id[iidx]
    
    # Organize indices for variants
    if(class(object) == "GHap.haplo"){
      if(max(vidx) > object$nalleles){
        stop("Some of the provided alleles are out of range")
      }else{
        names(vidx) <- paste(object$block[vidx],object$bp1[vidx],
                             object$bp2[vidx],object$allele[vidx], sep="_")
      }
    }else{
      if(max(vidx) > object$nmarkers){
        stop("Some of the provided markers are out of range")
      }else{
        names(vidx) <- object$marker[vidx]
      }
    }
    
    
  }else{
    
    # Organize indices for ids
    iidx <- which(object$id %in% ids)
    if(class(object) == "GHap.phase"){
      if(length(iidx) != 2*length(ids)){
        stop("Some of the provided ids were not found")
      }else{
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
      }
    }else{
      if(length(iidx) != length(ids)){
        stop("Some of the provided ids were not found")
      }else{
        names(iidx) <- object$id[iidx]
        iidx <- iidx[ids]
      }
    }
    
    # Organize indices for variants
    if(class(object) == "GHap.haplo"){
      vidx <- which(1:object$nalleles %in% variants)
      if(length(vidx) != length(alleles)){
        stop("Some of the provided marker names were not found")
      }else{
        vidx <- variants
        names(vidx) <- paste(object$block[vidx],object$bp1[vidx],
                             object$bp2[vidx],object$allele[vidx], sep="_")
      }
    }else{
      vidx <- which(object$marker %in% variants)
      if(length(vidx) != length(variants)){
        stop("Some of the provided marker names were not found")
      }else{
        names(vidx) <- object$marker[vidx]
        vidx <- vidx[variants]
      }
    }
  }
  
  # Define bit function --------------------------------------------------------
  if(class(object) == "GHap.phase"){
    getBitFun <- function(i){
      object.con <- file(object$phase, "rb")
      a <- seek(con = object.con, where = offset*(vidx[i]-1),
                origin = 'start',rw = 'r')
      geno <- readBin(object.con, what=integer(), size = 1,
                      n = offset, signed = FALSE)
      geno <- paste(lookup[geno+1], collapse = "")
      geno <- unlist(strsplit(x = geno, split = ""))
      geno <- as.integer(geno[1:(length(geno)-bitloss)])
      close.connection(object.con)
      return(geno[iidx])
    }
  }else if(class(object) == "GHap.haplo"){
    getBitFun <- function(i){
      object.con <- file(object$genotypes, "rb")
      a <- seek(con = object.con, where = 3 + offset*(vidx[i]-1),
                origin = 'start',rw = 'r')
      geno <- readBin(object.con, what=integer(), size = 1,
                      n = offset, signed = FALSE)
      geno <- paste(lookup[geno+1], collapse = "")
      nc <- nchar(geno)
      n <- seq(1, nc, by = 2)
      geno <- substring(geno, n, c(n[-1]-1, nc))
      geno[which(geno == "00")] <- "0"
      geno[which(geno == "01")] <- "1"
      geno[which(geno == "11")] <- "2"
      geno <- as.integer(geno[1:object$nsamples])
      close.connection(object.con)
      return(geno[iidx])
    }
  }else if(class(object) == "GHap.plink"){
    getBitFun <- function(i){
      object.con <- file(object$plink, "rb")
      a <- seek(con = object.con, where = 3 + offset*(vidx[i]-1),
                origin = 'start',rw = 'r')
      geno <- readBin(object.con, what=integer(), size = 1,
                      n = offset, signed = FALSE)
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
      geno <- as.integer(geno[1:object$nsamples])
      close.connection(object.con)
      return(geno[iidx])
    }
  }
  
  # Get data -------------------------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cl <- makeCluster(ncores)
    clusterExport(cl = cl, varlist = c("object"))
    X <- unlist(parLapply(cl = cl, fun = getBitFun, X = 1:length(vidx)))
    stopCluster(cl)
  }else{
    X <- unlist(mclapply(X = 1:length(vidx), FUN = getBitFun, mc.cores = ncores))
  }
  if(transposed == FALSE){
    X <- Matrix(data = X, nrow = length(vidx), ncol = length(iidx),
                byrow = TRUE, sparse = TRUE)
    colnames(X) <- names(iidx)
    rownames(X) <- names(vidx)
    if(class(object) == "GHap.phase" & unphase == TRUE){
      cols <- 1:ncol(X) %% 2
      X <- X[,which(cols == 1)] + X[,which(cols == 0)]
    }
  }else{
    X <- Matrix(data = X, ncol = length(vidx), nrow = length(iidx),
                byrow = FALSE, sparse = TRUE)
    rownames(X) <- names(iidx)
    colnames(X) <- names(vidx)
    if(class(object) == "GHap.phase" & unphase == TRUE){
      rows <- 1:nrow(X) %% 2
      X <- X[which(rows == 1),] + X[which(rows == 0),]
    }
  }
  if(sparse == FALSE){
    X <- as.matrix(X)
  }
  
  
  # Return X -------------------------------------------------------------------
  return(X)
  
}
