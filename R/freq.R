#Function: ghap.freq
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute marker allele frequencies

ghap.freq <- function(
  object,
  type="maf",
  only.active.samples=TRUE,
  only.active.markers=TRUE,
  ncores=1,
  verbose=TRUE
){
  
  # Check if input is a valid GHap object --------------------------------------
  obtype <- c("GHap.phase","GHap.plink")
  if(inherits(object, obtype) == FALSE){
    stop("\nInput must be a valid GHap object.")
  }
  fac <- c(2,1)
  names(fac) <- obtype
  
  # Check if inactive markers and samples should be reactived ------------------
  if(only.active.markers == FALSE){
    object$marker.in <- rep(TRUE,times=object$nmarkers)
    object$nmarkers.in <- length(which(object$marker.in))
  }
  if(only.active.samples == FALSE){
    object$id.in <- rep(TRUE,times=fac[class(object)]*object$nsamples)
    object$nsamples.in <- length(which(object$id.in))/fac[class(object)]
  }
  
  # Check if type is valid -----------------------------------------------------
  if(type %in% c("maf","A0","A1") == FALSE){
    stop("\nArgument type must be 'maf', 'A0' or 'A1'.")
  }
  
  # Calculate offset and bitloss -----------------------------------------------
  offset <- ceiling((2*object$nsamples)/8)
  bitloss <- 8 - ((2*object$nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  lookup <- rep(NA, times=offset*8)
  for(i in 1:offset){
    idx1 <- i*8
    idx2 <- idx1-7
    lookup[idx1:idx2] <- idx2:idx1
  }
  
  # Frequency function ---------------------------------------------------------
  if(inherits(object, "GHap.plink")){
    freqFun <- function(i){
      object.con <- file(object$plink, "rb")
      a <- seek(con = object.con, where = 3 + offset*(vidx[i]-1),
                origin = 'start',rw = 'r')
      geno <- readBin(object.con, what=raw(), size = 1,
                      n = offset, signed = FALSE)
      geno <- as.integer(rawToBits(geno))
      geno1 <- geno[1:length(geno) %% 2 == 1]
      geno2 <- geno[1:length(geno) %% 2 == 0]
      geno <- vector(mode = "integer", length = length(geno)/2)
      geno[which(geno1 == 0 & geno2 == 0)] <- 2
      geno[which(geno1 == 0 & geno2 == 1)] <- 1
      geno[which(geno1 == 1 & geno2 == 1)] <- 0
      geno[which(geno1 == 1 & geno2 == 0)] <- NA
      geno <- geno[1:object$nsamples]
      geno <- geno[iidx]
      geno <- geno[which(is.na(geno) == FALSE)]
      p <- sum(geno)/(2*length(geno))
      close.connection(object.con)
      return(p)
    }
  }else{
    freqFun <- function(i){
      object.con <- file(object$phase, "rb")
      a <- seek(con = object.con, where = offset*(vidx[i]-1),
                origin = 'start',rw = 'r')
      geno <- readBin(object.con, what=raw(), size = 1,
                      n = offset, signed = FALSE)
      geno <- as.integer(rawToBits(geno))
      geno <- geno[lookup]
      geno <- geno[1:(2*object$nsamples)]
      geno <- geno[iidx]
      geno1 <- geno[1:length(geno) %% 2 == 0]
      geno2 <- geno[1:length(geno) %% 2 == 1]
      geno <- geno1 + geno2
      p <- sum(geno)/(2*length(geno))
      close.connection(object.con)
      return(p)
    }
  }
  
  # Frequency calculation ------------------------------------------------------
  iidx <- which(object$id.in)
  vidx <- which(object$marker.in)
  freq <- rep(NA, times=object$nmarkers.in)
  ncores <- min(c(detectCores(), ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cl <- makeCluster(ncores)
    clusterExport(cl = cl, varlist = list("object","iidx","vidx"), envir=environment())
    freq <- unlist(parLapply(cl = cl, fun = freqFun, X = 1:length(vidx)))
    stopCluster(cl)
  }else{
    freq <- unlist(mclapply(X = 1:length(vidx), FUN = freqFun, mc.cores = ncores))
  }
  
  # Results --------------------------------------------------------------------
  names(freq) <- object$marker[vidx]
  if(type == "maf"){
    freq <- pmin(freq,1-freq)
  }else if(type == "A0"){
    freq <- 1-freq
  }
  return(freq)
  
}
