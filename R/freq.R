#Function: ghap.freq
#License: GPLv3 or later
#Modification date: 01 Oct 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute marker allele frequencies

ghap.freq <- function(
  object,
  type="maf",
  only.active.samples=TRUE,
  only.active.markers=TRUE,
  batchsize=NULL,
  ncores=1,
  verbose=TRUE
){
  
  # Check if input is a valid GHap object --------------------------------------
  obtype <- c("GHap.phase","GHap.plink")
  if(class(object) %in% obtype == FALSE){
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
  
  # Initialize lookup table ----------------------------------------------------
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
  
  # Generate batch index -------------------------------------------------------
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(object$nmarkers.in/10)
  }
  if(batchsize > object$nmarkers.in){
    batchsize <- object$nmarkers.in
  }
  id1 <- seq(1,object$nmarkers.in,by=batchsize)
  id2 <- (id1+batchsize)-1
  id1 <- id1[id2<=object$nmarkers.in]
  id2 <- id2[id2<=object$nmarkers.in]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,object$nmarkers.in)
  if(id1[length(id1)] > object$nmarkers.in){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  # Frequency calculation ------------------------------------------------------
  ids.in <- which(object$id.in)
  snps.in <- which(object$marker.in)
  freq <- rep(NA, times=object$nmarkers.in)
  ncores <- min(c(detectCores(), ncores))
  for(i in 1:length(id1)){
    X <- ghap.slice(object = object,
                    ids = ids.in,
                    variants = snps.in[id1[i]:id2[i]],
                    index = TRUE,
                    unphase = TRUE,
                    impute = FALSE,
                    lookup = lookup,
                    ncores = ncores)
    if(class(object) == "GHap.plink"){
      n <- apply(X = X, MARGIN = 1,
                 FUN = function(x){length(which(is.na(x) == FALSE))})
      freq[id1[i]:id2[i]] <- rowSums(X, na.rm = TRUE)/(2*n)
    }else{
      freq[id1[i]:id2[i]] <- rowSums(X)/(2*ncol(X))
    }
  }
  
  # Results --------------------------------------------------------------------
  names(freq) <- object$marker[snps.in]
  if(type == "maf"){
    freq <- pmin(freq,1-freq)
  }else if(type == "A0"){
    freq <- 1-freq
  }
  return(freq)
  
}
