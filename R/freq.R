#Function: ghap.freq
#License: GPLv3 or later
#Modification date: 23 May 2023
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute marker allele frequencies

ghap.freq <- function(
  object,
  type="maf",
  batchsize=NULL,
  only.active.samples=TRUE,
  only.active.markers=TRUE,
  verbose=TRUE
){
  
  # Check if input is a valid GHap object --------------------------------------
  obtype <- c("GHap.phase","GHap.plink","GHap.haplo")
  if(inherits(object, obtype) == FALSE){
    stop("\nInput must be a valid GHap object.")
  }
  fac <- c(2,1,1)
  names(fac) <- obtype
  
  # Check if inactive variants and samples should be reactived -----------------
  if(only.active.variants == FALSE){
    if(inherits(object, "GHap.haplo")){
      object$allele.in <- rep(TRUE,times=object$nalleles)
      object$nalleles.in <- length(which(object$allele.in))
    }else{
      object$marker.in <- rep(TRUE,times=object$nmarkers)
      object$nmarkers.in <- length(which(object$marker.in))
    }
  }
  if(only.active.samples == FALSE){
    object$id.in <- rep(TRUE,times=fac[class(object)]*object$nsamples)
    object$nsamples.in <- length(which(object$id.in))/fac[class(object)]
  }
  
  # Map individuals and variants -----------------------------------------------
  id.n <- object$nsamples.in
  id.in <- which(object$id.in)
  if(inherits(object, "GHap.haplo")){
    var.n <- object$nalleles.in
    var.in <- which(object$allele.in)
  }else{
    var.n <- object$nmarkers.in
    var.in <- which(object$marker.in)
  }
  
  # Check if type is valid -----------------------------------------------------
  if(type %in% c("maf","A0","A1") == FALSE){
    stop("\nArgument type must be 'maf', 'A0' or 'A1'.")
  }
  
  # Generate batch index -------------------------------------------------------
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(var.n/10)
  }
  if(batchsize > var.n){
    batchsize <- var.n
  }
  id1 <- seq(1,var.n,by=batchsize)
  id2 <- (id1+batchsize)-1
  id1 <- id1[id2<=var.n]
  id2 <- id2[id2<=var.n]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,var.n)
  if(id1[length(id1)] > var.n){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  #Log message -----------------------------------------------------------------
  if(verbose == TRUE){
    cat("Processing ", var.n, " variants in ", length(id1), " batches.\n", sep="")
    cat("Inactive variants will be ignored.\n")
  }
  
  #Initialize allele frequency function ----------------------------------------
  freqfun <- function(x){
    n <- length(which(is.na(x) == FALSE))
    p <- sum(x, na.rm = TRUE)/(2*n)
    return(p)
  }
  
  # Frequency iterate function -------------------------------------------------
  sumvariants <- 0
  freq <- rep(NA, times = var.n)
  names(freq) <- object$marker[var.in]
  for(i in 1:length(id1)){
    idx <- id1[i]:id2[i]
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = var.in[idx],
                       index = TRUE,
                       unphase = TRUE,
                       impute = TRUE)
    zids <- colnames(Ztmp)
    freq[idx] <- apply(X = Ztmp, MARGIN = 1, FUN = freqfun)
    if(verbose == TRUE){
      sumvariants <- sumvariants + length(idx)
      cat(sumvariants, "variants processed.\r")
    }
  }
  
  # Results --------------------------------------------------------------------
  if(type == "maf"){
    freq <- pmin(freq,1-freq)
  }else if(type == "A0"){
    freq <- 1-freq
  }
  return(freq)
  
}
