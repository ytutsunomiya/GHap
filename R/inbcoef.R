#Function: ghap.inbcoef
#License: GPLv3 or later
#Modification date: 06 Mar 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute inbreeding coeficients

ghap.inbcoef <- function(
  object,
  freq,
  type=1,
  batchsize=NULL,
  only.active.samples=TRUE,
  only.active.variants=TRUE,
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
  
  # Check if inactive variants and samples should be reactived -----------------
  if(only.active.variants == FALSE){
    object$marker.in <- rep(TRUE,times=object$nmarkers)
    object$nmarkers.in <- length(which(object$marker.in))
  }
  if(only.active.samples == FALSE){
    object$id.in <- rep(TRUE,times=fac[class(object)]*object$nsamples)
    object$nsamples.in <- length(which(object$id.in))/fac[class(object)]
  }
  
  # Map number of variants -----------------------------------------------------
  id.n <- object$nsamples.in
  id.in <- which(object$id.in)
  if(class(object) == "GHap.haplo"){
    var.n <- object$nalleles.in
    var.in <- which(object$allele.in)
  }else{
    var.n <- object$nmarkers.in
    var.in <- which(object$marker.in)
  }
  
  # Check if vector of allele frequencies is acceptable ------------------------
  if(sum(names(freq) %in% object$marker[var.in]) != length(freq)){
    stop("The vector of allele frequencies contains unknown markers")
  }
  freq <- freq[object$marker[var.in]]
  freq[which(freq == 0)] <- 1e-12

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
  
  #Initialize inbreeding functions ---------------------------------------------
  inbcoef <- vector(mode = "list", length = 3)
  inbcoef[[1]] <- function(x, p){
    mysum <- sum((x^2 - (1+2*p)*x + 2*p^2)/(2*p*(1-p)))
    return(mysum)
  }
  inbcoef[[2]] <- function(x, p){
    mysum <- sum(((x - 2*p)^2)/(2*p*(1-p)))
    return(mysum)
  }
  inbcoef[[3]] <- function(x, p){
    mysum <- sum((x*(2-x))/(2*p*(1-p)))
    return(mysum)
  }
  s <- c(1,-1,1)
  k <- c(0,1,-1)
  
  #Inbreeding iterate function -------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  sumvariants <- 0
  f <- rep(x = 0, times = length(id.in))
  for(i in 1:length(id1)){
    idx <- id1[i]:id2[i]
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = var.in[idx],
                       index = TRUE,
                       unphase = TRUE,
                       impute = TRUE,
                       ncores = ncores)
    tmp <- apply(X = Ztmp, MARGIN = 2, FUN = inbcoef[[type]],
                 p = freq[object$marker[var.in[idx]]])
    f <- f + tmp
    if(verbose == TRUE){
      sumvariants <- sumvariants + length(idx)
      cat(sumvariants, "variants processed.\r")
    }
  }
  f <- k[type] + s[type]*f/var.n
  if(type == 1){
    f[which(f < 0)] <- 0
    f[which(f > 1)] <- 1
  }
  tmp <- object$pop
  names(tmp) <- object$id
  tmp <- tmp[which(duplicated(names(tmp)) == FALSE)]
  out <- data.frame(POP = tmp[names(f)], ID = names(f), C = f)
  mycol <- c("UNI","HOM","GRM")
  colnames(out)[3] <- paste0("F",mycol[type])
  
  #Return output --------------------------------------------------------------
  return(out)
  
}