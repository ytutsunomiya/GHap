#Function: ghap.inbcoef
#License: GPLv3 or later
#Modification date: 06 Mar 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute inbreeding coeficients

ghap.inbcoef <- function(
  object,
  freq,
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
  var.n <- object$nmarkers.in
  var.in <- which(object$marker.in)
  
  # Check if vector of allele frequencies is acceptable ------------------------
  if(sum(names(freq) %in% object$marker[var.in]) != length(freq)){
    stop("The vector of allele frequencies contains unknown markers")
  }
  freq <- freq[object$marker[var.in]]
  
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
  inbcoef <- function(x, p){
    usable <- which(p > 0 & p < 1 & is.na(x) == FALSE)
    x <- x[usable]
    p <- p[usable]
    h <- 2*p*(1-p)
    mysum1 <- sum(((x - 2*p)^2/h))
    mysum2 <- sum(x*(2-x)/h)
    mysum3 <- sum((x^2 - (1+2*p)*x + 2*p^2)/h)
    return(c(mysum1,mysum2,mysum3,length(usable)))
  }
  
  #Inbreeding iterate function -------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  sumvariants <- 0
  fhat1 <- rep(x = 0, times = length(id.in))
  fhat2 <- fhat1; fhat3 <- fhat1
  n <- rep(x = 0, times = length(id.in))
  for(i in 1:length(id1)){
    idx <- id1[i]:id2[i]
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = var.in[idx],
                       index = TRUE,
                       unphase = TRUE,
                       impute = FALSE,
                       ncores = ncores)
    tmp <- apply(X = Ztmp, MARGIN = 2, FUN = inbcoef,
                 p = freq[object$marker[var.in[idx]]])
    tmp <- matrix(data = unlist(tmp), ncol = 4, byrow = TRUE)
    fhat1 <- fhat1 + tmp[,1]
    fhat2 <- fhat2 + tmp[,2]
    fhat3 <- fhat3 + tmp[,3]
    n <- n + tmp[,4]
    if(verbose == TRUE){
      sumvariants <- sumvariants + length(idx)
      cat(sumvariants, "variants processed.\r")
    }
  }
  fhat1 <- (fhat1/n)-1
  fhat2 <- 1-(fhat2/n)
  fhat3 <- fhat3/n
  tmp <- object$pop
  names(tmp) <- object$id
  tmp <- tmp[which(duplicated(names(tmp)) == FALSE)]
  out <- data.frame(POP = tmp[colnames(Ztmp)], ID = colnames(Ztmp), N = n,
                    fhat1, fhat2, fhat3, stringsAsFactors = FALSE)
  colnames(out)[4:6] <- c("Fhat1","Fhat2","Fhat3")
  
  #Return output --------------------------------------------------------------
  return(out)
  
}
