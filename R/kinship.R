#Function: ghap.kinship
#License: GPLv3 or later
#Modification date: 16 Dec 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute relationship matrix

ghap.kinship <- function(
    object,
    weights=NULL,
    sparsity=NULL,
    type=1,
    batchsize=NULL,
    freq=NULL,
    idrow=NULL,
    idcol=NULL,
    only.active.samples=TRUE,
    only.active.variants=TRUE,
    ncores=1,
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
  
  # Map individuals ------------------------------------------------------------
  if(is.null(idrow) & is.null(idcol)){
    id.n <- object$nsamples.in
    id.in <- which(object$id.in)
  }else if(is.null(idrow) == FALSE & is.null(idcol) == FALSE){
    ids <- c(idrow,idcol)
    id.n <- length(which(ids %in% object$id))
    if(id.n != length(ids)){
      stop("\nSome of the provided ids were not found")
    }
    id.in <- which(object$id %in% ids)
    if(is.null(freq) | type != 3){
      emsg <- paste0("\nSpecific ids for rows and columns are currently ",
                     "allowed only if freq is provided with type = 3")
      stop(emsg)
    }
  }else{
    stop("\nIDs must be provided for both rows and columns")
  }
  
  # Map number of variants -----------------------------------------------------
  if(inherits(object, "GHap.haplo")){
    var.n <- object$nalleles.in
    var.in <- which(object$allele.in)
  }else{
    var.n <- object$nmarkers.in
    var.in <- which(object$marker.in)
  }
  
  # Check frequency vector -----------------------------------------------------
  if(is.null(freq) == FALSE){
    freq <- freq[which(names(freq) %in% object$marker[var.in])]
    if(length(freq) != length(var.in)){
      stop("\nFrequency vector must include all active markers.")
    }else{
      freq <- freq[object$marker[var.in]]
    }
  }
  
  #Check weights
  if (is.null(weights) == FALSE & length(weights) != var.n) {
    stop("Vector of weights must have the same length as the number of variants.")
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
  
  #Initialize scaling function -------------------------------------------------
  scalefun <- vector(mode = "list", length = 6)
  scalefun[[1]] <- function(x){
    m <- mean(x, na.rm = TRUE)
    aa <- which(x == 0)
    ab <- which(x == 1)
    bb <- which(x == 2)
    x[aa] <- -m
    x[ab] <- 1 - m
    x[bb] <- 2 - m
    return(x)
  }
  scalefun[[2]] <- function(x){
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    aa <- which(x == 0)
    ab <- which(x == 1)
    bb <- which(x == 2)
    x[aa] <- -m/s
    x[ab] <- (1-m)/s
    x[bb] <- (2-m)/s
    return(x)
  }
  scalefun[[3]] <- function(x){
    n <- length(which(is.na(x) == FALSE))
    p <- sum(x, na.rm = TRUE)/(2*n)
    m <- 2*p
    aa <- which(x == 0)
    ab <- which(x == 1)
    bb <- which(x == 2)
    x[aa] <- -m
    x[ab] <- 1 - m
    x[bb] <- 2 - m
    return(x)
  }
  scalefun[[4]] <- function(x){
    n <- length(which(is.na(x) == FALSE))
    p <- sum(x, na.rm = TRUE)/(2*n)
    s <- sqrt(2*p*(1-p))
    m <- 2*p
    aa <- which(x == 0)
    ab <- which(x == 1)
    bb <- which(x == 2)
    x[aa] <- -m/s
    x[ab] <- (1-m)/s
    x[bb] <- (2-m)/s
    return(x)
  }
  scalefun[[5]] <- function(x){
    n <- length(which(is.na(x) == FALSE))
    p <- sum(x, na.rm = TRUE)/(2*n)
    aa <- which(x == 0)
    ab <- which(x == 1)
    bb <- which(x == 2)
    x[aa] <- -2*p^2
    x[ab] <- 2*p*(1-p)
    x[bb] <- -2*(1-p)^2
    return(x)
  }
  scalefun[[6]] <- scalefun[[5]]
  
  #Initialize denominators -----------------------------------------------------
  scaleval <- vector(mode = "list", length = 6)
  scaleval[[1]] <- function(){return(mean(diag(K)))}
  scaleval[[2]] <- function(){return(var.n)}
  scaleval[[3]] <- function(){return(2*sum(p*(1-p)))}
  scaleval[[4]] <- scaleval[[2]]
  scaleval[[5]] <- function(){return(4*sum(p^2*(1-p)^2))}
  scaleval[[6]] <- scaleval[[1]]
  
  #Initialize kinship matrix ---------------------------------------------------
  if(verbose == TRUE){
    if(is.null(idrow) & is.null(idcol)){
      cat("Preparing", id.n, "x", id.n, "kinship matrix.\n")
    }else{
      cat("Preparing", length(idrow), "x", length(idcol), "kinship matrix.\n")
    }
  }
  if(is.null(idrow) & is.null(idcol)){
    K <- Matrix(data = 0, nrow = id.n, ncol = id.n, doDiag = F)
  }else{
    K <- Matrix(data = 0, nrow = length(idrow), ncol = length(idcol), doDiag = F)
  }
  
  #Kinship iterate function ----------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  sumvariants <- 0
  q <- 0
  for(i in 1:length(id1)){
    idx <- id1[i]:id2[i]
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = var.in[idx],
                       index = TRUE,
                       unphase = TRUE,
                       impute = TRUE)
    zids <- colnames(Ztmp)
    if(is.null(freq)){
      p <- apply(X = Ztmp, MARGIN = 1, FUN = freqfun)
    }else{
      p <- freq[object$marker[var.in[idx]]]
    }
    exc <- which(pmin(p,1-p) == 0)
    if(length(exc) > 0){
      p <- p[-exc]
      Ztmp <- Ztmp[-exc,]
      var.n <- var.n - length(exc)
    }
    if(type %in% c(3,5)){
      q <- q + scaleval[[type]]() 
    }
    if(type == 3 & is.null(freq) == FALSE){
      Ztmp <- t(Ztmp - 2*p)
    }else{
      Ztmp <- apply(X = Ztmp, MARGIN = 1, FUN = scalefun[[type]])
    }
    if(is.null(weights) == FALSE){
      Ztmp*sqrt(weights[idx])
    }
    if(is.null(idrow) & is.null(idcol)){
      K <- K + tcrossprod(Ztmp)
    }else{
      K <- K + tcrossprod(Ztmp[idrow,],Ztmp[idcol,])
    }
    if(verbose == TRUE){
      sumvariants <- sumvariants + length(idx)
      cat(sumvariants, "variants processed.\r")
    }
  }
  
  #Scale kinship matrix -------------------------------------------------------
  if(type %in% c(3,5) == FALSE){
    q <- scaleval[[type]]() 
  }
  K <- K/q
  if(is.null(idrow) & is.null(idcol)){
    colnames(K) <- zids
    rownames(K) <- colnames(K)
  }else{
    rownames(K) <- idrow
    colnames(K) <- idcol
  }
  
  #Induce sparsity ------------------------------------------------------------
  if(is.null(sparsity) == FALSE){
    K <- drop0(K, tol = sparsity)
  }
  
  #Return output --------------------------------------------------------------
  return(K)
  
}
