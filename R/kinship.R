#Function: ghap.kinship
#License: GPLv3 or later
#Modification date: 23 Sep 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute relationship matrix

ghap.kinship <- function(
  object,
  weights=NULL,
  sparsity=NULL,
  batchsize=NULL,
  only.active.samples=TRUE,
  only.active.variants=TRUE,
  ncores=1,
  verbose=TRUE
){
  
  # Check if input is a valid GHap object --------------------------------------
  obtype <- c("GHap.phase","GHap.plink","GHap.haplo")
  if(class(object) %in% obtype == FALSE){
    stop("\nInput must be a valid GHap object.")
  }
  fac <- c(2,1,1)
  names(fac) <- obtype
  
  # Check if inactive variants and samples should be reactived -----------------
  if(only.active.variants == FALSE){
    if(class(object) == "GHap.haplo"){
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
  
  #Check weights
  if (is.null(weights) == TRUE) {
    weights <- rep(1,times=var.n)
  }
  if (length(weights) != var.n) {
    stop("Vector of weights must have the same length as the number of variants.")
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
  
  #Initialize kinship matrix ---------------------------------------------------
  if(verbose == TRUE){
    cat("Preparing", id.n, "x", id.n, "kinship matrix.\n")
  }
  K <- Matrix(data = 0, nrow = id.n, ncol = id.n, doDiag = F)
  K <- as(as(K,"dsyMatrix"),"dspMatrix")
  
  #Kinship iterate function ----------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  sumvariants <- 0
  for(i in 1:length(id1)){
    idx <- id1[i]:id2[i]
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = var.in[idx],
                       index = TRUE,
                       unphase = TRUE,
                       impute = TRUE,
                       lookup = lookup,
                       ncores = ncores)
    Ztmp.mean <- apply(X = Ztmp, MARGIN = 1, FUN = mean)
    Ztmp <- (Ztmp - Ztmp.mean)
    K <- K + crossprod(Ztmp*sqrt(weights[idx]))
    if(verbose == TRUE){
      sumvariants <- sumvariants + length(idx)
      cat(sumvariants, "variants processed.\r")
    }
  }
  
  #Scale kinship matrix -------------------------------------------------------
  q <- mean(diag(K))
  K <- K/q
  colnames(K) <- object$id[id.in]
  rownames(K) <- colnames(K)
  
  #Induce sparsity ------------------------------------------------------------
  if(is.null(sparsity) == FALSE){
    K[K < sparsity] <- 0
    K <- as(as(K,"dsyMatrix"),"dsCMatrix")
  }
  
  #Return output --------------------------------------------------------------
  return(K)
  
}
