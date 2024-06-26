#Function: ghap.varblup
#License: GPLv3 or later
#Modification date: 04 Jun 2024
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: convert blup of individuals into blup of variants

ghap.varblup <- function(
    object,
    gebv,
    covmat,
    type = 1,
    batchsize = NULL,
    only.active.variants = TRUE,
    weights = NULL,
    tol = 1e-12,
    vcp = NULL,
    errormat = NULL, 
    errorname = "",
    invcov = FALSE,
    nlambda = 1000,
    ncores = 1,
    verbose = TRUE
){
  
  # Check if input is a valid GHap object --------------------------------------
  obtype <- c("GHap.phase","GHap.plink","GHap.haplo")
  if(inherits(object, obtype) == FALSE){
    stop("\nInput must be a valid GHap object.")
  }
  fac <- c(2,1,1)
  names(fac) <- obtype
  
  # Sanity check for input objects ---------------------------------------------
  if(is.null(names(gebv))){
    stop("\nArgument 'gebv' must be a named vector.\n")
  }
  if(identical(rownames(covmat),colnames(covmat)) == FALSE){
    stop("\nRow and column names in 'covmat' are not identical.\n")
  }
  if(length(which(names(gebv) %in% colnames(covmat))) != length(gebv)){
    stop("\nSome 'gebv' names are missing in 'covmat'.\n")
  }
  if(is.null(errormat) == FALSE){
    if(identical(rownames(errormat),colnames(errormat)) == FALSE){
      stop("\nRow and column names in 'errormat' are not identical.\n")
    }
    mycol <- grep(pattern = errorname, x = colnames(errormat))
    errormat <- errormat[mycol,mycol]
    colnames(errormat) <- gsub(pattern = errorname,
                               replacement = "", x = colnames(errormat))
    rownames(errormat) <- colnames(errormat)
    if(identical(sort(rownames(errormat)),sort(names(gebv))) == FALSE){
      stop("\nNames in 'errormat' must match the names in 'gebv'.\n")
    }
    errormat <- errormat[names(gebv),names(gebv)]
  }
  if(type %in% 1:4 == FALSE){
    stop("Covariance matrix types currently supported are 1, 2, 3 or 4.")
  }
  
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
  
  # Map individuals and variants -----------------------------------------------
  id.in <- which(object$id %in% names(gebv))
  if(inherits(object, "GHap.haplo")){
    var.n <- object$nalleles.in
    var.in <- which(object$allele.in)
  }else{
    var.n <- object$nmarkers.in
    var.in <- which(object$marker.in)
  }  
  
  # Check weights --------------------------------------------------------------
  if(is.null(weights) == TRUE){
    weights <- rep(1,times=var.n)
  }
  if(length(weights) != var.n){
    stop("\nNumber of variant weights differs from the number of active variants.")
  }
  weights <- weights/mean(weights)
  
  # Invert reference matrix ----------------------------------------------------
  if(invcov == FALSE){
    covmat <- covmat[names(gebv),names(gebv)]
    icovmat <- try(solve(covmat), silent = TRUE)
    if(inherits(icovmat, "try-error")){
      icovmat <- try(solve(covmat + Diagonal(length(gebv))*tol), silent = TRUE)
      if(inherits(icovmat, "try-error")){
        emsg <- paste0("\nUnable to invert the covariance matrix",
                       " even after adding a tolerance of ",
                       tol)
        stop(emsg)
      }
    }
  }else{
    icovmat <- covmat
    gebv <- gebv[rownames(icovmat)]
  }
  
  # Compute rotated response ---------------------------------------------------
  k <- icovmat%*%gebv
  
  # Check if errors should be computed -----------------------------------------
  if(is.null(errormat) == FALSE){
    B <- icovmat%*%(vcp*covmat - errormat)%*%icovmat
  }else{
    B <- Diagonal(n = length(gebv))
  }
  
  # Initialize scaling function ------------------------------------------------
  scalefun <- vector(mode = "list", length = 6)
  scalefun[[1]] <- function(x){
    m <- mean(x)
    s <- sd(x)
    p <- sum(x)/(2*length(x))
    aa <- which(x == 0)
    ab <- which(x == 1)
    bb <- which(x == 2)
    x[aa] <- -m
    x[ab] <- 1 - m
    x[bb] <- 2 - m
    return(c(m,s,p,x))
  }
  scalefun[[2]] <- function(x){
    m <- mean(x)
    s <- sd(x)
    p <- sum(x)/(2*length(x))
    aa <- which(x == 0)
    ab <- which(x == 1)
    bb <- which(x == 2)
    x[aa] <- -m/s
    x[ab] <- (1-m)/s
    x[bb] <- (2-m)/s
    return(c(m,s,p,x))
  }
  scalefun[[3]] <- function(x){
    p <- sum(x)/(2*length(x))
    m <- 2*p
    s <- sqrt(2*p*(1-p))
    aa <- which(x == 0)
    ab <- which(x == 1)
    bb <- which(x == 2)
    x[aa] <- -m
    x[ab] <- 1 - m
    x[bb] <- 2 - m
    return(c(m,s,p,x))
  }
  scalefun[[4]] <- function(x){
    p <- sum(x)/(2*length(x))
    s <- sqrt(2*p*(1-p))
    m <- 2*p
    aa <- which(x == 0)
    ab <- which(x == 1)
    bb <- which(x == 2)
    x[aa] <- -m/s
    x[ab] <- (1-m)/s
    x[bb] <- (2-m)/s
    return(c(m,s,p,x))
  }
  
  # Initialize denominators ----------------------------------------------------
  scaleval <- vector(mode = "list", length = 4)
  scaleval[[1]] <- function(){return(sum(vareff[,5]^2))}
  scaleval[[2]] <- function(){return(length(var.in))}
  scaleval[[3]] <- scaleval[[1]]
  scaleval[[4]] <- scaleval[[2]]
  
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
  
  # Effect solver --------------------------------------------------------------
  varFun <- function(i){
    x <- scalefun[[type]](Ztmp[i,])
    m <- x[1]
    s <- x[2]
    p <- x[3]
    x <- x[-c(1:3)]
    b <- sum(weights[idx[i]]*x*k)
    varxb <- var(x*b)
    varx <- weights[i]*var(x)
    if(is.null(errormat) == FALSE){
      varb <- (weights[i]^2)*as.numeric(t(x)%*%B%*%x)
    }else{
      varb <- NA
    }
    return(c(p,b,varxb,m,s,varx,varb))
  }
  
  # Log message ----------------------------------------------------------------
  if(verbose == TRUE){
    cat("Processing ", var.n, " variants in ", length(id1), " batches.\n", sep="")
    cat("Inactive variants will be ignored.\n")
    if(is.null(errormat) == FALSE){
      if(verbose == TRUE){
        cat("Calculation of standard errors and test statistics activated.\n")
      }
    }
    cat("Converting breeding values into variant effects...\n")
  }
  
  # Batch iterate function -----------------------------------------------------
  sumvariants <- 0
  vareff <- rep(NA, times = 7*var.n)
  ncores <- min(c(detectCores(), ncores))
  for(i in 1:length(id1)){
    idx <- id1[i]:id2[i]
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = var.in[idx],
                       index = TRUE,
                       unphase = TRUE,
                       impute = TRUE)
    Ztmp <- Ztmp[,names(gebv)]
    ii <- ((id1[i]-1)*7) + 1
    fi <- ((id2[i]-1)*7) + 7
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      mylist <- list("Ztmp","idx","k","weights","B","scalefun")
      clusterExport(cl = cl, varlist = mylist, envir=environment())
      vareff <- unlist(parLapply(cl = cl, fun = varFun, X = 1:length(idx)))
      stopCluster(cl)
    }else{
      vareff[ii:fi] <- unlist(mclapply(X = 1:length(idx), FUN = varFun, mc.cores = ncores))
    }
    if(verbose == TRUE){
      sumvariants <- sumvariants + length(idx)
      cat(sumvariants, "variants processed.\r")
    }
  }
  
  # Build output ---------------------------------------------------------------
  if(verbose == TRUE){
    cat("Done! All variants processed.\n")
  }
  vareff <- matrix(data = vareff, ncol = 7, byrow = T)
  if(inherits(object, "GHap.haplo")){
    results <- matrix(data = NA, nrow = length(var.in), ncol = 17)
    results <- as.data.frame(results)
    colnames(results) <- c("CHR","BLOCK","BP1","BP2","ALLELE","FREQ",
                           "SCORE","VAR","pVAR","CENTER","SCALE",
                           "SE","CHISQ.EXP","CHISQ.OBS", "CHISQ.GC",
                           "LOGP","LOGP.GC")
    results$CHR <- object$chr[var.in]
    results$BLOCK <- object$block[var.in]
    results$BP1 <- object$bp1[var.in]
    results$BP2 <- object$bp2[var.in]
    results$ALLELE <- object$allele[var.in]
  }else{
    results <- matrix(data = NA, nrow = length(var.in), ncol = 16)
    results <- as.data.frame(results)
    colnames(results) <- c("CHR","MARKER","BP","ALLELE","FREQ",
                           "SCORE","VAR","pVAR","CENTER","SCALE",
                           "SE","CHISQ.EXP","CHISQ.OBS", "CHISQ.GC",
                           "LOGP","LOGP.GC")
    results$CHR <- object$chr[var.in]
    results$MARKER <- object$marker[var.in]
    results$BP <- object$bp[var.in]
    results$ALLELE <- object$A1[var.in]
  }
  results$FREQ <- vareff[,1]
  sumvar <- scaleval[[type]]()
  results$SCORE <- vareff[,2]/sumvar
  results$VAR <- vareff[,3]*(1/sumvar)^2
  results$pVAR <- results$VAR/sum(results$VAR)
  results$CENTER <- vareff[,4]
  if(type %in% c(1,3)){
    results$SCALE <- 1
  }else{
    results$SCALE <- vareff[,5]    
  }
  if(is.null(errormat) == FALSE){
    results$SE <- sqrt(vareff[,7]*(1/sumvar)^2)
    results$CHISQ.OBS <- (results$SCORE/results$SE)^2
    results$LOGP <- -1*pchisq(q = results$CHISQ.OBS, df = 1,
                              lower.tail = FALSE, log.p = TRUE)/log(10)
    poly <- which(results$FREQ > 0 & results$FREQ < 1)
    if(verbose == TRUE & length(poly) < nrow(results)){
      cat(nrow(results) - length(poly)," monomorphic variants in results.\n",
          "[NOTE] Subsetting polymorphic variants prior to the analysis is advised.\n", sep="")
    }
    results$CHISQ.EXP <- NA
    results$CHISQ.EXP[poly] <- qchisq(p = rank(results$CHISQ.OBS[poly])/(length(poly)+1), df = 1)
    chisq.mean <- mean(results$CHISQ.OBS, na.rm = TRUE)
    chisq.dev <- sd(results$CHISQ.OBS, na.rm = TRUE)
    chisq.sub <- which(is.na(results$CHISQ.EXP) == FALSE & 
                         results$CHISQ.OBS < chisq.mean + 3*chisq.dev)
    ranvars <- sample(x = chisq.sub, size = nlambda)
    lambda <- lm(formula = CHISQ.OBS ~ CHISQ.EXP, data = results[ranvars,])
    lambda <- as.numeric(lambda$coefficients[2])
    results$CHISQ.GC <- results$CHISQ.OBS/lambda
    results$LOGP.GC <- -1*pchisq(q = results$CHISQ.GC, df = 1,
                                 lower.tail = FALSE, log.p = TRUE)/log(10)
    if(verbose == TRUE){
      cat("Inflation factor estimated using ", nlambda,
          " variants = ", lambda, ".\n", sep = "")
    }
  }
  
  
  #Return output --------------------------------------------------------------
  return(results)
  
}
