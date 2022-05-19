#Function: ghap.varblup
#License: GPLv3 or later
#Modification date: 19 May 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: convert blup of individuals into blup of variants

ghap.varblup <- function(
  object,
  gebv,
  covmat,
  only.active.variants = TRUE,
  weights = NULL,
  tol = 1e-12,
  vcp = NULL,
  errormat = NULL, 
  errorname = "",
  nlambda = 1000,
  ncores = 1,
  verbose = TRUE
){
  
  # Sanity check for input objects ---------------------------------------------
  obtype <- c("GHap.phase","GHap.plink","GHap.haplo")
  if(class(object) %in% obtype == FALSE){
    stop("\nInput must be a valid GHap object.")
  }
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
  
  # Check if inactive variants should be reactivated ---------------------------
  if(only.active.variants == FALSE){
    if(class(object) == "GHap.haplo"){
      object$allele.in <- rep(TRUE,times=object$nalleles)
      object$nalleles.in <- length(which(object$allele.in))
    }else{
      object$marker.in <- rep(TRUE,times=object$nmarkers)
      object$nmarkers.in <- length(which(object$marker.in))
    }
  }
  if(class(object) == "GHap.haplo"){
    vidx <- which(object$allele.in)
  }else{
    vidx <- which(object$marker.in)
  }
  
  # Check weights --------------------------------------------------------------
  if(is.null(varweights) == TRUE){
    varweights <- rep(1,times=length(vidx))
  }
  if(length(varweights) != length(vidx)){
    stop("\nNumber of variant weights differs from the number of active variants.")
  }
  varweights <- varweights/mean(varweights)
  
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
  
  # Invert reference matrix ----------------------------------------------------
  covmat[names(gebv),names(gebv)]
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
  
  # Compute rotated response ---------------------------------------------------
  k <- icovmat%*%gebv
  
  # Check if errors should be computed -----------------------------------------
  if(is.null(errormat) == FALSE){
    B <- icovmat%*%(vcp*covmat - errormat)%*%icovmat
  }else{
    B <- Diagonal(n = length(gebv))
  }
  
  # Auxiliary functions --------------------------------------------------------
  if(class(object) %in% c("GHap.plink","GHap.haplo")){
    varFun <- function(i){
      object.con <- file(unlist(object[length(object)]), "rb")
      a <- seek(con = object.con, where = 3 + offset*(vidx[i]-1),
                origin = 'start',rw = 'r')
      x <- readBin(object.con, what=raw(), size = 1,
                   n = offset, signed = FALSE)
      x <- as.integer(rawToBits(x))
      x1 <- x[1:length(x) %% 2 == 1]
      x2 <- x[1:length(x) %% 2 == 0]
      x <- vector(mode = "integer", length = length(x)/2)
      x[which(x1 == 0 & x2 == 0)] <- 2
      x[which(x1 == 0 & x2 == 1)] <- 1
      x[which(x1 == 1 & x2 == 1)] <- 0
      x[which(x1 == 1 & x2 == 0)] <- 0
      x <- x[1:object$nsamples]
      names(x) <- object$id
      x <- x[names(gebv)]
      freq <- sum(x)/(2*length(x))
      cent <- mean(x)
      x <- x - cent
      b <- sum(varweights[i]*x*k)
      varxb <- var(x*b)
      varx <- varweights[i]*var(x)
      if(is.null(errormat) == FALSE){
        varb <- (varweights[i]^2)*as.numeric(t(x)%*%B%*%x)
      }else{
        varb <- NA
      }
      close.connection(object.con)
      return(c(freq,b,varxb,cent,varx,varb))
    }
  }else{
    varFun <- function(i){
      object.con <- file(object$phase, "rb")
      a <- seek(con = object.con, where = offset*(vidx[i]-1),
                origin = 'start',rw = 'r')
      x <- readBin(object.con, what=raw(), size = 1,
                   n = offset, signed = FALSE)
      x <- as.integer(rawToBits(x))
      x <- x[lookup]
      x <- x[1:(2*object$nsamples)]
      x1 <- x[1:length(x) %% 2 == 0]
      x2 <- x[1:length(x) %% 2 == 1]
      x <- x1 + x2
      names(x) <- object$id[1:length(object$id) %% 2 == 0]
      x <- x[names(gebv)]
      freq <- sum(x)/(2*length(x))
      cent <- mean(x)
      x <- x - cent
      b <- sum(varweights[i]*x*k)
      varxb <- var(x*b)
      varx <- varweights[i]*var(x)
      if(is.null(errormat) == FALSE){
        varb <- (varweights[i]^2)*as.numeric(t(x)%*%B%*%x)
      }else{
        varb <- NA
      }
      close.connection(object.con)
      return(c(freq,b,varxb,cent,varx,varb))
    }
  }
  
  # Association analysis -------------------------------------------------------
  if(is.null(errormat) == FALSE){
    if(verbose == TRUE){
      cat("Calculation of standard errors and test statistics activated.\n")
    }
  }
  if(verbose == TRUE){
    cat("Converting breeding values into variant effects... ")
  }
  if(Sys.info()["sysname"] == "Windows"){
    cl <- makeCluster(ncores)
    mylist <- list("object","vidx","k","varweights","B")
    clusterExport(cl = cl, varlist = mylist, envir=environment())
    vareff <- unlist(parLapply(cl = cl, fun = varFun, X = 1:length(vidx)))
    stopCluster(cl)
  }else{
    vareff <- unlist(mclapply(X = 1:length(vidx), FUN = varFun, mc.cores = ncores))
  }
  if(verbose == TRUE){
    cat("Done.\n")
  }
  vareff <- matrix(data = vareff, ncol = 6, byrow = T)
  if(class(object) == "GHap.haplo"){
    results <- matrix(data = NA, nrow = length(vidx), ncol = 17)
    results <- as.data.frame(results)
    colnames(results) <- c("CHR","BLOCK","BP1","BP2","ALLELE","FREQ",
                           "SCORE","VAR","pVAR","CENTER","SCALE",
                           "SE","CHISQ.EXP","CHISQ.OBS", "CHISQ.GC",
                           "LOGP","LOGP.GC")
    results$CHR <- object$chr[vidx]
    results$BLOCK <- object$block[vidx]
    results$BP1 <- object$bp1[vidx]
    results$BP2 <- object$bp2[vidx]
    results$ALLELE <- object$allele[vidx]
  }else{
    results <- matrix(data = NA, nrow = length(vidx), ncol = 16)
    results <- as.data.frame(results)
    colnames(results) <- c("CHR","MARKER","BP","ALLELE","FREQ",
                           "SCORE","VAR","pVAR","CENTER","SCALE",
                           "SE","CHISQ.EXP","CHISQ.OBS", "CHISQ.GC",
                           "LOGP","LOGP.GC")
    results$CHR <- object$chr[vidx]
    results$MARKER <- object$marker[vidx]
    results$BP <- object$bp[vidx]
    results$ALLELE <- object$A1[vidx]
  }
  sumvar <- sum(vareff[,5])
  results$FREQ <- vareff[,1]
  results$SCORE <- vareff[,2]/sumvar
  results$VAR <- vareff[,3]*(1/sumvar)^2
  results$pVAR <- results$VAR/sum(results$VAR)
  results$CENTER <- vareff[,4]
  results$SCALE <- 1
  if(is.null(errormat) == FALSE){
    results$SE <- sqrt(vareff[,6]*(1/sumvar)^2)
    results$CHISQ.OBS <- (results$SCORE/results$SE)^2
    results$LOGP <- -1*pchisq(q = results$CHISQ.OBS, df = 1,
                              lower.tail = FALSE, log.p = TRUE)/log(10)
    results$CHISQ.EXP <- qchisq(p = rank(results$CHISQ.OBS)/(nrow(results)+1),
                                df = 1)
    chisq.mean <- mean(results$CHISQ.OBS)
    chisq.dev <- sd(results$CHISQ.OBS)
    chisq.sub <- which(results$CHISQ.OBS < chisq.mean + 3*chisq.dev)
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