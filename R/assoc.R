#Function: ghap.assoc
#License: GPLv3 or later
#Modification date: 24 May 2023
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: phenotype-genotype association analysis

ghap.assoc <- function(
    object,
    formula,
    data,
    covmat,
    batchsize = NULL,
    ngamma = 100,
    nlambda = 1000,
    recalibrate = 0.01,
    only.active.variants=TRUE,
    tol = 1e-12,
    ncores=1,
    verbose=TRUE,
    ...
){
  
  # Check if input is a valid GHap object --------------------------------------
  obtype <- c("GHap.phase","GHap.plink","GHap.haplo")
  if(inherits(object, obtype) == FALSE){
    stop("\nInput must be a valid GHap object.")
  }
  fac <- c(2,1,1)
  names(fac) <- obtype
  
  # Check if inactive variants should be reactived -----------------------------
  if(only.active.variants == FALSE){
    if(inherits(object, "GHap.haplo")){
      object$allele.in <- rep(TRUE,times=object$nalleles)
      object$nalleles.in <- length(which(object$allele.in))
    }else{
      object$marker.in <- rep(TRUE,times=object$nmarkers)
      object$nmarkers.in <- length(which(object$marker.in))
    }
  }
  
  # Map variants ---------------------------------------------------------------
  if(inherits(object, "GHap.haplo")){
    var.n <- object$nalleles.in
    var.in <- which(object$allele.in)
  }else{
    var.n <- object$nmarkers.in
    var.in <- which(object$marker.in)
  }
  
  # Fit mixed model ------------------------------------------------------------
  model <- ghap.lmm(formula = formula, data = data, covmat = covmat,
                    verbose = verbose, extras = "V", errors = FALSE, ...)
  y <- model$residuals$Fixed
  names(y) <- data[,names(model$random[1])]
  Vi <- try(solve(model$extras$V), silent = TRUE)
  if(inherits(Vi, "try-error")){
    Vi <- try(solve(model$extras$V + Diagonal(n = length(y))*tol), silent = TRUE)
    if(inherits(Vi, "try-error")){
      emsg <- paste0("\nUnable to invert phenotypic (co)variance matrix",
                     " even after adding a tolerance of ", tol)
      stop(emsg)
    }
  }
  rm(model)
  
  # Map individuals ------------------------------------------------------------
  id.in <- which(object$id %in% names(y))
  id.n <- length(unique(names(y)))
  
  # Compute rotated response ---------------------------------------------------
  k <- as.numeric(Vi%*%y)
  
  # Auxiliary functions --------------------------------------------------------
  assocFun <- function(i){
    x <- Ztmp[i,unique(names(y))]
    ridx <- which(is.na(x) == FALSE)
    x <- x[ridx]
    freq <- sum(x)/(2*length(x))
    n <- length(x)
    x <- Ztmp[i,names(y)]
    ridx <- which(is.na(x) == FALSE)
    x <- x[ridx]
    x <- x - mean(x)
    varb <- as.numeric(1/(t(x)%*%Vi[ridx,ridx]%*%x))
    b <- varb*sum(x*k[ridx])
    return(c(length(x),n,freq,b,sqrt(varb)))
  }
  gammaFun <- function(i){
    x <- Ztmp[i,names(y)]
    ridx <- which(is.na(x) == FALSE)
    x <- x[ridx]
    x <- x - mean(x)
    g <- (t(x)%*%Vi[ridx,ridx]%*%x)/sum(x^2)
    return(as.numeric(g))
  }
  assocgammaFun <- function(i){
    x <- Ztmp[i,unique(names(y))]
    ridx <- which(is.na(x) == FALSE)
    x <- x[ridx]
    freq <- sum(x)/(2*length(x))
    n <- length(x)
    x <- Ztmp[i,names(y)]
    ridx <- which(is.na(x) == FALSE)
    x <- x[ridx]
    x <- x - mean(x)
    varb <- 1/sum(x^2)
    b <- varb*sum(x*k[ridx])
    return(return(c(length(x),n,freq,b/gamma,sqrt(varb/gamma))))
  }
  
  # Gamma factor calculation ---------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  if(ngamma > 0){
    ranvars <- sample(x = 1:var.n, size = ngamma)
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = var.in[ranvars],
                       index = TRUE,
                       unphase = TRUE,
                       impute = TRUE)
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      clusterExport(cl = cl, varlist = list("Ztmp","Vi","y"),
                    envir=environment())
      gamma <- unlist(parLapply(cl = cl, fun = gammaFun, X = 1:ngamma))
      stopCluster(cl)
    }else{
      gamma <- unlist(mclapply(X = 1:ngamma, FUN = gammaFun, mc.cores = ncores))
    }
    gamma <- gamma[which(is.na(gamma) == FALSE & is.nan(gamma) == FALSE)]
    if(verbose == TRUE){
      cat("Gamma factor estimated using ", ngamma,
          " variants:\n   mean = ", mean(gamma),
          "\n     sd = ", sd(gamma), ".\n", sep = "")
      if(length(gamma) < ngamma){
        cat(ngamma - length(gamma)," monomorphic variants ignored.\n", sep="")
      }
    }
    gamma <- mean(gamma)
  }else{
    if(verbose == TRUE){
      cat("Gamma approximation deactivated.\n")
    }
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
  
  # Batch iterate function -----------------------------------------------------
  sumvariants <- 0
  vareff <- rep(NA, times = 5*var.n)
  ncores <- min(c(detectCores(), ncores))
  if(verbose == TRUE){
    cat("Performing phenotype-genotype association analysis...\n")
  }
  for(i in 1:length(id1)){
    idx <- id1[i]:id2[i]
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = var.in[idx],
                       index = TRUE,
                       unphase = TRUE,
                       impute = TRUE)
    ii <- ((id1[i]-1)*5) + 1
    fi <- ((id2[i]-1)*5) + 5
    if(ngamma > 0){
      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(ncores)
        clusterExport(cl = cl, varlist = list("Ztmp","y","k","gamma"),
                      envir=environment())
        vareff[ii:fi] <- unlist(parLapply(cl = cl, fun = assocgammaFun, X = 1:length(idx)))
        stopCluster(cl)
      }else{
        vareff[ii:fi] <- unlist(mclapply(X = 1:length(idx), FUN = assocgammaFun, mc.cores = ncores))
      }
    }else{
      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(ncores)
        clusterExport(cl = cl, varlist =  list("Ztmp","y","k","Vi"),
                      envir=environment())
        vareff[ii:fi] <- unlist(parLapply(cl = cl, fun = assocFun, X = 1:length(idx)))
        stopCluster(cl)
      }else{
        vareff[ii:fi] <- unlist(mclapply(X = 1:length(idx), FUN = assocFun, mc.cores = ncores))
      }
    }
    if(verbose == TRUE){
      sumvariants <- sumvariants + length(idx)
      cat(sumvariants, "variants processed.\r")
    }
  }
  
  # Build output ---------------------------------------------------------------
  if(verbose == TRUE){
    cat("Done! All variantes processed.\n")
  }
  vareff <- matrix(data = vareff, ncol = 5, byrow = T)
  if(inherits(object, "GHap.haplo")){
    results <- matrix(data = NA, nrow = var.n, ncol = 15)
    results <- as.data.frame(results)
    colnames(results) <- c("CHR","BLOCK","BP1","BP2","ALLELE","FREQ",
                           "N.PHENO","N.GENO","BETA","SE","CHISQ.EXP",
                           "CHISQ.OBS", "CHISQ.GC","LOGP","LOGP.GC")
    results$CHR <- object$chr[var.in]
    results$BLOCK <- object$block[var.in]
    results$BP1 <- object$bp1[var.in]
    results$BP2 <- object$bp2[var.in]
    results$ALLELE <- object$allele[var.in]
  }else{
    results <- matrix(data = NA, nrow = var.n, ncol = 14)
    results <- as.data.frame(results)
    colnames(results) <- c("CHR","MARKER","BP","ALLELE","FREQ",
                           "N.PHENO","N.GENO","BETA","SE","CHISQ.EXP",
                           "CHISQ.OBS", "CHISQ.GC","LOGP","LOGP.GC")
    results$CHR <- object$chr[var.in]
    results$MARKER <- object$marker[var.in]
    results$BP <- object$bp[var.in]
    results$ALLELE <- object$A1[var.in]
  }
  results$N.PHENO <- vareff[,1]
  results$N.GENO <- vareff[,2]
  results$FREQ <- vareff[,3]
  results$BETA <- vareff[,4]
  results$SE <- vareff[,5]
  results$CHISQ.OBS <- (results$BETA/results$SE)^2
  results$LOGP <- -1*pchisq(q = results$CHISQ.OBS, df = 1,
                            lower.tail = FALSE, log.p = TRUE)/log(10)
  
  # Check if recalibration is required -----------------------------------------
  if(ngamma > 0  & recalibrate > 0){
    ntop <- ceiling(recalibrate*nrow(results))
    top <- order(results$LOGP, decreasing = TRUE)[1:ntop]
    if(verbose == TRUE){
      cat("Recalibrating statistics for the top ", recalibrate*100,
          "% (", ntop, ") variants... ", sep = "")
    }
    if(inherits(object, "GHap.haplo")){
      vidx <- 1:object$nalleles
      names(vidx) <- paste0(object$chr, object$block, object$bp1,
                            object$bp2, object$allele)
      vidx <- vidx[paste0(results$CHR,results$BLOCK,results$BP1,
                          results$BP2,results$ALLELE)[top]]
    }else{
      vidx <- 1:object$nmarkers
      names(vidx) <- object$marker
      vidx <- vidx[results$MARKER[top]]
    }
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = vidx,
                       index = TRUE,
                       unphase = TRUE,
                       impute = TRUE)
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      clusterExport(cl = cl, varlist = list("Ztmp","y","k","Vi"),
                    envir=environment())
      vareff <- unlist(parLapply(cl = cl, fun = assocFun, X = 1:ntop))
      stopCluster(cl)
    }else{
      vareff <- unlist(mclapply(X = 1:ntop, FUN = assocFun, mc.cores = ncores))
    }
    vareff <- matrix(data = vareff, ncol = 5, byrow = T)
    results$BETA[top] <- vareff[,4]
    results$SE[top] <- vareff[,5]
    results$CHISQ.OBS[top] <- (results$BETA[top]/results$SE[top])^2
    results$LOGP[top] <- -1*pchisq(q = results$CHISQ.OBS[top], df = 1,
                                   lower.tail = FALSE, log.p = TRUE)/log(10)
    if(verbose == TRUE){
      cat("Done.\n")
    }
  }
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
  
  
  #Return output --------------------------------------------------------------
  return(results)
  
}
