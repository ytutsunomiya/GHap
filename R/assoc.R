#Function: ghap.assoc
#License: GPLv3 or later
#Modification date: 5 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: phenotype-genotype association analysis

ghap.assoc <- function(
  object,
  formula,
  data,
  covmat,
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
  
  # Check if inactive variants should be reactivated ---------------------------
  if(only.active.variants == FALSE){
    if(inherits(object, "GHap.haplo")){
      object$allele.in <- rep(TRUE,times=object$nalleles)
      object$nalleles.in <- length(which(object$allele.in))
    }else{
      object$marker.in <- rep(TRUE,times=object$nmarkers)
      object$nmarkers.in <- length(which(object$marker.in))
    }
  }
  if(inherits(object, "GHap.haplo")){
    vidx <- which(object$allele.in)
  }else{
    vidx <- which(object$marker.in)
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
  
  # Compute rotated response ---------------------------------------------------
  k <- as.numeric(Vi%*%y)
  
  # Auxiliary functions --------------------------------------------------------
  if(inherits(object, c("GHap.plink","GHap.haplo"))){
    assocFun <- function(i){
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
      x[which(x1 == 1 & x2 == 0)] <- NA
      x <- x[1:object$nsamples]
      names(x) <- object$id
      x <- x[names(y)]
      ridx <- which(is.na(x) == FALSE)
      x <- x[ridx]
      freq <- sum(x)/(2*length(x))
      x <- x - mean(x)
      varb <- as.numeric(1/(t(x)%*%Vi[ridx,ridx]%*%x))
      b <- varb*sum(x*k[ridx])
      close.connection(object.con)
      return(c(length(x),freq,b,sqrt(varb)))
    }
    gammaFun <- function(i){
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
      x[which(x1 == 1 & x2 == 0)] <- NA
      x <- x[1:object$nsamples]
      names(x) <- object$id
      x <- x[names(y)]
      ridx <- which(is.na(x) == FALSE)
      x <- x[ridx]
      x <- x - mean(x)
      g <- (t(x)%*%Vi[ridx,ridx]%*%x)/sum(x^2)
      close.connection(object.con)
      return(as.numeric(g))
    }
    assocgammaFun <- function(i){
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
      x[which(x1 == 1 & x2 == 0)] <- NA
      x <- x[1:object$nsamples]
      names(x) <- object$id
      x <- x[names(y)]
      ridx <- which(is.na(x) == FALSE)
      x <- x[ridx]
      freq <- sum(x)/(2*length(x))
      x <- x - mean(x)
      varb <- 1/sum(x^2)
      b <- varb*sum(x*k[ridx])
      close.connection(object.con)
      return(return(c(length(x),freq,b/gamma,sqrt(varb/gamma))))
    }
  }else{
    assocFun <- function(i){
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
      x <- x[names(y)]
      ridx <- which(is.na(x) == FALSE)
      x <- x[ridx]
      freq <- sum(x)/(2*length(x))
      x <- x - mean(x)
      varb <- as.numeric(1/(t(x)%*%Vi[ridx,ridx]%*%x))
      b <- varb*sum(x*k[ridx])
      close.connection(object.con)
      return(c(length(x),freq,b,sqrt(varb)))
    }
    gammaFun <- function(i){
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
      x <- x[names(y)]
      ridx <- which(is.na(x) == FALSE)
      x <- x[ridx]
      x <- x - mean(x)
      g <- (t(x)%*%Vi[ridx,ridx]%*%x)/sum(x^2)
      close.connection(object.con)
      return(as.numeric(g))
    }
    assocgammaFun <- function(i){
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
      x <- x[names(y)]
      ridx <- which(is.na(x) == FALSE)
      x <- x[ridx]
      freq <- sum(x)/(2*length(x))
      x <- x - mean(x)
      varb <- 1/sum(x^2)
      b <- varb*sum(x*k[ridx])
      close.connection(object.con)
      return(return(c(length(x),freq,b/gamma,sqrt(varb/gamma))))
    }
  }
  
  # Gamma factor calculation ---------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  if(ngamma > 0){
    ranvars <- sample(x = 1:length(vidx), size = ngamma)
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      clusterExport(cl = cl, varlist = list("object","vidx","k","covmat"),
                    envir=environment())
      gamma <- unlist(parLapply(cl = cl, fun = gammaFun, X = ranvars))
      stopCluster(cl)
    }else{
      gamma <- unlist(mclapply(X = ranvars, FUN = gammaFun, mc.cores = ncores))
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
  
  # Association analysis -------------------------------------------------------
  if(verbose == TRUE){
    cat("Performing phenotype-genotype association analysis... ")
  }
  if(ngamma > 0){
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      clusterExport(cl = cl, varlist = list("object","vidx","k","gamma"),
                    envir=environment())
      assoc <- unlist(parLapply(cl = cl, fun = assocgammaFun, X = 1:length(vidx)))
      stopCluster(cl)
    }else{
      assoc <- unlist(mclapply(X = 1:length(vidx), FUN = assocgammaFun, mc.cores = ncores))
    }
  }else{
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      clusterExport(cl = cl, varlist = list("object","vidx","k","Vi"),
                    envir=environment())
      assoc <- unlist(parLapply(cl = cl, fun = assocFun, X = 1:length(vidx)))
      stopCluster(cl)
    }else{
      assoc <- unlist(mclapply(X = 1:length(vidx), FUN = assocFun, mc.cores = ncores))
    }
  }
  if(verbose == TRUE){
    cat("Done.\n")
  }
  assoc <- matrix(data = assoc, ncol = 4, byrow = T)
  if(inherits(object, "GHap.haplo")){
    results <- matrix(data = NA, nrow = length(vidx), ncol = 14)
    results <- as.data.frame(results)
    colnames(results) <- c("CHR","BLOCK","BP1","BP2","ALLELE","FREQ",
                           "N","BETA","SE","CHISQ.EXP","CHISQ.OBS",
                           "CHISQ.GC","LOGP","LOGP.GC")
    results$CHR <- object$chr[vidx]
    results$BLOCK <- object$block[vidx]
    results$BP1 <- object$bp1[vidx]
    results$BP2 <- object$bp2[vidx]
    results$ALLELE <- object$allele[vidx]
  }else{
    results <- matrix(data = NA, nrow = length(vidx), ncol = 13)
    results <- as.data.frame(results)
    colnames(results) <- c("CHR","MARKER","BP","ALLELE","FREQ",
                           "N","BETA","SE","CHISQ.EXP","CHISQ.OBS",
                           "CHISQ.GC","LOGP","LOGP.GC")
    results$CHR <- object$chr[vidx]
    results$MARKER <- object$marker[vidx]
    results$BP <- object$bp[vidx]
    results$ALLELE <- object$A1[vidx]
  }
  results$N <- assoc[,1]
  results$FREQ <- assoc[,2]
  results$BETA <- assoc[,3]
  results$SE <- assoc[,4]
  results$CHISQ.OBS <- (results$BETA/results$SE)^2
  results$LOGP <- -1*pchisq(q = results$CHISQ.OBS, df = 1,
                            lower.tail = FALSE, log.p = TRUE)/log(10)
  if(ngamma > 0  & recalibrate > 0){
    ntop <- ceiling(recalibrate*nrow(results))
    top <- order(results$LOGP, decreasing = TRUE)[1:ntop]
    if(verbose == TRUE){
      cat("Recalibrating statistics for the top ", recalibrate*100,
          "% (", ntop, ") variants... ", sep = "")
    }
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      clusterExport(cl = cl, varlist = list("object","vidx","k","Vi"),
                    envir=environment())
      assoc <- unlist(parLapply(cl = cl, fun = assocFun, X = top))
      stopCluster(cl)
    }else{
      assoc <- unlist(mclapply(X = top, FUN = assocFun, mc.cores = ncores))
    }
    assoc <- matrix(data = assoc, ncol = 4, byrow = T)
    results$BETA[top] <- assoc[,3]
    results$SE[top] <- assoc[,4]
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
