#Function: ghap.lmm
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Mixed linear model fitting

ghap.lmm<-function(
  fixed,
  random,
  covmat=NULL,       #List with covariance matrices for random effects
  data,              #Data frame containing model data
  weights = NULL,    #Weights of observations
  family = "gaussian",
  REML = TRUE,
  verbose=TRUE
){
  
  #Check family
  if(is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame(2))
  }
  if(is.function(family)){
    family <- family()
  }
  
  #Log message
  if(verbose==TRUE){
    cat("\nAssuming", family$family,"family with", family$link,"link function.\n")
  }
  
  #Make formula
  rand.labels <- attr(terms(random),"term.labels")
  form <- paste(format(fixed),paste("(1 | ",rand.labels,")",sep="",collapse=" + "),sep=" + ")
  form <- as.formula(form)
  
  #Log message
  if(verbose==TRUE){
    cat("\nAssembling design matrices... ")
  }
  
  #Sanity check for covariances
  cov.labels <- NULL
  if(is.null(covmat) == FALSE){
    cov.labels <- names(covmat)
    for(i in cov.labels){
      if(i %in% rand.labels == FALSE){
        stop("Random effect ", i, " was not specified in the formula.")
      }else{
        if(identical(colnames(covmat[[i]]),rownames(covmat[[i]])) != TRUE){
          stop("Factors do not match between columns and rows in the ", i," covariance matrix!")
        }
        if(any(data[,i] %in% colnames(covmat[[i]]))==F){
          stop("There are factors in ", i, " without covariance")
        }
        if(length(which(is.na(colnames(covmat[[i]])) == T)) != 0){
          stop("NAs as factors in the ", i," convariance matrix")
        }
        cvr.dup <- table(colnames(covmat[[i]]))
        cvr.dup <- length(which(cvr.dup > 1))
        if(cvr.dup != 0){
          stop("Duplicated factors declared in the ", i," covariance matrix!")
        }
        rm(cvr.dup)
        data[,i] <- factor(data[,i], levels = colnames(covmat[[i]]), ordered = TRUE)
        covmat[[i]] <- chol(nearPD(covmat[[i]])$mat)
      }
    }
  }
  
  #Assemble fixed effects for lme4
  mf <- model.frame(subbars(form),data=data)
  if(is.null(weights) == TRUE){
    w <- rep(1,times=nrow(mf))
  }else{
    w <- sqrt(weights/mean(weights))
  }
  mf[,1] <- w*mf[,1]
  X <- w*sparse.model.matrix(nobars(form),mf)
  
  #Assemble random effects for lme4
  rt <- mkReTrms(bars = findbars(form),fr = mf, drop.unused.levels = FALSE)
  for(i in 1:length(rt$Ztlist)){
    rt$Ztlist[[i]] <- rt$Ztlist[[i]]%*%Diagonal(x = w)
  }
  if(is.null(covmat) == FALSE){
    for(i in cov.labels){
      linklab <- paste("1 | ",i,sep="")
      rt$Ztlist[[linklab]] <- covmat[[i]]%*%rt$Ztlist[[linklab]]
      rt$Ztlist[[linklab]] <- as(rt$Ztlist[[linklab]], 'dgCMatrix')
    }
  }
  Zt <- rt$Ztlist[[1]]
  map <- rep(rand.labels[1],times=nrow(rt$Ztlist[[1]]))
  if(length(rt$Ztlist) > 1){
    for(i in 2:length(rt$Ztlist)){
      Zt <- rbind(Zt, rt$Ztlist[[i]])
      map <- rep(rand.labels[i],times=nrow(rt$Ztlist[[i]]))
    }
  }
  rt$Zt <- as(Zt,'dgCMatrix')
  rm(Zt); a<-gc()
  lmod <- list(fr = mf, X = X, reTrms = rt, start = rt$theta, REML = REML)
  
  #Number of parameters
  n <- nrow(mf)
  s <- ncol(X)
  q <- nrow(rt$Zt)
  
  #Log message
  if(verbose==TRUE){
    cat("Done.\n")
    cat(n,"records will be fitted to an intercept,",s-1,"fixed effects and",q,"random effects.\n")
    if(REML == TRUE & family$family == "gaussian" & family$link == "identity"){
      cat("\nMaximizing restricted likelihood... ")
    }else{
      cat("\nMaximizing likelihood... ")
    }
  }
  
  #Maximize
  if(family$family == "gaussian" & family$link == "identity"){
    devfun <- do.call(mkLmerDevfun,lmod)
    opt <- optimizeLmer(devfun, optimizer = "Nelder_Mead")
    lmout <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
  }else{
    lmod$family <- family
    devfun <- do.call(mkGlmerDevfun,lmod)
    opt <- optimizeGlmer(devfun, optimizer = "Nelder_Mead")
    lmout <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
  }
  
  #Assemble results
  results <- NULL
  results$fixed <- fixef(lmout)
  results$random <- ranef(lmout)
  results$residuals <- as.numeric(residuals(lmout))/w
  vcp <- as.data.frame(VarCorr(lmout))
  results$vcp <- vcp$vcov
  names(results$vcp) <- vcp$grp
  for(i in names(results$random)){
    if(i %in% cov.labels){
      rnames <- rownames(results$random[[i]])
      ref <- as.numeric(covmat[[i]]%*%unlist(results$random[[i]]))
      results$random[[i]] <- NULL
      results$random[[i]] <- ref
      names(results$random[[i]]) <- rnames
    }else{
      rnames <- rownames(results$random[[i]])
      ref <- as.numeric(unlist(results$random[[i]]))
      results$random[[i]] <- NULL
      results$random[[i]] <- ref
      names(results$random[[i]]) <- rnames
    }
  }
  results$lme4 <- lmout
  
  #Log message
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  #Results
  return(results)
  
  
}