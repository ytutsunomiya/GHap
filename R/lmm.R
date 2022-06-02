#Function: ghap.lmm
#License: GPLv3 or later
#Modification date: 2 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: mixed model fitting

ghap.lmm <- function(
  formula,
  data,
  covmat = NULL,
  weights = NULL,
  vcp.initial = NULL,
  vcp.estimate = TRUE,
  vcp.conv = 1e-12,
  errors = TRUE,
  invcov = FALSE,
  em.reml = 10,
  tol = 1e-12,
  extras = NULL,
  verbose = TRUE
){
  
  # Get formula information -------------------------------------------------
  allterms <- unlist(strsplit(x = as.character(formula)[-1], split = "\\+"))
  response <- gsub(pattern = "\\s+", replacement = "", x = allterms[1])
  ranterms <- grep(pattern = "\\(", x = allterms)
  fixterms <- gsub(pattern = "\\s+", replacement = "", allterms[-c(1,ranterms)])
  ranterms <- gsub(pattern = "\\s+|\\(|\\)", replacement = "", x = allterms[ranterms])
  if (verbose == TRUE) {
    cat("\nAssembling design matrices... ")
  }
  
  # Sanity check for random effects -----------------------------------------
  for(i in 1:length(ranterms)){
    ranname <- unlist(strsplit(x = ranterms[i], split = "\\|"))
    if(ranname[1] != "1"){
      stop("This function does not support random regression yet.\n")
    }else{
      ranterms[i] <- gsub(pattern = "^.+\\|", replacement = "", x = ranterms[i])
    }
  }
  
  # Assemble fixed effects -------------------------------------------------
  form <- paste0(response, "~", paste(fixterms, collapse = "+"))
  X <- sparse.model.matrix(as.formula(form), data=data)
  fixnames <- colnames(X)
  y <- data[,response]
  n <- length(y)
  s <- ncol(X)
  rankX <- as.numeric(rankMatrix(X, method = "qr"))
  
  # Sanity check for covariances --------------------------------------------
  if(is.null(covmat) == FALSE){
    for(i in names(covmat)){
      if(i %in% ranterms == FALSE){
        stop("Random effect ", i, " was not specified in the formula.")
      }else if(is.null(covmat[[i]]) == TRUE){
        lvl <- unique(data[,i])
        data[,i] <- factor(data[,i], levels = sort(lvl))
        covmat[[i]] <- Diagonal(nlevels(data[,i]))
      }else{
        if(identical(colnames(covmat[[i]]),rownames(covmat[[i]])) != TRUE){
          stop("Factors do not match between columns and rows in the ", i," covariance matrix!")
        }
        if(any(data[,i] %in% colnames(covmat[[i]]) == F)){
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
        data[,i] <- factor(data[,i], levels = colnames(covmat[[i]]))
      }
    }
  }else{
    covmat <- vector(mode = "list", length = length(ranterms))
    names(covmat) <- ranterms
    for(i in names(covmat)){
      lvl <- unique(data[,i])
      data[,i] <- factor(data[,i], levels = sort(lvl))
      covmat[[i]] <- Diagonal(nlevels(data[,i]))
    }
  }
  
  # Assemble random effects -------------------------------------------------
  q <- NULL
  Z <- vector(mode = "list", length = length(ranterms))
  rannames <- NULL
  for(i in 1:length(ranterms)){
    form <- paste0("~ 0 + ", ranterms[i])
    Z[[i]] <- sparse.model.matrix(as.formula(form), data=data, drop.unused.levels = FALSE)
    rannames <- c(rannames, colnames(Z[[i]]))
    colnames(Z[[i]]) <- gsub(pattern = ranterms[i], replacement = "", x = colnames(Z[[i]]))
    q <- c(q,ncol(Z[[i]]))
  }
  names(Z) <- ranterms
  names(q) <- ranterms
  vcp.n <- length(ranterms)+1
  
  # Sanity check for variance components ------------------------------------
  if(is.null(vcp.initial) == FALSE){
    if(length(vcp.initial) != vcp.n){
      emsg <- paste0("Number of initial values for variance components (",
                     length(vcp.initial),
                     ") does not match number of random effects (",
                     vcp.n,")\n")
      stop(emsg)
    }
    if(identical(sort(names(vcp.initial)), sort(c(ranterms,"Residual"))) == FALSE){
      emsg <- paste0("Labels of initial values for variance components (",
                     paste(names(vcp.initial), collapse=","),
                     ") do not match labels of random effects (",
                     paste(c(names(ranterms),"Residual"), collapse=","), ")\n")
      stop(emsg)
    }
    vcp.initial <- vcp.initial[c(ranterms,"Residual")]
  }
  
  # Check weights -----------------------------------------------------------
  if(is.null(weights) == FALSE){
    w <- sqrt(weights/mean(weights))
    y <- w*y
    X <- Diagonal(x = w)%*%X
    for(i in 1:length(Z)){
      Z[[i]] <- Diagonal(x = w)%*%Z[[i]]
    }
  }
  
  # Pre-assemble mixed model components -------------------------------------
  if(is.null(vcp.initial) == TRUE){
    vcp.new <- rep(var(y)/vcp.n, times = vcp.n)
  }else{
    vcp.new <- unlist(vcp.initial)
  }
  names(vcp.new) <- c(ranterms, "Residual")
  RHS <- crossprod(X, y)
  LHS <- crossprod(X)
  W <- X
  ZtX <- vector(mode = "list", length = vcp.n-1)
  ZtZ <- vector(mode = "list", length = vcp.n-1)
  names(ZtX) <- ranterms
  names(ZtZ) <- ranterms
  for(i in ranterms){
    RHS <- rbind(RHS, crossprod(Z[[i]], y))
    ZtX[[i]] <- crossprod(Z[[i]], X)
    LHS <- cbind(LHS, t(ZtX[[i]]))
    ZtZ[[i]] <- vector(mode = "list", length = length(ranterms))
    names(ZtZ[[i]]) <- ranterms
    for(j in ranterms){
      ZtZ[[i]][[j]] <- crossprod(Z[[i]],Z[[j]])
    }
    W <- cbind(W, Z[[i]])
  }
  for(i in ranterms){
    tmp <- ZtX[[i]]
    for(j in ranterms){
      tmp <- cbind(tmp,ZtZ[[i]][[j]])
    }
    LHS <- rbind(LHS, tmp)
  }
  
  # Initialize auxiliary functions ------------------------------------------
  sparsesolve <- function(X){
    Xchol <- cholPermute(X)
    X <- Takahashi_Davis(X, cholQp = Xchol$Qpermchol, P = Xchol$P)
    return(X)
  }
  pcgsolve <- function(A,b,criterion=1e-12){
    m = 0
    x <- rep(0, times=length(b))
    r <- b
    Minv <- (1/diag(A))
    z <- Minv*r
    p <- z
    rss.new <- sum(r*z)
    ss <- sum(b^2)
    conv <- 100
    while(conv > criterion){
      m <- m + 1
      g <- A%*%p
      alpha <- as.numeric(rss.new/sum(p*g))
      x <- x + alpha*p
      if(m %% 50 == 0){
        r <- b - A%*%x
      }else{
        r <- r - alpha*g
      }
      z <- Minv*r
      rss.old <- rss.new
      rss.new <- sum(r*z)
      beta <- rss.new/rss.old
      p <- z + beta*p
      conv <- rss.new/ss
    }
    return(as.numeric(x))
  }
  
  # Log message -------------------------------------------------------------
  if (verbose == TRUE) {
    cat("Done.\n")
    cat(n, " records will be fitted to\n ", s, " fixed effects\n ",
        sum(q), " random effects\n", sep="")
  }
  
  # Check if covariance matrices should be inverted -------------------------
  if(invcov == FALSE){
    if(verbose == TRUE){
      cat("Inverting covariance matrices...")
    }
    oricovmat <- covmat
    for(i in names(covmat)){
      icov <- try(solve(covmat[[i]]), silent = TRUE)
      if(inherits(icov, "try-error")){
        icov <- try(solve(covmat[[i]] + Diagonal(nrow(covmat[[i]]))*tol), silent = TRUE)
        if(inherits(icov, "try-error")){
          emsg <- paste0("\nUnable to invert covariance matrix for effect ",
                         ranterms[[i]], " even after adding a tolerance of ", tol)
          stop(emsg)
        }
      }
      covmat[[i]] <- icov
    }
    if(verbose == TRUE){
      cat(" Done.\n")
    }
  }
  
  # Variance components estimation ------------------------------------------
  if(vcp.estimate == TRUE){
    convall <- 100
    k <- 0
    k.reml <- 0
    if(verbose == TRUE){
      cat("\nVariance components:\n\n",
          "--------------------------------------------\n", sep="")
    }
    while(convall >= vcp.conv){
      
      # Increment iteration
      k <- k + 1
      k.reml <- k.reml + 1
      
      # Include variance components in mixed model equations
      for(i in 1:length(ranterms)){
        if(i == 1){
          idx <- (s+1):(s+q[i])
        }else{
          idx <- which(1:length(q) < i)
          idx <- (s+sum(q[idx])+1):(s+sum(q[1:i]))
        }
        LHS[idx,idx] <- ZtZ[[i]][[i]] + covmat[[ranterms[i]]]*(vcp.new["Residual"]/vcp.new[i])
      }
      LHS <- LHS/vcp.new["Residual"]
      RHS <- RHS/vcp.new["Residual"]
      
      # Solve mixed model equations
      coef <- pcgsolve(LHS,RHS)
      
      # Update fitted values
      fit <- as.numeric(W%*%coef)
      e <- y - fit
      
      # Get sparse solution of left hand side
      LHSi <- try(sparsesolve(LHS), silent = TRUE)
      if(inherits(LHSi, "try-error")){
        LHSi <- try(sparsesolve(LHS + Diagonal(nrow(LHS))*tol), silent = TRUE)
        if(inherits(LHSi, "try-error")){
          emsg <- paste0("\nUnable to solve LHS even after adding a tolerance of ",
                         tol)
          stop(emsg)
        }
      }
      
      # Get REML updates of variance components
      seconddev <- NULL
      vcp.old <- vcp.new
      if(k.reml <= em.reml){
        itermethod <- "EM-REML"
        vcp.new["Residual"] <- sum(y*e)/(n-rankX)
        for(i in 1:length(ranterms)){
          if(i == 1){
            idx <- (s+1):(s+q[i])      
          }else{
            idx <- which(1:length(q) < i)
            idx <- (s+sum(q[idx])+1):(s+sum(q[1:i]))
          }
          su <- as.numeric(t(coef[idx])%*%covmat[[ranterms[i]]]%*%coef[idx])
          tr <- sum(diag(covmat[[ranterms[i]]]%*%LHSi[idx,idx]))
          vcp.new[i] <- as.numeric(su + tr)/q[i]
        }
        convall <- 100
      }else{
        itermethod <- "AI-REML"
        Rinv <- Diagonal(length(y))/vcp.new["Residual"]
        Py <- Rinv%*%e
        firstdev <- NULL
        Q <- NULL
        D <- NULL
        for(i in 1:length(ranterms)){
          if(i == 1){
            idx <- (s+1):(s+q[i])
          }else{
            idx <- which(1:length(q) < i)
            idx <- (s+sum(q[idx])+1):(s+sum(q[1:i]))
          }
          su <- as.numeric(t(coef[idx])%*%covmat[[ranterms[i]]]%*%coef[idx])
          tr <- sum(diag(covmat[[ranterms[i]]]%*%LHSi[idx,idx]))
          firstdev <- c(firstdev, (q[i] - (su + tr)/vcp.new[i])/vcp.new[i])
          Q <- cbind(Q, Z[[i]]%*%as.numeric((coef[idx]/vcp.new[i])))
          RHSf <- crossprod(X, Q[,i])
          for(j in 1:length(ranterms)){
            RHSf <- rbind(RHSf, crossprod(Z[[ranterms[j]]], Q[,i]))
          }
          RHSf <- RHSf/vcp.new["Residual"]
          D <- cbind(D, pcgsolve(LHS, RHSf))
        }
        Q <- cbind(Q, Py)
        RHSf <- crossprod(X,Py)
        for(i in ranterms){
          RHSf <- rbind(RHSf, crossprod(Z[[i]], Py))
        }
        RHSf <- RHSf/vcp.new["Residual"]
        D <- cbind(D, pcgsolve(LHS, RHSf))
        trP <- sum(diag(Rinv)) - sum(diag(LHSi%*%crossprod(Rinv%*%W)))
        firstdev <- c(firstdev, as.numeric(trP - crossprod(Py)))
        RinvQ <- Rinv%*%Q
        WtRinvQ <- crossprod(W, RinvQ)
        seconddev <- solve(crossprod(Q, RinvQ) - crossprod(D,WtRinvQ))
        if(isSymmetric(seconddev) == FALSE){
          seconddevt <- t(seconddev)
          seconddev[upper.tri(seconddev)] <- seconddevt[upper.tri(seconddevt)]
        }
        delta <- as.numeric(seconddev%*%firstdev)
        vcp.new <- vcp.old - delta
        if(any(vcp.new < 0)){
          k.reml <- 0
          vcp.new <- vcp.old
        }
        convall <- sum((vcp.new - vcp.old)^2)/(sum(vcp.new)^2)
      }
      
      # Remove residual variance from mixed model equations
      LHS <- LHS*vcp.old["Residual"]
      RHS <- RHS*vcp.old["Residual"]
      
      # Log message for variance components update
      if(verbose == TRUE){
        if(itermethod == "EM-REML"){
          cat(
            "Variance components iteration: ", k, "\n",
            "Method: ", itermethod, "\n  ",
            paste(paste(c(ranterms,"Residual"), "=", sprintf("%.12f", vcp.new)), collapse = "\n  "),
            "\n  Convergence = ---",
            "\n--------------------------------------------\n", sep="")
        }else{
          cat(
            "Variance components iteration: ", k, "\n",
            "Method: ", itermethod, "\n  ",
            paste(paste(c(ranterms,"Residual"), "=", sprintf("%.12f", vcp.new)), collapse = "\n  "),
            "\n  Convergence = ", convall,
            "\n--------------------------------------------\n", sep="")
        }
      }
    }
    if(verbose == TRUE & vcp.estimate == TRUE){
      cat("\n\nVariance components converged.\n")
    }
  }else{
    if(verbose == TRUE){
      cat("Variance components assumed known\n  ",
          paste(paste(c(ranterms,"Residual"), "=", sprintf("%.12f", vcp.new)), collapse = "\n  "), "\n", sep="")
    }
    seconddev <- NULL
  }
  
  # Get solutions -----------------------------------------------------------
  for(i in 1:length(ranterms)){
    if(i == 1){
      idx <- (s+1):(s+q[i])
    }else{
      idx <- which(1:length(q) < i)
      idx <- (s+sum(q[idx])+1):(s+sum(q[1:i]))
    }
    LHS[idx,idx] <- ZtZ[[i]][[i]] + covmat[[ranterms[i]]]*(vcp.new["Residual"]/vcp.new[i])
  }
  LHS <- LHS/vcp.new["Residual"]
  RHS <- RHS/vcp.new["Residual"]
  coef <- pcgsolve(LHS,RHS)
  
  # Check extras ---------------------------------------------------------------
  
  # Full LHSi
  if("LHSi" %in% extras){
    if(verbose == TRUE){
      cat("Computing full inverse of coefficient matrix... ")
    }
    LHSi <- try(solve(LHS), silent = TRUE)
    if(inherits(LHSi, "try-error")) {
      LHSi <- solve(LHS + Diagonal(nrow(LHS))*tol)
    }
    colnames(LHSi) <- c(fixnames,rannames)
    rownames(LHSi) <- colnames(LHSi)
    if(verbose == TRUE){
      cat("Done.\n")     
    }
  }
  
  # Phenotypic variance-covariance matrix
  if("V" %in% extras){
    if(verbose == TRUE){
      cat("Computing phenotypic (co)variance matrix... ")
    }
    if(is.null(weights) == TRUE){
      w <- rep(1, times = n)
    }
    V <- Diagonal(n)
    diag(V) <- w*vcp.new["Residual"]
    for(i in ranterms){
      if(invcov == TRUE){
        V <- V + tcrossprod(Z[[i]]%*%solve(covmat[[i]])*vcp.new[i],Z[[i]])
      }else{
        V <- V + tcrossprod(Z[[i]]%*%oricovmat[[i]]*vcp.new[i],Z[[i]])
      }
    }
    if(verbose == TRUE){
      cat("Done.\n")
    }
  }
  
  # Calculate standard errors --------------------------------------------------
  if(errors == TRUE){
    if(verbose == TRUE){
      cat("Computing standard errors... ")      
    }
    if("LHSi" %in% extras == FALSE){
      LHSi <- try(sparsesolve(LHS), silent = TRUE)
      if(inherits(LHSi, "try-error")) {
        LHSi <- sparsesolve(LHS + Diagonal(nrow(LHS))*tol)
      }
    }
    stde <- sqrt(diag(LHSi))
    if(is.null(seconddev) == FALSE){
      vcp.stde <- sqrt(diag(seconddev))
    }
    if(verbose == TRUE){
      cat("Done.")     
    }
  }
  
  # Return results ----------------------------------------------------------
  results <- NULL
  tmp <- as.data.frame(matrix(data = NA, nrow = s, ncol = 4))
  colnames(tmp) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(tmp) <- colnames(X)
  tmp[,1] <- coef[1:s]
  if(errors == TRUE){
    tmp[,2] <- stde[1:s]
    tmp[,3] <- tmp[,1]/tmp[,2]
    a <- pt(q = -abs(tmp[,3]), df = n - s, lower.tail = T, log.p = TRUE)
    b <- pt(q = abs(tmp[,3]), df = n - s, lower.tail = F, log.p = TRUE)
    tmp[,4] <- 10^((a + log(1 + exp(b - a)))/log(10))
  }
  results$fixed <- tmp
  results$random <- vector(mode = "list", length = length(ranterms))
  names(results$random) <- ranterms
  for(i in 1:length(ranterms)){
    if(i == 1){
      idx <- (s+1):(s+q[i])
    }else{
      idx <- which(1:length(q) < i)
      idx <- (s+sum(q[idx])+1):(s+sum(q[1:i]))
    }
    tmp <- as.data.frame(matrix(data = NA, nrow = length(idx), ncol = 3))
    colnames(tmp) <- c("Estimate", "Std. Error", "Accuracy")
    rownames(tmp) <- colnames(Z[[ranterms[i]]])
    tmp[,1] <- coef[idx]
    if(errors == TRUE){
      tmp[,2] <- stde[idx]
      tmp[,3] <- 1 - (stde[idx]^2)/vcp.new[i]
      tmp[which(tmp[,3] <= 0),3] <- 1e-6
      tmp[,3] <- sqrt(tmp[,3])      
    }
    results$random[[ranterms[i]]] <- tmp
  }
  if(is.null(weights) == TRUE){
    w <- rep(1, times = length(y))
  }
  tmp <- as.data.frame(matrix(data = NA, nrow = length(y), ncol = 3))
  colnames(tmp) <- c("All", "Fixed", "Random")
  tmp[,1] <- as.numeric(W%*%coef)/w
  if(s == 1){
    tmp[,2] <- as.numeric(W[,1]*coef[1])/w
  }else{
    tmp[,2] <- as.numeric(W[,1:s]%*%coef[1:s])/w
  }
  
  tmp[,3] <- as.numeric(W[,-c(1:s)]%*%coef[-c(1:s)])/w
  results$fitted <- tmp
  tmp <- as.data.frame(matrix(data = NA, nrow = length(y), ncol = 3))
  colnames(tmp) <- c("All", "Fixed", "Random")
  tmp[,1] <- (y/w) - results$fitted[,1]
  tmp[,2] <- (y/w) - results$fitted[,2]
  tmp[,3] <- (y/w) - results$fitted[,3]
  results$residuals <- tmp
  tmp <- as.data.frame(matrix(data = NA, nrow = vcp.n, ncol = 2))
  colnames(tmp) <- c("Estimate", "Std. Error")
  rownames(tmp) <- names(vcp.new)
  tmp[,1] <- vcp.new
  if(is.null(seconddev) == FALSE & errors == TRUE){
    tmp[,2] <- vcp.stde
    results$vcp <- tmp
    results$AI <- seconddev
    colnames(results$AI) <- names(vcp.new)
    rownames(results$AI) <- names(vcp.new)
  }else{
    results$vcp <- tmp
  }
  results$extras <- NULL
  if("LHSi" %in% extras){
    results$extras$LHSi <- LHSi
  }
  if("V" %in% extras){
    results$extras$V <- V
  }
  return(results)
  
}
