#Function: ghap.anctrain
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Create prototype alleles for ancestry predictions

ghap.anctrain <- function(
  object,
  train = NULL,
  method = "unsupervised",
  K = 2,
  iter.max = 10,
  nstart = 10,
  nmarkers = 5000,
  tune = FALSE,
  only.active.samples = TRUE,
  only.active.markers = TRUE,
  batchsize=NULL,
  ncores = 1,
  verbose = TRUE
){
  
  # Check if phase is a GHap.phase object-------------------------------------------------------------
  if(inherits(object, "GHap.phase") == FALSE){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Check if method is valid--------------------------------------------------------------------------
  if(method %in% c("supervised","unsupervised") == FALSE){
    stop("Method should be either 'supervised' or 'unsupervised.")
  }
  
  # Check if inactive markers and samples should be reactived-----------------------------------------
  if(only.active.markers == FALSE){
    object$marker.in <- rep(TRUE,times=object$nmarkers)
    object$nmarkers.in <- length(which(object$marker.in))
  }
  if(only.active.samples == FALSE){
    object$id.in <- rep(TRUE,times=2*object$nsamples)
    object$nsamples.in <- length(which(object$id.in))/2
  }
  
  # Map training samples -----------------------------------------------------------------------------
  if(is.null(train) == TRUE){
    train.idx <- which(object$id.in == TRUE)
  }else{
    train.idx <- which(object$id %in% train & object$id.in == TRUE)
  }
  
  # Map population for supervised analysis -----------------------------------------------------------
  if(method == "supervised"){
    train.pop <- object$pop[train.idx]
    y <- object$pop[train.idx]
    y <- as.factor(y)
  }
  
  # Map parameters to use ----------------------------------------------------------------------------
  if(method == "unsupervised"){
    param <- list(K = K, iter.max = iter.max, nstart = nstart, nmarkers = nmarkers, tune = tune)
  }else{
    param <- table(y)
  }
  ncores <- min(c(detectCores(), ncores))
  
  # Log message of parameters-------------------------------------------------------------------------
  if(verbose == TRUE){
    if(method == "unsupervised"){
      printparams <- paste(names(param), "=", param, collapse=", ")
      cat("\nUsing method 'unsupervised' with parameters:\n[", printparams, "]\n", sep="")
    }else{
      printparams <- paste(names(param), " [n = ", param, "]", sep="", collapse="\n")
      cat("\nUsing method 'supervised' with reference haplotypes:\n", printparams, "\n", sep="")
    }
  }
  
  
  # Seed and tuning for kmeans------------------------------------------------------------------------
  if(method == "unsupervised"){
    mkr <- sample(x = which(object$marker.in), size = param$nmarkers, replace = FALSE)
    Mkm <- ghap.slice(object = object, ids = train.idx, variants = mkr, transposed = TRUE,
                      index = TRUE, ncores = ncores, verbose = FALSE)
    if(tune == TRUE){
      tune.FUN <- function(i){
        clk <- kmeans(x = Mkm, centers = i,
                      iter.max = param$iter.max, nstart = param$nstart)
        clout <- (clk$betweenss/clk$tot.withinss)*(sum(clk$size) - i)/(i-1)
        clout <- c(clout,clk$tot.withinss)
        return(clout)
      }
      if(verbose == TRUE){
        cat("\nQuantifying within-cluster dispersion from K = 1 to K = ", param$K, "... ", sep="")
      }
      if(ncores == 1){
        clout <- unlist(lapply(X = 1:param$K, FUN = tune.FUN))
        clout <- as.data.frame(matrix(data = clout, ncol = 2, byrow = T))
      }else{
        if(Sys.info()["sysname"] == "Windows"){
          cl <- makeCluster(ncores)
          clusterEvalQ(cl, library(Matrix))
          varlist <- list("Mkm","param")
          clusterExport(cl = cl, varlist = varlist, envir=environment())
          clout <- unlist(parLapply(cl = cl, fun = tune.FUN, X = 1:param$K))
          stopCluster(cl)
          clout <- as.data.frame(matrix(data = clout, ncol = 2, byrow = T))
        }else{
          clout <- unlist(mclapply(X = 1:param$K, FUN = tune.FUN, mc.cores = ncores))
          clout <- as.data.frame(matrix(data = clout, ncol = 2, byrow = T))
        }
      }
      colnames(clout) <- c("chi","sst")
      clout$chi[1] <- 0
      if(verbose == TRUE){
        cat("Done.\n")
      }
      sst <- clout$sst
      change <- 100*diff(sst)/sst[-param$K]
      names(change) <- paste("K", 2:param$K, " - K", 1:(param$K-1), sep="")
      # ymin <- min(sst)
      # ymax <- max(sst)
      # par(mfrow=c(1,2))
      # plot(x = 1:K, y = sst, ylim = c(ymin,ymax), type = "b", yaxt = "n", xaxt = "n",
      #      xlab="K value", ylab = "Total within-cluster sum of squares", col = "darkgrey", lwd=2,
      #      main = "Elbow method")
      # ssst <- seq(from = ymin, to = ymax, length.out = 5)
      # axis(side = 2, at = ssst, labels = sprintf("%.2g", ssst), las=3)
      # axis(side = 1, at = 1:param$K, labels = 1:K, las=1)
      # text(x = (2:param$K)-0.5, y = (sst[-length(sst)] + sst[-1])/2, 
      #      labels = paste(sprintf("%.1f", change),"%"), pos = 3)
      # plot(x = 1:K, y = clout$chi,type = "b", xlab="K value", ylab = "Calinski-Harabasz (CH) Index",
      #      col = "darkgrey", lwd=2, xaxt = "n",
      #      main = "Variance Ratio Criterion", las=1)
      # axis(side = 1, at = 1:param$K, labels = 1:K, las=1)
    }else{
      if(verbose == TRUE){
        cat("\nGrouping haplotypes into K = ", K," pseudo-lineages using K-means clustering... ", sep="")
      }
      clk <- kmeans(x = Mkm, centers = param$K,
                    iter.max = param$iter.max, nstart = param$nstart)
      y <- paste0("K",clk$cluster)
      if(verbose == TRUE){
        cat("Done.\n")
      }
    }
  }
  
  # Generate batch index -----------------------------------------------------------------------------
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(object$nmarkers.in/10)
  }
  if(batchsize > object$nmarkers.in){
    batchsize <- object$nmarkers.in
  }
  id1 <- seq(1,object$nmarkers.in,by=batchsize)
  id2 <- (id1+batchsize)-1
  id1 <- id1[id2<=object$nmarkers.in]
  id2 <- id2[id2<=object$nmarkers.in]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,object$nmarkers.in)
  if(id1[length(id1)] > object$nmarkers.in){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  # Prototype allele function ------------------------------------------------------------------------
  proto.fun <- function(k){
    x <- X[k,]
    dfp <- data.frame(geno = x, pop = y)
    res <- aggregate(formula = geno ~ pop, data = dfp, FUN = mean)
    return(res$geno)
  }
  
  # Prototype alleles calculation --------------------------------------------------------------------
  if(method == "unsupervised" & tune == TRUE){
    results <- NULL
    results$ssq <- clout$sst
    results$chindex <- clout$chi
    results$pchange <- change
  }else{
    snps.in <- which(object$marker.in)
    poplabs <- sort(unique(as.character(y)))
    results <- matrix(data = NA, nrow = length(snps.in), ncol = length(poplabs)+1)
    results <- as.data.frame(results)
    colnames(results) <- c("MARKER",poplabs)
    results$MARKER <- object$marker[snps.in]
    if(verbose == TRUE){
      cat("\nBuilding prototype alleles... ")
    }
    for(i in 1:length(id1)){
      X <- ghap.slice(object = object,
                      ids = train.idx,
                      variants = snps.in[id1[i]:id2[i]],
                      index = TRUE,
                      ncores = ncores)
      #Compute blocks
      if(ncores == 1){
        p <- unlist(lapply(FUN = proto.fun, X = 1:nrow(X)))
      }else{
        if(Sys.info()["sysname"] == "Windows"){
          cl <- makeCluster(ncores)
          clusterEvalQ(cl, library(Matrix))
          varlist <- list("X","y")
          clusterExport(cl = cl, varlist = varlist, envir=environment())
          p <- unlist(parLapply(cl = cl, fun = proto.fun, X = 1:nrow(X)))
          stopCluster(cl)
        }else{
          p <- unlist(mclapply(FUN = proto.fun, X = 1:nrow(X), mc.cores = ncores))
        }
      }
      p <- data.frame(matrix(p, ncol=length(poplabs), byrow=TRUE), stringsAsFactors = F)
      results[id1[i]:id2[i],-1] <- p
    }
    if(verbose == TRUE){
      cat("Done.\n")
    }
  }
  
  # Return results------------------------------------------------------------------------------------
  return(results)
  
}
