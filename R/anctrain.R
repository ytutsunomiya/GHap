#Function: ghap.anctrain
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Create prototype alleles for ancestry predictions

ghap.anctrain <- function(
  phase,
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
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Check if method is valid--------------------------------------------------------------------------
  if(method %in% c("supervised","unsupervised") == FALSE){
    stop("Method should be either 'supervised' or 'unsupervised.")
  }
  
  # Check if inactive markers and samples should be reactived-----------------------------------------
  if(only.active.markers == FALSE){
    phase$marker.in <- rep(TRUE,times=phase$nmarkers)
    phase$nmarkers.in <- length(which(phase$marker.in))
  }
  if(only.active.samples == FALSE){
    phase$id.in <- rep(TRUE,times=2*phase$nsamples)
    phase$nsamples.in <- length(which(phase$id.in))/2
  }
  
  # Map training samples -----------------------------------------------------------------------------
  if(is.null(train) == TRUE){
    train.idx <- which(phase$id.in == TRUE)
  }else{
    train.idx <- which(phase$id %in% train & phase$id.in == TRUE)
  }
  
  # Map population for supervised analysis -----------------------------------------------------------
  if(method == "supervised"){
    train.pop <- phase$pop[train.idx]
    y <- phase$pop[train.idx]
    y <- as.factor(y)
  }

  # Map parameters to use ----------------------------------------------------------------------------
  if(method == "unsupervised"){
    param <- list(K = K, iter.max = iter.max, nstart = nstart, nmarkers = nmarkers, tune = tune)
  }else{
    param <- table(y)
  }
  
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
  
  # Get number of cores-------------------------------------------------------------------------------
  if(Sys.info()["sysname"] == "Windows"){
    if(ncores > 1 & verbose == TRUE){
      cat("\nParallelization not supported yet under Windows (using a single core).")
    }
    ncores <- 1
  }else{
    ncores <- min(c(detectCores(), ncores))
  }
  
  # Initialize lookup table----------------------------------------------------------------------------
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
  
  # Seed and tuning for kmeans------------------------------------------------------------------------
  if(method == "unsupervised"){
    mkr <- sample(x = which(phase$marker.in), size = param$nmarkers, replace = FALSE)
    Mkm <- ghap.pslice(phase = phase, ids = train.idx, markers = mkr,
                       index = TRUE, lookup = lookup, ncores = ncores, verbose = FALSE)
    if(tune == TRUE){
      tune.FUN <- function(i){
        cl <- kmeans(x = t(Mkm), centers = i,
                     iter.max = param$iter.max, nstart = param$nstart)
        clout <- (cl$betweenss/cl$tot.withinss)*(sum(cl$size) - i)/(i-1)
        clout <- c(clout,cl$tot.withinss)
        return(clout)
      }
      if(verbose == TRUE){
        cat("\nQuantifying within-cluster dispersion from K = 1 to K = ", param$K, "... ", sep="")
      }
      if(Sys.info()["sysname"] == "Windows"){
        clout <- unlist(lapply(X = 1:param$K, FUN = tune.FUN))
        clout <- as.data.frame(matrix(data = clout, ncol = 2, byrow = T))
        colnames(clout) <- c("chi","sst")
        clout$chi[1] <- 0
      }else{
        clout <- unlist(mclapply(X = 1:param$K, FUN = tune.FUN, mc.cores = ncores))
        clout <- as.data.frame(matrix(data = clout, ncol = 2, byrow = T))
        colnames(clout) <- c("chi","sst")
        clout$chi[1] <- 0
      }
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
      cl <- kmeans(x = t(Mkm), centers = param$K,
                   iter.max = param$iter.max, nstart = param$nstart)
      y <- paste0("K",cl$cluster)
      if(verbose == TRUE){
        cat("Done.\n")
      }
    }
  }
  
  # Generate batch index -----------------------------------------------------------------------------
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(phase$nmarkers.in/10)
  }
  if(batchsize > phase$nmarkers.in){
    batchsize <- phase$nmarkers.in
  }
  id1 <- seq(1,phase$nmarkers.in,by=batchsize)
  id2 <- (id1+batchsize)-1
  id1 <- id1[id2<=phase$nmarkers.in]
  id2 <- id2[id2<=phase$nmarkers.in]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,phase$nmarkers.in)
  if(id1[length(id1)] > phase$nmarkers.in){
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
    snps.in <- which(phase$marker.in)
    poplabs <- sort(unique(as.character(y)))
    results <- matrix(data = NA, nrow = length(snps.in), ncol = length(poplabs)+1)
    results <- as.data.frame(results)
    colnames(results) <- c("MARKER",poplabs)
    results$MARKER <- phase$marker[snps.in]
    if(verbose == TRUE){
      cat("\nBuilding prototype alleles... ")
    }
    for(i in 1:length(id1)){
      X <- ghap.pslice(phase = phase,
                       ids = train.idx,
                       markers = snps.in[id1[i]:id2[i]],
                       index = TRUE,
                       lookup = lookup,
                       ncores = ncores)
      #Compute blocks
      if(Sys.info()["sysname"] == "Windows"){
        p <- unlist(lapply(FUN = proto.fun, X = 1:nrow(X)))
      }else{
        p <- unlist(mclapply(FUN = proto.fun, X = 1:nrow(X), mc.cores = ncores))
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