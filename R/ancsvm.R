#Function: ghap.ancsvm
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Predict ancestry of haplotypes using machine learning

ghap.ancsvm <- function(
  object,
  blocks = NULL,
  test = NULL,
  train = NULL,
  cost = 1,
  gamma = NULL,
  tune = FALSE,
  only.active.samples = TRUE,
  only.active.markers = TRUE,
  ncores = 1,
  verbose = TRUE
){
  
  # Check if object is a GHap.phase object-------------------------------------------------------------
  if(inherits(object, "GHap.phase") == FALSE){
    stop("Argument object must be a GHap.phase object.")
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
  
  # Map test samples----------------------------------------------------------------------------------
  if(is.null(test) == TRUE){
    test.idx <- which(object$id.in == TRUE)
  }else{
    test.idx <- which(object$id %in% test & object$id.in == TRUE)
  }
  
  # Map training samples -----------------------------------------------------------------------------
  if(is.null(train) == TRUE){
    train.idx <- which(object$id.in == TRUE)
  }else{
    train.idx <- which(object$id %in% train & object$id.in == TRUE)
  }
  y <- object$pop[train.idx]
  y <- as.factor(y)
  
  # Check if blocks exist ------------------------------------------------------
  if(is.null(blocks) == TRUE){
    
    # Calculate marker density
    mrkdist <- diff(object$bp)
    mrkdist <- mrkdist[which(mrkdist > 0)]
    density <- mean(mrkdist)
    
    # Generate blocks for admixture events
    g <- 10
    window <- (100e+6)/(2*g)
    window <- ceiling(window/density)
    step <- ceiling(window/4)
    blocks <- ghap.blockgen(object, windowsize = window,
                            slide = step, unit = "marker")
  }
  
  # Map parameters to use-----------------------------------------------------------------------------
  param <- list(cost = cost, gamma = gamma, tune = tune)
  if(is.null(param$gamma) == TRUE){
    param$gamma <- "1/blocksize"
  }
  ncores <- min(c(detectCores(), ncores))
  
  # Log message of parameters-------------------------------------------------------------------------
  if(verbose == TRUE){
    printparams <- paste(names(param), "=", param, collapse=", ")
    cat("\nUsing svm with parameters:\n[", printparams, "]\n", sep="")
  }
  
  # Initialize block iteration function---------------------------------------------------------------
  blockfun <- function(b){
    
    #Get block info
    block.info <- blocks[b, c("BLOCK","CHR","BP1","BP2")]
    
    #SNPs in the block
    snps <- which(object$chr == block.info$CHR &
                    object$bp >= block.info$BP1 &
                    object$bp <= block.info$BP2 &
                    object$marker.in == TRUE)
    blocksize <- length(snps)
    
    #Build model matrices
    Mtst <- ghap.slice(object = object, ids = test.idx, variants = snps,
                       index = TRUE, transposed = TRUE, verbose = FALSE)
    Mref <- ghap.slice(object = object, ids = train.idx, variants = snps,
                       index = TRUE, transposed = TRUE, verbose = FALSE)
    
    #Model training
    if(param$gamma == "1/blocksize"){
      gamma <- 1/blocksize
    }else{
      gamma <- param$gamma
    }
    model <- svm(y = y, x = Mref, kernel = "radial",
                 gamma = gamma, cost = param$cost)
    pred <- predict(model, Mtst)
    pred <- as.character(pred)
    ids <- object$id[test.idx]
    pops <- object$pop[test.idx]
    names(pred) <- ids
    
    #Make output
    ids <- ids[1:length(ids) %% 2 == 0]
    pops <- pops[1:length(pops) %% 2 == 0]
    out <- rep(NA, times=length(ids)*8)
    for(i in 1:length(ids)){
      haps <- pred[which(names(pred) == ids[i])]
      haps <- unlist(c(block.info,pops[i],ids[i],haps[1],haps[2]))
      haps <- as.vector(haps)
      out[(i*8 - 7):(i*8)] <- haps
    }
    
    #Return output
    return(out)
  }
  
  # Tuning for svm ----------------------------------------------------------------------------------
  if(tune == TRUE){
    if(verbose == TRUE){
      cat("\nPerforming 5-fold cross-validation... ")
    }
    test.old <- test
    train.old <- train
    param.old <- param
    groups <- as.character(1:5)
    train.group <- sample(x = groups, size = length(train.old), replace = TRUE)
    acc <- matrix(data = NA, nrow = length(param.old$cost)*length(param.old$gamma), ncol = 3)
    acc <- as.data.frame(acc)
    colnames(acc) <- c("cost", "gamma", "accuracy")
    l <- 1
    for(j in 1:length(param.old$cost)){
      for(k in 1:length(param.old$gamma)){
        perc <- rep(NA, times = 5)
        param$cost <- param.old$cost[j]
        param$gamma <- param.old$gamma[k]
        for(g in 1:length(groups)){
          train <- train.old[which(train.group != groups[g])]
          train.idx <- which(object$id %in% train)
          train.pop <- object$pop[train.idx]
          y <- object$pop[train.idx]
          y <- as.factor(y)
          test <- train.old[which(train.group == groups[g])]
          test.idx <- which(object$id %in% test)
          if(Sys.info()["sysname"] == "Windows"){
            cl <- makeCluster(ncores)
            results <- unlist(parLapply(cl = cl, fun = blockfun, X = 1:nrow(blocks)))
            stopCluster(cl)
          }else{
            results <- mclapply(FUN = blockfun, X = 1:nrow(blocks), mc.cores = ncores)
          }
          results <- data.frame(matrix(unlist(results), ncol=8, byrow=TRUE), stringsAsFactors = F)
          colnames(results) <- c("BLOCK","CHR","BP1","BP2","POP","ID","HAP1","HAP2")
          results$BP1 <- as.numeric(results$BP1)
          results$BP2 <- as.numeric(results$BP2)
          perc[g] <- length(which(results$POP == results$HAP1)) + length(which(results$POP == results$HAP2))
          perc[g] <- 100*perc[g]/(2*nrow(results))
        }
        acc$cost[l] <- param$cost
        acc$gamma[l] <- param$gamma
        acc$accuracy[l] <- mean(perc)
        l <- l + 1
      }
    }
  }
  
  
  # Check whether ancestries should be computed-------------------------------------------------------
  comp <- TRUE
  if(tune == TRUE){
    comp <- FALSE
    results <- acc
    if(verbose == TRUE){
      cat("Done.\n")
    }
  }
  
  # Compute ancestry----------------------------------------------------------------------------------
  if(comp == TRUE){
    if(verbose == TRUE){
      cat("\nPredicting ancestry of haplotypes... ")
    }
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      results <- unlist(parLapply(cl = cl, fun = blockfun, X = 1:nrow(blocks)))
      stopCluster(cl)
    }else{
      results <- mclapply(FUN = blockfun, X = 1:nrow(blocks), mc.cores = ncores)
    }
    if(verbose == TRUE){
      cat("Done.\n")
      cat("Assembling results... ")
    }
    results <- data.frame(matrix(unlist(results), ncol=8, byrow=TRUE), stringsAsFactors = F)
    colnames(results) <- c("BLOCK","CHR","BP1","BP2","POP","ID","HAP1","HAP2")
    results$BP1 <- as.numeric(results$BP1)
    results$BP2 <- as.numeric(results$BP2)
    if(verbose == TRUE){
      cat("Done.\n")
    }
  }
  
  # Return results------------------------------------------------------------------------------------
  return(results)
  
  
}
