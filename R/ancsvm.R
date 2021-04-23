#Function: ghap.ancsvm
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Predict ancestry of haplotypes using machine learning

ghap.ancsvm <- function(
  phase,
  blocks,
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
  
  # Check if phase is a GHap.phase object-------------------------------------------------------------
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
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
  
  # Map test samples----------------------------------------------------------------------------------
  if(is.null(test) == TRUE){
    test.idx <- which(phase$id.in == TRUE)
  }else{
    test.idx <- which(phase$id %in% test & phase$id.in == TRUE)
  }
  
  # Map training samples -----------------------------------------------------------------------------
  if(is.null(train) == TRUE){
    train.idx <- which(phase$id.in == TRUE)
  }else{
    train.idx <- which(phase$id %in% train & phase$id.in == TRUE)
  }
  y <- phase$pop[train.idx]
  y <- as.factor(y)
  
  # Map parameters to use-----------------------------------------------------------------------------
  param <- list(cost = cost, gamma = gamma, tune = tune)
  if(is.null(param$gamma) == TRUE){
    param$gamma <- "1/blocksize"
  }
  
  # Log message of parameters-------------------------------------------------------------------------
  if(verbose == TRUE){
    printparams <- paste(names(param), "=", param, collapse=", ")
    cat("\nUsing svm with parameters:\n[", printparams, "]\n", sep="")
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
  
  # Initialize block iteration function---------------------------------------------------------------
  blockfun <- function(b){
    
    #Get block info
    block.info <- blocks[b, c("BLOCK","CHR","BP1","BP2")]
    
    #SNPs in the block
    snps <- which(phase$chr == block.info$CHR &
                    phase$bp >= block.info$BP1 &
                    phase$bp <= block.info$BP2 &
                    phase$marker.in == TRUE)
    blocksize <- length(snps)
    
    #Build model matrices
    Mtst <- ghap.pslice(phase = phase, ids = test.idx, markers = snps,
                        index = TRUE, lookup = lookup, verbose = FALSE)
    Mref <- ghap.pslice(phase = phase, ids = train.idx, markers = snps,
                        index = TRUE, lookup = lookup, verbose = FALSE)
    Mtst <- t(Mtst)
    Mref <- t(Mref)
    
    #Model training
    if(param$gamma == "1/blocksize"){
      gamma <- 1/blocksize
    }else{
      gamma <- param$gamma
    }
    model <- svm(y = y, x = Mref, kernel = "radial", gamma = gamma, cost = param$cost)
    pred <- predict(model, Mtst)
    pred <- as.character(pred)
    ids <- phase$id[test.idx]
    pops <- phase$pop[test.idx]
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
            train.idx <- which(phase$id %in% train)
            train.pop <- phase$pop[train.idx]
            y <- phase$pop[train.idx]
            y <- as.factor(y)
            test <- train.old[which(train.group == groups[g])]
            test.idx <- which(phase$id %in% test)
            if(Sys.info()["sysname"] == "Windows"){
              results <- lapply(FUN = blockfun, X = 1:nrow(blocks))
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
      results <- lapply(FUN = blockfun, X = 1:nrow(blocks))
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