#Function: ghap.remlci
#License: GPLv3 or later
#Modification date: 09 Nov 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: confidence intervals for functions of
#             variance components obtained from AI-REML

ghap.remlci <- function(
  fun,
  vcp,
  ai,
  n = 10000,
  conf.level = 0.95,
  include.samples = FALSE,
  ncores = 1
){
  
  # Compute cholesky factor
  L <- t(chol(as.matrix(ai)))
  
  # Initialize sampling function
  samp <- function(i){
    d <- rnorm(n = nrow(L))
    x <- as.numeric(vcp[,1] + L%*%d)
    return(fun(x))
  }
  
  # Make sampling
  ncores <- min(c(detectCores(), ncores))
  if(ncores == 1){
    y <- lapply(FUN = samp, X = 1:n)
    y <- unlist(y)
  }else{
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      clusterEvalQ(cl, library(Matrix))
      varlist <- list("fun","ai","vcp","L")
      clusterExport(cl = cl, varlist = varlist, envir=environment())
      y <- unlist(parLapply(cl = cl, fun = samp, X = 1:n))
      stopCluster(cl)
    }else{
      y <- mclapply(FUN = samp, X = 1:n, mc.cores = ncores)
      y <- unlist(y)
    }
  }
  
  # Make output
  results <- NULL
  results$stde <- sd(y)
  alpha <- 1-conf.level
  results$ci <- quantile(x = y, probs = c(alpha/2, 1-(alpha/2)))
  if(include.samples == TRUE){
    results$samples <- y
  }
  return(results)
}
