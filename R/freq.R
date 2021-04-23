#Function: ghap.freq
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute marker allele frequencies

ghap.freq <- function(
  phase,
  type="maf",
  only.active.samples=TRUE,
  only.active.markers=TRUE,
  batchsize=NULL,
  ncores=1,
  verbose=TRUE
){
  
  #Check if phase is a GHap object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  #Check if inactive markers and samples should be reactived
  if(only.active.markers == FALSE){
    phase$marker.in <- rep(TRUE,times=phase$nmarkers)
    phase$nmarkers.in <- length(which(phase$marker.in))
  }
  if(only.active.samples == FALSE){
    phase$id.in <- rep(TRUE,times=2*phase$nsamples)
    phase$nsamples.in <- length(which(phase$id.in))/2
  }
  
  #Check if type is valid
  if(type %in% c("maf","A0","A1") == FALSE){
    stop("Argument type must be 'maf', 'A0' or 'A1'.")
  }
  
  # Initialize lookup table
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
  
  #Generate batch index
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
  
  #Windows warning
  ncores <- min(c(detectCores(), ncores))
  if(Sys.info()["sysname"] == "Windows" & ncores > 1 & verbose == TRUE){
    cat("\nParallelization not supported yet under Windows (using a single core).\n")
  }
  
  #Frequency function
  freq.fun <- function(i){
    x <- X[i,]
    p <- sum(x)/(2*phase$nsamples.in)
    return(p)
  }
  
  #Frequency calculation
  ids.in <- which(phase$id.in)
  snps.in <- which(phase$marker.in)
  freq <- rep(NA, times=phase$nmarkers.in)
  for(i in 1:length(id1)){
    X <- ghap.pslice(phase = phase,
                     ids = ids.in,
                     markers = snps.in[id1[i]:id2[i]],
                     index = TRUE,
                     lookup = lookup,
                     ncores = ncores)
    #Compute blocks
    if(Sys.info()["sysname"] == "Windows"){
      p <- unlist(lapply(FUN = freq.fun, X = 1:nrow(X)))
    }else{
      p <- unlist(mclapply(FUN = freq.fun, X = 1:nrow(X), mc.cores = ncores))
    } 
    freq[id1[i]:id2[i]] <- unlist(p)
  }
  
  #Results
  names(freq) <- phase$marker[snps.in]
  if(type == "maf"){
    freq <- pmin(freq,1-freq)
  }else if(type == "A0"){
    freq <- 1-freq
  }
  return(freq)
  
}