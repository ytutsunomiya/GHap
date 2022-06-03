#Function: ghap.froh
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Computation of genomic inbreeding from runs of homozygosity

ghap.froh<-function(
  object,
  roh,
  rohsizes=c(1,2,4,8,16),
  only.active.markers = TRUE,
  ncores=1
){
  
  # Check if input is a valid GHap object-------------------------------------------------------------
  obtype <- c("GHap.phase","GHap.plink")
  if(inherits(object, obtype) == FALSE){
    stop("\nInput must be a valid GHap object.")
  }
  
  # Check if inactive markers should be reactivated---------------------------------------
  if(only.active.markers == FALSE){
    object$marker.in <- rep(TRUE,times=object$nmarkers)
    object$nmarkers.in <- length(which(object$marker.in))
  }
  
  # Compute genome size-------------------------------------------------------------------------------
  mkrdist <- diff(object$bp[which(object$marker.in)])
  mkrdist <- mkrdist[which(mkrdist > 0)]
  genome <- sum(as.numeric(mkrdist))
  
  # ROH sum function----------------------------------------------------------------------------------
  rohsum <- function(i){
    sroh <- rep(NA,times=length(rohsizes))
    for(j in 1:length(rohsizes)){
      sroh[j] <- sum(roh$LENGTH[which(roh$ID == id[i,2] & roh$LENGTH >= rohsizes[j]*1e+6)])
    }
    out <- c(id[i,1],id[i,2],sroh/genome)
    return(out)
  }
  
  # Compute genomic inbreeding------------------------------------------------------------------------
  id <- unique(roh[,c("POP","ID")])
  ncores <- min(c(detectCores(), ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cl <- makeCluster(ncores)
    inb <- unlist(parLapply(cl = cl, fun = rohsum, X = 1:nrow(id)))
    stopCluster(cl)
  }else{
    inb <- mclapply(FUN = rohsum, X = 1:nrow(id), mc.cores = ncores)
  }
  froh <- matrix(data = unlist(inb), nrow = nrow(id), ncol = 2+length(rohsizes), byrow = TRUE)
  froh <- as.data.frame(froh, stringsAsFactors = FALSE)
  for(k in 3:ncol(froh)){
    froh[,k] <- as.numeric(froh[,k])
  }
  colnames(froh) <- c("POP","ID",paste0("FROH",rohsizes))
  return(froh)
  
}
