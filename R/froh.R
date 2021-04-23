#Function: ghap.froh
#License: GPLv3 or later
#Modification date: 07 Feb 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Computation of genomic inbreeding from runs of homozygosity

ghap.froh<-function(
  phase,
  roh,
  rohsizes=c(1,2,4,8,16),
  only.active.markers = TRUE,
  ncores=1
){

  # Check if phase is a GHap.phase object-------------------------------------------------------------
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Check if inactive markers should be reactivated---------------------------------------------------
  if(only.active.markers == FALSE){
    phase$marker.in <- rep(TRUE,times=phase$nmarkers)
    phase$nmarkers.in <- length(which(phase$marker.in))
  }
  
  # Compute genome size-------------------------------------------------------------------------------
  mkrdist <- diff(phase$bp[which(phase$marker.in)])
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
  if(Sys.info()["sysname"] == "Windows"){
    inb <- lapply(FUN = rohsum, X = 1:nrow(id))
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