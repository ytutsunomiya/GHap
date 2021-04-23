#Function: ghap.kinv
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Inverse kinship matrix

ghap.kinv<-function(
  kinship,           #Kinship matrix
  method="nearPD",   #Inversion method
  ncores=1,          #Number of cores to use
  proven=NULL,       #Names of proven subjects to use in APY
  verbose=TRUE
){
  
  #Check if kinship matrix is symmetrical
  if(identical(colnames(kinship),rownames(kinship)) == FALSE){
    stop("Names in rows and columns must be identical.")
  }
  
  #Check if method is valid
  if(method %in% c("common","nearPD","APY") == FALSE){
    stop('Method must be one of "common", "nearPD" and "APY" ')
  }
  
  #Inversion
  if(method == "common"){
    Kinv <- solve(kinship)
    colnames(Kinv) <- colnames(kinship)
    rownames(Kinv) <- rownames(kinship)
  }else if(method == "nearPD"){
    Kinv <- solve(nearPD(kinship)$mat)
  }else if(method == "APY"){
    if(is.null(proven)){
      stop('Names of proven subjects must be given for method "APY"')
    }
    if(length(which(proven %in% colnames(kinship))) != length(proven)){
      stop('All declared proven subjects must be present in the kinship matrix')
    }
    nonproven <- colnames(kinship)[colnames(kinship) %in% proven == FALSE]
    Kcci <- solve(kinship[proven,proven])
    colnames(Kcci) <- proven
    rownames(Kcci) <- proven
    Knc <- kinship[nonproven,proven]
    dFUN <- function(i){
      Kic <- rbind(kinship[i,proven])
      return(as.vector(kinship[i,i] - Kic%*%tcrossprod(Kcci,Kic)))
    }
    ncores <- min(c(detectCores(),ncores))
    if(Sys.info()["sysname"] == "Windows"){
      if(ncores > 1 & verbose == TRUE){
        cat("\nParallelization not supported yet under Windows (using a single core).\n")
      }
      Jnni <- Diagonal(x = 1/unlist(lapply(X = nonproven,FUN = dFUN)))
    }else{
      Jnni <- Diagonal(x = 1/unlist(mclapply(X = nonproven,FUN = dFUN, mc.cores = ncores)))
    }
    
    rownames(Jnni) <- nonproven
    colnames(Jnni) <- rownames(Jnni)
    KcciKcn <- tcrossprod(Kcci,Knc)
    JnniKncKcci <- Jnni%*%Knc%*%Kcci
    Kinv <- rbind(cbind(Kcci + KcciKcn%*%JnniKncKcci,-KcciKcn%*%Jnni),
                  cbind(-JnniKncKcci,Jnni))
    Kinv <- as(Kinv,"dgCMatrix")
    
  }
  return(Kinv[rownames(kinship),colnames(kinship)])
  
}