#Function: ghap.ibd
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: IBD estimates

ghap.ibd <- function(
  object,
  pairlist,
  freq,
  mafcut=0.05,
  refsize=10000,
  batchsize=NULL,
  ncores=1,
  verbose=TRUE
){
  
  # Check if input is a valid GHap object --------------------------------------
  obtype <- c("GHap.phase","GHap.plink")
  if(inherits(object, obtype) == FALSE){
    stop("\nInput must be a valid GHap object.")
  }
  
  # Check if pairs exist -------------------------------------------------------
  paircols <- length(which(colnames(pairlist) %in% c("ID1","ID2")))
  if(paircols != 2){
    stop("\nInvalid or missing columns in your pairs list.")
  }
  id.in <- unique(c(pairlist$ID1,pairlist$ID2))
  uniqids <- length(which(id.in %in% object$id))
  if(uniqids != length(id.in)){
    stop("\nSome of the provided ids were not found")
  }
  
  # Prune out monomorphic markers ----------------------------------------------
  maf <- pmin(freq,1-freq)
  freq <- freq[which(maf > mafcut)]
  rm(maf)
  
  # Check if markers exist -----------------------------------------------------
  var.n <- length(freq)
  var.in <- names(freq)
  if(sum(var.in %in% object$marker) != var.n){
    stop("The vector of allele frequencies contains unknown markers")
  }
  
  # Compute marker table -------------------------------------------------------
  N <- refsize
  A <- round(N*freq, digits <- 0)
  B <- N - A
  A1 <- (A-1)/A
  A2 <- (A-2)/A
  A3 <- (A-3)/A
  A12 <- A1*A2
  A123 <- A12*A3
  B1 <- (B-1)/B
  B2 <- (B-2)/B
  B3 <- (B-3)/B
  B12 <- B1*B2
  B123 <- B12*B3
  A1B1 <- A1*B1
  N1 <- N/(N-1)
  N2 <- N/(N-2)
  N3 <- N/(N-3)
  N12 <- N1*N2
  N123 <- N12*N3
  A1B1N123 <- A1*B1*N123
  p <- A/N
  q <- B/N
  p2 <- p^2
  p3 <- p^3
  p4 <- p^4
  q2 <- q^2
  q3 <- q^3
  q4 <- q^4
  P00 <- 2*p2*q2*A1B1N123
  P10 <- 4*p3*q*A12*N123 + 4*p*q3*B12*N123
  P20 <- p4*A123*N123 + q4*B123*N123 + 4*p2*q2*A1B1N123
  P01 <- rep(0, times = var.n)
  P11 <- 2*p2*q*A1*N12 + 2*p*q2*B1*N12
  P21 <- p3*A12*N12 + q3*B12*N12 + p2*q*A1*N12 + p*q2*B1*N12
  P22 <- rep(1, times = var.n)
  names(A) <- var.in
  names(B) <- var.in
  names(P00) <- var.in
  names(P10) <- var.in
  names(P20) <- var.in
  names(P01) <- var.in
  names(P11) <- var.in
  names(P21) <- var.in
  names(P22) <- var.in
  rm(A1,A2,A3,A12,A123,B1,B2,B3,B12,B123,A1B1,
     N1,N2,N3,A1B1N123,p2,p3,p4,q2,q3,q4)
  
  # Generate batch index -------------------------------------------------------
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(var.n/10)
  }
  if(batchsize > var.n){
    batchsize <- var.n
  }
  id1 <- seq(1,var.n,by=batchsize)
  id2 <- (id1+batchsize)-1
  id1 <- id1[id2<=var.n]
  id2 <- id2[id2<=var.n]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,var.n)
  if(id1[length(id1)] > var.n){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  # Log message ----------------------------------------------------------------
  if(verbose == TRUE){
    cat("Processing ", var.n, " variants in ", length(id1), " batches.\n", sep="")
    cat("Variants not present in 'freq' will be ignored.\n")
  }
  
  # IBD main function ----------------------------------------------------------
  ibdfun <- function(k){
    id1 <- Ztmp[,pairlist$ID1[k]]
    id2 <- Ztmp[,pairlist$ID2[k]]
    obs <- names(id1)[which(is.na(id1) == FALSE & is.na(id2) == FALSE)]
    id1 <- id1[obs]
    id2 <- id2[obs]
    ibs0 <- length(which(id1 == 2 & id2 == 0)) + length(which(id1 == 0 & id2 == 2))
    ibs1 <- length(which(id1 == 2 & id2 == 1)) + length(which(id1 == 0 & id2 == 1)) + 
      length(which(id1 == 1 & id2 == 2)) + length(which(id1 == 1 & id2 == 0))
    ibs2 <- length(which(id1 == 2 & id2 == 2)) + length(which(id1 == 1 & id2 == 1)) +
      length(which(id1 == 0 & id2 == 0))
    E <- c(sum(P00[obs]), sum(P10[obs]), sum(P20[obs]), sum(P11[obs]),
           sum(P21[obs]), sum(P22[obs]))
    return(c(ibs0,ibs1,ibs2,E))
  }
  
  # IBD iterate function -------------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  sumvariants <- 0
  pairlist$IBS0 <- 0
  pairlist$IBS1 <- 0
  pairlist$IBS2 <- 0
  pairlist$E00 <- 0
  pairlist$E10 <- 0
  pairlist$E20 <- 0
  pairlist$E11 <- 0
  pairlist$E21 <- 0
  pairlist$E22 <- 0
  for(i in 1:length(id1)){
    idx <- id1[i]:id2[i]
    Ztmp <- ghap.slice(object = object,
                       ids = id.in,
                       variants = var.in[idx],
                       index = FALSE,
                       unphase = TRUE,
                       impute = FALSE,
                       ncores = ncores)
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      clusterEvalQ(cl, library(Matrix))
      varlist <- list("Ztmp","pairlist","P00","P10","P20","P11","P21","P22")
      clusterExport(cl = cl, varlist = varlist, envir=environment())
      out <- unlist(parLapply(cl = cl, fun = ibdfun, X = 1:nrow(pairlist)))
      stopCluster(cl)
    }else{
      out <- mclapply(FUN = ibdfun, X = 1:nrow(pairlist), mc.cores = ncores)
      out <- matrix(data = unlist(out), nrow = nrow(pairlist), ncol = 9, byrow = TRUE)
    }
    pairlist$IBS0 <- pairlist$IBS0 + out[,1]
    pairlist$IBS1 <- pairlist$IBS1 + out[,2]
    pairlist$IBS2 <- pairlist$IBS2 + out[,3]
    pairlist$E00 <- pairlist$E00 + out[,4]
    pairlist$E10 <- pairlist$E10 + out[,5]
    pairlist$E20 <- pairlist$E20 + out[,6]
    pairlist$E11 <- pairlist$E11 + out[,7]
    pairlist$E21 <- pairlist$E21 + out[,8]
    pairlist$E22 <- pairlist$E22 + out[,9]
    if(verbose == TRUE){
      sumvariants <- sumvariants + length(idx)
      cat(sumvariants, "variants processed.\r")
    }
  }
  pairlist$Z0 <- pairlist$IBS0/pairlist$E00
  pairlist$Z1 <- (pairlist$IBS1 - pairlist$Z0*pairlist$E10)/pairlist$E11
  pairlist$Z2 <- (pairlist$IBS2 - pairlist$Z0*pairlist$E20 - pairlist$Z1*pairlist$E21)/pairlist$E22
  
  # Correct Z values -------------------------------------------------------------
  idx <- which(pairlist$Z0 < 0)
  pairlist$Z0[idx] <- 0
  idx <- which(pairlist$Z1 < 0)
  pairlist$Z1[idx] <- 0
  idx <- which(pairlist$Z2 < 0)
  pairlist$Z2[idx] <- 0
  idx <- which(pairlist$Z0 > 1)
  pairlist$Z0[idx] <- 1
  pairlist$Z1[idx] <- 0
  pairlist$Z2[idx] <- 0
  idx <- which(pairlist$Z1 > 1)
  pairlist$Z0[idx] <- 0
  pairlist$Z1[idx] <- 1
  pairlist$Z2[idx] <- 0
  idx <- which(pairlist$Z2 > 1)
  pairlist$Z0[idx] <- 0
  pairlist$Z1[idx] <- 0
  pairlist$Z2[idx] <- 1
  Zsum <- pairlist$Z0 + pairlist$Z1 + pairlist$Z2
  pairlist$Z0 <- pairlist$Z0/Zsum
  pairlist$Z1 <- pairlist$Z1/Zsum
  pairlist$Z2 <- pairlist$Z2/Zsum
  
  # Organize output -----------------------------------------------------------
  if(inherits(object, "GHap.phase")){
    myids <- 1:(2*object$nsamples)
    myids <- myids[myids %% 2 == 1]
    names(myids) <- object$id[myids]
  }else{
    myids <- 1:object$nsamples
    names(myids) <- object$id
  }
  m <- pairlist$IBS0 + pairlist$IBS1 + pairlist$IBS2
  pairlist$PERC <- pairlist$IBS2/m
  pairlist$DST <- (pairlist$IBS2 + 0.5*pairlist$IBS1)/m
  pairlist$PI_HAT <- pairlist$Z2 + 0.5*pairlist$Z1
  pairlist$POP1 <- object$pop[myids[pairlist$ID1]]
  pairlist$POP2 <- object$pop[myids[pairlist$ID2]]
  keep <- c('POP1','ID1','POP2','ID2',
            'IBS0','IBS1','IBS2','PERC','DST',
            'Z0','Z1','Z2','PI_HAT')
  pairlist <- pairlist[,keep]
  return(pairlist)
  
  
}
