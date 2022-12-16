#Function: ghap.slice
#License: GPLv3 or later
#Modification date: 16 Dec 2022
#Written by: Yuri Tani Utsunomiya, Adam Taiti Harth Utsunomiya
#Contact: ytutsunomiya@gmail.com, adamtaiti@gmail.com
#Description: Get a slice of a GHap object

ghap.slice <- function(
    object,
    ids,
    variants,
    index=FALSE,
    transposed=FALSE,
    sparse=TRUE,
    unphase=FALSE,
    impute=FALSE
){
  
  # Check object type ----------------------------------------------------------
  obtypes <- c("GHap.phase","GHap.haplo","GHap.plink")
  if(inherits(object, obtypes) == FALSE){
    stop("\nInput data must be a valid GHap object (phase, haplo or plink).")
  }
  
  # Build indices --------------------------------------------------------------
  # iidx = individual index (1 per individual)
  # hidx = haplotype index (2 per individual, only for mode 0 and 1)
  # oidx = output haplotype indices (1 or 2 per individual, only for index = T)
  # vidx = variant index
  if(inherits(object, "GHap.phase")){
    binfile <- object$phase
    mode <- object$mode
    nvars <- object$nmarkers
    nids <- object$nsamples
    phased <- 2
    imp <- 0
    if(index == TRUE){
      idlab <- object$id[ids]
      iidx <- 1:object$nsamples
      names(iidx) <- object$id[1:length(object$id) %% 2 == 0]
      iidx <- iidx[idlab[duplicated(idlab) == FALSE]]
      x <- rep((iidx*2)-1, each=2)
      x[1:length(x) %% 2 == 0] <- x[1:length(x) %% 2 == 0] + 1
      y <- 1:length(x)
      names(y) <- x
      oidx <- y[as.character(ids)]
      hidx <- rep(2*iidx, each=2)
      hidx[1:length(x) %% 2 == 1] <- hidx[1:length(x) %% 2 == 1] - 1
      vidx <- variants
      names(vidx) <- object$marker[vidx]
    }else{
      iidx <- 1:object$nsamples
      names(iidx) <- object$id[1:length(object$id) %% 2 == 0]
      iidx <- iidx[ids]
      hidx <- rep(2*iidx, each=2)
      hidx[1:length(hidx) %% 2 == 1] <- hidx[1:length(hidx) %% 2 == 0]-1
      vidx <- 1:object$nmarkers
      names(vidx) <- object$marker
      vidx <- vidx[variants]
    }
    nrows <- length(vidx)
    ncols <- 2*length(iidx)
    rowlab <- names(vidx)
    collab <- rep(names(iidx), each=2)
    brow <- mode != 2
  }
  if(inherits(object, "GHap.plink")){
    binfile <- object$plink
    mode <- 3
    nvars <- object$nmarkers
    nids <- object$nsamples
    phased <- 1
    imp <- as.integer(impute)
    if(index == TRUE){
      iidx <- ids
      names(iidx) <- object$id[ids]
      hidx <- iidx
      vidx <- variants
      names(vidx) <- object$marker[vidx]
    }else{
      iidx <- 1:object$nsamples
      names(iidx) <- object$id
      iidx <- iidx[ids]
      hidx <- iidx
      vidx <- 1:object$nmarkers
      names(vidx) <- object$marker
      vidx <- vidx[variants]
    }
    nrows <- length(vidx)
    ncols <- length(iidx)
    rowlab <- names(vidx)
    collab <- names(iidx)
    brow <- TRUE
  }
  if(inherits(object, "GHap.haplo")){
    binfile <- object$genotypes
    mode <- 3
    nvars <- object$nalleles
    nids <- object$nsamples
    phased <- 1
    imp <- 0
    if(index == TRUE){
      iidx <- ids
      names(iidx) <- object$id[ids]
      hidx <- iidx
      vidx <- variants
      names(vidx) <- object$allele[vidx]
    }else{
      iidx <- 1:object$nsamples
      names(iidx) <- object$id
      iidx <- iidx[ids]
      hidx <- iidx
      vidx <- variants
      names(vidx) <- paste(object$block[vidx],object$bp1[vidx],
                           object$bp2[vidx],object$allele[vidx], sep="_")
    }
    nrows <- length(vidx)
    ncols <- length(iidx)
    rowlab <- names(vidx)
    collab <- names(iidx)
    brow <- TRUE
  }
  
  # Get slice of data ----------------------------------------------------------
  X <- .sliceCpp(binfile, mode, nvars, nids, phased, imp, iidx, hidx, vidx)

  # Organize output matrix -----------------------------------------------------
  if(transposed == FALSE){
    X <- Matrix(data = X, nrow = nrows, ncol = ncols,
                byrow = brow, sparse = TRUE)
    rownames(X) <- rowlab
    colnames(X) <- collab
    if(unphase == TRUE & phased == 2){
      cols <- 1:ncol(X) %% 2
      X <- X[,which(cols == 1)] + X[,which(cols == 0)]
    }else if(index == TRUE & phased == 2){
      X <- X[,oidx]
    }
  }else{
    X <- Matrix(data = X, nrow = ncols, ncol = nrows, 
                byrow = (!brow), sparse = TRUE)
    rownames(X) <- collab
    colnames(X) <- rowlab
    if(unphase == TRUE & phased == 2){
      rows <- 1:nrow(X) %% 2
      X <- X[which(rows == 1),] + X[which(rows == 0),]
    }else if(index == TRUE & phased == 2){
      X <- X[oidx,]
    }
  }
  
  # Check if output should be sparse -------------------------------------------
  if(sparse == FALSE){
    X <- as.matrix(X)
  }

  # Return X -------------------------------------------------------------------
  return(X)
  
}
