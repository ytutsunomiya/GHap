#Function: ghap.subset
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Subset a GHap object

ghap.subset <- function(
  object,
  ids,
  variants,
  index = FALSE,
  verbose = TRUE
){
  
  # Check if object is a valid GHap object -------------------------------------
  obtypes <- c("GHap.phase","GHap.haplo","GHap.plink")
  if(inherits(object, obtypes) == FALSE){
    stop("\nInput data must be a valid GHap object (phase, haplo or plink).")
  }
  
  #Log message -----------------------------------------------------------------
  if(verbose == TRUE){
    cat("\nSubsetting ", length(ids), " individuals and ",
        length(variants), " variants... ",sep="")
  }
  
  # Retrieve indices -----------------------------------------------------------
  if(index == TRUE){
    iidx <- ids
    vidx <- variants
  }else{
    iidx <- which(object$id %in% ids)
    vidx <- which(object$marker %in% variants)
    if(identical(sort(unique(object$id[iidx])),sort(ids)) == FALSE){
      stop("Some of the provided ids were not found")
    }
    if(length(vidx) < length(variants)){
      stop("Some of the provided variants were not found")
    }
  }
  
  # Subset ids -----------------------------------------------------------------
  if(max(iidx) > 2*object$nsamples & inherits(object, "GHap.phase")){
    stop("Some of the provided ids are out of range")
  }else if(max(iidx) > object$nsamples & inherits(object, "GHap.phase") == FALSE){
    stop("Some of the provided ids are out of range")
  }else{
    ids <- rep(FALSE, times = length(object$id.in))
    ids[iidx] <- TRUE
    object$id.in <- ids
    object$nsamples.in <- length(unique(object$id[iidx]))
  }
  
  # Subset variants ------------------------------------------------------------
  if(inherits(object, "GHap.haplo")){
    if(max(vidx) > object$nalleles){
      stop("Some of the provided alleles are out of range")
    }else{
      variants <- rep(FALSE, times = length(object$allele.in))
      variants[vidx] <- TRUE
      object$allele.in <- variants
      object$nalleles.in <- length(vidx)
    }
  }else{
    if(max(vidx) > object$nmarkers){
      stop("Some of the provided markers are out of range")
    }else{
      variants <- rep(FALSE, times = length(object$marker.in))
      variants[vidx] <- TRUE
      object$marker.in <- variants
      object$nmarkers.in <- length(vidx)
    }
  }
  
  #Log message
  if(verbose == TRUE){
    cat("Done.\n")
    cat("Final data contains ", object$nsamples.in, " individuals and ",
        length(vidx), " variants.\n\n",sep="")
  }
  
  #Output object
  return(object)
  
}
