#Function: ghap.simpheno
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Simulate phenotypes

ghap.simpheno<-function(
  object,
  h2,
  r2=0,
  nrep=1,
  balanced=TRUE,
  seed=NULL,
  only.active.samples=TRUE,
  ncores=1,
  verbose=TRUE
){
  
  # Sanity check for input objects ---------------------------------------------
  obtype <- c("GHap.phase","GHap.plink")
  fac <- c(2,1)
  names(fac) <- obtype
  if(inherits(object, obtype) == FALSE){
    stop("\nInput must be a valid GHap.phase or GHap.plink object.")
  }
  if(nrep > 1 & r2 <= 0){
    stop("\nRepeatability must be greater than 0")
  }
  h2tot <- sum(h2)
  if(h2tot < 0 | h2tot > 1){
    stop("\nThe sum of variant heritabilities must be between 0 and 1.")
  }
  
  # Check if inactive samples should be reactived ------------------------------
  if(only.active.samples == FALSE){
    object$id.in <- rep(TRUE,times=fac[class(object)]*object$nsamples)
    object$nsamples.in <- length(which(object$id.in))/fac[class(object)]
  }
  ids <- unique(object$id[object$id.in])
  object <- ghap.subset(object = object, ids = ids,
                        variants = names(h2), index = FALSE, verbose = FALSE)
  vidx <- which(object$marker.in)
  h2 <- h2[object$marker[vidx]]
  
  # Compute allele frequencies -------------------------------------------------
  p <- ghap.freq(object = object, type = "A1", only.active.samples = TRUE,
                 only.active.markers = TRUE)
  p[which(p < 0.01)] <- 0.01
  p[which(p > 0.99)] <- 0.99
  
  # Make additive effects ------------------------------------------------------
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }
  a <- sqrt(h2/(2*p*(1-p)))
  a <- a*sample(x = c(-1,1), size = length(a), replace = TRUE)
  
  # Compute true breeding value ------------------------------------------------
  score <- data.frame(MARKER = object$marker[vidx], CHR = object$chr[vidx],
                      BP = object$bp[vidx], ALLELE = object$A1[vidx],
                      SCORE = a, CENTER = 0, SCALE = 1)
  scoreid <- ghap.profile(object = object, score = score, verbose = verbose,
                          ncores = ncores, only.active.samples = TRUE)
  u <- as.numeric(scoreid$SCORE)
  
  # Simulate repeatability -----------------------------------------------------
  if(nrep > 1){
    if(is.null(seed) == FALSE){
      set.seed(seed)
    }
    p <- rnorm(length(u),sd=sqrt(r2))
  }
  
  # Simulate residuals ---------------------------------------------------------
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }
  if(nrep > 1){
    e <- rnorm(nrep*length(u), sd=sqrt(1-h2tot-r2))
  }else{
    e <- rnorm(length(u), sd=sqrt(1-h2tot))
  }
  
  # Simulate phenotypes --------------------------------------------------------
  if(nrep > 1){
    y <- rep(u,each=nrep) + rep(p,each=nrep) + e
  }else{
    y <- u + e
  }
  
  # Assemble output ------------------------------------------------------------
  if(nrep > 1){
    df <- data.frame(POP = rep(scoreid$POP, each=nrep),
                     ID = rep(scoreid$ID, each=nrep),
                     PHENO = y,
                     TBV = rep(u, each=nrep),
                     REP = rep(p, each=nrep),
                     RESIDUAL = e)
    if(balanced == FALSE){
      keep <- NULL
      for(i in unique(df$ID)){
        nkeep <- ceiling(runif(n = 1, min = 0, max = nrep))
        idx <- which(df$ID == i)
        idx <- sample(x = idx, size = nkeep)
        keep <- c(keep, idx)
      }
      df <- df[sort(keep),]
    }
  }else{
    df <- data.frame(POP = scoreid$POP,
                     ID = scoreid$ID,
                     PHENO = y,
                     TBV = u,
                     RESIDUAL = e)
  }
  
  #Return output
  return(df)
  
}
