#Function: ghap.simpheno
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Simulate phenotypes

ghap.simpheno<-function(
  haplo,
  kinship,
  h2,
  g2,
  r2=0,
  nrep=1,
  balanced=TRUE,
  major=NULL,
  seed=NULL,
  ncores=1
){
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if kinship matrix is symmetrical
  if(identical(colnames(kinship),rownames(kinship)) == FALSE){
    stop("Names in rows and columns must be identical.")
  }
  
  #Check if names in the kinship matrix match with the GHap.haplo object
  if (length(which(colnames(kinship) %in% haplo$id)) != ncol(kinship)) {
    stop("All ids in the kinship matrix must be present in the GHap.haplo object.")
  }
  
  # Generate lookup table
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
  lookup <- sapply(lookup, function(i){intToUtf8(rev(utf8ToInt(i)))})
  
  # Simulate uncorrelated random effects
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }
  g <- rnorm(ncol(kinship),sd=sqrt(h2 - sum(g2*h2)))
  
  # Make random effects correlated by kinship
  # Now g = polygenic effect (sum of genome-wide haplotype effects) - major effect
  g <- crossprod(chol(nearPD(kinship)$mat),g)
  g <- as.vector(g)
  
  # Simulate major haplotypes
  if(is.null(major) == FALSE){
    cond <- which(g2 < 0 | g2 > 1)
    if(length(cond) > 0){
      stop("Argument g2 must be between zero and one.")
    }
    if(length(g2) != length(major)){
      stop("Vectors g2 and major must have equal length.")
    }
    if(sum(g2) > 1){
      stop("The sum of vector g2 must not exceed 1.")
    }
    ids <- 1:haplo$nsamples
    names(ids) <- haplo$id
    ids <- ids[colnames(kinship)]
    X <- ghap.hslice(haplo = haplo, ids = ids, alleles = major,
                     index = TRUE, lookup = lookup, ncores = ncores)
    if(length(major) > 1){
      X <- scale(t(X))
      X <- scale(X)
      b <- sqrt(g2*h2)
      Xb <- X%*%b
    }else{
      X <- scale(X)
      b <- sqrt(g2*h2)
      Xb <- X%*%b
    }
    if(is.null(seed) == FALSE){
      set.seed(seed)
    }
  }else{
    Xb <- 0
  }
  
  # Simulate breeding value
  u <- Xb + g
  
  # Simulate repeatability
  if(nrep > 1){
    if(is.null(seed) == FALSE){
      set.seed(seed)
    }
    p <- rnorm(ncol(kinship),sd=sqrt(r2))
  }
  
  # Simulate residuals
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }
  if(nrep > 1){
    e <- rnorm(nrep*ncol(kinship), sd=sqrt(1-h2-r2))
  }else{
    e <- rnorm(ncol(kinship), sd=sqrt(1-h2))
  }
  
  # Simulate phenotypes
  if(nrep > 1){
    y <- rep(u,each=nrep) + rep(p,each=nrep) + e
  }else{
    y <- u + e
  }
  
  #Assemble output
  sim <- NULL
  sim$h2 <- h2
  sim$g2 <- g2
  sim$major <- major
  sim$major.effect <- b
  sim$u <- as.vector(u)
  names(sim$u) <- colnames(kinship)
  if(nrep > 1){
    sim$p <- as.vector(p)
    names(sim$p) <- colnames(kinship)
  }
  sim$varu <- as.numeric(var(u))
  sim$vare <- var(e)
  if(nrep > 1){
    sim$varp <- var(p)
  }
  if(nrep > 1){
    sim$data <- data.frame(y,rep(colnames(kinship),each=nrep))
  }else{
    sim$data <- data.frame(y,colnames(kinship))
  }
  colnames(sim$data) <- c("phenotype","individual")
  sim$data$individual <- as.factor(sim$data$individual)
  
  #Generate unbalanced data
  if(balanced == FALSE){
    x <- rep(FALSE,times=nrow(sim$data))
    for(i in colnames(kinship)){
      if(is.null(seed) == FALSE){
        seed <- seed + 1
        set.seed(seed)
      }
      nrep.samp <- as.integer(runif(n = 1, min = 0, max = nrep))
      samp <- which(sim$data$individual == i)
      samp <- sample(samp, size = nrep.samp, replace = F)
      x[samp] <- TRUE
    }
    sim$data <- sim$data[x,]
  }
  
  #Return output
  return(sim)
  
}