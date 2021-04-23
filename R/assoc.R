#Function: ghap.assoc
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: fit ordinary least squares for each haplotype allele

ghap.assoc<-function(
  response,
  haplo,
  weights=NULL,
  gc=TRUE,
  only.active.alleles=TRUE,
  batchsize=NULL,
  ncores=1,
  verbose=TRUE
){
  
  #General data check
  if (class(haplo) != "GHap.haplo") {
    stop("Argument haplo must be a GHap.haplo object.")
  }
  if (only.active.alleles == FALSE) {
    haplo$allele.in <- rep(TRUE, times = haplo$nalleles)
    haplo$nalleles.in <- length(which(haplo$allele.in))
  }
  if (length(which(names(response) %in% haplo$id)) != length(response)) {
    stop("All ids in the response must be present in the GHap.haplo object.")
  }
  ids <- unique(names(response))
  
  #Check if response contains missing values
  if(length(na.omit(response)) != length(response)){
    stop("Missing values are not accepted.")
  }
  
  #Include weights
  if(is.null(weights)==FALSE){
    w <- sqrt(weights/mean(weights))
    response <- w*response
  }
  
  # Get number of cores
  if(Sys.info()["sysname"] == "Windows"){
    if(ncores > 1 & verbose == TRUE){
      cat("\nParallelization not supported yet under Windows (using a single core).")
    }
    ncores <- 1
  }else{
    ncores <- min(c(detectCores(), ncores))
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
  
  #Generate batch index
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(haplo$nalleles.in/10)
  }
  if(batchsize > haplo$nalleles.in){
    batchsize <- haplo$nalleles.in
  }
  activealleles <- which(haplo$allele.in)
  nbatches <- round(haplo$nalleles.in/(batchsize),digits=0) + 1
  mybatch <- paste("B",1:nbatches,sep="")
  batch <- rep(mybatch,each=batchsize)
  batch <- batch[1:haplo$nalleles.in]
  mybatch <- unique(batch)
  nbatches <- length(mybatch)
  
  #Log message
  if(verbose == TRUE){
    cat("Processing ", haplo$nalleles.in, " HapAlleles in ", nbatches, " batches.\n", sep="")
    cat("Inactive alleles will be ignored.\n")
  }
  
  #ols iterate function
  ols.FUN <- function(j){
    x <- hap.geno[j,names(response)]
    frq <- sum(x)/(2*length(x))
    x <- x - mean(x)
    if(is.null(weights)==F){
      x <- w*x
    }
    X <- cbind(1,x)
    Xt <- t(X)
    XtXi <- solve(t(X)%*%X)
    Xty <- Xt%*%response
    b <- XtXi%*%Xty
    vare <- sum((response - as.numeric(X%*%b))^2)/(length(response)-ncol(X))
    se <- sqrt(vare*diag(XtXi))
    b <- as.numeric(b)[2]
    se <- as.numeric(se)[2]
    return(c(b,se,frq))
  }
  
  #Iterate batches
  hapreg <- NULL
  hapreg$BLOCK <- haplo$block[haplo$allele.in]
  hapreg$CHR <- haplo$chr[haplo$allele.in]
  hapreg$BP1 <- haplo$bp1[haplo$allele.in]
  hapreg$BP2 <- haplo$bp2[haplo$allele.in]
  hapreg$ALLELE <- haplo$allele[haplo$allele.in]
  hapreg$BETA <- rep(NA, times=haplo$nalleles.in)
  hapreg$SE <- rep(NA, times=haplo$nalleles.in)
  hapreg$FREQ <- rep(NA, times=haplo$nalleles.in)
  sumalleles <- 0
  for(i in 1:nbatches){
    idx <- which(batch == mybatch[i])
    slice <- activealleles[idx]
    hap.geno <- ghap.hslice(haplo = haplo, ids = ids, alleles = slice,
                     index = FALSE, lookup = lookup, ncores = ncores)
    if(Sys.info()["sysname"] == "Windows"){
      a <- lapply(X = 1:nrow(hap.geno), FUN = ols.FUN)
    }else{
      a <- mclapply(FUN=ols.FUN, X=1:nrow(hap.geno), mc.cores = ncores)
    }
    a <- data.frame(matrix(unlist(a), nrow=nrow(hap.geno), byrow=TRUE))
    hapreg$BETA[idx] <- a[,1]
    hapreg$SE[idx] <- a[,2]
    hapreg$FREQ[idx] <- a[,3]
    if(verbose == TRUE){
      sumalleles <- sumalleles + length(idx)
      cat(sumalleles, "HapAlleles processed.\r")
    }
  }
  hapreg$CHISQ.OBS <- (hapreg$BETA^2)/(hapreg$SE^2)
  hapreg$CHISQ.EXP <- qchisq(p = rank(hapreg$CHISQ.OBS)/(haplo$nalleles.in+1), df=1)
  if(gc == TRUE){
    dev <- sd(hapreg$CHISQ.OBS)
    samp <- which(hapreg$CHISQ.OBS < 3*dev)
    lambda <- lm(hapreg$CHISQ.OBS[samp] ~ hapreg$CHISQ.EXP[samp])
    lambda <- lambda$coefficients[2]
    hapreg$CHISQ.OBS <- hapreg$CHISQ.OBS/lambda
  }
  hapreg$logP <- -1*pchisq(q = hapreg$CHISQ.OBS, df = 1, lower.tail=FALSE, log.p = TRUE)/log(10)
  hapreg <- data.frame(hapreg,stringsAsFactors = FALSE)

  #Return object
  return(hapreg)
  
}
