#Function: ghap.roh
#License: GPLv3 or later
#Modification date: 27 Apr 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Map streches of homozygous genotypes

ghap.roh <- function(
  phase,
  minroh=1e+6,
  method="hmm",
  freq=NULL,
  inbcoef=NULL,
  error=0.25/100,
  only.active.samples = TRUE,
  only.active.markers = TRUE,
  ncores = 1,
  verbose = TRUE
){
  
  # Check if phase is a GHap.phase object-------------------------------------------------------------
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Check if inactive markers and samples should be reactivated---------------------------------------
  if(only.active.markers == FALSE){
    phase$marker.in <- rep(TRUE,times=phase$nmarkers)
    phase$nmarkers.in <- length(which(phase$marker.in))
  }
  if(only.active.samples == FALSE){
    phase$id.in <- rep(TRUE,times=2*phase$nsamples)
    phase$nsamples.in <- length(which(phase$id.in))/2
  }
  
    # Get number of cores-------------------------------------------------------------------------------
  if(Sys.info()["sysname"] == "Windows"){
    if(ncores > 1 & verbose == TRUE){
      cat("\nParallelization not supported yet under Windows (using a single core).")
    }
    ncores <- 1
  }else{
    ncores <- min(c(detectCores(), ncores))
  }
  
  # Initialize lookup table----------------------------------------------------------------------------
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
  
  # ROH function---------------------------------------------------------------------------------------
  if(method == "hmm"){
    rohfun <- function(i){
      
      #Get vector of observations
      if(is.vector(geno)){
        x <- geno
      }else{
        x <- geno[,ids[i]]
      }
      x[which(x == 2)] <- 0
      x <- x + 1
      
      #Get inbreeding coefficient
      f <- inbcoef[ids[i]]
      
      #Starting state probabilities
      states <- c("ROH","N")
      start <- c(f,1-f)
      names(start) <- states
      
      #Emission probabilities
      emiss.roh <- c(1-error,error)
      emiss.n <- list(freqchr^2 + (1-freqchr)^2, 2*freqchr*(1-freqchr))

      #Transition probabilities
      expr <- exp(-2*(c(bps[1],diff(bps))/1e+6)/100)
      trans.roh2roh <- expr + (1-expr)*f
      trans.roh2n <- (1-expr)*(1-f)
      trans.n2n <- expr + (1-expr)*(1-f)
      trans.n2roh <- (1-expr)*f

      #Build states vector
      v <- array(NA, c(2, m))
      dimnames(v) = list(states = states, marker = mkrs)
      v[1,1] <- log(start[1]*emiss.roh[x[1]])
      v[2,1] <- log(start[2]*emiss.n[[x[1]]][1])

      #Get likelihood of states
      for (k in 2:m){
        maxi <- max(v[1,k-1] + log(trans.roh2roh[k]), v[2,k-1] + log(trans.n2roh[k]))
        v[1,k] <- log(emiss.roh[x[k]]) + maxi
        maxi <- max(v[1,k-1] + log(trans.roh2n[k]), v[2,k-1] + log(trans.n2n[k]))
        v[2,k] <- log(emiss.n[[x[k]]][k]) + maxi
      }

      #Guess Viterbi path
      viterbiPath <- rep(NA, m)
      viterbiPath[m] <- states[which(v[,m] == max(v[,m]))]
      for(k in (m - 1):1){
        if(viterbiPath[k+1] == "ROH"){
          L <- c(v[1,k] + log(trans.roh2roh[k]), v[2,k] + log(trans.n2roh[k]))
        }else{
          L <- c(v[1,k] + log(trans.roh2n[k]), v[2,k] + log(trans.n2n[k]))
        }
        viterbiPath[k] <- states[which(L == max(L))]
      }

      #Get runs
      runs <- rle(viterbiPath)
      runsum <- cumsum(runs$lengths)
      idx1 <- c(1,runsum[-length(runsum)]+1)
      idx2 <- c(runsum)
      bp1 <- bps[idx1]
      bp2 <- bps[idx2]
      runvals <- runs$values
      runs <- bp2-bp1+1
      keep <- which(runs > minroh & runvals == "ROH")
      bp1 <- bp1[keep]
      bp2 <- bp2[keep]
      runs <- runs[keep]
      pop <- unique(phase$pop[which(phase$id == ids[i])])
      if(length(runs) == 0){
        out <- NULL
      }else{
        out <- as.vector(rbind(pop,ids[i],chr,bp1,bp2,runs))
      }
      return(out)
      
    }
      
  }else if(method == "naive"){
    rohfun <- function(i){
      if(is.vector(geno)){
        x <- geno
      }else{
        x <- geno[,ids[i]]
      }
      x[which(x == 2)] <- 0
      runs <- rle(x)
      runsum <- cumsum(runs$lengths)
      idx1 <- c(1,runsum[-length(runsum)]+1)
      idx2 <- c(runsum)
      bp1 <- bps[idx1]
      bp2 <- bps[idx2]
      runvals <- runs$values
      runs <- bp2-bp1+1
      keep <- which(runs > minroh & runvals == "0")
      bp1 <- bp1[keep]
      bp2 <- bp2[keep]
      runs <- runs[keep]
      pop <- unique(phase$pop[which(phase$id == ids[i])])
      if(length(runs) == 0){
        out <- NULL
      }else{
        out <- as.vector(rbind(pop,ids[i],chr,bp1,bp2,runs))
      }
      return(out)
    }
  }else{
    stop('The method argument has to be either "hmm" or "naive"')
  }

  
  # Find runs of homozygosity--------------------------------------------------------------------------
  ids <- unique(phase$id[which(phase$id.in)])
  outruns <- NULL
  if(method == "hmm" & is.null(freq) == TRUE){
    stop('\nMethod "hmm" requires reference allele frequencies.\n')
  }
  if(method == "hmm" & is.null(freq) == TRUE){
    emsg <- '\nMethod "hmm" requires starting values for genomic inbreeding.\n'
    emsg <- paste0(emsg, "(i.e., proportion of the genome covered by ROH)")
    stop(emsg)
  }
  if(verbose == TRUE){
    cat('\n\nSearching for runs of homozygosity using the "', method, '" method.\n', sep="")
    cat("Number of individuals to search:", phase$nsamples.in, "\n")
    cat("Number of markers to search:", phase$nmarkers.in,"\n\n")
  }
  for(chr in unique(phase$chr)){
    if(verbose == TRUE){
      cat("Finding runs of homozygosity on chromosome", chr, "\r")
    }
    mkrs <- phase$marker[which(phase$chr == chr & phase$marker.in == TRUE)]
    m <- length(mkrs)
    bps <- phase$bp
    names(bps) <- phase$marker
    bps <- bps[mkrs]
    freqchr <- freq[mkrs]
    geno <- ghap.pslice(phase = phase, ids = ids, ncores = ncores,
                        markers = mkrs, unphase = TRUE, lookup = lookup)
    if(Sys.info()["sysname"] == "Windows"){
      #cl <- makeCluster(ncores)
      #segs <- parLapply(cl = cl, fun = rohfun, X = 1:length(ids))
      #stopCluster(cl)
      segs <- lapply(FUN = rohfun, X = 1:length(ids))
    }else{
      segs <- mclapply(FUN = rohfun, X = 1:length(ids), mc.cores = ncores)
    }
    outruns <- c(outruns,unlist(segs))
  }

  # Compile results------------------------------------------------------------------------------------
  if(verbose == TRUE){
    cat("\nCompiling results... ")
  }
  results <- matrix(data = outruns, ncol = 6, byrow = TRUE)
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  colnames(results) <- c("POP","ID","CHR","BP1","BP2","LENGTH")
  results$BP1 <- as.integer(results$BP1)
  results$BP2 <- as.integer(results$BP2)
  results$LENGTH <- as.integer(results$LENGTH)
  results <- results[order(results$POP, results$ID, results$CHR, results$BP1),]
  if(verbose == TRUE){
    cat("Done.\n\n")
  }
  return(results)

  
}
