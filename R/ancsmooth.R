#Function: ghap.ancsmooth
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Smoothing for predictions of haplotype ancestry

ghap.ancsmooth<-function(
  phase,
  admix,
  ncores=1,
  verbose=TRUE
){
  
  # Log-----------------------------------------------------------------------------------------------
  if(verbose == TRUE){
    cat("\nSmoothing haplotype ancestry predictions.\n")
  }
  
  # Smoothing function--------------------------------------------------------------------------------
  test <- unique(admix$ID)
  chr <- unique(admix$CHR)
  smooth.fun <- function(j){
    out <- NULL
    for(k in 1:length(chr)){
      pop <- unique(phase$pop[which(phase$id == test[j])])
      a <- admix[which(admix$ID == test[j] & admix$CHR == chr[k]),]
      start <- a$BP1[-nrow(a)]
      end <- a$BP1[-1]-1
      hap1 <- rep(NA, times=length(start))
      hap2 <- rep(NA, times=length(start))
      for(i in 1:length(start)){
        seg <- which(a$BP2 > start[i] & a$BP1 < end[i])
        h1 <- table(a$HAP1[seg])
        h1 <- names(h1)[which(h1 == max(h1))]
        h2 <- table(a$HAP2[seg])
        h2 <- names(h2)[which(h2 == max(h2))]
        if(length(h1) > 1){
          h1 <- NA
        }
        if(length(h2) > 1){
          h2 <- NA
        }
        hap1[i] <- h1
        hap2[i] <- h2
      }
      runs1 <- rle(hap1)
      hapout1 <- matrix(data = NA, nrow = length(runs1$lengths), ncol = 8)
      hapout1 <- as.data.frame(hapout1)
      colnames(hapout1) <- c("POP","ID","HAP","CHR","BP1","BP2","SIZE","ANCESTRY")
      hapout1$BP1 <- start[c(1,cumsum(runs1$lengths[-length(runs1$lengths)])+1)]
      hapout1$BP2 <- end[cumsum(runs1$lengths)]
      hapout1$ANCESTRY <- runs1$values
      hapout1$HAP <- 1
      hapout1$ID <- test[j]
      hapout1$POP <- pop
      hapout1$CHR <- chr[k]
      if(nrow(hapout1) > 1){
        hapout1$BP1[2:nrow(hapout1)] <- hapout1$BP2[1:(nrow(hapout1)-1)] + 1 
      }
      hapout1$BP1[1] <- 1
      hapout1$BP2[nrow(hapout1)] <- max(phase$bp[which(phase$chr == chr[k])])
      hapout1$SIZE <- 1 + hapout1$BP2 - hapout1$BP1
      runs2 <- rle(hap2)
      hapout2 <- matrix(data = NA, nrow = length(runs2$lengths), ncol = 8)
      hapout2 <- as.data.frame(hapout2)
      colnames(hapout2) <- c("POP","ID","HAP","CHR","BP1","BP2","SIZE","ANCESTRY")
      hapout2$BP1 <- start[c(1,cumsum(runs2$lengths[-length(runs2$lengths)])+1)]
      hapout2$BP2 <- end[cumsum(runs2$lengths)]
      hapout2$ANCESTRY <- runs2$values
      hapout2$HAP <- 2
      hapout2$ID <- test[j]
      hapout2$POP <- pop
      hapout2$CHR <- chr[k]
      if(nrow(hapout2) > 1){
        hapout2$BP1[2:nrow(hapout2)] <- hapout2$BP2[1:(nrow(hapout2)-1)] + 1 
      }
      hapout2$BP1[1] <- 1
      hapout2$BP2[nrow(hapout2)] <- max(phase$bp[which(phase$chr == chr[k])])
      hapout2$SIZE <- 1 + hapout2$BP2 - hapout2$BP1
      hapout <- rbind(hapout1,hapout2)
      out <- c(out, as.vector(t(hapout)))
    }
    return(out)
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
  
  # Compute runs of ancestry--------------------------------------------------------------------------
  if(Sys.info()["sysname"] == "Windows"){
    results <- lapply(FUN = smooth.fun, X = 1:length(test))
  }else{
    results <- mclapply(FUN = smooth.fun, X = 1:length(test), mc.cores = ncores)
  }
  hapout <- NULL
  results <- matrix(data = unlist(results), ncol = 8, byrow = TRUE)
  results <- as.data.frame(results, stringsAsFactors = F)
  colnames(results) <-   c("POP","ID","HAP","CHR","BP1","BP2","SIZE","ANCESTRY")
  results$HAP <- as.numeric(results$HAP)
  results$BP1 <- as.numeric(results$BP1)
  results$BP2 <- as.numeric(results$BP2)
  results$SIZE <- as.numeric(results$SIZE)
  results <- results[order(results$ID,results$HAP),]
  
  #Compute admixture proportions and generations------------------------------------------------------
  test <- unique(results$ID)
  pops <- unique(results$ANCESTRY)
  pops <- pops[which(is.na(pops) == FALSE)]
  pops <- pops[order(pops)]
  propout <- matrix(data = NA, ncol = 2+length(pops)+1, nrow = length(test))
  propout <- as.data.frame(propout)
  colnames(propout) <- c("POP","ID",pops,"UNK")
  # genout <- propout[,-ncol(propout)]
  for(j in 1:length(test)){
    propout$POP[j] <- unique(phase$pop[which(phase$id == test[j])])
    propout$ID[j] <- test[j]
    # genout$POP[j] <- propout$POP[j]
    # genout$ID[j] <- propout$ID[j]
    a <- results[which(results$ID == test[j]),]
    genome <- sum(a$SIZE)
    for(i in 1:length(pops)){
      tmp <- a$SIZE[which(a$ANCESTRY == pops[i])]
      propout[j,pops[i]] <- sum(tmp)/genome
      # L <- tmp/100e+6
      # g <- round(median(1/(2*L)),digits = 0)
      # genout[j,pops[i]] <- pmax(g,1)
    }
    propout$UNK[j] <- sum(a$SIZE[which(is.na(a$ANCESTRY) == TRUE)])/genome
  }
  
  #Re-calibration of ancestry proportions------------------------------------------------------------
  propout2 <- propout[,-ncol(propout)]
  for(i in 1:nrow(propout2)){
    props <- propout2[i,-c(1:2)]
    tot <- sum(props)
    propout2[i,-c(1:2)] <- props/tot
  }
  
  #Output--------------------------------------------------------------------------------------------
  if(verbose == TRUE){
    cat("\nProcessed a total of",j,"individuals.\n")
    cat("Analysis done.\n\n")
  }
  out <- NULL
  out$proportions1 <- propout
  out$proportions2 <- propout2
  #out$generations <- genout
  out$haplotypes <- results
  return(out)
  
}