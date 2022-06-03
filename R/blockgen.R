#Function: ghap.blockgen
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Generate blocks based on sliding windows

ghap.blockgen<-function(
  object,
  windowsize=10,
  slide=5,
  unit="marker",
  nsnp=2
){
  
  #Check if object is of class GHap.phase
  if(inherits(object, "GHap.phase") == FALSE){
    stop("Argument phase must be a GHap.phase object.")
  }
  if(unit %in% c("marker","kbp","ibd") == FALSE){
    stop("Unit must be specified as 'marker' or 'kbp'")
  }
  
  #Initialize vectors
  BLOCK <- rep(NA,times=object$nmarkers.in)
  CHR <- rep(NA,times=object$nmarkers.in)
  BP1 <- rep(NA,times=object$nmarkers.in)
  BP2 <- rep(NA,times=object$nmarkers.in)
  SIZE <- rep(NA,times=object$nmarkers.in)
  NSNP <- rep(NA,times=object$nmarkers.in)
  
  if(unit == "kbp"){
    
    windowsize <- windowsize*1e+3
    slide <- slide*1e+3
    offset <- 0
    for(k in unique(object$chr)){
      cmkr <- which(object$marker.in & object$chr == k)
      nmkr <- length(cmkr)
      bp <- object$bp[cmkr]
      minbp <- bp[1]
      maxbp <- bp[nmkr]
      id1<-seq(1,maxbp,by=slide)
      id2<-id1+(windowsize-1)
      id1<-id1[id2 <= maxbp]
      id2<-id2[id2 <= maxbp]
      for(i in 1:length(id1)){
        slice <- cmkr[which(bp >= id1[i] & bp <= id2[i])]
        CHR[i+offset] <- object$chr[slice[1]]
        BP1[i+offset] <- id1[i]
        BP2[i+offset] <- id2[i]
        NSNP[i+offset] <- length(slice)
        BLOCK[i+offset] <- paste("CHR",CHR[i+offset],"_B",i,sep="")
      }
      offset <- offset + length(id1)
    }
    
  }else if(unit == "marker"){
    
    offset <- 0
    for(k in unique(object$chr)){
      cmkr <- which(object$marker.in & object$chr == k)
      nmkr <- length(cmkr)
      id1<-seq(1,nmkr,by=slide)
      id2<-id1+(windowsize-1)
      id1<-id1[id2 <= nmkr]
      id2<-id2[id2 <= nmkr]
      for(i in 1:length(id1)){
        slice <- cmkr[id1[i]:id2[i]]
        CHR[i+offset] <- object$chr[slice[1]]
        BP1[i+offset] <- object$bp[slice[1]]
        BP2[i+offset] <- object$bp[slice[length(slice)]]
        NSNP[i+offset] <- length(slice)
        BLOCK[i+offset] <- paste("CHR",CHR[i+offset],"_B",i,sep="")
      }
      offset <- offset + length(id1)
    }
    
    
  }
  
  #Organize results
  results <- data.frame(BLOCK,CHR,BP1,BP2,NSNP,stringsAsFactors = FALSE)
  results <- unique(results)
  results$SIZE <- 1 + results$BP2 - results$BP1
  results$SIZE[results$NSNP == 1] <- 1
  results <- results[order(nchar(results$CHR),results$CHR,results$BP1,results$BP2),]
  results <- results[results$NSNP >= nsnp,]
  if(nrow(results) == 0){
    stop("No blocks could be generated with the specified options. Try setting the nsnp argument to a smaller value.")
  }
  results <- na.exclude(results)
  results <- results[,c("BLOCK","CHR","BP1","BP2","SIZE","NSNP")]
  return(results)
  
}
