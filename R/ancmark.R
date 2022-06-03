#Function: ghap.ancmark
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Per marker ancestry

ghap.ancmark <- function(object, ancsmooth, ids){
  
  # Check if object is of class GHap.phase -------------------------------------
  if(inherits(object, "GHap.phase") == FALSE){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Get populations ------------------------------------------------------------
  pops <- colnames(ancsmooth$proportions1)[-c(1:2)]
  npop <- length(pops)
  
  # Get chromosomes ------------------------------------------------------------
  chr <- unique(ancsmooth$haplotypes$CHR)
  
  # Get segments ---------------------------------------------------------------
  results <- NULL
  for(i in chr){
    tmp <- which(ancsmooth$haplotypes$CHR == i & ancsmooth$haplotypes$ID %in% ids)
    tmp <- ancsmooth$haplotypes[tmp,]
    bp1 <- min(tmp$BP1)
    bp2 <- min(tmp$BP2)
    maxbp <- max(tmp$BP2)+1
    while(bp2 < maxbp){
      segs <- which(tmp$BP1 <= bp1 & tmp$BP2 >= bp2)
      tbl <- table(tmp$ANCESTRY[segs], useNA = "always")
      names(tbl)[which(is.na(names(tbl)))] <- "UNK"
      tbl <- 100*tbl/sum(tbl)
      out <- rep(0, times = length(pops))
      names(out) <- pops
      out[names(tbl)] <- tbl
      out <- c(i, bp1, bp2, out)
      results <- c(results,out)
      bp1 <- bp2 + 1
      nxt1 <- tmp$BP1[which(tmp$BP1 > bp1)]
      nxt2 <- tmp$BP2[which(tmp$BP2 > bp1)]
      if(length(nxt1) != 0 | length(nxt2) != 0){
        bp2 <- min(c(nxt1,nxt2)) 
      }else{
        bp2 <- maxbp
      }
    }
  }
  
  # Organize results -----------------------------------------------------------
  a <- matrix(data = unlist(results), ncol = 3+length(pops), byrow = TRUE)
  results <- matrix(data = NA, nrow = nrow(a), ncol = ncol(a))
  results <- as.data.frame(results)
  colnames(results) <- c("CHR","BP1","BP2",pops)
  results[,1] <- a[,1]
  results[,-1] <- as.numeric(a[,-1])
  
  # Organize map ---------------------------------------------------------------
  map <- matrix(data = NA, ncol = 3+length(pops), nrow = object$nmarkers, byrow = TRUE)
  map <- as.data.frame(map)
  colnames(map) <- c("CHR","MARKER","BP", pops)
  map$CHR <- object$chr
  map$MARKER <- object$marker
  map$BP <- object$bp
  for(i in 1:nrow(results)){
    mkrs <- which(map$CHR == results$CHR[i] & map$BP >= results$BP1[i] & map$BP <= results$BP2[i])
    map[mkrs,-c(1:3)] <- results[i,-c(1:3)]
  }
  
  # #Set colors
  # if(is.null(colors) == TRUE){
  #   colors <- brewer.pal(7,"Set3")
  #   colors <- colors[c(5,4,1,3,2,6,7)]
  #   colors[npop] <- "grey"
  # }else if(length(colors) != npop){
  #   stop("Number of colors must match number of groups")
  # }
  
  # #Plot chromosomes
  # #layout(matrix(data = 1:10, nrow = 10, ncol = 1, byrow = TRUE), heights = rep(0.5, times=5))
  # par(mfrow = c(29,1), mai = c(0.02,0.82,0.02,0.42))
  # for(i in chr){
  #   tmp <- map[which(map$CHR == i),]
  #   plot(x = 1, y = 1, xlim = c(1,max(map$BP/1e+6)), ylim=c(0,100), pch="", las=1,
  #        xlab = "", bty = "n", xaxt = "n", yaxt = "n", ylab = "")
  #   # ticks <- seq(0,ceiling(max(tmp$BP/1e+6)), by = 20)
  #   # axis(side = 1, at = ticks, labels = ticks)
  #   summary <- rep(0, nrow(tmp))
  #   recent <- summary
  #   for(p in 1:length(pops)) {
  #     current <- tmp[,pops[p]]
  #     summary <- summary + current
  #     polygon(
  #       x=c(tmp$BP/1e+6, rev(tmp$BP/1e+6)),
  #       y=c(summary, rev(recent)),
  #       col=colors[p]
  #     )
  #     recent <- summary
  #   }
  #   
  # }
  
  return(map)
  
}
