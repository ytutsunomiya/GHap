#Function: ghap.manhattan
#License: GPLv3 or later
#Modification date: 15 May 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Generate a Manhattan plot

ghap.manhattan <- function(
  data,
  chr,
  bp, 
  y,
  colors=NULL,
  type="p",
  pch=20,
  cex=1,
  lwd=1,
  ylim=NULL,
  ylab="",
  xlab="",
  main="",
  backcolor="#F5EFE780",
  chr.ang=0,
  hlines = NULL,
  hcolors = NULL,
  hlty = 1,
  hlwd = 1
){
  
  # Check integrity of chromosome column--------------------------------------------------------------
  if(chr %in% colnames(data) == FALSE){
    emsg <- paste("There is no column '",chr,"' in your dataframe",sep="")
    stop(emsg)
  }
  if(is.factor(data[,chr]) == FALSE){
    chrs <- unique(unique(data[,chr]))
    nchr <- length(chrs)
    chrorder <- chrs[order(nchar(chrs),chrs)]
    data[,chr] <- factor(x = data[,chr], levels = chrorder, labels = chrorder)
  }
  
  
  # Check integrity of bp column----------------------------------------------------------------------
  if(bp %in% colnames(data) == FALSE){
    emsg <- paste("There is no column '",bp,"' in your dataframe",sep="")
    stop(emsg)
  }
  if(is.numeric(data[,bp]) == FALSE){
    emsg <- paste("Column '",bp,"' must be a numeric vector",sep="")
    stop(emsg)
  }
  if(length(which(data[,bp] < 0)) != 0){
    emsg <- paste("Negative values in the '",bp,"' column",sep="")
    stop(emsg)
  }
  
  # Check integrity of y column-----------------------------------------------------------------------
  if(y %in% colnames(data) == FALSE){
    emsg <- paste("There is no column '",y,"' in your dataframe",sep="")
    stop(emsg)
  }
  if(is.numeric(data[,y]) == FALSE){
    emsg <- paste("Column '",y,"' must be a numeric vector",sep="")
    stop(emsg)
  }
  
  # Order dataframe-----------------------------------------------------------------------------------
  data <- data[order(data[,chr], data[,bp]),]
  
  # Make relative coordinates-------------------------------------------------------------------------
  chrom <- levels(data[,chr])
  offset <- 0
  data$relpos <- NA
  ticks <- vector(mode = "numeric", length = length(chrom))
  for(i in 1:length(chrom)){
    snps <- which(data[,chr] == chrom[i])
    data$relpos[snps] <- data[snps,bp] + offset
    offset <- offset + max(data[snps,bp])
    ticks[i] <- mean(data$relpos[snps])
  }
  
  # Set colors----------------------------------------------------------------------------------------
  if(is.null(colors) == TRUE){
    colors <- hsv(h = c(0.6,0.35,1,0,0.12), s = c(0.7,0.5,0.7,0,1), v = c(0.8,0.7,0.9,0,0.95))
  }
  colors <- rep(x = colors, length.out = nlevels(data[,chr]))
  
  # Background plot-----------------------------------------------------------------------------------
  plot(x = data$relpos, y = data[,y], las=1, xaxt="n", pch = "",
       ylim = ylim, ylab = ylab, xlab = xlab, main = main)
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col=backcolor, border=T)
  text(x = ticks, y = u[3], labels = chrom, srt = chr.ang, pos = 1, xpd = TRUE)
  
  # Vertical lines------------------------------------------------------------------------------------
  if(is.null(hlines) == FALSE){
    
    #Check line types
    if(is.null(hlty) == TRUE){
      hlty <- rep(1, times=length(hlines))
    }else if(length(hlty) < length(hlines)){
      hlty <- rep(hlty, length.out = length(hlines))
    }
    
    #Check line widths
    if(is.null(hlwd) == TRUE){
      hlwd <- rep(1, times=length(hlines))
    }else if(length(hlwd) < length(hlines)){
      hlwd <- rep(hlwd, length.out = length(hlines))
    }
    
    #Plot lines
    for(i in 1:length(hlines)){
      abline(h = hlines[i], lty = hlty[i], lwd = hlwd[i], col = hcolors[i])
    }
    
  }
  
  # Plot markers--------------------------------------------------------------------------------------
  if(type == "l"){
    for(i in 1:length(chrom)){
      snps <- which(data[,chr] == chrom[i])
      points(x = data$relpos[snps], y = data[snps,y], col = colors[i],
             type = type, pch = pch, cex = cex, lwd=lwd)
    }
  }else{
    points(x = data$relpos, y = data[,y], col = colors[as.numeric(data[,chr])],
           type = type, pch = pch, cex = cex)
  }
  
  
  
}
