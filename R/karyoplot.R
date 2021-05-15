#Function: ghap.karyoplot
#License: GPLv3 or later
#Modification date: 15 May 2021
#Written by: Marco Milanesi, Yuri Tani Utsunomiya
#Contact: marco.milanesi.mm@gmail.com, ytutsunomiya@gmail.com
#Description: Individual chromosome painting

ghap.karyoplot <- function(
  ancsmooth, 
  ids=NULL, 
  colors=NULL,
  chr=NULL,
  chr.line=10, 
  plot.line=25, 
  chr.ang=45, 
  las=0){
  
  #Load admixture results
  admix <- ancsmooth$haplotypes
  
  #Check if chromosomes are factors
  if(is.factor(admix$CHR) == FALSE){
    chr <- unique(admix$CHR)
    nchr <- length(chr)
    chrorder <- chr[order(nchar(chr),chr)]
    admix$CHR <- factor(x = admix$CHR, levels = chrorder, labels = chrorder)
  }
  
  #Find number of Ks
  pop <- sort(unique(admix$ANCESTRY))
  pop <- c(pop,"?")
  npop <- length(pop)
  admix$ANCESTRY[is.na(admix$ANCESTRY)] <- "?"
  admix$BP1 <- admix$BP1/1e+6
  admix$BP2 <- admix$BP2/1e+6
  admix$SIZE <- admix$SIZE/1e+6
  
  #Extract the chromosome(s) to analyse
  if (!is.null(chr)){
    admix <- admix[which(admix$CHR %in% chr), ]
  }else{
    chr <- sort(unique(admix$CHR))
  }
  if (nrow(admix)==0){
    stop("The selected chromosome(s) is(are) not existing. Please verify the input file. ")
  }
  
  #Find number of chromosomes
  nchr <- length(unique(admix$CHR))
  if (nchr<chr.line){chr.line=nchr}
  
  #Find the max Mbp
  maxMb <- max(tapply(X = admix$BP2, INDEX = admix$CHR, FUN = max, na.rm=T))
  
  #Set colors
  if(is.null(colors) == TRUE){
    colors <- c("#80B1D3","#FB8072","#8DD3C7","#BEBADA","#FFFFB3","#FDB462","#B3DE69")
    colors <- colors[1:npop]
    colors[npop] <- "grey"
  }else if(length(colors) != npop){
    stop("Number of colors must match number of groups")
  }
  colors <- as.data.frame(colors, stringsAsFactors = F)
  colors$pop <- pop
  
  #Select ids
  if (is.null(ids) == TRUE){
    ids <- unique(admix$ID)
  }
  
  #Capture user's parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  #Plot
  for (id in ids){
    #Subset data
    id.admix <- admix[which(admix$ID == id), ]
    
    #Parameters
    par(mar=c(1, 4, 1, 2) + 0.1) 
    xlim <- chr.line*3
    chr.start <- 1
    chr.end <- chr.line
    if (chr.end > nchr){chr.end <- nchr}
    
    #Identify the number of plot line
    nline <- ceiling(nchr/chr.line)
    
    #Plot layout
    if (nline>=3){
      x.lay <- round(10/(nline+0.6),2)
      w.lay <- c(x.lay*1.3,rep(x.lay,nline-2),x.lay*1.3)
      layout(matrix(1:nline,ncol = 1), heights = w.lay)
      par(mar=c(5, 4, 4, 2) + 0.1)
    }else{
      layout(matrix(1:nline,ncol = 1))
      par(mar=c(5, 4, 4, 2) + 0.1)
    }
    
    #Plot
    for (line in 1:nline){
      
      #Margins
      par(mar=c(1, 4, 1, 2) + 0.1) 
      if (line == 1){
        par(mar=c(1, 4, 4, 2) + 0.1)
      }
      if (line == nline){
        par(mar=c(4, 4, 1, 2) + 0.1)
      }
      if (line == nline & line == 1){
        par(mar=c(4, 4, 4, 2) + 0.1)
      }
      
      #Plot
      if (line == 1){
        plot(0,0, ylim = c(0,maxMb), xlim = c(0,xlim), type = "n", 
             xlab = "", ylab = "Mbp", xaxt="n", bty="n", las=las, 
             main = paste(id.admix$POP[1]," - ",id.admix$ID[1],sep=""))
      }else{
        plot(0,0, ylim = c(0,maxMb), xlim = c(0,xlim), type = "n", 
             xlab = "", ylab = "Mbp", xaxt="n", bty="n", las=las)
      }
      
      #Abline
      if (is.numeric(plot.line)){
        abline(h = seq(0,(maxMb+plot.line),plot.line), col="grey", lwd=0.5, lty=2)
      }
      
      #Select chromosomes
      chr.list <- chr[c(chr.start:chr.end)]
      for (k in chr.list){
        k.admix <- id.admix[which(id.admix$CHR == k), ]
        pos <- (which(chr.list == k)-1)*3
        #Select haplotype 
        for (hap in c(1,2)){
          hap.admix <- k.admix[which(k.admix$HAP == hap), ]
          hap.admix <- merge(x = hap.admix, y = colors, by.x = "ANCESTRY", by.y = "pop", all.x = T,sort = F) 
          hap.admix <- hap.admix[order(hap.admix$BP1), ]
          #Plot Admixture
          for (seg in 1:nrow(hap.admix)){
            rect(xleft = pos+hap, ybottom = hap.admix$BP1[seg],
                 xright = pos+hap+0.70, ytop = hap.admix$BP2[seg],
                 col=hap.admix$colors[seg], border=NA)
          }
        }
      }
      
      #Labels
      u <- par("usr")
      text(x=((0:(length(chr.list)-1))*3)+1.90, y=u[3]-0.03*(u[4]-u[3]),
           labels=chr.list, srt=chr.ang, adj=1, xpd=TRUE)
      
      #Legend
      if (line == nline){
        par(xpd=T)
        ltext <- colors$pop
        legend(x =  median(((0:(chr.line-1))*3)+1.85), y = u[3]-(0.03*(u[4]-u[3]))*2.5, 
               legend = ltext, fill = colors$colors, 
               bty="n", ncol=npop, xjust=0.6)
        par(xpd=F)
      }
      
      #Change chr.start and chr.end
      if (line != nline){
        chr.start <- chr.start + chr.line
        chr.end <- chr.end + chr.line
        if (chr.end > nchr){chr.end <- nchr}
      }
    }
    par(mar=c(1, 4, 1, 2) + 0.1) 
  }
}
