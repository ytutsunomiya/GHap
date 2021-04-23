#Function: ghap.ancplot
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Plot ancestry proportions

ghap.ancplot <- function(
  ancsmooth,
  labels=TRUE,
  pop.ang=45, 
  group.ang=0,
  colors=NULL,
  pop.order=NULL,
  sortby=NULL,
  use.unk=TRUE,
  legend=TRUE
){
  
  #Load admixture results
  if(use.unk == TRUE){
    admix <- ancsmooth$proportions1
  }else{
    admix <- ancsmooth$proportions2
  }
  admix[,-c(1:2)] <- admix[,-c(1:2)]*100
  npop <- ncol(admix) - 2
  
  #Check population order
  if(is.null(pop.order) == TRUE){
    admix$POP <- factor(x = admix$POP)
  }else{
    admix$POP <- factor(x = admix$POP, levels = unlist(pop.order))
  }
  
  #Check sorting of individuals
  if(is.null(sortby) == TRUE){
    admix <- admix[order(admix$POP,admix$ID),]
  }else{
    admix <- admix[order(admix$POP,admix[,sortby]),]
  }
  
  #Set colors
  if(is.null(colors) == TRUE){
    colors <- c("#80B1D3","#FB8072","#8DD3C7","#BEBADA","#FFFFB3","#FDB462","#B3DE69")
    colors <- colors[1:npop]
    if(use.unk == TRUE){
      colors[npop] <- "grey"
    }
  }else if(length(colors) != npop){
    stop("Number of colors must match number of groups")
  }
  
  #Capture user's parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  #Plot
  if(legend == TRUE){
    par(mar=c(5, 4, 4, 2) + 0.1)
    par(xpd=T, mar=par()$mar+c(2,0,0,0))
  }
  p <- barplot(t(admix[,-c(1:2)]), space = 0, border = NA, col = colors, las=1, ylab = "Ancestry (%)",
               xaxt="n", ylim=c(0,125), yaxt = "n")
  axis(side = 2, at = seq(0,100,by=20), labels = seq(0,100,by=20), las =1)
  admix$POS <- p
  p <- aggregate(POS ~ POP, data = admix, FUN = median)
  xmin <- aggregate(POS ~ POP, data = admix, FUN = min)
  xmax <- aggregate(POS ~ POP, data = admix, FUN = max)
  u <- par("usr")
  if(labels == TRUE){
    text(x=p$POS, y=u[3]-0.1*(u[4]-u[3]),
         labels=p$POP, srt=pop.ang, adj=1, xpd=TRUE)
  }
  rect(xleft = min(admix$POS)-0.5, ybottom = 0,
       xright = max(admix$POS)+0.5, ytop = 100,
       col="#00000000", border=T)
  segments(x0 = xmax$POS+0.5, y0 = 0, x1 = xmax$POS+0.5, y1 = 100)
  if(is.list(pop.order) == TRUE){
    labcol <- rep(c("black","darkgrey"), times = length(pop.order))
    labcol <- labcol[1:length(pop.order)]
    for(i in 1:length(pop.order)){
      segments(x0 = min(xmin$POS[which(xmin$POP %in% pop.order[[i]])])-0.5, y0 = 102.5,
               x1 = max(xmax$POS[which(xmax$POP %in% pop.order[[i]])])+0.5, y1 = 102.5,
               col = labcol[i], lwd = 3)
      # arrows(x0 = min(xmin$POS[which(xmin$POP %in% pop.order[[i]])])-0.5, y0 = 103,
      #          x1 = max(xmax$POS[which(xmax$POP %in% pop.order[[i]])])+0.5, y1 = 103,
      #          col = "black", lwd = 2, angle = 90, code =3, length = 0.05)
      text(x = median(p$POS[which(p$POP %in% pop.order[[i]])]), y = 108,
           labels = names(pop.order)[i], srt=group.ang, col = "black")
    }
  }
  
  #Legend
  if(legend == TRUE){
    if(use.unk == TRUE){
      ltext <- colnames(admix[-c(1,2,ncol(admix),ncol(admix)-1)])
      ltext <- c(ltext, "?")
    }else{
      ltext <- colnames(admix[-c(1,2,ncol(admix),ncol(admix))])
    }
    legend(x =  median(p$POS), y = -40, legend = ltext,
           fill = colors, bty="n", ncol=npop, xjust=0.6)
    par(mar=c(5, 4, 4, 2) + 0.1)
  }
  
  
}