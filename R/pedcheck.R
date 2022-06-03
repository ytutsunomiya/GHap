#Function: ghap.pedcheck
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Pedigree check

ghap.pedcheck <- function(ped, depth.n.f = FALSE){
  
  # Sanity check for pedigree --------------------------------------------------
  if(inherits(ped, "data.frame") == FALSE){
    stop("A data frame must be passed to argument 'ped'")
  }
  if(identical(c("dam","id","sire"),sort(colnames(ped))) == FALSE){
    stop("Column names in the pedigree data frame must be 'id', 'sire' and 'dam'")
  }
  
  # Compute statistics ---------------------------------------------------------
  ped <- ped[,c("id","sire","dam")]
  ped[ped == 0] <- NA
  ids <- unique(c(ped$id,ped$sire,ped$dam))
  ids <- ids[which(is.na(ids) == FALSE)]
  sires <- unique(ped$sire)
  sires <- sires[which(is.na(sires) == FALSE)]
  dams <- unique(ped$dam)
  dams <- dams[which(is.na(dams) == FALSE)]
  nosire <- length(which(is.na(ped$sire) == TRUE & is.na(ped$dam) == FALSE))
  nodam <- length(which(is.na(ped$sire) == FALSE & is.na(ped$dam) == TRUE))
  noparents <- length(which(is.na(ped$sire) == TRUE & is.na(ped$dam) == TRUE))
  siredam <- table(c(sires,dams))
  siredam <- length(which(siredam > 1))
  self <- length(which(ped$sire == ped$dam))
  
  # Summary table --------------------------------------------------------------
  statistic = c("total records","unique individuals",
                "individuals listed as sires",
                "individuals listed as dams",
                "individuals listed as both sire and dam",
                "records with sire = dam (self-fertilization)",
                "records with missing sire",
                "records with missing dam",
                "records with both parents missing")
  values <- c(nrow(ped),length(ids),length(sires),length(dams),
              siredam,self,nosire,nodam,noparents)
  stats <- data.frame(statistic, values, stringsAsFactors = FALSE)
  
  # Inbreeding and generation depth --------------------------------------------
  if(depth.n.f == TRUE){
    
    # Inbreeding calculation
    tmp <- editPed(sire = ped$sire, dam = ped$dam, label = ped$id)
    tmp2 <- pedigree(sire = tmp$sire, dam = tmp$dam, label = tmp$label)
    f <- inbreeding(ped = tmp2)
    identical(tmp$label, tmp2@label)
    tmp$f <- f
    ped <- tmp
    colnames(ped) <- c("id","sire","dam","gen","f")
    
    # Generation depth per individual
    ped$ecg <- 0
    for(i in 1:nrow(ped)){
      myid <- ped[i,]
      parents <- c(myid$sire,myid$dam)
      parents <- parents[which(is.na(parents) == F)]
      newparents <- parents
      parentsg <- 0
      j <- 0
      while(length(newparents) > 0){
        j <- j + 1
        idx <- which(ped$id %in% newparents)
        sires <- ped$sire[idx]
        dams <- ped$dam[idx]
        sires <- sires[which(is.na(sires) == F)]
        dams <- dams[which(is.na(dams) == F)]
        parents <- c(parents, sires, dams)
        parentsg <- c(parentsg, rep(j, times = length(sires)+length(dams)))
        newparents <- c(sires,dams)
      }
      ped$ecg[i] <- sum((1/2)^parentsg)
    }
    
  }
  
  # Make output ----------------------------------------------------------------
  out <- NULL
  out$stats <- stats
  out$ped <- ped
  return(out)
  
  
}
