#Function: ghap.relfind
#License: GPLv3 or later
#Modification date: 12 Mar 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Find relatives in IBD estimates

ghap.relfind <- function(
  ibdpairs,
  v = 50,
  break3 = FALSE,
  ncores=1
){
  
  # Check if input is a valid IBD list -----------------------------------------
  ghapcols <- c("POP1","ID1","POP2","ID2","Z0","Z1","Z2","PI_HAT")
  plinkcols <- c("FID1","IID1","FID2","IID2","Z0","Z1","Z2","PI_HAT")
  ghapcols <- length(which(ghapcols %in% colnames(ibdpairs)))
  plinkcols <- length(which(plinkcols %in% colnames(ibdpairs)))
  if(ghapcols != 8 & plinkcols != 8){
    stop("\nInvalid or missing columns in your ibdpairs data.")
  }
  
  #Initialize main function ----------------------------------------------------
  ibdfun <- function(i){
    
    # Test pihat and Z values
    Z0 = ibdpairs$Z0[i]
    Z1 = ibdpairs$Z1[i]
    Z2 = ibdpairs$Z2[i]
    pihat <- ibdpairs$PI_HAT[i]
    s <- rep(NA, times = 7)
    
    # 1 - Duplicates or monozygotic twins (-1)
    m <- c(0.01,0.01,0.01,0.01)
    a <- m*v
    b <- (1-m)*v
    p <- pbeta(q = c(1-pihat,Z0,Z1,1-Z2), shape1 = a, shape2 = b, lower.tail = F, log.p = T)
    s[1] <- sum(p)
    print(s[1])
    
    # 2 - Parent-offspring with inbreeding (0)
    m <- c(0.75,0.01,0.50,0.50)
    a <- m*v
    b <- (1-m)*v
    p <- pbeta(q = c(pihat,Z0,Z1,Z2), shape1 = a, shape2 = b, lower.tail = F, log.p = T)
    s[2] <- sum(p)
    print(s[2])
    
    # 3 - Parent-offspring (1)
    m <- c(0.50,0.01,0.01,0.01)
    a <- m*v
    b <- (1-m)*v
    p <- pbeta(q = c(pihat,Z0,1-Z1,Z2), shape1 = a, shape2 = b, lower.tail = F, log.p = T)
    s[3] <- sum(p)
    print(s[3])
    
    # 4 - Full-siblings (2)
    m <- c(0.50,0.25,0.50,0.25)
    a <- m*v
    b <- (1-m)*v
    p <- pbeta(q = c(pihat,Z0,Z1,Z2), shape1 = a, shape2 = b, lower.tail = F, log.p = T)
    s[4] <- sum(p)
    print(s[4])
    
    # 5 - Half-siblings, Grandparent-grandchild or avuncular with inbreeding (3.1)
    m <- c(0.375,0.375,0.50,0.125)
    a <- m*v
    b <- (1-m)*v
    p <- pbeta(q = c(pihat,Z0,Z1,Z2), shape1 = a, shape2 = b, lower.tail = F, log.p = T)
    s[5] <- sum(p)
    print(s[5])
    
    # 6 - Half-siblings, Grandparent-grandchild or avuncular (3.2)
    m <- c(0.25,0.50,0.50,0.01)
    a <- m*v
    b <- (1-m)*v
    p <- pbeta(q = c(pihat,Z0,Z1,Z2), shape1 = a, shape2 = b, lower.tail = F, log.p = T)
    s[6] <- sum(p)
    print(s[6])
    
    # 7 - Cousin or Half-avuncular (3.3)
    m <- c(0.125,0.75,0.25,0.01)
    a <- m*v
    b <- (1-m)*v
    p <- pbeta(q = c(pihat,Z0,Z1,Z2), shape1 = a, shape2 = b, lower.tail = F, log.p = T)
    s[7] <- sum(p)
    print(s[7])
    
    # 8 - Half-cousin (3.4)
    m <- c(0.0625,0.875,0.125,0.01)
    a <- m*v
    b <- (1-m)*v
    p <- pbeta(q = c(pihat,Z0,Z1,Z2), shape1 = a, shape2 = b, lower.tail = F, log.p = T)
    s[8] <- sum(p)
    print(s[8])
    
    # 9 - Unrelated (4)
    m <- c(0.01,0.01,0.01,0.01)
    a <- m*v
    b <- (1-m)*v
    p <- pbeta(q = c(pihat,1-Z0,Z1,Z2), shape1 = a, shape2 = b, lower.tail = F, log.p = T)
    s[9] <- sum(p)
    print(s[9])
    
    # Class
    if(break3 == FALSE){
      scores <- c(-1,0,1,2,3,3,3,3,4)
      s <- scores[which(s == max(s))]
    }else{
      scores <- c(-1,0,1,2,3.1,3.2,3.3,3.4,4)
      s <- scores[which(s == max(s))]
    }
    return(s)
    
  }
  
  #iteration function ----------------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cl <- makeCluster(ncores)
    clusterExport(cl = cl, varlist = c("ibdpairs","inbreeding"), envir=environment())
    results <- unlist(parLapply(cl = cl, fun = ibdfun, X = 1:nrow(ibdpairs)))
    stopCluster(cl)
  }else{
    results <- mclapply(FUN = ibdfun, X = 1:nrow(ibdpairs), mc.cores = ncores)
  }
  
  #Return output ---------------------------------------------------------------
  ibdpairs$REL <- unlist(results)
  return(ibdpairs)
  
}
