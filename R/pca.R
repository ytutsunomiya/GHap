#Function: ghap.pca
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute principal components from a kinship matrix

ghap.pca<-function(haplo, kinship, npc=2){
  
  #Check if haplo is a GHap.haplo object
  if (class(haplo) != "GHap.haplo") {
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if kinship matrix is symmetrical
  if(identical(colnames(kinship),rownames(kinship)) == FALSE){
    stop("Names in rows and columns must be identical.")
  }
  
  #Check if names in the kinship matrix match with the GHap.haplo object
  if (length(which(colnames(kinship) %in% haplo$id)) != ncol(kinship)) {
    stop("All ids in the kinship matrix must be present in the GHap.haplo object.")
  }else{
    ids <- rep(NA, times = ncol(kinship))
    for (i in 1:length(ids)) {
      ids[i] <- which(haplo$id == colnames(kinship)[i])
    }
    pop <- haplo$pop[ids]
  }
  
  #Eigendecomposition
  eigK <- eigen(kinship, symmetric = TRUE)
  
  #Output
  results <- NULL
  results$eigenvec <- data.frame(pop,rownames(kinship),eigK$vectors[,1:npc],stringsAsFactors = FALSE)
  colnames(results$eigenvec) <- c("POP","ID",paste("PC",1:npc,sep=""))
  results$eigenval <- eigK$values[1:npc]
  results$propvar <- (eigK$values/sum(eigK$values))[1:npc]
  
  #Return output
  return(results)
  
}