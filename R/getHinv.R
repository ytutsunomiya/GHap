#Function: ghap.getHinv
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute the inverse of H

ghap.getHinv <- function(
  K,
  ped,
  include = NULL,
  depth = 3,
  alpha = 0.95,
  verbose=TRUE
){
  
  # Sanity check for pedigree ------------------------------------------------------
  if(inherits(ped, "data.frame") == FALSE){
    stop("A data frame must be passed to argument 'ped'")
  }
  if(identical(c("dam","id","sire"),sort(colnames(ped))) == FALSE){
    stop("Column names in the pedigree data frame must be 'id', 'sire' and 'dam'")
  }
  ped <- ped[,c("id","sire","dam")]
  ped[ped == 0] <- NA
  if(verbose == TRUE){
    ids <- unique(c(ped$id,ped$sire,ped$dam))
    ids <- ids[which(is.na(ids) == FALSE)]
    sires <- unique(ped$sire)
    sires <- sires[which(is.na(sires) == FALSE)]
    dams <- unique(ped$dam)
    dams <- dams[which(is.na(dams) == FALSE)]
    pad <- paste0("%",nchar(as.character(length(ids)))+2,"d")
    nosire <- length(which(is.na(ped$sire) == TRUE & is.na(ped$dam) == FALSE))
    nodam <- length(which(is.na(ped$sire) == FALSE & is.na(ped$dam) == TRUE))
    noparents <- length(which(is.na(ped$sire) == TRUE & is.na(ped$dam) == TRUE))
    siredam <- table(c(sires,dams))
    siredam <- length(which(siredam > 1))
    self <- length(which(ped$sire == ped$dam))
    out <- paste0("Raw pedigree summary:\n",
                  sprintf(pad,nrow(ped)), " total records\n",
                  sprintf(pad,length(ids)), " unique individuals in the pedigree\n",
                  sprintf(pad,length(sires)), " individuals listed as sires\n",
                  sprintf(pad,length(dams)), " individuals listed as dams\n",
                  sprintf(pad,siredam), " individuals listed as both sire and dam\n",
                  sprintf(pad,self), " records with sire = dam (self-fertilization)\n",
                  sprintf(pad,nosire), " records with missing sire\n",
                  sprintf(pad,nodam), " records with missing dam\n",
                  sprintf(pad,noparents), " records with both parents missing\n")
    cat(out)
  }
  
  # Sanity check for relationship matrix -------------------------------------------
  if(is.null(colnames(K)) | is.null(rownames(K))){
    stop("The genomic relationship matrix must contain row and column names")
  }
  if(identical(colnames(K),rownames(K)) == FALSE){
    stop("Row and column names in the genomic relationship matrix must be identical")
  }
  if(verbose == TRUE){
    g <- which(colnames(K) %in% ids)
    p <- which(ids %in% colnames(K) == FALSE)
    pad <- paste0("%",nchar(nrow(K))+2,"d")
    out <- paste0("\nGenotype summary:\n",
                  sprintf(pad,nrow(K)), " total genotyped individuals\n",
                  sprintf(pad,length(g)), " genotyped individuals included in the pedigree (",
                  sprintf("%.1f",100*length(g)/nrow(K)), "%)\n",
                  sprintf(pad,length(p)), " individuals from pedigree with no genotype (",
                  sprintf("%.1f",100*length(p)/length(ids)), "%)\n")
    cat(out)
  }
  
  # Crop pedigree with the desired generation depth --------------------------------
  ids <- unique(c(colnames(K),include))
  i <- 1
  while(i <= depth){
    tmp <- ped[which(ped$id %in% ids),]
    sire <- unique(tmp$sire[which(is.na(tmp$sire) == FALSE)])
    dam <- unique(tmp$dam[which(is.na(tmp$dam) == FALSE)])
    if(length(sire) == 0 & length(dam) == 0){
      i <- depth
    }else{
      ids <- unique(c(ids,sire,dam))
      i <- i + 1
    }
  }
  ped <- ped[which(ped$id %in% ids),]
  ids <- ids[which(ids %in% ped$id == FALSE)]
  ped <- rbind(ped, data.frame(id = ids, sire = NA, dam = NA, stringsAsFactors = FALSE))
  if(verbose == TRUE){
    ids <- unique(c(ped$id,ped$sire,ped$dam))
    ids <- ids[which(is.na(ids) == FALSE)]
    sires <- unique(ped$sire)
    sires <- sires[which(is.na(sires) == FALSE)]
    dams <- unique(ped$dam)
    dams <- dams[which(is.na(dams) == FALSE)]
    pad <- paste0("%",nchar(as.character(length(ids)))+2,"d")
    nosire <- length(which(is.na(ped$sire) == TRUE & is.na(ped$dam) == FALSE))
    nodam <- length(which(is.na(ped$sire) == FALSE & is.na(ped$dam) == TRUE))
    noparents <- length(which(is.na(ped$sire) == TRUE & is.na(ped$dam) == TRUE))
    siredam <- table(c(sires,dams))
    siredam <- length(which(siredam > 1))
    self <- length(which(ped$sire == ped$dam))
    out <- paste0("\nFiltered pedigree summary (", depth, " generations):\n",
                  sprintf(pad,length(ids)), " unique individuals in the pedigree\n",
                  sprintf(pad,length(sires)), " individuals listed as sires\n",
                  sprintf(pad,length(dams)), " individuals listed as dams\n",
                  sprintf(pad,siredam), " individuals listed as both sire and dam\n",
                  sprintf(pad,self), " records with sire = dam (self-fertilization)\n",
                  sprintf(pad,nosire), " records with missing sire\n",
                  sprintf(pad,nodam), " records with missing dam\n",
                  sprintf(pad,noparents), " records with both parents missing\n")
    cat(out)
  }
  
  # Compute numerator relationship matrix ------------------------------------------
  if(verbose == TRUE){
    cat("\nComputing numerator relationship matrix (A)... ")
  }
  ped <- editPed(sire = ped$sire, dam = ped$dam, label = ped$id)
  ped <- pedigree(sire = ped$sire, dam = ped$dam, label = ped$label)
  A <- drop0(x = getA(ped), tol = 1e-5)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  # Re-organize matrices -----------------------------------------------------------
  g <- colnames(K)
  n <- colnames(A)[which(colnames(A) %in% g == F)]
  A <- A[c(n,g),c(n,g)]
  K <- K[g,g]
  Agg <- A[g,g]
  
  # Check matrix regression --------------------------------------------------------
  if(verbose == TRUE){
    r <- cor(as.numeric(K), as.numeric(Agg))
    out <- paste0("Correlation between K and A: ", sprintf("%.3f", r), ".\n")
    cat(out)
  }
  
  # Compute inverse of relationship matrices ---------------------------------------
  if(verbose == TRUE){
    out <- paste0("Computing inverse of A... ")
    cat(out)
  }
  Ainv <- drop0(getAInv(ped), tol = 1e-5)
  Ainv <- Ainv[c(n,g),c(n,g)]
  if(verbose == TRUE){
    cat("Done.\n")
    out <- paste0("Computing inverse of ", alpha, "*K + ", 1-alpha, "*A... ")
    cat(out)
  }
  Kinv <- solve(alpha*K + (1-alpha)*Agg)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  # Blend inverses -----------------------------------------------------------------
  if(verbose == TRUE){
    cat("Computing inverse of blend matrix (H)... ")
  }
  #Hinv <- Ainv
  #Hinv[g,g] <- Hinv[g,g] + (Kinv - solve(Agg))
  #Hinv[g,g] <- Kinv + Ainv[g,n]%*%solve(Ainv[n,n])%*%Ainv[n,g]
  Hinv <- rbind(cbind(Ainv[n,n], Ainv[n,g]),
                cbind(Ainv[g,n], Kinv + Ainv[g,n]%*%solve(Ainv[n,n])%*%Ainv[n,g]))
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  # Return H inverse ---------------------------------------------------------------
  if(verbose == TRUE){
    Hinv <- drop0(Hinv)
    nz <- nnzero(Hinv)
    de <- prod(dim(Hinv))
    sparsity <- 100*(de - nz)/de
    out <- paste0("Final H inverse has a sparsity of ", sprintf("%.1f", sparsity), "%.\n",
                  length(g), " genotyped and ", length(n), " non-genotyped individuals.\n")
    cat(out)
  }
  return(Hinv)
  
}
