#Function: ghap.phase2plink
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Convert phase object into plink file

ghap.phase2plink <- function(
  object,
  out.file,
  only.active.samples=TRUE,
  only.active.markers=TRUE,
  batchsize=NULL,
  ncores=1,
  verbose=TRUE
){
  
  # Check if phase is a GHap.phase object --------------------------------------
  if(inherits(object, "GHap.phase") == FALSE){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Insert suffix to out.file --------------------------------------------------
  fam.file <- paste(out.file, "fam", sep=".")
  bim.file <- paste(out.file, "bim", sep=".")
  bed.file <- paste(out.file, "bed", sep=".")
  
  # Check if output will overwrite existing files ------------------------------
  if(file.exists(fam.file) == TRUE){
    stop(paste("File", fam.file, "already exists!"))
  }
  if(file.exists(bed.file) == TRUE){
    stop(paste("File", bed.file, "already exists!"))
  }
  if(file.exists(bim.file) == TRUE){
    stop(paste("File", bim.file, "already exists!"))
  }
  
  # Check if inactive markers and samples should be reactived ------------------
  if(only.active.markers == FALSE){
    object$marker.in <- rep(TRUE,times=object$nmarkers)
    object$nmarkers.in <- length(which(object$marker.in))
  }
  if(only.active.samples == FALSE){
    object$id.in <- rep(TRUE,times=2*object$nsamples)
    object$nsamples.in <- length(which(object$id.in))/2
  }
  nmarkers <- object$nmarkers.in
  
  # Identify activated samples -------------------------------------------------
  ids.in <- which(object$id.in)
  id <- object$id[ids.in]
  id <- id[1:length(id) %% 2 == 0]
  pop <- object$pop[ids.in]
  pop <- pop[1:length(pop) %% 2 == 0]
  ids.n <- length(id)
  sire <- object$sire
  dam <- object$dam
  sex <- object$sex
  sire[which(is.na(sire) == TRUE)] <- "0"
  dam[which(is.na(dam) == TRUE)] <- "0"
  sex[which(is.na(sex) == TRUE)] <- "0"
  sex[which(sex == "M")] <- "1"
  sex[which(sex == "F")] <- "2"
  
  # Output fam.file file -------------------------------------------------------
  fwrite(x = as.data.table(cbind(pop,id,sire,dam,sex,"-9")), file = fam.file,
         quote = FALSE, sep=" ", row.names = FALSE, col.names=FALSE)
  
  
  # Generate batch index -------------------------------------------------------
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(nmarkers/10)
  }
  if(batchsize > nmarkers){
    batchsize <- nmarkers
  }
  id1 <- seq(1,nmarkers,by=batchsize)
  id2 <- (id1+batchsize)-1
  id1 <- id1[id2<=nmarkers]
  id2 <- id2[id2<=nmarkers]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,nmarkers)
  if(id1[length(id1)] > nmarkers){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  # Log message ----------------------------------------------------------------
  if(verbose == TRUE){
    cat("Processing ", nmarkers, " markers in:\n", sep="")
    batch <- table((id2-id1)+1)
    for(i in 1:length(batch)){
      cat(batch[i]," batches of ",names(batch[i]),"\n",sep="")
    }
  }
  
  # Initialize lookup table for output -----------------------------------------
  lookup <- rep(NA,times=256)
  lookup[1:2] <- c(0,1)
  d <- 10
  i <- 3
  while(i <= 256){
    b <- d + lookup[1:(i-1)]
    lookup[i:(length(b)+i-1)] <- b
    i <- i + length(b)
    d <- d*10
  }
  lookup <- sprintf(fmt="%08d", lookup)
  lookup <- sapply(lookup, function(i){intToUtf8(rev(utf8ToInt(i)))})
  
  # Compute bitloss ------------------------------------------------------------
  bitloss <- 8 - ((2*ids.n) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  
  #Define binary conversion function
  toBitFUN <- function(k){
    line <- Ztmp[k,]
    line[which(line == "0")] <- "00"
    line[which(line == "1")] <- "01"
    line[which(line == "2")] <- "11"
    line <- c(line,rep("0",times=bitloss))
    line <- paste(line, collapse = "")
    nc <- nchar(line)
    n <- seq(1, nc, by = 8)
    line <- substring(line, n, c(n[-1]-1, nc))
    line <- strtoi(lookup[line], base=2)
    return(line)
  }
  
  # Print magic number [01101100 00011011] and mode [00000001] -----------------
  bed.file.con  <- file(bed.file, open = "ab")
  writeBin(object = as.integer(c(108,27,1)), con = bed.file.con, size = 1)
  close(con = bed.file.con)
  
  # Iterate batch --------------------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  sumvariants <- 0
  var.in <- which(object$marker.in)
  for(i in 1:length(id1)){
    idx <- id1[i]:id2[i]
    Ztmp <- ghap.slice(object = object,
                       ids = ids.in,
                       variants = var.in[idx],
                       index = TRUE,
                       unphase = TRUE,
                       impute = TRUE,
                       ncores = ncores)
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      bed.file.out <- unlist(parLapply(cl = cl, fun = toBitFUN, X = 1:nrow(Ztmp)))
      stopCluster(cl)
    }else{
      bed.file.out <- unlist(mclapply(FUN = toBitFUN, X = 1:nrow(Ztmp),
                                      mc.cores = ncores))
    }
    bed.file.con <- file(bed.file, open = "ab")
    writeBin(object = bed.file.out, con = bed.file.con, size = 1)
    close(con = bed.file.con)
    if(verbose == TRUE){
      sumvariants <- sumvariants + length(idx)
      cat(sumvariants, "variants processed.\r")
    }
  }
  
  # Export marker map ----------------------------------------------------------
  chr <- object$chr[var.in]
  marker <- object$marker[var.in]
  cm <- object$cm[var.in]
  bp <- object$bp[var.in]
  a1 <- object$A0[var.in]
  a2 <- object$A1[var.in]
  fwrite(x = as.data.table(cbind(chr,marker,cm,bp,a1,a2)), file = bim.file,
         quote = FALSE, sep=" ", row.names = FALSE, col.names=FALSE)
  
}
