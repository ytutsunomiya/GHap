#Function: ghap.anc2plink
#License: GPLv3 or later
#Modification date: 3 Jun 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Output ancestry genotype matrix

ghap.anc2plink <- function(
  object,
  ancsmooth,
  ancestry,
  outfile,
  freq=c(0,1),
  missingness=1,
  only.active.samples=TRUE,
  only.active.markers=TRUE,
  batchsize=NULL,
  binary=TRUE,
  ncores=1,
  verbose=TRUE
){
  
  # Check if object is a GHap.phase object
  if(inherits(object, "GHap.phase") == FALSE){
    stop("Argument object must be a GHap.phase object.")
  }
  
  # Insert suffix to outfile
  bed <- paste(outfile,"bed",sep=".")
  bim <- paste(outfile,"bim",sep=".")
  fam <- paste(outfile,"fam",sep=".")
  if(binary == FALSE){
    bed <- paste(outfile,"txt",sep=".")
  }
  
  # Check if output will overwrite existing files before opening connection
  if(file.exists(fam) == TRUE){
    stop(paste("File", fam, "already exists!"))
  }
  if(file.exists(bed) == TRUE){
    stop(paste("File", bed, "already exists!"))
  }
  if(file.exists(bim) == TRUE){
    stop(paste("File", bim, "already exists!"))
  }
  
  # Check if inactive markers should be reactived
  if(only.active.markers == FALSE){
    object$marker.in <- rep(TRUE,times=object$nmarkers)
    object$nmarkers.in <- length(which(object$marker.in))
  }
  if(only.active.samples == FALSE){
    object$id.in <- rep(TRUE,times=2*object$nsamples)
    object$nsamples.in <- length(which(object$id.in))/2
  }
  
  # Identify activated samples
  ids.in <- which(object$id.in)
  id <- object$id[ids.in]
  id <- id[1:length(id) %% 2 == 0]
  pop <- object$pop[ids.in]
  pop <- pop[1:length(pop) %% 2 == 0]
  ids.n <- length(id)
  anchaplotypes <- ancsmooth$haplotypes[which(ancsmooth$haplotypes$ID %in% id),]
  
  # Identify activated markers
  snp <- which(object$marker.in)
  nsnp <- length(snp)
  
  # Check if ancestry exists
  ancs <- unique(anchaplotypes$ANCESTRY)
  ancs <- ancs[which(is.na(ancs) == FALSE)]
  if(ancestry %in% ancs == FALSE | length(ancestry) != 1){
    emsg <- "\n\nSelected ancestry should be one in:"
    emsg <- paste(emsg, paste(sort(ancs), collapse = ", "), sep="\n")
    stop(emsg)
  }
  
  # Output fam file
  fwrite(x = as.data.table(cbind(pop,id,"0 0 0 -9")), file = fam, quote = FALSE,
         sep=" ", row.names = FALSE, col.names=FALSE)
  
  # Generate batch index
  if(is.null(batchsize) == TRUE){
    batchsize <- 1000
  }
  if(batchsize > nsnp){
    batchsize <- nsnp
  }
  id1 <- seq(1,nsnp,by=batchsize)
  id2 <- (id1 + batchsize) - 1
  id1 <- id1[id2 <= nsnp]
  id2 <- id2[id2 <= nsnp]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,nsnp)
  if(id1[length(id1)] > nsnp){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  # Log message
  if(verbose == TRUE){
    cat("Processing ", nsnp, " markers in:\n", sep="")
    batch <- table((id2-id1)+1)
    for(i in 1:length(batch)){
      cat(batch[i]," batches of ",names(batch[i]),"\n",sep="")
    }
  }
  
  # Initialize lookup table for output
  lookup2 <- rep(NA,times=256)
  lookup2[1:2] <- c(0,1)
  d <- 10
  i <- 3
  while(i <= 256){
    b <- d + lookup2[1:(i-1)]
    lookup2[i:(length(b)+i-1)] <- b
    i <- i + length(b)
    d <- d*10
  }
  lookup2 <- sprintf(fmt="%08d", lookup2)
  lookup2 <- sapply(lookup2, function(i){intToUtf8(rev(utf8ToInt(i)))})
  
  # Print magic number [01101100 00011011] and mode [00000001]
  # For compatibility with PLINK
  if(binary == TRUE){
    bed.con  <- file(bed, open = "ab")
    writeBin(object = as.integer(c(108,27,1)), con = bed.con, size = 1)
    close(con = bed.con)
  }
  
  # Define marker function
  marker.iter.FUN<-function(j){
    
    # Get marker info
    chr <- object$chr[j]
    bp <- object$bp[j]
    
    # Map overlapping ancestry tracks
    tmp <- which(anchaplotypes$CHR == chr & anchaplotypes$BP1 <= bp & anchaplotypes$BP2 >= bp)
    tmp <- anchaplotypes[tmp,]
    
    # Build ancestry counts
    X <- table(tmp$ID,tmp$ANCESTRY, useNA = "always")
    X <- X[-nrow(X),]
    miss <- which(X[,ncol(X)] != 0)
    x <- X[,ancestry]
    x[miss] <- NA
    
    # Order x vector
    x <- x[id]
    
    # Calculate frequency
    x2 <- x[which(is.na(x) == FALSE)]
    p <- sum(x2)/length(x2)
    
    # Calculate number of missing values
    nmiss <- length(miss)/length(x)
    
    # Output line
    x <- paste(x, collapse = " ")
    return(c(x,p,nmiss))
    
  }
  
  #Compute bitloss
  bitloss <- 8 - ((2*ids.n) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  
  #Define binary conversion function
  toBitFUN <- function(j){
    line <- scan(text=bed.out[j], what = "character", sep=" ", quiet = TRUE)
    line[which(line == "0")] <- "00"
    line[which(line == "1")] <- "01"
    line[which(line == "2")] <- "11"
    line[which(line == "NA")] <- "10"
    line <- c(line,rep("0",times=bitloss))
    line <- paste(line, collapse = "")
    nc <- nchar(line)
    n <- seq(1, nc, by = 8)
    line <- substring(line, n, c(n[-1]-1, nc))
    line <- strtoi(lookup2[line], base=2)
    return(line)
  }
  
  #Iterate markers
  nmarkers.done <- 0
  for(i in 1:length(id1)){
    
    #Get marker indices
    idx <- snp[id1[i]:id2[i]]
    
    #Compute markers
    ncores <- min(c(detectCores(), ncores))
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      mylines <- unlist(parLapply(cl = cl, fun = marker.iter.FUN, X =idx))
      stopCluster(cl)
    }else{
      mylines <- unlist(mclapply(FUN = marker.iter.FUN, X = idx, mc.cores = ncores))
    }
    mylines <- matrix(data = mylines, ncol = 3, byrow = T)
    bed.out <- mylines[,1]
    hapfreq <- as.numeric(mylines[,2])
    hapmiss <- as.numeric(mylines[,3])
    snpok <- which(hapfreq >= freq[1] & hapfreq <= freq[2] & hapmiss <= missingness)
    
    # Prepare output
    if(length(snpok) > 0){
      
      # Filter markers
      idx <- idx[snpok]
      bed.out <- bed.out[snpok]
      
      #Write bim
      bim.out <- paste(object$chr[idx], object$marker[idx], "0", object$bp[idx], "N H", sep=" ")
      bim.con  <- file(bim, open = "a") 
      writeLines(text = bim.out, con=bim.con)
      close(con = bim.con)
      
      #Write bed
      if(binary == TRUE){
        if(Sys.info()["sysname"] == "Windows"){
          bed.out <- unlist(lapply(FUN = toBitFUN, X = 1:length(bed.out)))
        }else{
          bed.out <- unlist(mclapply(FUN = toBitFUN, X = 1:length(bed.out), mc.cores = ncores))
        }
        bed.con  <- file(bed, open = "ab")
        writeBin(object = bed.out, con = bed.con, size = 1)
        close(con = bed.con)
      }else{
        bed.con  <- file(bed, open = "a")
        writeLines(text = bed.out, con=bed.con)
        close(con = bed.con)
      }
      
    }
    
    # Log message
    if(verbose == TRUE){
      nmarkers.done <- nmarkers.done + (id2[i]-id1[i]) + 1
      cat(nmarkers.done, "markers processed\r")
    }
    
  }
  
}
