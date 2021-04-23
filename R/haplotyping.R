#Function: ghap.haplotyping
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Output haplotype genotype matrix for user-defined haplotype blocks

ghap.haplotyping <- function(
  phase,
  blocks,
  outfile,
  freq=c(0,1),
  drop.minor=FALSE,
  only.active.samples=TRUE,
  only.active.markers=TRUE,
  batchsize=NULL,
  binary=TRUE,
  ncores=1,
  verbose=TRUE
){
  
  #Check if phase is a GHap.phase object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  #Insert suffix to outfile
  hapgenotypes <- paste(outfile,"hapgenotypes",sep=".")
  hapalleles <- paste(outfile,"hapalleles",sep=".")
  hapsamples <- paste(outfile,"hapsamples",sep=".")
  if(binary == TRUE){
    hapgenotypes <- paste(hapgenotypes,"b",sep="")
  }
  
  #Check if output will overwrite existing files before opening connection
  if(file.exists(hapsamples) == TRUE){
    stop(paste("File", hapsamples, "already exists!"))
  }
  if(file.exists(hapgenotypes) == TRUE){
    stop(paste("File", hapgenotypes, "already exists!"))
  }
  if(file.exists(hapalleles) == TRUE){
    stop(paste("File", hapalleles, "already exists!"))
  }
  
  #Check if inactive markers and samples should be reactived
  if(only.active.markers == FALSE){
    phase$marker.in <- rep(TRUE,times=phase$nmarkers)
    phase$nmarkers.in <- length(which(phase$marker.in))
  }
  if(only.active.samples == FALSE){
    phase$id.in <- rep(TRUE,times=2*phase$nsamples)
    phase$nsamples.in <- length(which(phase$id.in))/2
  }
  
  #Identify activated samples
  ids.in <- which(phase$id.in)
  id <- phase$id[ids.in]
  id <- id[1:length(id) %% 2 == 0]
  pop <- phase$pop[ids.in]
  pop <- pop[1:length(pop) %% 2 == 0]
  ids.n <- length(id)
  
  #Output hapsamples file
  fwrite(x = as.data.table(cbind(pop,id)), file = hapsamples, quote = FALSE,
         sep=" ", row.names = FALSE, col.names=FALSE)
  
  #Generate batch index
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(nrow(blocks)/10)
  }
  if(batchsize > nrow(blocks)){
    batchsize <- nrow(blocks)
  }
  id1<-seq(1,nrow(blocks),by=batchsize)
  id2<-(id1+batchsize)-1
  id1<-id1[id2<=nrow(blocks)]
  id2<-id2[id2<=nrow(blocks)]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,nrow(blocks))
  if(id1[length(id1)] > nrow(blocks)){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  #Log message
  if(verbose == TRUE){
    cat("Processing ", nrow(blocks), " blocks in:\n", sep="")
    batch <- table((id2-id1)+1)
    for(i in 1:length(batch)){
      cat(batch[i]," batches of ",names(batch[i]),"\n",sep="")
    }
  }
  
  #Windows warning
  ncores <- min(c(detectCores(), ncores))
  if(Sys.info()["sysname"] == "Windows" & ncores > 1 & verbose == TRUE){
    cat("\nParallelization not supported yet under Windows (using a single core).\n")
  }
  
  # Initialize lookup table for input
  lookup1 <- rep(NA,times=256)
  lookup1[1:2] <- c(0,1)
  d <- 10
  i <- 3
  while(i <= 256){
    b <- d + lookup1[1:(i-1)]
    lookup1[i:(length(b)+i-1)] <- b
    i <- i + length(b)
    d <- d*10
  }
  lookup1 <- sprintf(fmt="%08d", lookup1)
  
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
    hapgenotypes.con  <- file(hapgenotypes, open = "ab")
    writeBin(object = as.integer(c(108,27,1)), con = hapgenotypes.con, size = 1)
    close(con = hapgenotypes.con)
  }
  
  # Define block function
  block.iter.FUN<-function(i){
    
    outline <- NULL
    
    #Get block info
    block.info <- blocks[i,c("BLOCK","CHR","BP1","BP2")]
    
    #SNPs in the block
    snps <- which(phase$chr == block.info$CHR &
                    phase$bp >= block.info$BP1 &
                    phase$bp <= block.info$BP2 &
                    phase$marker.in == TRUE)
    phase.A0 <- phase$A0[snps]
    phase.A1 <- phase$A1[snps]
    
    if(length(snps) >= 1){
      
      #Subset block
      if(length(snps) == 1){
        block.subset <- ghap.pslice(phase = phase,
                                    ids = unique(phase$id[ids.in]),
                                    markers = phase$marker[snps],
                                    lookup = lookup1)
        haplotypes <- as.character(block.subset)
      }else{
        block.subset <- ghap.pslice(phase = phase,
                                    ids = unique(phase$id[ids.in]),
                                    markers = phase$marker[snps],
                                    lookup = lookup1)
        haplotypes <- apply(block.subset,MARGIN = 2, paste, collapse="")
      }
      
      #Haplotype library
      lib <- table(haplotypes)/(2*ids.n)
      lib <- sort(lib)
      if(drop.minor == TRUE & length(lib) > 1){
        lib <- lib[-1]
      }
      lib <- lib[lib >= freq[1] & lib <= freq[2]]
      
      #Output genotype matrix
      if(length(lib) >= 1){
        for(j in names(lib)){
          phase.geno <- as.integer(haplotypes[1:length(haplotypes) %% 2 == 0] == j)
          phase.geno <- phase.geno + as.integer(haplotypes[1:length(haplotypes) %% 2 == 1] == j)
          j <- unlist(strsplit(j,""))
          phase.allele <- rep(NA,times=length(j))
          phase.allele[j == "0"] <- phase.A0[j == "0"]
          phase.allele[j == "1"] <- phase.A1[j == "1"]
          phase.allele <- paste(phase.allele,collapse="")
          phase.info<-paste(paste(block.info,collapse=" "),phase.allele,sep=" ")
          phase.geno<-paste(phase.geno,collapse=" ")
          outline <- c(outline,phase.info,phase.geno)
        }
      }
      
    }
    return(outline)
  }
  
  #Compute bitloss
  bitloss <- 8 - ((2*ids.n) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  
  #Define binary conversion function
  toBitFUN <- function(i){
    line <- scan(text=hapgenotypes.out[i], what = "character", sep=" ", quiet = TRUE)
    line[which(line == "0")] <- "00"
    line[which(line == "1")] <- "01"
    line[which(line == "2")] <- "11"
    line <- c(line,rep("0",times=bitloss))
    line <- paste(line, collapse = "")
    nc <- nchar(line)
    n <- seq(1, nc, by = 8)
    line <- substring(line, n, c(n[-1]-1, nc))
    line <- strtoi(lookup2[line], base=2)
    return(line)
  }
  
  #Iterate blocks
  nblocks.done <- 0
  for(i in 1:length(id1)){
    
    #Compute blocks
    if(Sys.info()["sysname"] == "Windows"){
      mylines <- unlist(lapply(FUN = block.iter.FUN, X = id1[i]:id2[i]))
    }else{
      mylines <- unlist(mclapply(FUN = block.iter.FUN, X = id1[i]:id2[i], mc.cores = ncores))
    }
    
    #Write batch to files
    if(is.null(mylines) == F){
      
      #Write hapalleles
      hapalleles.out <- mylines[1:length(mylines) %% 2 == 1]
      hapalleles.con  <- file(hapalleles, open = "a") 
      writeLines(text = hapalleles.out, con=hapalleles.con)
      close(con = hapalleles.con)
      
      #Write hapgenotypes
      hapgenotypes.out <- mylines[1:length(mylines) %% 2 == 0]
      if(binary == TRUE){
        if(Sys.info()["sysname"] == "Windows"){
          hapgenotypes.out <- unlist(lapply(FUN = toBitFUN, X = 1:length(hapgenotypes.out)))
        }else{
          hapgenotypes.out <- unlist(mclapply(FUN = toBitFUN, X = 1:length(hapgenotypes.out), mc.cores = ncores))
        }
        hapgenotypes.con  <- file(hapgenotypes, open = "ab")
        writeBin(object = hapgenotypes.out, con = hapgenotypes.con, size = 1)
        close(con = hapgenotypes.con)
      }else{
        hapgenotypes.con  <- file(hapgenotypes, open = "a")
        writeLines(text = hapgenotypes.out, con=hapgenotypes.con)
        close(con = hapgenotypes.con)
      }
      
    }
    #Log message
    if(verbose == TRUE){
      nblocks.done <- nblocks.done + (id2[i]-id1[i]) + 1
      cat(nblocks.done, "blocks written to file\r")
    }
    
  }
  
  
}

