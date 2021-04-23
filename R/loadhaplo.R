#Function: ghap.loadhaplo
#License: GPLv3 or later
#Modification date: 18 Feb 2021
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Load haplotype genotypes

ghap.loadhaplo<-function(
  input.file=NULL,
  hapsamples.file=NULL,
  hapalleles.file=NULL,
  hapgenotypesb.file=NULL,
  verbose = TRUE
){
  
  # Check input file prefix
  if(is.null(input.file) == FALSE){
    hapsamples.file <- paste(input.file, "hapsamples", sep=".")
    hapalleles.file <- paste(input.file, "hapalleles", sep=".")
    hapgenotypesb.file <- paste(input.file, "hapgenotypesb", sep=".")
  }else if(is.null(hapgenotypesb.file)){
    stop("Please provide a binary haplotype genotypes file!")
  }else if(is.null(hapsamples.file)){
    stop("Please provide a samples file!")
  }else if(is.null(hapalleles.file)){
    stop("Please provide a haplotype alleles file!")
  }
  
  #Check files
  if(file.exists(hapgenotypesb.file) == FALSE){
    stop("Could not find binary haplotype genotypes file")
  }
  if(file.exists(hapalleles.file) == FALSE){
    stop("Could not find haplotype alleles file")
  }
  if(file.exists(hapsamples.file) == FALSE){
    stop("Could not find sample file")
  }
  
  #Load haplotype alleles file
  if(verbose == TRUE){
    cat("\nReading in haplotype allele information... ")
  }
  hapalleles <- fread(hapalleles.file, header=FALSE, colClasses = "character")
  
  #Check if the haplotype alleles file contains correct dimension
  if(ncol(hapalleles) != 5){
    stop("[ERROR]\n\nHaplotype alleles file contains wrong number of columns (expected 5)")
  }else{
    if(verbose == TRUE){
      cat("Done.\n")
      cat(paste("A total of ", nrow(hapalleles), " haplotype alleles were found.\n",sep=""))
    }
    hapalleles$V3 <- as.numeric(hapalleles$V3)
    hapalleles$V4 <- as.numeric(hapalleles$V4)
  }

  
  #Load sample file
  if(verbose == TRUE){
    cat("Reading in sample information... ")
  }
  sample <- fread(hapsamples.file, header=FALSE, colClasses = "character")
  
  #Check if the sample file contains correct dimension
  if(ncol(sample) != 2){
    stop("[ERROR]\n\nSample file contains wrong number of columns (expected 2)")
  }
  
  #Check for duplicated ids
  if(length(unique(sample$V2)) < nrow(sample)){
    stop("[ERROR]\n\nSample file contains duplicated ids!")
  }
  
  #Report individuals
  pop <- sample$V1
  ids <- sample$V2
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(sample), " individuals were found in ", length(unique(pop)), " populations.\n",sep=""))
  }
  
  #Create GHap.hap object
  hap<-NULL
  hap$nsamples <- nrow(sample)
  hap$nalleles <- nrow(hapalleles)
  hap$nsamples.in <- nrow(sample)
  hap$nalleles.in <- nrow(hapalleles)
  hap$pop <- pop
  hap$id <- ids
  hap$id.in <- rep(TRUE,times=length(hap$id))
  hap$chr <- hapalleles$V2
  hap$block <- hapalleles$V1
  hap$bp1 <- hapalleles$V3
  hap$bp2 <- hapalleles$V4
  hap$allele <- hapalleles$V5
  hap$allele.in <- rep(TRUE,times=length(hap$allele))
  hap$genotypes <- normalizePath(path = hapgenotypesb.file)
  
  # Check phaseb file
  if(verbose == TRUE){
    cat("Checking integrity of haplotype genotypes... ")
  }
  bitloss <- 8 - ((2*hap$nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  nbytes <- file.info(hapgenotypesb.file)$size
  ebytes <- hap$nalleles*(2*hap$nsamples + bitloss)/8
  ebytes <- ebytes + 3
  if(nbytes != ebytes){
    osize <- hap$nsamples*((nbytes-3)/ceiling((2*hap$nsamples)/8))
    osize <- floor(osize)
    esize <- hap$nsamples*hap$nalleles
    emsg <- "[ERROR]\n\n\nYour binary haplotype genotypes file contains wrong dimensions:\n"
    emsg <- paste(emsg, "Expected ",hap$nsamples,"*",hap$nalleles," = ", esize, " alleles",sep="")
    emsg <- paste(emsg, "but found", osize)
    stop(emsg)
  }else{
    if(verbose == TRUE){
      cat("Done.\n")
    }
  }
  
  #Return GHap.haplo object
  class(hap) <- "GHap.haplo"
  if(verbose == TRUE){
    cat("Your GHap.haplo object was successfully loaded without apparent errors.\n\n")
  }
  return(hap)
  
}