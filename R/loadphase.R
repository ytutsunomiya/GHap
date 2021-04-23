#Function: ghap.loadphase
#License: GPLv3 or later
#Modification date: 18 Feb 2021
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Load phased genotypes

ghap.loadphase <- function(
  input.file = NULL,
  samples.file = NULL,
  markers.file = NULL,
  phaseb.file = NULL,
  verbose = TRUE
){
  
  # Check stem file name
  if(is.null(input.file) == FALSE){
    samples.file <- paste(input.file, "samples", sep=".")
    markers.file <- paste(input.file, "markers", sep=".")
    phaseb.file <- paste(input.file, "phaseb", sep=".")
  }else if(is.null(phaseb.file)){
    stop("Please provide a binary phase file!")
  }else if(is.null(samples.file)){
    stop("Please provide a samples file!")
  }else if(is.null(markers.file)){
    stop("Please provide a markers file!")
  }
  
  #Check files
  if(file.exists(phaseb.file) == FALSE){
    stop("Could not find phased genotypes file")
  }
  if(file.exists(markers.file) == FALSE){
    stop("Could not find marker map file")
  }
  if(file.exists(samples.file) == FALSE){
    stop("Could not find sample file")
  }
  
  #Load marker map file
  if(verbose == TRUE){
    cat("\nReading in marker map information... ")
  }
  marker <- fread(markers.file, header=FALSE, colClasses = "character")
  
  #Check if the map file contains correct dimension
  if(ncol(marker) != 5){
    stop("[ERROR]\n\n Marker map contains wrong number of columns (expected 5)")
  }
  marker$V3 <- as.numeric(marker$V3)
  
  #Check if alleles are different
  equalalleles <- length(which(marker$V4 == marker$V5))
  if(equalalleles > 0){
    stop("[ERROR]\n\n The map contains markers with A0 = A1!")
  }
  
  #Check for duplicated marker ids
  dup <- which(duplicated(marker$V2))
  ndup <- length(dup)
  if(ndup > 0){
    emsg <- paste("[ERROR]\n\nYour marker map file contains", ndup, "duplicated ids")
    stop(emsg)
  }
  
  #Check if markers are sorted by bp
  chr <- unique(marker$V1)
  nchr <- length(chr)
  chrorder <- chr[order(nchar(chr),chr)]
  negpos <- diff(marker$V3)
  negpos <- length(which(negpos < 0)) + 1
  if(identical(chr,chrorder) == FALSE | negpos != nchr){
    stop("[ERROR]\n\nMarkers are not sorted by chromosome and base pair position")
  }
  
  #Check for duplicated bp
  dup <- paste(marker$V1,marker$V3)
  dup <- which(duplicated(dup))
  ndup <- length(dup)
  if(ndup > 0){
    emsg <- paste("[ERROR]\n\nYour marker map file contains", ndup, "duplicated positions!")
    stop(emsg)
  }
  
  #Map passed checks
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(marker),
              " markers were found in ", nchr," chromosomes.\n",sep=""))
  }
  
  #Load sample file
  if(verbose == TRUE){
    cat("Reading in sample information... ")
  }
  sample <- fread(samples.file, header=FALSE, colClasses = "character")
  
  #Check if the sample file contains correct dimension
  if(ncol(sample) != 2){
    stop("[ERROR]\n\nSample file contains wrong number of columns (expected 2)")
  }

  #Check for duplicated ids
  dup <- which(duplicated(sample$V2))
  ndup <- length(dup)
  if(ndup > 0){
    emsg <- paste("[ERROR]\n\nSample file contains", ndup, "duplicated ids!")
    stop(emsg)
  }
  
  # Samples passed check
  pop <- rep(sample$V1,each=2)
  ids <- rep(sample$V2,each=2)
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(sample), " individuals were found in ",
              length(unique(pop)), " populations.\n",sep=""))
  }
  
  #Create GHap2.phase object
  phase <- NULL
  phase$chr <- marker$V1
  phase$nsamples <- nrow(sample)
  phase$nmarkers <- nrow(marker)
  phase$nsamples.in <- nrow(sample)
  phase$nmarkers.in <- nrow(marker)
  phase$pop <- pop
  phase$id <- ids
  phase$id.in <- rep(TRUE,times=length(phase$id))
  phase$marker <- marker$V2
  phase$marker.in <- rep(TRUE,times=length(phase$marker))
  phase$bp <- marker$V3
  phase$A0 <- marker$V4
  phase$A1 <- marker$V5
  phase$phase <- normalizePath(path = phaseb.file)
  
  # Check phaseb file
  if(verbose == TRUE){
    cat("Checking integrity of phased genotypes... ")
  }
  bitloss <- 8 - ((2*phase$nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  nbytes <- file.info(phaseb.file)$size
  ebytes <- phase$nmarkers*(2*phase$nsamples + bitloss)/8
  if(nbytes != ebytes){
    osize <- 2*phase$nsamples*(nbytes/ceiling((2*phase$nsamples)/8))
    osize <- floor(osize)
    esize <- 2*phase$nsamples*phase$nmarkers
    emsg <- "[ERROR]\n\n\nYour binary phased genotypes file contains wrong dimensions:\n"
    emsg <- paste(emsg, "Expected 2*",phase$nsamples,"*",phase$nmarkers," = ", esize, " alleles",sep="")
    emsg <- paste(emsg, "but found", osize)
    stop(emsg)
  }else{
    if(verbose == TRUE){
      cat("Done.\n") 
    }
  }
  
  #Return GHap object
  class(phase) <- "GHap.phase"
  if(verbose == TRUE){
    cat("Your GHap.phase object was successfully created without apparent errors.\n\n")
  }
  return(phase)
  
}