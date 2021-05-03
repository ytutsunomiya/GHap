#Function: ghap.loadplink
#License: GPLv3 or later
#Modification date: 03 May 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Build a GHap.plink object from PLINK binary files

ghap.loadplink<-function(
  input.file = NULL,
  bed.file = NULL,
  bim.file = NULL,
  fam.file = NULL,
  ncores = 1,
  verbose = TRUE
){
  
  # Check stem file name -------------------------------------------------------
  if(is.null(input.file) == FALSE){
    fam.file <- paste(input.file, "fam", sep=".")
    bim.file <- paste(input.file, "bim", sep=".")
    bed.file <- paste(input.file, "bed", sep=".")
  }else if(is.null(bed.file) & is.null(bim.file) & is.null(fam.file)){
    stop("Please provide a valid set of binary PLINK (bed/bim/fam) files!")
  }else if(is.null(bed.file)){
    stop("Please provide a plink bed file!")
  }else if(is.null(bim.file)){
    stop("Please provide a plink bim file!")
  }else if(is.null(fam.file)){
    stop("Please provide a plink fam file!")
  }
  
  # Check if files exist -------------------------------------------------------
  if(file.exists(bed.file) == FALSE){
    stop("Could not find bed file")
  }
  if(file.exists(bim.file) == FALSE){
    stop("Could not find bim file")
  }
  if(file.exists(fam.file) == FALSE){
    stop("Could not find fam file")
  }
  
  # Load bim file --------------------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  if(verbose == TRUE){
    cat("\nReading in marker information... ")
  }
  bim <- fread(bim.file, header = FALSE,
               colClasses = "character", nThread = ncores)
  
  #Check if the bim file contains correct dimension
  if(ncol(bim) != 6){
    emsg <- "\n\nThe bim file contains wrong number of columns (expected 6)"
    stop(emsg)
  }else{
    bim$V3 <- as.numeric(bim$V3)
    bim$V4 <- as.numeric(bim$V4)
  }
  
  # Check if alleles are different ---------------------------------------------
  equalalleles <- length(which(bim$V5 == bim$V6 & bim$V6 != "0"))
  if(equalalleles > 0){
    emsg <- "\n\nThe bim file contains markers with A0 = A1!"
    stop(emsg)
  }
  
  # Check for duplicated marker ids --------------------------------------------
  dup <- which(duplicated(bim$V2))
  ndup <- length(dup)
  if(ndup > 0){
    emsg <- paste("\n\nThe bim file contains", ndup, "duplicated ids")
    stop(emsg)
  }
  
  # Check if markers are sorted by bp ------------------------------------------
  chr <- unique(bim$V1)
  nchr <- length(chr)
  chrorder <- chr[order(nchar(chr),chr)]
  negpos <- diff(bim$V4)
  if("0" %in% chr){
    negpos <- length(which(negpos < 0)) + 2
  }else{
    negpos <- length(which(negpos < 0)) + 1
  }
  if(identical(chr,chrorder) == FALSE | negpos != nchr){
    emsg <- "\n\nMarkers are not sorted by chromosome and base pair position"
    stop(emsg)
  }
  
  # Set to missing chromosome labels and positions with code "0" ---------------
  idx1 <- which(bim$V1 == "0")
  idx2 <- which(bim$V3 == "0") 
  idx3 <- which(bim$V4 == "0")
  note <- NULL
  if(length(idx1) > 0){
    note <- paste(note, "\n[NOTE] Found", length(idx1),
                  "markers with missing chromosome labels!")
  }
  if(length(idx2) > 0){
    note <- paste(note, "\n[NOTE] Found", length(idx2),
                  "markers with missing genetic positions!")
  }
  if(length(idx3) > 0){
    note <- paste(note, "\n[NOTE] Found", length(idx3),
                  "markers with missing physical positions!")
  }
  bim$V1[idx1] <- NA
  bim$V3[idx2] <- NA
  bim$V4[idx3] <- NA
  bim$V3[idx2] <- as.numeric(bim$V4[idx2])/1e+6
  
  # Get report of potentially monomorphic markers ------------------------------
  idx1 <- which(bim$V5 == "0" | bim$V6 == "0")
  if(length(idx1) > 0){
    note <- paste(note, "\n[NOTE] Found", length(idx1),
                  'markers that seem monomorphic (one of the allele codes is "0")!')
  }
  
  # Check for duplicated bp ----------------------------------------------------
  dup <- paste(bim$V1,bim$V4)
  dup <- which(duplicated(dup[which(dup != "NA NA")]))
  ndup <- length(dup)
  if(ndup > 0){
    note <- paste(note, "\n[NOTE] Found", ndup,
                  "duplicated physical positions!")
  }
  
  # Map passed checks ----------------------------------------------------------
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(bim),
              " markers were found in ", nchr," chromosomes.\n",sep=""))
  }
  
  # Load fam file --------------------------------------------------------------
  if(verbose == TRUE){
    cat("Reading in sample information... ")
  }
  fam <- fread(fam.file, header = FALSE,
               colClasses = "character", nThread = ncores)
  
  # Check if the fam file contains correct dimension ---------------------------
  if(ncol(fam) != 6){
    emsg <- "\n\nThe fam file contains wrong number of columns (expected 6)"
    stop(emsg)
  }
  
  # Check for duplicated ids ---------------------------------------------------
  if(length(unique(fam$V2)) < nrow(fam)){
    emsg <- "\n\nThe fam file contains duplicated ids!"
    stop(emsg)
  }
  
  # Report individuals ---------------------------------------------------------
  fids <- fam$V1
  iids <- fam$V2
  fam$V3[which(fam$V3 == "0")] <- NA
  fam$V4[which(fam$V4 == "0")] <- NA
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(fam), " individuals were found in ",
              length(unique(fids)),
              " families (treated as populations by GHap).\n",sep=""))
  }
  
  # Create GHap.plink object ---------------------------------------------------
  sexcode <- c(NA,"M","F")
  names(sexcode) <- c("0","1","2")
  plink <- NULL
  plink$nsamples <- nrow(fam)
  plink$nmarkers <- nrow(bim)
  plink$nsamples.in <- nrow(fam)
  plink$nmarkers.in <- nrow(bim)
  plink$pop <- fids
  plink$id <- iids
  plink$id.in <- rep(TRUE,times=length(plink$id))
  plink$sire <- fam$V3
  plink$dam <- fam$V4
  plink$sex <- as.character(sexcode[fam$V5])
  plink$chr <- bim$V1
  plink$marker <- bim$V2
  plink$marker.in <- rep(TRUE,times=length(plink$marker))
  plink$cm <- bim$V3
  plink$bp <- bim$V4
  plink$A0 <- bim$V6
  plink$A1 <- bim$V5
  plink$plink <- normalizePath(path = bed.file)
  
  # Get report of missing sex and parent information ---------------------------
  idx1 <- which(is.na(plink$sire))
  idx2 <- which(is.na(plink$dam)) 
  idx3 <- which(is.na(plink$sex))
  if(length(idx1) > 0){
    note <- paste(note, "\n[NOTE] Found", length(idx1),
                  "individuals with missing sire id!")
  }
  if(length(idx2) > 0){
    note <- paste(note, "\n[NOTE] Found", length(idx2),
                  "individuals with missing dam id!")
  }
  if(length(idx3) > 0){
    note <- paste(note, "\n[NOTE] Found", length(idx3),
                  "individuals with missing sex code!")
  }
  
  # Check phaseb file ----------------------------------------------------------
  if(verbose == TRUE){
    cat("Checking integrity of bed file... ")
  }
  bitloss <- 8 - ((2*plink$nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  nbytes <- file.info(bed.file)$size
  ebytes <- plink$nmarkers*(2*plink$nsamples + bitloss)/8
  ebytes <- ebytes + 3
  if(nbytes != ebytes){
    osize <- plink$nsamples*((nbytes-3)/ceiling((2*plink$nsamples)/8))
    osize <- floor(osize)
    esize <- plink$nsamples*plink$nmarkers
    emsg <- "\n\nYour bed file contains wrong dimensions:\n"
    emsg <- paste(emsg, "Expected ",plink$nsamples,"*",plink$nmarkers," = ",
                  esize, " markers",sep="")
    emsg <- paste(emsg, "but found", osize)
    stop(emsg)
  }else{
    if(verbose == TRUE){
      cat("Done.\n")
    }
  }
  
  # Return GHap.plink object ---------------------------------------------------
  class(plink) <- "GHap.plink"
  cmsg <- "Your GHap.plink object was successfully loaded"
  if(verbose == TRUE){
    if(is.null(note)){
      cmsg <- paste(cmsg, "without apparent errors.\n\n")
    }else{
      cmsg <- paste0(cmsg, ", but with the following warnings:\n",
                     note)
    }
    cat(cmsg)
  }
  return(plink)
  
}
