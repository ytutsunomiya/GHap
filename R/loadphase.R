#Function: ghap.loadphase
#License: GPLv3 or later
#Modification date: 29 Aug 2022
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Load phased genotypes

ghap.loadphase <- function(
  input.file = NULL,
  samples.file = NULL,
  markers.file = NULL,
  phaseb.file = NULL,
  ncores = 1,
  verbose = TRUE
){
  
  # Check stem file name -------------------------------------------------------
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
  
  # Check files ----------------------------------------------------------------
  if(file.exists(phaseb.file) == FALSE){
    stop("Could not find phased genotypes file")
  }
  if(file.exists(markers.file) == FALSE){
    stop("Could not find marker map file")
  }
  if(file.exists(samples.file) == FALSE){
    stop("Could not find sample file")
  }
  
  # Load marker map file -------------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  if(verbose == TRUE){
    cat("\nReading in marker map information... ")
  }
  marker <- fread(markers.file, header=FALSE,
                  colClasses = "character", nThread = ncores)
  
  # Check if the map file contains correct dimension ---------------------------
  if(ncol(marker) %in% c(5,6) == FALSE){
    stop("\n\nMarker map contains wrong number of columns (expected 5 or 6)")
  }
  marker$V3 <- as.numeric(marker$V3)
  if(ncol(marker) == 5){
    tmp <- as.data.frame(matrix(data = NA, nrow = nrow(marker),
                                ncol = 6))
    colnames(tmp) <- paste0("V",1:6)
    tmp[,1:3] <- marker[,1:3]
    idx <- which(is.na(tmp$V4))
    tmp$V4[idx] <- as.numeric(tmp$V3[idx])/1e+6
    tmp[,5:6] <- marker[,4:5]
    marker <- tmp
  }else{
    marker$V4 <- as.numeric(marker$V4)
  }
  
  # Check if alleles are different ---------------------------------------------
  equalalleles <- length(which(marker$V5 == marker$V6))
  if(equalalleles > 0){
    stop("\n\nThe map contains markers with A0 = A1!")
  }
  
  # Check for duplicated marker ids --------------------------------------------
  dup <- which(duplicated(marker$V2))
  ndup <- length(dup)
  if(ndup > 0){
    emsg <- paste("\n\nYour marker map file contains", ndup, "duplicated ids")
    stop(emsg)
  }
  
  # Check if markers are sorted by bp ------------------------------------------
  chr <- unique(marker$V1)
  nchr <- length(chr)
  chrorder <- chr[order(nchar(chr),chr)]
  negpos <- diff(marker$V3)
  negpos <- length(which(negpos < 0)) + 1
  if(identical(chr,chrorder) == FALSE | negpos != nchr){
    stop("\n\nMarkers are not sorted by chromosome and base pair position")
  }
  
  # Check for duplicated bp ----------------------------------------------------
  dup <- paste(marker$V1,marker$V3)
  dup <- which(duplicated(dup))
  ndup <- length(dup)
  note <- NULL
  if(ndup > 0){
    note <- paste(note, "\n[NOTE] Found", ndup,
                  "duplicated physical positions!")
  }
  
  # Map passed checks ----------------------------------------------------------
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(marker),
              " markers were found in ", nchr," chromosomes.\n",sep=""))
  }
  
  # Load sample file -----------------------------------------------------------
  if(verbose == TRUE){
    cat("Reading in sample information... ")
  }
  sample <- fread(samples.file, header=FALSE,
                  colClasses = "character", nThread = ncores)
  sample <- as.data.frame(sample)
  
  # Check if the sample file contains correct dimension ------------------------
  if(ncol(sample) %in% 2:5 == FALSE){
    stop("\n\nSample file contains wrong number of columns (expected 2 to 5)")
  }
  if(ncol(sample) < 5){
    tmp <- as.data.frame(matrix(data = NA, nrow = nrow(sample), ncol = 5))
    for(i in 1:ncol(sample)){
      tmp[,i] <- sample[,i]
    }
    colnames(tmp) <- paste0("V",1:5)
    sample <- tmp
  }
  sample$V3[which(sample$V3 == "0")] <- NA
  sample$V4[which(sample$V4 == "0")] <- NA
  sample$V5[which(is.na(sample$V5))] <- "0"
  
  # Check for duplicated ids ---------------------------------------------------
  dup <- which(duplicated(sample$V2))
  ndup <- length(dup)
  if(ndup > 0){
    emsg <- paste("\n\nSample file contains", ndup, "duplicated ids!")
    stop(emsg)
  }
  
  # Samples passed check -------------------------------------------------------
  pop <- rep(sample$V1,each=2)
  ids <- rep(sample$V2,each=2)
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(sample), " individuals were found in ",
              length(unique(pop)), " populations.\n",sep=""))
  }
  
  # Create GHap.phase object ---------------------------------------------------
  sexcode <- c(NA,"M","F")
  names(sexcode) <- c("0","1","2")
  phase <- NULL
  phase$nsamples <- nrow(sample)
  phase$nmarkers <- nrow(marker)
  phase$nsamples.in <- nrow(sample)
  phase$nmarkers.in <- nrow(marker)
  phase$pop <- pop
  phase$id <- ids
  phase$id.in <- rep(TRUE,times=length(phase$id))
  phase$sire <- sample$V3
  phase$dam <- sample$V4
  phase$sex <- as.character(sexcode[sample$V5])
  phase$chr <- marker$V1
  phase$marker <- marker$V2
  phase$marker.in <- rep(TRUE,times=length(phase$marker))
  phase$bp <- marker$V3
  phase$cm <- marker$V4
  phase$A0 <- marker$V5
  phase$A1 <- marker$V6
  phase$phase <- normalizePath(path = phaseb.file)
  phase$mode <- NA
  
  # Check phaseb file ----------------------------------------------------------
  if(verbose == TRUE){
    cat("Checking integrity of phased genotypes... ")
  }
  bitloss <- rep(NA, times = 2)
  ebytes <- rep(NA, times = 2)
  bitloss[1] <- 8 - ((2*phase$nsamples) %% 8)
  if(bitloss[1] == 8){
    bitloss[1] <- 0
  }
  bitloss[2] <- 8 - (phase$nmarkers %% 8)
  if(bitloss[2] == 8){
    bitloss[2] <- 0
  }
  nbytes <- file.info(phase$phase)$size
  ebytes[1] <- phase$nmarkers*(2*phase$nsamples + bitloss[1])/8
  ebytes[2] <- 2*phase$nsamples*(phase$nmarkers + bitloss[2])/8
  if(nbytes == ebytes[1]){
    phase$mode <- 0
    if(verbose == TRUE){
      cat("Done.\n")
      cat("Phase object format is mode 0 (variants x individuals).\n")
    }
  }else if(nbytes == ebytes[1] + 1 | nbytes == ebytes[2] + 1){
    phase.con <- file(phase$phase, "rb")
    a <- seek(con = phase.con, where = 0,
              origin = 'start',rw = 'r')
    fmt <- readBin(phase.con, what=raw(), size = 1,
                   n = 1, signed = FALSE)
    fmt <- paste(as.integer(rawToBits(fmt)), collapse="")
    close.connection(phase.con)
    if(fmt == "00000001"){
      phase$mode <- 1
      if(verbose == TRUE){
        cat("Done.\n")
        cat("Phase object format is mode 1 (variants x individuals).\n")
      }
    }else if(fmt == "0000010"){
      phase$mode <- 2
      if(verbose == TRUE){
        cat("Done.\n")
        cat("Phase object format is mode 2 (individuals x variants).\n")
      }
    }else{
      emsg <- "\n\nYour binary phased genotypes file has an unknown format!\n"
      stop(emsg)
    }
  }else{
    emsg <- "\n\nYour binary phased genotypes file contains wrong dimensions:\n"
    emsg <- paste(emsg,
                  "\nExpected", ebytes[1]+1, "bytes for variants x individuals",
                  "\nExpected", ebytes[2]+1, "bytes for individuals x variants",
                  "\nObserved", nbytes, "bytes instead\n\n")
    stop(emsg)
  }
  
  
  # Return GHap object ---------------------------------------------------------
  class(phase) <- "GHap.phase"
  if(verbose == TRUE){
    cat("Your GHap.phase object was successfully created without apparent errors.\n\n")
  }
  return(phase)
  
}
