#Function: ghap.compress
#License: GPLv3 or later
#Modification date: 18 Nov 2022
#Written by: Yuri Tani Utsunomiya, Adam Taiti Harth Utsunomiya
#Contact: ytutsunomiya@gmail.com, adamtaiti@gmail.com
#Description: Compress phased data into GHap binary

ghap.compress <- function(
  input.file=NULL,
  out.file,
  samples.file=NULL,
  markers.file=NULL,
  phase.file=NULL,
  mode=1,
  ncores=1,
  verbose=TRUE
){
  
  # Check input file prefix ----------------------------------------------------
  if(is.null(input.file) == FALSE){
    samples.file <- paste(input.file, "samples", sep=".")
    markers.file <- paste(input.file, "markers", sep=".")
    phase.file <- paste(input.file, "phase", sep=".")
  }else if(is.null(phase.file)){
    stop("Please provide a phase file!")
  }else if(is.null(samples.file)){
    stop("Please provide a samples file!")
  }else if(is.null(markers.file)){
    stop("Please provide a markers file!")
  }
  
  # Check if phase file exist
  if(file.exists(phase.file) == FALSE){
    stop("Could not find the phase file!")
  }
  
  # Check if samples file exist
  if(file.exists(samples.file) == FALSE){
    stop("Could not find the samples file!")
  }
  
  # Check if markers file exist
  if(file.exists(markers.file) == FALSE){
    stop("Could not find the markers file!")
  }
  
  # Check if out file exist
  if(file.exists(paste(out.file,".phaseb",sep="")) == TRUE){
    stop("Output file already exists!")
  }else{
    rnumb <- runif(n = 1, min = 1, max = 1e+6)
    rnumb <- ceiling(rnumb)
    tmp.file <- paste(tempdir(),"/tmp",rnumb,sep="")
  }
  
  # Check file mode
  if(mode %in% 0:2 == FALSE){
    emsg <- paste0("\n\nUnrecognized file mode. Please use one of:\n",
                   "\nmode = 0 (variants x individuals, backward compatibility)",
                   "\nmode = 1 (variants x individuals, default)",
                   "\nmode = 2 (individuals x variants)\n\n")
    stop(emsg)
  }else{
    pmode <- c("mode = 0 (variants x individuals, backward compatibility)",
               "mode = 1 (variants x individuals, default)",
               "mode = 2 (individuals x variants)")
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
  nmarkers <- nrow(marker)
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nmarkers,
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
  nsamples <- nrow(sample)
  pop <- rep(sample$V1,each=2)
  ids <- rep(sample$V2,each=2)
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nsamples, " individuals were found in ",
              length(unique(pop)), " populations.\n",sep=""))
    cat(paste0("Using file ", pmode[mode+1], ".\n"))
  }
  
  # Compute bit loss -----------------------------------------------------------
  if(mode %in% c(0,1)){
    bitloss <- 8 - ((2*nsamples) %% 8)
    if(bitloss == 8){
      bitloss <- 0
    }
    linelen <- 2*nsamples
    nlines <- nmarkers
  }else if(mode == 2){
    bitloss <- 8 - ((2*nmarkers) %% 8)
    if(bitloss == 8){
      bitloss <- 0
    }
    linelen <- 2*nmarkers
    nlines <- nsamples
  }
  
  # Process line function ------------------------------------------------------
  status <- .compressCpp(infile = phase.file, outfile = tmp.file,
                         nunits = nlines, nbits = linelen,
                         tbits = bitloss, fmode = mode)
  if(status == 1){
    emsg <- "Please check the format of your phase file.\n\n"
    stop(emsg)
  }
  if(verbose == TRUE){
    if(mode %in% c(0,1)){
      cat(nmarkers, "markers written to file\n\n")
    }else if(mode == 2){
      cat(nsamples, "individuals written to file\n\n")
    }
  }
  
  # Output results -------------------------------------------------------------
  sup <- file.copy(from = tmp.file, to = paste(out.file,".phaseb",sep=""))
  sup <- file.remove(tmp.file)
  if(verbose == TRUE){
    cat("Phase file succesfully compressed.\n\n")
  }
  
  
}
