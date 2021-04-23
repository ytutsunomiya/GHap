#Function: ghap.compress
#License: GPLv3 or later
#Modification date: 18 Feb 2021
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Compress Oxford phased data into GHap binary

ghap.compress <- function(
  input.file=NULL,
  out.file,
  samples.file=NULL,
  markers.file=NULL,
  phase.file=NULL,
  verbose=TRUE
){
  
  # Check input file prefix
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
  
  #Load marker map file
  if(verbose == TRUE){
    cat("\nReading in marker map information... ")
  }
  marker <- fread(markers.file, header=FALSE, colClasses = "character")
  
  #Check if the map file contains correct dimension
  if(ncol(marker) != 5){
    stop("[ERROR]\n\nMarker map contains wrong number of columns (expected 5)")
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
  nmarkers <- nrow(marker)
  percent <- round(nmarkers/20, digits = 0)
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nmarkers,
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
  nsamples <- nrow(sample)
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nsamples, " individuals were found in ",
              length(unique(pop)), " populations.\n\n",sep=""))
  }
  
  # Open connection with phase file
  phase.con <- file(phase.file,"r")
  
  # Set parameters for connection
  i <- 0
  bitloss <- 8 - ((2*nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  linelen <- 2*nsamples
  p <- 0
  
  # Read haps file
  while(length(line <- readLines(con = phase.con, n = 1)) > 0){
    
    # Update line number
    i <- i + 1
    
    # Break line elements
    line <- scan(text=line, what = "character", sep=" ", quiet = TRUE)
    
    # Check line integrity
    if(length(line) != linelen){
      emsg <- paste("\n\nExpected", linelen, "columns in line",i,"of",
                    phase.file,"but found",length(line),"\n\n")
      stop(emsg)
    }
    
    # Check integrity of alleles
    strings <- unique(line)
    if(length(which(strings %in% c("0","1") == FALSE)) > 0){
      stop("Phased genotypes should be coded as 0 and 1")
    }
    
    # Convert alleles to bits
    line <- c(line,rep("0",times=bitloss))
    line <- paste(line, collapse = "")
    nc <- nchar(line)
    n <- seq(1, nc, by = 8)
    line <- substring(line, n, c(n[-1]-1, nc))
    line <- strtoi(line, base=2)
    
    # Write to output file
    out.con <- file(tmp.file, "ab")
    writeBin(object = line, con = out.con, size = 1)
    close.connection(out.con)
    
    # Report progress
    if(verbose == TRUE){
      if(i == 1 | i == p*percent){
        prog <- paste("Compression status: [",
                      paste(rep("=",times=p),collapse="", sep=""),
                      paste(rep(" ", times=20-p), collapse="", sep=""),
                      "] ", 5*p, "%", collpase="", sep="")
        cat(prog, "\r")
        p <- p + 1
      }
    }
    
    
  }
  
  # Close connection with phase file
  close.connection(phase.con)
  
  # Report progress
  if(verbose == TRUE){
    p <- 20
    prog <- paste("Compression status: [",
                  paste(rep("=",times=p),collapse="", sep=""),
                  paste(rep(" ", times=20-p), collapse="", sep=""),
                  "] ", 5*p, "%", collpase="", sep="")
    cat(prog, "\n")
  }
  
  # Check if number of markers is ok
  if(nmarkers != i){
    file.remove(tmp.file)
    emsg <- paste("\n\nExpected", nmarkers, "lines in",
                  phase.file,"but found",i,"\n\n")
    stop(emsg)
  }else{
    file.copy(from = tmp.file, to = paste(out.file,".phaseb",sep=""))
    file.remove(tmp.file)
  }
  
}

