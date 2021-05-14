#Function: ghap.compress
#License: GPLv3 or later
#Modification date: 13 May 2021
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Compress phased data into GHap binary

ghap.compress <- function(
  input.file=NULL,
  out.file,
  samples.file=NULL,
  markers.file=NULL,
  phase.file=NULL,
  batchsize=NULL,
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
  }
  
  # Compute bit loss -----------------------------------------------------------
  bitloss <- 8 - ((2*nsamples) %% 8)
  if(bitloss == 8){
    bitloss <- 0
  }
  linelen <- 2*nsamples
  
  # Generate batch index -------------------------------------------------------
  if(is.null(batchsize) == TRUE){
    batchsize <- ceiling(nmarkers/10)
  }
  if(batchsize > nmarkers){
    batchsize <- nmarkers
  }
  id1<-seq(1,nmarkers,by=batchsize)
  id2<-(id1+batchsize)-1
  id1<-id1[id2<=nmarkers]
  id2<-id2[id2<=nmarkers]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,nmarkers)
  if(id1[length(id1)] > nmarkers){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  if(verbose == TRUE){
    cat("Processing ", nmarkers, " markers in:\n", sep="")
    batch <- table((id2-id1)+1)
    for(i in 1:length(batch)){
      cat(batch[i]," batches of ",names(batch[i]),"\n",sep="")
    }
  }
  
  # Process line function ------------------------------------------------------
  lineprocess <- function(i){
    line <- scan(text=batchline[i], what = "character", sep=" ", quiet = TRUE)
    if(length(line) != linelen){
      emsg <- paste("\n\nExpected", linelen, "columns in line",i+nlines.skip[b],"of",
                    phase.file,"but found",length(line),"\n\n")
      stop(emsg)
    }
    line <- c(line,rep("0",times=bitloss))
    strings <- unique(line)
    if(length(which(strings %in% c("0","1") == FALSE)) > 0){
      stop("Phased genotypes should be coded as 0 and 1")
    }else{
      line <- paste(line, collapse = "")
      nc <- nchar(line)
      n <- seq(1, nc, by = 8)
      line <- substring(line, n, c(n[-1]-1, nc))
      line <- strtoi(line, base=2)
    }
    return(line)
  }
  
  # Iterate batches ------------------------------------------------------------
  phase.con <- file(phase.file,"r")
  nmarkers.done <- 0
  nlines.read <- id2-id1+1
  nlines.skip <- c(0,cumsum(nlines.read)[-length(nlines.read)])
  for(b in 1:length(id1)){
    
    # Load batch
    batchline <- readLines(con = phase.con, n = nlines.read[b])
    
    # Check if batch is ok
    if(length(batchline) != nlines.read[b]){
      emsg <- paste("\n\nExpected", nmarkers, "lines in",
                    phase.file,"but found", length(batchline)+nlines.skip[b],"\n\n")
      stop(emsg)
    }
    
    # Transform lines
    ncores <- min(c(detectCores(), ncores))
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      results <- unlist(parLapply(cl = cl, fun = lineprocess, X = 1:length(batchline)))
      stopCluster(cl)
    }else{
      results <- unlist(mclapply(FUN = lineprocess, X = 1:length(batchline), mc.cores = ncores))
    }
    
    # Write to output file
    out.con <- file(tmp.file, "ab")
    writeBin(object = results, con = out.con, size = 1)
    close.connection(out.con)
    
    # Log message
    if(verbose == TRUE){
      nmarkers.done <- nmarkers.done + (id2[i]-id1[i]) + 1
      cat(nmarkers.done, "markers written to file\r")
    }
    
  }
  
  # Last integrity check -------------------------------------------------------
  batchline <- readLines(con = phase.con, n = 1)
  if(length(batchline) != 0){
    emsg <- paste("\n\nExpected", nmarkers, "lines in",
                  phase.file,"but found more\n\n")
    stop(emsg)
  }
  
  # Close connection with phase file -------------------------------------------
  close.connection(phase.con)
  if(verbose == TRUE){
    cat(nmarkers, "markers written to file\n\n")
  }
  
  # Output results -------------------------------------------------------------
  sup <- file.copy(from = tmp.file, to = paste(out.file,".phaseb",sep=""))
  sup <- file.remove(tmp.file)
  if(verbose == TRUE){
    cat("Phase file succesfully compressed.\n\n")
  }
  
  
}
