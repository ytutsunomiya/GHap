#Function: ghap.fast2phase
#License: GPLv3 or later
#Modification date: 2 Jun 2022
#Written by: Mario Barbato, Yuri Tani Utsunomiya
#Contact: mario.barbato@unicatt.it, ytutsunomiya@gmail.com
#Description: Converts fastPHASE to GHap phase format

ghap.fast2phase <- function(
  input.files = NULL, 
  switchout.files = NULL,
  map.files = NULL,
  fam.file = NULL,
  out.file = NULL,
  overwrite = FALSE,
  ncores = 1,
  verbose = TRUE
){
  
  # Check input file prefix
  if(is.null(input.files) == FALSE){
    map.files <- paste(input.files, "map", sep=".")
    switchout.files <- paste(input.files, "_switch.out", sep="")
  }else if(is.null(switchout.files)){
    stop("Please provide a vector of _switch.out files!")
  }else if(is.null(map.files)){
    stop("Please provide a vector of map files!")
  }else if(is.null(fam.file)){
    stop("Please provide a fam file!")
  }else if(is.null(out.file)){
    stop("Please provide a valid output file name!")
  }
  
  # Check if switchout files exist
  for(i in 1:length(switchout.files)){
    if(file.exists(switchout.files[i]) == FALSE){
      switchout.files[i] <- paste(switchout.files[i], "gz", sep=".")
      if(file.exists(switchout.files[i]) == FALSE){
        switchout.files[i] <- gsub(pattern = "\\.gz$", replacement = "", x = switchout.files[i])
        emsg <- paste("Could not find file ",switchout.files[i],"{.gz}",sep="")
        stop(emsg)
      }
    }
  }
  
  # Check if fam files exist
  if(file.exists(fam.file) == FALSE){
    fam.file <- paste(fam.file, "gz", sep=".")
    if(file.exists(fam.file) == FALSE){
      fam.file <- gsub(pattern = "\\.gz$", replacement = "", x = fam.file)
      emsg <- paste("Could not find file ", fam.file,"{.gz}", sep="")
      stop(emsg)
    }
  }
  
  # Check if map files exist
  for(i in 1:length(map.files)){
    if(file.exists(map.files[i]) == FALSE){
      map.files[i] <- paste(map.files[i], "gz", sep=".")
      if(file.exists(map.files[i]) == FALSE){
        map.files[i] <- gsub(pattern = "\\.gz$", replacement = "", x = map.files[i])
        emsg <- paste("Could not find file ",map.files[i],"{.gz}",sep="")
        stop(emsg)
      }
    }
  }
  
  # Check matching number of files
  if(length(switchout.files) != length(map.files)){
    stop("Number of _switch.out files and map files should be equal!")
  }
  
  # Check if out file exist
  samples.file <- paste(out.file,"samples",sep=".")
  phase.file <- paste(out.file,"phase",sep=".")
  markers.file <- paste(out.file,"markers",sep=".")
  if(overwrite == FALSE & any(file.exists(samples.file) == TRUE, file.exists(markers.file) == TRUE, file.exists(markers.file) == TRUE)){
    stop("Output file already exists!")
  }else{
    rnumb <- runif(n = 1, min = 1, max = 1e+6)
    rnumb <- ceiling(rnumb)
    tmp.file <- paste(tempdir(),"/tmp",rnumb,sep="")
    tmp.samples.file <- paste(tmp.file,"samples",sep=".")
    tmp.phase.file <- paste(tmp.file,"phase",sep=".")
    tmp.markers.file <- paste(tmp.file,"markers",sep=".")
  }
  
  # Process files
  if(verbose == TRUE){
    cat("\nConverting fastPHASE files to the GHap phase format.\n")
  }
  
  # Get fam info
  ncores <- min(c(detectCores(), ncores))
  samples.fam <- fread(fam.file, select = 1:2, colClasses = "character", nThread = ncores)
  
  for(i in 1:length(switchout.files)){
    
    # Load marker map file
    marker <- fread(map.files[i], header=FALSE, colClasses = "character", nThread = ncores)
    
    # Check if the map file contains correct dimension
    if(ncol(marker) != 5){
      emsg <- paste("File", map.files[i], "contains wrong number of columns (expected 5)")
      stop(emsg)
    }
    marker$V3 <- as.numeric(marker$V3)
    
    #Check if alleles are different
    equalalleles <- length(which(marker$V4 == marker$V5))
    if(equalalleles > 0){
      emsg <- paste("File", map.files[i], "countains markers with A0 = A1!")
      stop(emsg)
    }
    
    # Check for duplicated marker ids
    dup <- which(duplicated(marker[,2]))
    ndup <- length(dup)
    if(ndup > 0){
      emsg <- paste("File", map.files[i], "contains", ndup, "duplicated marker ids")
      stop(emsg)
    }
    
    # Check if markers are sorted by bp
    chr <- unique(marker[,1])$V1
    nchr <- length(chr)
    chrorder <- chr[order(nchar(chr),chr)]
    negpos <- diff(marker[,3])
    negpos <- length(which(negpos < 0)) + 1
    if(identical(chr,chrorder) == FALSE | negpos != nchr){
      emsg <- paste("Markers are not sorted by chromosome and base pair position in file", map.files[i])
      stop(emsg)
    }
    
    # Check for duplicated bp
    dup <- paste(marker[,1],marker[,3])
    dup <- which(duplicated(dup))
    ndup <- length(dup)
    if(ndup > 0){
      emsg <- paste("File", map.files[i], "contains", ndup, "duplicated positions")
      stop(emsg)
    }
    
    # Load fastPHASE file
    fp <- fread(switchout.files[i], fill = TRUE, header = FALSE, nThread = ncores)
    
    # Map samples and haplotypes
    for(j in 1:nrow(fp)){if(fp[j] == "BEGIN GENOTYPES"){break()}}
    idx <- c(j + 1, seq(j + 4, to = nrow(fp) - 1, by = 3))
    hapidx <- which(1:(nrow(fp)-1) %in% c(1:j,idx) == FALSE)
    
    # Get ID info from files
    samples <- stri_match(str = fp$V1[idx], regex = "(?<=\\bID\\s)\\S+\\b")
    
    # Check if samples are the same across files
    if(is.null(samples.fam) == FALSE){
      if(identical(samples.fam$V2, as.character(samples)) == FALSE){
        stop("Samples differ between files!")
      }
    }
    
    # Isolate haplotypes
    fp_split <- strsplit(fp$V1[hapidx], split = " ")
    fp_split <- t(matrix(unlist(fp_split), ncol = length(fp_split[[1]]), byrow = T))
    
    # Check number of markers
    if(nrow(fp_split) != nrow(marker)){
      emsg <- paste("Expected",nrow(marker),"markers in file", switchout.files[i],
                    "but found",nrow(fp_split),"!")
      stop(emsg)
    }
    
    # Check number of samples
    if(ncol(fp_split) != 2*nrow(samples)){
      emsg <- paste("Expected",2*nrow(samples),"haplotypes in file", switchout.files[i],
                    "but found",ncol(fp_split),"!")
      stop(emsg)
    }
    
    # Recode alleles
    ccx <- 0
    for(k in seq(1, ncol(fp_split), by=2)){
      ccx = ccx + 1
      fp_split[fp_split[, k] == marker$V5[ccx], k] <- 0
      fp_split[fp_split[, k+1] == marker$V5[ccx], k+1] <- 0
    }
    fp_split[fp_split != "0"] <- 1
    
    # Write to file
    fwrite(x = marker, file = tmp.markers.file, nThread = ncores,
           col.names = FALSE, row.names = FALSE, sep = " ", append = TRUE)
    fwrite(x = as.data.table(fp_split), file = tmp.phase.file, nThread = ncores,
           col.names = FALSE, row.names = FALSE, sep = " ", append = TRUE, verbose = FALSE)
    
    # Log message
    if(verbose == TRUE){
      cat("Processed ",i," files of",length(switchout.files), "\r")
    }
    
  }
  
  # Output samples file
  fwrite(x = as.data.table(samples.fam), file = tmp.samples.file, col.names = FALSE,
         row.names = FALSE, sep = " ", nThread = ncores)
  
  # Get files
  if(verbose == TRUE){
    cat("Copying files to the working directory...")
  }
  ok <- file.copy(from = tmp.phase.file, to = phase.file, overwrite = overwrite)
  ok <- file.remove(tmp.phase.file)
  ok <- file.copy(from = tmp.markers.file, to = markers.file, overwrite = overwrite)
  ok <- file.remove(tmp.markers.file)
  ok <- file.copy(from = tmp.samples.file, to = samples.file, overwrite = overwrite)
  ok <- file.remove(tmp.samples.file)
  if(verbose == TRUE){
    cat("Done.\n\n")
  } 
}
