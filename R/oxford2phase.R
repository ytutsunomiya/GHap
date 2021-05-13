#Function: ghap.oxford2phase
#License: GPLv3 or later
#Modification date: 13 May 2021
#Written by: Mario Barbato & Yuri Utsunomiya
#Contact: mario.barbato@unicatt.it, ytutsunomiya@gmail.com
#Description: Convert an Oxford HAPS/SAMPLE file to the GHap phase format

ghap.oxford2phase <- function(
  input.files = NULL, 
  haps.files = NULL, 
  sample.files = NULL,
  out.file,
  ncores = 1,
  verbose = TRUE
){
  
  # Check input file prefix
  if(is.null(input.files) == FALSE){
    sample.files <- paste(input.files, "sample", sep=".")
    haps.files <- paste(input.files, "haps", sep=".")
  }else if(is.null(haps.files)){
    stop("Please provide a vector of haps files!")
  }else if(is.null(sample.files)){
    stop("Please provide a vector of sample files!")
  }
  
  # Check if haps files exist
  for(i in 1:length(haps.files)){
    if(file.exists(haps.files[i]) == FALSE){
      haps.files[i] <- paste(haps.files[i], "gz", sep=".")
      if(file.exists(haps.files[i]) == FALSE){
        haps.files[i] <- gsub(pattern = "\\.gz$", replacement = "", x = haps.files[i])
        emsg <- paste("Could not find file ",haps.files[i],"{.gz}",sep="")
        stop(emsg)
      }
    }
  }
  
  # Check if sample files exist
  for(i in 1:length(sample.files)){
    if(file.exists(sample.files[i]) == FALSE){
      sample.files[i] <- paste(sample.files[i], "gz", sep=".")
      if(file.exists(sample.files[i]) == FALSE){
        sample.files[i] <- gsub(pattern = "\\.gz$", replacement = "", x = sample.files[i])
        emsg <- paste("Could not find file ",sample.files[i],"{.gz}",sep="")
        stop(emsg)
      }
    }
  }
  
  # Check matching number of files
  if(length(haps.files) != length(sample.files)){
    stop("Number of haps files and sample files should be equal!")
  }
  
  # Check if out file exist
  samples.file <- paste(out.file,"samples",sep=".")
  phase.file <- paste(out.file,"phase",sep=".")
  markers.file <- paste(out.file,"markers",sep=".")
  if(file.exists(samples.file) == TRUE | file.exists(markers.file) == TRUE | file.exists(markers.file) == TRUE){
    stop("Output file already exists!")
  }else{
    rnumb <- runif(n = 1, min = 1, max = 1e+6)
    rnumb <- ceiling(rnumb)
    tmp.file <- paste(tempdir(),"/tmp",rnumb,sep="")
    tmp.samples.file <- paste(tmp.file,"samples",sep=".")
    tmp.phase.file <- paste(tmp.file,"phase",sep=".")
    tmp.markers.file <- paste(tmp.file,"markers",sep=".")
  }
  
  # Convert sample files
  ncores <- min(c(detectCores(), ncores))
  mysamp <- fread(file = sample.files[1], header = FALSE, skip = 2, drop = 3, nThread = ncores)
  if(length(sample.files) > 1){
    for(i in 2:length(sample.files)){
      tmp.mysamp <- fread(file = sample.files[i], header = FALSE, skip = 2, drop = 3, nThread = ncores)
      if(identical(tmp.mysamp,mysamp) == FALSE){
        stop("Sample files differ!")
      }
    }
  }
  fwrite(x = mysamp[,1:2], file = tmp.samples.file, col.names = FALSE, row.names = FALSE, sep = " ", nThread = ncores)
  
  # Convert haps files
  nids <- nrow(mysamp)
  expcols <- 2*nids + 5
  if(verbose == TRUE){
    cat("\nConverting Oxford files to the GHap phase format.\n")
  }
  for(i in 1:length(haps.files)){
    myfile <- fread(file = haps.files[i], header = FALSE, verbose = FALSE, showProgress = FALSE, nThread = ncores)
    if(ncol(myfile) != expcols){
      emsg <- paste("Expected 5 + 2*",nids," = ",expcols,
                    " columns in file ",haps.files[i]," but found ", ncol(myfile), "!", sep="")
      stop(emsg)
    }
    fwrite(x = myfile[,1:5], file = tmp.markers.file, nThread = ncores,
           col.names = FALSE, row.names = FALSE, sep = " ", append = TRUE)
    fwrite(x = myfile[,-c(1:5)], file = tmp.phase.file, nThread = ncores,
           col.names = FALSE, row.names = FALSE, sep = " ", append = TRUE)
    if(verbose == TRUE){
      cat("Processed ",i," files of",length(haps.files), "\r")
    }
  }
  
  # Get files
  if(verbose == TRUE){
    cat("Copying files to the working directory... ")
  }
  ok <- file.copy(from = tmp.phase.file, to = phase.file)
  ok <- file.remove(tmp.phase.file)
  ok <- file.copy(from = tmp.markers.file, to = markers.file)
  ok <- file.remove(tmp.markers.file)
  ok <- file.copy(from = tmp.samples.file, to = samples.file)
  ok <- file.remove(tmp.samples.file)
  if(verbose == TRUE){
    cat("Done.\n\n")
  }
  
}
