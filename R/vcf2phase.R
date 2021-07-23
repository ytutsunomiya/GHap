#Function: ghap.vcf2phase
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Utsunomiya
#Contact: mario.barbato@unicatt.it, ytutsunomiya@gmail.com
#Description: Convert VCF files to the GHap phase format

ghap.vcf2phase <- function(
  input.files = NULL, 
  vcf.files = NULL, 
  sample.files = NULL,
  out.file,
  ncores = 1,
  verbose = TRUE
){
  
  # Check input file prefix
  if(is.null(input.files) == FALSE){
    sample.files <- paste(input.files, "sample", sep=".")
    vcf.files <- paste(input.files, "vcf", sep=".")
  }else if(is.null(vcf.files)){
    stop("Please provide a vector of vcf files!")
  }else if(is.null(sample.files)){
    stop("Please provide a vector of sample files!")
  }
  
  # Check if vcf files exist
  for(i in 1:length(vcf.files)){
    if(file.exists(vcf.files[i]) == FALSE){
      vcf.files[i] <- paste(vcf.files[i], "gz", sep=".")
      if(file.exists(vcf.files[i]) == FALSE){
        vcf.files[i] <- gsub(pattern = "\\.gz$", replacement = "", x = vcf.files[i])
        emsg <- paste("Could not find file ",vcf.files[i],"{.gz}",sep="")
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
  if(length(vcf.files) != length(sample.files)){
    stop("Number of vcf files and sample files should be equal!")
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
  mysamp <- fread(sample.files[1], header=FALSE, nThread = ncores)
  if(length(sample.files) > 1){
    for(i in 2:length(sample.files)){
      tmp.mysamp <- fread(sample.files[i], header=FALSE, nThread = ncores)
      if(identical(tmp.mysamp,mysamp) == FALSE){
        stop("Sample files differ!")
      }
    }
  }
  mysamp$V2 <- as.character(mysamp$V2)
  fwrite(x = mysamp, file = tmp.samples.file, col.names = FALSE,
         row.names = FALSE, sep = " ", nThread = ncores)
  
  # Convert vcf files
  nids <- nrow(mysamp)
  expcols <- nids + 9
  if(verbose == TRUE){
    cat("\nConverting vcf files to the GHap phase format.\n")
    cat("Note: only bi-allelic variants will be processed.\n")
  }
  for(i in 1:length(vcf.files)){
    head <- scan(file = vcf.files[i], what = "character", sep = "\n", n = 5000, quiet = TRUE)
    head <- grep(pattern = "^#CHROM", x = head)
    if(length(head) != 1){
      emsg <- paste("Could not find unique header in file ",vcf.files[i],"!", sep="")
      stop(emsg)
    }
    myfile <- fread(file = vcf.files[i], header = TRUE, skip = head-1, nThread = ncores)
    if(ncol(myfile) != expcols){
      emsg <- paste("Expected 9 + ",nids," = ",expcols,
                    " columns in file ",vcf.files[i]," but found ", ncol(myfile), "!", sep="")
      stop(emsg)
    }
    expcolnames <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
    obscolnames <- colnames(myfile)[1:9]
    if(identical(obscolnames, expcolnames) == FALSE){
      expcolnames <- c("Expected",expcolnames)
      obscolnames <- c("Observed",obscolnames)
      explen <- max(nchar(expcolnames))
      obslen <- max(nchar(obscolnames))
      explen <- paste("%",explen+1,"s",sep="")
      obslen <- paste("%",obslen+1,"s",sep="")
      emsg <- paste("Unmatching column names in file ",vcf.files[i],":\n\n", sep="")
      emsg <- paste(emsg, paste(sprintf(explen,expcolnames), sprintf(obslen,obscolnames), sep="  ",collapse="\n"))
      stop(emsg)
    }
    if(identical(sort(colnames(myfile)[-c(1:9)]), sort(mysamp$V2)) == FALSE){
      emsg <- paste("Unmatching sample names between files",vcf.files[i],"and",sample.files[i])
      stop(emsg)
    }
    fgt <- grep(x = myfile$FORMAT, pattern = "^GT", invert = TRUE)
    if(length(fgt) > 0){
      emsg <- paste("\nUnsupported FORMAT filed in ",vcf.files[i],"!\n", sep="")
      emsg <- paste(emsg, 'Genotypes should be provided first, as in the example "GT:DS"',sep = "")
      stop(emsg)
    }
    missnpid <- which(myfile$ID == ".")
    if(length(missnpid) > 0){
      newsnpid <- paste(myfile$`#CHROM`[missnpid],myfile$POS[missnpid],sep="_")
      myfile$ID[missnpid] <- newsnpid
    }
    biallelic <- grep(pattern = ",", x = myfile$ALT, invert = TRUE)
    myfile <- myfile[biallelic,]
    fwrite(x = myfile[,c(1,3,2,4,5)], file = tmp.markers.file, nThread = ncores,
           col.names = FALSE, row.names = FALSE, sep = " ", append = TRUE)
    myfile <- myfile[, mysamp$V2, with=FALSE]
    myfile[ , (colnames(myfile)) := lapply(.SD, function(myfile) {gsub("\\:.+", "", myfile)}), .SDcols = colnames(myfile)]
    myfile[ , (colnames(myfile)) := lapply(.SD, function(myfile) {gsub("\\|", " ", myfile)}), .SDcols = colnames(myfile)]
    unphased <- lapply(X = colnames(myfile),
                       FUN = function(selcol){x <- grep("0 0|0 1|1 0|1 1",unlist(myfile[,selcol,with=FALSE]),invert = TRUE); length(x)})
    unphased <- sum(unlist(unphased))
    if(unphased > 0){
      emsg <- paste("Found unexpected characters in ", unphased," genotypes out of ",nrow(myfile),"*",ncol(myfile),
                    " = ",nrow(myfile)*ncol(myfile)," bi-allelic genotypes in file ",vcf.files[i],
                    "!\nAre these genotypes really phased?", sep="")
      stop(emsg)
    }
    fwrite(x = myfile, file = tmp.phase.file, quote = FALSE, nThread = ncores,
           col.names = FALSE, row.names = FALSE, sep = " ", append = TRUE)
    if(verbose == TRUE){
      cat("Processed ",i," files of",length(vcf.files), "\r")
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
