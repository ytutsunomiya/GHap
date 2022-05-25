#Function: ghap.makefile
#License: GPLv3 or later
#Modification date: 25 May 2022
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Create a copy of an example file in a temporary directory

ghap.makefile <- function(
  dataset = "example",
  format = "phase",
  verbose = TRUE
){
  
  # Check if format is supported -----------------------------------------------
  extension <- NULL
  extension$phase <- c(".markers",".samples",".phaseb")
  extension$plink <- c(".bed",".bim",".fam")
  extension$haplo <- c(".hapgenotypesb",".hapsamples",".hapalleles")
  extension$raw <- c(".markers",".samples",".phase")
  extension$vcf <- c(".vcf",".sample")
  extension$oxford <- c(".haps",".sample")
  extension$meta <- c(".phenotypes",".pedigree")
  if(format %in% names(extension) == FALSE){
    emsg <- paste("\nFormat", format, "is not supported")
    stop(emsg)
  }
  
  # Check if requested dataset exists ------------------------------------------
  # Code borrowed from stackoverflow
  urlfile <- paste0("https://github.com/ytutsunomiya/GHap/blob/main/datasets/",
                    dataset, "_", format, ".zip?raw=true")
  con <- url(urlfile)
  check <- suppressWarnings(try(open.connection(con,open="rt",timeout=2),silent=T)[1])
  suppressWarnings(try(close.connection(con),silent=T))
  check <- ifelse(is.null(check),TRUE,FALSE)
  if(check == FALSE){
    emsg <- paste("\nCould not find dataset", dataset, "with format", format)
    stop(emsg)
  }
  
  # Download files to temporary directory ---------------------------------------
  tmpdir <- tempdir()
  tmpfile <- tempfile(dataset, fileext = ".zip")
  if(Sys.info()["sysname"] == "Windows"){
    cpfile <- download.file(url = urlfile, destfile = tmpfile, quiet = TRUE, mode = "wb")
  }else{
    cpfile <- download.file(url = urlfile, destfile = tmpfile, quiet = TRUE)
  }
  unzip(zipfile = tmpfile, exdir = tmpdir)
  rfile <- file.remove(tmpfile)
  
  # Check whether files have been copied and extracted ---------------------------
  patterns <- paste(dataset,".+\\", unlist(extension[format]), sep="", collapse="|")
  patterns <- c(patterns, paste(dataset, "\\", unlist(extension[format]), sep="", collapse="|"))
  patterns <- paste(patterns, collapse = "|")
  outfiles <- paste0(tmpdir, "/", list.files(path = tmpdir, pattern = patterns))
  if(cpfile == 0 & verbose == TRUE){
    cat("\nFiles successfully created!\n")
    cat("\nImportant NOTE:")
    cat("\nIn compliance to CRAN policies, the files have been downloaded to:\n\n")
    cat(paste(outfiles, sep="", collapse="\n"))
    cat("\n\nThese files will be deleted when the R session is finished.\n\n")
  }else{
    stop("Error when copying the files!\nThis function needs super-user privilege and a working untar engine.")
  }
  
  # Return file names ------------------------------------------------------------
  return(outfiles)
  
}
