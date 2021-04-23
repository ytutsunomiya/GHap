#Function: ghap.makefile
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Create a copy of the example file in a temporary directory

ghap.makefile<-function(verbose = TRUE){
  
  #Download files to temporary directory
  tmpdir <- tempdir()
  tmpfile <- tempfile("HapMap3_chr2", fileext = ".zip")
  urlfile <- "https://bitbucket.org/marcomilanesi/ghap/raw/92c0e2bd4370d9531048bb47b4f57d520be6501e/HapMap3_chr2.zip"
  cpfile <- download.file(url = urlfile, destfile = tmpfile, quiet = TRUE)
  unzip(zipfile = tmpfile, exdir = tmpdir)
  rfile <- file.remove(tmpfile)
  
  #Check whether files have been copied and extracted
  outfiles <- paste(tmpdir,"/human",c(".markers",".samples",".phase"),sep="")
  exfile <- sum(file.exists(outfiles))
  if(cpfile == 0 & exfile == 3 & verbose == TRUE){
    cat("\nFiles successfully created!\n")
    cat("\nImportant NOTE:")
    cat("\nIn compliance to CRAN policies, the files have been downloaded to:\n\n")
    cat(paste(outfiles, sep="", collapse="\n"))
    cat("\n\nThese files will be deleted when the R session is finished.\n\n")
  }else{
    stop("Error when copying the files!\nThis function needs super-user privilege and a working untar engine.")
  }
  
  #Return file names
  return(outfiles)
  
}