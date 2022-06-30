.onAttach <- function(libname, pkgname) {
  mymsg <- "Loading required package: GHap\n\n\n"
  mymsg <- paste(mymsg,"Thanks for using GHap v3.0.0!\n")
  mymsg <- paste(mymsg,"For more information use: help(package = 'GHap')\n")
  mymsg <- paste(mymsg,"                          citation(package = 'GHap')\n")
  mymsg <- paste(mymsg,"                          browseVignettes(package = 'GHap')\n\n")
  mymsg <- paste(mymsg,"Version date: 30 Jun 2022\n\n")
  packageStartupMessage(mymsg)
}
