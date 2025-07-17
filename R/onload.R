.onAttach <- function(libname, pkgname) {
  mymsg <- "Loading required package: GHap\n\n\n"
  mymsg <- paste(mymsg,"Thanks for using GHap v4.0.0.13 (beta/dev)!\n")
  mymsg <- paste(mymsg,"For more information use: help(package = 'GHap')\n")
  mymsg <- paste(mymsg,"                          citation(package = 'GHap')\n")
  mymsg <- paste(mymsg,"                          browseVignettes(package = 'GHap')\n\n")
  mymsg <- paste(mymsg,"Version date: 17 Jul 2025\n\n")
  packageStartupMessage(mymsg)
}
