#Function: ghap.exfiles
#License: GPLv3 or later
#Modification date: 13 May 2021
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: list example files available

ghap.exfiles <- function(){
  
  # Get file list
  out <- fread("https://raw.githubusercontent.com/ytutsunomiya/GHap/main/datasets/exlist.txt", showProgress = FALSE)
 
  # Return file names
  return(out)
  
}
