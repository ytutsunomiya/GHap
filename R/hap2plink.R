#Function: ghap.hap2plink
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Convert haplotype allele counts to plink bed/bim/fam

ghap.hap2plink <- function(
  object,
  outfile
){
  
  #Check if haplo is a GHap.haplo object
  if(inherits(object, "GHap.haplo") == FALSE){
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if output will ovewrite existing files
  bed <- paste(outfile,"bed",sep=".")
  bim <- paste(outfile,"bim",sep=".")
  fam <- paste(outfile,"fam",sep=".")
  if(file.exists(bed) == TRUE){
    stop("The bed file already exists!")
  }
  if(file.exists(bim) == TRUE){
    stop("The bim file already exists!")
  }
  if(file.exists(fam) == TRUE){
    stop("The fam file already exists!")
  }
  
  #Generate fam file
  famfile <- cbind(object$pop,object$id,"0 0 0 -9")
  fwrite(x = as.data.table(famfile), file = fam,
         quote = FALSE, sep=" ", row.names = FALSE, col.names = FALSE)
  
  #Generate bim file
  bimfile <- data.frame(CHR = object$chr,
                        SNP = paste(object$block, object$bp1, object$bp2, object$allele, sep="_"),
                        CM = (object$bp1+object$bp2)/2e+6,
                        BP = (object$bp1+object$bp2)/2,
                        A1 = "N",
                        A2 = "H")
  bimfile$BP <- as.integer(bimfile$BP)
  fwrite(x = as.data.table(bimfile), file = bim, quote = FALSE, sep=" ",
         row.names = FALSE, col.names = FALSE)
  
  #Generate bed file
  file.copy(from = object$genotypes, to = bed)
  
}
