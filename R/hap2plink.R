#Function: ghap.hap2plink
#License: GPLv3 or later
#Modification date: 11 Sep 2020
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Convert haplotype allele counts to plink bed/bim/fam

ghap.hap2plink <- function(
  haplo,
  outfile
){
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
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
  famfile <- cbind(haplo$pop,haplo$id,"0 0 0 -9")
  fwrite(x = as.data.table(famfile), file = fam,
         quote = FALSE, sep=" ", row.names = FALSE, col.names = FALSE)
  
  #Generate bim file
  bimfile <- data.frame(CHR = haplo$chr,
                        SNP = paste(haplo$block, haplo$bp1, haplo$bp2, haplo$allele, sep="_"),
                        CM = (haplo$bp1+haplo$bp2)/2e+6,
                        BP = (haplo$bp1+haplo$bp2)/2,
                        A1 = "N",
                        A2 = "H")
  bimfile$BP <- as.integer(bimfile$BP)
  fwrite(x = as.data.table(bimfile), file = bim, quote = FALSE, sep=" ",
         row.names = FALSE, col.names = FALSE)
  
  #Generate bed file
  file.copy(from = haplo$genotypes, to = bed)
  
}
