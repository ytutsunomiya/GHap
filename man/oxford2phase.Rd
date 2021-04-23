\name{ghap.oxford2phase}
\alias{ghap.oxford2phase}

\title{
  Convert Oxford data into GHap phase
}
\description{
  This function takes phased genotype data in Oxford HAPS/SAMPLES format and converts them into the GHap phase format.
}
\usage{
  ghap.oxford2phase(input.files = NULL, haps.files = NULL, sample.files = NULL,
                    out.file, verbose = TRUE)
}

\arguments{

  If all input files share the same prefix, the user can use the following shortcut options:

  \item{input.files}{
    Character vector with the list of prefixes for input files.
  }
  \item{out.file}{
    Character value for the output file name.
  }
  
  The user can also opt to point to input files separately:
  
  \item{haps.files}{
    Character vector containing the list of Oxford HAPS files.
  }
  \item{sample.files}{
    Character vector containing the list of Oxford SAMPLES files.
  }

  To turn conversion progress-tracking on or off please use:

  \item{verbose}{
  	A logical value specfying whether log messages should be printed (default = TRUE).
  }
}
\details{
The Oxford HAPS/SAMPLE format output of widely used phasing software such as SHAPEIT2 (O'Connell et al., 2014) or Eagle (Loh et al., 2106) is here manipulated to obtain the GHap phase format.
}

\references{
R-R. Loh P-R et al. Reference-based phasing using the Haplotype Reference Consortium panel. Nat Genet. 2016. 48(11):1443-1448.

J. O'Connell et al. A general approach for haplotype phasing across the full spectrum of relatedness. PLOS Genet. 2014. 10:e1004234.
}
\author{
Mario Barbato <mario.barbato@unicatt.it>, Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}
\seealso{
  \code{\link{ghap.compress}}, \code{\link{ghap.loadphase}}, \code{\link{ghap.fast2phase}}, \code{\link{ghap.vcf2phase}}
}
\examples{

# # Build toy example data (2 chromosomes, 5 samples and 5 markers per chromosome)
# samples <- data.frame(V1 = c("ID_1","0","POP1","POP1","POP1","POP2","POP2"),
#                       V2 = c("ID_2","0","ID1","ID2","ID3","ID4","ID5"),
#                       V3 = c("missing","0","0","0","0","0","0"),
#                       stringsAsFactors = FALSE)
# bases <- c("A","T","C","G")
# for(i in 1:2){
#   map <- data.frame(V1 = rep(i, times = 5),
#                     V2 = paste("C",i,"P",1:5,sep=""),
#                     V3 = 1:5,
#                     V4 = sample(x = bases, size = 5, replace = TRUE),
#                     V5 = sample(x = bases, size = 5, replace = TRUE))
#   phase <- matrix(data = sample(x = c(0,1), size = 2*5*5, replace = TRUE), nrow = 5, ncol = 2*5)
#   haps <- cbind(map,phase)
#   write.table(file = paste("example_chr",i,".sample",sep=""), x = samples,
#               quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
#   write.table(file = paste("example_chr",i,".haps",sep=""), x = haps,
#               quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
# }
# 
# # Convert Oxford HAPS/SAMPLES to GHap phase using prefix
# ghap.oxford2phase(input.files = paste("example_chr",1:2,sep=""),
#                   out.file = "example")
# 
# # Convert Oxford HAPS/SAMPLES to GHap phase using file names
# ghap.oxford2phase(haps.files = paste("example_chr",1:2,".haps",sep=""),
#                   sample.files = paste("example_chr",1:2,".sample",sep=""),
#                   out.file = "example")
# 
# # A more efficient alternative for *nix system users
# # Note: replace "cat" by "zcat" if files are gzipped
# haps.files = paste("example_chr",1:2,".haps",sep="")
# command <- "tail -n+3 example_chr1.sample | cut -d' ' -f1,2 > example.samples"
# system(command)
# for(i in 1:2){
#   command <- paste("cat",haps.files[i],"| cut -d' ' -f1-5 >> example.markers")
#   system(command)
#   command <- paste("cat",haps.files[i],"| cut -d' ' -f1-5 --complement >> example.phase")
#   system(command)
# }

}
