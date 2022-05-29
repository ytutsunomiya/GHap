\name{ghap.oxford2phase}
\alias{ghap.oxford2phase}

\title{
  Convert Oxford data into GHap phase
}
\description{
  This function takes phased genotype data in Oxford HAPS/SAMPLES format and converts them into the GHap phase format.
}
\usage{
  ghap.oxford2phase(input.files = NULL, haps.files = NULL,
                    sample.files = NULL, out.file,
                    ncores = 1, verbose = TRUE)
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

  To turn conversion progress-tracking on or off or use multiple cores please use:
  
  \item{ncores}{
  	A numeric value specfying the number of cores to use (default = 1).
  }

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

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy the example data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "oxford",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# ### RUN ###
# 
# # Convert from a single genome-wide file
# ghap.oxford2phase(input.files = "example",
#                   out.file = "example")
# 
# # Convert from a list of chromosome files
# ghap.oxford2phase(input.files = paste0("example_chr",1:10),
#                   out.file = "example")
# 
# # Convert using separate lists for file extensions
# ghap.oxford2phase(haps.files = paste0("example_chr",1:10,".haps"),
#                   sample.files = paste0("example_chr",1:10,".sample"),
#                   out.file = "example")
# 
# # A more efficient alternative for *nix system users
# # Note: replace "cat" by "zcat" if files are gzipped
# haps.files = paste("example_chr",1:10,".haps",sep="")
# command <- "tail -n+3 example_chr1.sample | cut -d' ' -f1,2 > example.samples"
# system(command)
# for(i in 1:10){
#   command <- paste("cat",haps.files[i],"| cut -d' ' -f1-5 >> example.markers")
#   system(command)
#   command <- paste("cat",haps.files[i],"| cut -d' ' -f1-5 --complement >> example.phase")
#   system(command)
# }

}
