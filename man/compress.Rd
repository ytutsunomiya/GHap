\name{ghap.compress}
\alias{ghap.compress}

\title{
  Compress phased genotype data
}
\description{
  This function takes phased genotype data and converts them into a compressed binary format.
}
\usage{
  ghap.compress(input.file = NULL, out.file, samples.file = NULL,
                markers.file = NULL, phase.file = NULL, verbose = TRUE)
}

\arguments{

  If all input files share the same prefix, the user can use the following shortcut options:

  \item{input.file}{
    Prefix for input files.
  }
  \item{out.file}{
    Output file name.
  }
  
  For backward compatibility, the user can still point to input files separately:
  
  \item{samples.file}{
    Individual information.
  }
  \item{markers.file}{
    Variant map information.
  }
  \item{phase.file}{
    Phased genotype matrix.
  }
  
  To turn compression progress-tracking on or off please use:
  
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\details{
  The supported input format is composed of three files with suffix:
    \itemize{
      \item \strong{.samples}: space-delimited file without header containing two columns: Population and ID. Please notice that the Population column serves solely for the purpose of grouping samples, so the user can define any arbitrary family/cluster/subgroup and use as a "population" tag. 
      \item \strong{.markers}: space-delimited file without header containing five columns: Chromosome, Marker, Position (in bp), Reference Allele (A0) and Alternative Allele (A1). Markers should be sorted by chromosome and position. Repeated positions are not tolerated.
      \item \strong{.phase}: space-delimited file without header containing the phased genotype matrix. The dimension of the matrix is expected to be \emph{m x 2n}, where \emph{m} is the number of markers and \emph{n} is the number of individuals (i.e., two columns per individual, representing the two phased chromosome alleles). Alleles must be coded as 0 or 1. No missing values are allowed, since imputation is assumed to be part of the phasing procedure.
    }
  The function outputs a binary file with suffix \strong{.phaseb}. Each allele is stored as a bit in that file. Bits for any given marker are arranged in a sequence of bytes. Since each marker requires storage of 2*nsamples bits, the number of bytes consumed by a single marker in the output file is ceiling(2*nsamples). If the number of alleles is not a multiple of 8, bits in the remainder of the last byte are filled with 0. All functions in GHap were carefully designed to decode the bytes of a marker in such a way that trailing bits are ignored if present.
}

\author{
  Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy the example data in the current working directory
# exfiles <- ghap.makefile()
# file.copy(from = exfiles, to = "./")
# 
# ### RUN ###
# 
# # Compress phase data using prefix
# ghap.compress(input.file = "human", out.file = "human")
#
# # Compress phase data using file names
# ghap.compress(samples.file = "human.samples",
#               markers.file = "human.markers",
#               phase.file = "human.phase",
#               out.file = "human")

}