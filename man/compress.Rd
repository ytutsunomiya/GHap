\name{ghap.compress}
\alias{ghap.compress}

\title{
  Compress phased genotype data
}
\description{
  This function takes phased genotype data and converts them into a compressed binary format.
}
\usage{
  ghap.compress(input.file = NULL, out.file,
                samples.file = NULL, markers.file = NULL,
                phase.file = NULL, batchsize = NULL,
                ncores = 1, verbose = TRUE)
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
  
  To turn compression progress-tracking on or off, or to control parallelization of the task please use:
  
  \item{batchsize}{
   A numeric value controlling the number of markers to be compressed and written to output at a time (default = nmarkers/10).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\details{
  The supported input format is composed of three files with suffix:
    \itemize{
      \item \strong{.samples}: space-delimited file without header containing two mandatory columns: Population and ID. Please notice that the Population column serves solely for the purpose of grouping samples, so the user can define any arbitrary family/cluster/subgroup and use as a "population" tag. This file may further contain three additional columns, which are optional: Sire, Dam and Sex (with code 1 = M and 2 = F). Values "0" and "NA" in these additional columns are treated as missing values. 
      \item \strong{.markers}: space-delimited file without header containing five mandatory columns: Chromosome, Marker, Position (in bp), Reference Allele (A0) and Alternative Allele (A1). Markers should be sorted by chromosome and position. Repeated positions are tolerated, but the user is warned of their presence in the data. Optionally, the user may provide a file containing an additional column with genetic positions (in cM), which has to be placed between the base pair position and the reference allele columns.
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
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "raw",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# ### RUN ###
# 
# # Compress phase data using prefix
# ghap.compress(input.file = "example",
#               out.file = "example")
# 
# # Compress phase data using file names
# ghap.compress(samples.file = "example.samples",
#               markers.file = "example.markers",
#               phase.file = "example.phase",
#               out.file = "example")

}
