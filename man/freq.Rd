\name{ghap.freq}
\alias{ghap.freq}

\title{
  Compute marker allele frequencies
}
\description{
  This function takes a GHap.phase object and computes the allele frequency for each marker.
}
\usage{
  ghap.freq(object, type = "maf",
            only.active.samples = TRUE,
            only.active.markers = TRUE,
            ncores = 1, verbose = TRUE)
}

\arguments{
  \item{object}{
    A GHap object of type phase or plink.
  }
  \item{type}{
    A character value indicating which allele frequency to compute. Valid options are minor allele frequency ('maf', default), frequency of allele 0 ('A0') and frequency of allele 1 ('A1').
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be used for calculations (default = TRUE).
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be included in the output (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
  The function outputs a numeric vector of the same length of active markers containing allele frequencies based on the active samples.
}

\author{
  Yuri Tani Utsunomiya <ytutsunomiya@gmail.com> \cr
  Marco Milanesi <marco.milanesi.mm@gmail.com>
}

\examples{
  
# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy phase data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "phase",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Copy plink data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "plink",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# ### RUN ###
# 
# # Calculate allele frequency for phase data
# phase <- ghap.loadphase("example")
# q <- ghap.freq(phase, type = 'A0')
# p <- ghap.freq(phase, type = 'A1')
# maf <- ghap.freq(phase, type = 'maf')
# 
# # Calculate allele frequency for plink data
# plink <- ghap.loadplink("example")
# q <- ghap.freq(plink, type = 'A0')
# p <- ghap.freq(plink, type = 'A1')
# maf <- ghap.freq(plink, type = 'maf')

}
