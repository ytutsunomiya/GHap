\name{ghap.freq}
\alias{ghap.freq}

\title{
  Compute marker allele frequencies
}
\description{
  This function takes a GHap.phase object and computes the allele frequency for each marker.
}
\usage{
  ghap.freq(phase, type = "maf", only.active.samples = TRUE, only.active.markers = TRUE,
           batchsize = NULL, ncores = 1, verbose = TRUE)
}

\arguments{
  \item{phase}{
    A GHap.phase object.
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
  \item{batchsize}{
   A numeric value controlling the number of markers to be processed at a time (default = nmarkers/10).
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
# # Copy the example data in the current working directory
# exfiles <- ghap.makefile()
# file.copy(from = exfiles, to = "./")
# 
# # Compress phase data
# ghap.compress(input.file = "human", out.file = "human")
# 
# # Load phase data
# phase <- ghap.loadphase(input.file = "human")
#
#
# ### RUN ###
# 
# # Calculate minor allele frequency
# maf <- ghap.freq(phase, type = 'maf', ncores = 1)

}