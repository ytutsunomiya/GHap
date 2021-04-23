\name{ghap.loadphase}
\alias{ghap.loadphase}

\title{
  Load binary phased genotype data
}
\description{
  This function loads binary phased genotype data and converts them into a native GHap.phase object.
}
\usage{
  ghap.loadphase(input.file = NULL, samples.file = NULL, markers.file = NULL,
                 phaseb.file = NULL, verbose = TRUE)
}

\arguments{

  If all input files share the same prefix, the user can use the following shortcut option:

  \item{input.file}{
    Prefix for input files.
  }
  
  For backward compatibility, the user can still point to input files separately:
  
  \item{samples.file}{
    Individual information.
  }
  \item{markers.file}{
    Variant map information.
  }
  \item{phaseb.file}{
    Binary phased genotype matrix, such as supplied by the \code{\link{ghap.compress}} function.
  }
  
  To turn loading progress-tracking on or off please use:
  
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
  The returned GHap.phase object is a list with components:
    \item{chr}{
      A character vector indicating chromosome identity for each marker.
    }
  \item{nsamples}{
    An integer value for the sample size.
  }
  \item{nmarkers}{
    An integer value for the number of markers.
  }
  \item{nsamples.in}{
    An integer value for the number of active samples.
  }
  \item{nmarkers.in}{
    An integer value for the number of active markers.
  }
  \item{pop}{
    A character vector relating chromosome alleles to populations. This information is obtained from the first column of the sample file.
  }
  \item{id}{
    A character vector mapping chromosome alleles to samples. This information is obtained from the second column of the sample file.
  }
  \item{id.in}{
    A logical vector indicating active chromosome alleles. By default, all chromosomes are set to TRUE.
  }
  \item{marker}{
    A character vector containing marker names. This information is obtained from the second column of the marker map file.
  }
  \item{marker.in}{
    A logical vector indicating active markers. By default, all markers are set to TRUE.
  }
  \item{bp}{
    A numeric vector with marker positions. This information is obtained from the third column of the marker map file.
  }
  \item{A0}{
    A character vector with reference alleles. This information is obtained from the fourth column of the marker map file.
  }
  \item{A1}{
    A character vector with alternative alleles. This information is obtained from the fifth column of the marker map file.
  }
  \item{phase}{
    A character value giving the pathway to the binary phased genotype matrix.
  }
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
# ### RUN ###
# 
# # Load data using prefix
# phase <- ghap.loadphase(input.file = "human")
#
# # Load data using file names
# phase <- ghap.loadphase(samples.file = "human.samples",
#                         markers.file = "human.markers",
#                         phaseb.file = "human.phaseb")

}