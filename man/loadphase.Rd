\name{ghap.loadphase}
\alias{ghap.loadphase}

\title{
  Load binary phased genotype data
}
\description{
  This function loads binary phased genotype data and converts them into a native GHap.phase object.
}
\usage{
  ghap.loadphase(input.file = NULL,
                 samples.file = NULL,
                 markers.file = NULL,
                 phaseb.file = NULL,
                 ncores = 1, verbose = TRUE)
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
  
  To turn loading progress-tracking on or off, or use multiple cores, please use:
  \item{ncores}{
    A numerical value specfying the number of cores to use while loading the input files (default = 1).
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
  The returned GHap.phase object is a list with components:
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
  \item{sire}{
    A character vector indicating sire names, as provided in the third column of the sample file (optional).
  }
  \item{dam}{
    A character vector indicating dam names, as provided in the fourth column of the sample file (optional).
  }
  \item{sex}{
    A character vector indicating individual sex, as provided in the fifth column of the sample file (optional). Codes are converted as follows: 0 = NA, 1 = Male and 2 = Female.
  }
  \item{chr}{
    A character vector indicating chromosome identity for each marker.
  }
  \item{marker}{
    A character vector containing marker names. This information is obtained from the second column of the marker map file.
  }
  \item{marker.in}{
    A logical vector indicating active markers. By default, all markers are set to TRUE.
  }
  \item{cm}{
    A numeric vector with genetic positions for markers. This information is obtained from the third column of the marker map file if it contains 6 columns. Otherwise, if the map file contains only 5 columns, genetic positions are considered absent and approximated from physical positions (in this case assumed to be the third column) as 1 Mb ~ 1 cM.
  }
  \item{bp}{
    A numeric vector with marker positions. This information is obtained from the third column of the marker map file if it contains 5 columns, or from the fourth column if it contains 6 columns.
  }
  \item{A0}{
    A character vector with reference alleles. This information is obtained from the fourth column of the marker map file in case it contains 5 columns, or from the fifth column if it contains 6 columns.
  }
  \item{A1}{
    A character vector with alternative alleles. This information is obtained from the last column of the marker map file.
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
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "phase",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# ### RUN ###
# 
# # Load data using prefix
# phase <- ghap.loadphase(input.file = "example")
# 
# # Load data using file names
# phase <- ghap.loadphase(samples.file = "example.samples",
#                         markers.file = "example.markers",
#                         phaseb.file = "example.phaseb")

}
