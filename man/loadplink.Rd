\name{ghap.loadplink}
\alias{ghap.loadplink}

\title{
  Load binary PLINK data
}
\description{
  This function loads binary PLINK files (bed/bim/fam) and converts them into a native GHap.plink object.
}
\usage{
  ghap.loadplink(input.file = NULL, bed.file = NULL,
                 bim.file = NULL, fam.file = NULL,
                 ncores = 1, verbose = TRUE)
}

\arguments{

  If all input files share the same prefix, the user can use the following shortcut option:

  \item{input.file}{
    Prefix for input files.
  }
  
  For backward compatibility, the user can still point to input files separately:
  
  \item{bed.file}{
    The binary genotype matrix (in SNP-major format).
  }
  \item{bim.file}{
    Variant map file.
  }
  \item{fam.file}{
    Pedigree (family) file.
  }
  
  To turn loading progress-tracking on or off, or engage multiple cores, please use:
  
  \item{ncores}{
    A numerical value specfying the number of cores to use while loading the input files (default = 1).
  }
  
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
  The returned GHap.plink object is a list with components:
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
    A character vector relating genotypes to populations. This information is obtained from the FID (1st) column in the fam file.
  }
  \item{id}{
    A character vector mapping genotypes to samples. This information is obtained from the IID (2nd) column in the fam file.
  }
  \item{id.in}{
    A logical vector indicating active chromosome alleles. By default, all chromosomes are set to TRUE.
  }
  \item{sire}{
      A character vector indicating sire names, as provided in the SID (3rd) column of the fam file.
  }
  \item{dam}{
      A character vector indicating dam names, as provided in the DID (4th) column of the fam file.
  }
  \item{sex}{
      A character vector indicating individual sex, as provided in the SEX (5th) column of the fam file. Codes are converted as follows: 0 = NA, 1 = Male and 2 = Female.
  }
  \item{chr}{
      A character vector indicating chromosome identity for each marker.
  }
  \item{marker}{
    A character vector containing marker names.
  }
  \item{marker.in}{
    A logical vector indicating active markers. By default, all markers are set to TRUE.
  }
  \item{cm}{
    A numeric vector with genetic positions for markers. This information is obtained from the third column of the bim file. If genetic positions are absent (coded as "0"), they are approximated from physical positions assuming 1 Mb ~ 1 cM.
  }
  \item{bp}{
    A numeric vector with physical positions for markers.
  }
  \item{A0}{
    A character vector with reference alleles. For convenience, this information is obtained from the 6th column of the bim file. If "--keep-allele-order" is not used while generating the PLINK binary file, A0 will correspond to the major allele.
  }
  \item{A1}{
    A character vector with alternative alleles. As for A0, if "--keep-allele-order" is not used A1 will correspond to the minor allele.
  }
  \item{plink}{
    A character value giving the pathway to the binary genotype matrix.
  }
}

\author{
  Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy phase data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "plink",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# ### RUN ###
# 
# # Load data using prefix
# plink <- ghap.loadplink("example")
# 
# # Load data using file names
# plink <- ghap.loadplink(bed.file = "example.bed",
#                         bim.file = "example.bim",
#                         fam.file = "example.fam")


}
