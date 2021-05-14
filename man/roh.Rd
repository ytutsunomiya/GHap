\name{ghap.roh}
\alias{ghap.roh}

\title{
 Detection of runs of homozygosity (ROH)
}
\description{
  Map haplotype segments that are likely identical-by-descent.
}
\usage{
 ghap.roh(object, minroh = 1e+6, method = "hmm", freq = NULL,
          inbcoef = NULL, error = 0.25/100, only.active.samples = TRUE,
          only.active.markers = TRUE, ncores = 1, verbose = TRUE)
}

\arguments{

  The following arguments are used by both the 'hmm' and 'naive' methods:

  \item{object}{
   A valid GHap object (phase or plink).
  }
  \item{minroh}{
   Minimum ROH length to output.
  } 
  \item{method}{
   Character value indicating which method to use: 'naive' or 'hmm' (default).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the search (default = TRUE).
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be used in the search (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
  
  The following arguments are only used by the 'hmm' method:
  
  \item{freq}{
   Named numeric vector of allele frequencies, such as provided by function \code{\link{ghap.freq}}.
  }
  \item{inbcoef}{
   Named numeric vector of starting values for genomic inbreeding (i.e., guess for the proportion of the genome covered by ROH).
  }
  \item{error}{
   Numeric value representing the expected genotyping error rate (default = 0.001).
  }


}

\value{
The function returns a dataframe with the following columns:

\item{POP}{
  Original population label.
}

\item{ID}{
  Individual name.
}

\item{CHR}{
  Chromosome name.
}

\item{BP1}{
  Segment start position.
}

\item{BP2}{
  Segment end position.
}

\item{LENGTH}{
  Length of run of homozygosity.
}


}

\details{
This function searchs for runs of homozygosity (ROH) via two different methods: \cr

The 'naive' method simply finds streches of homozygous genotypes in the observed haplotypes that are larger then a user-definied minimum size (default is 1 Mbp). The 'hmm' method uses a Hidden Markov Model that takes genotyping error and recombination into account while detecting ROHs. The 'hmm' model in GHap is similar to the ones described by Narasimhan et al (2016) and Druet & Gautier (2017), differing slightly in model fitting and definition of transition and emission probabilities (details are covered in our vignette). \cr

The 'hmm' method requires allele frequencies for each marker, as well as starting values for the expected proportion of the genome covered by ROH (genomic inbreeding) for each individual. Estimates of allele frequencies can be either based on a reference or estimated from the data with the \code{\link{ghap.freq}} function. Starting values for genomic inbreeding can be obtained by running the function with the 'naive' method first and then computing starting values with \code{\link{ghap.froh}} (see the examples).\cr
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
V. Narasimhan et al. BCFtools/RoH: a hidden Markov model approach for detecting autozygosity from next-generation sequencing data. Bioinformatics. 2016. 32:1749-1751.

T. Druet & M. Gautier. A model-based approach to characterize individual inbreeding at both global and local genomic scales. Molecular Ecology. 2017. 26:5820-5841.

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
# ### RUN ###
# 
# # ROH via the 'naive' method
# roh <- ghap.roh(phase, method = "naive", ncores = 1)
# froh <- ghap.froh(phase, roh, ncores = 1)
# 
# # ROH via the 'hmm' method
# freq <- ghap.freq(phase, ncores = 1)
# inbcoef <- froh$FROH1; names(inbcoef) <- froh$ID
# roh <- ghap.roh(phase, method = "hmm", freq = freq,
#                 inbcoef = inbcoef, ncores = 1)
# froh <- ghap.froh(phase, roh, ncores = 1)
}

\seealso{
 \code{\link{ghap.freq}}, \code{\link{ghap.froh}}
}
