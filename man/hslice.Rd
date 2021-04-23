\name{ghap.hslice}
\alias{ghap.hslice}

\title{
  Get a slice of the haplo object
}
\description{
  This function parses the binary HapGenotypes matrix and returns the slice as an R matrix.
}
\usage{
ghap.hslice(haplo, ids, alleles, index=FALSE, lookup=NULL,
            ncores=1, verbose=TRUE)
}

\arguments{
  \item{haplo}{
    A GHap.haplo object.
  }
  \item{ids}{
    A character or numeric vector indicating individuals to parse.
  }
  \item{alleles}{
    A numeric vector indicating alleles to parse.
  }
  \item{index}{
    A logical value specfying if values provided for ids are indices (see details).
  }
  \item{lookup}{
    A character vector containing the look up table to decode the binary file. Only meaningful for developers.
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
  An R matrix with haplotype alleles in rows and individuals in columns.
}

\details{
  This function parses the file with suffix \emph{.hapgenotypesb} and returns an R matrix with the requested list of haplotype alleles and individuals. The argument index allows the user to specify individuals either by name (index = FALSE) or by indices as stored in the GHap.haplo object (index = TRUE). HapAlleles can only by parsed via indices. The lookup argument also allows the user to provide a vector containing the decoding scheme to translate bits into integers in R. However, this argument should only be used by developers.
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
# # Compress phase data
# ghap.compress(input.file = "human", out.file = "human")
# 
# # Load phase data
# phase <- ghap.loadphase(input.file = "human")
# 
# # Generate blocks of 5 markers sliding 5 markers at a time
# blocks <- ghap.blockgen(phase, windowsize = 5, slide = 5, unit = "marker")
# 
# # Haplotyping
# ghap.haplotyping(phase = phase, blocks = blocks, outfile = "human",
#                  binary = T, ncores = 1)
#
# # Load haplotype genotypes
# haplo <- ghap.loadhaplo(input.file = "human")
# 
# ### RUN ###
#
# ids <- sample(x = haplo$id, size = 10, replace = F)
# alleles <- sample(x = 1:haplo$nalleles, size = 10, replace = F)
# X <- ghap.hslice(haplo = haplo, ids = ids, alleles = alleles, index = F, ncores = 1)
  
}
