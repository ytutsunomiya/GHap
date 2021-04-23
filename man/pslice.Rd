\name{ghap.pslice}
\alias{ghap.pslice}

\title{
 Get a slice of the phase object
}
\description{
  This function parses the binary phased genotype matrix and returns the slice as an R matrix.
}
\usage{
 ghap.pslice(phase, ids, markers, index=FALSE, unphase=FALSE, 
             lookup=NULL, ncores=1, verbose=TRUE)
}

\arguments{
  \item{phase}{
   A GHap.phase object.
  }
  \item{ids}{
   A character or numeric vector indicating individuals to parse.
  }
  \item{markers}{
   A character or numeric vector indicating markers to parse.
  }
  \item{index}{
   A logical value specfying if values provided for ids and markers are indices (see details).
  }
  \item{unphase}{
   A logical value specfying if phased genotypes should be retrieved as unphased allele counts.
  }
  \item{lookup}{
   A character vector containing the look up table to decode the binary phase file. Only meaningful for developers.
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
 An R matrix with markers in rows and haplotypes in columns.
}

\details{
 This function parses the file with suffix \emph{.phaseb} and returns an R matrix with the requested list of markers and individuals. The argument index allows the user to specify markers and individuals either by name (index = FALSE) or by indices as stored in the GHap.phase object (index = TRUE). The lookup argument also allows the user to provide a vector containing the decoding scheme to translate bits into integers in R. However, this argument should only be used by developers.
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
# ### RUN ###
# 
# # Get 1000 random markers for CEU individuals
# ids <- unique(phase$id[which(phase$pop == "CEU")])
# markers <- sample(phase$marker, size = 1000, replace = FALSE)
# X <- ghap.pslice(phase = phase, ids = ids, markers = markers, index = FALSE)
#
# # Get the same data as before but unphased
# X <- ghap.pslice(phase = phase, ids = ids, markers = markers,
#                  index = FALSE, unphase = TRUE)

}
