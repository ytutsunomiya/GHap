\name{ghap.phase2plink}
\alias{ghap.phase2plink}

\title{
 Export phase object to PLINK binary
}
\description{
 This function takes a phase object and converts it to the PLINK binary (bed/bim/fam) format.
}
\usage{
ghap.phase2plink(object, out.file,
                 only.active.samples=TRUE,
                 only.active.markers=TRUE,
                 batchsize=NULL, ncores=1,
                 verbose=TRUE)
}

\arguments{
  \item{object}{
    A GHap.phase object.
  }
  \item{out.file}{
   A character value specifying the name used for the .bed, .bim and .fam output files.
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be included in the output (default = TRUE).
  }
  \item{batchsize}{
   A numeric value controlling the number of markers to be processed and written to output at a time (default = nmarkers/10).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}
\details{
 The returned output is a standard set of PLINK (Purcell et al., 2007; Chang et al., 2015) binary file (bed/bim/fam), meaning that phase information will be lost during conversion.
}
\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}
\references{
C. C. Chang et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015. 4, 7.

S. Purcell et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am. J. Hum. Genet. 2007. 81, 559-575.
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
# # Load data
# phase <- ghap.loadphase("example")
# 
# ### RUN ###
# 
# # Convert to plink
# ghap.phase2plink(object = phase, out.file = "example")

}
