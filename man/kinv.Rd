\name{ghap.kinv}
\alias{ghap.kinv}

\title{
Inverse of kinship matrix
}
\description{
Inversion of the haplotype covariance matrix
}
\usage{
ghap.kinv(kinship,  method="nearPD", ncores=1, proven=NULL, verbose=TRUE)
}

\arguments{
  \item{kinship}{
   A HapAllele-based kinship matrix, as supplied by \code{\link{ghap.kinship}}.
  }
  \item{method}{
   Inversion method: common inverse ("common"), inverse of the nearest positive definite matrix ("nearPD", default) or approximate inverse based on the Algorithm for Proven and Young animals ("APY").
  }
  \item{ncores}{
  Number of cores to be used for computations (default = 1). Only relevant for method "APY".
 }
  \item{proven}{
   Character vector with the names of proven subjects to be used as reference for method "APY".
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
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
# # Subset
# ids <- haplo$id[which(haplo$pop == "CEU")]
# haplo <- ghap.subsethaplo(haplo = haplo, ids = ids, alleles = haplo$allele.in)
# 
# # Exclude minor alleles and singletons
# hapstats <- ghap.hapstats(haplo, ncores = 1)
# haplo <- ghap.subsethaplo(haplo,ids=ids, alleles = hapstats$TYPE %in% c("REGULAR","MAJOR"))
# 
# # Compute Kinship matrix
# K <- ghap.kinship(haplo, ncores = 1)
#
#
# ### RUN ###
# 
# # Common inverse
# Kinv <- ghap.kinv(K, method="common")
#
# # Inverse of nearest positive definite matrix
# Kinv <- ghap.kinv(K, method="nearPD")
#
# # APY inverse
# Kinv <- ghap.kinv(K, method="APY", proven=colnames(K)[1:500], ncores=1)

}
