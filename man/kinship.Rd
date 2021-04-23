\name{ghap.kinship}
\alias{ghap.kinship}

\title{
Kinship matrix from haplotypes
}
\description{
This function computes a HapAllele-based kinship matrix from a GHap.haplo object.
}
\usage{
ghap.kinship(haplo, weights=NULL, batchsize=NULL, only.active.samples=TRUE,
             only.active.alleles=TRUE, ncores=1, verbose=TRUE)
}

\arguments{
  \item{haplo}{
    A GHap.haplo object.
  }
  \item{weights}{
    A numeric vector providing HapAllele-specific weights.
  }
  \item{batchsize}{
   A numeric value controlling the number of HapAlleles to be processed at a time (default = nalleles/10).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.alleles}{
    A logical value specifying whether only active haplotype alleles should be included in the output (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
The function returns a \emph{n} x \emph{n} matrix of HapAllele-based kinships, where \emph{n} is the number of individuals.
}

\details{

Let \eqn{\mathbf{M}} be the centered \emph{N} x \emph{H} matrix of HapGenotypes, where \emph{N} is the number of individuals and \emph{H} is the number of HapAlleles. The HapAllele covariance among individuals is computed as:

\deqn{\mathbf{K} = q\mathbf{MDM}'}

where \eqn{\mathbf{D} = diag(d_i)}, \eqn{d_i} is the weight of HapAllele \emph{i} (default \eqn{d_i = 1}), and \eqn{q} is a scaling factor defined as \eqn{tr(\mathbf{MDM}')^{-1}M}. This is a generalization of the SNP-based genomic relationship matrix (VanRaden, 2008).

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com> \cr
Marco Milanesi <marco.milanesi.mm@gmail.com>
}


\references{
P. M. VanRaden. Efficient methods to compute genomic predictions. J. Dairy. Sci. 2008. 91:4414-4423.
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
#
# ### RUN ###
# 
# # Exclude minor alleles and singletons
# hapstats <- ghap.hapstats(haplo, ncores = 1)
# haplo <- ghap.subsethaplo(haplo,ids=haplo$id, alleles = hapstats$TYPE \%in\% c("REGULAR","MAJOR"))
#
# # Compute Kinship matrix
# K <- ghap.kinship(haplo, ncores = 1)

}
