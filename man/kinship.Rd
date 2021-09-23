\name{ghap.kinship}
\alias{ghap.kinship}

\title{
Relationship matrix based on genomic data
}
\description{
This function computes marker-based and HapAllele-based relationship matrices.
}
\usage{
ghap.kinship(object, weights = NULL, sparsity = NULL,
             batchsize = NULL, only.active.samples = TRUE,
             only.active.variants = TRUE, ncores = 1,
             verbose = TRUE)
}

\arguments{
  \item{object}{
    A valid GHap object (phase, haplo or plink).
  }
  \item{weights}{
    A numeric vector providing variant-specific weights.
  }
  \item{sparsity}{
    A numeric value specifying a relationship cut-off (default = NULL). All relationships below the specified cut-off will be set to zero, inducing sparsity into the relationship matrix.
  }
  \item{batchsize}{
   A numeric value controlling the number of variants to be processed at a time (default = nalleles/10).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.alleles}{
    A logical value specifying whether only active variants should be included in the output (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
The function returns a \emph{n} x \emph{n} relationship matrix, where \emph{n} is the number of individuals.
}

\details{

Let \eqn{\mathbf{M}} be the centered \emph{n} x \emph{m} matrix of genotypes, where \emph{n} is the number of individuals and \emph{m} is the number of variants (i.e, markers or HapAlleles). The relationship matrix is computed as:

\deqn{\mathbf{K} = q\mathbf{MDM}'}
where \eqn{\mathbf{D} = diag(d_i)}, \eqn{d_i} is the weight of variant \emph{i} (default \eqn{d_i = 1}), and \eqn{q} is a scaling factor defined as \eqn{tr(\mathbf{MDM}')^{-1}m}. This is a generalization of the SNP-based genomic relationship matrix (VanRaden, 2008).

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
# # Copy phase data to the current working directory
# exfiles <- ghap.makefile()
# file.copy(from = exfiles, to = "./")
# 
# # Copy plink data to the current working directory
# exfiles <- ghap.makefile(format = "plink")
# file.copy(from = exfiles, to = "./")
# 
# # Load phase data
# phase <- ghap.loadphase(input.file = "human_hapmap3_chr2")
# 
# # Load plink data
# plink <- ghap.loadplink(input.file = "human_hapmap3_chr2")
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
# # Exclude minor alleles and singletons
# hapstats <- ghap.hapstats(haplo, ncores = 1)
# haplo <- ghap.subsethaplo(haplo,ids=haplo$id,
#                           alleles = hapstats$TYPE \%in\% c("REGULAR","MAJOR"))
# 
# 
# ### RUN ###
# 
# # Relationship matrix with phase object
# K1 <- ghap.kinship(phase, ncores = 1)
# 
# # Relationship matrix with plink object
# K2 <- ghap.kinship(plink, ncores = 1)
# 
# # Relationship matrix with haplo object
# K3 <- ghap.kinship(haplo, ncores = 1)

}
