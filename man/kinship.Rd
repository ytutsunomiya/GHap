\name{ghap.kinship}
\alias{ghap.kinship}

\title{
Relationship matrix based on genomic data
}
\description{
This function computes marker-based and HapAllele-based relationship matrices.
}
\usage{
ghap.kinship(object, weights = NULL,
             sparsity = NULL, type = 1,
             batchsize = NULL, freq = NULL,
             idrow = NULL, idcol = NULL,
             only.active.samples = TRUE,
             only.active.variants = TRUE,
             ncores = 1, verbose = TRUE)
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
  \item{type}{
    A numeric value indicating the type of relationship matrix (see details).
  }
  \item{batchsize}{
   A numeric value controlling the number of variants to be processed at a time (default = nalleles/10).
  }
  \item{freq}{
    A named numeric vector with (A1) allele frequencies computed in a reference sample. If not provided, frequencies will be computed using the active individuals.
  }
  \item{idrow}{
    A character vector providing ids of individuals to be included in rows (default = NULL, see details).
  }
  \item{idcol}{
    A character vector providing ids of individuals to be included in columns (default = NULL, see details).
  }  
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.variants}{
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

Let \eqn{\mathbf{M}} be the \emph{n} x \emph{m} matrix of genotypes, where \emph{n} is the number of individuals and \emph{m} is the number of variants (i.e, markers or HapAlleles). Prior to computation, genotypes in matrix \eqn{\mathbf{M}} are transformed according to the desired relationship type. After that transformation, the relationship matrix is computed as:

\deqn{\mathbf{K} = q^{-1}\mathbf{MDM}'}
where \eqn{\mathbf{D} = diag(d_i)}, \eqn{d_i} is the weight of variant \emph{i} (default \eqn{d_i = 1}), and \eqn{q} is a scaling factor. The argument type controls the genotype transformation and the scaling factor, and includes the following option:
\cr
Type = 1 (General additive 1)\cr
Genotype transformation: m - mean(m)\cr
Scaling factor: \eqn{tr(\mathbf{MDM}')^{-1}n}

Type = 2 (General additive 2)\cr
Genotype transformation: (m - mean(m))/sd(m)\cr
Scaling factor: m

Type = 3 (VanRaden, 2008)\cr
Genotype transformation: m - 2*p[j]\cr
Scaling factor: 2*sum(p*(1-p))

Type = 4 (GCTA)\cr
Genotype transformation: (m - 2*p[j])/sqrt(2*p[j]*(1-p[j]))\cr
Scaling factor: m

Type = 5 (Dominance 1)\cr
Genotype transformation:\cr
0 = -2*p[j]^2\cr
1 = 2*p[j]*(1-p[j])\cr
2 = -2*(1-p[j])^2\cr
Scaling factor: 4*sum(p^2*(1-p)^2).

Type = 6 (Dominance 2)\cr
Genotype transformation:\cr
0 = -2*p[j]^2\cr
1 = 2*p[j]*(1-p[j])\cr
2 = -2*(1-p[j])^2\cr
Scaling factor:\cr
\eqn{tr(\mathbf{MDM}')^{-1}n}.

The idrow and idcol arguments are intended for faster calculations when the user wishes to obtain a non-symmetric relationship matrix with a set of individuals in rows and a different set of individuals in columns. For example, training x testing samples, males x females, etc. Currently only supported if freq != NULL and for type 3 matrices.

}
\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com> \cr
Marco Milanesi <marco.milanesi.mm@gmail.com>
}
\references{
J. F. Taylor. Implementation and accuracy of genomic selection. Aquaculture 2014. 420, S8-S14.
P. M. VanRaden. Efficient methods to compute genomic predictions. J. Dairy. Sci. 2008. 91:4414-4423.
J. Yang et al. GCTA: A Tool for Genome-wide Complex Trait Analysis. Am. J. Hum. Genet. 2011. 88:76–82.
Z. G. Vitezica et al. On the additive and dominant variance and covariance of individuals within the genomic selection scope. Genetics. 2013. 195:1223-1230.
}
\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy plink data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "plink",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Copy metadata in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "meta",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load plink data
# plink <- ghap.loadplink("example")
# 
# # Load phenotype data
# df <- read.table(file = "example.phenotypes", header=T)
# 
# ### RUN ###
# 
# # Subset pure1 population
# pure1 <- plink$id[which(plink$pop == "Pure1")]
# plink <- ghap.subset(object = plink, ids = pure1, variants = plink$marker)
# 
# # Compute different types of relationship matrices
# K1 <- ghap.kinship(plink, type = 1) # General additive 1
# K2 <- ghap.kinship(plink, type = 2) # General additive 2
# K3 <- ghap.kinship(plink, type = 3) # VanRaden 2008
# K4 <- ghap.kinship(plink, type = 4) # GCTA GRM
# K5 <- ghap.kinship(plink, type = 5) # Dominance 1
# K6 <- ghap.kinship(plink, type = 6) # Dominance 2


}
