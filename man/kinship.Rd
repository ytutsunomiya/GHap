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
             ancestry=NULL, idrow = NULL,
             idcol = NULL, outputM = FALSE,
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
    A named numeric vector with (A1) allele frequencies computed in a reference sample. If not provided, frequencies will be computed using the active individuals. If type = 7, this argument must be a data.frame with rownames indicating marker names and one column per ancestral population.
  }
  \item{ancestry}{
    A data.frame containing individual ancestry proportions, with rownames indicating individual names and one column per ancestral population.
  }
  \item{idrow}{
    A character vector providing ids of individuals to be included in rows (default = NULL, see details).
  }
  \item{idcol}{
    A character vector providing ids of individuals to be included in columns (default = NULL, see details).
  }
  \item{outputM}{
    A logical value indicating if the transformed genotype matrix should be returned (default = FALSE, see details).
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

\deqn{\mathbf{K} = q^{-1}\mathbf{MDM}^{T}}
where \eqn{\mathbf{D} = diag(d_j)}, \eqn{d_j} is the weight of variant \emph{j} (default \eqn{d_j = 1}), and \eqn{q} is a scaling factor. The argument type controls the genotype transformation and the scaling factor, and includes the following options:

Type = 1 (General additive 1)\cr
Genotype transformation: \eqn{\mathbf{M}_{ij} - \text{mean}(\mathbf{M}_j)}\cr
Scaling factor: \eqn{tr(\mathbf{MDM}')^{-1}n}

Type = 2 (General additive 2)\cr
Genotype transformation: \eqn{(\mathbf{M}_{ij} - \text{mean}(\mathbf{M}_j))/\text{sd}(\mathbf{M}_j)}\cr
Scaling factor: \emph{m}

Type = 3 (VanRaden, 2008)\cr
Genotype transformation: \eqn{\mathbf{M}_{ij} - 2p_j}\cr
Scaling factor: \eqn{2\sum p_j(1-p_j)}

Type = 4 (GCTA)\cr
Genotype transformation: \eqn{(\mathbf{M}_{ij} - 2p_j)/\sqrt{2p_j(1-p_j)}}\cr
Scaling factor: \emph{m}

Type = 5 (Dominance 1)\cr
Genotype transformation:\cr
\eqn{\mathbf{M}_{ij} = 0 \rightarrow -2p_j^2}\cr
\eqn{\mathbf{M}_{ij} = 1 \rightarrow 2p_j(1-p)}\cr
\eqn{\mathbf{M}_{ij} = 2 \rightarrow -2(1-p_j)^2}\cr
Scaling factor: \eqn{4\sum p^2(1-p)^2}

Type = 6 (Dominance 2)\cr
Genotype transformation:\cr
\eqn{\mathbf{M}_{ij} = 0 \rightarrow -2p_j^2}\cr
\eqn{\mathbf{M}_{ij} = 1 \rightarrow 2p_j(1-p)}\cr
\eqn{\mathbf{M}_{ij} = 2 \rightarrow -2(1-p_j)^2}\cr
Scaling factor: \eqn{tr(\mathbf{MDM}')^{-1}n}

Type = 7 (Ancestry-weighted additive)\cr
Genotype transformation:\eqn{(\mathbf{M}_{ij} - 2p_{ij})/\sqrt{2p_{ij}(1-p_{ij})}}\cr
Scaling factor: \emph{m}

Where subscript \emph{i} refers to an individual. For relationship types 3 through 6, \eqn{p_j} is the counted allele frequency at variant site \emph{j}. For relationship type 7, we have \eqn{p_{ij} = \sum_{k=1}^{K} p_{jk}\pi_{ik}}, where \emph{K} is the number of ancestral populations, \eqn{p_{jk}} is the counted allele frequency at variant site \emph{j} in population \emph{k}, and \eqn{\pi_{ik}} is the ancestry proportion of population \emph{k} in individual \emph{i}. The idrow and idcol arguments are intended for faster calculations when the user wishes to obtain a non-symmetric relationship matrix with a set of individuals in rows and a different set of individuals in columns. For example, training x testing samples, males x females, etc. Currently only supported if freq != NULL and for type 3 matrices.

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
