\name{ghap.simpheno}
\alias{ghap.simpheno}

\title{
Quantitative trait simulation using real HapGenotype data
}
\description{
Simulates phenotypes from a quantitative trait with arbitrary major alleles.
}
\usage{
ghap.simpheno(haplo, kinship, h2, g2, r2=0, nrep=1,
              balanced=TRUE, major=NULL, seed=NULL, ncores=1)
}

\arguments{
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{kinship}{
    Covariance matrix to be used in polygenic effects simulation.
  }
  \item{h2}{
   A numeric value specifying the heritability.
  }
  \item{g2}{
   A numeric vector specifying the proportion of genetic variance explained by each major allele. The sum of the proportions must not exceed 1. In cases where the sum is less than 1, polygenic effects are simulated such that the remaining variance is uniformly distributed throughout the genome.
  }
  \item{r2}{
   A numeric value specifying the repeatability (default = 0). Only relevant if nrep > 1.
  }
  \item{nrep}{
   A numeric value specifying the number of repeated measures per subject.
  }
  \item{balanced}{
   A logical value specifying whether the output data should be balanced (default = TRUE). If balanced = FALSE, the number of repeated measures per subject will be heterogeneous, following a uniform distribution with minimum zero and maximum nrep. Only relevant if nrep > 1.
  }
  \item{major}{
   A numeric vector specifying the indices of the alleles to be considered as major.
  }
  \item{seed}{
    A numeric value used to set the random number generation state (default = NULL). This is useful for reproducibility of the results.
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
}

\value{
The function returns a list with items:
  \item{h2}{
   A numeric value specifying the heritability.
  }
  \item{g2}{
   A numeric vector specifying the proportion of genetic variance explained by each major HapAllele.
  }
  \item{major}{
   A numeric vector specifying the indices of the HapAlleles to be considered as major.
  }
  \item{major.effect}{
   A numeric vector containing the simulated major HapAllele effects.
  }
  \item{u}{
   A numeric vector containing breeding values.
  }
  \item{p}{
   A numeric vector containing permanent environmental effects. Suppressed if r2 = NULL and nrep = 1.
  }
  \item{varu}{
   A numeric value corresponding to the genetic variance.
  }
  \item{varp}{
   A numeric vector corresponding to the variance in permanent environmental effects. Suppressed if nrep = 1.
  }
  \item{vare}{
   A numeric value corresponding to the residual variance.
  }
  \item{data}{
   A data.frame containing columns: phenotype = a numeric vector containing the simulated phenotypes; individual = a character vector containing the IDs of the corresponding phenotypes.
  }

  
}

\details{

The simulation considers the model:

\deqn{\mathbf{y} = \mathbf{Zu} + \mathbf{Zp} + \mathbf{e}}

where \eqn{\mathbf{u}} is a vector of breeding values, \eqn{\mathbf{p}} is a vector of permanent environmental effects, \eqn{\mathbf{Z}} is an incidence matrix mapping \eqn{\mathbf{y}} to \eqn{\mathbf{u}} and \eqn{\mathbf{p}}, and \eqn{\mathbf{e}} is the vector of residuals. Breeding values are assumed:

\deqn{\mathbf{u} = \mathbf{H}\mathbf{a} + \mathbf{g}}

where \eqn{\mathbf{a}} is the vector of major allele effects, \eqn{\mathbf{H}} is the centered matrix of major allele counts, and \eqn{\mathbf{g}} is a vector of polygenic effects.

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
#
# ### RUN ###
#
# # Subset common haplotypes in Europeans
# EUR.ids <- haplo$id[which(haplo$pop \%in\% c("TSI","CEU"))]
# haplo <- ghap.subsethaplo(haplo,EUR.ids,haplo$allele.in)
# hapstats <- ghap.hapstats(haplo, ncores = 1)
# common <- hapstats$TYPE \%in\% c("REGULAR","MAJOR") &
#  hapstats$FREQ > 0.05 &
#  hapstats$FREQ < 0.95
# haplo <- ghap.subsethaplo(haplo,EUR.ids,common)
# 
# #Compute relationship matrix
# K <- ghap.kinship(haplo, ncores = 1)
# 
# # Quantitative trait with 50\% heritability
# # Unbalanced repeated measurements (0 to 30)
# # Two major haplotypes accounting for 50\% of the genetic variance
# myseed <- 123456789
# set.seed(myseed)
# major <- sample(which(haplo$allele.in == TRUE), size = 2)
# g2 <- runif(n = 2, min = 0, max = 1)
# g2 <- (g2/sum(g2))*0.5
# sim <- ghap.simpheno(haplo, kinship = K, h2 = 0.5, g2 = g2, nrep = 30,
#                      balanced = FALSE, major = major, seed = myseed)


}
