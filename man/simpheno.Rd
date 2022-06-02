\name{ghap.simpheno}
\alias{ghap.simpheno}

\title{
Quantitative trait simulation using real genotype data
}
\description{
Simulates phenotypes from a quantitative trait with arbitrary variant-specific heritabilities.
}
\usage{
ghap.simpheno(object, h2, r2 = 0, nrep = 1,
              balanced = TRUE, seed = NULL, 
              only.active.samples = TRUE,
              ncores = 1, verbose = TRUE)
}

\arguments{
  \item{object}{
   A valid GHap object (phase or plink).
  }
  \item{h2}{
   A named numeric value specifying the heritability per variant. The sum of variant-specific heritabilities will be set as the narrow-sense heritability, and must not exceed 1.
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
  \item{seed}{
    A numeric value used to set the random number generation state (default = NULL). This is useful for reproducibility of the results.
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be used for calculations (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
The function returns a data frame with items:
  \item{POP}{
   Original population label.
  }
  \item{ID}{
   Individual name.
  }
  \item{PHENO}{
   Phenotypic observation.
  }
  \item{TBV}{
   True breeding value.
  }
  \item{REP}{
   Permanent environmental effect (only present if nrep > 1).
  }
  \item{RESIDUAL}{
   Residual value.
  }
  
}

\details{

The simulation considers the model:

\deqn{\mathbf{y} = \mathbf{Zu} + \mathbf{Zp} + \mathbf{e}}

where \eqn{\mathbf{u}} is a vector of breeding values, \eqn{\mathbf{p}} is a vector of permanent environmental effects, \eqn{\mathbf{Z}} is an incidence matrix mapping \eqn{\mathbf{y}} to \eqn{\mathbf{u}} and \eqn{\mathbf{p}}, and \eqn{\mathbf{e}} is the vector of residuals. True breeding values are computed from the sum of causal variant effects specified in the 'h2' argument. Both the residual and permanent environmental effects are sampled from normal distributions.
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy the example data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "phase",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load phase data
# phase <- ghap.loadphase("example")
# 
# ### RUN ###
# 
# # Subset Pure1 population
# pure1 <- unique(phase$id[which(phase$pop == "Pure1")])
# phase <- ghap.subset(object = phase, ids = pure1,
#                      variants = phase$marker)
# freq <- ghap.freq(object = phase, type = "maf")
# 
# # Heritability   = 0.3
# # Number of QTLs = 1000
# # Major QTLs     = 0
# # Records per id = 1
# nqtl <- 1000
# h2 <- 0.3
# mkr <- sample(names(freq[which(freq > 0.05)]), size = nqtl)
# eff <- runif(n = nqtl, min = 0, max = 1)
# eff <- h2*eff/sum(eff)
# names(eff) <- mkr
# df1 <- ghap.simpheno(object = phase, h2 = eff)
# 
# # Heritability   = 0.5
# # Number of QTLs = 100
# # Major QTLs     = 1
# # Records per id = 5 (balanced)
# # Repeatability = 0.2
# nqtl <- 100
# h2 <- 0.4
# r2 <- 0.2
# reps <- 5
# mkr <- sample(names(freq[which(freq > 0.05)]), size = nqtl)
# eff <- runif(n = nqtl, min = 0, max = 1)
# eff <- h2*eff/sum(eff)
# eff[which(eff == min(eff))] <- 0.1
# names(eff) <- mkr
# df2 <- ghap.simpheno(object = phase, h2 = eff,
#                      r2 = r2, nrep = reps)
# 
# # Heritability   = 0.5
# # Number of QTLs = 100
# # Major QTLs     = 1
# # Records per id = 5 (unbalanced)
# # Repeatability = 0.2
# nqtl <- 100
# h2 <- 0.4
# r2 <- 0.2
# reps <- 5
# mkr <- sample(names(freq[which(freq > 0.05)]), size = nqtl)
# eff <- runif(n = nqtl, min = 0, max = 1)
# eff <- h2*eff/sum(eff)
# eff[which(eff == min(eff))] <- 0.1
# names(eff) <- mkr
# df3 <- ghap.simpheno(object = phase, h2 = eff, r2 = r2,
#                      nrep = reps, balanced = FALSE)


}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}
