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

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}
