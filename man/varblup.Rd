\name{ghap.varblup}
\alias{ghap.varblup}

\title{
Convert BLUP of individuals into BLUP of variants
}
\description{
Given genomic estimated breeding values (GEBVs), compute Best Linear Unbiased Predictor (BLUP) solutions for variant effects.
}
\usage{
ghap.varblup(object, gebv, covmat, type = 2,
             only.active.variants = TRUE,
             weights = NULL, tol = 1e-12,
             vcp = NULL, errormat = NULL, 
             errorname = "", nlambda = 1000,
             ncores = 1, verbose = TRUE)
}

\arguments{
  \item{object}{
    A valid GHap object (phase, haplo or plink).
  }
  \item{gebv}{
   A named vector of genomic estimated breeding values.
  }
  \item{covmat}{
    An additive genomic relationship matrix, such as obtained with type=1 or type=2 in the \code{\link{ghap.kinship}} function. 
  }
  \item{type}{
    A numeric value indicating the type of relationship matrix (see details in the \code{\link{ghap.kinship}} function).
  }
  \item{only.active.variants}{
    A logical value specifying whether only active variants should be included in the calculations (default = TRUE).
  }
  \item{weights}{
    A numeric vector providing variant-specific weights.
  }
  \item{tol}{
  A numeric value specifying the scalar to add to the diagonal of the relationship matrix it is not inversible (default = 1e-10).
 }
 \item{vcp}{
    A numeric value for the variance in GEBVs.
 }
 \item{errormat}{
    A square error matrix for GEBVs. This matrix can be obtained with argument extras = "LHSi" in the \code{\link{ghap.lmm}} function. If provided, calculation of standard errors and test statistics for the variants is activated.
  }
  \item{errorname}{
    The name used for the random effect representing GEBVs in the \code{\link{ghap.lmm}} function. If the error matrix was imported from somewhere else, this argument can be ignored provided that the names in the error matrix match the ones in the relationship matrix.
  }
  \item{nlambda}{
  A numeric value for the number of variants to be used in the estimation of the inflation factor (default = 1000).
 }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
 \item{verbose}{
   A logical value specifying whether log messages should be printed (default = TRUE).
 }
}

\value{
The function returns a data frame with results from the genome-wide conversion of BLUP of individuals into BLUP of variants. If a GHap.haplo object is used, the first columns of the data frame will be:

 \item{CHR}{
  Chromosome name.
 }
 \item{BLOCK}{
  Block alias.
 }
 \item{BP1}{
  Block start position.
 }
 \item{BP2}{
  Block end position.
 }
 
 For GHap.phase and GHap.plink objects, the first columns will be:
 
 \item{CHR}{
  Chromosome name.
 }
 \item{MARKER}{
  Block start position.
 }
 \item{BP}{
  Block end position.
 }
 
 The remaining columns of the data frame will be equal for any class of GHap objects:
 
  \item{ALLELE}{
  Identity of the counted (A1 or haplotype) allele.
 }
 \item{FREQ}{
  Frequency of the allele.
 }
 \item{SCORE}{
  Estimated BLUP of the allele.
 }
 \item{VAR}{
  Variance in allele-specific breeding values.
 }
  \item{pVAR}{
  Proportion of variance explained by the allele.
 }
 \item{CENTER}{
  Average genotype (meaningful only for predictions with \code{\link{ghap.profile}}).
 }
 \item{SCALE}{
  A constant set to 1 (meaningful only for predictions with \code{\link{ghap.profile}}).
 }
 
 If an error matrix for GEBVs is provided through the 'errormat' argument, the following additional columns are included in the data frame:
 
 \item{SE}{
  Standard error for the BLUP of the allele.
 }
 \item{CHISQ.EXP}{
  Expected values for the test statistics.
 }
 \item{CHISQ.OBS}{
  Observed value for the test statistics.
 }
 \item{CHISQ.GC}{
  Test statistics scaled by the inflation factor (Genomic Control). Inflation is computed through regression of observed quantiles onto expected quantiles. In order to avoid overestimation by variants rejecting the null hypothesis, a random sample of variants (with size controled via the nlambda argument) is taken within three standard deviations from the mean of the distribution of test statistics.
 }
 \item{LOGP}{
  log10(1/P) or -log10(P) for the BLUP of the allele.
 }
 \item{LOGP.GC}{
  log10(1/P) or -log10(P) for the BLUP of the allele (scaled by the inflation factor).
 }

}

\details{

The function uses the equation:

\deqn{\mathbf{\hat{a}} = q\mathbf{DM}^T\mathbf{K}^{-1}\mathbf{\hat{u}}}

where \eqn{\mathbf{M}} is the \emph{N} x \emph{H} centered matrix of variants observed for \emph{N} individuals and \emph{H} variants, \eqn{\mathbf{D} = diag(d_i)}, \eqn{d_i} is the weight of variant \emph{i} (default \eqn{d_i = 1}), \eqn{q} is the inverse weighted sum of variances in the columns of \eqn{\mathbf{M}}, \eqn{\mathbf{K}} is the additive genomic relationship matrix and \eqn{\hat{u}} is the vector of GEBVs.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
I. Stranden and D.J. Garrick. Technical note: derivation of equivalent computing algorithms for genomic predictions and reliabilities of animal merit. J Dairy Sci. 2009. 92:2971-2975.

J.L.G. Duarte et al. Rapid screening for phenotype-genotype associations by linear transformations of genomic evaluations. BMC Bioinformatics. 2014, 15:246.
}
