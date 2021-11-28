\name{ghap.remlci}
\alias{ghap.remlci}

\title{
Confidence intervals for functions of variance components
}
\description{
Approximation of standard errors and confidence intervals of arbitrary functions of variance components using a sampling-based method.
}
\usage{
ghap.remlci(fun, vcp, ai, n = 10000,
            conf.level = 0.95,
            include.samples = FALSE,
            ncores = 1)
}

\arguments{
  \item{fun}{
   Function of variance components. See details.
  }
  \item{vcp}{
   A matrix or dataframe containing variance components, such as supplied by the \code{\link{ghap.lmm}} function.
 }
  \item{ai}{
   The inverse of the average information matrix, such as supplied by the \code{\link{ghap.lmm}} function.
  }
  \item{conf.level}{
   A numeric value informing the confidence level (default = 0.95).
 }
 \item{include.samples}{
  A logical value indicating if samples of the likelihood function should be kept.
 }
 \item{n}{
  A numerical value giving the number of samples to draw from the likelihood function (default = 10000).
 }
 \item{ncores}{
   A numerical value specifying the number of cores to use in parallel computations (default = 1).
 }
}

\value{
The returned object is a list with the following items:
 \item{stde}{
  The standard error estimate.
 }
 \item{ci}{
  The lower and upper limits of the confidence interval.
 }
 \item{samples}{
  A vector containing the samples from the multivariate normal distribution. Only present in the output if the argument include.samples is set to TRUE.
 }

}

\details{

The function implements the method of Meyer & Houle (2013) to approximate standard errors and confidence intervals of arbitrary functions of variance components. The method consists in assuming that REML estimates of variance components asymptotically follow a multivariate normal distribution with mean equals to the REML estimates themselves and covariance matrix equals the inverse of the average information matrix. Following that assumption, samples are drawn from that distribution and the user-defined function is applied to the samples. Standard errors are obtained as the standard deviation of the resulting vector, and confidence intervals are derived from vector quantiles.

The user must provide a function of variance components to the fun argument. The function has to be set assuming that the components are stored in a vector named x. For example, if variance components come from an animal model and the user wishes to obtain standard errors and confidence intervals for the heritability,  simply use fun = function(x)\{x[1]/sum(x)\}.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{

K. Meyer, D. Houle. Sampling based approximation of confidence intervals for functions of genetic covariance matrices. Proc. Assoc. Adv. Anim. Breed. 2013. 20, 523–527
         
J. Jensen et al. Residual maximum likelihood estimation of (Co)variance components in multivariate mixed linear models using average information. J. Ind. Soc. Ag. Statistics 1997. 49, 215-236.

}
