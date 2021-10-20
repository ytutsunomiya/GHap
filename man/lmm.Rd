\name{ghap.lmm}
\alias{ghap.lmm}

\title{
Linear mixed model
}
\description{
Linear mixed model fitting for fixed effects, random effects and variance components.
}
\usage{
ghap.lmm(formula, data, covmat = NULL, weights = NULL,
         vcp.initial = NULL, vcp.estimate = TRUE,
         vcp.conv = 1e-12, errors = TRUE, em.reml = 10,
         tol = 1e-10, verbose = TRUE)
}

\arguments{
  \item{formula}{
   Formula describing the model. The synthax is consistent with lme4. The response is declared first, followed by the ~ operator. Predictors are then separated by + operators. Currently only random intercepts are supported, which are distinguished from fixed effects by the notation (1|x).
  }
  \item{covmat}{
   A list of inverse covariance matrices for each group of random effects. If a matrix is not defined for a given group, an identity matrix will be used.
 }
  \item{data}{
   A dataframe containing the data.
  }
  \item{weights}{
   A numeric vector with weights for observations. These weights are treated as diagonal elements of the inverse weight matrix. If not supplied, the analysis is carried out assuming all observations are equally important.
 }
 \item{vcp.initial}{
  A list of initial values for variance components. If not provided, the sample variance will be equally divided across the variance components.
 }
 \item{vcp.estimate}{
  A logical value specifying whether variance components should be estimated (default = TRUE). If FALSE, values passed to vcp.initial will be regarded as known variances.
 }
 \item{vcp.conv}{
  A numeric value specifying the convergence criterion for variance components (default = 1e-12).
 }
 \item{errors}{
  A logical value specifying whether standard errors should be computed (default = TRUE).
 }
 \item{em.reml}{
  A numeric value specifying the number of EM-REML iterations to carry out before switching to AI-REML (default = 10).
 }
 \item{tol}{
  A numeric value specifying the scalar to add to the diagonal of the left hand side of mixed models equations if it is not inversible (default = 1e-12).
 }
 \item{verbose}{
   A logical value specifying whether log messages should be printed (default = TRUE).
 }
}

\value{
The returned GHap.lmm object is a list with the following items:
 \item{fixed}{
  A dataframe containing estimates, standard errors, t values and p-values of fixed effects.
 }
 \item{random}{
  A list of dataframes, one for each group of random effects, containing estimates, standard errors and accuracy of random effects.
 }
 \item{fitted}{
  A dataframe with fitted values using all effects, only fixed effects and only random effects.
 }
 \item{residuals}{
  A dataframe with residuals from subtracting fitted values using all effects, only fixed effects and only random effects.
 }
 \item{vcp}{
  A dataframe with estimates and standard errors of variance components.
 }


}

\details{

The function fits mixed models using a combination of EM-REML, AI-REML and PCG. Random regression and multivariate analyses are currently not supported.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
D. Bates et al. Fitting Linear Mixed-Effects Models Using lme4. J. Stat. Soft., 67:1-48.
         
J. Jensen et al. Residual maximum likelihood estimation of (Co)variance components in multivariate mixed linear models using average information. J. Ind. Soc. Ag. Statistics 1997. 49, 215-236.

A. I. Vazquez. Technical note: An R package for fitting generalized linear mixed models in animal breeding. J. Anim. Sci. 2010. 88, 497-504.
}
