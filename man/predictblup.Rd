\name{ghap.predictblup}
\alias{ghap.predictblup}

\title{
  Predict BLUP from reference
}
\description{
  Prediction of BLUP values in test individuals based on reference individuals.
}
\usage{
  ghap.predictblup(refblup, vcp, covmat,
                   errormat = NULL,
                   errorname = "",
                   include.ref = TRUE,
                   diagonals = FALSE,
                   tol = 1e-12)
}

\arguments{
  \item{refblup}{
    A named numeric vector of reference BLUP values.
  }
  \item{vcp}{
    A numeric value for the variance in BLUP values.
  }
  \item{covmat}{
    A square matrix containing correlations among individuals. Both test and reference indiviudals must be present in the matrix.
  }
  \item{errormat}{
    A square error matrix for reference individuals. This matrix can be obtained with argument extras = "LHSi" in the \code{\link{ghap.lmm}} function.
  }
  \item{errorname}{
    The name used for the random effect in the \code{\link{ghap.lmm}} function. If the error matrix was imported from somewhere else, this argument can be ignored provided that the names in the error matrix match the ones in the covariance matrix.
  }
  \item{include.ref}{
    A logical value indicating if reference individuals should be included in the output (default = TRUE).
  }
  \item{diagonals}{
    A logical value indicating if diagonals of the covariance matrix should be used in calculations of accuracy and standard errors (default = FALSE). The default is to set diagonals to 1. For genomic estimated breeding values, using TRUE will account for inbreeding in the computation of accuracies and standard errors.
  }
  \item{tol}{
    A numeric value specifying the scalar to add to the diagonal of the covariance matrix if it is not inversible (default = 1e-12).
  }
}

\value{
  A data frame with predictions of BLUP values. If an error matrix is provided, standard errors and accuracies are also included.
}


\author{
  Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
  J.F. Taylor. Implementation and accuracy of genomic selection. Aquaculture 2014. 420, S8-S14.
}
