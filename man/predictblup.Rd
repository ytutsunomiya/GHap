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
# # Load phenotype and pedigree data
# df <- read.table(file = "example.phenotypes", header=T)
# 
# ### RUN ###
# 
# # Subset individuals from the pure1 population
# pure1 <- plink$id[which(plink$pop == "Pure1")]
# plink <- ghap.subset(object = plink, ids = pure1, variants = plink$marker)
# 
# # Subset markers with MAF > 0.05
# freq <- ghap.freq(plink)
# mkr <- names(freq)[which(freq > 0.05)]
# plink <- ghap.subset(object = plink, ids = pure1, variants = mkr)
# 
# # Compute genomic relationship matrix
# # Induce sparsity to help with matrix inversion
# K <- ghap.kinship(plink, sparsity = 0.01)
# 
# # Fit mixed model
# df$rep <- df$id
# model <- ghap.lmm(formula = pheno ~ 1 + (1|id) + (1|rep),
#                   data = df,
#                   covmat = list(id = K, rep = NULL),
#                   extras = "LHSi")
# refblup <- model$random$id$Estimate
# names(refblup) <- rownames(model$random$id)
# 
# # Predict blup of reference and test individuals
# blup <- ghap.predictblup(refblup, vcp = model$vcp$Estimate[1],
#                          covmat = as.matrix(K),
#                          errormat = model$extras$LHSi,
#                          errorname = "id")
# 
# # Compare predictions
# plot(blup$Estimate, model$random$id$Estimate)
# abline(0,1)
}
