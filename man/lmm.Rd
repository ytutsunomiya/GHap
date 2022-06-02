\name{ghap.lmm}
\alias{ghap.lmm}

\title{
Linear mixed model
}
\description{
Linear mixed model fitting for fixed effects, random effects and variance components.
}
\usage{
ghap.lmm(formula, data, covmat = NULL,
         weights = NULL, vcp.initial = NULL,
         vcp.estimate = TRUE, vcp.conv = 1e-12,
         errors = TRUE, invcov = FALSE,
         em.reml = 10, tol = 1e-12, extras = NULL,
         verbose = TRUE)
}

\arguments{
  \item{formula}{
   Formula describing the model. The synthax is consistent with lme4. The response is declared first, followed by the ~ operator. Predictors are then separated by + operators. Currently only random intercepts are supported, which are distinguished from fixed effects by the notation (1|x).
  }
  \item{covmat}{
   A list of covariance matrices for each group of random effects. If a matrix is not defined for a given group, an identity matrix will be used. Inverse covariance matrices can also be provided, as long as argument invcov = TRUE is used.
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
 \item{invcov}{
  A logical value specifying whether the provided covariance matrices are already inverted (default = FALSE).
 }
 \item{em.reml}{
  A numeric value specifying the number of EM-REML iterations to carry out before switching to AI-REML (default = 10).
 }
 \item{tol}{
  A numeric value specifying the scalar to add to the diagonal of the left hand side of mixed models equations if it is not inversible (default = 1e-12).
 }
 \item{extras}{
  A character vector indicating extra output to be included in the results list. Currently supported extras are the full inverse of the coefficient matrix ("LHSi") and the phenotypic (co)variance matrix ("V").
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
 \item{extras}{
  A list of extra outputs requested by the user (see arguments for possible extras).
 }

}

\details{

The function fits mixed models using a combination of EM-REML, AI-REML and PCG. Random regression and multivariate analyses are currently not supported.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
J. Jensen et al. Residual maximum likelihood estimation of (Co)variance components in multivariate mixed linear models using average information. J. Ind. Soc. Ag. Statistics 1997. 49, 215-236.
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
# # Fit mixed model with variance components estimation
# df$rep <- df$id
# model1 <- ghap.lmm(formula = pheno ~ 1 + (1|id) + (1|rep),
#                    data = df,
#                    covmat = list(id = K, rep = NULL))
# 
# # Fit mixed model with fixed variance components
# df$rep <- df$id
# model2 <- ghap.lmm(formula = pheno ~ 1 + (1|id) + (1|rep),
#                    data = df,
#                    covmat = list(id = K, rep = NULL),
#                    vcp.initial = list(id = 0.4, rep = 0.2, Residual = 0.4),
#                    vcp.estimate = F)

}
