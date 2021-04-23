\name{ghap.ancsvm}
\alias{ghap.ancsvm}

\title{
 SVM-based predictions of haplotype ancestry
}
\description{
  This function uses Support Vector Machines (SVM) to predict ancestry of haplotype alleles in test samples.
}
\usage{
 ghap.ancsvm(phase, blocks, test = NULL, train = NULL,
             cost = 1, gamma = NULL, tune = FALSE,
             only.active.samples = TRUE, only.active.markers = TRUE,
             ncores = 1, verbose = TRUE)
}

\arguments{

  \item{phase}{
   A GHap.phase object.
  }
  \item{blocks}{
   A data frame containing block boundaries, such as supplied by the \code{\link{ghap.blockgen}} function.
  }
  \item{test}{
   Character vector of individuals to test.
  }
  \item{train}{
   Character vector of individuals to use as reference samples.
  } 
  \item{cost}{
   A numeric value specifying the C constant of the regularization term in the Lagrange formulation.
  }
  \item{gamma}{
   A numeric value specifying the gamma parameter of the RBF kernel (default = 1/blocksize).
  }
  \item{tune}{
   A logical value specfying if a grid search is to be performed for parameters (default = FALSE).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in predictions (default = TRUE).
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be used for predictions (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
  
}

\value{
If ran with tune = FALSE, the function returns a dataframe with the following columns:
 \item{BLOCK}{
  Block alias.
 }
 \item{CHR}{
  Chromosome name.
 }
 \item{BP1}{
  Block start position.
 }
 \item{BP2}{
  Block end position.
 }
 \item{POP}{
  Original population label.
 }
 \item{ID}{
  Individual name.
 }
 \item{HAP1}{
  Predicted ancestry of haplotype 1.
 }
 \item{HAP2}{
  Predicted ancestry of haplotype 2.
 }
 
If tune = TRUE, the function returns a dataframe with the following columns:

\item{cost}{
  The candidate value of the C constant.
}

\item{gamma}{
  The canidate value of the gamma parameter.
}

\item{accuracy}{
  The percentage of correctly assigned ancestries.
}

}

\details{
This function predicts haplotype allele ancestry using Support Vector Machines (SVM) together with a Gaussian Radial Basis Function (RBF) kernel. The user is required to specify the C constant of the regularization term in the Lagrange formulation (default cost = 1) and the gamma parameter (default gamma = 1/blocksize) of the RBF kernel. \cr
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
# ### RUN ###
# 
# # Calculate marker density
# mrkdist <- diff(phase$bp)
# mrkdist <- mrkdist[which(mrkdist > 0)]
# density <- mean(mrkdist)
# 
# # Generate blocks for admixture events up to g = 10 generations in the past
# # Assuming mean block size in Morgans of 1/(2*g)
# # Approximating 1 Morgan ~ 100 Mbp
# g <- 10
# window <- (100e+6)/(2*g)
# window <- ceiling(window/density)
# step <- ceiling(window/4)
# blocks <- ghap.blockgen(phase, windowsize = window, slide = step, unit = "marker")
# 
# # User-defined labels
# pops <- phase$pop
# labels <- pops
# labels[which(labels \%in\% c("CEU","TSI"))] <- "Europe"
# labels[which(labels \%in\% c("CHB","JPT","CHD"))] <- "Asia"
# labels[which(labels \%in\% c("YRI","LWK"))] <- "Africa"
# labels[which(labels \%in\% c("MKK","GIH","ASW","MEX"))] <- "Test"
# phase$pop <- labels
# train <- unique(phase$id[which(phase$pop != "Test")])
# test <- unique(phase$id)
# 
# # Tune supervised analysis
# ranblocks <- sample(x = 1:nrow(blocks), size = 5, replace = FALSE)
# tunesvm <- ghap.ancsvm(phase = phase, blocks = blocks[ranblocks,], train = train,
#                        gamma = 1/window*c(0.1,1,10), tune = TRUE, ncores = 1)
# 
# # Supervised analysis with default parameters
# hapadmix <- ghap.ancsvm(phase = phase, blocks = blocks, test = test,
#                         train = train, ncores = 1)
# phase$pop <- pops

}

\seealso{
 \code{\link{svm}}, \code{\link{ghap.ancsmooth}}, \code{\link{ghap.ancplot}}, \code{\link{ghap.ancmark}}
}

\references{
R. J. Haasl  et al. Genetic ancestry inference using support vector machines, and the active emergence of a unique American population. Eur J Hum Genet. 2013. 21(5):554-62.

D. Meyer et al. e1071: Misc Functions of the Department of Statistics, Probability Theory Group (e1071). TU Wien. 2019 R Package Version 1.7-0.1. http://cran.r-project.org/web/packages/e1071/index.html.
}