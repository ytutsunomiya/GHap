\name{ghap.anctest}
\alias{ghap.anctest}

\title{
 Prediction of haplotype ancestry
}
\description{
  This function uses prototype alleles to predict ancestry of haplotypes in test samples.
}
\usage{
 ghap.anctest(object, blocks = NULL,
              prototypes, test = NULL,
              only.active.samples = TRUE,
              only.active.markers = TRUE,
              ncores = 1, verbose = TRUE)
}

\arguments{

  \item{object}{
   A GHap.phase object.
  }
  \item{blocks}{
   A data frame containing block boundaries, such as supplied by the \code{\link{ghap.blockgen}} function.
  }
  \item{prototypes}{
   A data frame containing prototype alleles, such as supplied by the \code{\link{ghap.anctrain}} function.
  }
  \item{test}{
   Character vector of individuals to test. All active individuals are used if this vector is not provided.
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
The function returns a dataframe with the following columns:

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
}

\details{
For each interrogated block, tested haplotypes are assigned to their nearest centroids (i.e., the pseudo-lineages with the smallest Euclidean distances). If no blocks are supplied, the function automatically builds blocks compatible with admixture up to 10 generations in the past based on intermarker distances. This has been chosen according to simulation results, where the use of haplotype blocks compatible with recent admixture (~10 generations) retained reasonable accuracy across most scenarios, regardless of the age of admixture.
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy phase data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "phase",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load phase data
# 
# phase <- ghap.loadphase("example")
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
# blocks <- ghap.blockgen(phase, windowsize = window,
#                         slide = step, unit = "marker")
# 
# # BestK analysis
# bestK <- ghap.anctrain(object = phase, K = 5, tune = TRUE)
# plot(bestK$ssq, type = "b", xlab = "K",
#      ylab = "Within-cluster sum of squares")
# 
# # Unsupervised analysis with best K
# prototypes <- ghap.anctrain(object = phase, K = 2)
# hapadmix <- ghap.anctest(object = phase,
#                          blocks = blocks,
#                          prototypes = prototypes,
#                          test = unique(phase$id))
# anctracks <- ghap.ancsmooth(object = phase, admix = hapadmix)
# ghap.ancplot(ancsmooth = anctracks)
# 
# # Supervised analysis
# train <- unique(phase$id[which(phase$pop != "Cross")])
# prototypes <- ghap.anctrain(object = phase, train = train,
#                             method = "supervised")
# hapadmix <- ghap.anctest(object = phase,
#                          blocks = blocks,
#                          prototypes = prototypes,
#                          test = unique(phase$id))
# anctracks <- ghap.ancsmooth(object = phase, admix = hapadmix)
# ghap.ancplot(ancsmooth = anctracks)

}

\seealso{
 \code{\link{ghap.anctrain}}, \code{\link{ghap.ancsmooth}}, \code{\link{ghap.ancplot}}, \code{\link{ghap.ancmark}}
}

\references{
Y.T. Utsunomiya et al. Unsupervised detection of ancestry tracks with the GHap R package. Methods in Ecology and Evolution. 2020. 11:1448–54.
}
