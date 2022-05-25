\name{ghap.anctrain}
\alias{ghap.anctrain}

\title{
 Construction of prototype alleles
}
\description{
  This function builds prototype alleles to be used in ancestry predictions.
}
\usage{
 ghap.anctrain(object, train = NULL,
               method = "unsupervised",
               K = 2, iter.max = 10, nstart = 10,
               nmarkers = 5000, tune = FALSE,
               only.active.samples = TRUE,
               only.active.markers = TRUE,
               batchsize = NULL, ncores = 1,
               verbose = TRUE)
}

\arguments{

  The following arguments are used by both the 'supervised' and 'unsupervised' methods:

  \item{object}{
   A GHap.phase object.
  }
  \item{train}{
   Character vector of individuals to use as reference samples. All active individuals are used if this    vector is not provided.
  } 
  \item{method}{
   Character value indicating which method to use: 'supervised' or 'unsupervised' (default).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in predictions (default = TRUE).
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be used for predictions (default = TRUE).
  }
  \item{batchsize}{
   A numeric value controlling the number of markers to be processed at a time (default = nmarkers/10).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
  
  The following arguments are only used by the 'unsupervised' method:
  
  \item{K}{
   A numeric value specifying the number of clusters in K-means (default = 2). Proxy for the number of ancestral populations.
  }
  \item{iter.max}{
   A numeric value specifying the maximum number of iterations of the K-means clustering (default = 10).
  }
  \item{nstart}{
   A numeric value specifying the number of independent runs of the K-means clustering (default = 10).
  }
  \item{nmarkers}{
   A numeric value specifying the number of seeding markers to be used by the K-means clustering (default = 10).
  }
  \item{tune}{
   A logical value specfying if a Best K analysis should be performed (default = FALSE).
  }

}

\value{
The function returns a dataframe with the first column giving marker names and remaining columns containing prototype alleles for each pseudo-lineage. If method 'unsupervised' is ran with tune = TRUE, the function returns the following list:

\item{ssq}{
  Within-cluster sum of squares for each value of K.
}

\item{chindex}{
  Calinski Harabasz Index for consecutive values of K.
}

\item{pchange}{
  Percent change in ssq for consecutive values of K.
}

}

\details{
This function builds prototype alleles (i.e., cluster centroids, representing lineage-specific allele frequencies) through two methods: \cr

The 'unsupervised' method uses the K-means clustering algorithm to group haplotypes into K pseudo-lineages. A random sample of seeding markers (default value of nmarkers = 5000) is used to group all 2*nsamples haplotypes in a user-specified number of clusters (default value of K = 2). Then, for each interrogated block, prototype alleles are built for every cluster using the arithmetic mean of observed haplotypes initially assigned to that cluster. If train = NULL, the function uses all active haplotypes to build prototype alleles. If the user is working with a severely unbalanced data set (ex. one population with a large number of individuals and others with few individuals), it is recommended that a vector of individual names is provided via the train argument such that prototype alleles are built using a more balanced subset of the data.

The 'supervised' method works in a similar way, but skips the K-means algorithm and uses population labels present in the GHap.phase object as clusters. \cr
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
 \code{\link{ghap.anctest}}, \code{\link{ghap.ancsmooth}}, \code{\link{ghap.ancplot}}, \code{\link{ghap.ancmark}}
}

\references{
Y.T. Utsunomiya et al. Unsupervised detection of ancestry tracks with the GHap R package. Methods in Ecology and Evolution. 2020. 11:1448–54.
}
