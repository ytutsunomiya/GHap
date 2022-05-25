\name{ghap.ancplot}
\alias{ghap.ancplot}

\title{
 Barplot of predictions of ancestry proportions
}
\description{
  Given smoothed ancestry predictions obtained with the \code{\link{ghap.ancsmooth}} function, an admixture barplot is generated.
}
\usage{
 ghap.ancplot(ancsmooth, labels = TRUE,
              pop.ang = 45, group.ang = 0,
              colors = NULL, pop.order = NULL,
              sortby = NULL, use.unk = FALSE,
              legend = TRUE)
}

\arguments{
  \item{ancsmooth}{
   A list containing smoothed ancestry classifications, such as supplied by the \code{\link{ghap.ancsmooth}} function.
  }
  \item{labels}{
   A logic value indicating if population labels should be plotted (default = TRUE).
  }
  \item{pop.ang}{
   A numeric value representing the rotation of population labels in degrees (default = 45).
  }
  \item{group.ang}{
   A numeric value representing the rotation of group labels in degrees (default = 0).
  }
  \item{colors}{
   A character vector of colors to use for each ancestry label.
  }
  \item{pop.order}{
   A single character vector or a list of character vectors specifying the order of populations to plot (see details).
  }
  \item{sortby}{
   A character value indicating the ancestry label to use for sorting individuals within populations.
  }
  \item{use.unk}{
   A logical value indicating if the plot should be generated using missing values (default = TRUE).
  }
  \item{legend}{
   A logical value indicating if a legend should be included to the plot (default = TRUE).
  }
}

\details{
 This function takes smoothed ancestry classifications provided by the \code{\link{ghap.ancsmooth}} function and generates a traditional admixture/structure barplot. The argument pop.order allows the user to organize the displaying order of populations through a vector or a list. If a vector is provided, the populations are plotted following the order of elements within the vector. Otherwise, if a named list of vectors is provided, populations are first grouped by list elements and then displayed in the order they appear within their respective group vector.
 
 The same data could be used to generate a circular barplot, for example with the \pkg{BITE} package. 
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
 \code{\link{ghap.ancsmooth}}, \code{\link{ghap.ancmark}}}

\references{
Milanesi, M., Capomaccio, S., Vajana, E., Bomba, L., Garcia, J.F., Ajmone-Marsan, P., Colli, L., 2017. BITE: an R package for biodiversity analyses. bioRxiv 181610. https://doi.org/10.1101/181610

Y.T. Utsunomiya et al. Unsupervised detection of ancestry tracks with the GHap R package. Methods in Ecology and Evolution. 2020. 11:1448–54.

}
