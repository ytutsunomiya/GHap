\name{ghap.karyoplot}
\alias{ghap.karyoplot}

\title{
 Individual chromosome painting
}
\description{
  Given smoothed ancestry predictions obtained with the \code{\link{ghap.ancsmooth}} function and the name of the indivual, an individual karyotype plot is generated.
}
\usage{
 ghap.karyoplot(ancsmooth, ids=NULL, colors=NULL, chr=NULL, 
                chr.line=10, plot.line=25, chr.ang=45, las=0)
}

\arguments{
  \item{ancsmooth}{
   A list containing smoothed ancestry classifications, such as supplied by the \code{\link{ghap.ancsmooth}} function.
  }
  \item{ids}{
   A character vector of individual(s) to plot. If NULL all the individuals will be plotted (default = NULL).
  }
  \item{colors}{
   A character vector of colors to use for each ancestry label (default = NULL).
  }
  \item{chr}{
   A vector with the chromosome(s) to plot. If NULL, all the chromosomes will be plotted (default = NULL).
  }
  \item{chr.line}{
   A numeric value representing the number of chromosomes per plot line (default = 10).
  }
  \item{plot.line}{
   A numeric value representing the distance of horizontal guide (default = 25).
  }
  \item{chr.ang}{
   A numeric value representing the rotation of chromosome labels in degrees (default = 45).
  }
  \item{las}{
   A numeric value representing the las of y-axes (default = 0).
  }
}

\details{
 This function takes smoothed ancestry classifications provided by the \code{\link{ghap.ancsmooth}} function and "paint" the chomosomes of one individual using the ancestry proportions. One or more individuals could be plotted in separated graph. 
}

\author{
Marco Milanesi <marco.milanesi.mm@gmail.com>
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy the example data to the current working directory
# exfiles <- ghap.makefile()
# file.copy(from = exfiles, to = "./")
# 
# # Compress phase data
# ghap.compress(input.file = "human", out.file = "human")
# 
# # Load phase data
# phase <- ghap.loadphase(input.file = "human")
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
# # Unsupervised analysis with best K
# admix.kmeans <- ghap.ancestry(phase = phase, blocks = blocks, test = unique(phase$id),
#                               method = "unsupervised", K = 3, ncores = 1)
# out.kmeans <- ghap.ancsmooth(phase = phase, admix = admix.kmeans, ncores = 1)
# 
# # Training and test set for supervised analysis
# pops <- phase$pop
# names(pops) <- phase$id
# labels <- pops
# labels[which(labels \%in\% c("CEU","TSI"))] <- "Europe"
# labels[which(labels \%in\% c("CHB","JPT","CHD"))] <- "Asia"
# labels[which(labels \%in\% c("YRI","LWK"))] <- "Africa"
# labels[which(labels \%in\% c("MKK","GIH","ASW","MEX"))] <- "Test"
# phase$pop <- labels
# train <- unique(phase$id[which(phase$pop != "Test")])
# test <- unique(phase$id[which(phase$pop == "Test")])
# 
# # Supervised analysis with best parameters
# admix.svm <- ghap.ancestry(phase = phase, blocks = blocks, train = train, test = test,
#                            method = "supervised", fit.train = TRUE, ncores = 1)
# phase$pop <- pops
# admix.svm$POP <- pops[admix.svm$ID]
# out.svm <- ghap.ancsmooth(phase = phase, admix = admix.svm, ncores = 1)
# 
# ### RUN ###
# 
# # Select some individuals
# ids <- c("NA19835", "NA12003", "NA18749", "NA17969", "NA21108",
#          "NA18998", "NA19394", "NA19678", "NA21453", "NA20529", "NA19223")
# 
# # Plot karyoplot 
# ghap.karyoplot(ancsmooth = out.kmeans, ids = ids[1], chr.line = 11,
#                plot.line = 50, las=1, chr=NULL)
# ghap.karyoplot(ancsmooth = out.svm, ids = ids[1], chr.line = 11,
#                plot.line = 50, las=1)

}

\seealso{
 \code{\link{ghap.anctrain}}, \code{\link{ghap.anctest}}, \code{\link{ghap.ancsvm}}, \code{\link{ghap.ancsmooth}}, \code{\link{ghap.ancmark}}, \code{\link{ghap.ancplot}}
}
