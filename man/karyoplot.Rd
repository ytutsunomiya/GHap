\name{ghap.karyoplot}
\alias{ghap.karyoplot}

\title{
 Individual chromosome painting
}
\description{
  Given smoothed ancestry predictions obtained with the \code{\link{ghap.ancsmooth}} function and the name of the indivual, an individual karyotype plot is generated.
}
\usage{
 ghap.karyoplot(ancsmooth, ids = NULL,
                colors = NULL, chr = NULL, 
                chr.line = 10, plot.line = 25,
                chr.ang = 45, las = 0)
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
# 
# ### RUN ###
# 
# # Plot karyoplot
# pure1 <- unique(phase$id[which(phase$pop == "Pure1")])
# pure2 <- unique(phase$id[which(phase$pop == "Pure2")])
# cross <- unique(phase$id[which(phase$pop == "Cross")])
# ghap.karyoplot(ancsmooth = anctracks, ids = pure1[1],
#                chr.line = 11, plot.line = 50, las=1, chr=NULL)
# ghap.karyoplot(ancsmooth = anctracks, ids = pure2[1],
#                chr.line = 11, plot.line = 50, las=1, chr=NULL)
# ghap.karyoplot(ancsmooth = anctracks, ids = cross[1],
#                chr.line = 11, plot.line = 50, las=1, chr=NULL)

}

\seealso{
 \code{\link{ghap.anctrain}}, \code{\link{ghap.anctest}}, \code{\link{ghap.ancsvm}}, \code{\link{ghap.ancsmooth}}, \code{\link{ghap.ancmark}}, \code{\link{ghap.ancplot}}
}
