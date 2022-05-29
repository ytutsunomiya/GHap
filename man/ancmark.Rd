\name{ghap.ancmark}
\alias{ghap.ancmark}

\title{
 Per marker ancestry proportions
}
\description{
   Given smoothed ancestry predictions obtained with the \code{\link{ghap.ancsmooth}} function, per marker ancestry proportions are calculated across selected individuals.
}
\usage{
 ghap.ancmark(object, ancsmooth, ids)
}

\arguments{
  \item{object}{
   A GHap.phase object.
  }
  \item{ancsmooth}{
   A list containing smoothed ancestry classifications, such as supplied by the \code{\link{ghap.ancsmooth}} function.
  }
  \item{ids}{
   A character vector specifying which individuals to use for the calculations.
  }
}

\value{
The function returns a dataframes containing the following columns:

 \item{CHR}{
  Chromosome name.
 }
 \item{MARKER}{
  Marker name.
 }
 \item{BP}{
  Marker position.
 }
 \item{...}{
  A number of columns (one for each ancestry label) giving the proportion of haplotypes carrying the respective ancestry label.
 }
}

\details{
 This function takes smoothed ancestry classifications provided by the \code{\link{ghap.ancsmooth}} function and calculates, for each marker, the proportion of haplotypes carrying each ancestry label. The resulting output serve as a proxy for locus-specific ancestry proportions.
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
# # Get per marker ancestry proportions for 'Pure1'
# pure1 <- unique(phase$id[which(phase$pop == "Pure1")])
# ancmark <- ghap.ancmark(object = phase,
#                         ancsmooth = anctracks,
#                         ids = pure1)
# 
# # Plot 'Pure2' introgression into 'Pure1'
# ghap.manhattan(data = ancmark, chr = "CHR",
#                bp = "BP", y = "Pure2", type = "h")

}

\seealso{
 \code{\link{ghap.ancsmooth}}
}
