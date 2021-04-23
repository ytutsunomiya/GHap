\name{ghap.ancmark}
\alias{ghap.ancmark}

\title{
 Per marker ancestry proportions
}
\description{
   Given smoothed ancestry predictions obtained with the \code{\link{ghap.ancsmooth}} function, per marker ancestry proportions are calculated across selected individuals.
}
\usage{
 ghap.ancmark(phase, ancsmooth, ids)
}

\arguments{
  \item{phase}{
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
# # Supervised analysis with user-defined labels
# pops <- phase$pop
# labels <- pops
# labels[which(labels \%in\% c("CEU","TSI"))] <- "Europe"
# labels[which(labels \%in\% c("CHB","JPT","CHD"))] <- "Asia"
# labels[which(labels \%in\% c("YRI","LWK"))] <- "Africa"
# labels[which(labels \%in\% c("MKK","GIH","ASW","MEX"))] <- "Test"
# phase$pop <- labels
# train <- unique(phase$id[which(phase$pop != "Test")])
# prototypes <- ghap.anctrain(phase = phase, train = train,
#                             method = "supervised", ncores = 1)
# hapadmix <- ghap.anctest(phase = phase, blocks = blocks, prototypes = prototypes,
#                          test = unique(phase$id), ncores = 1)
# phase$pop <- pops
# anctracks <- ghap.ancsmooth(phase = phase, admix = hapadmix, ncores = 1)
# 
# ### RUN ###
# 
# # Get per marker ancestry proportions for MKK
# ancmark <- ghap.ancmark(phase = phase, ancsmooth = anctracks,
#                         ids = unique(phase$id[which(phase$pop == "MKK")]))
# plot(ancmark$BP/1e+6, ancmark$Africa, ylim = c(0,100), las=1, type="l", col="blue",
#      xlab = "Chromosome 2 (Mbp)", ylab = "Ancestry (\%)", main = "Unsupervised")
# points(ancmark$BP/1e+6, ancmark$Asia, type="l", col="red")
# points(ancmark$BP/1e+6, ancmark$Europe, type="l", col="green")
# points(ancmark$BP/1e+6, ancmark$UNK, type="l", col="grey")

}

\seealso{
 \code{\link{ghap.ancsmooth}}
}
