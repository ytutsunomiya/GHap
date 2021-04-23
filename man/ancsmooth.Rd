\name{ghap.ancsmooth}
\alias{ghap.ancsmooth}

\title{
 Smoothing of haplotype ancestry predictions
}
\description{
  Given ancestry predictions obtained with the \code{\link{ghap.anctest}} or \code{\link{ghap.ancsvm}} functions, overlapping classifications are smoothed to refine the boundaries of recombination breakpoints.
}
\usage{
 ghap.ancsmooth(phase, admix, ncores=1, verbose=TRUE)
}

\arguments{
  \item{phase}{
   A GHap.phase object.
  }
  \item{admix}{
   A data frame containing ancestry classifications, such as supplied by the \code{\link{ghap.anctest}} or \code{\link{ghap.ancsvm}} functions.
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
The function returns three dataframes: 'proportions1', 'proportions2' and 'haplotypes'. The 'proportions1' dataframe contains the following columns:
 \item{POP}{
  Original population label.
 }
 \item{ID}{
  Individual name.
 }
 \item{...}{
  A number of columns giving the predicted ancestry proportions.
 }
 \item{UNK}{
  The proportion of the genome without ancestry assignment.
 }

The 'proportions2' dataframe is similar to 'proportions1', expect that ancestry contributions are re-calibrated using only genome segments with ancestry assignments (therefore does not include the 'UNK' column). The 'haplotypes' dataframe contains the following columns:

 \item{POP}{
  Original population label.
 }
 \item{ID}{
  Individual name.
 }
 \item{HAP}{
  Haplotype number.
 }
 \item{CHR}{
  Chromosome name.
 }
 \item{BP1}{
  Segment start position.
 }
 \item{BP2}{
  Segment end position.
 }
 \item{SIZE}{
  Segment size.
 }
 \item{ANCESTRY}{
  Predicted ancestry of the segment.
 }
}

\details{
 This function takes results from ancestry classifications provided by the \code{\link{ghap.anctest}} or \code{\link{ghap.ancsvm}} functions and converts them into runs of ancestry. Since the classifiers assume exactly one ancestry per HapBlock, segments encompassing breakpoints are miss-classified as pertaining to a single origin, as opposed to a recombinant mixture of hybrid ancestry. When \code{\link{ghap.anctest}}/\code{\link{ghap.ancsvm}} are ran with overlapping HapBlocks, the smoothing function interrogates the ancestry of each overlapped segment by majority voting of all blocks containing it. After the ancestry of all segments have been resolved, contiguous sites sharing the same classification are converted into runs or segments of ancestry (i.e., ancestry tracks), which comprise the final output ('haplotypes' dataframe). These segments are then used to predict ancestry contributions ('proportions1' and 'proportions2' dataframes).
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
# # Unsupervised analysis with best K
# prototypes <- ghap.anctrain(phase = phase, K = 3, ncores = 1)
# hapadmix <- ghap.anctest(phase = phase, blocks = blocks, prototypes = prototypes,
#                          test = unique(phase$id), ncores = 1)
# anctracks <- ghap.ancsmooth(phase = phase, admix = hapadmix, ncores = 1)
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
}

\seealso{
 \code{\link{ghap.anctrain}}, \code{\link{ghap.ancsvm}}, \code{\link{ghap.ancplot}}, \code{\link{ghap.ancmark}}
}
