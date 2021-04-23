\name{ghap.anc2plink}
\alias{ghap.anc2plink}

\title{
 Convert ancestry tracks to PLINK binary
}
\description{
 This function takes smoothed ancestry predictions obtained with the \code{\link{ghap.ancsmooth}} function and converts them to PLINK binary (bed/bim/fam) format.
}
\usage{
ghap.anc2plink(phase, ancsmooth, ancestry, outfile, freq=c(0,1),
               missingness=1, only.active.samples=TRUE,
               only.active.markers=TRUE, batchsize=NULL,
               binary=TRUE, ncores=1, verbose=TRUE)
}

\arguments{
  \item{phase}{
   A GHap.phase object.
  }
  \item{ancsmooth}{
   A list containing smoothed ancestry classifications, such as supplied by the \code{\link{ghap.ancsmooth}} function.
  }
  \item{ancestry}{
  Character value indicating which ancestry to count at each observed marked site.
  }
  \item{outfile}{
  Character value for the output file name.
  }
  \item{freq}{
   A numeric vector of length 2 specifying the range of ancestry frequency to be included in the output. Default is c(0,1), which includes all marked sites.
  }
  \item{missingness}{
   A numeric value providing the missingness threshold to exclude marked sites with poor ancestry assignments (default = 1, with all sites retained).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be used for haplotyping (default = TRUE).
  }
  \item{batchsize}{
   A numeric value controlling the number of haplotype blocks to be processed and written to output at a time (default = 500).
  }
  \item{binary}{
   A logical value specfying whether the output file should be binary (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}


\details{
 The returned file mimics a standard PLINK (Purcell et al., 2007; Chang et al., 2015) binary file (bed/bim/fam), where counts 0, 1 and 2 represent the number of alleles assigned to the selected ancestry. For compatibility with PLINK, counts are coded as NN, NH and HH genotypes (N = NULL and H = haplotype allele), as if ancestry counts were bi-alelic markers. This codification is acceptable for any given analysis relying on SNP genotype counts, as long as the user specifies that the analysis should be done using the H character as reference for counts. You can specify reference alleles using the .tref file in PLINK with the \emph{--reference-allele} command. This is desired for very large datasets, as softwares such as PLINK and GCTA (Yang et al., 2011) have faster implementations for regression, principal components and kinship matrix analyses. Optionally, the user can use binary = FALSE to replace the bed file with a plain txt with ancestry counts.
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
C. C. Chang et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015. 4, 7.

S. Purcell et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am. J. Hum. Genet. 2007. 81, 559-575.

J. Yang et al. GCTA: A tool for genome-wide complex trait analysis. Am. J. Hum. Genet. 2011. 88, 76-82.
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
# # Export African ancestry in MEX to PLINK binary
# phase <- ghap.subsetphase(phase = phase,
#                           ids = unique(phase$id[which(phase$pop == "MEX")]),
#                           markers = phase$marker)
# ghap.anc2plink(phase = phase, ancsmooth = anctracks, ancestry = "Africa",
#                outfile = "mex_africa", ncores = 1, batchsize = 1000)


}
