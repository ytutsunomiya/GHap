\name{ghap.blup}
\alias{ghap.blup}

\title{
Convert breeding values into BLUP solutions of HapAllele effects
}
\description{
Given genomic estimated breeding values (GEBVs), compute Best Linear Unbiased Predictor (BLUP) solutions for HapAllele effects.
}
\usage{
ghap.blup(gebvs,haplo,invcov, gebvsweights = NULL, haploweights = NULL,
          only.active.alleles = TRUE, batchsize = NULL, ncores = 1, verbose = TRUE)
}

\arguments{
  \item{gebvs}{
   A vector of GEBVs. The vector must be named and all names must be present in the GHap.haplo object.
  }
  \item{haplo}{
   A GHap.haplo object.
  }
   \item{invcov}{
  The inverse covariance (i.e., inverse genomic kinship) matrix for GEBVs.
 }
  \item{gebvsweights}{
    A numeric vector providing individual-specific weights.
  }
  \item{haploweights}{
    A numeric vector providing HapAllele-specific weights.
  }
  \item{only.active.alleles}{
    A logical value specifying whether only active HapAlleles should be included in the calculations (default = TRUE).
  }
  \item{batchsize}{
   A numeric value controlling the number of HapAlleles to be processed at a time (default = nalleles/10).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
The function returns a dataframe with columns:
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
 \item{ALLELE}{
  Haplotype allele identity.
 }
 \item{SCORE}{
  BLUP for the random effect of the haplotype allele.
 }
 \item{FREQ}{
  Frequency of the haplotype allele.
 }
 \item{VAR}{
  Variance in allele-specific breeding values.
 }
  \item{pVAR}{
  Proportion of variance explained by the haplotype allele.
 }
 \item{CENTER}{
  Average genotype (meaningful only for predictions with \code{\link{ghap.profile}}).
 }
 \item{SCALE}{
  A constant set to 1 (meaningful only for predictions with \code{\link{ghap.profile}}).
 }

}

\details{

The function uses the equation:

\deqn{\mathbf{\hat{a}} = q\mathbf{DM}^T\mathbf{K}^{-1}\mathbf{\hat{u}}}

where \eqn{\mathbf{M}} is the \emph{N} x \emph{H} centered matrix of HapGenotypes observed for \emph{N} individuals and \emph{H} HapAlleles, \eqn{\mathbf{D} = diag(d_i)}, \eqn{d_i} is the weight of HapAllele \emph{i} (default \eqn{d_i = 1}), \eqn{q} is the inverse weighted sum of variances in the columns of \eqn{\mathbf{M}}, \eqn{\mathbf{K}} is the haplotype-based kinship matrix and \eqn{\hat{u}} is the vector of GEBVs. The permutation procedure consists in randomizing the vector \eqn{\hat{\mathbf{u}}} and computing the null statistic \eqn{max(\mathbf{\hat{a}})}. The permutation p-value is computed as the number of times the HapAllele effect was smaller than the null statistic.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
I. Stranden and D.J. Garrick. Technical note: derivation of equivalent computing algorithms for genomic predictions and reliabilities of animal merit. J Dairy Sci. 2009. 92:2971-2975.
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
# # Generate blocks of 5 markers sliding 5 markers at a time
# blocks <- ghap.blockgen(phase, windowsize = 5, slide = 5, unit = "marker")
# 
# # Haplotyping
# ghap.haplotyping(phase = phase, blocks = blocks, outfile = "human",
#                  binary = T, ncores = 1)
#
# # Load haplotype genotypes
# haplo <- ghap.loadhaplo(input.file = "human")
#
# # Subset common haplotypes in Europeans
# EUR.ids <- haplo$id[which(haplo$pop \%in\% c("TSI","CEU"))]
# haplo <- ghap.subsethaplo(haplo,EUR.ids,haplo$allele.in)
# hapstats <- ghap.hapstats(haplo, ncores = 1)
# common <- hapstats$TYPE \%in\% c("REGULAR","MAJOR") &
#  hapstats$FREQ > 0.05 &
#  hapstats$FREQ < 0.95
# haplo <- ghap.subsethaplo(haplo,EUR.ids,common)
# 
# # Compute relationship matrix
# K <- ghap.kinship(haplo, ncores = 1)
# 
# # Quantitative trait with 50\% heritability
# # Unbalanced repeated measurements (0 to 30)
# # Two major haplotypes accounting for 50\% of the genetic variance
# myseed <- 123456789
# set.seed(myseed)
# major <- sample(which(haplo$allele.in == TRUE), size = 2)
# g2 <- runif(n = 2, min = 0, max = 1)
# g2 <- (g2/sum(g2))*0.5
# sim <- ghap.simpheno(haplo, kinship = K, h2 = 0.5, g2 = g2, nrep = 30,
#                      balanced = FALSE, major = major, seed = myseed)
# 
# # Fit model using REML
# model <- ghap.lmm(fixed = phenotype ~ 1, random = ~ individual,
#                   covmat = list(individual = K), data = sim$data)
# 
#
# ### RUN ###
# 
# # BLUP GWAS
# gebvs <- model$random$individual
# gebvsw <- table(sim$data$individual)
# gebvsw <- gebvsw + mean(gebvsw)
# gebvsw <- gebvsw[names(gebvs)]
# Kinv <- ghap.kinv(K)
# gwas.blup <- ghap.blup(gebvs = gebvs, haplo = haplo, gebvsweights = gebvsw,
#                        ncores = 1, invcov = Kinv)
# plot(gwas.blup$BP1/1e+6,gwas.blup$pVAR*100,pch=20,
#      xlab="Position (in Mb)",ylab="Variance explained (%)")
# abline(v=haplo$bp1[major]/1e+6)
# 
# # BLUP with one update
# w <- gwas.blup$VAR*nrow(gwas.blup)
# K2 <- ghap.kinship(haplo=haplo, weights = w, ncores=4)
# Kinv2 <- ghap.kinv(K2)
# gwas.blup2 <- ghap.blup(gebvs = gebvs, haplo = haplo, invcov = Kinv2, ncores = 1,
#                         gebvsweights = gebvsw, haploweights = w)
# plot(gwas.blup2$BP1/1e+6,gwas.blup2$pVAR*100,pch=20,
#      xlab="Position (in Mb)",ylab="Variance explained (%)")
# abline(v=haplo$bp1[major]/1e+6)


}
