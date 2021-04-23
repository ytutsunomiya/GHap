\name{ghap.profile}
\alias{ghap.profile}

\title{
Haplotype allele profile
}
\description{
Given a data.frame of user-defined haplotype allele scores, compute individual profiles.
}
\usage{
ghap.profile(score, haplo, only.active.samples = TRUE, batchsize = NULL,
             ncores = 1, verbose = TRUE)
}

\arguments{
  \item{score}{
   A data.frame containing columns: BLOCK, CHR, BP1, BP2, ALLELE, SCORE, CENTER and SCALE.
  }
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{only.active.samples}{
    A logical value specifying whether calculations should be reported only for active samples (default = TRUE).
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
The function returns a data.frame with the following columns:
 \item{POP}{
 Population ID.
 }
 \item{ID}{
 Individual name.
 }
 \item{PROFILE}{
 Individual profile.
 }
}

\details{
The profile for each individual is calculated as sum(b*(x-c)/s), where x is a vector of number of copies of each haplotype allele, c is a constant to center the genotypes (taken from the CENTER column of the score dataframe), s is a constant to scale the genotypes (taken from the SCALE column of the score dataframe), and b is a vector of user-defined scores for each haplotype allele (taken from the SCORE column of the score dataframe). If no centering or scaling is required, the user can set the CENTER and SCALE columns to 0 and 1, respectively. By default, if scores are provided for only a subset of the haplotype alleles, the missing alleles scores will be set to zero. This function has the same spirit as the profiling routine implemented in the \emph{score} option in PLINK (Purcell et al., 2007; Chang et al., 2015).
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>\cr
Marco Milanesi <marco.milanesi.mm@gmail.com>
}

\references{
C. C. Chang et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015. 4, 7.

S. Purcell et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am. J. Hum. Genet. 2007. 81, 559-575.
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
# # BLUP GWAS
# gebvs <- model$random$individual
# Kinv <- ghap.kinv(K)
# gwas.blup <- ghap.blup(gebvs = gebvs, haplo = haplo,
#                        ncores = 1, invcov = Kinv)
# 
# 
# ### RUN ###
# 
# # Scoring
# score <- ghap.profile(score = gwas.blup, haplo = haplo, ncores = 1)
# plot(score$SCORE, model$random$individual); abline(0,1)

}
