\name{ghap.profile}
\alias{ghap.profile}

\title{
Genomic profile
}
\description{
Given a data.frame of user-defined marker or haplotype allele scores, compute individual genomic profiles.
}
\usage{
ghap.profile(object, score, only.active.samples = TRUE,
             batchsize = NULL, ncores = 1, verbose = TRUE)
}

\arguments{

  \item{object}{
   A valid GHap object (phase, haplo or plink).
  }
  \item{score}{
   For HapAlleles (to use with GHap.haplo objects), the columns should be: BLOCK, CHR, BP1, BP2, ALLELE, SCORE, CENTER and SCALE. In the case of markers, the columns are: MARKER, CHR, BP, ALLELE, SCORE, CENTER and SCALE.
  }
  \item{only.active.samples}{
    A logical value specifying whether calculations should be reported only for active samples (default = TRUE).
  }
  \item{batchsize}{
   A numeric value controlling the number of HapAlleles or markers to be processed at a time (default = n/10).
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
The profile for each individual is calculated as sum(b*(x-c)/s), where x is a vector of number of copies of a reference allele, c is a constant to center the genotypes (taken from the CENTER column of the score dataframe), s is a constant to scale the genotypes (taken from the SCALE column of the score dataframe), and b is a vector of user-defined scores for each reference allele (taken from the SCORE column of the score dataframe). If no centering or scaling is required, the user can set the CENTER and SCALE columns to 0 and 1, respectively. By default, if scores are provided for only a subset of the HapAlleles or markers, the absent scores will be set to zero. If the input genomic data is a GHap.phase or GHap.plink object, the ALLELE column in the data frame is used as a guide to keep track of the direction of the scores. The default coding in GHap is A0/A0 = 0, A0|A1 = A1|A0 = 1 and A1/A1 = 2. For each marker where the declared allele is A0 instead of A1, that coding is flipped to A0/A0 = 2, A0|A1 = A1|A0 = 1 and A1/A1 = 1 so the profile is computed in the correct direction. Therefore, the user must be careful regarding the allele order, as well as the centering and scaling for each marker, since profiling is allele-oriented. This function has the same spirit as the profiling routine implemented in the \emph{score} option in PLINK (Purcell et al., 2007; Chang et al., 2015).
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
# # Copy plink data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "plink",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Copy metadata in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "meta",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load plink data
# plink <- ghap.loadplink("example")
# 
# # Load phenotype and pedigree data
# df <- read.table(file = "example.phenotypes", header=T)
# 
# ### RUN ###
# 
# # Subset individuals from the pure1 population
# pure1 <- plink$id[which(plink$pop == "Pure1")]
# plink <- ghap.subset(object = plink, ids = pure1, variants = plink$marker)
# 
# # Subset markers with MAF > 0.05
# freq <- ghap.freq(plink)
# mkr <- names(freq)[which(freq > 0.05)]
# plink <- ghap.subset(object = plink, ids = pure1, variants = mkr)
# 
# # Compute genomic relationship matrix
# # Induce sparsity to help with matrix inversion
# K <- ghap.kinship(plink, sparsity = 0.01)
# 
# # Fit mixed model
# df$rep <- df$id
# model <- ghap.lmm(formula = pheno ~ 1 + (1|id) + (1|rep),
#                   data = df,
#                   covmat = list(id = K, rep = NULL))
# refblup <- model$random$id$Estimate
# names(refblup) <- rownames(model$random$id)
# 
# # Convert blup of individuals into blup of variants
# mkrblup <- ghap.varblup(object = plink,
#                         gebv = refblup,
#                         covmat = K)
# 
# # Build GEBVs from variant effects and compare predictions
# gebv <- ghap.profile(object = plink, score = mkrblup)
# plot(gebv$SCORE, refblup); abline(0,1)

}
