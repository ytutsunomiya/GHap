\name{ghap.getHinv}
\alias{ghap.getHinv}

\title{
Compute the inverse of H
}
\description{
This function combines an additive genomic relationship matrix with pedigree data to form the inverse of the H relationship matrix used in single-step GBLUP.
}
\usage{
ghap.getHinv(K, ped, include = NULL, depth = 3,
             alpha = 0.95, verbose=TRUE)
}

\arguments{
  \item{K}{
   An additive genomic relationship matrix, as provided for example by type=3 in \code{\link{ghap.kinship}}.
  }
  \item{ped}{
   A dataframe with columns "id", "sire" and "dam" containing pedigree data.
 }
  \item{include}{
   An optional vector of ids to be forced into the output (default = NULL). See details.
  }
  \item{depth}{
   A numeric value specifying the pedigree depth in number of generations to be used (default = 3). See details.
 }
 \item{alpha}{
  A numeric value between 0 and 1 specifying the weight of the genomic relationship matrix in the blend alpha*K + (1-alpha)*K. Default = 0.95.
 }
 \item{verbose}{
   A logical value specifying whether log messages should be printed (default = TRUE).
 }
}

\value{
  A matrix consisting of the inverse of H.
}

\details{

The pedigree is pruned to include only the number of generations specified by argument "depth". This prunning starts by seeding the genealogy with all individuals included in the genomic relationship matrix plus all individuals listed in the "include" argument. Then, the genealogy is increased by advancing one generation back at a time up to "depth". For example, if depth = 3, the genealogical tree will include all individuals listed in K and "include", plus their parents (depth = 1), grandparents (depth = 2) and great-grandparents (depth = 3). After the pedigree has been pruned, pedigree and genomic relationships are blended to form the inverse of H.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
A. I. Vazquez. Technical note: An R package for fitting generalized linear mixed models in animal breeding. J. Anim. Sci. 2010. 88, 497-504.
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
# ped <- read.table(file = "example.pedigree", header=T)
# 
# ### RUN ###
# 
# # This analysis emulates a scenario of
# # 100 individuals with genotypes and phenotypes
# # 200 individuals with only phenotypes
# # 400 individuals from pedigree with no data
# 
# # Subset 100 individuals from the pure1 population
# pure1 <- plink$id[which(plink$pop == "Pure1")]
# pure1 <- sample(x = pure1, size = 100)
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
# # Exclude pedigree records with missing sire
# ped <- ped[which(is.na(ped$sire) == F),]
# 
# # Make inverse of blended pedigree/genomic matrix
# ids <- unique(c(ped$id, ped$sire, ped$dam, colnames(K)))
# Hinv <- ghap.getHinv(K = K, ped = ped[,-1], include = ids)
# 
# # Run single-step GBLUP
# df$rep <- df$id
# model <- ghap.lmm(formula = pheno ~ 1 + (1|id) + (1|rep),
#                   data = df,
#                   covmat = list(id = Hinv, rep = NULL),
#                   invcov = T)

}
