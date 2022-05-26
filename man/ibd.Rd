\name{ghap.ibd}
\alias{ghap.ibd}

\title{
Estimation of IBD sharing
}
\description{
This function estimates the same IBD statistics computed by plink's 'genome' option.
}
\usage{
ghap.ibd(object, pairlist, freq, mafcut=0.05,
         refsize=10000, batchsize=NULL,
         ncores=1, verbose=TRUE)
}
\arguments{
  \item{object}{
    A valid GHap object (phase or plink).
  }
  \item{pairlist}{
    A dataframe containing columns ID1 and ID2 specifying the pairs of individual ids to compare.
  }
  \item{freq}{
    A named numeric vector with (A1) allele frequencies computed in a reference sample.
  }
  \item{mafcut}{
   A numeric value specifying the minor allele frequency threshold for IBD calculations (default = 0.05).
  }
  \item{refsize}{
   A numeric value representing the reference sample size used in allele frequncy calculations. If not specified, a large reference sample is assumed (default = 10000)
  }
  \item{batchsize}{
   A numeric value controlling the number of variants to be processed at a time (default = nalleles/10).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}
\value{
A dataframe with columns:
POP1 = Population of individual 1\cr
ID1 = Name of individual 1\cr
POP2 = Population of individual 2\cr
ID2 = Name of individual 2\cr
IBS0 = Variant sites where ID1 and ID2 share no identical alleles\cr
IBS1 = Variant sites where ID1 and ID2 share 1 identical allele\cr
IBS2 = Variant sites where ID1 and ID2 share 2 identical alleles\cr
PERC = IBS2/(IBS0+IBS1+IBS2) [proportion of identical genotypes]\cr
DST = (IBS2 + 0.5*IBS1)/(IBS0+IBS1+IBS2) [proportion of shared alleles]\cr
Z0 = Proportion of the genome where ID1 and ID2 share no alleles IBD\cr
Z1 = Proportion of the genome where ID1 and ID2 share 1 allele IBD\cr
Z2 = Proportion of the genome where ID1 and ID2 share 2 alleles IBD\cr
PI_HAT = Z2 + 0.5*Z1 (proportion of IBD alleles shared between ID1 and ID2)
}
\details{
This function implements plink's method-of-moments for IBD estimation. Although not as efficient as plink's implementation, our function allows the user to restrict calculations to specific individual pairs, as well as ground all computations on allele frequencies obtained from a reference population. This is useful for routine pedigree confirmation, since a smaller set of indviduals and comparisons are typically targeted in these situations. The original --genome flag in plink not only performs all possible comparisons given a set of individuals, but also estimates allele frequencies on-the-fly, which may be unreliable if the number of individuals is small. We still recommend using plink for large problems, such as all pairwise comparisons from thousands of individuals, because it is more efficient. Nevertheless, we offer a more convenient alternative for the validation of smaller pedigrees in routine analyses where a lookup table of allele frequencies is available and maintained from a large reference population.
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
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
# # Load pedigree data
# ped <- read.table(file = "example.pedigree", header=T)
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
# # Compute A1 allele frequencies
# p <- ghap.freq(plink, type = "A1")
# 
# # Compute IBD statistics for individual 1
# pairlist <- data.frame(ID1 = pure1[1], ID2 = pure1[-1])
# ibd <- ghap.ibd(object = plink, pairlist = pairlist, freq = p,
#                 refsize = length(pure1))
# 
# # Predict relationships for individual 1
# # 1 = parent-offspring
# # 3 = other types of relationship
# # 4 = unrelated
# rel <- ghap.relfind(ibdpairs = ibd)
# table(rel$REL)
# 
# # Confirm with pedigree
# toprel <- rel$ID2[which(rel$REL == 1)]
# ped[which(ped$id \%in\% toprel & ped$dam == pure1[1]),]

}
