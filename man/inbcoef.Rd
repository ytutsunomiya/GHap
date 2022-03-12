\name{ghap.inbcoef}
\alias{ghap.inbcoef}

\title{
Compute measures of inbreeding
}
\description{
This function computes genomic measures of inbreeding.
}
\usage{
ghap.inbcoef(object, freq, batchsize=NULL,
             only.active.samples=TRUE,
             only.active.variants=TRUE,
             ncores=1, verbose=TRUE)
}

\arguments{
  \item{object}{
    A valid GHap object (phase, haplo or plink).
  }
  \item{freq}{
    A named numeric vector providing allele frequencies.
  }
  \item{batchsize}{
   A numeric value controlling the number of variants to be processed at a time (default = nalleles/10).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.alleles}{
    A logical value specifying whether only active variants should be included in the output (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
The function returns a dataframe containing population name, id and inbreeding measures for each individual.
}

\details{

The inbreeding measures are computed as k + s*sum(f)/n, where k is a constant, s is a sign shifting scalar, f is a per-variant function and n is the number of variants (only non-monomorphic ones are considered). Three different measures of inbreeding are currently available:

Type = 1 (based on genomic relationship)\cr
k = -1\cr
s = 1\cr
f = ((m - 2p)^2)/(2p*(1-p))\cr

Type = 2 (excess homozygosity)\cr
k = 1\cr
s = -1\cr
f = m*(2-m)/(2p*(1-p))\cr

Type = 3 (correlation between uniting gametes)\cr
k = 0\cr
s = 1\cr
f = (m^2 - (1+2p) + 2p^2)/(2p*(1-p))\cr
}
\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}
\references{
J. Yang et al. GCTA: A Tool for Genome-wide Complex Trait Analysis. Am. J. Hum. Genet. 2011. 88:76–82.
}
\examples{
# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy plink data to the current working directory
# exfiles <- ghap.makefile(format = "plink")
# file.copy(from = exfiles, to = "./")
#
# ### RUN ###
# 
# # Load plink data
# plink <- ghap.loadplink(input.file = "human_hapmap3_chr2")
# 
# # Subset population
# pop <- "CEU"
# ids <- plink$id[which(plink$pop == pop)]
# mkrs <- plink$marker
# plink <- ghap.subset(object = plink, ids = ids, variants = mkrs)
# 
# # Calculate allele frequencies (has to be A1)
# p <- ghap.freq(plink, type = 'A1', ncores = 1)
# 
# # Compute inbreeding coefficients
# ibc <- ghap.inbcoef(object = plink, freq = p, ncores = 1)

}
