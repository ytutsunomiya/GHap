\name{ghap.hapstats}
\alias{ghap.hapstats}

\title{
Haplotype allele statistics
}
\description{
Summary statistics for HapAlleles.
}
\usage{
ghap.hapstats(haplo, alpha=c(1,1), batchsize=NULL, only.active.samples=TRUE,
              only.active.alleles=TRUE, ncores=1, verbose=TRUE)
}

\arguments{
  \item{haplo}{
    A GHap.haplo object.
  }
  \item{alpha}{
   A numeric vector of size 2 specifying the shrinkage parameters for the expected-to-observed homozygotes ratio. Default is c(1,1).
  }
  \item{batchsize}{
   A numeric value controlling the number of HapAlleles to be processed at a time (default = nalleles/10).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.alleles}{
    A logical value specifying whether only active haplotype alleles should be included in the output (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
 A data frame with columns:
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
 \item{N}{
  Number of observations for the haplotype.
 }
 \item{FREQ}{
  Haplotype frequency.
 }
 \item{O.HOM}{
  Observed number of homozygotes.
 }
 \item{O.HET}{
  Observed number of heterozygotes.
 }
 \item{E.HOM}{
  Expected number of homozygotes.
 }
 \item{RATIO}{
  Shrinkage expected-to-observed ratio for the number of homozygotes.
 }
 \item{BIN.logP}{
  log10(1/P) or -log10(P) for Hardy-Weinberg equilibrium assuming number of homozygotes follows a Binomial distribution.
 }
 \item{POI.logP}{
  log10(1/P) or -log10(P) for Hardy-Weinberg equilibrium assuming number of homozygotes follows a Poisson distribution.
 }
 \item{TYPE}{
  Category of the HapAllele: "SINGLETON" = single allele of its block; "ABSENT" = the frequency of the allele is 0; "MINOR" = the least frequent allele of its block (in the case of ties, only the first allele is marked); "MAJOR" = the most frequent allele of its block (ties are also resolved by marking the first allele); "REGULAR" = the allele does not fall in any of the previous categories. Categories "SINGLETON", "MINOR" and "MAJOR" only apply for blocks where frequencies sum to 1.
 }
 
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com> \cr
Marco Milanesi <marco.milanesi.mm@gmail.com>
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
#
# ### RUN ###
#
# Subset
# ids <- haplo$id[which(haplo$pop == "CEU")]
# haplo <- ghap.subsethaplo(haplo = haplo, ids = ids, alleles = haplo$allele.in)
#
# #Compute haplotype statistics
# hapstats <- ghap.hapstats(haplo, ncores = 1)

}
