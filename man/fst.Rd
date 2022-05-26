
   
\name{ghap.fst}
\alias{ghap.fst}

\title{
 Haplotype-based Fst
}

\description{
Multi-alleleic Fst computed using block summary statistics generated from \code{\link{ghap.blockstats}}.
}

\usage{
 ghap.fst(blockstats.pop1,
          blockstats.pop2,
          blockstats.tot)
}

\arguments{
  \item{blockstats.pop1}{
   A data.frame containing block statistics computed on population 1.
  }
  \item{blockstats.pop2}{
   A data.frame containing block statistics computed on population 2.
  }
  \item{blockstats.tot}{
   A data.frame containing block statistics computed on population 1 + population 2.
  }
}

\value{
The function returns a data.frame with the following columns:
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
  \item{EXP.H.pop1}{
  Expected heterozygosity in population 1.
 }
  \item{EXP.H.pop2}{
  Expected heterozygosity in population 2.
 }
  \item{EXP.H.tot}{
  Expected heterozygosity in the total population.
 }
 \item{FST}{
  Fst value.
 }

}

\details{
 This function calculates Fst (Nei, 1973) based on the formula for multi-allelic markers:
 
 \deqn{Fst = (Ht - Hs) / Ht}
 
 where \emph{Ht} is the total gene diversity (i.e., expected heterozygosity in the population) and \emph{Hs} is the subpopulation gene diversity (i.e., the average expected heterozygosity in the subpopulations).
  
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com> \cr
Marco Milanesi <marco.milanesi.mm@gmail.com>
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy phase data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "phase",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load data
# phase <- ghap.loadphase("example")
# 
# # Generate blocks of 5 markers
# blocks <- ghap.blockgen(phase, windowsize = 5,
#                         slide = 5, unit = "marker")
# 
# # Haplotyping
# ghap.haplotyping(phase = phase, blocks = blocks, outfile = "example",
#                  binary = T, ncores = 1)
# 
# # Load haplotype genotypes using prefix
# haplo <- ghap.loadhaplo("example")
# 
# ### RUN ###
# 
# # Compute block statistics for population 1
# ids <- which(haplo$pop == "Pure1")
# haplo <- ghap.subset(haplo, ids = ids,
#                      variants = haplo$allele.in,
#                      index = TRUE)
# hapstats1 <- ghap.hapstats(haplo)
# blockstats1 <- ghap.blockstats(hapstats1)
# 
# # Compute block statistics for population 2
# ids <- which(haplo$pop == "Pure2")
# haplo <- ghap.subset(haplo, ids = ids,
#                      variants = haplo$allele.in,
#                      index = TRUE)
# hapstats2 <- ghap.hapstats(haplo)
# blockstats2 <- ghap.blockstats(hapstats2)
# 
# # Compute block statistics for combined populations
# ids <- which(haplo$pop %in% c("Pure1","Pure2"))
# haplo <- ghap.subset(haplo, ids = ids,
#                      variants = haplo$allele.in,
#                      index = TRUE)
# hapstats12 <- ghap.hapstats(haplo)
# blockstats12 <- ghap.blockstats(hapstats12)
# 
# # Compute FST
# fst <- ghap.fst(blockstats1, blockstats2, blockstats12)
# ghap.manhattan(data = fst, chr = "CHR", bp = "BP1",
#                y = "FST", type = "h")

}

\references{
M. Nei. Analysis of Gene Diversity in Subdivided Populations. PNAS. 1973. 70, 3321-3323.
}
