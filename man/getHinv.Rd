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
