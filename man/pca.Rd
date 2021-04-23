\name{ghap.pca}
\alias{ghap.pca}

\title{
Principal Components Analysis
}
\description{
PCA from a HapAllele-based kinship matrix.
}
\usage{
ghap.pca(haplo, kinship, npc = 2)
}

\arguments{
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{kinship}{
    A HapAllele-based kinship matrix, as supplied by \code{\link{ghap.kinship}}.
  }
  \item{npc}{
    Number of principal components to be retrieved (default = 2).
  }
}

\value{
The returned object is a list with items:
 \item{eigenvec}{
  A data.frame containing the principal components of the kinship matrix.
 }
 \item{eigenval}{
  Vector with eigenvalues of the kinship matrix.
 }
  \item{propvar}{
  Vector with the proportion of variance explained by each principal component.
 }
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
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
# # Exclude minor alleles and singletons
# hapstats <- ghap.hapstats(haplo, ncores = 1)
# haplo <- ghap.subsethaplo(haplo,ids=haplo$id, alleles = hapstats$TYPE \%in\% c("REGULAR","MAJOR"))
#
# # Compute Kinship matrix
# K <- ghap.kinship(haplo, ncores = 1)
#
#
# ### RUN ###
# 
# # PCA analysis
# pca <- ghap.pca(haplo,K)
# 
# # Plot
# plot(x=pca$eigenvec$PC1, y=pca$eigenvec$PC2, xlab="PC1", ylab="PC2", pch="")
# pop <- pca$eigenvec$POP
# pop.col <- as.numeric(as.factor(pop))
# pop <- sort(unique(pop))
# legend("bottomleft", legend = pop, col = 1:length(pop), pch = 1:length(pop), ncol = 3)
# points(x=pca$eigenvec$PC1, y=pca$eigenvec$PC2, pch = pop.col, col = pop.col, cex = 1.2)

}
