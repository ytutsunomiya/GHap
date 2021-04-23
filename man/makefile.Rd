\name{ghap.makefile}
\alias{ghap.makefile}

\title{
 Create example input files
}
\description{
 Create example files to test the package.
}
\usage{
 ghap.makefile(verbose = TRUE) 
}
\arguments{
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}
\details{
This function downloads the following example files to the R temporary directory (requires internet connection):

\emph{human.phase} \cr
\emph{human.markers} \cr
\emph{human.samples} \cr

For details about the format of these files, see \code{\link{ghap.compress}}. The dataset was extracted from the reference phased data available at the IMPUTE2 (Howie et al., 2009) software website (\url{https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference}).

The genotypes derive from the International HapMap Project Phase 3 (The International HapMap 3 Consortium, 2010), and comprise 1,011 subjects (from 11 populations) and 20,000 SNPs (randomly sampled from chromosome 2) mapped to the NCBI build 36 (hg18) assembly. The source link is (\url{https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference})

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com> \cr
Marco Milanesi <marco.milanesi.mm@gmail.com>
}

\examples{
# # Copy the example data in the current working directory
# exfiles <- ghap.makefile()
# file.copy(from = exfiles, to = "./")
}

\references{
B. N. Howie, P. Donnelly, and J. Marchini. A flexible and accurate genotype imputation method for the next generation of genome-wide association studies. PLOS Genet. 2009. 5, e1000529.

The International HapMap 3 Consortium. Integrating common and rare genetic variation in diverse human populations. Nature. 2010. 467, 52-58.
}
