\name{ghap.makefile}
\alias{ghap.makefile}

\title{
 Create example input files
}
\description{
 Create example files to test the package.
}
\usage{
 ghap.makefile(dataset = "example",
               format = "phase",
               verbose = TRUE) 
}
\arguments{
  \item{dataset}{
    A character value specfying the name of the dataset.
  }
  \item{format}{
    A character value specfying the format of the dataset.
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}
\details{
This function downloads example files to the R temporary directory (requires internet connection). The default dataset comprises the following group of files:

\emph{example.phaseb} \cr
\emph{example.markers} \cr
\emph{example.samples} \cr

For details about the format of these files, see \code{\link{ghap.compress}}. The dataset was simulated using the QMSim v1.10 software and contains 450 individuals genotyped for 15,000 markers. These markers were randomly distributed along 10 chromosomes of 100 Mbp each (i.e., 1,500 markers per chromosome). Two divergent lineages were created, namely 'Pure1' (n = 300) and 'Pure2' (n = 100), and gene flow between these two lineages was allowed to produce low levels of admixture. An additional set of 50 crossbred individuals was also included. The same dataset is available in the following formats: raw (equal to the phase format, except that the genotype matrix is not compressed), vcf and oxford. By using format = 'meta', metadata including pedigree and phenotypes for the 'Pure1' population can be downloaded. The pedigree contains 700 records, spanning 5 generations. The records in the phenotypes file are unbalanced repeated measurements (1 to 5 records per individual, with an average of 3) of a trait with heritability of 0.4, repeatability of 0.2, and a major QTL located on chromosome 3.

Since version 2.1.0, GHap maintains additional example files in its github repository (\url{https://github.com/ytutsunomiya/GHap}). In order to see which files are available, please see \code{\link{ghap.exfiles}}.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com> \cr
Marco Milanesi <marco.milanesi.mm@gmail.com>
}

\examples{
# # See list of example files
# exlist <- ghap.exfiles()
# View(exlist)
# 
# # Copy example data in phase format
# exfiles <- ghap.makefile()
# file.copy(from = exfiles, to = "./")
# 
# # Copy example data in plink format
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "plink", verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Copy phenotypes and pedigree data
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "meta", verbose = TRUE)
# file.copy(from = exfiles, to = "./")
}

\references{
M. Sargolzaei and F. S. Schenkel. QMSim: A large-scale genome simulator for livestock. Bioinformatics. 2009. 25, 680–681.
}
