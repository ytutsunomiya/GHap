\name{ghap.haplotyping}
\alias{ghap.haplotyping}

\title{
Haplotype genotypes
}
\description{
Generate matrix of HapGenotypes for user-defined blocks.
}
\usage{
ghap.haplotyping(phase, blocks, outfile, freq=c(0,1), drop.minor=FALSE,
                 only.active.samples=TRUE, only.active.markers=TRUE,
                 batchsize=NULL, binary=TRUE, ncores=1, verbose=TRUE)

}

\arguments{
  \item{phase}{
   A GHap.phase object.
  }
  \item{blocks}{
   A data frame containing block boundaries, such as supplied by the \code{\link{ghap.blockgen}} function.
  }
  \item{outfile}{
   A character value specifying the name for the output files.
  }
  \item{freq}{
   A numeric vector of length 2 specifying the range of haplotype allele frequency to be included in the output. Default is c(0,1), which includes all alleles.
  }
  \item{drop.minor}{
   A logical value specfying whether the minor allele should be excluded from the output (default = FALSE).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be used for haplotyping (default = TRUE).
  }
  \item{batchsize}{
   A numeric value controlling the number of haplotype blocks to be processed and written to output at a time (default = nblocks/10).
  }
  \item{binary}{
   A logical value specfying whether the output file should be binary (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{

  The function outputs three files with suffix:
  \itemize{
    \item \strong{.hapsamples}: space-delimited file without header containing two columns: Population and Individual ID.
    \item \strong{.hapalleles}: space-delimited file without header containing five columns: Block Name, Chromosome, Start and End Position (in bp), and HapAllele.
    \item \strong{.hapgenotypes}: if binary = FALSE, a space-delimited file without header containing the HapGenotype matrix (coded as 0, 1 or 2 copies of the HapAllele). The dimension of the matrix is \emph{m x n}, where \emph{m} is the number of HapAlleles and \emph{n} is the number of individuals.
    \item \strong{.hapgenotypesb}: if binary = TRUE (default), the same matrix as described above compressed into bits. For seamless compatibility with softwares that use PLINK binary files, the compression is performed using the SNP-major bed format.
}
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>

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
# 
# ### RUN ###
# 
# # Haplotyping
# ghap.haplotyping(phase = phase, blocks = blocks, outfile = "human",
#                  binary = T, ncores = 1)

}
