\name{ghap.slice}
\alias{ghap.slice}

\title{
  Get a slice of a GHap object
}
\description{
  This function parses a binary PLINK, phased or HapGenotypes matrix and returns the slice as an R matrix.
}
\usage{
ghap.slice(object, ids, variants, index=FALSE,
           transposed=FALSE, sparse=TRUE, unphase=FALSE,
           impute=FALSE, lookup=NULL, ncores=1, verbose=TRUE)
}

\arguments{
  \item{object}{
    A valid GHap object (phase, haplo or plink).
  }
  \item{ids}{
    A character or numeric vector indicating individuals to parse.
  }
  \item{variants}{
    A character or numeric vector indicating variants to parse. If a "GHap.haplo" object is provided, the vector must be numeric.
  }
  \item{index}{
    A logical value specfying if values provided for ids and variants are indices (see details).
  }
  \item{transposed}{
   A logical value specfying if genotypes should be transposed. If FALSE (default), the matrix is returned as variants by individuals. Otherwise, the retrieved matrix is organized as individuals by variants.
  }
  \item{sparse}{
   A logical value specfying if the returned matrix should be formatted as a sparse matrix (default TRUE).
  }
  \item{unphase}{
   A logical value specfying if phased genotypes should be retrieved as unphased allele counts. Only meaningful for "GHap.phase" objects.
  }
  \item{impute}{
   A logical value specfying if missing genotypes should be replaced by 0 (i.e., A0/A0 genotypes, default = FALSE). Only meaningful for "GHap.plink" objects.
  }
  \item{lookup}{
    A character vector containing the look up table to decode the binary file. Only meaningful for developers.
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
  An R matrix with variants in rows and individuals in columns (this is inverted if transposed = TRUE).
}

\details{
  This function parses the binary input file and returns an R matrix with the requested list of variants (markers or alleles) and individuals. The argument index allows the user to specify individuals either by name (index = FALSE) or by indices as stored in the GHap object (index = TRUE). In the case of "GHap.haplo" objects, HapAlleles can only by parsed via indices. The lookup argument also allows the user to provide a vector containing the decoding scheme to translate bits into integers in R. However, this argument should only be used by developers.
}

\author{
  Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}