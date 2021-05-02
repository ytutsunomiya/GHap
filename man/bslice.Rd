\name{ghap.bslice}
\alias{ghap.bslice}

\title{
 Get a slice of the plink object
}
\description{
  This function parses the binary genotype matrix from PLINK and returns the slice as an R matrix.
}
\usage{
 ghap.bslice(plink, ids, markers, index=FALSE,
             transposed=FALSE, sparse=TRUE, impute=FALSE,
             lookup=NULL, ncores=1, verbose=TRUE)
}

\arguments{
  \item{plink}{
   A GHap.plink object.
  }
  \item{ids}{
   A character or numeric vector indicating individuals to parse.
  }
  \item{markers}{
   A character or numeric vector indicating markers to parse.
  }
  \item{index}{
   A logical value specfying if values provided for ids and markers are indices (see details).
  }
  \item{transposed}{
   A logical value specfying if genotypes should be transposed. If FALSE (default), the matrix is returned as markers x individuals. Otherwise, the retrieved matrix is organized as individuals x markers.
  }
  \item{sparse}{
   A logical value specfying if the returned matrix should be formatted as a sparse matrix (default TRUE).
  }
  \item{impute}{
   A logical value specfying if missing genotypes should be replaced by 0 (i.e., A0/A0 genotypes, default = FALSE).
  }
  \item{lookup}{
   A character vector containing the look up table to decode the binary plink file. Only meaningful for developers.
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel parsing (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
 An R matrix with markers in rows and individuals in columns. If transposed = TRUE, then individuals are placed in rows and markers in columns.
}

\details{
 This function parses the PLINK binary file and returns an R matrix with the requested list of markers and individuals. The argument index allows the user to specify markers and individuals either by name (index = FALSE) or by indices as stored in the GHap.plink object (index = TRUE). The lookup argument also allows the user to provide a vector containing the decoding scheme to translate bits into integers in R. However, this argument should only be used by developers.
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}
