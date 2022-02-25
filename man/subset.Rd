\name{ghap.subset}
\alias{ghap.subset}

\title{
 Subset GHap objects
}

\description{
 This function takes a list of variants and individuals and subsets a GHap object.
}

\usage{
ghap.subset(object, ids, variants,
            index=FALSE, verbose=TRUE)
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
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
  A GHap object of the same type as the one used in the object argument.
}

\details{
  This function sets to FALSE (i.e., inactivates) all individuals and variants not included in the provided arguments. This procedure avoids expensive subsetting operations by simply flagging which variants and individuals should be used in downstream analyses. The argument index allows the user to specify individuals either by name (index = FALSE) or by indices as stored in the GHap object (index = TRUE). In the case of "GHap.haplo" objects, HapAlleles can only by parsed via indices.
}

\author{
  Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}
