\name{ghap.fast2phase}
\alias{ghap.fast2phase}

\title{
  Convert fastPHASE data into the GHap phase format
}
\description{
  This function takes phased genotype data in fastPHASE format and converts them into a GHap plan and compressed binary format.
}
\usage{
  ghap.fast2phase(input.files = NULL, switchout.files = NULL,
                  map.files = NULL, fam.file = NULL,
                  out.file = NULL, overwrite = FALSE,
                  ncores = 1, verbose = TRUE)
}

\arguments{

  If all input files share the same prefix, the user can use the following shortcut options:

  \item{input.files}{
    Character vector with the list of prefixes for input files.
  }
  \item{out.file}{
    Character value for the output file name.
  }

  The user can also opt to point to input files separately:
  
  \item{switchout.files}{
    Character vector containing the list of fastPHASE files.
  }
  \item{map.files}{
    Character vector containing the list of map files.
  }
  \item{fam.file}{
    Name of the file containing the population and individual ids.
  }
  
  The function avoids overwritting existing files by default. To change this behavior, please use:
  
  \item{overwrite}{
  	A logical value controling if existing files with the same name as the selected output should be overwritten (default = FALSE).
  }
  
  To turn conversion progress-tracking on or off or set the number of cores please use:

  \item{ncores}{
  	A numeric value specfying the number of cores (default = 1).
  }

  \item{verbose}{
  	A logical value specfying whether log messages should be printed (default = TRUE).
  }
  


}

\details{
Currently this function handles _switch.out files from fastPHASE v1.4.0 or later. The map files should contain the following 5 space-delimited columns: chromosome, marker, position, allele 0 and allele 1 (no header). The fam file can have an arbitrary number of columns (with no header), but by default only the first two are read and should be space-delimited with the following data: population and individual name.
}

\references{
Scheet, P., Stephens, M., 2006. A Fast and Flexible Statistical Model for Large-Scale Population Genotype Data: Applications to Inferring Missing Genotypes and Haplotypic Phase. Am. J. Hum. Genet. 78, 629-644. https://doi.org/10.1086/502802.
}

\author{
Mario Barbato <mario.barbato@unicatt.it>
}

\seealso{
  \code{\link{ghap.compress}}, \code{\link{ghap.loadphase}}, \code{\link{ghap.oxford2phase}}, \code{\link{ghap.vcf2phase}}
}
