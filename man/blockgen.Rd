\name{ghap.blockgen}
\alias{ghap.blockgen}

\title{
 Haplotype block generator
}
\description{
 This function generates HapBlocks based on sliding windows. The window and the step size can be specified in markers or kbp. For each window, block coordinates are generated.
}
\usage{
 ghap.blockgen(object, windowsize = 10, slide = 5,
               unit = "marker", nsnp = 2)
}

\arguments{
  \item{object}{
   A GHap.phase object.
  }
  \item{windowsize}{
   A numeric value for the size of the window (default = 10).
  }
  \item{slide}{
   A numeric value for the step size (default = 5).
  }
  \item{unit}{
   A character value for the size unit used for the window and the step. It can be either "marker" or "kbp" (default = "marker").
  }
  \item{nsnp}{
   A numeric value for the minimum number of markers per block.
  }
}

\value{
 A data frame with columns:
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
 \item{SIZE}{
  Haplotype size.
 }
 \item{NSNP}{
  Number of marker.
 }
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
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
# ### RUN ###
# 
# # Generate blocks of 5 markers sliding 5 markers at a time
# blocks.mkr <- ghap.blockgen(phase, windowsize = 5,
#                             slide = 5, unit = "marker")
# 
# # Generate blocks of 100 kbp sliding 100 kbp at a time
# blocks.kb <- ghap.blockgen(phase, windowsize = 100,
#                            slide = 100, unit = "kbp")


}
