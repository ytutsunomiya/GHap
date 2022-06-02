\name{ghap.froh}
\alias{ghap.froh}

\title{
 Calculation of genomic inbreeding (FROH)
}
\description{
  Given runs of homozygosity (ROH) obtained with the \code{\link{ghap.roh}} function, this function computes the proportion of the genome covered by ROHs (FROH) of certain lengths.
}
\usage{
 ghap.froh(object, roh, rohsizes = c(1, 2, 4, 8, 16),
           only.active.markers = TRUE, ncores = 1)
}

\arguments{
  \item{object}{
   A valid GHap object (phase or plink).
  }
  \item{roh}{
   A data frame containing runs of homozygosity, such as supplied by the \code{\link{ghap.roh}} function.
  }
  \item{rohsizes}{
   A numeric vector providing the minimum ROH length (in Mbp) to use in the calculation of ROH (default is 1, 2, 4, 8 and 16). 
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be used in the calculation of genome size (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
}

\value{
The function returns a dataframe the following columns:
 \item{POP}{
  Original population label.
 }
 \item{ID}{
  Individual name.
 }
 \item{FROH...}{
  A number of columns giving FROH calculated over runs of length greater than each of the sizes informed by the rohsize argument. Default values will return FROH1, FROH2, FROH4, FROH8 and FROH16.
 }

}

\details{
 This function takes runs of homozygosity obtained with \code{\link{ghap.roh}} and returns estimates of genomic inbreeding (FROH). The user can specify the minimum ROH length considered in the calculation using the rohsize argument. A vector of values will cause the function to add an extra column for each specified ROH size. Since the average size (measured in Morgans) of identical-by-descent segments after g generations of the inbreeding event is 1/2g, the default lengths 1, 2, 4, 8 adn 16 are proxies for inbreeding that occurred 50, 25, 13 6 and 3 generations in the past, respectively (assuming an average recombination rate of 1 Mbp ~ cM).
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy plink data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "plink",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load plink data
# plink <- ghap.loadplink("example")
# 
# ### RUN ###
# 
# # Subset pure1 population
# pure1 <- plink$id[which(plink$pop == "Pure1")]
# plink <- ghap.subset(object = plink, ids = pure1, variants = plink$marker)
# 
# # ROH via the 'naive' method
# roh1 <- ghap.roh(plink, method = "naive")
# froh1 <- ghap.froh(plink, roh1)
# 
# # ROH via the 'hmm' method
# freq <- ghap.freq(plink, type = 'A1')
# inbcoef <- froh1$FROH1; names(inbcoef) <- froh1$ID
# roh2 <- ghap.roh(plink, method = "hmm", freq = freq,
#                 inbcoef = inbcoef)
# froh2 <- ghap.froh(plink, roh2)
#
# # Method 'hmm' using Fhat3 as starting values
# inbcoef <- ibc$Fhat3; names(inbcoef) <- ibc$ID
# inbcoef[which(inbcoef < 0)] <- 0.01
# roh3 <- ghap.roh(plink, method = "hmm", freq = freq,
#                  inbcoef = inbcoef)
# froh3 <- ghap.froh(plink, roh3)

}

\seealso{
 \code{\link{ghap.roh}}
}
