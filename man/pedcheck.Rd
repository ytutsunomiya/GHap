\name{ghap.pedcheck}
\alias{ghap.pedcheck}

\title{
Summary statistics for pedigree
}
\description{
This function summarizes pedigree data and calculates inbreeding coefficients and equivalent complete generations.
}
\usage{
ghap.pedcheck(ped, depth.n.f = FALSE)
}

\arguments{
 \item{ped}{
   A dataframe with columns "id", "sire" and "dam" containing pedigree data.
 }
 \item{depth.n.f}{
   A logical indicating if equivalent complete generations (depth) and inbreeding coefficients (f) should be calculated.
 }
}

\value{
  A list containing two data frames: stats, which includes the pedigree summary; and ped, consisting of the original pedigree. If depth.n.f = TRUE, three columns are added to the ped data frame: gen (generation number), f (inbreding coefficient) and ecg (equivalent complete generations). The generation number and the inbreeding coefficient are computed with the help of the pedigreemm package.
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
A. I. Vazquez. Technical note: An R package for fitting generalized linear mixed models in animal breeding. J. Anim. Sci. 2010. 88, 497-504.
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy metadata in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "meta",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load pedigree data
# ped <- read.table(file = "example.pedigree", header=T)
# 
# ### RUN ###
# 
# # Descriptive statistics for the pedigree
# pedstat <- ghap.pedcheck(ped[,-1])
# print(pedstat$stats)
# 
# # Retrieve inbreeding and pedigree depth
# pedstat <- ghap.pedcheck(ped[,-1], depth.n.f = TRUE)
# print(pedstat$ped)
# hist(pedstat$ped$f)
# hist(pedstat$ped$ecg)

}
