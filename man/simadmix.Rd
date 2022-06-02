\name{ghap.simadmix}
\alias{ghap.simadmix}

\title{
  Simulate individuals from specified admixture proportions
}
\description{
  Generation of simulated haplotypes based on a list of user-defined ancestral populations.
}
\usage{
  ghap.simadmix(object, n.individuals,
                n.generations, ancestors,
                proportions = NULL,
                alpha = NULL,
                out.file,
                only.active.markers = TRUE,
                ncores = 1,verbose = TRUE)
}

\arguments{

  \item{object}{
   A GHap.phase object.
  }
  \item{n.individuals}{
   Number of individuals to simulate.
  }
  \item{n.generations}{
   Number of generations past the admixture event.
  }
  \item{ancestors}{
   List of ancestral populations. Each ancestral population is given as a vector of ids of the ancestors belonging to that population.
  }
  \item{proportions}{
   A dataframe containing the ancestry proportions in the simulated individuals. The number of columns has to be equal to the number of ancestral populations defined in the ancestors list, and the number of row has to be equal to the number of simulated individuals. See argument 'alpha' if you want the function to generate random samples for admixture proportions.
  }
  \item{alpha}{
   A list with same size of the 'ancestors' list, with each value representing the vector of parameters to be used for sampling ancestry proportions from a Direchelet distribution.
  }
  \item{out.file}{
    Output file name.
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be used in simulations (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\details{
  Given a list of ancestral populations, this function simulates haplotypes of individuals descending from admixture among these populations after a defined number of generations. The ancestry proportions can be sampled from a Direchelet distribution with a vector of parameters alpha. For example, if three ancestral populations are considered and the parameters are set as alpha = list(pop1 = 1, pop2 = 1, pop3 = 1), the resulting simulated individuals will have a highly diverse configuration of ancestry proportions. On the other hand, by setting alpha = list(pop1 = 0, pop2 = 1, pop3 = 0) for example, all individuals will descend entirely from population 2 without admixture. The user can play with these values to fine tune the desired sampling distribution. Alternatively, exact proportions for each simulated individual can be set through the 'proportions' argument by providing a dataframe containing the desired ancestry values.
}

\value{
The function outputs the following files:

    \itemize{
      \item \strong{.samples}: space-delimited file without header containing two columns: Population and ID of simulated individuals (numbered from 1 to n). The first column is filled with "SIM". 
      \item \strong{.markers}: space-delimited file without header containing five columns: Chromosome, Marker, Position (in bp), Reference Allele (A0) and Alternative Allele (A1).
      \item \strong{.phase}: space-delimited file without header containing the phased genotype matrix of the simulated progeny. The dimension of the matrix is \emph{m x 2n}, where \emph{m} is the number of markers and \emph{n} is the number of simulated individuals (i.e., two columns per individual, representing the two phased chromosome alleles).
      \item \strong{.proportions}: space-delimited file with header containing the following columns: Population, ID and K columns of ancestry proportions, one for each ancestral population.
      \item \strong{.haplotypes}: space-delimited file with header containing ancestry tracks with the following columns: Population, ID, haplotype number, chromosome name, starting position, ending position, track size and ancestry of the segment. 
    }
    
  The user can then treat these files as a regular GHap input, or use them to build the Oxford HAPS/SAMPLES format for analysis with other software. The simulated data can be useful in the evaluation of accuracy of different algorithms designed for admixture analysis, as well as of other methods in the field of population genomics.
}

\author{
  Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{

# # Copy the example data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "phase",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load phase data
# phase <- ghap.loadphase("example")
# 
# # Make vectors of ancestors
# pure1 <- unique(phase$id[which(phase$pop == "Pure1")])
# pure2 <- unique(phase$id[which(phase$pop == "Pure2")])
# 
# # Simulate proportions
# ngroup <- 30
# pop1 <- c(rep(1, times = ngroup), # purebred from population 1
#           runif(n = ngroup, min = 0.1, max = 0.9), # admixed individuals
#           rep(0, times = ngroup)) # purebred from population 2
# pop2 <- 1-pop1
# prop <- data.frame(pop1,pop2)
# 
# # Simulate individuals
# set.seed(1988)
# ghap.simadmix(object = phase, n.individuals = nrow(prop),
#               n.generations = 10,
#               ancestors = list(pop1 = pure1, pop2 = pure2),
#               proportions = prop, out.file = "sim")
# ghap.compress(input.file = "sim", out.file = "sim")
# sim <- ghap.loadphase("sim")
# 
# # Unsupervised analysis with K = 2
# prototypes <- ghap.anctrain(object = sim, K = 2)
# hapadmix <- ghap.anctest(object = sim,
#                          prototypes = prototypes,
#                          test = unique(sim$id))
# anctracks <- ghap.ancsmooth(object = sim, admix = hapadmix)
# 
# # Load simulated ancestry proportions
# ancsim <- NULL
# ancsim$proportions2 <- read.table(file = "sim.proportions", header=T)
# ancsim$haplotypes <- read.table(file = "sim.haplotypes", header=T)
# 
# # Compare estimates with real values
# # Obs: the original populations had introgression from each other
# # Admixture seen in the estimates of purebreds reflect that introgression
# ghap.ancplot(ancsmooth = ancsim)
# ghap.ancplot(ancsmooth = anctracks)
# cor(anctracks$proportions2$K1, ancsim$proportions2$pop1)
# ghap.karyoplot(ancsmooth = ancsim, ids = sim$id[66])
# ghap.karyoplot(ancsmooth = anctracks, ids = sim$id[66])

}
