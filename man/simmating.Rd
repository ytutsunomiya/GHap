\name{ghap.simmating}
\alias{ghap.simmating}

\title{
  Simulate individuals from specified matings
}
\description{
  Generation of simulated haplotypes based on in silico matings.
}
\usage{
  ghap.simmating(object, n.individuals = 1,
                 parent1 = NULL, parent2 = NULL,
                 model = "proportional", out.file,
                 only.active.markers = TRUE,
                 ncores = 1, verbose = TRUE)
}

\arguments{

  \item{object}{
   A GHap.phase object.
  }
  \item{n.individuals}{
   Number of individuals to simulate.
  }
  \item{parent1}{
   Vector containing the ids of candidate sires. When a character vector is provided, all candidates are considered equally likely to be selected. Alternatively, a named numeric vector of probabilities (with sire ids used as names) can be used in order to make specific sires more likely to be selected.
  }
  \item{parent2}{
   Vector containing the ids of candidate dams. When a character vector is provided, all candidates are considered equally likely to be selected. Alternatively, a named numeric vector of probabilities (with dam ids used as names) can be used in order to make specific dams more likely to be selected.
  }
  \item{model}{
   The model used for sampling the number of recombinations per chromosome (default = "proportional"). See details for more options.
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
  Given a list of candidate parents, this function samples sire i with probability p_i and dam j with probability p_j (these probabilities should be provided by the user, otherwise matings will occur at random). Once sire and dam are sampled, gametes are simulated by creating recombinant parental haplotypes. The progeny is then obtained by uniting the simulated gametes.
  
The default model ("proportional") assumes that the number of recombinations per meiosis across all chromosomes follows a Poisson distribution with mean equal to nchr (the number of chromosome pairs), such that the number of recombinations for a given chromosome is sampled from a Poisson distribution with mean prop*nchr, where prop is the proportion of the genome size covered by that chromosome. Therefore, this option takes chromosome size in consideration while sampling the number of recombination events. In option "uniform", the number of recombinations is sampled from a Poisson distribution with mean 1 for each chromosome instead. Alternatively, the model argument can also take a named vector of chromosome-specific recombination rates in cM/Mb. In this option, the number of recombinations for a given chromosome is sampled from a Poisson distribution with mean chrsize*chrrate/100, where chrsize is the size of the chromosome in Mb and chrrate is the recombination rate in cM/Mb. A last option is offered where the model argument is a named vector of marker-specific recombination rates in cM/Mb. The number of recombinations for a given chromosome is also sampled from a Poisson distribution with mean chrsize*chrrate/100, but with chrrate calculated as the average across markers within the same chromosome. In this model, instead of placing the recombination breakpoint randomly within a chromosome for each gamete, marker-specific recombination rates are taken into account, making regions in the chromosome with higher recombination rates more susceptible to breaks.
}

\value{
The function outputs the following files:

    \itemize{
      \item \strong{.samples}: space-delimited file without header containing two columns: Population and ID of simulated individuals (numbered from 1 to n). The first column is filled with "SIM". 
      \item \strong{.markers}: space-delimited file without header containing five columns: Chromosome, Marker, Position (in bp), Reference Allele (A0) and Alternative Allele (A1).
      \item \strong{.phase}: space-delimited file without header containing the phased genotype matrix of the simulated progeny. The dimension of the matrix is \emph{m x 2n}, where \emph{m} is the number of markers and \emph{n} is the number of simulated individuals (i.e., two columns per individual, representing the two phased chromosome alleles).
      \item \strong{.pedrigree}: space-delimited file without header containing three columns: individual number, sire ID and dam ID. 
    }
    
  The user can then treat these files as a regular GHap input, or use them to build the Oxford HAPS/SAMPLES format for analysis with other software. The simulated data can be useful in the evaluation of mating plans, since the in silico progeny can help in the characterization of the expected distribution of genomic inbreeding coefficients, ancestry proportions and EBVs in the progeny. Since the function does not check for sex differences between parents, the user can perform simulations of self-fertilization (relevant for plant breeders) and same sex matings.
  
}

\author{
  Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy the example data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "phase",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load phase data
# phase <- ghap.loadphase("example")
# 
# ### RUN ###
# 
# # Simulation using only two specific parents
# parent1 <- phase$id[1]
# parent2 <- phase$id[3]
# ghap.simmating(phase, n.individuals = 100,
#                parent1 = parent1, parent2 = parent2,
#                out.file = "sim1", ncores = 1)
# 
# # Simulation using candidates with unequal probabilities
# parent1 <- c(0.5,0.25,0.25)
# names(parent1) <- phase$id[c(1,3,5)]
# parent2 <- c(0.7,0.2,0.1)
# names(parent2) <- phase$id[c(7,9,11)]
# ghap.simmating(phase, n.individuals = 100,
#                parent1 = parent1, parent2 = parent2,
#                out.file = "sim2", ncores = 1)

}
