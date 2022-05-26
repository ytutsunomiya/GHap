\name{ghap.relfind}
\alias{ghap.relfind}

\title{
Find relatives in IBD estimates
}
\description{
This function infers relationships from IBD sharing.
}
\usage{
ghap.relfind(ibdpairs, v = 50,
             breakclass = FALSE, ncores=1)
}

\arguments{
  \item{ibdpairs}{
    A dataframe containing IBD estimates, such as those provided by the output of --genome in plink.
  }
  \item{v}{
    A hyperparameter controling the variance of the likelihood functions (smaller values lead to more variance).
  }
  \item{breakclass}{
    A logical value indicating if relationship types 2 and 3 should be reported using their subclasses (default = FALSE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
}

\value{
The function returns the original dataframe with the extra column 'REL', containing the relationship predictions.
}

\details{

The input dataframe must contain columns POP1 (or FID1), ID1 (or IID1), POP2 (or FID2), ID2 (or IID2), Z0, Z1, Z2 and PI_HAT. Columns Z0, Z1 and Z2 are the IBD sharing estimates representing the proportions of the genome where the two individuals being compared share exactly 0, 1 and 2 alleles identically by descent, respectively. The last column is Z2 + Z1/2, namely the proportion of the genome shared identically by descent.

This function implements a method based on composite likelihood scores to infer relationships based on the values of Z0, Z1, Z2 and PI_HAT. The predicted values are:

-1 = duplicates or monozygotic twins \cr
0 = parent-offspring with inbreeding or self-fertilization \cr
1 = parent-offspring \cr
2 = full-siblings \cr
3 = other types of relationships \cr
4 = unrelated \cr

Briefly, for each relationship type, the likelihood of each of the four IBD values is computed from a beta distribution with parameters a = m*v and b = (1-m)*v, where m is the expected value according to the relationship type and v is a hyperparameter controling the variance around the expected value (default v = 50). The composite likelihood is computed by summing the log-likelihoods for the four IBD values. The prediction is made by adopting the relationship type with the highest composite likelihood score. Details of the expected values of each relationship type is found in our vignette. This classification strategy is inspired by the method reported by Staples et al. (2014), albeit it is a different method. While Staples and collaborators infer relationships through Gaussian Kernel Density Estimation using only Z0 and Z1 values, our strategy uses the composite score formed by the sum of beta log-likelihoods for all IBD values.

An important detail is that relationship types 2 and 3 are in fact modelled through two and four different composite likelihoods, respectively, which are combined in order to achieve more robust and stable predictions. The user can choose to break down the subclasses via the argument 'breakclass = TRUE'. Reporting the subclasses is not default due to the reduced accuracy in distinguishing them. In addition, other relationship types not implicitly modelled in this version of the package may be confused with one of those subclasses. However, breaking these subclasses down may be useful for users seeking specific relationships in the data, as well as for those willing to perform pedigree simulations. The additional subclasses are:

2.1 = full-siblings\cr
2.2 = full-siblings from a self-fertilized parent (or related parents)\cr
3.1 = half-siblings, grandparent-grandchild or avuncular with inbreeding \cr
3.2 = half-siblings, grandparent-grandchild or avuncular \cr
3.3 = first-cousin or half-avuncular \cr
3.4 = half-cousin and distant relatives \cr

If the user chooses to run the analysis with 'breakclass = TRUE', beware that other types of cryptic relationships will be misclassified as pertaining to one of those four subclasses.
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
C. C. Chang et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015. 4, 7.

J. Staples et al. PRIMUS: Rapid Reconstruction of Pedigrees from Genome-wide Estimates of Identity by Descent. 2014. 95, 553-564.

S. Purcell et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am. J. Hum. Genet. 2007. 81, 559-575.
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
# # Copy metadata in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "meta",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# # Load plink data
# plink <- ghap.loadplink("example")
# 
# # Load pedigree data
# ped <- read.table(file = "example.pedigree", header=T)
# 
# ### RUN ###
# 
# # Subset individuals from the pure1 population
# pure1 <- plink$id[which(plink$pop == "Pure1")]
# plink <- ghap.subset(object = plink, ids = pure1, variants = plink$marker)
# 
# # Subset markers with MAF > 0.05
# freq <- ghap.freq(plink)
# mkr <- names(freq)[which(freq > 0.05)]
# plink <- ghap.subset(object = plink, ids = pure1, variants = mkr)
# 
# # Compute A1 allele frequencies
# p <- ghap.freq(plink, type = "A1")
# 
# # Compute IBD statistics for individual 1
# pairlist <- data.frame(ID1 = pure1[1], ID2 = pure1[-1])
# ibd <- ghap.ibd(object = plink, pairlist = pairlist, freq = p,
#                 refsize = length(pure1))
# 
# # Predict relationships for individual 1
# # 1 = parent-offspring
# # 3 = other types of relationship
# # 4 = unrelated
# rel <- ghap.relfind(ibdpairs = ibd)
# table(rel$REL)
# 
# # Confirm with pedigree
# toprel <- rel$ID2[which(rel$REL == 1)]
# ped[which(ped$id \%in\% toprel & ped$dam == pure1[1]),]

}
