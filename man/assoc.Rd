\name{ghap.assoc}
\alias{ghap.assoc}

\title{
Genome-wide association analysis
}
\description{
This function performs phenotype-genotype association analysis based on mixed models.
}
\usage{
ghap.assoc(object, formula, data, covmat,
           ngamma = 100, nlambda = 1000,
           recalibrate = 0.01,
           only.active.variants=TRUE,
           tol = 1e-12, ncores=1,
           verbose=TRUE, ...)
}

\arguments{
  \item{object}{
    A valid GHap object (phase, haplo or plink).
  }
  \item{formula}{
   Formula describing the model. The synthax is consistent with lme4. The response is declared first, followed by the ~ operator. Predictors are then separated by + operators. Currently only random intercepts are supported, which are distinguished from fixed effects by the notation (1|x). If multiple random effects are specified, the first declared in the formula will be assumed to be the genetic (polygenic) effects.
  }
  \item{data}{
   A dataframe containing the data.
  }
  \item{covmat}{
   A list of covariance matrices for each group of random effects. If a matrix is not defined for a given group, an identity matrix will be used. Inverse covariance matrices can also be provided, as long as argument invcov = TRUE is used.
 }
  \item{ngamma}{
   A numeric value for the number of variants to be used in the estimation of the gamma factor (default = 100). The grammar-gamma approximation is turned off if ngamma = 0.
 }
 \item{nlambda}{
  A numeric value for the number of variants to be used in the estimation of the inflation factor (default = 1000).
 }
 \item{recalibrate}{
  A numeric value for the proportion of top scoring variants to re-analyze without the grammar-gamma approximation (default = 0.01). Not relevant if ngamma = 0.
 }
 \item{only.active.variants}{
  A logical value specifying whether only active variants should be included in the output (default = TRUE).
 }
 \item{tol}{
  A numeric value specifying the scalar to add to the diagonal of the phenotypic (co)variance matrix if it is not inversible (default = 1e-12).
 }
 \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
 }
 \item{verbose}{
   A logical value specifying whether log messages should be printed (default = TRUE).
 }
  \item{...}{
   Additional arguments to be passed to the \code{\link{ghap.lmm}} function.
 }
}

\value{
The function returns a data frame with results from the genome-wide association analysis. If a GHap.haplo object is used, the first columns of the data frame will be:
 \item{CHR}{
  Chromosome name.
 }
 \item{BLOCK}{
  Block alias.
 }
 \item{BP1}{
  Block start position.
 }
 \item{BP2}{
  Block end position.
 }
 For GHap.phase and GHap.plink objects, the first columns will be:
 \item{CHR}{
  Chromosome name.
 }
 \item{MARKER}{
  Block start position.
 }
 \item{BP}{
  Block end position.
 }
 The remaining columns of the data frame will be equal for any class of GHap objects:
  \item{ALLELE}{
  Identity of the counted (A1 or haplotype) allele.
 }
 \item{FREQ}{
  Frequency of the allele.
 }
 \item{N}{
  Number of non-missing observations.
 }
 \item{BETA}{
  Estimated allele substitution effect.
 }
 \item{SE}{
  Standard error for the allele substitution effect.
 }
 \item{CHISQ.EXP}{
  Expected values for the test statistics.
 }
 \item{CHISQ.OBS}{
  Observed value for the test statistics.
 }
 \item{CHISQ.GC}{
  Test statistics scaled by the inflation factor (Genomic Control). Inflation is computed through regression of observed quantiles onto expected quantiles. In order to avoid overestimation by variants rejecting the null hypothesis, a random sample of variants (with size controled via the nlambda argument) is taken within three standard deviations from the mean of the distribution of test statistics.
 }
 \item{LOGP}{
  log10(1/P) or -log10(P) for the allele substitution effect.
 }
 \item{LOGP.GC}{
  log10(1/P) or -log10(P) for the allele substitution effect (scaled by the inflation factor).
 }

}

\details{
This function uses mixed models and the grammar-gamma approximation for fast genome-wide association analysis. Since mixed models are fit using the \code{\link{ghap.lmm}} function, the association analysis can be performed using more flexible models than those offered by alternative software, including the use of repeated measurements and other random effects apart from polygenic effects.
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
N. Amin et al. A Genomic Background Based Method for Association Analysis in Related Individuals. PLoS ONE. 2007. 2:e1274.

Y. Da. Multi-allelic haplotype model based on genetic partition for genomic prediction and variance component estimation using SNP markers. BMC Genet. 2015. 16:144.

B. Devlin and K. Roeder. Genomic control for association studies. Biometrics. 1999. 55:997-1004.

C. C. Ekine et al. Why breeding values estimated using familial data should not be used for genome-wide association studies. G3. 2014. 4:341-347.

L. Jiang et al. A resource-efficient tool for mixed model association analysis of large-scale data. Nat. Genet. 2019. 51:1749-1755.

J. Listgarten et al. Improved linear mixed models for genome-wide association studies. Nat. Methods. 2012. 9:525-526.

G. R. Svishcheva et al. Rapid variance components-based method for whole-genome association analysis. Nat Genet. 2012. 44:1166-1170.

J. Yang et al. Advantages and pitfalls in the application of mixed-model association methods. Nat. Genet. 2014. 46: 100-106.
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
# # Load phenotype data
# df <- read.table(file = "example.phenotypes", header=T)
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
# # Compute genomic relationship matrix
# # Induce sparsity to help with matrix inversion
# K <- ghap.kinship(plink, sparsity = 0.01)
# 
# # Perform GWAS on repeated measures
# # Use grammar-gama approximation
# # Recalibrate top 1 percent variants
# df$rep <- df$id
# gwas1 <- ghap.assoc(object = plink,
#                     formula = pheno ~ 1 + (1|id) + (1|rep),
#                     data = df,
#                     covmat = list(id = K, rep = NULL),
#                     ngamma = 100, nlambda = 1000, recalibrate = 0.01)
# ghap.manhattan(data = gwas1, chr = "CHR", bp = "BP", y = "LOGP")
# 
# # GWAS with no approximaion (slow)
# gwas2 <- ghap.assoc(object = plink,
#                     formula = pheno ~ 1 + (1|id) + (1|rep),
#                     data = df,
#                     covmat = list(id = K, rep = NULL),
#                     ngamma = 0, nlambda = 1000)
# ghap.manhattan(data = gwas2, chr = "CHR", bp = "BP", y = "LOGP")
# 
# # Correlation between methods
# cor(gwas1$LOGP, gwas2$LOGP)
# plot(gwas1$LOGP, gwas2$LOGP); abline(0,1)

}
