\name{ghap.assoc}
\alias{ghap.assoc}

\title{
Association analysis for HapAlleles
}
\description{
Given a GHap.lmm object (as supplied by the  \code{\link{ghap.lmm}} function), this function returns association statistics for HapAlleles.
}
\usage{
ghap.assoc(response, haplo, weights=NULL, gc = TRUE, only.active.alleles = TRUE,
           batchsize = NULL, ncores = 1, verbose = TRUE)
}

\arguments{
  \item{response}{
   A vector of phenotypes. The vector must be named and all names must be present in the GHap.haplo object.
  }
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{weights}{
   A numeric vector with weights for phenotypes. These weights are treated as diagonal elements of the inverse weight matrix. If not supplied, the analysis is carried out assuming all observations are equally important.
  }
  \item{gc}{
   A logical value specifying whether genomic control should be performed (default = TRUE). Currently, this option does not take effect on HapBlocks.
  }
  \item{only.active.alleles}{
    A logical value specifying whether calculations should be reported only for active haplotype alleles (default = TRUE).
  }
  \item{batchsize}{
   A numeric value controlling the number of HapAlleles to be processed at a time (default = nalleles/10).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
The function returns a data.frame with the following columns:
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
 \item{ALLELE}{
  Haplotype allele identity.
 }
 \item{BETA}{
  BLUE for the fixed effect of the haplotype allele.
 }
 \item{SE}{
  Standard error for the fixed effect.
 }
 \item{FREQ}{
  Frequency of the haplotype allele.
 }
 \item{CHISQ.OBS}{
  Observed value for the test statistics. If gc = TRUE (default), these values are scaled by the inflation factor. Inflation is computed through regression of observed quantiles onto expected quantiles. In order to avoid overestimation, only HapAlleles with test statistics within three standard deviations from the mean are used to compute the inflation factor.
 }
 \item{CHISQ.EXP}{
  Expected values for the test statistics.
 }
 \item{logP}{
  log10(1/P) or -log10(P) for the fixed effect.
 }

}

\details{
This function uses least squares regression to test each HapAllele at a time for association with phenotypes. The fixed effect, error variance and test statistic of a given HapAllele are estimated as:

\deqn{\hat{a}_i = (\mathbf{x}_i^T\mathbf{x}_i)^{-1}\mathbf{x}_i^T\mathbf{y}}
\deqn{VAR(\hat{a}_i) = (\mathbf{x}_i^T\mathbf{x}_i)^{-1}\hat{\sigma}_e^2}
\deqn{t_i^2 = \frac{\hat{a}_i^2}{VAR(\hat{a}_i)}}

Under the null hypothesis that the regression coefficient is zero \eqn{t_i^2 \sim \chi^2(\nu = 1)}. This function supports repeated measures, and records are dully mapped against the vectors of HapGenotypes.

If the vector of responses comprises adjusted records (i.e., residuals) from a linear mixed model, the regression analysis approximates the model implemented in other GWAS tools. However, the user must be aware of two known caveats associated with this approach. First, by pre-adjusting records instead of estimating HapAllele effects based on generalized least squares equations we ignore covariance structure and therefore bias the estimates downwards (Svishcheva et al., 2012). Second, each HapAllele being tested was also potentially included in the kinship matrix in the mixed model analysis, such that the HapAllele is included twice in the model: as fixed and random effect. This problem is known as proximal contamination (Listgarten et al., 2012). In the first case, we can use genomic control to recover p-values to an unbiased scale (Devlin and Roeder, 1999; Amin et al., 2007). However, not much can be done regarding the estimates of the effects. As a general recommendation, if the user is only interested in p-values, the ghap.assoc analysis should be sufficient. When effect estimates are of interest, the user can include the candidate HapAllele as a fixed effect in the full model in ghap.lmm. For the second case, a leave-one-chromosome-out (LOCO analysis) procedure can mitigate proximal contamination (Yang et al., 2014). An alternative to these methods is to use polygenic effects as response instead of residuals. However, this can lead to a higher false-positive rate (Ekine et al., 2014). Finally, if the analysis is intended to screen the genome for recovering effects that are not captured by single markers, a more traditional kinship based on single markers instead of HapAlleles could be used to fit the linear mixed model.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
N. Amin et al. A Genomic Background Based Method for Association Analysis in Related Individuals. PLoS ONE. 2007. 2:e1274.

Y. Da. Multi-allelic haplotype model based on genetic partition for genomic prediction and variance component estimation using SNP markers. BMC Genet. 2015. 16:144.

B. Devlin and K. Roeder. Genomic control for association studies. Biometrics. 1999. 55:997-1004.

C. C. Ekine et al. Why breeding values estimated using familial data should not be used for genome-wide association studies. G3. 2014. 4:341-347.

J. Listgarten et al. Improved linear mixed models for genome-wide association studies. Nat. Methods. 2012. 9:525-526.

G. R. Svishcheva et al. Rapid variance components-based method for whole-genome association analysis. Nat Genet. 2012. 44:1166-1170.

J. Yang et al. Advantages and pitfalls in the application of mixed-model association methods. Nat. Genet. 2014. 46: 100-106.
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy the example data in the current working directory
# exfiles <- ghap.makefile()
# file.copy(from = exfiles, to = "./")
# 
# # Compress phase data
# ghap.compress(input.file = "human", out.file = "human")
# 
# # Load phase data
# phase <- ghap.loadphase(input.file = "human")
# 
# # Generate blocks of 5 markers sliding 5 markers at a time
# blocks <- ghap.blockgen(phase, windowsize = 5, slide = 5, unit = "marker")
# 
# # Haplotyping
# ghap.haplotyping(phase = phase, blocks = blocks, outfile = "human",
#                  binary = T, ncores = 1)
#
# # Load haplotype genotypes
# haplo <- ghap.loadhaplo(input.file = "human")
#
# # Subset common haplotypes in Europeans
# EUR.ids <- haplo$id[which(haplo$pop \%in\% c("TSI","CEU"))]
# haplo <- ghap.subsethaplo(haplo,EUR.ids,haplo$allele.in)
# hapstats <- ghap.hapstats(haplo, ncores = 1)
# common <- hapstats$TYPE \%in\% c("REGULAR","MAJOR") &
#  hapstats$FREQ > 0.05 &
#  hapstats$FREQ < 0.95
# haplo <- ghap.subsethaplo(haplo,EUR.ids,common)
# 
# #Compute relationship matrix
# K <- ghap.kinship(haplo, ncores = 1)
# 
# # Quantitative trait with 50\% heritability
# # Unbalanced repeated measurements (0 to 30)
# # Two major haplotypes accounting for 50\% of the genetic variance
# myseed <- 123456789
# set.seed(myseed)
# major <- sample(which(haplo$allele.in == TRUE), size = 2)
# g2 <- runif(n = 2, min = 0, max = 1)
# g2 <- (g2/sum(g2))*0.5
# sim <- ghap.simpheno(haplo, kinship = K, h2 = 0.5, g2 = g2, nrep = 30,
#                      balanced = FALSE, major = major, seed = myseed)
#
# #Fit model using REML
# model <- ghap.lmm(fixed = phenotype ~ 1, random = ~ individual,
#                   covmat = list(individual = K), data = sim$data)
# 
#
# ### RUN ###
# 
# #HapAllele GWAS using GEBVs as response
# pheno <- model$random$individual
# gwas1 <- ghap.assoc(response = pheno, haplo = haplo, ncores = 1)
# 
# #HapAllele GWAS using GEBVs as response
# #Weight observations by number of repeated measurements
# pheno <- model$random$individual
# w <- table(sim$data$individual)
# w <- w + mean(w)
# w <- w[names(pheno)]
# gwas2 <- ghap.assoc(response = pheno, haplo = haplo, ncores = 1, weights = w)
# 
# #HapAllele GWAS using residuals as response
# pheno <- model$residuals
# names(pheno) <- sim$data$individual
# gwas3 <- ghap.assoc(response = pheno, haplo = haplo, ncores = 1)
# 
# #Plot results
# plot(gwas1$BP1/1e+6,gwas1$logP,pch=20,col="darkgreen",ylim=c(0,20),
#      xlab="Position (in Mb)",ylab=expression(-log[10](p)))
# points(gwas2$BP1/1e+6,gwas2$logP,pch=20,col="gray")
# points(gwas3$BP1/1e+6,gwas3$logP,pch=20,col="blue")
# abline(v=haplo$bp1[major]/1e+6,lty=3)
# abline(h=-log10(0.05/nrow(gwas1)),lty=3)
# legend("topleft",legend = c("GEBVs","weighted GEBVs","residuals"),
#        pch = 20,col=c("darkgreen","gray","blue"))

}
