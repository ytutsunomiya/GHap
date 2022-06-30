# GHap

Haplotype calling from phased marker data. Given user-defined haplotype blocks (HapBlock), the package identifies the different haplotype alleles (HapAllele) present in the data and scores sample haplotype allele genotypes (HapGenotype) based on HapAllele dose (i.e. 0, 1 or 2 copies). The output is not only useful for analyses that can handle multi-allelic markers, but is also conveniently formatted for existing pipelines intended for bi-allelic markers.

The package was first described in Bioinformatics by Utsunomiya et al. (2016, <doi:10.1093/bioinformatics/btw356>). Since the v2 release, the package provides functions for unsupervised and supervised detection of ancestry tracks. The methods implemented in these functions were described in an article published in Methods in Ecology and Evolution by Utsunomiya et al. (2020, <doi:10.1111/2041-210X.13467>).

The source code for v3 has been modified for improved performance and inclusion of new functionality, including analysis of unphased data, runs of homozygosity, sampling methods for virtual gamete mating, mixed model fitting and GWAS.
