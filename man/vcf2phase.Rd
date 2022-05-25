\name{ghap.vcf2phase}
\alias{ghap.vcf2phase}

\title{
  Convert VCF data into GHap phase
}
\description{
  This function takes phased genotype data in the Variant Call Format (VCF) and converts them into the GHap phase format.
}
\usage{
  ghap.vcf2phase(input.files = NULL, vcf.files = NULL,
                 sample.files = NULL, out.file,
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
  
  \item{vcf.files}{
    Character vector containing the list of VCF files.
  }
  \item{sample.files}{
    Character vector containing the list of SAMPLE files.
  }

  To turn conversion progress-tracking on or off or set the number of cores please use:

  \item{ncores}{
  	A numeric value specfying the number of cores to use (default = 1).
  }
  
  \item{verbose}{
  	A logical value specfying whether log messages should be printed (default = TRUE).
  }
}
\details{
The Variant Call Format (VCF) - as described in https://github.com/samtools/hts-specs - is here manipulated to obtain the GHap phase format. Important: the function does not apply filters to the data, except for skipping multi-allelic variants. Should variants be filtered, the user is advised to pre-process the VCF files with third-party software (such as BCFTools). The FORMAT field should also follow the "GT:..." specification, with genotypes placed first in each sample column. Finally, all genotypes should be phased and take one of the following values: "0|0", "0|1", "1|0" or "1|1".
Warning: this function is not optimized for very large datasets.
}

\references{

H. Li et al. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics. 2009. 25:2078-2079.

H. Li. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011. 27(21):2987-2993.

}
\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\seealso{
  \code{\link{ghap.compress}}, \code{\link{ghap.loadphase}}, \code{\link{ghap.fast2phase}}, \code{\link{ghap.oxford2phase}}
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy the example data in the current working directory
# exfiles <- ghap.makefile(dataset = "example",
#                          format = "vcf",
#                          verbose = TRUE)
# file.copy(from = exfiles, to = "./")
# 
# ### RUN ###
# 
# # Convert from a single genome-wide file
# ghap.vcf2phase(input.files = "example",
#                out.file = "example")
# 
# # Convert from a list of chromosome files
# ghap.vcf2phase(input.files = paste0("example_chr",1:10),
#                out.file = "example")
# 
# # Convert using separate lists for file extensions
# ghap.vcf2phase(vcf.files = paste0("example_chr",1:10,".vcf"),
#                sample.files = paste0("example_chr",1:10,".sample"),
#                out.file = "example")

}
