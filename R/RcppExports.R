# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.compressCpp <- function(infile, outfile, nunits, nbits, tbits, fmode) {
    .Call(`_GHap_compressCpp`, infile, outfile, nunits, nbits, tbits, fmode)
}

.sliceCpp <- function(binfile, mode, nvars, nids, phased, imp, iidx, hidx, vidx) {
    .Call(`_GHap_sliceCpp`, binfile, mode, nvars, nids, phased, imp, iidx, hidx, vidx)
}
