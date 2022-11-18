// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compressCpp
int compressCpp(const char* infile, const char* outfile, const int nunits, const int nbits, const int tbits, const int fmode);
RcppExport SEXP _GHap_compressCpp(SEXP infileSEXP, SEXP outfileSEXP, SEXP nunitsSEXP, SEXP nbitsSEXP, SEXP tbitsSEXP, SEXP fmodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type infile(infileSEXP);
    Rcpp::traits::input_parameter< const char* >::type outfile(outfileSEXP);
    Rcpp::traits::input_parameter< const int >::type nunits(nunitsSEXP);
    Rcpp::traits::input_parameter< const int >::type nbits(nbitsSEXP);
    Rcpp::traits::input_parameter< const int >::type tbits(tbitsSEXP);
    Rcpp::traits::input_parameter< const int >::type fmode(fmodeSEXP);
    rcpp_result_gen = Rcpp::wrap(compressCpp(infile, outfile, nunits, nbits, tbits, fmode));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GHap_compressCpp", (DL_FUNC) &_GHap_compressCpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_GHap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}