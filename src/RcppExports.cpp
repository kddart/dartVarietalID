// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dartDistanceMatrixCpp
NumericMatrix dartDistanceMatrixCpp(const IntegerMatrix& snpData, const LogicalVector& sampleWiseAnalysis, const IntegerVector& nThreads);
RcppExport SEXP _dartVarietalID_dartDistanceMatrixCpp(SEXP snpDataSEXP, SEXP sampleWiseAnalysisSEXP, SEXP nThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type snpData(snpDataSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type sampleWiseAnalysis(sampleWiseAnalysisSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type nThreads(nThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(dartDistanceMatrixCpp(snpData, sampleWiseAnalysis, nThreads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dartVarietalID_dartDistanceMatrixCpp", (DL_FUNC) &_dartVarietalID_dartDistanceMatrixCpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_dartVarietalID(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
