// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// citbincvr
void citbincvr(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, Rcpp::NumericVector CG, int& maxit, int& nrow, int& ncol, int& ncolc, int& ncolct, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, int& rseed);
RcppExport SEXP _CITpkg_citbincvr(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP CGSEXP, SEXP maxitSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP ncolcSEXP, SEXP ncolctSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval3SEXP, SEXP pval4SEXP, SEXP pval3ncSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolc(ncolcSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolct(ncolctSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3(pval3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3nc(pval3ncSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citbincvr(L, G, T, C, CG, maxit, nrow, ncol, ncolc, ncolct, pval1, pval2, pval3, pval4, pval3nc, rseed);
    return R_NilValue;
END_RCPP
}
// citbincvr_linear
void citbincvr_linear(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, Rcpp::NumericVector CG, int& maxit, int& nrow, int& ncol, int& ncolc, int& ncolct, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, int& rseed);
RcppExport SEXP _CITpkg_citbincvr_linear(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP CGSEXP, SEXP maxitSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP ncolcSEXP, SEXP ncolctSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval3SEXP, SEXP pval4SEXP, SEXP pval3ncSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolc(ncolcSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolct(ncolctSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3(pval3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3nc(pval3ncSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citbincvr_linear(L, G, T, C, CG, maxit, nrow, ncol, ncolc, ncolct, pval1, pval2, pval3, pval4, pval3nc, rseed);
    return R_NilValue;
END_RCPP
}
// citbinmcvr
void citbinmcvr(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, int& maxit, int& nrow, int& dfz, int& dfx, int& dfc, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, int& rseed);
RcppExport SEXP _CITpkg_citbinmcvr(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP maxitSEXP, SEXP nrowSEXP, SEXP dfzSEXP, SEXP dfxSEXP, SEXP dfcSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval4SEXP, SEXP pval3ncSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type dfz(dfzSEXP);
    Rcpp::traits::input_parameter< int& >::type dfx(dfxSEXP);
    Rcpp::traits::input_parameter< int& >::type dfc(dfcSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3nc(pval3ncSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citbinmcvr(L, G, T, C, maxit, nrow, dfz, dfx, dfc, pval1, pval2, pval4, pval3nc, rseed);
    return R_NilValue;
END_RCPP
}
// citbinmcvr_linear
void citbinmcvr_linear(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, int& maxit, int& nrow, int& dfz, int& dfx, int& dfc, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, int& rseed);
RcppExport SEXP _CITpkg_citbinmcvr_linear(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP maxitSEXP, SEXP nrowSEXP, SEXP dfzSEXP, SEXP dfxSEXP, SEXP dfcSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval4SEXP, SEXP pval3ncSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type dfz(dfzSEXP);
    Rcpp::traits::input_parameter< int& >::type dfx(dfxSEXP);
    Rcpp::traits::input_parameter< int& >::type dfc(dfcSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3nc(pval3ncSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citbinmcvr_linear(L, G, T, C, maxit, nrow, dfz, dfx, dfc, pval1, pval2, pval4, pval3nc, rseed);
    return R_NilValue;
END_RCPP
}
// citbinmpcvr
void citbinmpcvr(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, int& maxit, int& permit, int& permnum, int& nobs, int& dfz, int& dfx, int& dfc, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, Rcpp::IntegerVector perm_index, int& rseed);
RcppExport SEXP _CITpkg_citbinmpcvr(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP maxitSEXP, SEXP permitSEXP, SEXP permnumSEXP, SEXP nobsSEXP, SEXP dfzSEXP, SEXP dfxSEXP, SEXP dfcSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval4SEXP, SEXP pval3ncSEXP, SEXP perm_indexSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type permit(permitSEXP);
    Rcpp::traits::input_parameter< int& >::type permnum(permnumSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int& >::type dfz(dfzSEXP);
    Rcpp::traits::input_parameter< int& >::type dfx(dfxSEXP);
    Rcpp::traits::input_parameter< int& >::type dfc(dfcSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3nc(pval3ncSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type perm_index(perm_indexSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citbinmpcvr(L, G, T, C, maxit, permit, permnum, nobs, dfz, dfx, dfc, pval1, pval2, pval4, pval3nc, perm_index, rseed);
    return R_NilValue;
END_RCPP
}
// citbinmpcvr_linear
void citbinmpcvr_linear(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, int& maxit, int& permit, int& permnum, int& nobs, int& dfz, int& dfx, int& dfc, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, Rcpp::IntegerVector perm_index, int& rseed);
RcppExport SEXP _CITpkg_citbinmpcvr_linear(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP maxitSEXP, SEXP permitSEXP, SEXP permnumSEXP, SEXP nobsSEXP, SEXP dfzSEXP, SEXP dfxSEXP, SEXP dfcSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval4SEXP, SEXP pval3ncSEXP, SEXP perm_indexSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type permit(permitSEXP);
    Rcpp::traits::input_parameter< int& >::type permnum(permnumSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int& >::type dfz(dfzSEXP);
    Rcpp::traits::input_parameter< int& >::type dfx(dfxSEXP);
    Rcpp::traits::input_parameter< int& >::type dfc(dfcSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3nc(pval3ncSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type perm_index(perm_indexSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citbinmpcvr_linear(L, G, T, C, maxit, permit, permnum, nobs, dfz, dfx, dfc, pval1, pval2, pval4, pval3nc, perm_index, rseed);
    return R_NilValue;
END_RCPP
}
// citbinpcvr
void citbinpcvr(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, Rcpp::NumericVector CG, int& maxit, int& permit, int& boots, int& nrow, int& ncol, int& ncolc, int& ncolct, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, Rcpp::IntegerVector perm_index, int& rseed);
RcppExport SEXP _CITpkg_citbinpcvr(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP CGSEXP, SEXP maxitSEXP, SEXP permitSEXP, SEXP bootsSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP ncolcSEXP, SEXP ncolctSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval3SEXP, SEXP pval4SEXP, SEXP pval3ncSEXP, SEXP perm_indexSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type permit(permitSEXP);
    Rcpp::traits::input_parameter< int& >::type boots(bootsSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolc(ncolcSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolct(ncolctSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3(pval3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3nc(pval3ncSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type perm_index(perm_indexSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citbinpcvr(L, G, T, C, CG, maxit, permit, boots, nrow, ncol, ncolc, ncolct, pval1, pval2, pval3, pval4, pval3nc, perm_index, rseed);
    return R_NilValue;
END_RCPP
}
// citbinpcvr_linear
void citbinpcvr_linear(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, Rcpp::NumericVector CG, int& maxit, int& permit, int& boots, int& nrow, int& ncol, int& ncolc, int& ncolct, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, Rcpp::IntegerVector perm_index, int& rseed);
RcppExport SEXP _CITpkg_citbinpcvr_linear(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP CGSEXP, SEXP maxitSEXP, SEXP permitSEXP, SEXP bootsSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP ncolcSEXP, SEXP ncolctSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval3SEXP, SEXP pval4SEXP, SEXP pval3ncSEXP, SEXP perm_indexSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type permit(permitSEXP);
    Rcpp::traits::input_parameter< int& >::type boots(bootsSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolc(ncolcSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolct(ncolctSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3(pval3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3nc(pval3ncSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type perm_index(perm_indexSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citbinpcvr_linear(L, G, T, C, CG, maxit, permit, boots, nrow, ncol, ncolc, ncolct, pval1, pval2, pval3, pval4, pval3nc, perm_index, rseed);
    return R_NilValue;
END_RCPP
}
// citconlog2cvr
void citconlog2cvr(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, Rcpp::NumericVector CG, int& nrow, int& ncol, int& ncolc, int& ncolct, Rcpp::NumericVector pval, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, int& maxit, int& rseed);
RcppExport SEXP _CITpkg_citconlog2cvr(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP CGSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP ncolcSEXP, SEXP ncolctSEXP, SEXP pvalSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval3SEXP, SEXP pval4SEXP, SEXP maxitSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolc(ncolcSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolct(ncolctSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3(pval3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citconlog2cvr(L, G, T, C, CG, nrow, ncol, ncolc, ncolct, pval, pval1, pval2, pval3, pval4, maxit, rseed);
    return R_NilValue;
END_RCPP
}
// citconlog2cvr_linear
void citconlog2cvr_linear(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, Rcpp::NumericVector CG, int& nrow, int& ncol, int& ncolc, int& ncolct, Rcpp::NumericVector pval, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, int& maxit, int& rseed);
RcppExport SEXP _CITpkg_citconlog2cvr_linear(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP CGSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP ncolcSEXP, SEXP ncolctSEXP, SEXP pvalSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval3SEXP, SEXP pval4SEXP, SEXP maxitSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolc(ncolcSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolct(ncolctSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3(pval3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citconlog2cvr_linear(L, G, T, C, CG, nrow, ncol, ncolc, ncolct, pval, pval1, pval2, pval3, pval4, maxit, rseed);
    return R_NilValue;
END_RCPP
}
// citconlog3pcvr
void citconlog3pcvr(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, Rcpp::NumericVector CG, int& nrow, int& ncol, int& ncolc, int& ncolct, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, int& maxit, int& permit, int& boots, Rcpp::NumericVector Pind, int& rseed);
RcppExport SEXP _CITpkg_citconlog3pcvr(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP CGSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP ncolcSEXP, SEXP ncolctSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval3SEXP, SEXP pval4SEXP, SEXP maxitSEXP, SEXP permitSEXP, SEXP bootsSEXP, SEXP PindSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolc(ncolcSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolct(ncolctSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3(pval3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type permit(permitSEXP);
    Rcpp::traits::input_parameter< int& >::type boots(bootsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Pind(PindSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citconlog3pcvr(L, G, T, C, CG, nrow, ncol, ncolc, ncolct, pval1, pval2, pval3, pval4, maxit, permit, boots, Pind, rseed);
    return R_NilValue;
END_RCPP
}
// citconlog3pcvr_linear
void citconlog3pcvr_linear(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, Rcpp::NumericVector CG, int& nrow, int& ncol, int& ncolc, int& ncolct, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, int& maxit, int& permit, int& boots, Rcpp::NumericVector Pind, int& rseed);
RcppExport SEXP _CITpkg_citconlog3pcvr_linear(SEXP LSEXP, SEXP GSEXP, SEXP TSEXP, SEXP CSEXP, SEXP CGSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP ncolcSEXP, SEXP ncolctSEXP, SEXP pval1SEXP, SEXP pval2SEXP, SEXP pval3SEXP, SEXP pval4SEXP, SEXP maxitSEXP, SEXP permitSEXP, SEXP bootsSEXP, SEXP PindSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< int& >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int& >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolc(ncolcSEXP);
    Rcpp::traits::input_parameter< int& >::type ncolct(ncolctSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval1(pval1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval2(pval2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval3(pval3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pval4(pval4SEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int& >::type permit(permitSEXP);
    Rcpp::traits::input_parameter< int& >::type boots(bootsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Pind(PindSEXP);
    Rcpp::traits::input_parameter< int& >::type rseed(rseedSEXP);
    citconlog3pcvr_linear(L, G, T, C, CG, nrow, ncol, ncolc, ncolct, pval1, pval2, pval3, pval4, maxit, permit, boots, Pind, rseed);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CITpkg_citbincvr", (DL_FUNC) &_CITpkg_citbincvr, 16},
    {"_CITpkg_citbincvr_linear", (DL_FUNC) &_CITpkg_citbincvr_linear, 16},
    {"_CITpkg_citbinmcvr", (DL_FUNC) &_CITpkg_citbinmcvr, 14},
    {"_CITpkg_citbinmcvr_linear", (DL_FUNC) &_CITpkg_citbinmcvr_linear, 14},
    {"_CITpkg_citbinmpcvr", (DL_FUNC) &_CITpkg_citbinmpcvr, 17},
    {"_CITpkg_citbinmpcvr_linear", (DL_FUNC) &_CITpkg_citbinmpcvr_linear, 17},
    {"_CITpkg_citbinpcvr", (DL_FUNC) &_CITpkg_citbinpcvr, 19},
    {"_CITpkg_citbinpcvr_linear", (DL_FUNC) &_CITpkg_citbinpcvr_linear, 19},
    {"_CITpkg_citconlog2cvr", (DL_FUNC) &_CITpkg_citconlog2cvr, 16},
    {"_CITpkg_citconlog2cvr_linear", (DL_FUNC) &_CITpkg_citconlog2cvr_linear, 16},
    {"_CITpkg_citconlog3pcvr", (DL_FUNC) &_CITpkg_citconlog3pcvr, 18},
    {"_CITpkg_citconlog3pcvr_linear", (DL_FUNC) &_CITpkg_citconlog3pcvr_linear, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_CITpkg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
