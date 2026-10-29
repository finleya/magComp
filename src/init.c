#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP magcomp_cov_models(void);
extern SEXP magcomp_cor_matrix(SEXP, SEXP, SEXP);
extern SEXP magcomp_sampler(SEXP config);
extern SEXP magcomp_predict_magnitude_cpp(SEXP XS, SEXP magnitudeCoefficientsS,
                                          SEXP effectS, SEXP covariateS);
extern SEXP magcomp_predict_iid_scalar_cpp(SEXP effectS,
                                           SEXP observedIndexS,
                                           SEXP newIndexS,
                                           SEXP nNewLevelS,
                                           SEXP priorVarS);
extern SEXP magcomp_predict_iid_factor_cpp(SEXP effectS,
                                           SEXP observedIndexS,
                                           SEXP newIndexS,
                                           SEXP nNewLevelS);
extern SEXP magcomp_spatial_predict_cpp(SEXP effectS, SEXP coordsObsS,
                                        SEXP coordsNewS, SEXP phiS,
                                        SEXP sigmaS, SEXP mS,
                                        SEXP methodS, SEXP isFactorS, SEXP modelS, SEXP stScaleS);
extern SEXP magcomp_predict_block_proportion_cpp(SEXP XS, SEXP coefS,
                                                 SEXP loadingS,
                                                 SEXP factorS,
                                                 SEXP covariateS,
                                                 SEXP parentIndexS,
                                                 SEXP JS,
                                                 SEXP linkS);
extern SEXP magcomp_combine_conditional_cpp(SEXP parentS, SEXP childS);
extern SEXP magcomp_predict_category_cpp(SEXP magnitudeS,
                                         SEXP proportionS,
                                         SEXP areaS,
                                         SEXP typeS);

static const R_CallMethodDef CallEntries[] = {
  {"magcomp_cov_models", (DL_FUNC) &magcomp_cov_models, 0},
  {"magcomp_cor_matrix", (DL_FUNC) &magcomp_cor_matrix, 3},
  {"magcomp_sampler", (DL_FUNC) &magcomp_sampler, 1},
  {"magcomp_predict_magnitude_cpp",
   (DL_FUNC) &magcomp_predict_magnitude_cpp, 4},
  {"magcomp_predict_iid_scalar_cpp",
   (DL_FUNC) &magcomp_predict_iid_scalar_cpp, 5},
  {"magcomp_predict_iid_factor_cpp",
   (DL_FUNC) &magcomp_predict_iid_factor_cpp, 4},
  {"magcomp_spatial_predict_cpp",
   (DL_FUNC) &magcomp_spatial_predict_cpp, 10},
  {"magcomp_predict_block_proportion_cpp",
   (DL_FUNC) &magcomp_predict_block_proportion_cpp, 8},
  {"magcomp_combine_conditional_cpp",
   (DL_FUNC) &magcomp_combine_conditional_cpp, 2},
  {"magcomp_predict_category_cpp",
   (DL_FUNC) &magcomp_predict_category_cpp, 4},
  {NULL, NULL, 0}
};

void R_init_magComp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
