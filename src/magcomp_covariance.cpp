#define R_NO_REMAP
#include "magcomp_covariance.h"

extern "C" SEXP magcomp_cov_models()
{
  return magcomp::nativeCall([]() {
    constexpr int count = sizeof(magcomp::covarianceModels) / sizeof(magcomp::CovModel);
    SEXP out = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, count));
    SEXP names = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, count));
    int i = 0;
    for(const auto &model : magcomp::covarianceModels){
      SEXP entry = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, 3));
      SEXP fields = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, 3));
      const char *keys[] = {"parameters", "domains", "space_time"};
      for(int j = 0; j < 3; ++j) SET_STRING_ELT(fields, j, magcomp::guarded_Rf_mkChar(keys[j]));
      SEXP pars = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, model.size));
      SEXP domains = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, model.size));
      for(int p = 0; p < model.size; ++p){
        SET_STRING_ELT(pars, p, magcomp::guarded_Rf_mkChar(model.parameters[p]));
        const char *domain = model.domains[p] == magcomp::CovDomain::positive ? "positive" :
          model.domains[p] == magcomp::CovDomain::positiveUnit ? "positive_unit" :
          model.domains[p] == magcomp::CovDomain::unit ? "unit" : "nonnegative";
        SET_STRING_ELT(domains, p, magcomp::guarded_Rf_mkChar(domain));
      }
      SET_VECTOR_ELT(entry, 0, pars);
      SET_VECTOR_ELT(entry, 1, domains);
      SET_VECTOR_ELT(entry, 2, magcomp::rCall([&]() { return Rf_ScalarLogical(model.spaceTime); }));
      magcomp::guarded_Rf_setAttrib(entry, R_NamesSymbol, fields);
      SET_VECTOR_ELT(out, i, entry);
      SET_STRING_ELT(names, i++, magcomp::guarded_Rf_mkChar(model.name));
      UNPROTECT(4);
    }
    magcomp::guarded_Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
  });
}

// Shared evaluator used by numerical tests and development simulations.
extern "C" SEXP magcomp_cor_matrix(SEXP coordsS, SEXP modelS, SEXP thetaS)
{
  return magcomp::nativeCall([&]() {
    if(TYPEOF(modelS) != STRSXP || Rf_length(modelS) != 1 || STRING_ELT(modelS, 0) == NA_STRING)
      magcomp::error("covariance model must be a single name");
    const auto &model = magcomp::covarianceModel(CHAR(STRING_ELT(modelS, 0)));
    if(TYPEOF(coordsS) != REALSXP || !Rf_isMatrix(coordsS) || TYPEOF(thetaS) != REALSXP ||
       Rf_length(thetaS) != model.size) magcomp::error("invalid covariance inputs");
    int n = Rf_nrows(coordsS), dim = Rf_ncols(coordsS);
    if(dim < (model.spaceTime ? 2 : 1)) magcomp::error("insufficient coordinate columns");
    for(R_xlen_t i = 0; i < XLENGTH(coordsS); ++i)
      if(!R_FINITE(REAL(coordsS)[i])) magcomp::error("coordinates must be finite");
    magcomp::validateCovariance(model, REAL(thetaS));
    SEXP out = PROTECT(magcomp::guarded_Rf_allocMatrix(REALSXP, n, n));
    for(int j = 0; j < n; ++j){
      magcomp::rCall([]() { R_CheckUserInterrupt(); });
      for(int i = j; i < n; ++i){
        double value = magcomp::correlationAt(model, REAL(thetaS), REAL(coordsS), n, dim, i, j);
        REAL(out)[i + j * n] = REAL(out)[j + i * n] = value;
      }
    }
    UNPROTECT(1);
    return out;
  });
}
