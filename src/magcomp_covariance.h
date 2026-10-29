#ifndef MAGCOMP_COVARIANCE_H
#define MAGCOMP_COVARIANCE_H

#include <cmath>
#include <cstring>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "magcomp_runtime.h"

namespace magcomp {

enum class CovDomain { positive, unit, positiveUnit, nonnegative };
struct CovModel {
  const char *name;
  bool spaceTime;
  int size;
  const char *parameters[6];
  CovDomain domains[6];
};

inline const CovModel covarianceModels[] = {
  {"exp", false, 1, {"phi"}, {CovDomain::positive}},
  {"matern", false, 2, {"phi", "nu"},
    {CovDomain::positive, CovDomain::positive}},
  {"sep_exp", true, 2, {"phi", "lambda"},
    {CovDomain::positive, CovDomain::positive}},
  {"multi_res_sep_exp", true, 5,
    {"alpha", "phi_1", "lambda_1", "phi_2", "lambda_2"},
    {CovDomain::unit, CovDomain::positive, CovDomain::positive,
     CovDomain::positive, CovDomain::positive}},
  {"gneiting", true, 6, {"a", "c", "alpha", "beta", "gamma", "delta"},
    {CovDomain::positive, CovDomain::positive, CovDomain::positiveUnit,
     CovDomain::unit, CovDomain::positiveUnit, CovDomain::nonnegative}}
};

inline const CovModel &covarianceModel(const char *name)
{
  for(const auto &model : covarianceModels)
    if(std::strcmp(name, model.name) == 0) return model;
  error("unknown covariance model '%s'", name);
}

inline void validateCovariance(const CovModel &model, const double *theta)
{
  for(int p = 0; p < model.size; ++p){
    double x = theta[p];
    auto domain = model.domains[p];
    bool valid = R_FINITE(x);
    if(domain == CovDomain::positive) valid = valid && x > 0;
    if(domain == CovDomain::positiveUnit) valid = valid && x > 0 && x <= 1;
    if(domain == CovDomain::unit) valid = valid && x >= 0 && x <= 1;
    if(domain == CovDomain::nonnegative) valid = valid && x >= 0;
    if(!valid) error("invalid %s covariance parameter '%s'", model.name, model.parameters[p]);
  }
}

struct CovLag { double space, time; };

inline CovLag covarianceLag(const double *a, int na, const double *b, int nb,
                           int dim, int i, int j, bool spaceTime)
{
  double h = 0;
  int spaceDim = dim - (spaceTime ? 1 : 0);
  for(int d = 0; d < spaceDim; ++d){
    double delta = a[i + d * na] - b[j + d * nb];
    h = std::hypot(h, delta);
  }
  return {h, spaceTime ? std::fabs(a[i + spaceDim * na] - b[j + spaceDim * nb]) : 0.0};
}

inline double correlation(const CovModel &model, const double *theta,
                          CovLag lag, int spaceDim)
{
  const double h = lag.space, u = lag.time;
  if(h == 0 && (!model.spaceTime || u == 0)) return 1.0;
  if(std::strcmp(model.name, "exp") == 0) return std::exp(-theta[0] * h);
  if(std::strcmp(model.name, "matern") == 0){
    double x = theta[0] * h, nu = theta[1];
    if(x == 0) return 1.0;
    if(nu == 0.5) return std::exp(-x);
    // Exponentially scaled Bessel K avoids underflow at large distances.
    // Guard Rmath warnings/errors so native allocations are released on exit.
    double bk = rCall([&]() { return Rf_bessel_k(x, nu, 2.0); });
    double logCor = nu * std::log(x) - (nu - 1) * std::log(2.0) -
      std::lgamma(nu) + std::log(bk) - x;
    if(!R_FINITE(logCor)) error("non-finite Matern correlation (phi*h=%g, nu=%g)", x, nu);
    return std::exp(logCor);
  }
  if(std::strcmp(model.name, "sep_exp") == 0)
    return std::exp(-theta[0] * h) * std::exp(-theta[1] * u);
  if(std::strcmp(model.name, "multi_res_sep_exp") == 0)
    return theta[0] * std::exp(-theta[1] * h) * std::exp(-theta[2] * u) +
      (1 - theta[0]) * std::exp(-theta[3] * h) * std::exp(-theta[4] * u);
  if(std::strcmp(model.name, "gneiting") != 0)
    error("no correlation evaluator registered for %s", model.name);
  double logTime = 0;
  if(u > 0){
    double x = std::log(theta[0]) + 2 * theta[2] * std::log(u);
    logTime = x > 0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
  }
  double spatial = h == 0 ? 0 : std::exp(std::log(theta[1]) +
    2 * theta[4] * std::log(h) - theta[3] * theta[4] * logTime);
  return std::exp(-(theta[5] + 0.5 * spaceDim) * logTime - spatial);
}

inline double correlationAt(const CovModel &model, const double *theta,
                            const double *coords, int n, int dim, int i, int j)
{
  return correlation(model, theta, covarianceLag(coords, n, coords, n, dim, i, j, model.spaceTime),
                     dim - (model.spaceTime ? 1 : 0));
}

} // namespace magcomp
#endif
