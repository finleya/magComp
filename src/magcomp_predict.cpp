#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

#define USE_FC_LEN_T
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "magcomp_runtime.h"
#include "magcomp_covariance.h"

#ifndef FCONE
# define FCONE
#endif

static std::vector<int> getDims(SEXP x, const char *name)
{
  SEXP dimS = magcomp::guarded_Rf_getAttrib(x, R_DimSymbol);
  if(Rf_isNull(dimS))
    magcomp::error("%s must have a dim attribute", name);

  std::vector<int> out(magcomp::guarded_Rf_length(dimS));
  for(int i = 0; i < magcomp::guarded_Rf_length(dimS); i++)
    out[i] = magcomp::guarded_INTEGER(dimS)[i];
  return out;
}

static const double *realData(SEXP x, const char *name)
{
  if(TYPEOF(x) != REALSXP)
    magcomp::error("%s must be numeric", name);
  return magcomp::guarded_REAL(x);
}

static const int *intData(SEXP x, const char *name)
{
  if(TYPEOF(x) != INTSXP)
    magcomp::error("%s must be integer", name);
  return magcomp::guarded_INTEGER(x);
}

static SEXP allocRealArray(const std::vector<int> &dim)
{
  SEXP dimS = PROTECT(magcomp::guarded_Rf_allocVector(INTSXP, (int) dim.size()));
  for(int i = 0; i < (int) dim.size(); i++)
    magcomp::guarded_INTEGER(dimS)[i] = dim[i];

  SEXP out = magcomp::guarded_Rf_allocArray(REALSXP, dimS);
  UNPROTECT(1);
  return out;
}

static double invLogit(double x)
{
  if(x >= 0.0){
    double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }

  double z = std::exp(x);
  return z / (1.0 + z);
}

enum CompositionLink {
  COMPOSITION_LINK_LOGIT = 0,
  COMPOSITION_LINK_PROBIT = 1
};

static CompositionLink compositionLinkFromSexp(SEXP linkS)
{
  if(magcomp::guarded_Rf_length(linkS) < 1)
    magcomp::error("composition link must have length one");
  const char *link = CHAR(magcomp::guarded_Rf_asChar(linkS));
  if(std::strcmp(link, "logit") == 0)
    return COMPOSITION_LINK_LOGIT;
  if(std::strcmp(link, "probit") == 0)
    return COMPOSITION_LINK_PROBIT;
  magcomp::error("unsupported composition link '%s'", link);
  return COMPOSITION_LINK_LOGIT;
}

static double invCompositionLink(double x, CompositionLink link)
{
  if(link == COMPOSITION_LINK_PROBIT)
    return Rf_pnorm5(x, 0.0, 1.0, 1, 0);
  return invLogit(x);
}

static double coordDist(const double *coordsA, int nA,
                        const double *coordsB, int nB,
                        int dim, int i, int j)
{
  double out = 0.0;
  for(int k = 0; k < dim; k++){
    double d = coordsA[i + k * nA] - coordsB[j + k * nB];
    out += d * d;
  }
  return std::sqrt(out);
}

static void fillCorrelationSym(double *A, const double *coords, int n,
                           int coordDim, const magcomp::CovModel &model, const double *theta)
{
  for(int j = 0; j < n; j++){
    A[j + j * n] = 1.0;
    for(int i = j + 1; i < n; i++){
      double value = magcomp::correlationAt(model, theta, coords, n, coordDim, i, j);
      A[i + j * n] = value;
      A[j + i * n] = value;
    }
  }
}

static void fillCorrelationCross(double *A, const double *coordsA, int nA,
                             const double *coordsB, int nB, int coordDim,
                             const magcomp::CovModel &model, const double *theta)
{
  for(int j = 0; j < nB; j++){
    for(int i = 0; i < nA; i++){
      A[i + j * nA] = magcomp::correlation(model, theta, magcomp::covarianceLag(coordsA, nA, coordsB, nB, coordDim, i, j, model.spaceTime), coordDim - model.spaceTime);
    }
  }
}

static void choleskyLower(const std::vector<double> &A,
                          std::vector<double> &chol, int n,
                          const char *context)
{
  int info = 0;
  const char *lower = "L";
  chol = A;
  magcomp::rCall([&]() {
    F77_CALL(dpotrf)(lower, &n, chol.data(), &n, &info FCONE);
  });
  if(info != 0)
    magcomp::error("dpotrf failed in %s (info = %i)", context, info);
}

static void solveCholLower(std::vector<double> &chol, int n,
                           double *rhs, int nrhs, const char *context)
{
  int info = 0;
  const char *lower = "L";
  magcomp::rCall([&]() {
    F77_CALL(dpotrs)(lower, &n, &nrhs, chol.data(), &n, rhs, &n, &info FCONE);
  });
  if(info != 0)
    magcomp::error("dpotrs failed in %s", context);
}

static double getSampleScalar(SEXP x, int s, int ch, int nSamples,
                              double defaultValue)
{
  if(Rf_isNull(x))
    return defaultValue;

  const double *v = realData(x, "sample scalar");
  if(magcomp::guarded_Rf_length(x) == 1)
    return v[0];

  std::vector<int> d = getDims(x, "sample scalar");
  if(d.size() == 2)
    return v[s];
  if(d.size() == 3)
    return v[s + ch * nSamples];

  magcomp::error("sample scalar must be scalar, matrix, or array");
  return defaultValue;
}

static double getSampleFactorScalar(SEXP x,
                                    int k,
                                    int s,
                                    int ch,
                                    int K,
                                    int nSamples,
                                    double defaultValue)
{
  if(Rf_isNull(x))
    return defaultValue;

  const double *v = realData(x, "factor sample scalar");
  if(magcomp::guarded_Rf_length(x) == 1)
    return v[0];
  if(magcomp::guarded_Rf_length(x) == K && Rf_isNull(magcomp::guarded_Rf_getAttrib(x, R_DimSymbol)))
    return v[k];

  std::vector<int> d = getDims(x, "factor sample scalar");
  if(d.size() == 2){
    if(d[0] == 1)
      return v[s];
    if(d[0] == K)
      return v[k + s * K];
  }
  if(d.size() == 3){
    if(d[0] == 1)
      return v[s + ch * nSamples];
    if(d[0] == K)
      return v[k + s * K + (size_t) ch * K * nSamples];
  }

  magcomp::error("factor sample scalar must be scalar, length-K vector, matrix, "
           "or array");
  return defaultValue;
}

typedef struct {
  const double *loading;
  const double *factor;
  const double *covariate;
  int K;
} FactorPredictionTerm;

static void addFactorPredictionTerm(std::vector<FactorPredictionTerm>& terms,
                                    SEXP loadingS,
                                    SEXP factorS,
                                    SEXP covariateS,
                                    int n,
                                    int coefStick,
                                    int nSamples,
                                    int nChains,
                                    const char *context)
{
  FactorPredictionTerm term;

  if(Rf_isNull(loadingS) || Rf_isNull(factorS))
    magcomp::error("%s loading and factor samples must be supplied together",
             context);

  term.loading = realData(loadingS, "loading_samples");
  term.factor = realData(factorS, "factor_samples");
  term.covariate = NULL;

  std::vector<int> lDim = getDims(loadingS, "loading_samples");
  std::vector<int> fDim = getDims(factorS, "factor_samples");
  if((lDim.size() != 3 && lDim.size() != 4) ||
     (fDim.size() != 3 && fDim.size() != 4) ||
     lDim[0] != coefStick || fDim[0] != n ||
     lDim[2] != nSamples || fDim[2] != nSamples ||
     (int) lDim.size() != (int) fDim.size()){
    magcomp::error("%s loading and factor samples have incompatible dimensions",
             context);
  }
  term.K = lDim[1];
  if(fDim[1] != term.K)
    magcomp::error("%s loading and factor samples have incompatible K", context);
  if(nChains == 1){
    if(lDim.size() != 3 || fDim.size() != 3)
      magcomp::error("%s factor samples must match single-chain coef_samples",
               context);
  }else{
    if(lDim.size() != 4 || fDim.size() != 4 ||
       lDim[3] != nChains || fDim[3] != nChains){
      magcomp::error("%s factor samples must match multi-chain coef_samples",
               context);
    }
  }

  if(!Rf_isNull(covariateS)){
    term.covariate = realData(covariateS, "factor_covariate");
    if(magcomp::guarded_Rf_length(covariateS) != n)
      magcomp::error("%s factor_covariate must have one value per prediction row",
               context);
  }

  terms.push_back(term);
}

static void gpPredictDraw(const double *wObs, const double *coordsObs,
                          const double *coordsNew, int nObs, int nNew,
                          int coordDim, const magcomp::CovModel &model, const double *theta, double sigmaSq,
                          double *out)
{
  std::vector<double> RObs((size_t) nObs * nObs);
  std::vector<double> cholObs;
  std::vector<double> RNewObs((size_t) nNew * nObs);
  std::vector<double> alpha(nObs);
  std::vector<double> B((size_t) nObs * nNew);
  std::vector<double> RNew((size_t) nNew * nNew);
  std::vector<double> cond((size_t) nNew * nNew);
  std::vector<double> cholNew;

  fillCorrelationSym(RObs.data(), coordsObs, nObs, coordDim, model, theta);
  choleskyLower(RObs, cholObs, nObs, "GP prediction observed covariance");

  std::memcpy(alpha.data(), wObs, sizeof(double) * (size_t) nObs);
  solveCholLower(cholObs, nObs, alpha.data(), 1,
                 "GP prediction observed solve");

  fillCorrelationCross(RNewObs.data(), coordsNew, nNew, coordsObs, nObs,
                   coordDim, model, theta);
  std::memcpy(B.data(), RNewObs.data(), sizeof(double) * (size_t) nNew * nObs);

  for(int j = 0; j < nNew; j++){
    for(int i = 0; i < nObs; i++)
      B[i + j * nObs] = RNewObs[j + i * nNew];
  }
  solveCholLower(cholObs, nObs, B.data(), nNew,
                 "GP prediction cross-covariance solve");

  for(int i = 0; i < nNew; i++){
    out[i] = 0.0;
    for(int h = 0; h < nObs; h++)
      out[i] += RNewObs[i + h * nNew] * alpha[h];
  }

  fillCorrelationSym(RNew.data(), coordsNew, nNew, coordDim, model, theta);
  for(int j = 0; j < nNew; j++){
    for(int i = 0; i < nNew; i++){
      double value = RNew[i + j * nNew];
      for(int h = 0; h < nObs; h++)
        value -= RNewObs[i + h * nNew] * B[h + j * nObs];
      cond[i + j * nNew] = value;
    }
  }
  for(int j = 0; j < nNew; j++){
    for(int i = j; i < nNew; i++){
      double value = sigmaSq * 0.5 * (cond[i+j*nNew] + cond[j+i*nNew]);
      cond[i+j*nNew] = cond[j+i*nNew] = value;
    }
  }

  choleskyLower(cond, cholNew, nNew, "GP prediction covariance");

  std::vector<double> z(nNew);
  for(int i = 0; i < nNew; i++)
    z[i] = Rf_rnorm(0.0, 1.0);

  for(int i = 0; i < nNew; i++){
    for(int j = 0; j <= i; j++)
      out[i] += cholNew[i + j * nNew] * z[j];
  }
}

struct NNGPPredictionLocation {
  int newId;
  std::vector<int> neighbor;
  std::vector<magcomp::CovLag> distanceToNew;
  // Strict lower triangle, by column. The diagonal is zero and distances
  // are symmetric, so only q * (q - 1) / 2 entries need to be retained.
  std::vector<magcomp::CovLag> neighborDistance;
};

// Coordinates and m determine the directed prediction graph. Build it once
// per spatial-term call, independently of chains, draws, and factor covariance parameters.
static std::vector<NNGPPredictionLocation> buildNNGPPredictionGeometry(
    const double *coordsObs, const double *coordsNew,
    int nObs, int nNew, int coordDim, int m, bool spaceTime, double stScale)
{
  std::vector<int> order(nNew);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    for(int k = 0; k < coordDim; k++){
      double av = coordsNew[a + k * nNew];
      double bv = coordsNew[b + k * nNew];
      if(av < bv)
        return true;
      if(av > bv)
        return false;
    }
    return a < b;
  });

  std::vector<NNGPPredictionLocation> geometry(nNew);

  for(int pos = 0; pos < nNew; pos++){
    int newId = order[pos];
    int nCandidate = nObs + pos;
    std::vector<std::pair<double, int> > dist;
    dist.reserve(nCandidate);

    for(int i = 0; i < nObs; i++){
      auto lag = magcomp::covarianceLag(coordsObs, nObs, coordsNew, nNew, coordDim, i, newId, spaceTime);
      double d = std::hypot(lag.space, stScale * lag.time);
      dist.push_back(std::make_pair(d, i));
    }
    for(int h = 0; h < pos; h++){
      int prevId = order[h];
      auto lag = magcomp::covarianceLag(coordsNew, nNew, coordsNew, nNew, coordDim, prevId, newId, spaceTime);
      double d = std::hypot(lag.space, stScale * lag.time);
      dist.push_back(std::make_pair(d, nObs + prevId));
    }

    int q = std::min(m, (int) dist.size());
    if(q < (int) dist.size()){
      std::nth_element(dist.begin(), dist.begin() + q, dist.end(),
                       [](const std::pair<double, int> &a,
                          const std::pair<double, int> &b) {
                         return a.first == b.first ? a.second < b.second : a.first < b.first;
                       });
    }
    std::sort(dist.begin(), dist.begin() + q,
              [](const std::pair<double, int> &a,
                 const std::pair<double, int> &b) {
                return a.first == b.first ? a.second < b.second : a.first < b.first;
              });

    NNGPPredictionLocation &location = geometry[pos];
    location.newId = newId;
    location.neighbor.resize(q);
    location.distanceToNew.resize(q);
    for(int h = 0; h < q; h++){
      location.neighbor[h] = dist[h].second;
      int id = location.neighbor[h];
      location.distanceToNew[h] = magcomp::covarianceLag(
        id < nObs ? coordsObs : coordsNew, id < nObs ? nObs : nNew,
        coordsNew, nNew, coordDim, id < nObs ? id : id-nObs, newId, spaceTime);
    }
    location.neighborDistance.reserve((size_t) q * (q > 0 ? q - 1 : 0) / 2);
    for(int c = 0; c < q; c++){
      for(int row = c + 1; row < q; row++){
        int ia = location.neighbor[row];
        int ib = location.neighbor[c];
        const double *ca = ia < nObs ? coordsObs : coordsNew;
        const double *cb = ib < nObs ? coordsObs : coordsNew;
        int na = ia < nObs ? nObs : nNew;
        int nb = ib < nObs ? nObs : nNew;
        int la = ia < nObs ? ia : ia - nObs;
        int lb = ib < nObs ? ib : ib - nObs;
        location.neighborDistance.push_back(magcomp::covarianceLag(ca, na, cb, nb, coordDim, la, lb, spaceTime));
      }
    }
  }
  return geometry;
}

static void nngpPredictDraw(const double *wObs, int nObs,
                            const std::vector<NNGPPredictionLocation> &geometry,
                            const magcomp::CovModel &model, const double *theta, int spaceDim,
                            double sigmaSq, double *out)
{
  std::vector<double> values((size_t) nObs + geometry.size(), NA_REAL);
  std::copy(wObs, wObs + nObs, values.begin());
  std::vector<double> R, chol, r, rhs;

  // Locations remain sequential: a neighbor may be an earlier new location
  // whose effect was sampled in this same posterior draw.
  for(const NNGPPredictionLocation &location : geometry){
    int newId = location.newId;
    const std::vector<int> &neighbor = location.neighbor;
    int q = (int) neighbor.size();

    double mean = 0.0;
    double var = 1.0;
    if(q > 0){
      R.resize((size_t) q * q);
      r.resize(q);

      size_t offset = 0;
      for(int c = 0; c < q; c++){
        R[c + c * q] = 1.0;
        for(int row = c + 1; row < q; row++){
          double value = magcomp::correlation(model, theta, location.neighborDistance[offset++], spaceDim);
          R[row + c * q] = value;
          R[c + row * q] = value;
        }
      }
      for(int h = 0; h < q; h++)
        r[h] = magcomp::correlation(model, theta, location.distanceToNew[h], spaceDim);

      rhs = r;
      choleskyLower(R, chol, q, "NNGP prediction neighbor covariance");
      solveCholLower(chol, q, rhs.data(), 1,
                     "NNGP prediction neighbor solve");
      for(int h = 0; h < q; h++){
        mean += rhs[h] * values[neighbor[h]];
        var -= rhs[h] * r[h];
      }
      if(!R_FINITE(var) || var <= 0.0)
        magcomp::error("NNGP prediction conditional variance is non-positive or "
                 "non-finite at new location %i (q = %i)",
                 newId + 1, q);
    }

    values[nObs + newId] = mean + std::sqrt(sigmaSq * var) *
      Rf_rnorm(0.0, 1.0);
    out[newId] = values[nObs + newId];
  }
}

static SEXP predictMagnitude(SEXP XS, SEXP magnitudeCoefficientsS,
                                              SEXP effectS, SEXP covariateS)
{
  const double *X = realData(XS, "X");
  const double *magnitudeCoefficients = realData(magnitudeCoefficientsS, "magnitude_coefficient_samples");
  std::vector<int> xDim = getDims(XS, "X");
  std::vector<int> gDim = getDims(magnitudeCoefficientsS, "magnitude_coefficient_samples");
  if(xDim.size() != 2 || (gDim.size() != 2 && gDim.size() != 3))
    magcomp::error("invalid dimensions in magnitude prediction");

  int n = xDim[0];
  int p = xDim[1];
  int nSamples = gDim[1];
  int nChains = gDim.size() == 3 ? gDim[2] : 1;
  if(gDim[0] != p)
    magcomp::error("X and magnitude_coefficient_samples have incompatible dimensions");

  const double *effect = NULL;
  std::vector<int> eDim;
  if(!Rf_isNull(effectS)){
    effect = realData(effectS, "effect_samples");
    eDim = getDims(effectS, "effect_samples");
    if((nChains == 1 && (eDim.size() != 2 || eDim[0] != n ||
                         eDim[1] != nSamples)) ||
       (nChains > 1 && (eDim.size() != 3 || eDim[0] != n ||
                        eDim[1] != nSamples || eDim[2] != nChains)))
      magcomp::error("effect_samples has incompatible dimensions");
  }

  const double *covariate = NULL;
  if(!Rf_isNull(covariateS)){
    covariate = realData(covariateS, "covariate");
    if(magcomp::guarded_Rf_length(covariateS) != n)
      magcomp::error("covariate must have one value per prediction row");
  }

  SEXP outS = PROTECT(nChains == 1 ?
    magcomp::guarded_Rf_allocMatrix(REALSXP, n, nSamples) :
    allocRealArray(std::vector<int>{n, nSamples, nChains}));
  double *out = magcomp::guarded_REAL(outS);

  const char *transN = "N";
  double one = 1.0, zero = 0.0;
  for(int ch = 0; ch < nChains; ch++){
    const double *coefficientChain = magnitudeCoefficients + (size_t) ch * p * nSamples;
    double *outCh = out + (size_t) ch * n * nSamples;
    magcomp::rCall([&]() {
      F77_CALL(dgemm)(transN, transN, &n, &nSamples, &p,
                      &one, X, &n, coefficientChain, &p, &zero, outCh, &n
                      FCONE FCONE);
    });
    for(int s = 0; s < nSamples; s++){
      for(int i = 0; i < n; i++){
        double eta = outCh[i + s * n];
        if(effect != NULL){
          double x = covariate == NULL ? 1.0 : covariate[i];
          eta += x * effect[i + s * n + (size_t) ch * n * nSamples];
        }
        outCh[i + s * n] = std::exp(eta);
      }
    }
  }

  UNPROTECT(1);
  return outS;
}

static SEXP predictIidScalar(SEXP effectS,
                                               SEXP observedIndexS,
                                               SEXP newIndexS,
                                               SEXP nNewLevelS,
                                               SEXP priorVarS, bool &rngActive)
{
  const double *effect = realData(effectS, "effect_samples");
  const int *observedIndex = intData(observedIndexS, "observed_index");
  const int *newIndex = intData(newIndexS, "new_index");
  std::vector<int> eDim = getDims(effectS, "effect_samples");
  if(eDim.size() != 2 && eDim.size() != 3)
    magcomp::error("effect_samples must be a matrix or array");

  int nGroup = eDim[0];
  int nSamples = eDim[1];
  int nChains = eDim.size() == 3 ? eDim[2] : 1;
  int n = magcomp::guarded_Rf_length(observedIndexS);
  int nNewLevel = magcomp::guarded_Rf_asInteger(nNewLevelS);
  if(magcomp::guarded_Rf_length(newIndexS) != n || nNewLevel < 0)
    magcomp::error("invalid iid scalar prediction inputs");
  const double *variance = realData(priorVarS, "iid variance samples");
  if(magcomp::guarded_Rf_length(priorVarS) != 1){
    std::vector<int> vDim = getDims(priorVarS, "iid variance samples");
    if((vDim.size() != 2 && vDim.size() != 3) ||
       vDim[0] != 1 || vDim[1] != nSamples ||
       (vDim.size() == 3 && vDim[2] != nChains))
      magcomp::error("iid variance sample dimensions do not match effect samples");
  }
  for(int j = 0; j < magcomp::guarded_Rf_length(priorVarS); j++){
    if(!R_FINITE(variance[j]) || variance[j] <= 0.0)
      magcomp::error("iid variance samples must be positive and finite");
  }

  SEXP outS = PROTECT(nChains == 1 ?
    magcomp::guarded_Rf_allocMatrix(REALSXP, n, nSamples) :
    allocRealArray(std::vector<int>{n, nSamples, nChains}));
  double *out = magcomp::guarded_REAL(outS);
  std::vector<double> drawNew(std::max(1, nNewLevel));

  magcomp::rCall([]() { GetRNGstate(); });
  rngActive = true;
  for(int ch = 0; ch < nChains; ch++){
    for(int s = 0; s < nSamples; s++){
      double sd = std::sqrt(getSampleScalar(priorVarS, s, ch, nSamples, 1.0));
      for(int g = 0; g < nNewLevel; g++)
        drawNew[g] = Rf_rnorm(0.0, sd);
      for(int i = 0; i < n; i++){
        int obs = observedIndex[i];
        int nw = newIndex[i];
        double value;
        if(obs == NA_INTEGER || obs <= 0){
          if(nw == NA_INTEGER || nw <= 0 || nw > nNewLevel)
            magcomp::error("invalid new iid scalar group index");
          value = drawNew[nw - 1];
        }else{
          if(obs > nGroup)
            magcomp::error("invalid observed iid scalar group index");
          value = effect[(obs - 1) + s * nGroup +
                         (size_t) ch * nGroup * nSamples];
        }
        out[i + s * n + (size_t) ch * n * nSamples] = value;
      }
    }
  }

  UNPROTECT(1);
  return outS;
}

static SEXP predictIidFactor(SEXP effectS,
                                               SEXP observedIndexS,
                                               SEXP newIndexS,
                                               SEXP nNewLevelS, bool &rngActive)
{
  const double *effect = realData(effectS, "effect_samples");
  const int *observedIndex = intData(observedIndexS, "observed_index");
  const int *newIndex = intData(newIndexS, "new_index");
  std::vector<int> eDim = getDims(effectS, "effect_samples");
  if(eDim.size() != 3 && eDim.size() != 4)
    magcomp::error("factor effect_samples must be a three- or four-way array");

  int nGroup = eDim[0];
  int K = eDim[1];
  int nSamples = eDim[2];
  int nChains = eDim.size() == 4 ? eDim[3] : 1;
  int n = magcomp::guarded_Rf_length(observedIndexS);
  int nNewLevel = magcomp::guarded_Rf_asInteger(nNewLevelS);
  if(magcomp::guarded_Rf_length(newIndexS) != n || nNewLevel < 0)
    magcomp::error("invalid iid factor prediction inputs");

  SEXP outS = PROTECT(nChains == 1 ?
    allocRealArray(std::vector<int>{n, K, nSamples}) :
    allocRealArray(std::vector<int>{n, K, nSamples, nChains}));
  double *out = magcomp::guarded_REAL(outS);
  std::vector<double> drawNew((size_t) std::max(1, nNewLevel) * K);

  magcomp::rCall([]() { GetRNGstate(); });
  rngActive = true;
  for(int ch = 0; ch < nChains; ch++){
    for(int s = 0; s < nSamples; s++){
      for(int h = 0; h < nNewLevel * K; h++)
        drawNew[h] = Rf_rnorm(0.0, 1.0);
      for(int i = 0; i < n; i++){
        int obs = observedIndex[i];
        int nw = newIndex[i];
        for(int k = 0; k < K; k++){
          double value;
          if(obs == NA_INTEGER || obs <= 0){
            if(nw == NA_INTEGER || nw <= 0 || nw > nNewLevel)
              magcomp::error("invalid new iid factor group index");
            value = drawNew[(nw - 1) + k * nNewLevel];
          }else{
            if(obs > nGroup)
              magcomp::error("invalid observed iid factor group index");
            value = effect[(obs - 1) + k * nGroup +
                           (size_t) s * nGroup * K +
                           (size_t) ch * nGroup * K * nSamples];
          }
          out[i + k * n + (size_t) s * n * K +
              (size_t) ch * n * K * nSamples] = value;
        }
      }
    }
  }

  UNPROTECT(1);
  return outS;
}

static SEXP predictSpatial(SEXP effectS, SEXP coordsObsS,
                                            SEXP coordsNewS, SEXP phiS,
                                            SEXP sigmaS, SEXP mS,
                                            SEXP methodS, SEXP isFactorS, SEXP modelS, SEXP stScaleS, bool &rngActive)
{
  const double *effect = realData(effectS, "effect_samples");
  const double *coordsObs = realData(coordsObsS, "coords_obs");
  const double *coordsNew = realData(coordsNewS, "coords_new");
  std::vector<int> eDim = getDims(effectS, "effect_samples");
  std::vector<int> obsDim = getDims(coordsObsS, "coords_obs");
  std::vector<int> newDim = getDims(coordsNewS, "coords_new");
  if(obsDim.size() != 2 || newDim.size() != 2 ||
     obsDim[1] != newDim[1])
    magcomp::error("prediction coordinates have incompatible dimensions");

  bool isFactor = magcomp::guarded_Rf_asLogical(isFactorS) == TRUE;
  const char *method = CHAR(magcomp::guarded_STRING_ELT(methodS, 0));
  bool denseGp = std::strcmp(method, "gp") == 0;
  bool nngp = std::strcmp(method, "nngp") == 0;
  if(!denseGp && !nngp)
    magcomp::error("spatial prediction method must be 'gp' or 'nngp'");

  const auto &model = magcomp::covarianceModel(CHAR(STRING_ELT(modelS, 0)));
  double stScale = Rf_asReal(stScaleS);
  if(!R_FINITE(stScale) || stScale <= 0) magcomp::error("invalid prediction st_scale");
  int nObs = obsDim[0];
  int nNew = newDim[0];
  int coordDim = obsDim[1];
  int K = 1;
  int nSamples, nChains;
  if(isFactor){
    if(eDim.size() != 3 && eDim.size() != 4)
      magcomp::error("factor spatial effect samples have invalid dimensions");
    K = eDim[1];
    nSamples = eDim[2];
    nChains = eDim.size() == 4 ? eDim[3] : 1;
  }else{
    if(eDim.size() != 2 && eDim.size() != 3)
      magcomp::error("scalar spatial effect samples have invalid dimensions");
    nSamples = eDim[1];
    nChains = eDim.size() == 3 ? eDim[2] : 1;
  }
  if(eDim[0] != nObs)
    magcomp::error("effect_samples and observed coordinates have incompatible rows");

  int m = magcomp::guarded_Rf_asInteger(mS);
  if(m < 1)
    m = nObs + nNew - 1;

  SEXP outS = PROTECT(isFactor ?
    (nChains == 1 ?
       allocRealArray(std::vector<int>{nNew, K, nSamples}) :
       allocRealArray(std::vector<int>{nNew, K, nSamples, nChains})) :
    (nChains == 1 ?
       magcomp::guarded_Rf_allocMatrix(REALSXP, nNew, nSamples) :
       allocRealArray(std::vector<int>{nNew, nSamples, nChains})));
  double *out = magcomp::guarded_REAL(outS);
  std::vector<double> w(nObs);
  std::vector<NNGPPredictionLocation> geometry;
  if(nngp)
    geometry = buildNNGPPredictionGeometry(coordsObs, coordsNew,
                                           nObs, nNew, coordDim, m, model.spaceTime, stScale);

  magcomp::rCall([]() { GetRNGstate(); });
  rngActive = true;
  for(int ch = 0; ch < nChains; ch++){
    for(int s = 0; s < nSamples; s++){
      double sigmaSq = getSampleScalar(sigmaS, s, ch, nSamples, 1.0);
      if(sigmaSq <= 0.0 || !R_FINITE(sigmaSq))
        magcomp::error("spatial prediction sigma_sq must be positive");

      for(int k = 0; k < K; k++){
        std::vector<double> theta(model.size);
        for(int p = 0; p < model.size; ++p)
          theta[p] = getSampleFactorScalar(phiS, k+K*p, s, ch, K*model.size, nSamples, NA_REAL);
        magcomp::validateCovariance(model, theta.data());

        for(int i = 0; i < nObs; i++){
          if(isFactor){
            w[i] = effect[i + k * nObs + (size_t) s * nObs * K +
                          (size_t) ch * nObs * K * nSamples];
          }else{
            w[i] = effect[i + s * nObs +
                          (size_t) ch * nObs * nSamples];
          }
        }

        double *target = isFactor ?
          out + k * nNew + (size_t) s * nNew * K +
            (size_t) ch * nNew * K * nSamples :
          out + s * nNew + (size_t) ch * nNew * nSamples;

        if(denseGp)
          gpPredictDraw(w.data(), coordsObs, coordsNew, nObs, nNew,
                        coordDim, model, theta.data(), sigmaSq, target);
        else
          nngpPredictDraw(w.data(), nObs, geometry, model, theta.data(), coordDim-model.spaceTime, sigmaSq, target);
      }
    }
  }

  UNPROTECT(1);
  return outS;
}

static SEXP predictBlockProportion(SEXP XS, SEXP coefS,
                                                     SEXP loadingS,
                                                     SEXP factorS,
                                                     SEXP covariateS,
                                                     SEXP parentIndexS,
                                                     SEXP JS,
                                                     SEXP linkS)
{
  const double *X = realData(XS, "X");
  const double *coef = realData(coefS, "coef_samples");
  std::vector<int> xDim = getDims(XS, "X");
  std::vector<int> cDim = getDims(coefS, "coef_samples");
  if(xDim.size() != 2 || (cDim.size() != 3 && cDim.size() != 4))
    magcomp::error("invalid dimensions in block proportion prediction");

  int n = xDim[0];
  int p = xDim[1];
  int coefStick = cDim[0];
  int nSamples = cDim[2];
  int nChains = cDim.size() == 4 ? cDim[3] : 1;
  int J = magcomp::guarded_Rf_asInteger(JS);
  int Jstick = J - 1;
  CompositionLink compositionLink = compositionLinkFromSexp(linkS);
  if(cDim[1] != p || J < 2)
    magcomp::error("X, coef_samples, and J have incompatible dimensions");

  bool hasParent = !Rf_isNull(parentIndexS) && magcomp::guarded_Rf_length(parentIndexS) > 0;
  const int *parentIndex = hasParent ? intData(parentIndexS,
                                              "parent_index") : NULL;
  if(hasParent && magcomp::guarded_Rf_length(parentIndexS) != n)
    magcomp::error("parent_index must have one value per prediction row");

  std::vector<FactorPredictionTerm> factorTerms;
  bool hasFactor = !Rf_isNull(loadingS) || !Rf_isNull(factorS);
  if(hasFactor){
    if(TYPEOF(loadingS) == VECSXP || TYPEOF(factorS) == VECSXP){
      if(TYPEOF(loadingS) != VECSXP || TYPEOF(factorS) != VECSXP)
        magcomp::error("loading_samples and factor_samples lists must be supplied together");
      if(magcomp::guarded_Rf_length(loadingS) != magcomp::guarded_Rf_length(factorS))
        magcomp::error("loading_samples and factor_samples lists must have the same length");
      bool covariateList = !Rf_isNull(covariateS) && TYPEOF(covariateS) == VECSXP;
      if(!Rf_isNull(covariateS) && !covariateList)
        magcomp::error("factor_covariate must be a list when factor samples are lists");
      if(covariateList && magcomp::guarded_Rf_length(covariateS) != magcomp::guarded_Rf_length(loadingS))
        magcomp::error("factor_covariate list must match factor sample list length");
      for(int f = 0; f < magcomp::guarded_Rf_length(loadingS); f++){
        SEXP covS = covariateList ? magcomp::guarded_VECTOR_ELT(covariateS, f) : R_NilValue;
        addFactorPredictionTerm(factorTerms,
                                magcomp::guarded_VECTOR_ELT(loadingS, f),
                                magcomp::guarded_VECTOR_ELT(factorS, f),
                                covS, n, coefStick, nSamples, nChains,
                                "composition factor");
      }
    }else{
      addFactorPredictionTerm(factorTerms, loadingS, factorS, covariateS,
                              n, coefStick, nSamples, nChains,
                              "composition factor");
    }
  }

  SEXP outS = PROTECT(nChains == 1 ?
    allocRealArray(std::vector<int>{n, J, nSamples}) :
    allocRealArray(std::vector<int>{n, J, nSamples, nChains}));
  double *out = magcomp::guarded_REAL(outS);

  for(int ch = 0; ch < nChains; ch++){
    for(int s = 0; s < nSamples; s++){
      for(int i = 0; i < n; i++){
        int parent = 0;
        if(hasParent){
          if(parentIndex[i] == NA_INTEGER || parentIndex[i] <= 0)
            magcomp::error("parent_index contains invalid values");
          parent = parentIndex[i] - 1;
        }

        double remaining = 1.0;
        for(int j = 0; j < Jstick; j++){
          int row = hasParent ? parent * Jstick + j : j;
          if(row < 0 || row >= coefStick)
            magcomp::error("parent_index maps outside coefficient rows");

          double eta = 0.0;
          for(int l = 0; l < p; l++){
            eta += X[i + l * n] *
              coef[row + l * coefStick +
                   (size_t) s * coefStick * p +
                   (size_t) ch * coefStick * p * nSamples];
          }
          if(hasFactor){
            for(int f = 0; f < (int) factorTerms.size(); f++){
              const FactorPredictionTerm& term = factorTerms[f];
              double x = term.covariate == NULL ? 1.0 : term.covariate[i];
              double fw = 0.0;
              for(int k = 0; k < term.K; k++){
                fw += term.loading[row + k * coefStick +
                                   (size_t) s * coefStick * term.K +
                                   (size_t) ch * coefStick * term.K *
                                   nSamples] *
                  term.factor[i + k * n +
                              (size_t) s * n * term.K +
                              (size_t) ch * n * term.K * nSamples];
              }
              eta += x * fw;
            }
          }

          double pi = invCompositionLink(eta, compositionLink);
          double value = remaining * pi;
          out[i + j * n + (size_t) s * n * J +
              (size_t) ch * n * J * nSamples] = value;
          remaining *= 1.0 - pi;
        }
        out[i + Jstick * n + (size_t) s * n * J +
            (size_t) ch * n * J * nSamples] = remaining;
      }
    }
  }

  UNPROTECT(1);
  return outS;
}

static SEXP combineConditional(SEXP parentS, SEXP childS)
{
  const double *parent = realData(parentS, "parent_proportion");
  const double *child = realData(childS, "child_proportion");
  std::vector<int> pDim = getDims(parentS, "parent_proportion");
  std::vector<int> cDim = getDims(childS, "child_proportion");
  if((pDim.size() != 3 && pDim.size() != 4) ||
     cDim.size() != pDim.size())
    magcomp::error("parent and child proportions must both be arrays");

  int n = pDim[0];
  int JA = pDim[1];
  int nSamples = pDim[2];
  int nChains = pDim.size() == 4 ? pDim[3] : 1;
  if(cDim[0] != n * JA || cDim[2] != nSamples ||
     (nChains > 1 && cDim[3] != nChains))
    magcomp::error("parent and child proportions have incompatible dimensions");

  int JB = cDim[1];
  SEXP outS = PROTECT(nChains == 1 ?
    allocRealArray(std::vector<int>{n, JA * JB, nSamples}) :
    allocRealArray(std::vector<int>{n, JA * JB, nSamples, nChains}));
  double *out = magcomp::guarded_REAL(outS);

  for(int ch = 0; ch < nChains; ch++){
    for(int s = 0; s < nSamples; s++){
      for(int a = 0; a < JA; a++){
        for(int j = 0; j < JB; j++){
          int outCol = a * JB + j;
          for(int i = 0; i < n; i++){
            int childRow = a * n + i;
            out[i + outCol * n + (size_t) s * n * JA * JB +
                (size_t) ch * n * JA * JB * nSamples] =
              parent[i + a * n + (size_t) s * n * JA +
                     (size_t) ch * n * JA * nSamples] *
              child[childRow + j * n * JA +
                    (size_t) s * n * JA * JB +
                    (size_t) ch * n * JA * JB * nSamples];
          }
        }
      }
    }
  }

  UNPROTECT(1);
  return outS;
}

static SEXP predictCategory(SEXP magnitudeS,
                                             SEXP proportionS,
                                             SEXP areaS,
                                             SEXP typeS)
{
  const double *magnitude = realData(magnitudeS, "magnitude_samples");
  const double *proportion = realData(proportionS, "proportion_samples");
  std::vector<int> mDim = getDims(magnitudeS, "magnitude_samples");
  std::vector<int> pDim = getDims(proportionS, "proportion_samples");
  if((mDim.size() != 2 && mDim.size() != 3) ||
     (pDim.size() != 3 && pDim.size() != 4))
    magcomp::error("invalid category prediction sample dimensions");

  int n = pDim[0];
  int J = pDim[1];
  int nSamples = pDim[2];
  int nChains = pDim.size() == 4 ? pDim[3] : 1;
  if(mDim[0] != n || mDim[1] != nSamples ||
     (nChains > 1 && (mDim.size() != 3 || mDim[2] != nChains)) ||
     (nChains == 1 && mDim.size() != 2))
    magcomp::error("magnitude and proportion samples have incompatible dimensions");

  int type = magcomp::guarded_Rf_asInteger(typeS);
  if(type < 0 || type > 2)
    magcomp::error("prediction type must be 0, 1, or 2");
  const double *area = NULL;
  if(type == 2){
    area = realData(areaS, "area");
    if(magcomp::guarded_Rf_length(areaS) != n)
      magcomp::error("area must have one value per prediction row");
  }

  SEXP outS = PROTECT(nChains == 1 ?
    allocRealArray(std::vector<int>{n, J, nSamples}) :
    allocRealArray(std::vector<int>{n, J, nSamples, nChains}));
  double *out = magcomp::guarded_REAL(outS);

  for(int ch = 0; ch < nChains; ch++){
    for(int s = 0; s < nSamples; s++){
      for(int j = 0; j < J; j++){
        for(int i = 0; i < n; i++){
          double value = proportion[i + j * n + (size_t) s * n * J +
                                    (size_t) ch * n * J * nSamples];
          if(type > 0){
            double scale = magnitude[i + s * n +
                                     (size_t) ch * n * nSamples];
            if(type == 2)
              scale *= area[i];
            value *= scale;
          }
          out[i + j * n + (size_t) s * n * J +
              (size_t) ch * n * J * nSamples] = value;
        }
      }
    }
  }

  UNPROTECT(1);
  return outS;
}

// Keep R conditions inside the unwind bridge so every C++ workspace is freed.
extern "C" SEXP magcomp_predict_magnitude_cpp(SEXP XS, SEXP magnitudeCoefficientsS,
                                              SEXP effectS, SEXP covariateS)
{
  return magcomp::nativeCall([&]() {
    return predictMagnitude(XS, magnitudeCoefficientsS, effectS, covariateS);
  });
}

extern "C" SEXP magcomp_predict_iid_scalar_cpp(SEXP effectS,
                                               SEXP observedIndexS,
                                               SEXP newIndexS,
                                               SEXP nNewLevelS,
                                               SEXP priorVarS)
{
  return magcomp::nativeCallWithRng([&](bool &rngActive) {
    return predictIidScalar(effectS, observedIndexS, newIndexS, nNewLevelS, priorVarS,
      rngActive);
  });
}

extern "C" SEXP magcomp_predict_iid_factor_cpp(SEXP effectS,
                                               SEXP observedIndexS,
                                               SEXP newIndexS,
                                               SEXP nNewLevelS)
{
  return magcomp::nativeCallWithRng([&](bool &rngActive) {
    return predictIidFactor(effectS, observedIndexS, newIndexS, nNewLevelS, rngActive);
  });
}

extern "C" SEXP magcomp_spatial_predict_cpp(SEXP effectS, SEXP coordsObsS,
                                            SEXP coordsNewS, SEXP phiS,
                                            SEXP sigmaS, SEXP mS,
                                            SEXP methodS, SEXP isFactorS, SEXP modelS, SEXP stScaleS)
{
  return magcomp::nativeCallWithRng([&](bool &rngActive) {
    return predictSpatial(effectS, coordsObsS, coordsNewS, phiS, sigmaS, mS, methodS,
      isFactorS, modelS, stScaleS, rngActive);
  });
}

extern "C" SEXP magcomp_predict_block_proportion_cpp(SEXP XS, SEXP coefS,
                                                     SEXP loadingS,
                                                     SEXP factorS,
                                                     SEXP covariateS,
                                                     SEXP parentIndexS,
                                                     SEXP JS,
                                                     SEXP linkS)
{
  return magcomp::nativeCall([&]() {
    return predictBlockProportion(XS, coefS, loadingS, factorS, covariateS,
      parentIndexS, JS, linkS);
  });
}

extern "C" SEXP magcomp_combine_conditional_cpp(SEXP parentS, SEXP childS)
{
  return magcomp::nativeCall([&]() {
    return combineConditional(parentS, childS);
  });
}

extern "C" SEXP magcomp_predict_category_cpp(SEXP magnitudeS,
                                             SEXP proportionS,
                                             SEXP areaS,
                                             SEXP typeS)
{
  return magcomp::nativeCall([&]() {
    return predictCategory(magnitudeS, proportionS, areaS, typeS);
  });
}
