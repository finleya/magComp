#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

#define USE_FC_LEN_T
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>


#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include <Matrix.h>
#include <cholmod.h>
#include <BayesLogit.h>
#include "magcomp_runtime.h"
#include "magcomp_covariance.h"

#ifndef FCONE
# define FCONE
#endif

typedef struct {
  int n;
  int maxNeighbor;
  bool active;
  cholmod_common cm;
  cholmod_sparse *M;
  cholmod_factor *MFactor;
  std::vector<int> diagCache;
  std::vector<int> nngpCacheStart;
  std::vector<int> nngpCacheIndex;
  std::vector<double> rhsWork;
  std::vector<double> zWork;
} SparseNngpPrecision;

typedef struct {
  int n;
  int maxNeighbor;
  bool active;
  std::vector<int> orderIndex;
  std::vector<int> neighborCount;
  std::vector<int> neighborIndex;
  SparseNngpPrecision sparse;
} NngpGraphBackend;

typedef struct {
  int nNngp;
  int reservedNngp;
  std::vector<NngpGraphBackend> nngp;
} LatentBackendRegistry;

typedef struct {
  double samplerStart;
  double setupSeconds;
  double loopSeconds;
  double nbSizeSeconds;
  double magnitudeCoefficientSeconds;
  double magnitudeProcessSeconds;
  double allocationBlockSeconds;
  double allocationFactorSeconds;
  double saveSeconds;
  double pgMagnitudeSeconds;
	  double pgCompositionSeconds;
	  double pgAggregatedCompositionSeconds;
	  double pgOtherSeconds;
	  double pgTotalSeconds;
	  double probitCompositionSeconds;
	  long long pgMagnitudeDraws;
	  long long pgCompositionDraws;
	  long long pgAggregatedCompositionDraws;
	  long long pgOtherDraws;
	  long long pgFallbackDraws;
	  long long probitCompositionDraws;
	  long long iterations;
	  long long retained;
	} SamplerProfile;

static SamplerProfile *gSamplerProfile = NULL;

static double profileNowSeconds()
{
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(
    clock::now().time_since_epoch()
  ).count();
}

typedef struct {
  bool adapt;
  int batchLength;
  int batchPos;
  int batchIndex;
  int batchAccept;
  double targetAccept;
  double lastBatchAccept;
  std::vector<double> batchAcceptance;
  std::vector<double> tuning;
} ScalarMetropolisState;

#include "magcomp_covariance_metropolis.h"

typedef struct {
  std::string name;
  std::string termType;
  bool isIid;
  bool isSpatial;
  bool denseGp;
  bool varying;
  int n;
  int nGroup;
  int coordDim;
  int nNeighbor;
  const double *coords;
  NngpGraphBackend *nngpGraph;
  SEXP orderS;
  SEXP neighborCountS;
  SEXP neighborS;
  std::vector<int> groupIndex;
  std::vector<double> covariate;
  std::vector<double> effect;
  std::vector<double> gpPrecision;
  std::vector<double> gpPrecisionProposal;
  std::vector<double> nngpB;
  std::vector<double> nngpF;
  std::vector<double> nngpBProposal;
  std::vector<double> nngpFProposal;
  double gpLogDet;
  double gpLogDetProposal;
  double sigmaSq;
  double varianceShape;
  double varianceScale;
  double covScale;
  bool estimateCovariance;
  int covAccept;
  int covAttempt;
  CovarianceMetropolis covMetrop;
} MagnitudeProcessState;

typedef struct {
  int n;
  int p;
  const double *X;
  std::vector<int> total;
  std::vector<double> area;
  std::vector<double> coef;
  std::vector<double> omega;
  int nIid;
  int nSpatial;
  std::vector<MagnitudeProcessState> process;
} MagnitudeState;

typedef struct {
  bool estimate;
  double value;
  double shape;
  double rate;
  double tuning;
  int accept;
  int attempt;
  ScalarMetropolisState metrop;
} MagnitudeNbSizeState;

typedef struct {
  std::string name;
  std::string termType;
  bool isIid;
  bool isSpatial;
  bool denseGp;
  bool varying;
  int K;
  int nGroup;
  std::vector<int> groupIndex;
  std::vector<int> anchorIndex;
  std::vector<double> covariate;
  std::vector<double> coordStorage;
  double loadingPriorPrec;
  std::vector<double> loading;
  std::vector<double> factorScore;
  const double *coords;
  int coordDim;
  int nNeighbor;
  NngpGraphBackend *nngpGraph;
  SEXP orderS;
  SEXP neighborCountS;
  SEXP neighborS;
  std::vector< std::vector<double> > gpPrecision;
  std::vector< std::vector<double> > gpPrecisionProposal;
  std::vector< std::vector<double> > nngpB;
  std::vector< std::vector<double> > nngpF;
  std::vector< std::vector<double> > nngpBProposal;
  std::vector< std::vector<double> > nngpFProposal;
  std::vector<double> gpLogDet;
  std::vector<double> gpLogDetProposal;
  std::vector<double> covScale;
  bool estimateCovariance;
  std::vector<int> covAccept;
  std::vector<int> covAttempt;
  std::vector<CovarianceMetropolis> covMetrop;
} CompositionFactorState;

typedef struct {
  int n;
  int J;
  int Jstick;
  int coefStick;
  int p;
  const double *X;
  std::vector<int> count;
  std::vector<int> trials;
  std::vector<int> unitIndex;
  std::vector<int> parentIndex;
	  std::vector<double> coef;
	  std::vector<double> omega;
	  std::vector<double> workingResponse;
	  std::vector<double> proportion;
  std::vector<int> aggregateCount;
  std::vector<int> aggregateTrials;
  bool hasFactor;
  int nFactor;
  int nSpatialFactor;
  bool aggregateInterceptOnly;
  std::vector<CompositionFactorState> factor;
  bool parentSpecificCoef;
  int interceptSmoothOrder;
  int slopeSmoothOrder;
  double interceptPriorPrec;
  double betaPriorPrec;
  double interceptSmoothPrec;
  double slopeSmoothPrec;
} AllocationBlockState;

static SEXP getListElement(SEXP list, const char *name)
{
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);

  if(Rf_isNull(names))
    magcomp::error("expected a named list while looking for '%s'", name);

  for(int i = 0; i < Rf_length(list); i++){
    if(std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
      return VECTOR_ELT(list, i);
  }

  magcomp::error("missing list element '%s'", name);
  return R_NilValue;
}

static SEXP getOptionalListElement(SEXP list, const char *name)
{
  if(Rf_isNull(list))
    return R_NilValue;

  SEXP names = Rf_getAttrib(list, R_NamesSymbol);
  if(Rf_isNull(names))
    return R_NilValue;

  for(int i = 0; i < Rf_length(list); i++){
    if(std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
      return VECTOR_ELT(list, i);
  }

  return R_NilValue;
}

static SEXP getOptionalBlockStart(SEXP startingS, const char *listName,
                                  const char *blockName)
{
  SEXP blockListS = getOptionalListElement(startingS, listName);
  if(Rf_isNull(blockListS))
    return R_NilValue;
  if(TYPEOF(blockListS) != VECSXP || Rf_isNull(Rf_getAttrib(blockListS,
                                                            R_NamesSymbol))){
    magcomp::error("starting$%s must be a named list", listName);
  }
  if(blockName == NULL)
    return R_NilValue;
  return getOptionalListElement(blockListS, blockName);
}

static double getListDouble(SEXP list, const char *name)
{
  SEXP value = getListElement(list, name);
  if(Rf_length(value) < 1)
    magcomp::error("list element '%s' must have length at least one", name);

  return magcomp::guarded_Rf_asReal(value);
}

static void getListDoubleVector(SEXP list,
                                const char *name,
                                int length,
                                const char *context,
                                std::vector<double>& out,
                                bool positive)
{
  SEXP value = getListElement(list, name);
  int valueLength = Rf_length(value);
  if(valueLength != 1 && valueLength != length){
    magcomp::error("%s %s must have length one or K", context, name);
  }

  SEXP tmpS = PROTECT(magcomp::guarded_Rf_coerceVector(value, REALSXP));
  out.assign(length, 0.0);
  for(int i = 0; i < length; i++){
    out[i] = REAL(tmpS)[valueLength == 1 ? 0 : i];
    if(!R_FINITE(out[i]) || (positive && out[i] <= 0.0)){
      UNPROTECT(1);
      magcomp::error("%s %s must contain positive finite values", context, name);
    }
  }
  UNPROTECT(1);
}

static int getListInteger(SEXP list, const char *name)
{
  SEXP value = getListElement(list, name);
  if(Rf_length(value) < 1)
    magcomp::error("list element '%s' must have length at least one", name);

  return magcomp::guarded_Rf_asInteger(value);
}

static bool getListLogical(SEXP list, const char *name)
{
  SEXP value = getListElement(list, name);
  if(Rf_length(value) < 1)
    magcomp::error("list element '%s' must have length at least one", name);

  return magcomp::guarded_Rf_asLogical(value) == TRUE;
}

static const char *getListString(SEXP list, const char *name)
{
  SEXP value = getListElement(list, name);
  if(TYPEOF(value) != STRSXP || Rf_length(value) < 1)
    magcomp::error("list element '%s' must be a character value", name);

  return CHAR(STRING_ELT(value, 0));
}

static void checkArraySize(int d1, int d2, int d3 = 1)
{
  if(d1 < 0 || d2 < 0 || d3 < 0 || (double) d1 * d2 * d3 > INT_MAX)
    magcomp::error("sampler array exceeds the supported integer indexing limit; reduce dimensions or retained samples");
}

static SEXP allocMatrixReal(int nrow, int ncol)
{
  checkArraySize(nrow, ncol);
  SEXP out = magcomp::guarded_Rf_allocMatrix(REALSXP, nrow, ncol);
  double *x = REAL(out);

  for(int i = 0; i < nrow * ncol; i++)
    x[i] = 0.0;

  return out;
}

static SEXP allocArrayReal3(int d1, int d2, int d3)
{
  checkArraySize(d1, d2, d3);
  SEXP dimS = PROTECT(magcomp::guarded_Rf_allocVector(INTSXP, 3));
  INTEGER(dimS)[0] = d1;
  INTEGER(dimS)[1] = d2;
  INTEGER(dimS)[2] = d3;

  SEXP out = magcomp::guarded_Rf_allocArray(REALSXP, dimS);
  double *x = REAL(out);
  for(int i = 0; i < d1 * d2 * d3; i++)
    x[i] = 0.0;

  UNPROTECT(1);
  return out;
}

static double logistic(double x)
{
  double ex;

  if(x >= 0.0)
    return 1.0 / (1.0 + std::exp(-x));

  ex = std::exp(x);
  return ex / (1.0 + ex);
}

enum CompositionLink {
  COMPOSITION_LINK_LOGIT = 0,
  COMPOSITION_LINK_PROBIT = 1
};

static const char *compositionLinkName(CompositionLink link)
{
  return link == COMPOSITION_LINK_PROBIT ? "probit" : "logit";
}

static CompositionLink getCompositionLink(SEXP controlS)
{
  SEXP linkS = getOptionalListElement(controlS, "composition_link");
  if(Rf_isNull(linkS))
    return COMPOSITION_LINK_PROBIT;
  if(TYPEOF(linkS) != STRSXP || Rf_length(linkS) < 1)
    magcomp::error("control$composition_link must be 'logit' or 'probit'");

  const char *link = CHAR(STRING_ELT(linkS, 0));
  if(std::strcmp(link, "logit") == 0)
    return COMPOSITION_LINK_LOGIT;
  if(std::strcmp(link, "probit") == 0)
    return COMPOSITION_LINK_PROBIT;

  magcomp::error("unsupported composition_link '%s'", link);
  return COMPOSITION_LINK_LOGIT;
}

static double compositionInverseLink(double x, CompositionLink link)
{
  if(link == COMPOSITION_LINK_PROBIT)
    return Rf_pnorm5(x, 0.0, 1.0, 1, 0);
  return logistic(x);
}

static double compositionLinkQuantile(double p, CompositionLink link)
{
  if(link == COMPOSITION_LINK_PROBIT)
    return Rf_qnorm5(p, 0.0, 1.0, 1, 0);
  return std::log(p / (1.0 - p));
}

static double rstdnormLowerTrunc(double lower)
{
  double x, alpha, u, logu;

  if(!R_FINITE(lower))
    magcomp::error("non-finite probit truncation point");

  if(lower <= 0.0){
    do {
      x = norm_rand();
    } while(x < lower);
    return x;
  }

  alpha = 0.5 * (lower + std::sqrt(lower * lower + 4.0));
  for(;;){
    u = unif_rand();
    if(u <= 0.0)
      continue;
    x = lower - std::log(u) / alpha;
    logu = std::log(unif_rand());
    if(logu <= -0.5 * (x - alpha) * (x - alpha))
      return x;
  }
}

static double rnormProbitLatent(double mean, int y)
{
  double lower;

  if(!R_FINITE(mean))
    magcomp::error("non-finite probit linear predictor");

  lower = -mean;
  if(y)
    return mean + rstdnormLowerTrunc(lower);

  return mean - rstdnormLowerTrunc(-lower);
}

// N(mean, sd^2) restricted to (0, infinity). In the positive tail, draw
// the excess above the truncation point directly to avoid cancellation.
static double rnormPositive(double mean, double sd)
{
  double lower = -mean / sd;
  if(!R_FINITE(mean) || !R_FINITE(sd) || sd <= 0.0 || !R_FINITE(lower))
    magcomp::error("invalid positive-normal parameters");

  if(lower <= 0.0){
    double value;
    do {
      value = mean + sd * norm_rand();
    } while(value <= 0.0);
    return value;
  }

  double alpha = 0.5 * lower + 0.5 * std::hypot(lower, 2.0);
  for(;;){
    double u = unif_rand();
    if(u <= 0.0)
      continue;
    double excess = -std::log(u) / alpha;
    double difference = excess - 1.0 / alpha;
    if(std::log(unif_rand()) <= -0.5 * difference * difference){
      double value = sd * excess;
      if(!R_FINITE(value) || value <= 0.0)
        magcomp::error("positive-normal draw is outside the representable range");
      return value;
    }
  }
}

static double drawProbitBinomialLatentSum(double eta,
                                          int successes,
                                          int trials,
                                          const char *context,
                                          int iter)
{
  if(trials <= 0)
    return 0.0;
  if(successes < 0 || successes > trials){
    magcomp::error("invalid probit binomial counts in %s at iteration %i",
             context, iter);
  }

  double out = 0.0;
  for(int k = 0; k < successes; k++)
    out += rnormProbitLatent(eta, 1);
  for(int k = successes; k < trials; k++)
    out += rnormProbitLatent(eta, 0);

  return out;
}

static double drawPolyaGamma(BayesLogit_rpg_hybrid_t rpg,
                             BayesLogit_rpg_gamma_t rpgGamma,
                             double h,
                             double z,
                             const char *context,
                             int iter)
{
  if(!R_FINITE(h) || h <= 0.0 || !R_FINITE(z)){
    magcomp::error("invalid Polya-Gamma parameters in %s at iteration %i "
             "(h = %.17g, z = %.17g)",
	             context, iter, h, z);
	  }

  double profileStart = profileNowSeconds();
	  double omega = magcomp::rCall([&]() { return rpg(h, z); });
  bool usedFallback = false;
	  if(!R_FINITE(omega) || omega <= 0.0){
    usedFallback = true;
	    omega = magcomp::rCall([&]() { return rpgGamma(h, z, 1000); });
  }

  if(gSamplerProfile != NULL){
    double elapsed = profileNowSeconds() - profileStart;
    gSamplerProfile->pgTotalSeconds += elapsed;
    if(usedFallback)
      gSamplerProfile->pgFallbackDraws++;
    if(std::strcmp(context, "magnitude update") == 0){
      gSamplerProfile->pgMagnitudeDraws++;
      gSamplerProfile->pgMagnitudeSeconds += elapsed;
    }else if(std::strcmp(context, "composition update") == 0){
      gSamplerProfile->pgCompositionDraws++;
      gSamplerProfile->pgCompositionSeconds += elapsed;
    }else if(std::strcmp(context, "aggregated composition update") == 0){
      gSamplerProfile->pgAggregatedCompositionDraws++;
      gSamplerProfile->pgAggregatedCompositionSeconds += elapsed;
    }else{
      gSamplerProfile->pgOtherDraws++;
      gSamplerProfile->pgOtherSeconds += elapsed;
    }
  }

	  if(!R_FINITE(omega) || omega <= 0.0){
	    magcomp::error("invalid Polya-Gamma draw in %s at iteration %i "
	             "(h = %.17g, z = %.17g)",
             context, iter, h, z);
  }

  return omega;
}

static double logit(double p)
{
  return std::log(p / (1.0 - p));
}

static double clipProbability(double p)
{
  if(p < 1.0e-8)
    return 1.0e-8;
  if(p > 1.0 - 1.0e-8)
    return 1.0 - 1.0e-8;

  return p;
}

static int getMatrixCount(SEXP y, int idx)
{
  double value;
  int out;

  if(TYPEOF(y) == INTSXP){
    out = INTEGER(y)[idx];
  }else if(TYPEOF(y) == REALSXP){
    value = REAL(y)[idx];
    out = (int) std::floor(value + 0.5);
  }else{
    magcomp::error("y must be an integer or numeric matrix");
    out = 0;
  }

  if(out < 0)
    magcomp::error("y must contain non-negative counts");

  return out;
}

static double rowDot(const double *X, const double *beta, int n, int p, int i)
{
  double out = 0.0;

  for(int l = 0; l < p; l++)
    out += X[i + l * n] * beta[l];

  return out;
}

static int blockCoefIndex(const AllocationBlockState& block, int i, int j)
{
  if(block.parentSpecificCoef)
    return block.parentIndex[i] * block.Jstick + j;

  return j;
}

static bool loadingIsConstrained(int j,
                                 int k,
                                 const std::vector<int>& anchorIndex,
                                 double& value)
{
  for(int h = 0; h < k; h++){
    if(j == anchorIndex[h]){
      value = 0.0;
      return true;
    }
  }

  value = 0.0;
  return false;
}

static bool loadingIsFree(int j, int k, const std::vector<int>& anchorIndex)
{
  double value = 0.0;
  return !loadingIsConstrained(j, k, anchorIndex, value);
}

static void enforceLoadingConstraints(std::vector<double>& loading,
                                      int Jstick,
                                      int K,
                                      const std::vector<int>& anchorIndex)
{
  for(int k = 0; k < K; k++){
    for(int h = 0; h < k; h++)
      loading[anchorIndex[h] + k * Jstick] = 0.0;
    double diagonal = loading[anchorIndex[k] + k * Jstick];
    if(!R_FINITE(diagonal) || diagonal <= 0.0)
      magcomp::error("anchor diagonal loadings must be positive and finite");
  }
}

static double compositionFactorEta(const AllocationBlockState& block,
                                   const CompositionFactorState& factor,
                                   int i,
                                   int j)
{
  double out = 0.0;

  int g = factor.groupIndex[i];
  int coefIdx = blockCoefIndex(block, i, j);
  for(int k = 0; k < factor.K; k++){
    out += factor.loading[coefIdx + k * block.coefStick] *
      factor.factorScore[g + k * factor.nGroup];
  }

  return factor.covariate[i] * out;
}

static double blockFactorEtaExcept(const AllocationBlockState& block,
                                   int i,
                                   int j,
                                   int excludeFactor)
{
  double out = 0.0;

  if(!block.hasFactor)
    return out;

  for(int f = 0; f < (int) block.factor.size(); f++){
    if(f == excludeFactor)
      continue;
    out += compositionFactorEta(block, block.factor[f], i, j);
  }

  return out;
}

static double blockFactorEta(const AllocationBlockState& block,
                             int i,
                             int j)
{
  return blockFactorEtaExcept(block, i, j, -1);
}

static double blockLinearPredictor(const AllocationBlockState& block,
                                   int i,
                                   int j)
{
  double out = 0.0;

  for(int l = 0; l < block.p; l++)
    out += block.X[i + l * block.n] *
      block.coef[blockCoefIndex(block, i, j) + l * block.coefStick];

  out += blockFactorEta(block, i, j);

  return out;
}

static int smoothOrder(const char *method, const char *name)
{
  if(std::strcmp(method, "none") == 0)
    return 0;
  if(std::strcmp(method, "rw1") == 0)
    return 1;
  if(std::strcmp(method, "rw2") == 0)
    return 2;

  magcomp::error("unsupported %s smoothing method '%s'", name, method);
  return 0;
}

static void addCurvePrior(double *Q,
                          int dim,
                          int order,
                          double basePrec,
                          double smoothPrec)
{
  for(int j = 0; j < dim; j++)
    Q[j + j * dim] += basePrec;

  if(order == 0)
    return;

  if(order == 1){
    for(int d = 0; d < dim - 1; d++){
      Q[d + d * dim] += smoothPrec;
      Q[(d + 1) + (d + 1) * dim] += smoothPrec;
      Q[d + (d + 1) * dim] -= smoothPrec;
      Q[(d + 1) + d * dim] -= smoothPrec;
    }
    return;
  }

  for(int d = 0; d < dim - 2; d++){
    int id[3] = {d, d + 1, d + 2};
    double val[3] = {1.0, -2.0, 1.0};

    for(int a = 0; a < 3; a++){
      for(int b = 0; b < 3; b++)
        Q[id[a] + id[b] * dim] += smoothPrec * val[a] * val[b];
    }
  }
}

static void addCurvePriorBlock(double *Q,
                               int offset,
                               int dim,
                               int leadingDim,
                               int order,
                               double basePrec,
                               double smoothPrec)
{
  for(int j = 0; j < dim; j++)
    Q[(offset + j) + (offset + j) * leadingDim] += basePrec;

  if(order == 0)
    return;

  if(order == 1){
    for(int d = 0; d < dim - 1; d++){
      int a = offset + d;
      int b = offset + d + 1;
      Q[a + a * leadingDim] += smoothPrec;
      Q[b + b * leadingDim] += smoothPrec;
      Q[a + b * leadingDim] -= smoothPrec;
      Q[b + a * leadingDim] -= smoothPrec;
    }
    return;
  }

  for(int d = 0; d < dim - 2; d++){
    int id[3] = {offset + d, offset + d + 1, offset + d + 2};
    double val[3] = {1.0, -2.0, 1.0};

    for(int a = 0; a < 3; a++){
      for(int b = 0; b < 3; b++)
        Q[id[a] + id[b] * leadingDim] += smoothPrec * val[a] * val[b];
    }
  }
}

static void sampleGaussianPrecision(double *draw,
                                    double *precision,
                                    double *rhs,
                                    double *chol,
                                    double *z,
                                    int dim,
                                    const char *context,
                                    int iter)
{
  int info, nrhs = 1;
  const int inc = 1;
  const double one = 1.0;
  char const *lower = "L";
  char const *trans = "T";
  char const *nonUnit = "N";

  for(int i = 0; i < dim * dim; i++){
    if(!R_FINITE(precision[i]))
      magcomp::error("non-finite precision in %s update at iteration %i",
               context, iter);
    chol[i] = precision[i];
  }

  F77_CALL(dpotrf)(lower, &dim, chol, &dim, &info FCONE);
  if(info != 0)
    magcomp::error("dpotrf failed in %s update at iteration %i", context, iter);

  for(int i = 0; i < dim; i++){
    if(!R_FINITE(rhs[i]))
      magcomp::error("non-finite rhs in %s update at iteration %i", context, iter);
    draw[i] = rhs[i];
  }

  F77_CALL(dpotrs)(lower, &dim, &nrhs, chol, &dim, draw, &dim, &info FCONE);
  if(info != 0)
    magcomp::error("dpotrs failed in %s update at iteration %i", context, iter);

  for(int i = 0; i < dim; i++)
    z[i] = Rf_rnorm(0.0, 1.0);

  F77_CALL(dtrsv)(lower, trans, nonUnit, &dim, chol, &dim, z, &inc
                  FCONE FCONE FCONE);
  F77_CALL(daxpy)(&dim, &one, z, &inc, draw, &inc);
}

static void buildGpPrecision(const double *coords,
                             int n,
                             int coordDim,
                             const CovarianceState& covariance,
                             std::vector<double>& precision,
                             double *logDet)
{
  int info;
  char const *lower = "L";

  if(n < 1)
    magcomp::error("dense GP precision requires at least one location");
  if(coordDim < 1)
    magcomp::error("dense GP precision requires at least one coordinate dimension");
  if((int) precision.size() < n * n)
    magcomp::error("dense GP precision workspace is too small");

  std::fill(precision.begin(), precision.end(), 0.0);
  for(int j = 0; j < n; j++){
    precision[j + j * n] = 1.0;
    for(int i = j + 1; i < n; i++){
      double value = magcomp::correlationAt(*covariance.model, covariance.theta.data(), coords, n, coordDim, i, j);
      precision[i + j * n] = value;
      precision[j + i * n] = value;
    }
  }

  F77_CALL(dpotrf)(lower, &n, precision.data(), &n, &info FCONE);
  if(info != 0)
    magcomp::error("dense GP correlation matrix is not positive definite");

  *logDet = 0.0;
  for(int i = 0; i < n; i++){
    double diag = precision[i + i * n];
    if(!R_FINITE(diag) || diag <= 0.0)
      magcomp::error("invalid dense GP Cholesky diagonal");
    *logDet += 2.0 * std::log(diag);
  }

  F77_CALL(dpotri)(lower, &n, precision.data(), &n, &info FCONE);
  if(info != 0)
    magcomp::error("dense GP precision inversion failed");

  for(int j = 0; j < n; j++){
    for(int i = j + 1; i < n; i++)
      precision[j + i * n] = precision[i + j * n];
  }
}

static void addDensePrecision(double *Q,
                              int n,
                              const std::vector<double>& precision,
                              double scale)
{
  for(int j = 0; j < n; j++){
    for(int i = 0; i < n; i++)
      Q[i + j * n] += scale * precision[i + j * n];
  }
}

static double denseVectorScaledSsq(const std::vector<double>& value,
                                   int n,
                                   const std::vector<double>& precision)
{
  double out = 0.0;

  for(int j = 0; j < n; j++){
    double tmp = 0.0;
    for(int i = 0; i < n; i++)
      tmp += precision[i + j * n] * value[i];
    out += value[j] * tmp;
  }

  return out;
}

static double gpVectorLogPrior(const std::vector<double>& value,
                               int n,
                               const std::vector<double>& precision,
                               double logDet,
                               double sigmaSq)
{
  if(!R_FINITE(sigmaSq) || sigmaSq <= 0.0)
    magcomp::error("invalid dense GP marginal variance");

  return -0.5 * ((double) n * std::log(sigmaSq) + logDet +
                 denseVectorScaledSsq(value, n, precision) / sigmaSq);
}

static double gpFactorColumnLogPrior(const std::vector<double>& factorScore,
                                     int n,
                                     int k,
                                     const std::vector<double>& precision,
                                     double logDet)
{
  double ssq = 0.0;

  for(int j = 0; j < n; j++){
    double tmp = 0.0;
    for(int i = 0; i < n; i++)
      tmp += precision[i + j * n] * factorScore[i + k * n];
    ssq += factorScore[j + k * n] * tmp;
  }

  return -0.5 * (logDet + ssq);
}

static void initScalarMetropolisState(ScalarMetropolisState& metrop)
{
  metrop.adapt = false;
  metrop.batchLength = 25;
  metrop.batchPos = 0;
  metrop.batchIndex = 0;
  metrop.batchAccept = 0;
  metrop.targetAccept = 0.44;
  metrop.lastBatchAccept = NA_REAL;
}

static void initMagnitudeNbSizeState(MagnitudeNbSizeState& nb)
{
  nb.estimate = false;
  nb.value = 1000.0;
  nb.shape = 1.0;
  nb.rate = 1.0;
  nb.tuning = 0.25;
  nb.accept = 0;
  nb.attempt = 0;
  initScalarMetropolisState(nb.metrop);
}

static void prepareScalarMetropolisState(ScalarMetropolisState& metrop,
                                         bool adapt,
                                         int nBurn,
                                         int batchLength,
                                         double targetAccept)
{
  int maxBatchHistory;

  metrop.adapt = adapt && nBurn > 0;
  metrop.batchLength = batchLength;
  metrop.batchPos = 0;
  metrop.batchIndex = 0;
  metrop.batchAccept = 0;
  metrop.targetAccept = targetAccept > 0 ? targetAccept : 0.44;
  metrop.lastBatchAccept = NA_REAL;

  maxBatchHistory = metrop.adapt ? nBurn / batchLength : 0;
  metrop.batchAcceptance.assign(std::max(1, maxBatchHistory), NA_REAL);
  metrop.tuning.assign(std::max(1, maxBatchHistory), NA_REAL);
}

static void updateScalarMetropolisTuning(ScalarMetropolisState& metrop,
                                         int accepted,
                                         double *tuning)
{
  double acceptRate, step;

  if(!metrop.adapt)
    return;

  metrop.batchAccept += accepted;
  metrop.batchPos++;
  if(metrop.batchPos < metrop.batchLength)
    return;

  metrop.batchIndex++;
  acceptRate = (double) metrop.batchAccept / (double) metrop.batchLength;
  metrop.lastBatchAccept = acceptRate;

  step = std::min(0.01, 1.0 / std::sqrt((double) metrop.batchIndex));
  if(acceptRate > metrop.targetAccept)
    *tuning = std::exp(std::log(*tuning) + step);
  else
    *tuning = std::exp(std::log(*tuning) - step);
  if(!R_FINITE(*tuning) || *tuning <= 0.0)
    magcomp::error("invalid adaptive Metropolis tuning value");

  if(metrop.batchIndex <= (int) metrop.batchAcceptance.size()){
    metrop.batchAcceptance[metrop.batchIndex - 1] = acceptRate;
    metrop.tuning[metrop.batchIndex - 1] = *tuning;
  }

  metrop.batchPos = 0;
  metrop.batchAccept = 0;
}

static SEXP scalarMetropolisHistorySexp(const ScalarMetropolisState& metrop,
                                        const std::vector<double>& history)
{
  int n = metrop.batchIndex;
  SEXP out = magcomp::guarded_Rf_allocVector(REALSXP, n);

  for(int i = 0; i < n; i++)
    REAL(out)[i] = history[i];

  return out;
}

static void appendSparseLowerPattern(std::vector<std::vector<int> >& cols,
                                     int row,
                                     int col)
{
  if(row < col)
    std::swap(row, col);

  cols[col].push_back(row);
}

static int findSparseLowerEntry(cholmod_sparse *M, int row, int col)
{
  int *Mp, *Mi;
  int left, right, mid;

  if(row < col)
    std::swap(row, col);

  Mp = (int*) M->p;
  Mi = (int*) M->i;

  left = Mp[col];
  right = Mp[col + 1] - 1;
  while(left <= right){
    mid = left + (right - left) / 2;
    if(Mi[mid] == row)
      return mid;
    if(Mi[mid] < row)
      left = mid + 1;
    else
      right = mid - 1;
  }

  magcomp::error("sparse precision entry (%i,%i) is missing from the pattern",
           row + 1, col + 1);
  return -1;
}

static void sparseCholmodErrorHandler(int, const char *, int, const char *)
{
  // Inspect CHOLMOD's status after it returns, so its own cleanup completes
  // before raising a C++ exception in the caller.
}

static void initSparseNngpPrecisionState(SparseNngpPrecision& sp)
{
  sp.n = 0;
  sp.maxNeighbor = 0;
  sp.active = false;
  sp.M = NULL;
  sp.MFactor = NULL;
}

static void initNngpGraphBackend(NngpGraphBackend& graph)
{
  graph.n = 0;
  graph.maxNeighbor = 0;
  graph.active = false;
  initSparseNngpPrecisionState(graph.sparse);
}

static void initLatentBackendRegistry(LatentBackendRegistry& registry)
{
  registry.nNngp = 0;
  registry.reservedNngp = 0;
  registry.nngp.clear();
}

static void reserveLatentBackendRegistry(LatentBackendRegistry& registry,
                                         int nNngp)
{
  if(nNngp < 0)
    magcomp::error("internal error: invalid NNGP graph reserve count");
  registry.reservedNngp = nNngp;
  registry.nngp.reserve(std::max(1, nNngp));
}

static void freeSparseNngpPrecision(SparseNngpPrecision& sp)
{
  if(!sp.active)
    return;

  if(sp.MFactor != NULL)
    M_cholmod_free_factor(&sp.MFactor, &sp.cm);
  if(sp.M != NULL)
    M_cholmod_free_sparse(&sp.M, &sp.cm);
  M_cholmod_finish(&sp.cm);

  sp.M = NULL;
  sp.MFactor = NULL;
  sp.active = false;
}

static void freeNngpGraphBackend(NngpGraphBackend& graph)
{
  freeSparseNngpPrecision(graph.sparse);
  graph.active = false;
}

static void freeLatentBackendRegistry(LatentBackendRegistry& registry)
{
  for(auto &graph : registry.nngp)
    freeNngpGraphBackend(graph);
  registry.nNngp = 0;
  registry.reservedNngp = 0;
  registry.nngp.clear();
}

struct LatentBackendGuard {
  LatentBackendRegistry &registry;
  ~LatentBackendGuard() { freeLatentBackendRegistry(registry); }
};

static void loadNngpGraph(const char *termName,
                          int n,
                          int maxNeighbor,
                          SEXP orderS,
                          SEXP neighborCountS,
                          SEXP neighborS,
                          std::vector<int>& orderIndex,
                          std::vector<int>& neighborCount,
                          std::vector<int>& neighborIndex)
{
  SEXP neighborDimS;
  std::vector<int> orderSeen(n, 0);
  std::vector<int> orderPosition(n, -1);
  int pos, i, q, h;

  if(!Rf_isMatrix(neighborS))
    magcomp::error("%s neighbor_index must be a matrix", termName);
  neighborDimS = Rf_getAttrib(neighborS, R_DimSymbol);
  if(INTEGER(neighborDimS)[0] != n)
    magcomp::error("%s neighbor_index must have one row per ordered location",
             termName);
  if(INTEGER(neighborDimS)[1] != maxNeighbor)
    magcomp::error("%s neighbor_index has an inconsistent neighbor dimension",
             termName);
  if(Rf_length(orderS) != n || Rf_length(neighborCountS) != n)
    magcomp::error("%s order and neighbor_count must have length n", termName);

  for(pos = 0; pos < n; pos++){
    i = INTEGER(orderS)[pos] - 1;
    if(i < 0 || i >= n)
      magcomp::error("%s order index is out of range", termName);
    if(orderSeen[i])
      magcomp::error("%s order indices must be unique", termName);
    orderSeen[i] = 1;
    orderPosition[i] = pos;
    orderIndex[pos] = i;
  }

  for(pos = 0; pos < n; pos++){
    q = INTEGER(neighborCountS)[pos];
    if(q < 0 || q > maxNeighbor || q > pos)
      magcomp::error("%s neighbor_count is invalid", termName);
    neighborCount[pos] = q;
    for(h = 0; h < q; h++){
      i = INTEGER(neighborS)[pos + h * n] - 1;
      if(i < 0 || i >= n)
        magcomp::error("%s neighbor index is out of range", termName);
      if(orderPosition[i] >= pos)
        magcomp::error("%s neighbor must precede current location", termName);
      neighborIndex[pos + h * n] = i;
    }
  }
}

static void initSparseNngpPrecision(SparseNngpPrecision& sp,
                                    int n,
                                    int maxNeighbor,
                                    const std::vector<int>& orderIndex,
                                    const std::vector<int>& neighborCount,
                                    const std::vector<int>& neighborIndex)
{
  int i, pos, q, h, a, b, col, row, nz, cachePos;
  int supportSize, nodeA, nodeB;
  std::vector<std::vector<int> > cols(n);
  std::vector<int> nodes(std::max(1, maxNeighbor + 1), 0);
  int *Mp, *Mi;
  double *Mx;

  if(n < 1)
    magcomp::error("sparse NNGP precision requires at least one location");
  if(maxNeighbor < 0)
    magcomp::error("sparse NNGP precision received an invalid neighbor count");

  freeSparseNngpPrecision(sp);

  for(i = 0; i < n; i++)
    appendSparseLowerPattern(cols, i, i);

  for(pos = 0; pos < n; pos++){
    q = neighborCount[pos];
    if(q < 0 || q > maxNeighbor)
      magcomp::error("invalid neighbor count while building sparse NNGP pattern");

    nodes[0] = orderIndex[pos];
    for(h = 0; h < q; h++)
      nodes[h + 1] = neighborIndex[pos + h * n];

    supportSize = q + 1;
    for(a = 0; a < supportSize; a++){
      nodeA = nodes[a];
      for(b = 0; b <= a; b++){
        nodeB = nodes[b];
        appendSparseLowerPattern(cols, nodeA, nodeB);
      }
    }
  }

  nz = 0;
  for(col = 0; col < n; col++){
    std::sort(cols[col].begin(), cols[col].end());
    cols[col].erase(std::unique(cols[col].begin(), cols[col].end()),
                    cols[col].end());
    nz += (int) cols[col].size();
  }

  M_R_cholmod_start(&sp.cm);
  sp.cm.final_ll = TRUE;
  sp.cm.error_handler = sparseCholmodErrorHandler;
  sp.active = true;
  sp.M = M_cholmod_allocate_sparse(n, n, nz, 1, 1, -1, CHOLMOD_REAL,
                                   &sp.cm);
  if(sp.M == NULL || sp.cm.status != CHOLMOD_OK)
    magcomp::error("failed to allocate sparse NNGP precision matrix");

  Mp = (int*) sp.M->p;
  Mi = (int*) sp.M->i;
  Mx = (double*) sp.M->x;

  nz = 0;
  Mp[0] = 0;
  for(col = 0; col < n; col++){
    for(row = 0; row < (int) cols[col].size(); row++){
      Mi[nz] = cols[col][row];
      Mx[nz] = 0.0;
      nz++;
    }
    Mp[col + 1] = nz;
  }

  sp.n = n;
  sp.maxNeighbor = maxNeighbor;
  sp.diagCache.assign(n, 0);
  for(i = 0; i < n; i++)
    sp.diagCache[i] = findSparseLowerEntry(sp.M, i, i);

  sp.nngpCacheStart.assign(n + 1, 0);
  for(pos = 0; pos < n; pos++){
    q = neighborCount[pos];
    sp.nngpCacheStart[pos + 1] =
      sp.nngpCacheStart[pos] + 1 + q + q * (q + 1) / 2;
  }
  sp.nngpCacheIndex.assign(sp.nngpCacheStart[n], 0);
  cachePos = 0;
  for(pos = 0; pos < n; pos++){
    q = neighborCount[pos];
    i = orderIndex[pos];

    sp.nngpCacheIndex[cachePos++] = findSparseLowerEntry(sp.M, i, i);
    for(a = 0; a < q; a++){
      nodeA = neighborIndex[pos + a * n];
      sp.nngpCacheIndex[cachePos++] = findSparseLowerEntry(sp.M, i, nodeA);
    }
    for(a = 0; a < q; a++){
      nodeA = neighborIndex[pos + a * n];
      for(b = 0; b <= a; b++){
        nodeB = neighborIndex[pos + b * n];
        sp.nngpCacheIndex[cachePos++] =
          findSparseLowerEntry(sp.M, nodeA, nodeB);
      }
    }
  }

  sp.rhsWork.assign(n, 0.0);
  sp.zWork.assign(n, 0.0);

  sp.MFactor = M_cholmod_analyze(sp.M, &sp.cm);
  if(sp.MFactor == NULL || sp.cm.status != CHOLMOD_OK)
    magcomp::error("CHOLMOD analyze failed for sparse NNGP precision");
}

static bool nngpGraphMatches(const NngpGraphBackend& graph,
                             int n,
                             int maxNeighbor,
                             const std::vector<int>& orderIndex,
                             const std::vector<int>& neighborCount,
                             const std::vector<int>& neighborIndex)
{
  if(!graph.active)
    return false;
  if(graph.n != n || graph.maxNeighbor != maxNeighbor)
    return false;
  if(graph.orderIndex != orderIndex)
    return false;
  if(graph.neighborCount != neighborCount)
    return false;
  if(graph.neighborIndex != neighborIndex)
    return false;

  return true;
}

static NngpGraphBackend *registerNngpGraphBackend(
  LatentBackendRegistry& registry,
  const char *termName,
  int n,
  int maxNeighbor,
  SEXP orderS,
  SEXP neighborCountS,
  SEXP neighborS)
{
  std::vector<int> orderIndex(n, 0);
  std::vector<int> neighborCount(n, 0);
  std::vector<int> neighborIndex(std::max(1, n * std::max(1, maxNeighbor)),
                                 0);
  NngpGraphBackend *graph;

  loadNngpGraph(termName, n, maxNeighbor, orderS, neighborCountS, neighborS,
                orderIndex, neighborCount, neighborIndex);

  for(int g = 0; g < registry.nNngp; g++){
    if(nngpGraphMatches(registry.nngp[g], n, maxNeighbor, orderIndex,
                        neighborCount, neighborIndex))
      return &registry.nngp[g];
  }

  if(registry.nNngp >= registry.reservedNngp)
    magcomp::error("internal error: NNGP graph registry reserve exceeded");

  registry.nngp.push_back(NngpGraphBackend());
  graph = &registry.nngp[registry.nNngp];
  initNngpGraphBackend(*graph);
  graph->n = n;
  graph->maxNeighbor = maxNeighbor;
  graph->orderIndex.swap(orderIndex);
  graph->neighborCount.swap(neighborCount);
  graph->neighborIndex.swap(neighborIndex);
  initSparseNngpPrecision(graph->sparse, n, maxNeighbor, graph->orderIndex,
                          graph->neighborCount, graph->neighborIndex);
  graph->active = true;
  registry.nNngp++;

  return graph;
}

static int choleskySolveSmallInPlace(double *A, double *rhs, int n)
{
  for(int j = 0; j < n; j++){
    double diag = A[j + j * n];
    for(int k = 0; k < j; k++){
      double ljk = A[j + k * n];
      diag -= ljk * ljk;
    }
    if(!R_FINITE(diag) || diag <= 0.0)
      return j + 1;
    A[j + j * n] = std::sqrt(diag);

    for(int i = j + 1; i < n; i++){
      double value = A[i + j * n];
      for(int k = 0; k < j; k++)
        value -= A[i + k * n] * A[j + k * n];
      value /= A[j + j * n];
      if(!R_FINITE(value))
        return j + 1;
      A[i + j * n] = value;
    }
  }

  for(int i = 0; i < n; i++){
    double value = rhs[i];
    for(int k = 0; k < i; k++)
      value -= A[i + k * n] * rhs[k];
    value /= A[i + i * n];
    if(!R_FINITE(value))
      return i + 1;
    rhs[i] = value;
  }

  for(int i = n - 1; i >= 0; i--){
    double value = rhs[i];
    for(int k = i + 1; k < n; k++)
      value -= A[k + i * n] * rhs[k];
    value /= A[i + i * n];
    if(!R_FINITE(value))
      return i + 1;
    rhs[i] = value;
  }

  return 0;
}

static int buildNngpConditionalOne(const double *coords,
                                   int n,
                                   int coordDim,
                                   const CovarianceState& covariance,
                                   int pos,
                                   const std::vector<int>& orderIndex,
                                   const std::vector<int>& neighborCount,
                                   const std::vector<int>& neighborIndex,
                                   std::vector<double>& nngpB,
                                   std::vector<double>& nngpF,
                                   double *C,
                                   double *rhs,
                                   double *c)
{
  int info = 1;

  int q = neighborCount[pos];
  int i = orderIndex[pos];
  if(q <= 0){
    nngpF[pos] = 1.0;
    return 0;
  }

  for(int a = 0; a < q; a++){
    int na = neighborIndex[pos + a * n];
    c[a] = magcomp::correlationAt(*covariance.model, covariance.theta.data(), coords, n, coordDim, i, na);
  }

  for(int a = 0; a < q; a++){
    int na = neighborIndex[pos + a * n];
    rhs[a] = c[a];
    for(int b = 0; b < q; b++){
      int nb = neighborIndex[pos + b * n];
      C[a + b * q] = magcomp::correlationAt(*covariance.model, covariance.theta.data(), coords, n, coordDim, na, nb);
    }
  }
  info = choleskySolveSmallInPlace(C, rhs, q);
  if(info != 0)
    return pos + 1;

  double cSolve = 0.0;
  for(int a = 0; a < q; a++){
    nngpB[pos + a * n] = rhs[a];
    cSolve += c[a] * rhs[a];
  }

  nngpF[pos] = 1.0 - cSolve;
  if(!R_FINITE(nngpF[pos]) || nngpF[pos] <= 0.0)
    return -(pos + 1);

  return 0;
}

static void buildNngpConditionals(const double *coords,
                                  int n,
                                  int coordDim,
                                  int m,
                                  const CovarianceState& covariance,
                                  const std::vector<int>& orderIndex,
                                  const std::vector<int>& neighborCount,
                                  const std::vector<int>& neighborIndex,
                                  std::vector<double>& nngpB,
                                  std::vector<double>& nngpF,
                                  double *C,
                                  double *rhs,
                                  double *c)
{
  std::fill(nngpB.begin(), nngpB.end(), 0.0);
  std::fill(nngpF.begin(), nngpF.end(), 1.0);
  for(int pos = 0; pos < n; pos++){
    int fail = buildNngpConditionalOne(
      coords, n, coordDim, covariance, pos, orderIndex, neighborCount,
      neighborIndex, nngpB, nngpF, C, rhs, c
    );
    if(fail < 0)
      magcomp::error("NNGP conditional variance is non-positive or non-finite at ordered location %i", pos + 1);
    if(fail > 0)
      magcomp::error("NNGP neighbor covariance solve failed at ordered location %i", pos + 1);
  }
}

static void assembleSparseNngpPrecision(SparseNngpPrecision& sp,
                                        const std::vector<int>& neighborCount,
                                        const std::vector<double>& nngpB,
                                        const std::vector<double>& nngpF,
                                        const double *diagPrecision,
                                        double priorScale)
{
  int i, pos, q, a, b, cachePos, nz, idx;
  double *Mx;
  double invF, ba, bb;

  if(!sp.active || sp.M == NULL)
    magcomp::error("sparse NNGP precision workspace is not initialized");

  Mx = (double*) sp.M->x;
  nz = (int) sp.M->nzmax;
  for(i = 0; i < nz; i++)
    Mx[i] = 0.0;

  for(pos = 0; pos < sp.n; pos++){
    q = neighborCount[pos];
    invF = priorScale / nngpF[pos];
    cachePos = sp.nngpCacheStart[pos];

    Mx[sp.nngpCacheIndex[cachePos++]] += invF;
    for(a = 0; a < q; a++){
      ba = nngpB[pos + a * sp.n];
      Mx[sp.nngpCacheIndex[cachePos++]] -= ba * invF;
    }
    for(a = 0; a < q; a++){
      ba = nngpB[pos + a * sp.n];
      for(b = 0; b <= a; b++){
        bb = nngpB[pos + b * sp.n];
        Mx[sp.nngpCacheIndex[cachePos++]] += ba * bb * invF;
      }
    }
  }

  for(i = 0; i < sp.n; i++){
    idx = sp.diagCache[i];
    Mx[idx] += diagPrecision[i];
  }
}

static void sampleSparseGaussianPrecision(SparseNngpPrecision& sp,
                                          double *draw,
                                          const double *rhs,
                                          const char *context,
                                          int iter)
{
  cholmod_dense rhsView, zView, *sol = NULL, *noisePerm = NULL, *noise = NULL;
  struct DenseGuard {
    cholmod_common &cm;
    cholmod_dense *&sol, *&noisePerm, *&noise;
    ~DenseGuard() {
      if(noise) M_cholmod_free_dense(&noise, &cm);
      if(noisePerm) M_cholmod_free_dense(&noisePerm, &cm);
      if(sol) M_cholmod_free_dense(&sol, &cm);
    }
  } denseGuard{sp.cm, sol, noisePerm, noise};
  double *solx, *noisex;
  int i;

  if(!sp.active || sp.M == NULL || sp.MFactor == NULL)
    magcomp::error("sparse precision workspace is not initialized");

  sp.cm.status = CHOLMOD_OK;
  if(!M_cholmod_factorize(sp.M, sp.MFactor, &sp.cm) ||
     sp.cm.status != CHOLMOD_OK){
    magcomp::error("CHOLMOD factorize failed in %s update at iteration %i",
             context, iter);
  }

  for(i = 0; i < sp.n; i++){
    if(!R_FINITE(rhs[i]))
      magcomp::error("non-finite rhs in %s update at iteration %i", context, iter);
    sp.rhsWork[i] = rhs[i];
  }

  std::memset(&rhsView, 0, sizeof(cholmod_dense));
  rhsView.nrow = sp.n;
  rhsView.ncol = 1;
  rhsView.nzmax = sp.n;
  rhsView.d = sp.n;
  rhsView.xtype = CHOLMOD_REAL;
  rhsView.dtype = CHOLMOD_DOUBLE;
  rhsView.x = (void*) sp.rhsWork.data();
  rhsView.z = NULL;

  sol = M_cholmod_solve(CHOLMOD_A, sp.MFactor, &rhsView, &sp.cm);
  if(sol == NULL || sp.cm.status != CHOLMOD_OK)
    magcomp::error("CHOLMOD solve(A) failed in %s update at iteration %i",
             context, iter);

  solx = (double*) sol->x;
  for(i = 0; i < sp.n; i++)
    draw[i] = solx[i];

  for(i = 0; i < sp.n; i++)
    sp.zWork[i] = Rf_rnorm(0.0, 1.0);

  std::memset(&zView, 0, sizeof(cholmod_dense));
  zView.nrow = sp.n;
  zView.ncol = 1;
  zView.nzmax = sp.n;
  zView.d = sp.n;
  zView.xtype = CHOLMOD_REAL;
  zView.dtype = CHOLMOD_DOUBLE;
  zView.x = (void*) sp.zWork.data();
  zView.z = NULL;

  noisePerm = M_cholmod_solve(CHOLMOD_Lt, sp.MFactor, &zView, &sp.cm);
  if(noisePerm == NULL || sp.cm.status != CHOLMOD_OK){
    M_cholmod_free_dense(&sol, &sp.cm);
    magcomp::error("CHOLMOD solve(Lt) failed in %s update at iteration %i",
             context, iter);
  }

  noise = M_cholmod_solve(CHOLMOD_Pt, sp.MFactor, noisePerm, &sp.cm);
  if(noise == NULL || sp.cm.status != CHOLMOD_OK){
    M_cholmod_free_dense(&noisePerm, &sp.cm);
    M_cholmod_free_dense(&sol, &sp.cm);
    magcomp::error("CHOLMOD solve(Pt) failed in %s update at iteration %i",
             context, iter);
  }

  noisex = (double*) noise->x;
  for(i = 0; i < sp.n; i++)
    draw[i] += noisex[i];

  M_cholmod_free_dense(&noise, &sp.cm);
  M_cholmod_free_dense(&noisePerm, &sp.cm);
  M_cholmod_free_dense(&sol, &sp.cm);
}

static double nngpVectorScaledSsq(const std::vector<double>& value,
                                  int n,
                                  const std::vector<int>& orderIndex,
                                  const std::vector<int>& neighborCount,
                                  const std::vector<int>& neighborIndex,
                                  const std::vector<double>& nngpB,
                                  const std::vector<double>& nngpF)
{
  double out = 0.0;

  for(int pos = 0; pos < n; pos++){
    int i = orderIndex[pos];
    int q = neighborCount[pos];
    double mean = 0.0;

    for(int a = 0; a < q; a++){
      int nb = neighborIndex[pos + a * n];
      mean += nngpB[pos + a * n] * value[nb];
    }

    double res = value[i] - mean;
    double f = nngpF[pos];
    if(!R_FINITE(f) || f <= 0.0)
      magcomp::error("invalid NNGP conditional variance");
    out += res * res / f;
  }

  return out;
}

static double nngpVectorLogPrior(const std::vector<double>& value,
                                 int n,
                                 const std::vector<int>& orderIndex,
                                 const std::vector<int>& neighborCount,
                                 const std::vector<int>& neighborIndex,
                                 const std::vector<double>& nngpB,
                                 const std::vector<double>& nngpF,
                                 double sigmaSq)
{
  double out = 0.0;

  if(!R_FINITE(sigmaSq) || sigmaSq <= 0.0)
    magcomp::error("invalid NNGP marginal variance");

  for(int pos = 0; pos < n; pos++){
    double f = nngpF[pos];
    if(!R_FINITE(f) || f <= 0.0)
      magcomp::error("invalid NNGP conditional variance");
  }

  out = -0.5 * (double) n * std::log(sigmaSq);
  out += -0.5 * nngpVectorScaledSsq(value, n, orderIndex, neighborCount,
                                    neighborIndex, nngpB, nngpF) / sigmaSq;
  for(int pos = 0; pos < n; pos++)
    out += -0.5 * std::log(nngpF[pos]);

  return out;
}

static double nngpFactorColumnLogPrior(
  const std::vector<double>& factorScore,
  int n,
  int k,
  const std::vector<int>& orderIndex,
  const std::vector<int>& neighborCount,
  const std::vector<int>& neighborIndex,
  const std::vector<double>& nngpB,
  const std::vector<double>& nngpF)
{
  double out = 0.0;

  for(int pos = 0; pos < n; pos++){
    double f = nngpF[pos];
    if(!R_FINITE(f) || f <= 0.0)
      magcomp::error("invalid NNGP conditional variance");
    out += -0.5 * std::log(f);

    int i = orderIndex[pos];
    int q = neighborCount[pos];
    double mean = 0.0;
    for(int a = 0; a < q; a++){
      int nb = neighborIndex[pos + a * n];
      mean += nngpB[pos + a * n] * factorScore[nb + k * n];
    }
    double res = factorScore[i + k * n] - mean;
    out += -0.5 * res * res / f;
  }

  return out;
}

static std::string getCompositionFactorName(SEXP termS,
                                            const char *termType,
                                            int termIndex)
{
  SEXP nameS = getOptionalListElement(termS, "name");
  if(!Rf_isNull(nameS)){
    if(TYPEOF(nameS) != STRSXP || Rf_length(nameS) < 1)
      magcomp::error("composition factor term name must be a character value");
    const char *name = CHAR(STRING_ELT(nameS, 0));
    if(std::strlen(name) > 0)
      return std::string(name);
  }

  char tmp[64];
  std::snprintf(tmp, sizeof(tmp), "%s_%i", termType, termIndex + 1);
  return std::string(tmp);
}

static void initCompositionFactorState(CompositionFactorState& factor, int n)
{
  factor.name.clear();
  factor.termType.clear();
  factor.isIid = false;
  factor.isSpatial = false;
  factor.denseGp = false;
  factor.varying = false;
  factor.K = 0;
  factor.nGroup = 0;
  factor.groupIndex.assign(n, 0);
  factor.anchorIndex.clear();
  factor.covariate.assign(n, 1.0);
  factor.coordStorage.clear();
  factor.loading.clear();
  factor.factorScore.clear();
  factor.coords = NULL;
  factor.coordDim = 0;
  factor.nNeighbor = 0;
  factor.nngpGraph = NULL;
  factor.orderS = R_NilValue;
  factor.neighborCountS = R_NilValue;
  factor.neighborS = R_NilValue;
  factor.gpPrecision.clear();
  factor.gpPrecisionProposal.clear();
  factor.nngpB.clear();
  factor.nngpF.clear();
  factor.nngpBProposal.clear();
  factor.nngpFProposal.clear();
  factor.gpLogDet.clear();
  factor.gpLogDetProposal.clear();
  factor.covScale.clear();
  factor.estimateCovariance = false;
  factor.covAccept.clear();
  factor.covAttempt.clear();
  factor.covMetrop.clear();
}

static std::string getMagnitudeProcessName(SEXP termS,
                                           const char *termType,
                                           int termIndex)
{
  SEXP nameS = getOptionalListElement(termS, "name");
  if(!Rf_isNull(nameS)){
    if(TYPEOF(nameS) != STRSXP || Rf_length(nameS) < 1)
      magcomp::error("magnitude process term name must be a character value");
    const char *name = CHAR(STRING_ELT(nameS, 0));
    if(std::strlen(name) > 0)
      return std::string(name);
  }

  char tmp[64];
  std::snprintf(tmp, sizeof(tmp), "%s_%i", termType, termIndex + 1);
  return std::string(tmp);
}

static void initMagnitudeProcessState(MagnitudeProcessState& proc, int n)
{
  proc.name.clear();
  proc.termType.clear();
  proc.isIid = false;
  proc.isSpatial = false;
  proc.denseGp = false;
  proc.varying = false;
  proc.n = n;
  proc.nGroup = 0;
  proc.coordDim = 0;
  proc.nNeighbor = 0;
  proc.coords = NULL;
  proc.nngpGraph = NULL;
  proc.orderS = R_NilValue;
  proc.neighborCountS = R_NilValue;
  proc.neighborS = R_NilValue;
  proc.groupIndex.assign(n, 0);
  proc.covariate.assign(n, 1.0);
  proc.effect.clear();
  proc.gpPrecision.clear();
  proc.gpPrecisionProposal.clear();
  proc.nngpB.clear();
  proc.nngpF.clear();
  proc.nngpBProposal.clear();
  proc.nngpFProposal.clear();
  proc.gpLogDet = 0.0;
  proc.gpLogDetProposal = 0.0;
  proc.sigmaSq = 1.0;
  proc.covScale = 0.0;
  proc.estimateCovariance = false;
  proc.covAccept = 0;
  proc.covAttempt = 0;
  initScalarMetropolisState(proc.covMetrop);
}

static void loadMagnitudeState(MagnitudeState& mag,
                               SEXP yS,
                               SEXP areaS,
                               SEXP xS)
{
  SEXP yDimS = Rf_getAttrib(yS, R_DimSymbol);
  SEXP xDimS = Rf_getAttrib(xS, R_DimSymbol);

  if(Rf_length(yDimS) != 2)
    magcomp::error("y must be a matrix");
  if(Rf_length(xDimS) != 2)
    magcomp::error("X_magnitude must be a matrix");

  mag.n = INTEGER(yDimS)[0];
  mag.p = INTEGER(xDimS)[1];
  if(INTEGER(xDimS)[0] != mag.n)
    magcomp::error("X_magnitude must have one row per sampling unit");
  if(Rf_length(areaS) != mag.n)
    magcomp::error("area must have one value per sampling unit");

  mag.X = REAL(xS);
  mag.total.assign(mag.n, 0);
  mag.area.assign(mag.n, 0.0);
  mag.coef.assign(std::max(1, mag.p), 0.0);
  mag.omega.assign(mag.n, 0.0);
  mag.nIid = 0;
  mag.nSpatial = 0;
  mag.process.clear();

  for(int i = 0; i < mag.n; i++){
    mag.area[i] = REAL(areaS)[i];
    if(!R_FINITE(mag.area[i]) || mag.area[i] <= 0.0)
      magcomp::error("area must contain positive finite values");
  }
}

static void loadCovariance(CovarianceMetropolis& block, SEXP term, int k, int K)
{
  block.model = &magcomp::covarianceModel(getListString(term, "cov_model"));
  int P = block.model->size;
  SEXP values = getListElement(term, "covariance_values");
  SEXP lower = getListElement(term, "covariance_lower");
  SEXP upper = getListElement(term, "covariance_upper");
  SEXP tuning = getListElement(term, "covariance_tuning");
  SEXP fixed = getListElement(term, "covariance_fixed");
  for(SEXP x : {values, lower, upper, tuning, fixed})
    if(Rf_length(x) != K*P) magcomp::error("covariance setting dimensions disagree with registry");
  block.theta.resize(P); block.lower.resize(P); block.upper.resize(P); block.free.clear();
  for(int p = 0; p < P; ++p){
    int i = k + K*p;
    block.theta[p] = REAL(values)[i];
    block.lower[p] = REAL(lower)[i]; block.upper[p] = REAL(upper)[i];
    if(!LOGICAL(fixed)[i]){
      if(!(block.theta[p] > block.lower[p] && block.theta[p] < block.upper[p]))
        magcomp::error("covariance starting value outside prior support");
      block.free.push_back(p);
    }
  }
  magcomp::validateCovariance(*block.model, block.theta.data());
  int q = block.free.size();
  block.shape.assign(q*q, 0); block.chol.assign(q*q, 0);
  for(int j = 0; j < q; ++j){
    double sd = REAL(tuning)[k + K*block.free[j]];
    if(!(R_FINITE(sd) && sd > 0)) magcomp::error("invalid covariance proposal tuning");
    block.shape[j+j*q] = sd*sd;
    block.chol[j+j*q] = sd;
  }
}

static void loadMagnitudeProcessTerm(MagnitudeState& mag, SEXP processS)
{
  if(Rf_isNull(processS) || Rf_length(processS) == 0)
    return;
  if(TYPEOF(processS) != VECSXP)
    magcomp::error("magnitude_process must be a list");

  for(int termIndex = 0; termIndex < Rf_length(processS); termIndex++){
    SEXP termS = VECTOR_ELT(processS, termIndex);
    const char *termType = getListString(termS, "term_type");
    MagnitudeProcessState proc;
    initMagnitudeProcessState(proc, mag.n);
    proc.name = getMagnitudeProcessName(termS, termType, termIndex);
    proc.termType = std::string(termType);
    SEXP variancePriorS = getListElement(termS, "variance_prior");
    if(Rf_length(variancePriorS) != 2)
      magcomp::error("magnitude variance prior must have shape and scale");
    proc.varianceShape = REAL(variancePriorS)[0];
    proc.varianceScale = REAL(variancePriorS)[1];
    if(!R_FINITE(proc.varianceShape) || proc.varianceShape <= 0.0 ||
       !R_FINITE(proc.varianceScale) || proc.varianceScale <= 0.0)
      magcomp::error("magnitude variance prior must be positive and finite");

    if(std::strcmp(termType, "iid") != 0 &&
       std::strcmp(termType, "gp") != 0 &&
       std::strcmp(termType, "nngp") != 0){
      magcomp::error("sampler currently supports iid(), gp(), and nngp() magnitude terms only");
    }

    SEXP covariateS = getOptionalListElement(termS, "covariate");

    if(std::strcmp(termType, "iid") == 0){
      proc.isIid = true;
      proc.nGroup = getListInteger(termS, "n_group");
      if(proc.nGroup < 1)
        magcomp::error("iid() must have at least one group");
      proc.effect.assign(proc.nGroup, 0.0);

      SEXP groupS = getListElement(termS, "group");
      if(Rf_length(groupS) != mag.n)
        magcomp::error("iid() group index must have one value per sampling unit");
      for(int i = 0; i < mag.n; i++){
        proc.groupIndex[i] = INTEGER(groupS)[i] - 1;
        if(proc.groupIndex[i] < 0 || proc.groupIndex[i] >= proc.nGroup)
          magcomp::error("iid() group index is out of range");
      }

      proc.varying = !Rf_isNull(covariateS);
      if(proc.varying){
        if(Rf_length(covariateS) != mag.n)
          magcomp::error("iid() covariate must have one value per sampling unit");
        for(int i = 0; i < mag.n; i++){
          proc.covariate[i] = REAL(covariateS)[i];
          if(!R_FINITE(proc.covariate[i]))
            magcomp::error("iid() covariate must contain finite values");
        }
      }
      mag.nIid++;
      mag.process.push_back(proc);
      continue;
    }

    proc.isSpatial = true;
    proc.denseGp = std::strcmp(termType, "gp") == 0;
    proc.nGroup = mag.n;
    proc.effect.assign(mag.n, 0.0);
    if(proc.denseGp){
      checkArraySize(mag.n, mag.n);
      proc.gpPrecision.assign(std::max(1, mag.n * mag.n), 0.0);
      proc.gpPrecisionProposal.assign(std::max(1, mag.n * mag.n), 0.0);
    }

    SEXP coordsS = getListElement(termS, "coords");
    if(!Rf_isMatrix(coordsS))
      magcomp::error("%s() coordinates must be a matrix", termType);
    SEXP coordsDimS = Rf_getAttrib(coordsS, R_DimSymbol);
    if(INTEGER(coordsDimS)[0] != mag.n)
      magcomp::error("%s() coordinates must have one row per sampling unit", termType);
    proc.coordDim = INTEGER(coordsDimS)[1];
    if(proc.coordDim < 1)
      magcomp::error("%s() coordinate dimension is invalid", termType);
    proc.coords = REAL(coordsS);

    loadCovariance(proc.covMetrop, termS, 0, 1);
    proc.estimateCovariance = !proc.covMetrop.free.empty();
    proc.covScale = 1.0;

    if(!proc.denseGp){
      proc.nNeighbor = getListInteger(termS, "m");
      if(proc.nNeighbor < 1)
        magcomp::error("nngp() m must be positive");
      proc.orderS = getListElement(termS, "order");
      proc.neighborCountS = getListElement(termS, "neighbor_count");
      proc.neighborS = getListElement(termS, "neighbor_index");
      SEXP neighborDimS = Rf_getAttrib(proc.neighborS, R_DimSymbol);
      if(Rf_length(neighborDimS) != 2)
        magcomp::error("nngp() neighbor_index must be a matrix");
      if(INTEGER(neighborDimS)[0] != mag.n)
        magcomp::error("nngp() neighbor_index must have one row per ordered location");
      proc.nNeighbor = INTEGER(neighborDimS)[1];
      if(proc.nNeighbor < 0)
        magcomp::error("nngp() neighbor dimension is invalid");
      proc.nngpB.assign(std::max(1, mag.n * std::max(1, proc.nNeighbor)), 0.0);
      proc.nngpF.assign(std::max(1, mag.n), 1.0);
      proc.nngpBProposal.assign(std::max(1, mag.n *
                                         std::max(1, proc.nNeighbor)), 0.0);
      proc.nngpFProposal.assign(std::max(1, mag.n), 1.0);
    }

    proc.varying = !Rf_isNull(covariateS);
    if(proc.varying){
      if(Rf_length(covariateS) != mag.n)
        magcomp::error("%s() covariate must have one value per sampling unit",
                 termType);
      for(int i = 0; i < mag.n; i++){
        proc.covariate[i] = REAL(covariateS)[i];
        if(!R_FINITE(proc.covariate[i]))
          magcomp::error("%s() covariate must contain finite values", termType);
      }
    }

    mag.nSpatial++;
    mag.process.push_back(proc);
  }
}

static void loadAllocationBlock(AllocationBlockState& block,
                                SEXP blockS,
                                SEXP yS,
                                SEXP xS,
                                SEXP smoothS,
                                SEXP priorsS)
{
  SEXP yDimS = Rf_getAttrib(yS, R_DimSymbol);
  SEXP xDimS = Rf_getAttrib(xS, R_DimSymbol);

  if(Rf_length(yDimS) != 2)
    magcomp::error("y must be a matrix");
  if(Rf_length(xDimS) != 2)
    magcomp::error("X_composition must be a matrix");

  block.n = INTEGER(yDimS)[0];
  block.J = INTEGER(yDimS)[1];
  block.Jstick = block.J - 1;
  block.coefStick = block.Jstick;
  block.p = INTEGER(xDimS)[1];
  if(block.J < 2)
    magcomp::error("fixed-area stem counts require at least two categories");
  if(INTEGER(xDimS)[0] != block.n)
    magcomp::error("X_composition must have one row per sampling unit");

  block.parentSpecificCoef = false;
  if(!Rf_isNull(blockS)){
    SEXP coefStickS = getOptionalListElement(blockS, "coef_stick_count");
    if(!Rf_isNull(coefStickS))
      block.coefStick = magcomp::guarded_Rf_asInteger(coefStickS);
    SEXP parentSpecificS =
      getOptionalListElement(blockS, "parent_specific_coef");
    if(!Rf_isNull(parentSpecificS))
      block.parentSpecificCoef = magcomp::guarded_Rf_asLogical(parentSpecificS) == TRUE;
  }
  if(block.coefStick < block.Jstick)
    magcomp::error("allocation block coef_stick_count is too small");
  if(block.parentSpecificCoef && block.coefStick % block.Jstick != 0)
    magcomp::error("parent-specific coef_stick_count must be a multiple of J - 1");

  block.X = REAL(xS);
  block.count.assign(block.n * block.J, 0);
  block.trials.assign(block.n * block.Jstick, 0);
  block.unitIndex.assign(block.n, 0);
  block.parentIndex.assign(block.n, -1);
	  block.coef.assign(block.coefStick * block.p, 0.0);
	  block.omega.assign(block.n * block.Jstick, 0.0);
	  block.workingResponse.assign(block.n * block.Jstick, 0.0);
	  block.proportion.assign(block.n * block.J, 0.0);
  block.aggregateCount.clear();
  block.aggregateTrials.clear();
  block.hasFactor = false;
  block.nFactor = 0;
  block.nSpatialFactor = 0;
  block.aggregateInterceptOnly = false;
  block.factor.clear();

  SEXP blockSmoothS = R_NilValue;
  if(!Rf_isNull(blockS))
    blockSmoothS = getOptionalListElement(blockS, "smooth");
  if(Rf_isNull(blockSmoothS))
    blockSmoothS = smoothS;

  block.interceptSmoothOrder =
    smoothOrder(getListString(blockSmoothS, "intercept"), "intercept");
  block.slopeSmoothOrder =
    smoothOrder(getListString(blockSmoothS, "slope"), "slope");
  block.interceptPriorPrec = 1.0 / getListDouble(priorsS, "intercept_variance");
  block.betaPriorPrec = 1.0 / getListDouble(priorsS, "slope_variance");
  block.interceptSmoothPrec = 1.0 / getListDouble(priorsS,
                                              "intercept_smoothing_variance");
  block.slopeSmoothPrec = 1.0 / getListDouble(priorsS,
                                             "slope_smoothing_variance");
}

static void loadBlockIndex(AllocationBlockState& block,
                           SEXP blockS,
                           int nUnit)
{
  SEXP unitIndexS = getOptionalListElement(blockS, "unit_index");
  SEXP parentIndexS = getOptionalListElement(blockS, "parent_index");

  if(Rf_isNull(unitIndexS)){
    for(int i = 0; i < block.n; i++)
      block.unitIndex[i] = i;
  }else{
    if(Rf_length(unitIndexS) != block.n)
      magcomp::error("allocation block unit_index has wrong length");
    for(int i = 0; i < block.n; i++){
      block.unitIndex[i] = INTEGER(unitIndexS)[i] - 1;
      if(block.unitIndex[i] < 0 || block.unitIndex[i] >= nUnit)
        magcomp::error("allocation block unit_index is out of range");
    }
  }

  if(!Rf_isNull(parentIndexS)){
    if(Rf_length(parentIndexS) != block.n)
      magcomp::error("allocation block parent_index has wrong length");
    for(int i = 0; i < block.n; i++){
      block.parentIndex[i] = INTEGER(parentIndexS)[i] - 1;
      if(block.parentIndex[i] < 0)
        magcomp::error("allocation block parent_index must be positive");
    }
  }
}

static void loadAllocationFactorTerms(AllocationBlockState& block,
                                      SEXP compositionBlockS,
                                      int nUnit)
{
  if(Rf_isNull(compositionBlockS))
    return;

  SEXP processS = getOptionalListElement(compositionBlockS, "process_terms");
  if(Rf_isNull(processS) || Rf_length(processS) == 0)
    return;
  if(TYPEOF(processS) != VECSXP)
    magcomp::error("composition process_terms must be a list");
  block.factor.reserve(block.factor.size() + Rf_length(processS));

  for(int termIdx = 0; termIdx < Rf_length(processS); termIdx++){
    SEXP termS = VECTOR_ELT(processS, termIdx);
    const char *termType = getListString(termS, "term_type");
    bool isIid = std::strcmp(termType, "factor_iid") == 0;
    bool isGp = std::strcmp(termType, "factor_gp") == 0;
    bool isNngp = std::strcmp(termType, "factor_nngp") == 0;
    if(!isIid && !isGp && !isNngp){
      magcomp::error("sampler currently supports factor_iid(), "
               "factor_gp(), and factor_nngp() composition terms only");
    }

    CompositionFactorState factor;
    initCompositionFactorState(factor, block.n);
    factor.loadingPriorPrec = 1.0 / getListDouble(termS, "loading_variance");
    factor.name = getCompositionFactorName(termS, termType, termIdx);
    factor.termType = std::string(termType);
    factor.isIid = isIid;
    factor.isSpatial = !isIid;
    factor.denseGp = isGp;
    factor.K = getListInteger(termS, "K");
    if(factor.K < 1 || factor.K > block.coefStick){
      magcomp::error("%s() K must be between one and the number of coefficient sticks",
               termType);
    }

    if(isIid){
      factor.nGroup = getListInteger(termS, "n_group");
      if(factor.nGroup < 1)
        magcomp::error("factor_iid() must have at least one group");

      SEXP groupS = getListElement(termS, "group");
      if(Rf_length(groupS) != block.n){
        magcomp::error("factor_iid() group index must have one value per "
                 "allocation row");
      }
      for(int i = 0; i < block.n; i++){
        factor.groupIndex[i] = INTEGER(groupS)[i] - 1;
        if(factor.groupIndex[i] < 0 ||
           factor.groupIndex[i] >= factor.nGroup){
          magcomp::error("factor_iid() group index is out of range");
        }
      }
    }else{
      SEXP coordsS = getListElement(termS, "coords");
      if(!Rf_isMatrix(coordsS))
        magcomp::error("%s() coordinates must be a matrix", termType);
      SEXP coordsDimS = Rf_getAttrib(coordsS, R_DimSymbol);
      if(INTEGER(coordsDimS)[0] != block.n)
        magcomp::error("%s() coordinates must have one row per allocation row",
                 termType);
      factor.coordDim = INTEGER(coordsDimS)[1];
      if(factor.coordDim < 1)
        magcomp::error("%s() coordinate dimension is invalid", termType);
      if(nUnit < 1)
        magcomp::error("%s() requires at least one sampling unit", termType);

      factor.nGroup = nUnit;
      factor.coordStorage.assign(factor.nGroup * factor.coordDim, NA_REAL);
      std::vector<int> seen(factor.nGroup, 0);
      for(int i = 0; i < block.n; i++){
        int g = block.unitIndex[i];
        if(g < 0 || g >= factor.nGroup)
          magcomp::error("%s() unit index is out of range", termType);
        factor.groupIndex[i] = g;
        for(int d = 0; d < factor.coordDim; d++){
          double value = REAL(coordsS)[i + d * block.n];
          if(!R_FINITE(value))
            magcomp::error("%s() coordinates must contain finite values", termType);
          int outIdx = g + d * factor.nGroup;
          if(!seen[g]){
            factor.coordStorage[outIdx] = value;
          }else if(std::fabs(factor.coordStorage[outIdx] - value) >
                   1.0e-10){
            magcomp::error("%s() coordinates for repeated sampling units must match",
                     termType);
          }
        }
        seen[g] = 1;
      }
      for(int g = 0; g < factor.nGroup; g++){
        if(!seen[g])
          magcomp::error("%s() spatial factor is missing a sampling unit", termType);
      }
      factor.coords = factor.coordStorage.data();
      if(isNngp){
        factor.nNeighbor = getListInteger(termS, "m");
        if(factor.nNeighbor < 1)
          magcomp::error("factor_nngp() m must be positive");
        factor.orderS = getListElement(termS, "order");
        factor.neighborCountS = getListElement(termS, "neighbor_count");
        factor.neighborS = getListElement(termS, "neighbor_index");
        SEXP neighborDimS = Rf_getAttrib(factor.neighborS, R_DimSymbol);
        if(Rf_length(neighborDimS) != 2)
          magcomp::error("factor_nngp() neighbor_index must be a matrix");
        if(INTEGER(neighborDimS)[0] != factor.nGroup){
          magcomp::error("factor_nngp() neighbor_index must have one row per "
                   "spatial factor location");
        }
        factor.nNeighbor = INTEGER(neighborDimS)[1];
        if(factor.nNeighbor < 0)
          magcomp::error("factor_nngp() neighbor dimension is invalid");
        factor.nngpB.assign(
          factor.K,
          std::vector<double>(
            std::max(1, factor.nGroup * std::max(1, factor.nNeighbor)), 0.0)
        );
        factor.nngpF.assign(factor.K,
                            std::vector<double>(std::max(1, factor.nGroup),
                                                1.0));
        factor.nngpBProposal.assign(
          factor.K,
          std::vector<double>(
            std::max(1, factor.nGroup * std::max(1, factor.nNeighbor)), 0.0)
        );
        factor.nngpFProposal.assign(
          factor.K,
          std::vector<double>(std::max(1, factor.nGroup), 1.0)
        );
      }

      factor.covMetrop.resize(factor.K);
      factor.covAccept.assign(factor.K, 0);
      factor.covAttempt.assign(factor.K, 0);
      factor.covScale.assign(factor.K, 1.0);
      factor.estimateCovariance = false;
      for(int k = 0; k < factor.K; ++k){
        loadCovariance(factor.covMetrop[k], termS, k, factor.K);
        factor.estimateCovariance = factor.estimateCovariance || !factor.covMetrop[k].free.empty();
      }
      if(factor.denseGp){
        checkArraySize(factor.nGroup, factor.nGroup);
        factor.gpPrecision.assign(
          factor.K,
          std::vector<double>(std::max(1, factor.nGroup * factor.nGroup),
                              0.0)
        );
        factor.gpPrecisionProposal.assign(
          factor.K,
          std::vector<double>(std::max(1, factor.nGroup * factor.nGroup),
                              0.0)
        );
        factor.gpLogDet.assign(factor.K, 0.0);
        factor.gpLogDetProposal.assign(factor.K, 0.0);
      }
    }

    SEXP covariateS = getOptionalListElement(termS, "covariate");
    factor.varying = !Rf_isNull(covariateS);
    if(factor.varying){
      if(Rf_length(covariateS) != block.n)
        magcomp::error("%s() covariate must have one value per allocation row",
                 termType);
      for(int i = 0; i < block.n; i++){
        factor.covariate[i] = REAL(covariateS)[i];
        if(!R_FINITE(factor.covariate[i]))
          magcomp::error("%s() covariate must contain finite values", termType);
      }
    }

    SEXP anchorS = getListElement(termS, "anchor_index");
    if(Rf_length(anchorS) != factor.K)
      magcomp::error("%s() anchor_index must have length K", termType);
    factor.anchorIndex.assign(factor.K, 0);
    for(int k = 0; k < factor.K; k++){
      factor.anchorIndex[k] = INTEGER(anchorS)[k] - 1;
      if(factor.anchorIndex[k] < 0 ||
         factor.anchorIndex[k] >= block.coefStick)
        magcomp::error("%s() anchor_index is out of range", termType);
      for(int h = 0; h < k; h++){
        if(factor.anchorIndex[h] == factor.anchorIndex[k])
          magcomp::error("%s() anchor_index values must be unique", termType);
      }
    }

    factor.loading.assign(block.coefStick * factor.K, 0.0);
    factor.factorScore.assign(factor.nGroup * factor.K, 0.0);
    for(int k = 0; k < factor.K; k++)
      factor.loading[factor.anchorIndex[k] + k * block.coefStick] = 1.0;

    if(factor.isSpatial)
      block.nSpatialFactor++;
    block.factor.push_back(factor);
    if(block.factor.back().isSpatial)
      block.factor.back().coords = block.factor.back().coordStorage.data();
  }

  block.nFactor = (int) block.factor.size();
  block.hasFactor = block.nFactor > 0;
}

static void loadMagnitudeTotals(MagnitudeState& mag, SEXP yS)
{
  int totalY = 0;
  SEXP yDimS = Rf_getAttrib(yS, R_DimSymbol);
  int n = INTEGER(yDimS)[0];
  int J = INTEGER(yDimS)[1];

  if(n != mag.n)
    magcomp::error("y must have one row per sampling unit");

  for(int j = 0; j < J; j++){
    for(int i = 0; i < n; i++){
      int y = getMatrixCount(yS, i + j * n);
      mag.total[i] += y;
      totalY += y;
    }
  }

  if(totalY <= 0)
    magcomp::error("y must contain at least one observed count");
}

static void loadAllocationCounts(AllocationBlockState& block, SEXP yS)
{
  int totalY = 0;

  for(int j = 0; j < block.J; j++){
    for(int i = 0; i < block.n; i++){
      int idx = i + j * block.n;
      int y = getMatrixCount(yS, idx);
      block.count[idx] = y;
      totalY += y;
    }
  }

  if(totalY <= 0)
    magcomp::error("y must contain at least one observed count");

  for(int i = 0; i < block.n; i++){
    int remaining = 0;
    for(int j = 0; j < block.J; j++)
      remaining += block.count[i + j * block.n];
    for(int j = 0; j < block.Jstick; j++){
      block.trials[i + j * block.n] = remaining;
      remaining -= block.count[i + j * block.n];
    }
  }
}

static void configureAllocationAggregation(AllocationBlockState& block)
{
  bool interceptOnly = block.p == 1 && !block.hasFactor;

  if(interceptOnly){
    for(int i = 0; i < block.n; i++){
      if(!R_FINITE(block.X[i]) || std::fabs(block.X[i] - 1.0) > 1e-12){
        interceptOnly = false;
        break;
      }
    }
  }

  block.aggregateInterceptOnly = interceptOnly;
  block.aggregateCount.clear();
  block.aggregateTrials.clear();
  if(!interceptOnly)
    return;

  block.aggregateCount.assign(block.coefStick, 0);
  block.aggregateTrials.assign(block.coefStick, 0);

  for(int j = 0; j < block.Jstick; j++){
    for(int i = 0; i < block.n; i++){
      int idx = i + j * block.n;
      int coefIdx = blockCoefIndex(block, i, j);
      block.aggregateCount[coefIdx] += block.count[idx];
      block.aggregateTrials[coefIdx] += block.trials[idx];
    }
  }
}

static void initializeMagnitude(MagnitudeState& mag, SEXP startingS)
{
  double totalY = 0.0;
  double totalArea = 0.0;

  for(int i = 0; i < mag.n; i++){
    totalY += (double) mag.total[i];
    totalArea += mag.area[i];
  }

  if(mag.p > 0)
    mag.coef[0] = std::log((totalY + 0.5) / totalArea);

  SEXP magnitudeCoefficientStartS = getOptionalListElement(startingS, "magnitude_coefficients");
  if(!Rf_isNull(magnitudeCoefficientStartS)){
    if(Rf_length(magnitudeCoefficientStartS) == 1){
      mag.coef[0] = magcomp::guarded_Rf_asReal(magnitudeCoefficientStartS);
    }else if(Rf_length(magnitudeCoefficientStartS) == mag.p){
      SEXP tmpS = PROTECT(magcomp::guarded_Rf_coerceVector(magnitudeCoefficientStartS, REALSXP));
      for(int l = 0; l < mag.p; l++)
        mag.coef[l] = REAL(tmpS)[l];
      UNPROTECT(1);
    }else{
      magcomp::error("starting$magnitude_coefficients must have length one or ncol(X_magnitude)");
    }
  }

  for(int l = 0; l < mag.p; l++){
    if(!R_FINITE(mag.coef[l]))
      magcomp::error("starting$magnitude_coefficients must contain finite values");
  }

  for(int termIndex = 0; termIndex < (int) mag.process.size(); termIndex++){
    MagnitudeProcessState& proc = mag.process[termIndex];
    SEXP effectStartS =
      getOptionalBlockStart(startingS, "magnitude_process",
                            proc.name.c_str());

    if(!Rf_isNull(effectStartS)){
      if(Rf_length(effectStartS) != (int) proc.effect.size()){
        magcomp::error("starting$magnitude_process$%s has incompatible length",
                 proc.name.c_str());
      }
      SEXP tmpS = PROTECT(magcomp::guarded_Rf_coerceVector(effectStartS, REALSXP));
      for(int i = 0; i < (int) proc.effect.size(); i++){
        proc.effect[i] = REAL(tmpS)[i];
        if(!R_FINITE(proc.effect[i])){
          magcomp::error("starting$magnitude_process$%s must contain finite values",
                   proc.name.c_str());
        }
      }
      UNPROTECT(1);
    }

    SEXP sigmaStartS =
      getOptionalBlockStart(startingS, "magnitude_process_variance",
                            proc.name.c_str());
    if(!Rf_isNull(sigmaStartS)){
      if(Rf_length(sigmaStartS) != 1)
        magcomp::error("starting magnitude process sigma_sq must have length one");
      proc.sigmaSq = magcomp::guarded_Rf_asReal(sigmaStartS);
      if(!R_FINITE(proc.sigmaSq) || proc.sigmaSq <= 0.0){
        magcomp::error("starting magnitude process sigma_sq must be positive and finite");
      }
    }


  }
}

static void initializeAllocationBlock(AllocationBlockState& block,
                                      SEXP startingS,
                                      const char *blockName,
                                      CompositionLink compositionLink)
{
  double successTotal, totalRemaining, pi;
  int coefIdx;

  for(int c = 0; c < block.coefStick; c++){
    int parent = block.parentSpecificCoef ? c / block.Jstick : -1;
    int j = block.parentSpecificCoef ? c % block.Jstick : c;
    successTotal = 0.0;
    totalRemaining = 0.0;
    for(int i = 0; i < block.n; i++){
      if(block.parentSpecificCoef && block.parentIndex[i] != parent)
        continue;
      successTotal += (double) block.count[i + j * block.n];
      totalRemaining += (double) block.trials[i + j * block.n];
    }
    pi = clipProbability((successTotal + 0.5) / (totalRemaining + 1.0));
    block.coef[c] = compositionLinkQuantile(pi, compositionLink);
  }

  SEXP compStartS =
    getOptionalBlockStart(startingS, "composition_block_coefficient", blockName);
  if(!Rf_isNull(compStartS)){
    if(Rf_length(compStartS) != block.coefStick * block.p){
      magcomp::error("starting composition coefficients must have length "
               "coef_stick_count * ncol(X_composition)");
    }
    SEXP tmpS = PROTECT(magcomp::guarded_Rf_coerceVector(compStartS, REALSXP));
    for(int l = 0; l < block.p; l++){
      for(int c = 0; c < block.coefStick; c++)
        block.coef[c + l * block.coefStick] =
          REAL(tmpS)[c + l * block.coefStick];
    }
    UNPROTECT(1);
  }

  for(int l = 0; l < block.p; l++){
    for(int c = 0; c < block.coefStick; c++){
      coefIdx = c + l * block.coefStick;
      if(!R_FINITE(block.coef[coefIdx])){
        magcomp::error("starting composition coefficients must contain finite "
                 "values");
      }
    }
  }
}

static void initializeAllocationFactor(AllocationBlockState& block,
                                       SEXP startingS,
                                       const char *blockName)
{
  if(!block.hasFactor)
    return;

  SEXP blockLoadingStartS =
    getOptionalBlockStart(startingS, "composition_block_loading", blockName);
  SEXP blockFactorStartS =
    getOptionalBlockStart(startingS, "composition_block_score", blockName);


  for(int f = 0; f < (int) block.factor.size(); f++){
    CompositionFactorState& factor = block.factor[f];
    SEXP loadingStartS = R_NilValue;
    SEXP factorStartS = R_NilValue;


    if(!Rf_isNull(blockLoadingStartS)){
      loadingStartS = getOptionalListElement(blockLoadingStartS,
                                             factor.name.c_str());
      if(Rf_isNull(loadingStartS)){
        magcomp::error("starting composition_block_loading for block '%s' must be "
                 "a named list with one entry per factor term",
                 blockName == NULL ? "" : blockName);
      }
    }
    if(!Rf_isNull(loadingStartS)){
      if(Rf_length(loadingStartS) != block.coefStick * factor.K){
        magcomp::error("starting loading values for %s must have length "
                 "coef_stick_count * K", factor.name.c_str());
      }
      SEXP tmpS = PROTECT(magcomp::guarded_Rf_coerceVector(loadingStartS, REALSXP));
      for(int k = 0; k < factor.K; k++){
        for(int j = 0; j < block.coefStick; j++){
          double constrainedValue = 0.0;
          if(!loadingIsConstrained(j, k, factor.anchorIndex,
                                   constrainedValue)){
            factor.loading[j + k * block.coefStick] =
              REAL(tmpS)[j + k * block.coefStick];
            if(!R_FINITE(factor.loading[j + k * block.coefStick])){
              magcomp::error("starting loading values for %s must contain finite "
                       "values", factor.name.c_str());
            }
          }else{
            factor.loading[j + k * block.coefStick] = constrainedValue;
          }
        }
      }
      UNPROTECT(1);
    }

    if(!Rf_isNull(blockFactorStartS)){
      factorStartS = getOptionalListElement(blockFactorStartS,
                                            factor.name.c_str());
      if(Rf_isNull(factorStartS)){
        magcomp::error("starting composition_block_score for block '%s' must be "
                 "a named list with one entry per factor term",
                 blockName == NULL ? "" : blockName);
      }
    }
    if(!Rf_isNull(factorStartS)){
      if(Rf_length(factorStartS) != factor.nGroup * factor.K)
        magcomp::error("starting factor values for %s must have length n_group * K",
                 factor.name.c_str());
      SEXP tmpS = PROTECT(magcomp::guarded_Rf_coerceVector(factorStartS, REALSXP));
      for(int k = 0; k < factor.K; k++){
        for(int g = 0; g < factor.nGroup; g++){
          factor.factorScore[g + k * factor.nGroup] =
            REAL(tmpS)[g + k * factor.nGroup];
          if(!R_FINITE(factor.factorScore[g + k * factor.nGroup])){
            magcomp::error("starting factor values for %s must contain finite values",
                     factor.name.c_str());
          }
        }
      }
      UNPROTECT(1);
    }

    enforceLoadingConstraints(factor.loading, block.coefStick, factor.K,
                              factor.anchorIndex);
  }
}

static double magnitudeProcessContribution(const MagnitudeProcessState& proc,
                                           int i)
{
  int idx = proc.isIid ? proc.groupIndex[i] : i;
  if(idx < 0 || idx >= (int) proc.effect.size())
    magcomp::error("internal error: magnitude process effect index out of range");

  return proc.covariate[i] * proc.effect[idx];
}

static double magnitudeProcessEta(const MagnitudeState& mag, int i)
{
  double eta = 0.0;
  for(int termIndex = 0; termIndex < (int) mag.process.size(); termIndex++)
    eta += magnitudeProcessContribution(mag.process[termIndex], i);

  return eta;
}

static double magnitudeProcessEtaExcept(const MagnitudeState& mag,
                                        int i,
                                        int skipTerm)
{
  double eta = 0.0;
  for(int termIndex = 0; termIndex < (int) mag.process.size(); termIndex++){
    if(termIndex == skipTerm)
      continue;
    eta += magnitudeProcessContribution(mag.process[termIndex], i);
  }

  return eta;
}

static double logSumExp2(double a, double b)
{
  if(a > b)
    return a + std::log1p(std::exp(b - a));

  return b + std::log1p(std::exp(a - b));
}

static double magnitudeNbLogLikelihood(const MagnitudeState& mag,
                                       double nbSize)
{
  if(!R_FINITE(nbSize) || nbSize <= 0.0)
    return R_NegInf;

  double out = 0.0;
  double logNbSize = std::log(nbSize);
  for(int i = 0; i < mag.n; i++){
    double theta = rowDot(mag.X, mag.coef.data(), mag.n, mag.p, i);
    theta += magnitudeProcessEta(mag, i);
    double logMu = std::log(mag.area[i]) + theta;
    if(!R_FINITE(logMu))
      return R_NegInf;

    double y = (double) mag.total[i];
    double logDenom = logSumExp2(logNbSize, logMu);
    double value =
      Rf_lgammafn(y + nbSize) -
      Rf_lgammafn(nbSize) -
      Rf_lgammafn(y + 1.0) +
      nbSize * (logNbSize - logDenom) +
      y * (logMu - logDenom);
    if(!R_FINITE(value))
      return R_NegInf;
    out += value;
  }

  return out;
}

static double magnitudeNbSizeLogTarget(const MagnitudeState& mag,
                                       const MagnitudeNbSizeState& nb,
                                       double logNbSize)
{
  if(!R_FINITE(logNbSize))
    return R_NegInf;

  double nbSize = std::exp(logNbSize);
  if(!R_FINITE(nbSize) || nbSize <= 0.0)
    return R_NegInf;

  double logLik = magnitudeNbLogLikelihood(mag, nbSize);
  if(!R_FINITE(logLik))
    return R_NegInf;

  double logPrior = (nb.shape - 1.0) * logNbSize - nb.rate * nbSize;
  return logLik + logPrior + logNbSize;
}

static void updateMagnitudeNbSize(MagnitudeState& mag,
                                  MagnitudeNbSizeState& nb,
                                  int iter,
                                  int nBurn)
{
  if(!nb.estimate)
    return;

  nb.attempt++;
  int accepted = 0;
  double logCurrent = std::log(nb.value);
  double logProposal = logCurrent + nb.tuning * Rf_rnorm(0.0, 1.0);
  double logPostCurrent = magnitudeNbSizeLogTarget(mag, nb, logCurrent);
  double logPostProposal = magnitudeNbSizeLogTarget(mag, nb, logProposal);
  double logAccept = logPostProposal - logPostCurrent;

  if(R_FINITE(logAccept) &&
     std::log(Rf_runif(0.0, 1.0)) < logAccept){
    nb.value = std::exp(logProposal);
    nb.accept++;
    accepted = 1;
  }

  if(iter <= nBurn)
    updateScalarMetropolisTuning(nb.metrop, accepted, &nb.tuning);
}

static void updateMagnitudeCoefficients(MagnitudeState& mag,
                                        double nbSize,
                                        double betaPriorPrec,
                                        double *Q,
                                        double *rhs,
                                        double *draw,
                                        double *chol,
		                                  double *z,
		                                  BayesLogit_rpg_hybrid_t rpg,
		                                  BayesLogit_rpg_gamma_t rpgGamma,
		                                  int iter)
{
  double theta, eta, omega, kappa, adj, xValue;
  double logNbSize = std::log(nbSize);
  double processEta;

  for(int l = 0; l < mag.p; l++){
    rhs[l] = 0.0;
    for(int h = 0; h < mag.p; h++)
      Q[l + h * mag.p] = 0.0;
    Q[l + l * mag.p] = betaPriorPrec;
  }

  for(int i = 0; i < mag.n; i++){
    theta = rowDot(mag.X, mag.coef.data(), mag.n, mag.p, i);
    processEta = magnitudeProcessEta(mag, i);
    theta += processEta;
    eta = std::log(mag.area[i]) + theta - logNbSize;
    omega = drawPolyaGamma(rpg, rpgGamma, (double) mag.total[i] + nbSize,
                           eta, "magnitude update", iter);
    mag.omega[i] = omega;

    kappa = 0.5 * ((double) mag.total[i] - nbSize);
    adj = kappa + omega * (logNbSize - std::log(mag.area[i]) -
                            processEta);

    for(int l = 0; l < mag.p; l++){
      xValue = mag.X[i + l * mag.n];
      rhs[l] += xValue * adj;
      for(int h = 0; h < mag.p; h++)
        Q[l + h * mag.p] += xValue * omega * mag.X[i + h * mag.n];
    }
  }

  sampleGaussianPrecision(draw, Q, rhs, chol, z, mag.p,
                          "magnitude coefficient", iter);
  for(int l = 0; l < mag.p; l++)
    mag.coef[l] = draw[l];
}

static void updateMagnitudeSpatialEffect(MagnitudeState& mag,
                                         int termIndex,
                                         double nbSize,
                                         double spatialShape,
                                         double spatialScale,
                                         double *Q,
                                         double *rhs,
                                         double *draw,
                                         double *chol,
                                         double *z,
                                         double *latentDiag,
                                         double *nngpWorkC,
                                         double *nngpWorkRhs,
                                         double *nngpWorkCvec,
                                         int iter,
                                         int nBurn)
{
  double logNbSize = std::log(nbSize);
  MagnitudeProcessState& proc = mag.process[termIndex];

  if(!proc.isSpatial)
    return;

  for(int i = 0; i < mag.n; i++){
    rhs[i] = 0.0;
    latentDiag[i] = 0.0;
    if(proc.denseGp){
      for(int h = 0; h < mag.n; h++)
        Q[i + h * mag.n] = 0.0;
    }
  }
  if(proc.denseGp)
    addDensePrecision(Q, mag.n, proc.gpPrecision, 1.0 / proc.sigmaSq);
  else if(proc.nngpGraph == NULL)
    magcomp::error("magnitude NNGP graph is not initialized");

  for(int i = 0; i < mag.n; i++){
    double kappa = 0.5 * ((double) mag.total[i] - nbSize);
    double fixedEta = rowDot(mag.X, mag.coef.data(), mag.n, mag.p, i);
    fixedEta += magnitudeProcessEtaExcept(mag, i, termIndex);
    double x = proc.covariate[i];

    rhs[i] += x * (kappa + mag.omega[i] *
                   (logNbSize - std::log(mag.area[i]) - fixedEta));
    if(proc.denseGp)
      Q[i + i * mag.n] += mag.omega[i] * x * x;
    else
      latentDiag[i] = mag.omega[i] * x * x;
  }

  if(proc.denseGp){
    sampleGaussianPrecision(draw, Q, rhs, chol, z, mag.n,
                            "magnitude GP spatial effect", iter);
  }else{
    assembleSparseNngpPrecision(proc.nngpGraph->sparse,
                                proc.nngpGraph->neighborCount,
                                proc.nngpB, proc.nngpF,
                                latentDiag, 1.0 / proc.sigmaSq);
    sampleSparseGaussianPrecision(proc.nngpGraph->sparse, draw, rhs,
                                  "magnitude NNGP spatial effect", iter);
  }
  for(int i = 0; i < mag.n; i++)
    proc.effect[i] = draw[i];

  double ssq;
  if(proc.denseGp){
    ssq = denseVectorScaledSsq(proc.effect, mag.n,
                               proc.gpPrecision);
  }else{
    ssq = nngpVectorScaledSsq(proc.effect, mag.n,
                              proc.nngpGraph->orderIndex,
                              proc.nngpGraph->neighborCount,
                              proc.nngpGraph->neighborIndex,
                              proc.nngpB, proc.nngpF);
  }
  proc.sigmaSq = 1.0 / Rf_rgamma(spatialShape + 0.5 * (double) mag.n,
                                 1.0 / (spatialScale + 0.5 * ssq));
  if(!R_FINITE(proc.sigmaSq) || proc.sigmaSq <= 0.0)
    magcomp::error("invalid magnitude spatial sigma_sq draw");

  if(proc.estimateCovariance){
    proc.covAttempt++;
    int accepted = 0;
    CovarianceState covarianceProposal;
    double jacobian = 0;
    if(proc.covMetrop.propose(proc.covScale, covarianceProposal, jacobian)){
      double logPostCurrent, logPostProposal;
      if(proc.denseGp){
        buildGpPrecision(proc.coords, mag.n, proc.coordDim, covarianceProposal,
                         proc.gpPrecisionProposal, &proc.gpLogDetProposal);
        logPostCurrent = gpVectorLogPrior(proc.effect, mag.n,
                                          proc.gpPrecision,
                                          proc.gpLogDet,
                                          proc.sigmaSq);
        logPostProposal = gpVectorLogPrior(proc.effect, mag.n,
                                           proc.gpPrecisionProposal,
                                           proc.gpLogDetProposal,
                                           proc.sigmaSq);
      }else{
        buildNngpConditionals(
          proc.coords, mag.n, proc.coordDim, proc.nNeighbor, covarianceProposal,
          proc.nngpGraph->orderIndex, proc.nngpGraph->neighborCount,
          proc.nngpGraph->neighborIndex, proc.nngpBProposal,
          proc.nngpFProposal, nngpWorkC, nngpWorkRhs, nngpWorkCvec
        );
        logPostCurrent = nngpVectorLogPrior(proc.effect, mag.n,
                                              proc.nngpGraph->orderIndex,
                                              proc.nngpGraph->neighborCount,
                                              proc.nngpGraph->neighborIndex,
                                              proc.nngpB, proc.nngpF,
                                              proc.sigmaSq);
        logPostProposal = nngpVectorLogPrior(proc.effect, mag.n,
                                               proc.nngpGraph->orderIndex,
                                               proc.nngpGraph->neighborCount,
                                               proc.nngpGraph->neighborIndex,
                                               proc.nngpBProposal,
                                               proc.nngpFProposal,
                                               proc.sigmaSq);
      }
      double logAccept = logPostProposal - logPostCurrent + jacobian;

      if(R_FINITE(logAccept) && std::log(Rf_runif(0.0, 1.0)) < logAccept){
        proc.covMetrop.theta.swap(covarianceProposal.theta);
        if(proc.denseGp){
          proc.gpPrecision.swap(proc.gpPrecisionProposal);
          proc.gpLogDet = proc.gpLogDetProposal;
        }else{
          proc.nngpB.swap(proc.nngpBProposal);
          proc.nngpF.swap(proc.nngpFProposal);
        }
        proc.covAccept++;
        accepted = 1;
      }
    }

    proc.covMetrop.update(accepted, proc.covScale);
  }
}

static void updateMagnitudeIidEffects(MagnitudeState& mag,
                                      int termIndex,
                                      double nbSize,
                                      double iidShape,
                                      double iidScale)
{
  double logNbSize = std::log(nbSize);
  MagnitudeProcessState& proc = mag.process[termIndex];

  if(!proc.isIid)
    return;

  for(int g = 0; g < proc.nGroup; g++){
    double prec = 1.0 / proc.sigmaSq;
    double rhs = 0.0;

    for(int i = 0; i < mag.n; i++){
      if(proc.groupIndex[i] != g)
        continue;

      double kappa = 0.5 * ((double) mag.total[i] - nbSize);
      double fixedEta = rowDot(mag.X, mag.coef.data(), mag.n, mag.p, i);
      fixedEta += magnitudeProcessEtaExcept(mag, i, termIndex);
      double x = proc.covariate[i];
      rhs += x * (kappa + mag.omega[i] *
                  (logNbSize - std::log(mag.area[i]) - fixedEta));
      prec += mag.omega[i] * x * x;
    }

    if(!R_FINITE(prec) || prec <= 0.0 || !R_FINITE(rhs)){
      magcomp::error("invalid scalar Gaussian precision update for magnitude iid "
               "term '%s' group %i", proc.name.c_str(), g + 1);
    }
    double mean = rhs / prec;
    double sd = std::sqrt(1.0 / prec);
    if(!R_FINITE(mean) || !R_FINITE(sd) || sd < 0.0){
      magcomp::error("invalid scalar Gaussian draw for magnitude iid term '%s' "
               "group %i", proc.name.c_str(), g + 1);
    }
    proc.effect[g] = mean + sd * Rf_rnorm(0.0, 1.0);
    if(!R_FINITE(proc.effect[g])){
      magcomp::error("non-finite magnitude iid draw for term '%s' group %i",
               proc.name.c_str(), g + 1);
    }
  }
  double ssq = 0.0;
  for(int g = 0; g < proc.nGroup; g++)
    ssq += proc.effect[g] * proc.effect[g];
  proc.sigmaSq = 1.0 / Rf_rgamma(iidShape + 0.5 * (double) proc.nGroup,
                                 1.0 / (iidScale + 0.5 * ssq));
  if(!R_FINITE(proc.sigmaSq) || proc.sigmaSq <= 0.0)
    magcomp::error("invalid magnitude iid sigma_sq draw for term '%s'", proc.name.c_str());
}

static void updateAllocationBlock(AllocationBlockState& block,
                                  double *Q,
                                  double *rhs,
                                  double *draw,
                                  double *chol,
	                                  double *z,
	                                  BayesLogit_rpg_hybrid_t rpg,
	                                  BayesLogit_rpg_gamma_t rpgGamma,
                                  CompositionLink compositionLink,
	                                  int iter)
{
  double theta, omega, kappa, xValue;
  int idx, trials, curveOrder, coefIdx;
  double curvePriorPrec, curveSmoothPrec;

	  if(block.aggregateInterceptOnly){
    double probitProfileStart = profileNowSeconds();
    long long probitDraws = 0;
	    curveOrder = block.interceptSmoothOrder;
    curvePriorPrec = block.interceptPriorPrec;
    curveSmoothPrec = block.interceptSmoothPrec;

    if(curveOrder == 0){
      for(int c = 0; c < block.coefStick; c++){
        double rhsValue = 0.0;
        double prec = curvePriorPrec;
        trials = block.aggregateTrials[c];

	        if(trials > 0){
	          theta = block.coef[c];
          if(compositionLink == COMPOSITION_LINK_PROBIT){
            rhsValue += drawProbitBinomialLatentSum(
              theta, block.aggregateCount[c], trials,
              "aggregated probit composition update", iter);
            probitDraws += trials;
            omega = (double) trials;
          }else{
	            omega = drawPolyaGamma(rpg, rpgGamma, (double) trials, theta,
	                                   "aggregated composition update", iter);
	            kappa = (double) block.aggregateCount[c] -
	              0.5 * (double) trials;
	            rhsValue += kappa;
          }
	          prec += omega;
	        }

        if(!R_FINITE(prec) || prec <= 0.0 || !R_FINITE(rhsValue)){
          magcomp::error("invalid scalar Gaussian precision update for aggregated "
                   "composition intercept row %i", c + 1);
        }
        double mean = rhsValue / prec;
        double sd = std::sqrt(1.0 / prec);
        if(!R_FINITE(mean) || !R_FINITE(sd) || sd < 0.0){
          magcomp::error("invalid scalar Gaussian draw for aggregated composition "
                   "intercept row %i", c + 1);
        }
        block.coef[c] = mean + sd * Rf_rnorm(0.0, 1.0);
        if(!R_FINITE(block.coef[c])){
          magcomp::error("non-finite aggregated composition intercept draw for "
                   "row %i", c + 1);
        }
      }
      if(compositionLink == COMPOSITION_LINK_PROBIT && gSamplerProfile != NULL){
        gSamplerProfile->probitCompositionSeconds +=
          profileNowSeconds() - probitProfileStart;
        gSamplerProfile->probitCompositionDraws += probitDraws;
      }
	      return;
    }

    for(int j = 0; j < block.coefStick; j++){
      rhs[j] = 0.0;
      for(int h = 0; h < block.coefStick; h++)
        Q[j + h * block.coefStick] = 0.0;
    }
    if(block.parentSpecificCoef){
      int nParent = block.coefStick / block.Jstick;
      for(int a = 0; a < nParent; a++){
        addCurvePriorBlock(Q, a * block.Jstick, block.Jstick,
                           block.coefStick, curveOrder,
                           curvePriorPrec, curveSmoothPrec);
      }
    }else{
      addCurvePrior(Q, block.coefStick, curveOrder, curvePriorPrec,
                    curveSmoothPrec);
    }

    for(int c = 0; c < block.coefStick; c++){
      trials = block.aggregateTrials[c];
      if(trials <= 0)
        continue;

	      theta = block.coef[c];
      if(compositionLink == COMPOSITION_LINK_PROBIT){
        omega = (double) trials;
        rhs[c] += drawProbitBinomialLatentSum(
          theta, block.aggregateCount[c], trials,
          "aggregated probit composition update", iter);
        probitDraws += trials;
      }else{
	        omega = drawPolyaGamma(rpg, rpgGamma, (double) trials, theta,
	                               "aggregated composition update", iter);
	        kappa = (double) block.aggregateCount[c] -
          0.5 * (double) trials;
	        rhs[c] += kappa;
      }
	      Q[c + c * block.coefStick] += omega;
	    }

    sampleGaussianPrecision(draw, Q, rhs, chol, z, block.coefStick,
                            "aggregated composition intercept curve", iter);
    for(int j = 0; j < block.coefStick; j++)
      block.coef[j] = draw[j];
    if(compositionLink == COMPOSITION_LINK_PROBIT && gSamplerProfile != NULL){
      gSamplerProfile->probitCompositionSeconds +=
        profileNowSeconds() - probitProfileStart;
      gSamplerProfile->probitCompositionDraws += probitDraws;
    }
	    return;
	  }

  {
    double probitProfileStart = profileNowSeconds();
    long long probitDraws = 0;
	  for(int j = 0; j < block.Jstick; j++){
	    for(int i = 0; i < block.n; i++){
	      idx = i + j * block.n;
	      trials = block.trials[idx];
	      if(trials <= 0){
	        block.omega[idx] = 0.0;
          block.workingResponse[idx] = 0.0;
	        continue;
	      }

	      theta = blockLinearPredictor(block, i, j);
      if(compositionLink == COMPOSITION_LINK_PROBIT){
        block.omega[idx] = (double) trials;
        block.workingResponse[idx] = drawProbitBinomialLatentSum(
          theta, block.count[idx], trials, "probit composition update", iter);
        probitDraws += trials;
      }else{
	        omega = drawPolyaGamma(rpg, rpgGamma, (double) trials, theta,
	                               "composition update", iter);
	        block.omega[idx] = omega;
        block.workingResponse[idx] =
          (double) block.count[idx] - 0.5 * (double) trials;
      }
	    }
	  }
    if(compositionLink == COMPOSITION_LINK_PROBIT && gSamplerProfile != NULL){
      gSamplerProfile->probitCompositionSeconds +=
        profileNowSeconds() - probitProfileStart;
      gSamplerProfile->probitCompositionDraws += probitDraws;
    }
  }

  for(int l = 0; l < block.p; l++){
    curveOrder = (l == 0) ? block.interceptSmoothOrder :
      block.slopeSmoothOrder;
    curvePriorPrec = (l == 0) ? block.interceptPriorPrec :
      block.betaPriorPrec;
    curveSmoothPrec = (l == 0) ? block.interceptSmoothPrec :
      block.slopeSmoothPrec;

    for(int j = 0; j < block.coefStick; j++){
      rhs[j] = 0.0;
      for(int h = 0; h < block.coefStick; h++)
        Q[j + h * block.coefStick] = 0.0;
    }
    if(block.parentSpecificCoef){
      int nParent = block.coefStick / block.Jstick;
      for(int a = 0; a < nParent; a++){
        addCurvePriorBlock(Q, a * block.Jstick, block.Jstick,
                           block.coefStick, curveOrder,
                           curvePriorPrec, curveSmoothPrec);
      }
    }else{
      addCurvePrior(Q, block.coefStick, curveOrder, curvePriorPrec,
                    curveSmoothPrec);
    }

    for(int j = 0; j < block.Jstick; j++){
      for(int i = 0; i < block.n; i++){
        idx = i + j * block.n;
        coefIdx = blockCoefIndex(block, i, j);
        trials = block.trials[idx];
        if(trials <= 0)
          continue;

        theta = blockFactorEta(block, i, j);
        for(int h = 0; h < block.p; h++){
          if(h != l)
            theta += block.X[i + h * block.n] *
              block.coef[coefIdx + h * block.coefStick];
        }

	        xValue = block.X[i + l * block.n];
        rhs[coefIdx] += xValue *
          (block.workingResponse[idx] - block.omega[idx] * theta);
	        Q[coefIdx + coefIdx * block.coefStick] +=
	          block.omega[idx] * xValue * xValue;
      }
    }

    sampleGaussianPrecision(draw, Q, rhs, chol, z, block.coefStick,
                            "composition coefficient curve", iter);
    for(int j = 0; j < block.coefStick; j++)
      block.coef[j + l * block.coefStick] = draw[j];
  }
}

static void updateAllocationFactor(AllocationBlockState& block,
                                   int factorIndex,
                                   double loadingPriorPrec,
                                   double *Q,
                                   double *rhs,
                                   double *draw,
                                   double *chol,
                                   double *z,
                                   double *latentDiag,
                                   double *nngpWorkC,
                                   double *nngpWorkRhs,
                                   double *nngpWorkCvec,
                                   int iter,
                                   int nBurn)
{
  if(!block.hasFactor)
    return;
  if(factorIndex < 0 || factorIndex >= (int) block.factor.size())
    magcomp::error("internal error: composition factor index out of range");

  CompositionFactorState& factor = block.factor[factorIndex];

  for(int c = 0; c < block.coefStick; c++){
    int j = block.parentSpecificCoef ? c % block.Jstick : c;
    for(int k = 0; k < factor.K; k++){
      if(!loadingIsFree(c, k, factor.anchorIndex)){
        double constrainedValue = 0.0;
        loadingIsConstrained(c, k, factor.anchorIndex, constrainedValue);
        factor.loading[c + k * block.coefStick] = constrainedValue;
        continue;
      }

      double prec = loadingPriorPrec;
      double rhsValue = 0.0;
      for(int i = 0; i < block.n; i++){
        int idx = i + j * block.n;
        int coefIdx = blockCoefIndex(block, i, j);
        if(coefIdx != c)
          continue;
        int trials = block.trials[idx];
        if(trials <= 0)
          continue;

        int g = factor.groupIndex[i];
        double fixedEta = 0.0;
        for(int l = 0; l < block.p; l++)
          fixedEta += block.X[i + l * block.n] *
            block.coef[coefIdx + l * block.coefStick];

        double otherFactorEta = 0.0;
        otherFactorEta += blockFactorEtaExcept(block, i, j, factorIndex);
        for(int kk = 0; kk < factor.K; kk++){
          if(kk != k){
            otherFactorEta += factor.covariate[i] *
              factor.loading[coefIdx + kk * block.coefStick] *
              factor.factorScore[g + kk * factor.nGroup];
          }
        }

	        double x = factor.covariate[i] *
	          factor.factorScore[g + k * factor.nGroup];
        rhsValue += x * (block.workingResponse[idx] -
                         block.omega[idx] * (fixedEta + otherFactorEta));
	        prec += block.omega[idx] * x * x;
      }

      if(!R_FINITE(prec) || prec <= 0.0 || !R_FINITE(rhsValue)){
        magcomp::error("invalid scalar Gaussian precision update for composition "
                 "factor loading term '%s', row %i, factor %i",
                 factor.name.c_str(), c + 1, k + 1);
      }
      double mean = rhsValue / prec;
      double sd = std::sqrt(1.0 / prec);
      if(!R_FINITE(mean) || !R_FINITE(sd) || sd < 0.0){
        magcomp::error("invalid scalar Gaussian draw for composition factor loading "
                 "term '%s', row %i, factor %i",
                 factor.name.c_str(), c + 1, k + 1);
      }
      factor.loading[c + k * block.coefStick] =
        c == factor.anchorIndex[k] ? rnormPositive(mean, sd) :
        mean + sd * Rf_rnorm(0.0, 1.0);
      if(!R_FINITE(factor.loading[c + k * block.coefStick])){
        magcomp::error("non-finite composition factor loading draw for term '%s', "
                 "row %i, factor %i", factor.name.c_str(), c + 1, k + 1);
      }
    }
  }

  enforceLoadingConstraints(factor.loading, block.coefStick, factor.K,
                            factor.anchorIndex);

  if(factor.isSpatial){
    if(!factor.denseGp && factor.nngpGraph == NULL)
      magcomp::error("composition factor NNGP graph is not initialized");

    for(int k = 0; k < factor.K; k++){
      for(int g = 0; g < factor.nGroup; g++){
        rhs[g] = 0.0;
        latentDiag[g] = 0.0;
        if(factor.denseGp){
          for(int gg = 0; gg < factor.nGroup; gg++)
            Q[g + gg * factor.nGroup] = 0.0;
        }
      }
      if(factor.denseGp)
        addDensePrecision(Q, factor.nGroup, factor.gpPrecision[k], 1.0);

      for(int i = 0; i < block.n; i++){
        int g = factor.groupIndex[i];
        for(int j = 0; j < block.Jstick; j++){
          int idx = i + j * block.n;
          int coefIdx = blockCoefIndex(block, i, j);
          int trials = block.trials[idx];
          if(trials <= 0)
            continue;

          double fixedEta = 0.0;
          for(int l = 0; l < block.p; l++)
            fixedEta += block.X[i + l * block.n] *
              block.coef[coefIdx + l * block.coefStick];

          double otherFactorEta = 0.0;
          otherFactorEta += blockFactorEtaExcept(block, i, j, factorIndex);
          for(int kk = 0; kk < factor.K; kk++){
            if(kk != k){
              otherFactorEta += factor.covariate[i] *
                factor.loading[coefIdx + kk * block.coefStick] *
                factor.factorScore[g + kk * factor.nGroup];
            }
          }

	          double x = factor.covariate[i] *
	            factor.loading[coefIdx + k * block.coefStick];
          rhs[g] += x * (block.workingResponse[idx] -
                         block.omega[idx] * (fixedEta + otherFactorEta));
	          if(factor.denseGp)
	            Q[g + g * factor.nGroup] += block.omega[idx] * x * x;
          else
            latentDiag[g] += block.omega[idx] * x * x;
        }
      }

      if(factor.denseGp){
        sampleGaussianPrecision(draw, Q, rhs, chol, z, factor.nGroup,
                                "composition GP factor score", iter);
      }else{
        assembleSparseNngpPrecision(factor.nngpGraph->sparse,
                                    factor.nngpGraph->neighborCount,
                                    factor.nngpB[k],
                                    factor.nngpF[k],
                                    latentDiag, 1.0);
        sampleSparseGaussianPrecision(factor.nngpGraph->sparse,
                                      draw, rhs,
                                      "composition NNGP factor score",
                                      iter);
      }
      for(int g = 0; g < factor.nGroup; g++)
        factor.factorScore[g + k * factor.nGroup] = draw[g];
    }

    enforceLoadingConstraints(factor.loading, block.coefStick, factor.K,
                              factor.anchorIndex);

    if(factor.estimateCovariance){
      for(int k = 0; k < factor.K; k++){
        if(factor.covMetrop[k].free.empty()) continue;
        factor.covAttempt[k]++;
        int accepted = 0;
        CovarianceState covarianceProposal;
        double jacobian = 0;
        if(factor.covMetrop[k].propose(factor.covScale[k], covarianceProposal, jacobian)){
          double logPostCurrent, logPostProposal;
          if(factor.denseGp){
            buildGpPrecision(factor.coords, factor.nGroup,
                             factor.coordDim, covarianceProposal,
                             factor.gpPrecisionProposal[k],
                             &factor.gpLogDetProposal[k]);
            logPostCurrent = gpFactorColumnLogPrior(
              factor.factorScore, factor.nGroup, k,
              factor.gpPrecision[k], factor.gpLogDet[k]);
            logPostProposal = gpFactorColumnLogPrior(
              factor.factorScore, factor.nGroup, k,
              factor.gpPrecisionProposal[k], factor.gpLogDetProposal[k]);
          }else{
            buildNngpConditionals(
              factor.coords, factor.nGroup, factor.coordDim,
              factor.nNeighbor, covarianceProposal, factor.nngpGraph->orderIndex,
              factor.nngpGraph->neighborCount, factor.nngpGraph->neighborIndex,
              factor.nngpBProposal[k], factor.nngpFProposal[k],
              nngpWorkC, nngpWorkRhs, nngpWorkCvec
            );
            logPostCurrent = nngpFactorColumnLogPrior(
                factor.factorScore, factor.nGroup, k,
                factor.nngpGraph->orderIndex,
                factor.nngpGraph->neighborCount,
                factor.nngpGraph->neighborIndex,
                factor.nngpB[k], factor.nngpF[k]);
            logPostProposal = nngpFactorColumnLogPrior(
                factor.factorScore, factor.nGroup, k,
                factor.nngpGraph->orderIndex,
                factor.nngpGraph->neighborCount,
                factor.nngpGraph->neighborIndex,
                factor.nngpBProposal[k], factor.nngpFProposal[k]);
          }
          double logAccept = logPostProposal - logPostCurrent + jacobian;

          if(R_FINITE(logAccept) &&
             std::log(Rf_runif(0.0, 1.0)) < logAccept){
            factor.covMetrop[k].theta.swap(covarianceProposal.theta);
            if(factor.denseGp){
              factor.gpPrecision[k].swap(factor.gpPrecisionProposal[k]);
              factor.gpLogDet[k] = factor.gpLogDetProposal[k];
            }else{
              factor.nngpB[k].swap(factor.nngpBProposal[k]);
              factor.nngpF[k].swap(factor.nngpFProposal[k]);
            }
            factor.covAccept[k]++;
            accepted = 1;
          }
        }

        factor.covMetrop[k].update(accepted, factor.covScale[k]);
      }
    }
    return;
  }

  for(int g = 0; g < factor.nGroup; g++){
    for(int k = 0; k < factor.K; k++){
      rhs[k] = 0.0;
      for(int kk = 0; kk < factor.K; kk++)
        Q[k + kk * factor.K] = 0.0;
      Q[k + k * factor.K] = 1.0;
    }

    for(int i = 0; i < block.n; i++){
      if(factor.groupIndex[i] != g)
        continue;

      for(int j = 0; j < block.Jstick; j++){
        int idx = i + j * block.n;
        int coefIdx = blockCoefIndex(block, i, j);
        int trials = block.trials[idx];
        if(trials <= 0)
          continue;

        double fixedEta = 0.0;
        for(int l = 0; l < block.p; l++)
          fixedEta += block.X[i + l * block.n] *
            block.coef[coefIdx + l * block.coefStick];
        fixedEta += blockFactorEtaExcept(block, i, j, factorIndex);
	        for(int k = 0; k < factor.K; k++){
	          double xk = factor.covariate[i] *
	            factor.loading[coefIdx + k * block.coefStick];
          rhs[k] += xk * (block.workingResponse[idx] -
                          block.omega[idx] * fixedEta);
	          for(int kk = 0; kk < factor.K; kk++){
            double xkk = factor.covariate[i] *
              factor.loading[coefIdx + kk * block.coefStick];
            Q[k + kk * factor.K] += block.omega[idx] * xk * xkk;
          }
        }
      }
    }

    sampleGaussianPrecision(draw, Q, rhs, chol, z, factor.K,
                            "composition iid factor score", iter);
    for(int k = 0; k < factor.K; k++)
      factor.factorScore[g + k * factor.nGroup] = draw[k];
  }

  enforceLoadingConstraints(factor.loading, block.coefStick, factor.K,
                            factor.anchorIndex);
}

static void fillBlockProportions(AllocationBlockState& block,
                                 CompositionLink compositionLink)
{
  double remaining, eta, pi;

  for(int i = 0; i < block.n; i++){
    remaining = 1.0;
	    for(int j = 0; j < block.Jstick; j++){
	      eta = blockLinearPredictor(block, i, j);
      pi = compositionInverseLink(eta, compositionLink);
	      block.proportion[i + j * block.n] = remaining * pi;
	      remaining *= 1.0 - pi;
    }
    block.proportion[i + (block.J - 1) * block.n] = remaining;
  }
}

static SEXP runSampler(SEXP config, bool &rngActive)
{
  SEXP yS, yDimS, areaS, xMagS, xCompS, priorsS, smoothS;
  SEXP controlS, startingS, magnitudeProcessS, metropolisS;
  SEXP allocationBlocksS, compositionBlocksS, allocationBlockS;
  SEXP compositionBlockS, blockCountsS, blockNamesS;
  SEXP magnitudeCoefficientSamplesS, piSamplesS;
  SEXP pSamplesS, magnitudeSamplesS, countMeanS, densityMeanS;
  SEXP pMeanS, totalCountMeanS, totalDensityMeanS;
  SEXP nbSizeSamplesS, nbSizeAcceptanceS, nbSizeTuningS;
  SEXP nbSizeBatchS, nbSizeTuningHistoryS;
  SEXP blockCoefSamplesS, blockCoefNamesS;
  SEXP magProcessSamplesS, magProcessNamesS;
  SEXP magProcessSigmaSamplesS, magProcessSigmaNamesS;
  SEXP magProcessPhiSamplesS, magProcessPhiNamesS;
  SEXP magProcessPhiAcceptanceS, magProcessPhiTuningS;
  SEXP magProcessPhiBatchS, magProcessPhiTuningHistoryS;
  SEXP blockLoadingSamplesS, blockFactorSamplesS, blockFactorPhiSamplesS;
  SEXP blockLoadingNamesS, blockFactorNamesS, blockFactorPhiNamesS;
  SEXP blockFactorPhiScalarNamesS;
  SEXP blockFactorPhiAcceptanceS, blockFactorPhiTuningS;
  SEXP blockFactorPhiBatchS, blockFactorPhiTuningHistoryS;
	  SEXP outS, namesS, profileStageS, profileStageNamesS;
  SEXP profilePgS, profilePgNamesS;
  MagnitudeState mag = {};
  MagnitudeNbSizeState nb;
  std::vector<AllocationBlockState> blocks;
  std::vector<SEXP> blockCoefArrayS;
  std::vector<SEXP> magProcessArrayS;
  std::vector<SEXP> magProcessSigmaArrayS;
  std::vector<SEXP> magProcessPhiArrayS;
  std::vector<SEXP> blockLoadingArrayS;
  std::vector<SEXP> blockFactorArrayS;
  std::vector<SEXP> blockFactorPhiArrayS;
  std::vector<int> spatialProcessIndex;
  std::vector<int> factorBlockIndex;
  std::vector<int> factorTermIndex;
  std::vector<int> spatialFactorBlockIndex;
  std::vector<int> spatialFactorTermIndex;
  std::vector<int> childRowByUnitParent;
  LatentBackendRegistry latentRegistry;
  int nSamples, nBurn, nThin, nReport, nIter, saveIdx;
  int nUnits, nBlocks, outputJ, stickTotal, maxDim, nProtect, idx;
  int parentJ, childJ, row, stickCol;
  int nFactorBlocks, nSpatialFactorBlocks, nSpatialFactorPhi;
  int metropolisBatchLength, maxNeighbor, nRegisteredNngpGraphs;
  int expectedNngpGraphs;
  double nbSize, betaPriorPrec;
	  double nbSizeTuning, metropolisTargetAccept;
		  double theta, magnitude, density, count, profileStageStart, profileLoopStart;
	  double pValue, piValue;
	  bool verbose, metropolisAdapt;
  CompositionLink compositionLink;
		  BayesLogit_rpg_hybrid_t rpg;
		  BayesLogit_rpg_gamma_t rpgGamma;
  SamplerProfile profile = {};
  profile.samplerStart = profileNowSeconds();

  yS = getListElement(config, "y");
  allocationBlocksS = getOptionalListElement(config, "allocation_blocks");
  compositionBlocksS = getOptionalListElement(config, "composition_blocks");
  areaS = getOptionalListElement(config, "area");
  xMagS = getOptionalListElement(config, "X_magnitude");
  const bool fitMagnitude = !Rf_isNull(xMagS);
  yDimS = Rf_getAttrib(yS, R_DimSymbol);
  if(Rf_length(yDimS) != 2 || INTEGER(yDimS)[0] < 1 || INTEGER(yDimS)[1] < 1)
    magcomp::error("y must be a nonempty count matrix");
  nUnits = INTEGER(yDimS)[0];
  if(TYPEOF(allocationBlocksS) != VECSXP)
    magcomp::error("allocation_blocks must be a list");
  nBlocks = Rf_length(allocationBlocksS);
  if(nBlocks > 2 || (!fitMagnitude && nBlocks == 0))
    magcomp::error("specify magnitude, composition with one or two allocation blocks, or both");
  const bool fitComposition = nBlocks > 0;
  priorsS = getListElement(config, "priors");
  smoothS = getListElement(config, "smooth");
  startingS = getListElement(config, "starting");
  controlS = getListElement(config, "control");
  metropolisS = getListElement(controlS, "metropolis");
  magnitudeProcessS = getOptionalListElement(config, "magnitude_process");

  nIter = getListInteger(controlS, "n_samples");
  nBurn = getListInteger(controlS, "n_burn");
  nThin = getListInteger(controlS, "n_thin");
  nReport = getListInteger(controlS, "n_report");
  nbSize = getListDouble(controlS, "nb_size");
  nbSizeTuning = getListDouble(controlS, "nb_size_tuning");
  metropolisAdapt = getListLogical(metropolisS, "adapt");
	  metropolisBatchLength = getListInteger(metropolisS, "batch_length");
	  metropolisTargetAccept = getListDouble(metropolisS, "target_accept");
	  verbose = getListLogical(controlS, "verbose");
  compositionLink = getCompositionLink(controlS);
	  betaPriorPrec = fitMagnitude ? 1.0 / getListDouble(priorsS, "coefficient_variance") : 0.0;

  if(nIter < 1 || nBurn < 0 || nThin < 1 || nIter <= nBurn)
    magcomp::error("invalid MCMC control values");
  nSamples = (nIter - nBurn) / nThin;
  if(nSamples < 1)
    magcomp::error("n_thin must leave at least one retained sample after burn-in");
  if(!R_FINITE(nbSize) || nbSize <= 0.0)
    magcomp::error("nb_size must be positive and finite");
  if(!R_FINITE(nbSizeTuning) || nbSizeTuning <= 0.0)
    magcomp::error("nb_size_tuning must be positive and finite");
  if(metropolisBatchLength < 1)
    magcomp::error("metropolis$batch_length must be positive");
  if(!R_FINITE(metropolisTargetAccept) ||
     metropolisTargetAccept < 0.0 || metropolisTargetAccept >= 1.0)
    magcomp::error("metropolis$target_accept must be between zero and one");


  if(nIter >= INT_MAX)
    magcomp::error("total MCMC iterations exceed the supported integer indexing limit");
  nProtect = 0;
  nRegisteredNngpGraphs = 0;
  initLatentBackendRegistry(latentRegistry);
  LatentBackendGuard backendGuard{latentRegistry};
  initMagnitudeNbSizeState(nb);
  nb.value = nbSize;
  nb.tuning = nbSizeTuning;

  if(fitMagnitude) {
    {
      SEXP nbPriorS = getOptionalListElement(priorsS, "nb_size");
      if(!Rf_isNull(nbPriorS)){
        if(Rf_length(nbPriorS) != 2)
          magcomp::error("priors$nb_size must have length two");
        nb.estimate = true;
        nb.shape = REAL(nbPriorS)[0];
        nb.rate = REAL(nbPriorS)[1];
        if(!R_FINITE(nb.shape) || !R_FINITE(nb.rate) ||
           nb.shape <= 0.0 || nb.rate <= 0.0){
          magcomp::error("priors$nb_size must contain positive finite Gamma shape "
                   "and rate values");
        }
        prepareScalarMetropolisState(nb.metrop, metropolisAdapt, nBurn,
                                     metropolisBatchLength,
                                     metropolisTargetAccept);
      }
    }
    {
      SEXP nbStartS = getOptionalListElement(startingS, "nb_size");
      if(!Rf_isNull(nbStartS)){
        if(!nb.estimate){
          magcomp::error("starting$nb_size is only used when priors$nb_size "
                   "requests size estimation; set control$nb_size for a "
                   "fixed size");
        }
        nb.value = magcomp::guarded_Rf_asReal(nbStartS);
      }
      if(!R_FINITE(nb.value) || nb.value <= 0.0)
        magcomp::error("starting$nb_size must be positive and finite");
    }

    loadMagnitudeState(mag, yS, areaS, xMagS);
    loadMagnitudeTotals(mag, yS);
    loadMagnitudeProcessTerm(mag, magnitudeProcessS);
    if(nb.estimate && !mag.process.empty()){
      magcomp::error("estimating nb_size is currently supported only for magnitude "
               "models without random-effect process terms");
    }

    initializeMagnitude(mag, startingS);
    for(int mp = 0; mp < (int) mag.process.size(); mp++){
      if(mag.process[mp].isSpatial){
        spatialProcessIndex.push_back(mp);
        mag.process[mp].covMetrop.prepare(metropolisAdapt, nBurn, nIter,
                                         metropolisBatchLength, metropolisTargetAccept);
      }
    }
  }
  if(fitComposition && (TYPEOF(compositionBlocksS) != VECSXP ||
      Rf_length(compositionBlocksS) != nBlocks))
    magcomp::error("composition_blocks must match allocation_blocks");
  blockNamesS = Rf_getAttrib(allocationBlocksS, R_NamesSymbol);
  blocks.resize(nBlocks);
  blockCoefArrayS.resize(nBlocks);
  for(int b = 0; b < nBlocks; b++){
    allocationBlockS = VECTOR_ELT(allocationBlocksS, b);
    blockCountsS = getListElement(allocationBlockS, "counts");
    const char *blockType = getListString(allocationBlockS, "type");
    if(b == 0 && std::strcmp(blockType, "marginal") != 0)
      magcomp::error("first allocation block must be marginal");
    if(b == 1 && std::strcmp(blockType, "conditional") != 0)
      magcomp::error("second allocation block must be conditional");
    compositionBlockS = VECTOR_ELT(compositionBlocksS, b);
    xCompS = getListElement(compositionBlockS, "X");
    loadAllocationBlock(blocks[b], allocationBlockS, blockCountsS, xCompS,
                        smoothS, getListElement(compositionBlockS, "coefficient_priors"));
    if(b == 0 && blocks[b].n != nUnits)
      magcomp::error("marginal allocation block must have one row per unit");
    loadBlockIndex(blocks[b], allocationBlockS, nUnits);
    loadAllocationFactorTerms(blocks[b], compositionBlockS, nUnits);
    loadAllocationCounts(blocks[b], blockCountsS);
    configureAllocationAggregation(blocks[b]);
    const char *blockName = NULL;
    if(!Rf_isNull(blockNamesS) && Rf_length(blockNamesS) == nBlocks)
      blockName = CHAR(STRING_ELT(blockNamesS, b));
    initializeAllocationBlock(blocks[b], startingS, blockName,
                              compositionLink);
    initializeAllocationFactor(blocks[b], startingS, blockName);
  }

  nFactorBlocks = 0;
  nSpatialFactorBlocks = 0;
  nSpatialFactorPhi = 0;
  for(int b = 0; b < nBlocks; b++){
    for(int f = 0; f < (int) blocks[b].factor.size(); f++){
      factorBlockIndex.push_back(b);
      factorTermIndex.push_back(f);
      nFactorBlocks++;
      if(blocks[b].factor[f].isSpatial){
        spatialFactorBlockIndex.push_back(b);
        spatialFactorTermIndex.push_back(f);
        nSpatialFactorBlocks++;
        nSpatialFactorPhi += blocks[b].factor[f].K;
        for(int k = 0; k < blocks[b].factor[f].K; k++){
          blocks[b].factor[f].covMetrop[k].prepare(metropolisAdapt, nBurn, nIter,
                                                 metropolisBatchLength, metropolisTargetAccept);
        }
      }
    }
  }

  yDimS = Rf_getAttrib(yS, R_DimSymbol);
  outputJ = fitComposition ? INTEGER(yDimS)[1] : 1;
  stickTotal = fitComposition ? blocks[0].coefStick : 0;
  if(nBlocks == 1){
    if(outputJ != blocks[0].J)
      magcomp::error("output count matrix does not match allocation block");
  }else if(nBlocks == 2){
    parentJ = blocks[0].J;
    childJ = blocks[1].J;
    if(outputJ != parentJ * childJ)
      magcomp::error("output count matrix does not match conditional blocks");
    stickTotal += blocks[1].coefStick;
    childRowByUnitParent.assign(nUnits * parentJ, -1);
    for(int i = 0; i < blocks[1].n; i++){
      if(blocks[1].parentIndex[i] < 0 ||
         blocks[1].parentIndex[i] >= parentJ)
        magcomp::error("conditional block parent_index is out of range");
      idx = blocks[1].unitIndex[i] +
        blocks[1].parentIndex[i] * nUnits;
      childRowByUnitParent[idx] = i;
    }
    for(int i = 0; i < nUnits * parentJ; i++){
      if(childRowByUnitParent[i] < 0)
        magcomp::error("conditional block is missing a unit-parent row");
    }
  }

  rpg = magcomp::rCall([]() { return BayesLogit_rpg_hybrid(); });
  rpgGamma = magcomp::rCall([]() { return BayesLogit_rpg_gamma(); });

  maxDim = mag.p;
  for(int mp = 0; mp < (int) mag.process.size(); mp++){
    if(mag.process[mp].isSpatial)
      maxDim = std::max(maxDim, nUnits);
  }
  for(int b = 0; b < nBlocks; b++){
    maxDim = std::max(maxDim, blocks[b].coefStick);
    for(int f = 0; f < (int) blocks[b].factor.size(); f++){
      CompositionFactorState& factor = blocks[b].factor[f];
      maxDim = std::max(maxDim, factor.K);
      if(factor.isSpatial)
        maxDim = std::max(maxDim, factor.nGroup);
    }
  }
  if((double) maxDim * maxDim > INT_MAX)
    magcomp::error("sampler workspace exceeds the supported integer indexing limit");
  std::vector<double> Q(maxDim * maxDim, 0.0);
  std::vector<double> rhs(maxDim, 0.0);
  std::vector<double> draw(maxDim, 0.0);
  std::vector<double> chol(maxDim * maxDim, 0.0);
  std::vector<double> z(maxDim, 0.0);
  std::vector<double> latentDiag(std::max(1, maxDim), 0.0);
  maxNeighbor = 0;
  expectedNngpGraphs = 0;
  for(int mp = 0; mp < (int) mag.process.size(); mp++){
    MagnitudeProcessState& proc = mag.process[mp];
    if(proc.isSpatial && !proc.denseGp){
      maxNeighbor = std::max(maxNeighbor, proc.nNeighbor);
      expectedNngpGraphs++;
    }
  }
  for(int b = 0; b < nBlocks; b++){
    for(int f = 0; f < (int) blocks[b].factor.size(); f++){
      CompositionFactorState& factor = blocks[b].factor[f];
      if(factor.isSpatial && !factor.denseGp){
        maxNeighbor = std::max(maxNeighbor, factor.nNeighbor);
        expectedNngpGraphs++;
      }
    }
  }
  std::vector<double> nngpWorkC(std::max(1, maxNeighbor * maxNeighbor), 0.0);
  std::vector<double> nngpWorkRhs(std::max(1, maxNeighbor), 0.0);
  std::vector<double> nngpWorkCvec(std::max(1, maxNeighbor), 0.0);

  reserveLatentBackendRegistry(latentRegistry, expectedNngpGraphs);

  for(int mp = 0; mp < (int) mag.process.size(); mp++){
    MagnitudeProcessState& proc = mag.process[mp];
    if(!proc.isSpatial)
      continue;
    if(proc.denseGp){
      buildGpPrecision(proc.coords, nUnits, proc.coordDim, proc.covMetrop,
                       proc.gpPrecision, &proc.gpLogDet);
    }else{
      proc.nngpGraph = registerNngpGraphBackend(
        latentRegistry, proc.name.c_str(), nUnits, proc.nNeighbor,
        proc.orderS, proc.neighborCountS, proc.neighborS);
      buildNngpConditionals(proc.coords, nUnits, proc.coordDim,
                            proc.nNeighbor, proc.covMetrop,
                            proc.nngpGraph->orderIndex,
                            proc.nngpGraph->neighborCount,
                            proc.nngpGraph->neighborIndex,
                            proc.nngpB, proc.nngpF,
                            nngpWorkC.data(), nngpWorkRhs.data(),
                            nngpWorkCvec.data());
    }
  }
  for(int b = 0; b < nBlocks; b++){
    for(int f = 0; f < (int) blocks[b].factor.size(); f++){
      CompositionFactorState& factor = blocks[b].factor[f];
      if(!factor.isSpatial)
        continue;
      if(factor.denseGp){
        for(int k = 0; k < factor.K; k++){
          buildGpPrecision(factor.coords, factor.nGroup,
                           factor.coordDim, factor.covMetrop[k],
                           factor.gpPrecision[k],
                           &factor.gpLogDet[k]);
        }
      }else{
        factor.nngpGraph = registerNngpGraphBackend(
          latentRegistry, factor.name.c_str(), factor.nGroup,
          factor.nNeighbor, factor.orderS,
          factor.neighborCountS, factor.neighborS);
        for(int k = 0; k < factor.K; k++){
          buildNngpConditionals(factor.coords, factor.nGroup,
                                factor.coordDim,
                                factor.nNeighbor,
                                factor.covMetrop[k],
                                factor.nngpGraph->orderIndex,
                                factor.nngpGraph->neighborCount,
                                factor.nngpGraph->neighborIndex,
                                factor.nngpB[k],
                                factor.nngpF[k],
                                nngpWorkC.data(), nngpWorkRhs.data(),
                                nngpWorkCvec.data());
        }
      }
    }
  }

  // Check every draw-array size before allocating any potentially large one.
  checkArraySize(mag.p, nSamples);
  checkArraySize(nUnits, stickTotal, nSamples);
  checkArraySize(nUnits, outputJ, nSamples);
  for(const auto &proc : mag.process)
    checkArraySize(proc.nGroup, nSamples);
  for(const auto &block : blocks){
    checkArraySize(block.coefStick, block.p, nSamples);
    for(const auto &factor : block.factor){
      checkArraySize(block.coefStick, factor.K, nSamples);
      checkArraySize(factor.nGroup, factor.K, nSamples);
    }
  }

  magnitudeCoefficientSamplesS = R_NilValue;
  if(fitMagnitude){
    magnitudeCoefficientSamplesS = PROTECT(allocMatrixReal(mag.p, nSamples));
    nProtect++;
  }
  blockCoefSamplesS = R_NilValue;
  if(fitComposition) {
    blockCoefSamplesS = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, nBlocks));
    nProtect++;
    blockCoefNamesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, nBlocks));
    nProtect++;
    for(int b = 0; b < nBlocks; b++){
      blockCoefArrayS[b] = PROTECT(allocArrayReal3(blocks[b].coefStick,
                                                   blocks[b].p,
                                                   nSamples));
      nProtect++;
      if(!Rf_isNull(blockNamesS) && Rf_length(blockNamesS) == nBlocks){
        SET_STRING_ELT(blockCoefNamesS, b, STRING_ELT(blockNamesS, b));
      }else{
        char tmp[32];
        std::snprintf(tmp, sizeof(tmp), "block%i", b + 1);
        SET_STRING_ELT(blockCoefNamesS, b, magcomp::guarded_Rf_mkChar(tmp));
      }
      SET_VECTOR_ELT(blockCoefSamplesS, b, blockCoefArrayS[b]);
    }
    magcomp::guarded_Rf_setAttrib(blockCoefSamplesS, R_NamesSymbol, blockCoefNamesS);
  }
  piSamplesS = R_NilValue;
  if(fitComposition){
    piSamplesS = PROTECT(allocArrayReal3(nUnits, stickTotal, nSamples));
    nProtect++;
  }
  pSamplesS = R_NilValue;
  if(fitComposition){
    pSamplesS = PROTECT(allocArrayReal3(nUnits, outputJ, nSamples));
    nProtect++;
  }
  magnitudeSamplesS = R_NilValue;
  if(fitMagnitude){
    magnitudeSamplesS = PROTECT(allocMatrixReal(nUnits, nSamples));
    nProtect++;
  }
  nbSizeSamplesS = R_NilValue;
  nbSizeAcceptanceS = R_NilValue;
  nbSizeTuningS = R_NilValue;
  nbSizeBatchS = R_NilValue;
  nbSizeTuningHistoryS = R_NilValue;
  if(nb.estimate){
    nbSizeSamplesS = PROTECT(allocMatrixReal(1, nSamples));
    nProtect++;
  }

  magProcessSamplesS = R_NilValue;
  magProcessNamesS = R_NilValue;
  magProcessSigmaSamplesS = R_NilValue;
  magProcessSigmaNamesS = R_NilValue;
  magProcessPhiSamplesS = R_NilValue;
  magProcessPhiNamesS = R_NilValue;
  magProcessPhiAcceptanceS = R_NilValue;
  magProcessPhiTuningS = R_NilValue;
  magProcessPhiBatchS = R_NilValue;
  magProcessPhiTuningHistoryS = R_NilValue;
  magProcessArrayS.assign(mag.process.size(), R_NilValue);
  magProcessSigmaArrayS.assign(mag.process.size(), R_NilValue);
  magProcessPhiArrayS.assign(spatialProcessIndex.size(), R_NilValue);
  if(!mag.process.empty()){
    magProcessSamplesS = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP,
                                                mag.process.size()));
    magProcessNamesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP,
                                              mag.process.size()));
    magProcessSigmaSamplesS = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, mag.process.size()));
    magProcessSigmaNamesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, mag.process.size()));
    nProtect += 4;
    for(int mp = 0; mp < (int) mag.process.size(); mp++){
      int nRow = (int) mag.process[mp].effect.size();
      magProcessArrayS[mp] = PROTECT(allocMatrixReal(nRow, nSamples));
      magProcessSigmaArrayS[mp] = PROTECT(allocMatrixReal(1, nSamples));
      nProtect += 2;
      SET_STRING_ELT(magProcessNamesS, mp,
                     magcomp::guarded_Rf_mkChar(mag.process[mp].name.c_str()));
      SET_VECTOR_ELT(magProcessSamplesS, mp, magProcessArrayS[mp]);
      SET_STRING_ELT(magProcessSigmaNamesS, mp,
                     magcomp::guarded_Rf_mkChar(mag.process[mp].name.c_str()));
      SET_VECTOR_ELT(magProcessSigmaSamplesS, mp, magProcessSigmaArrayS[mp]);
    }
    magcomp::guarded_Rf_setAttrib(magProcessSamplesS, R_NamesSymbol, magProcessNamesS);
    magcomp::guarded_Rf_setAttrib(magProcessSigmaSamplesS, R_NamesSymbol, magProcessSigmaNamesS);
  }
  if(!spatialProcessIndex.empty()){
    magProcessPhiSamplesS =
      PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, spatialProcessIndex.size()));
    magProcessPhiNamesS =
      PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, spatialProcessIndex.size()));
    nProtect += 2;
    for(int sf = 0; sf < (int) spatialProcessIndex.size(); sf++){
      int mp = spatialProcessIndex[sf];
      magProcessPhiArrayS[sf] = PROTECT(allocMatrixReal(mag.process[mp].covMetrop.model->size, nSamples));
      nProtect++;
      SET_STRING_ELT(magProcessPhiNamesS, sf,
                     magcomp::guarded_Rf_mkChar(mag.process[mp].name.c_str()));
      SET_VECTOR_ELT(magProcessPhiSamplesS, sf,
                     magProcessPhiArrayS[sf]);
    }
    magcomp::guarded_Rf_setAttrib(magProcessPhiSamplesS, R_NamesSymbol,
                 magProcessPhiNamesS);
  }

  blockLoadingSamplesS = R_NilValue;
  blockFactorSamplesS = R_NilValue;
  blockFactorPhiSamplesS = R_NilValue;
  blockFactorPhiNamesS = R_NilValue;
  blockFactorPhiScalarNamesS = R_NilValue;
  blockFactorPhiAcceptanceS = R_NilValue;
  blockFactorPhiTuningS = R_NilValue;
  blockFactorPhiBatchS = R_NilValue;
  blockFactorPhiTuningHistoryS = R_NilValue;
  if(nFactorBlocks > 0){
    blockLoadingSamplesS = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, nFactorBlocks));
    blockFactorSamplesS = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, nFactorBlocks));
    blockLoadingNamesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, nFactorBlocks));
    blockFactorNamesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, nFactorBlocks));
    nProtect += 4;

    for(int fb = 0; fb < nFactorBlocks; fb++){
      int b = factorBlockIndex[fb];
      int f = factorTermIndex[fb];
      CompositionFactorState& factor = blocks[b].factor[f];
      blockLoadingArrayS.push_back(PROTECT(
        allocArrayReal3(blocks[b].coefStick, factor.K, nSamples)));
      blockFactorArrayS.push_back(PROTECT(
        allocArrayReal3(factor.nGroup, factor.K, nSamples)));
      nProtect += 2;
      std::string sampleName;
      if(!Rf_isNull(blockNamesS) && Rf_length(blockNamesS) == nBlocks){
        sampleName = CHAR(STRING_ELT(blockNamesS, b));
      }else{
        char tmp[32];
        std::snprintf(tmp, sizeof(tmp), "block%i", b + 1);
        sampleName = tmp;
      }
      if(blocks[b].factor.size() > 1)
        sampleName += ":" + factor.name;
      SET_STRING_ELT(blockLoadingNamesS, fb, magcomp::guarded_Rf_mkChar(sampleName.c_str()));
      SET_STRING_ELT(blockFactorNamesS, fb, magcomp::guarded_Rf_mkChar(sampleName.c_str()));
      SET_VECTOR_ELT(blockLoadingSamplesS, fb, blockLoadingArrayS[fb]);
      SET_VECTOR_ELT(blockFactorSamplesS, fb, blockFactorArrayS[fb]);
    }
    magcomp::guarded_Rf_setAttrib(blockLoadingSamplesS, R_NamesSymbol, blockLoadingNamesS);
    magcomp::guarded_Rf_setAttrib(blockFactorSamplesS, R_NamesSymbol, blockFactorNamesS);
  }
  if(nSpatialFactorBlocks > 0){
    blockFactorPhiSamplesS =
      PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, nSpatialFactorBlocks));
    blockFactorPhiNamesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP,
                                                  nSpatialFactorBlocks));
    nProtect += 2;
    for(int sf = 0; sf < nSpatialFactorBlocks; sf++){
      int b = spatialFactorBlockIndex[sf];
      int f = spatialFactorTermIndex[sf];
      CompositionFactorState& factor = blocks[b].factor[f];
      blockFactorPhiArrayS.push_back(PROTECT(
        allocMatrixReal(factor.K * factor.covMetrop[0].model->size, nSamples)));
      nProtect++;
      std::string sampleName;
      if(!Rf_isNull(blockNamesS) && Rf_length(blockNamesS) == nBlocks){
        sampleName = CHAR(STRING_ELT(blockNamesS, b));
      }else{
        char tmp[32];
        std::snprintf(tmp, sizeof(tmp), "block%i", b + 1);
        sampleName = tmp;
      }
      if(blocks[b].factor.size() > 1)
        sampleName += ":" + factor.name;
      SET_STRING_ELT(blockFactorPhiNamesS, sf, magcomp::guarded_Rf_mkChar(sampleName.c_str()));
      SET_VECTOR_ELT(blockFactorPhiSamplesS, sf, blockFactorPhiArrayS[sf]);
    }
    magcomp::guarded_Rf_setAttrib(blockFactorPhiSamplesS, R_NamesSymbol,
                 blockFactorPhiNamesS);
  }
  pMeanS = R_NilValue;
  if(fitComposition){
    pMeanS = PROTECT(allocMatrixReal(nUnits, outputJ));
    nProtect++;
  }
  densityMeanS = R_NilValue;
  if(fitMagnitude && fitComposition){
    densityMeanS = PROTECT(allocMatrixReal(nUnits, outputJ));
    nProtect++;
  }
  countMeanS = R_NilValue;
  if(fitMagnitude){
    countMeanS = PROTECT(allocMatrixReal(nUnits, outputJ));
    nProtect++;
  }
  totalDensityMeanS = R_NilValue;
  if(fitMagnitude){
    totalDensityMeanS = PROTECT(allocMatrixReal(nUnits, 1));
    nProtect++;
  }
  totalCountMeanS = R_NilValue;
  if(fitMagnitude){
    totalCountMeanS = PROTECT(allocMatrixReal(nUnits, 1));
    nProtect++;
  }

  if(verbose){
    magcomp::guarded_Rprintf("magComp sampler\n");
    magcomp::guarded_Rprintf("  components: %s%s%s\n", fitMagnitude ? "magnitude" : "",
      fitMagnitude && fitComposition ? " and " : "", fitComposition ? "composition" : "");
    magcomp::guarded_Rprintf("  sampling units: %i, output cells: %i\n", nUnits, outputJ);
    magcomp::guarded_Rprintf("  allocation blocks: %i\n", nBlocks);
    if(fitComposition) magcomp::guarded_Rprintf("  composition link: %s\n",
            compositionLinkName(compositionLink));
    if(fitMagnitude) magcomp::guarded_Rprintf("  magnitude coefficients: %i\n", mag.p);
    if(nb.estimate){
      magcomp::guarded_Rprintf("  magnitude NB size start: %.6f, prior Gamma(%.6f, %.6f), tuning: %.6f\n",
              nb.value, nb.shape, nb.rate, nb.tuning);
    }else if(fitMagnitude){
      magcomp::guarded_Rprintf("  magnitude NB size fixed: %.6f\n", nb.value);
    }
    for(int mp = 0; mp < (int) mag.process.size(); mp++){
      MagnitudeProcessState& proc = mag.process[mp];
      if(proc.isIid){
        magcomp::guarded_Rprintf("  magnitude %s iid%s groups: %i\n",
                proc.name.c_str(),
                proc.varying ? " varying-coefficient" : "",
                proc.nGroup);
        magcomp::guarded_Rprintf("  magnitude %s sigma_sq start: %.6f, prior IG(%.6f, %.6f)\n",
                proc.name.c_str(), proc.sigmaSq, proc.varianceShape, proc.varianceScale);
      }else if(proc.isSpatial){
        magcomp::guarded_Rprintf("  magnitude %s %s, covariance %s, sigma_sq %.6f\n",
          proc.name.c_str(), proc.denseGp ? "GP" : "NNGP", proc.covMetrop.model->name, proc.sigmaSq);
      }
    }
    for(int fb = 0; fb < nFactorBlocks; fb++){
      int b = factorBlockIndex[fb];
      int f = factorTermIndex[fb];
      CompositionFactorState& factor = blocks[b].factor[f];
      magcomp::guarded_Rprintf("  composition block %i %s %s%s: K %i, groups %i\n",
              b + 1, factor.name.c_str(), factor.termType.c_str(),
              factor.varying ? " varying-coefficient" : "",
              factor.K, factor.nGroup);
      if(factor.isSpatial){
        magcomp::guarded_Rprintf("    covariance %s\n", factor.covMetrop[0].model->name);
      }
    }
    if(fitComposition) magcomp::guarded_Rprintf("  composition smoothing: intercept %s, slope %s\n",
            getListString(smoothS, "intercept"), getListString(smoothS, "slope"));
    magcomp::guarded_Rprintf("  total iterations: %i, burn-in: %i, thin: %i, retained samples: %i\n",
            nIter, nBurn, nThin, nSamples);
    magcomp::guarded_R_FlushConsole();
  }

	  magcomp::rCall([]() { GetRNGstate(); });
  rngActive = true;
  profile.setupSeconds = profileNowSeconds() - profile.samplerStart;
  profileLoopStart = profileNowSeconds();
  gSamplerProfile = &profile;
	  saveIdx = 0;
	  for(int iter = 1; iter <= nIter; iter++){
    magcomp::rCall([]() { R_CheckUserInterrupt(); });
    profile.iterations++;
    if(fitMagnitude){
      profileStageStart = profileNowSeconds();
      updateMagnitudeNbSize(mag, nb, iter, nBurn);
      profile.nbSizeSeconds += profileNowSeconds() - profileStageStart;

      profileStageStart = profileNowSeconds();
      updateMagnitudeCoefficients(mag, nb.value, betaPriorPrec, Q.data(),
                                  rhs.data(), draw.data(), chol.data(),
                                  z.data(), rpg, rpgGamma, iter);
      profile.magnitudeCoefficientSeconds +=
        profileNowSeconds() - profileStageStart;
    }
    for(int mp = 0; mp < (int) mag.process.size(); mp++){
      profileStageStart = profileNowSeconds();
	      if(mag.process[mp].isIid){
	        updateMagnitudeIidEffects(mag, mp, nb.value,
                                  mag.process[mp].varianceShape, mag.process[mp].varianceScale);
	      }else if(mag.process[mp].isSpatial){
	        updateMagnitudeSpatialEffect(mag, mp, nb.value, mag.process[mp].varianceShape,
                                     mag.process[mp].varianceScale, Q.data(), rhs.data(),
                                     draw.data(), chol.data(), z.data(),
	                                     latentDiag.data(), nngpWorkC.data(),
	                                     nngpWorkRhs.data(),
	                                     nngpWorkCvec.data(), iter, nBurn);
	      }
      profile.magnitudeProcessSeconds +=
        profileNowSeconds() - profileStageStart;
	    }
	    for(int b = 0; b < nBlocks; b++){
      profileStageStart = profileNowSeconds();
	      updateAllocationBlock(blocks[b], Q.data(), rhs.data(), draw.data(),
	                            chol.data(), z.data(), rpg, rpgGamma,
                              compositionLink, iter);
      profile.allocationBlockSeconds +=
        profileNowSeconds() - profileStageStart;
	      for(int f = 0; f < (int) blocks[b].factor.size(); f++){
        profileStageStart = profileNowSeconds();
	        updateAllocationFactor(blocks[b], f, blocks[b].factor[f].loadingPriorPrec, Q.data(),
	                               rhs.data(), draw.data(), chol.data(),
	                               z.data(), latentDiag.data(),
	                               nngpWorkC.data(), nngpWorkRhs.data(),
	                               nngpWorkCvec.data(), iter, nBurn);
        profile.allocationFactorSeconds +=
          profileNowSeconds() - profileStageStart;
	      }
	    }

	    if(iter > nBurn && ((iter - nBurn) % nThin == 0)){
      profileStageStart = profileNowSeconds();
	      if(saveIdx >= nSamples)
	        magcomp::error("internal save index exceeded retained samples");

      if(nb.estimate)
        REAL(nbSizeSamplesS)[saveIdx] = nb.value;
      for(int l = 0; l < mag.p; l++)
        REAL(magnitudeCoefficientSamplesS)[l + saveIdx * mag.p] = mag.coef[l];
      for(int mp = 0; mp < (int) mag.process.size(); mp++){
        int nRow = (int) mag.process[mp].effect.size();
        REAL(magProcessSigmaArrayS[mp])[saveIdx] = mag.process[mp].sigmaSq;
        for(int r = 0; r < nRow; r++){
          REAL(magProcessArrayS[mp])[r + saveIdx * nRow] =
            mag.process[mp].effect[r];
        }
      }
      for(int sf = 0; sf < (int) spatialProcessIndex.size(); sf++){
        int mp = spatialProcessIndex[sf];
        const auto &cov = mag.process[mp].covMetrop;
        for(int p = 0; p < cov.model->size; ++p)
          REAL(magProcessPhiArrayS[sf])[p + saveIdx * cov.model->size] = cov.theta[p];
      }
      for(int fb = 0; fb < nFactorBlocks; fb++){
        int b = factorBlockIndex[fb];
        int f = factorTermIndex[fb];
        CompositionFactorState& factor = blocks[b].factor[f];
        for(int k = 0; k < factor.K; k++){
          for(int j = 0; j < blocks[b].coefStick; j++){
            REAL(blockLoadingArrayS[fb])[j + k * blocks[b].coefStick +
                                         saveIdx * blocks[b].coefStick *
                                         factor.K] =
              factor.loading[j + k * blocks[b].coefStick];
          }
          for(int g = 0; g < factor.nGroup; g++){
            REAL(blockFactorArrayS[fb])[g + k * factor.nGroup +
                                        saveIdx * factor.nGroup *
                                        factor.K] =
              factor.factorScore[g + k * factor.nGroup];
          }
        }
      }
      for(int sf = 0; sf < nSpatialFactorBlocks; sf++){
        int b = spatialFactorBlockIndex[sf];
        int f = spatialFactorTermIndex[sf];
        CompositionFactorState& factor = blocks[b].factor[f];
        for(int k = 0; k < factor.K; k++){
          const auto &cov = factor.covMetrop[k];
          for(int p = 0; p < cov.model->size; ++p)
            REAL(blockFactorPhiArrayS[sf])[k + factor.K*p + saveIdx*factor.K*cov.model->size] = cov.theta[p];
        }
      }

      for(int b = 0; b < nBlocks; b++){
        for(int j = 0; j < blocks[b].coefStick; j++){
          for(int l = 0; l < blocks[b].p; l++){
            REAL(blockCoefArrayS[b])[j + l * blocks[b].coefStick +
                                     saveIdx * blocks[b].coefStick *
                                     blocks[b].p] =
              blocks[b].coef[j + l * blocks[b].coefStick];
          }
        }
      }
      for(int b = 0; b < nBlocks; b++)
        fillBlockProportions(blocks[b], compositionLink);

      if(fitComposition) for(int i = 0; i < nUnits; i++){
        for(int j = 0; j < blocks[0].Jstick; j++){
          piValue = compositionInverseLink(
            blockLinearPredictor(blocks[0], i, j), compositionLink);
          REAL(piSamplesS)[i + j * nUnits +
                           saveIdx * nUnits * stickTotal] = piValue;
        }
        if(nBlocks == 2){
          stickCol = blocks[0].Jstick;
          for(int a = 0; a < blocks[0].J; a++){
            row = childRowByUnitParent[i + a * nUnits];
            for(int j = 0; j < blocks[1].Jstick; j++){
              piValue = compositionInverseLink(
                blockLinearPredictor(blocks[1], row, j), compositionLink);
              REAL(piSamplesS)[i + (stickCol + a * blocks[1].Jstick + j) *
                               nUnits + saveIdx * nUnits * stickTotal] =
                piValue;
            }
          }
        }
      }

      for(int i = 0; i < nUnits; i++){
        magnitude = 0.0;
        if(fitMagnitude){
          theta = rowDot(mag.X, mag.coef.data(), nUnits, mag.p, i);
          theta += magnitudeProcessEta(mag, i);
          magnitude = std::exp(theta);
          REAL(magnitudeSamplesS)[i + saveIdx * nUnits] = magnitude;
          REAL(totalDensityMeanS)[i] += magnitude / (double) nSamples;
          REAL(totalCountMeanS)[i] += mag.area[i] * magnitude / (double) nSamples;
          if(!fitComposition)
            REAL(countMeanS)[i] += mag.area[i] * magnitude / (double) nSamples;
        }
        if(!fitComposition)
          continue;
        for(int j = 0; j < outputJ; j++){
          idx = i + j * nUnits;
          if(nBlocks == 1){
            pValue = blocks[0].proportion[idx];
          }else{
            const int parent = j / childJ;
            const int child = j % childJ;
            row = childRowByUnitParent[i + parent * nUnits];
            pValue = blocks[0].proportion[i + parent * nUnits] *
              blocks[1].proportion[row + child * blocks[1].n];
          }
          REAL(pSamplesS)[idx + saveIdx * nUnits * outputJ] = pValue;
          REAL(pMeanS)[idx] += pValue / (double) nSamples;
          if(fitMagnitude){
            density = magnitude * pValue;
            count = mag.area[i] * density;
            REAL(densityMeanS)[idx] += density / (double) nSamples;
            REAL(countMeanS)[idx] += count / (double) nSamples;
          }
        }
      }

	      saveIdx++;
      profile.retained++;
      profile.saveSeconds += profileNowSeconds() - profileStageStart;
	    }

    if(verbose && nReport > 0 && iter % nReport == 0){
      magcomp::guarded_Rprintf("  iteration %i of %i\n", iter, nIter);
      magcomp::guarded_R_FlushConsole();
    }
	  }
  profile.loopSeconds = profileNowSeconds() - profileLoopStart;
  gSamplerProfile = NULL;
	  rngActive = false;
  magcomp::rCall([]() { PutRNGstate(); });

  if(saveIdx != nSamples)
    magcomp::error("internal save count mismatch");

  nRegisteredNngpGraphs = latentRegistry.nNngp;
  freeLatentBackendRegistry(latentRegistry);

  profileStageS = PROTECT(magcomp::guarded_Rf_allocVector(REALSXP, 9));
  profileStageNamesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, 9));
  nProtect += 2;
  SET_STRING_ELT(profileStageNamesS, 0, magcomp::guarded_Rf_mkChar("setup_seconds"));
  SET_STRING_ELT(profileStageNamesS, 1, magcomp::guarded_Rf_mkChar("loop_seconds"));
  SET_STRING_ELT(profileStageNamesS, 2, magcomp::guarded_Rf_mkChar("nb_size_seconds"));
  SET_STRING_ELT(profileStageNamesS, 3,
                 magcomp::guarded_Rf_mkChar("magnitude_coefficient_seconds"));
  SET_STRING_ELT(profileStageNamesS, 4,
                 magcomp::guarded_Rf_mkChar("magnitude_process_seconds"));
  SET_STRING_ELT(profileStageNamesS, 5,
                 magcomp::guarded_Rf_mkChar("allocation_block_seconds"));
  SET_STRING_ELT(profileStageNamesS, 6,
                 magcomp::guarded_Rf_mkChar("allocation_factor_seconds"));
  SET_STRING_ELT(profileStageNamesS, 7, magcomp::guarded_Rf_mkChar("save_seconds"));
  SET_STRING_ELT(profileStageNamesS, 8, magcomp::guarded_Rf_mkChar("iterations"));
  REAL(profileStageS)[0] = profile.setupSeconds;
  REAL(profileStageS)[1] = profile.loopSeconds;
  REAL(profileStageS)[2] = profile.nbSizeSeconds;
  REAL(profileStageS)[3] = profile.magnitudeCoefficientSeconds;
  REAL(profileStageS)[4] = profile.magnitudeProcessSeconds;
  REAL(profileStageS)[5] = profile.allocationBlockSeconds;
  REAL(profileStageS)[6] = profile.allocationFactorSeconds;
  REAL(profileStageS)[7] = profile.saveSeconds;
  REAL(profileStageS)[8] = (double) profile.iterations;
  magcomp::guarded_Rf_setAttrib(profileStageS, R_NamesSymbol, profileStageNamesS);

  long long profilePgTotalDraws = profile.pgMagnitudeDraws +
    profile.pgCompositionDraws + profile.pgAggregatedCompositionDraws +
    profile.pgOtherDraws;
  profilePgS = PROTECT(magcomp::guarded_Rf_allocVector(REALSXP, 16));
  profilePgNamesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, 16));
  nProtect += 2;
  SET_STRING_ELT(profilePgNamesS, 0, magcomp::guarded_Rf_mkChar("total_draws"));
  SET_STRING_ELT(profilePgNamesS, 1, magcomp::guarded_Rf_mkChar("magnitude_draws"));
  SET_STRING_ELT(profilePgNamesS, 2, magcomp::guarded_Rf_mkChar("composition_draws"));
  SET_STRING_ELT(profilePgNamesS, 3,
                 magcomp::guarded_Rf_mkChar("aggregated_composition_draws"));
  SET_STRING_ELT(profilePgNamesS, 4, magcomp::guarded_Rf_mkChar("other_draws"));
  SET_STRING_ELT(profilePgNamesS, 5, magcomp::guarded_Rf_mkChar("fallback_draws"));
  SET_STRING_ELT(profilePgNamesS, 6, magcomp::guarded_Rf_mkChar("magnitude_seconds"));
  SET_STRING_ELT(profilePgNamesS, 7, magcomp::guarded_Rf_mkChar("composition_seconds"));
  SET_STRING_ELT(profilePgNamesS, 8,
                 magcomp::guarded_Rf_mkChar("aggregated_composition_seconds"));
  SET_STRING_ELT(profilePgNamesS, 9, magcomp::guarded_Rf_mkChar("other_seconds"));
  SET_STRING_ELT(profilePgNamesS, 10, magcomp::guarded_Rf_mkChar("total_seconds"));
  SET_STRING_ELT(profilePgNamesS, 11, magcomp::guarded_Rf_mkChar("seconds_per_draw"));
  SET_STRING_ELT(profilePgNamesS, 12, magcomp::guarded_Rf_mkChar("fallback_fraction"));
  SET_STRING_ELT(profilePgNamesS, 13,
                 magcomp::guarded_Rf_mkChar("probit_composition_draws"));
  SET_STRING_ELT(profilePgNamesS, 14,
                 magcomp::guarded_Rf_mkChar("probit_composition_seconds"));
  SET_STRING_ELT(profilePgNamesS, 15,
                 magcomp::guarded_Rf_mkChar("probit_composition_seconds_per_draw"));
  REAL(profilePgS)[0] = (double) profilePgTotalDraws;
  REAL(profilePgS)[1] = (double) profile.pgMagnitudeDraws;
  REAL(profilePgS)[2] = (double) profile.pgCompositionDraws;
  REAL(profilePgS)[3] = (double) profile.pgAggregatedCompositionDraws;
  REAL(profilePgS)[4] = (double) profile.pgOtherDraws;
  REAL(profilePgS)[5] = (double) profile.pgFallbackDraws;
  REAL(profilePgS)[6] = profile.pgMagnitudeSeconds;
  REAL(profilePgS)[7] = profile.pgCompositionSeconds;
  REAL(profilePgS)[8] = profile.pgAggregatedCompositionSeconds;
  REAL(profilePgS)[9] = profile.pgOtherSeconds;
  REAL(profilePgS)[10] = profile.pgTotalSeconds;
  REAL(profilePgS)[11] = profilePgTotalDraws > 0 ?
    profile.pgTotalSeconds / (double) profilePgTotalDraws : NA_REAL;
  REAL(profilePgS)[12] = profilePgTotalDraws > 0 ?
    (double) profile.pgFallbackDraws / (double) profilePgTotalDraws : NA_REAL;
  REAL(profilePgS)[13] = (double) profile.probitCompositionDraws;
  REAL(profilePgS)[14] = profile.probitCompositionSeconds;
  REAL(profilePgS)[15] = profile.probitCompositionDraws > 0 ?
    profile.probitCompositionSeconds /
    (double) profile.probitCompositionDraws : NA_REAL;
  magcomp::guarded_Rf_setAttrib(profilePgS, R_NamesSymbol, profilePgNamesS);

	  outS = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, 34));
	  namesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP, 34));
	  nProtect += 2;

  if(nb.estimate){
    nbSizeAcceptanceS = PROTECT(magcomp::guarded_Rf_ScalarReal(
      nb.attempt > 0 ? (double) nb.accept / (double) nb.attempt : NA_REAL
    ));
    nbSizeTuningS = PROTECT(magcomp::guarded_Rf_ScalarReal(nb.tuning));
    nbSizeBatchS = PROTECT(scalarMetropolisHistorySexp(
      nb.metrop, nb.metrop.batchAcceptance
    ));
    nbSizeTuningHistoryS = PROTECT(scalarMetropolisHistorySexp(
      nb.metrop, nb.metrop.tuning
    ));
    nProtect += 4;
  }

  if(!spatialProcessIndex.empty()){
    magProcessPhiAcceptanceS =
      PROTECT(magcomp::guarded_Rf_allocVector(REALSXP, spatialProcessIndex.size()));
    magProcessPhiTuningS =
      PROTECT(magcomp::guarded_Rf_allocVector(REALSXP, spatialProcessIndex.size()));
    magProcessPhiBatchS =
      PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, spatialProcessIndex.size()));
    magProcessPhiTuningHistoryS =
      PROTECT(magcomp::guarded_Rf_allocVector(VECSXP, spatialProcessIndex.size()));
    nProtect += 4;
    for(int sf = 0; sf < (int) spatialProcessIndex.size(); sf++){
      MagnitudeProcessState& proc = mag.process[spatialProcessIndex[sf]];
      if(proc.estimateCovariance && proc.covAttempt > 0){
        REAL(magProcessPhiAcceptanceS)[sf] =
          (double) proc.covAccept / (double) proc.covAttempt;
      }else{
        REAL(magProcessPhiAcceptanceS)[sf] = NA_REAL;
      }
      REAL(magProcessPhiTuningS)[sf] =
        proc.estimateCovariance ? proc.covScale : NA_REAL;
      if(proc.estimateCovariance){
        SET_VECTOR_ELT(
          magProcessPhiBatchS, sf,
          scalarMetropolisHistorySexp(proc.covMetrop,
                                      proc.covMetrop.batchAcceptance)
        );
        SET_VECTOR_ELT(
          magProcessPhiTuningHistoryS, sf,
          scalarMetropolisHistorySexp(proc.covMetrop,
                                      proc.covMetrop.tuning)
        );
      }else{
        SET_VECTOR_ELT(magProcessPhiBatchS, sf,
                       magcomp::guarded_Rf_allocVector(REALSXP, 0));
        SET_VECTOR_ELT(magProcessPhiTuningHistoryS, sf,
                       magcomp::guarded_Rf_allocVector(REALSXP, 0));
      }
    }
    magcomp::guarded_Rf_setAttrib(magProcessPhiAcceptanceS, R_NamesSymbol,
                 magProcessPhiNamesS);
    magcomp::guarded_Rf_setAttrib(magProcessPhiTuningS, R_NamesSymbol,
                 magProcessPhiNamesS);
    magcomp::guarded_Rf_setAttrib(magProcessPhiBatchS, R_NamesSymbol,
                 magProcessPhiNamesS);
    magcomp::guarded_Rf_setAttrib(magProcessPhiTuningHistoryS, R_NamesSymbol,
                 magProcessPhiNamesS);
  }

  if(nSpatialFactorBlocks > 0){
    blockFactorPhiAcceptanceS = PROTECT(magcomp::guarded_Rf_allocVector(REALSXP,
                                                       nSpatialFactorPhi));
    blockFactorPhiTuningS = PROTECT(magcomp::guarded_Rf_allocVector(REALSXP,
                                                   nSpatialFactorPhi));
    blockFactorPhiBatchS = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP,
                                                  nSpatialFactorPhi));
    blockFactorPhiTuningHistoryS = PROTECT(magcomp::guarded_Rf_allocVector(VECSXP,
                                                          nSpatialFactorPhi));
    blockFactorPhiScalarNamesS = PROTECT(magcomp::guarded_Rf_allocVector(STRSXP,
                                                        nSpatialFactorPhi));
    nProtect += 5;
    int scalarPhiIndex = 0;
    for(int sf = 0; sf < nSpatialFactorBlocks; sf++){
      int b = spatialFactorBlockIndex[sf];
      int f = spatialFactorTermIndex[sf];
      CompositionFactorState& factor = blocks[b].factor[f];
      const char *sampleName = CHAR(STRING_ELT(blockFactorPhiNamesS, sf));
      for(int k = 0; k < factor.K; k++){
        char phiName[256];
        if(factor.K == 1){
          std::snprintf(phiName, sizeof(phiName), "%s", sampleName);
        }else{
          std::snprintf(phiName, sizeof(phiName), "%s:factor%i",
                        sampleName, k + 1);
        }
        SET_STRING_ELT(blockFactorPhiScalarNamesS, scalarPhiIndex,
                       magcomp::guarded_Rf_mkChar(phiName));

        if(!factor.covMetrop[k].free.empty() && factor.covAttempt[k] > 0){
          REAL(blockFactorPhiAcceptanceS)[scalarPhiIndex] =
            (double) factor.covAccept[k] / (double) factor.covAttempt[k];
        }else{
          REAL(blockFactorPhiAcceptanceS)[scalarPhiIndex] = NA_REAL;
        }
        REAL(blockFactorPhiTuningS)[scalarPhiIndex] =
          !factor.covMetrop[k].free.empty() ? factor.covScale[k] : NA_REAL;
        if(!factor.covMetrop[k].free.empty()){
          SET_VECTOR_ELT(
            blockFactorPhiBatchS, scalarPhiIndex,
            scalarMetropolisHistorySexp(factor.covMetrop[k],
                                        factor.covMetrop[k].batchAcceptance)
          );
          SET_VECTOR_ELT(
            blockFactorPhiTuningHistoryS, scalarPhiIndex,
            scalarMetropolisHistorySexp(factor.covMetrop[k],
                                        factor.covMetrop[k].tuning)
          );
        }else{
          SET_VECTOR_ELT(blockFactorPhiBatchS, scalarPhiIndex,
                         magcomp::guarded_Rf_allocVector(REALSXP, 0));
          SET_VECTOR_ELT(blockFactorPhiTuningHistoryS, scalarPhiIndex,
                         magcomp::guarded_Rf_allocVector(REALSXP, 0));
        }
        scalarPhiIndex++;
      }
    }
    magcomp::guarded_Rf_setAttrib(blockFactorPhiAcceptanceS, R_NamesSymbol,
                 blockFactorPhiScalarNamesS);
    magcomp::guarded_Rf_setAttrib(blockFactorPhiTuningS, R_NamesSymbol,
                 blockFactorPhiScalarNamesS);
    magcomp::guarded_Rf_setAttrib(blockFactorPhiBatchS, R_NamesSymbol,
                 blockFactorPhiScalarNamesS);
    magcomp::guarded_Rf_setAttrib(blockFactorPhiTuningHistoryS, R_NamesSymbol,
                 blockFactorPhiScalarNamesS);
  }

  SET_STRING_ELT(namesS, 0, magcomp::guarded_Rf_mkChar("magnitude_coefficient_samples"));
  SET_STRING_ELT(namesS, 1, magcomp::guarded_Rf_mkChar("composition_block_coefficient_samples"));
  SET_STRING_ELT(namesS, 2, magcomp::guarded_Rf_mkChar("pi_samples"));
  SET_STRING_ELT(namesS, 3, magcomp::guarded_Rf_mkChar("proportion_samples"));
  SET_STRING_ELT(namesS, 4, magcomp::guarded_Rf_mkChar("magnitude_samples"));
  SET_STRING_ELT(namesS, 5, magcomp::guarded_Rf_mkChar("proportion_fitted_mean"));
  SET_STRING_ELT(namesS, 6, magcomp::guarded_Rf_mkChar("density_fitted_mean"));
  SET_STRING_ELT(namesS, 7, magcomp::guarded_Rf_mkChar("count_fitted_mean"));
  SET_STRING_ELT(namesS, 8, magcomp::guarded_Rf_mkChar("total_density_fitted_mean"));
  SET_STRING_ELT(namesS, 9, magcomp::guarded_Rf_mkChar("total_count_fitted_mean"));
  SET_STRING_ELT(namesS, 10, magcomp::guarded_Rf_mkChar("nb_size"));
  SET_STRING_ELT(namesS, 11, magcomp::guarded_Rf_mkChar("nngp_graph_count"));
  SET_STRING_ELT(namesS, 12, magcomp::guarded_Rf_mkChar("composition_block_loading_samples"));
  SET_STRING_ELT(namesS, 13, magcomp::guarded_Rf_mkChar("composition_block_score_samples"));
  SET_STRING_ELT(namesS, 14, magcomp::guarded_Rf_mkChar("composition_block_covariance_samples"));
  SET_STRING_ELT(namesS, 15,
                 magcomp::guarded_Rf_mkChar("composition_block_covariance_acceptance"));
  SET_STRING_ELT(namesS, 16,
                 magcomp::guarded_Rf_mkChar("composition_block_covariance_tuning"));
  SET_STRING_ELT(namesS, 17,
                 magcomp::guarded_Rf_mkChar("composition_block_covariance_batch_acceptance"));
  SET_STRING_ELT(namesS, 18,
                 magcomp::guarded_Rf_mkChar("composition_block_covariance_tuning_history"));
  SET_STRING_ELT(namesS, 19, magcomp::guarded_Rf_mkChar("magnitude_process_effect_samples"));
  SET_STRING_ELT(namesS, 20,
                 magcomp::guarded_Rf_mkChar("magnitude_process_variance_samples"));
  SET_STRING_ELT(namesS, 21, magcomp::guarded_Rf_mkChar("magnitude_process_covariance_samples"));
  SET_STRING_ELT(namesS, 22,
                 magcomp::guarded_Rf_mkChar("magnitude_process_covariance_acceptance"));
  SET_STRING_ELT(namesS, 23, magcomp::guarded_Rf_mkChar("magnitude_process_covariance_tuning"));
  SET_STRING_ELT(namesS, 24,
                 magcomp::guarded_Rf_mkChar("magnitude_process_covariance_batch_acceptance"));
  SET_STRING_ELT(namesS, 25,
                 magcomp::guarded_Rf_mkChar("magnitude_process_covariance_tuning_history"));
  SET_STRING_ELT(namesS, 26, magcomp::guarded_Rf_mkChar("nb_size_samples"));
  SET_STRING_ELT(namesS, 27, magcomp::guarded_Rf_mkChar("nb_size_acceptance"));
  SET_STRING_ELT(namesS, 28, magcomp::guarded_Rf_mkChar("nb_size_tuning"));
	  SET_STRING_ELT(namesS, 29, magcomp::guarded_Rf_mkChar("nb_size_batch_acceptance"));
	  SET_STRING_ELT(namesS, 30, magcomp::guarded_Rf_mkChar("nb_size_tuning_history"));
  SET_STRING_ELT(namesS, 31, magcomp::guarded_Rf_mkChar("profile_stage_seconds"));
  SET_STRING_ELT(namesS, 32, magcomp::guarded_Rf_mkChar("profile_pg"));
  SET_STRING_ELT(namesS, 33, magcomp::guarded_Rf_mkChar("composition_link"));

  SET_VECTOR_ELT(outS, 0, magnitudeCoefficientSamplesS);
  SET_VECTOR_ELT(outS, 1, blockCoefSamplesS);
  SET_VECTOR_ELT(outS, 2, piSamplesS);
  SET_VECTOR_ELT(outS, 3, pSamplesS);
  SET_VECTOR_ELT(outS, 4, magnitudeSamplesS);
  SET_VECTOR_ELT(outS, 5, pMeanS);
  SET_VECTOR_ELT(outS, 6, densityMeanS);
  SET_VECTOR_ELT(outS, 7, countMeanS);
  SET_VECTOR_ELT(outS, 8, totalDensityMeanS);
  SET_VECTOR_ELT(outS, 9, totalCountMeanS);
  SET_VECTOR_ELT(outS, 10, fitMagnitude ? magcomp::guarded_Rf_ScalarReal(nb.value) : R_NilValue);
  SET_VECTOR_ELT(outS, 11, magcomp::guarded_Rf_ScalarInteger(nRegisteredNngpGraphs));
  SET_VECTOR_ELT(outS, 12, blockLoadingSamplesS);
  SET_VECTOR_ELT(outS, 13, blockFactorSamplesS);
  SET_VECTOR_ELT(outS, 14, blockFactorPhiSamplesS);
  SET_VECTOR_ELT(outS, 15, blockFactorPhiAcceptanceS);
  SET_VECTOR_ELT(outS, 16, blockFactorPhiTuningS);
  SET_VECTOR_ELT(outS, 17, blockFactorPhiBatchS);
  SET_VECTOR_ELT(outS, 18, blockFactorPhiTuningHistoryS);
  SET_VECTOR_ELT(outS, 19, magProcessSamplesS);
  SET_VECTOR_ELT(outS, 20, magProcessSigmaSamplesS);
  SET_VECTOR_ELT(outS, 21, magProcessPhiSamplesS);
  SET_VECTOR_ELT(outS, 22, magProcessPhiAcceptanceS);
  SET_VECTOR_ELT(outS, 23, magProcessPhiTuningS);
  SET_VECTOR_ELT(outS, 24, magProcessPhiBatchS);
  SET_VECTOR_ELT(outS, 25, magProcessPhiTuningHistoryS);
  SET_VECTOR_ELT(outS, 26, nbSizeSamplesS);
  SET_VECTOR_ELT(outS, 27, nbSizeAcceptanceS);
	  SET_VECTOR_ELT(outS, 28, nbSizeTuningS);
	  SET_VECTOR_ELT(outS, 29, nbSizeBatchS);
	  SET_VECTOR_ELT(outS, 30, nbSizeTuningHistoryS);
  SET_VECTOR_ELT(outS, 31, profileStageS);
  SET_VECTOR_ELT(outS, 32, profilePgS);
  SET_VECTOR_ELT(outS, 33,
                 fitComposition ? magcomp::guarded_Rf_mkString(compositionLinkName(compositionLink)) : R_NilValue);
	  magcomp::guarded_Rf_setAttrib(outS, R_NamesSymbol, namesS);

  UNPROTECT(nProtect);
  return outS;
}


extern "C" SEXP magcomp_sampler(SEXP config)
{
  struct Call { SEXP config; bool rngActive; } call{config, false};
  return R_UnwindProtect([](void *data) -> SEXP {
    Call &call = *static_cast<Call *>(data);
    return magcomp::nativeCall([&]() { return runSampler(call.config, call.rngActive); });
  }, &call, [](void *data, Rboolean) {
    Call &call = *static_cast<Call *>(data);
    gSamplerProfile = NULL;
    if(call.rngActive){
      call.rngActive = false;
      PutRNGstate();
    }
  }, &call, NULL);
}
