#ifndef MAGCOMP_COVARIANCE_METROPOLIS_H
#define MAGCOMP_COVARIANCE_METROPOLIS_H

// Included after ScalarMetropolisState: its small history interface is shared
// with the separate negative-binomial size update.
struct CovarianceState {
  const magcomp::CovModel *model = nullptr;
  std::vector<double> theta;
};

struct CovarianceMetropolis : ScalarMetropolisState, CovarianceState {
  std::vector<int> free;
  std::vector<double> lower, upper, shape, chol, states;
  int searchLimit = 0, searchBatches = 0, adaptiveBatches = 0;
  bool searching = false;

  void prepare(bool enabled, int burn, int iterations, int batch, double target)
  {
    adapt = enabled && !free.empty();
    batchLength = batch;
    batchPos = batchIndex = batchAccept = 0;
    targetAccept = target > 0 ? target : (free.size() == 1 ? 0.44 : 0.234);
    lastBatchAccept = NA_REAL;
    batchAcceptance.clear(); tuning.clear(); states.clear();
    batchAcceptance.reserve(iterations / batch);
    tuning.reserve(iterations / batch);
    searchLimit = std::min(20, burn / batch);
    searching = adapt && searchLimit > 0;
    searchBatches = adaptiveBatches = 0;
  }

  double transformed(int p) const
  { return std::log(theta[p] - lower[p]) - std::log(upper[p] - theta[p]); }

  bool propose(double scale, CovarianceState &proposal, double &jacobian)
  {
    proposal.model = model;
    proposal.theta = theta;
    const int q = free.size();
    std::vector<double> z(q);
    for(double &x : z) x = Rf_rnorm(0, 1);
    jacobian = 0;
    for(int i = 0; i < q; ++i){
      int p = free[i];
      double eta = transformed(p);
      for(int j = 0; j <= i; ++j) eta += scale * chol[i + j*q] * z[j];
      double x = lower[p] + (upper[p] - lower[p]) * Rf_plogis(eta, 0, 1, 1, 0);
      if(!(x > lower[p] && x < upper[p])) return false;
      jacobian += std::log(x - lower[p]) + std::log(upper[p] - x) -
        std::log(theta[p] - lower[p]) - std::log(upper[p] - theta[p]);
      proposal.theta[p] = x;
    }
    return true;
  }

  void update(int accepted, double &scale)
  {
    if(!adapt) return;
    const int q = free.size();
    for(int p : free) states.push_back(transformed(p));
    batchAccept += accepted;
    if(++batchPos < batchLength) return;
    double rate = double(batchAccept) / batchPos;
    lastBatchAccept = rate;
    if(searching){
      const double low = q == 1 ? 0.30 : 0.15;
      const double high = q == 1 ? 0.60 : 0.45;
      if(rate < 0.02) scale *= 0.1;
      else if(rate < 0.10) scale *= 0.25;
      else if(rate < low) scale *= 0.5;
      else if(rate > 0.65) scale *= 3;
      else if(rate > high) scale *= 1.5;
      if(++searchBatches >= searchLimit || (rate >= low && rate <= high)) searching = false;
    }else{
      const double gain = std::pow(double(++adaptiveBatches) + 1, -0.8);
      scale = std::exp(std::log(scale) + gain * (rate - targetAccept));
      if(batchPos > 1){
        std::vector<double> mean(q, 0);
        for(int t = 0; t < batchPos; ++t)
          for(int i = 0; i < q; ++i) mean[i] += states[t*q+i] / batchPos;
        for(int j = 0; j < q; ++j) for(int i = 0; i < q; ++i){
          double v = 0;
          for(int t = 0; t < batchPos; ++t)
            v += (states[t*q+i]-mean[i]) * (states[t*q+j]-mean[j]) / (batchPos-1);
          shape[i+j*q] += gain * (v - shape[i+j*q]);
        }
        // Convex updates retain positive definiteness of the initial proposal
        // shape. No diagonal noise is added to any model correlation matrix.
        chol = shape;
        int info = 0;
        F77_CALL(dpotrf)("L", &q, chol.data(), &q, &info FCONE);
        if(info) magcomp::error("adaptive covariance proposal lost positive definiteness");
      }
    }
    scale = std::max(1e-6, std::min(1e6, scale));
    ++batchIndex;
    batchAcceptance.push_back(rate); tuning.push_back(scale);
    batchPos = batchAccept = 0;
    states.clear();
  }
};

#endif
