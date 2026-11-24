# magComp 0.0.8

* Allow separate observations to share GP or NNGP effects at identical inputs.
  Each observation retains its own count, sampled area, and covariates. Shared
  effects work in magnitude and composition models, including conditional
  allocation, varying slopes, and space-time terms.
* Support one-dimensional and multivariate covariate inputs, such as
  `gp(height)`, `gp(height, cover)`, and `x:nngp(height)`, for nonlinear and
  varying-coefficient relationships. Input scaling remains under user control.
* Map observations to distinct process inputs throughout fitting, starting
  values, prediction, combined-coefficient summaries, and cross-validation.
  Existing models with distinct inputs retain their ordering and calculations.
* Allow an explicit Gamma prior to request negative-binomial size estimation
  alongside magnitude random effects. Automatic defaults remain unchanged:
  estimate size without magnitude random effects and fix it otherwise.

* Support data and cross-validation training folds with all-zero counts.
  Magnitude uses the zero-count likelihood; composition with no allocation
  information is governed by its priors.

* Reject fractional chain counts, worker counts, and Metropolis batch lengths
  instead of silently truncating them.

* Add an article on nonlinear covariate effects and varying coefficients.
  Update the space-time article to distinguish plot and measurement IDs and
  explain shared effects, and clarify size-estimation options in the
  random-effect articles.

# magComp 0.0.7

* Improve exploration of composition factor signs by sampling unrestricted
  working loadings and jointly aligning loading columns and factor scores
  to positive anchors. This preserves the existing model, half-normal
  diagonal priors, and upper-anchor zeros for group, spatial, and space-time
  factors in leading and conditional allocation blocks.

* Add a Harvard Forest inventory article with a mapped census, five model
  specifications, and a 1,000-layout repeated-sampling assessment. The data,
  analysis code, and compact results are available with the article.

# magComp 0.0.6

Initial public release.

* Fit magnitude and composition models together or separately, using a
  negative-binomial model for total counts and multinomial allocation among
  one or two categorical variables.
* Specify separate covariates for magnitude, leading allocation, and
  conditional allocation. Include group, spatial, or space-time effects
  in intercepts and slopes.
* Represent dependence among composition categories with shared latent
  factors and smooth coefficients across ordered classes.
* Use GP or NNGP terms with exponential, Matérn, separable space-time,
  multiresolution, or Gneiting covariance functions. Covariance parameters
  have automatic defaults and can be estimated or fixed individually.
* Inspect posterior draws and chain diagnostics through `samples()`,
  `summary()`, and `plot()`, including combined coefficients and effects.
* Obtain fitted values and predictions for expected quantities or simulated
  observations. Aggregate density draws to reporting areas with
  `standSummary()`, and assess models with posterior predictive checks and
  cross-validation.
* Follow worked Quarto articles with downloadable data, from individual
  components to models with group, spatial, and space-time effects.
