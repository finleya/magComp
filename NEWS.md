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
