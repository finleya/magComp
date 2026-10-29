# Conditional DBH effects example

This directory supplies the population for
[`conditional-dbh-effects-plot-sampling.qmd`](../../../vignettes/articles/conditional-dbh-effects-plot-sampling.qmd).
It reuses the shared forest geometry, 151 plots, 4,812 prediction cells,
canopy heights, standardization, and cell–stand intersection areas from
`../plot-sampling/data/`. The earlier populations are unchanged.

## Generate the population

Run from this directory to reproduce the supplied CSV files:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript generate_population.R
```

The script retains only spatial inputs and covariates from the supplied
location tables before generating new counts with seed 1. It draws one
unit-variance exponential GP jointly over plots and prediction cells, with
decay 3 and distances in kilometers. All 12 conditional DBH allocation rows
share this factor, with their own loadings. `dbh_parameters.csv` gives the
intercepts, slopes, and loadings. `pine:10-20` is the positive anchor.

Magnitude and leading species allocation have canopy-height regressions
and no random effects. Magnitude coefficients are log(350) and -0.2;
negative-binomial size is 50. The generating magnitude and species
coefficients and conditional factor decay are in `true_parameters.csv`.
The location tables include the true factor `w`, expected density, total
counts, and the 15 species-by-DBH counts.

## Fit and preview

From the package root, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 quarto render vignettes/articles/conditional-dbh-effects-plot-sampling.qmd
```

The article fits a 15-neighbor NNGP in `dbh_class | species`, with three
50,000-iteration chains and 20,000 burn-in iterations per chain, retaining
every tenth draw. Size is estimated because magnitude has no random effect.
Joint prediction uses
100 retained draws per chain over the whole grid. Combined-intercept maps
use those same stored latent draws. Parameter and observed-location
diagnostics use all 3,000 retained draws per chain.

The download includes this population and its optional generation script,
not fitted models or saved predictions. Rendering runs the analysis afresh.
Rebuild the archive with `Rscript build_downloads.R` from `dev/article-data/`.
