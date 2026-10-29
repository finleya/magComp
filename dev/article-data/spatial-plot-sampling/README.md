# Spatial plot-sampling example

This directory supplies the generated population for
[`spatial-effects-plot-sampling.qmd`](../../../vignettes/articles/spatial-effects-plot-sampling.qmd).
It reuses the forest, 151 plots, 4,812 prediction cells, canopy heights,
standardization, and cell–stand areas from `../plot-sampling/data/`.
The earlier articles' population is not modified.

## Generation and downloads

Run from this directory to reproduce the supplied CSV files:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript generate_population.R
```

The script reads the locations from the supplied CSV files, retaining only
spatial inputs and predictors before replacing the simulated columns. It
uses `set.seed(1)` and full exponential GP draws jointly over plots and
prediction cells. The independent magnitude field has variance 0.25 and
decay 3; the unit-variance species field has decay 5. Coordinates for these
calculations are in kilometers. Magnitude coefficients are log(350) and
-0.2; species loadings are 0.7 and -0.5. Species and DBH regressions follow
the shared population, and NB size is fixed at 1,000.

`plots.csv` and `prediction.csv` contain counts, true fields `v` and `w`,
combined intercepts, and expected density. `true_parameters.csv` provides
scalar generating values for the article's comparison table. The complete
DBH coefficients and fixed size are explicit in the generator. The supplied
polygons and overlap areas are unchanged copies of the shared inputs.

From `dev/article-data/`, run `Rscript build_downloads.R` to refresh the ZIP.
It includes `generate_population.R` alongside `data/`. Article rendering
reads these generated files; it does not regenerate the population. The
optional displayed generation chunks come directly from the script.
No fitted models or prediction arrays are included in the archive.

## Fit and preview

Install the current package, including `predict(keep_effects = TRUE)` and
`samples.magCompPrediction()`. From the package root, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 quarto render vignettes/articles/spatial-effects-plot-sampling.qmd
```

The article runs without cached fits. It fits 15-neighbor NNGPs with three
6,000-iteration chains, discards 2,000 iterations per chain, and uses 100
retained draws per chain for joint prediction. `n_cores = 3` assigns one
worker to each chain on systems supporting forked processes.
The generated full GP and fitted NNGP are deliberately distinguished.

The new-unit combined intercepts come from the same latent draws as the
posterior predictive densities. Parameter and observed-location combined
intercept diagnostics use all retained fitting draws; map summaries use
the smaller prediction subset. Forest and stand comparisons use realized
simulated cell counts and the same area weights as prediction summaries.
