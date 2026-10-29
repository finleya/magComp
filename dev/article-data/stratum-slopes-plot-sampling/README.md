# Stratum-specific intercepts and slopes

This population supports
[`stratum-slopes-plot-sampling.qmd`](../../../vignettes/articles/stratum-slopes-plot-sampling.qmd).
It reads the prepared portions and standardized canopy height from
`../plot-sampling/data/stratum_locations.csv`, just as the stratum-intercept
article does. It generates a population in memory using those prepared inputs.

## Run

From the package root, render the article:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 quarto render vignettes/articles/stratum-slopes-plot-sampling.qmd
```

The article displays and runs the generator's named chunks. In an interactive
session, run `source("../stratum-slopes-plot-sampling/simulate_population.R")`
from `dev/article-data/plot-sampling/`. Displayed code reads `data/...`, so
it also runs from the folder containing the extracted `stratum-plot-sampling.zip`.
No new geographic dataset or download is needed. There are no cached fits,
saved predictions, or saved simulated population files.

## Model

All four sets of stratum effects are independent draws after `set.seed(1)`.
Magnitude uses coefficients log(350) and -0.2, intercept standard deviation
0.5, slope standard deviation 0.2, and fixed negative-binomial size 1,000.
Species regression coefficients and DBH regressions follow the shared
population. Species intercept loadings are 0.7 for pine and -0.5 for oak;
slope loadings are 0.4 and -0.3, respectively. Intercept and slope scores
have separate unit-variance Normal distributions. Pine anchors each term.
DBH allocation has no random effects or smoothing.

The same realized stratum effects apply to fitted and prediction units.
The prepared geography has 206 observed portions from 151 original plots
and 6,575 prediction portions, with areas conserved as described in the
[stratum-intercept dataset notes](../stratum-plot-sampling/README.md).

`true_coefficients` records 18 combined coefficients: an intercept and slope
for magnitude, pine, and oak, in each of three strata. These are compared
with `summary(fit, type = "combined", group = "stratum_id")`. Separate
intercept and slope terms do not estimate an intercept–slope covariance.

## Fitting and predictions

The fit uses default priors and starts, three 4,000-iteration chains,
1,000 burn-in iterations per chain, and three fitting workers. Parameter
and combined-coefficient diagnostics use all 3,000 retained draws per chain.
Predictions use 100 draws per chain, with one BLAS/OpenMP thread per worker
in the development checks.

Relationship figures evaluate 50 canopy-height values within each stratum's
observed range. They use the expected target, including stratum effects.
Species probabilities sum the five DBH probabilities within each draw before
summarizing. Forest densities instead use observation-target counts and
allocations, divided by unit area and then area-weighted with `standSummary()`.
