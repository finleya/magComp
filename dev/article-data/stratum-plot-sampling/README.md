# Population with stratum effects

This population supports the [stratum-effects article](../../../vignettes/articles/stratum-effects-plot-sampling.qmd).
It reads the prepared stratum portions and standardized canopy-height predictor
from the shared `plot-sampling/data/` directory, then generates random effects and new counts.
It does not change the shared introductory population or the manuscript workflow.

## Run

Run `prepare_data.R` in `../plot-sampling/` to prepare the geographic inputs.
Then render the complete article from the package root:

```sh
quarto render vignettes/articles/stratum-effects-plot-sampling.qmd
```

The article displays and executes the simulation script's named chunks,
then fits the model and predicts. There is no fit or prediction cache.
The source and preview are tracked in Git but excluded from the built
R package. The simulated population stays in memory. To generate it in an
interactive R session, run `source("../stratum-plot-sampling/simulate_population.R")`
from `dev/article-data/plot-sampling/`. The displayed code uses `data/` paths,
so it also runs from a folder containing the extracted article download.
Prepared inputs and their upstream provenance are documented in
[the shared dataset README](../plot-sampling/README.md).

## Spatial construction

The shared `prepare_data.R` assigns stands to three strata, intersects
plots and cells with their boundaries, and combines their attributes in
`../plot-sampling/data/stratum_locations.csv`. Its `dataset` column identifies
the 206 observed portions and 6,575 prediction portions. The same preparation
saves `plot_locations.gpkg` for the map. This setup stays in development
code; neither article needs to construct the strata or combine input tables.

Each portion gets one `stratum_id`, a unique numeric `plot_id`, its polygon
area in hectares, and the original plot or cell's p95 and standardized `x`.
`source_id` links it to the original plot or prediction cell. The plot and
prediction ID sequences are separate. IDs are sorted by original unit and
stratum before being numbered. A multi-part intersection within one stratum
remains one unit.

Prediction portions conserve the original cell areas and total covered area,
512.770705 ha. Forest densities average with these
areas, as in the other articles. A forest total would extend that mean to
the original forest polygon area, following the shared dataset convention.

Counts are simulated separately on the portions. No spatial process is
included, and conditional count draws are independent given the parameters
and shared stratum effects. The same group effect is used for observed and
prediction portions in each stratum. For real boundary-crossing plots,
this construction would require counts by stratum portion.

## Generating model

- Magnitude: `beta_0 = log(350)`, `beta_1 = -0.2`, stratum-effect standard
  deviation 0.5, negative-binomial size 1,000.
- Species: the same pine and oak intercepts and slopes as the shared
  population, plus one unit-variance iid score per stratum. Pine's loading
  is 0.7 and oak's is -0.5. Birch receives the remainder.
- DBH: the same four sequential probit regressions within each species
  as the shared population; no stratum effects or smoothing.

Simulation starts with `set.seed(1)` before drawing the stratum effects,
independently of the seed used for geographic preparation. Effects are drawn
once per stratum. Their realized mean is not forced to zero. Their
generating distributions have mean zero, while the combined
intercepts record the particular effects realized in this population.

## Simulation objects

- `plots`, `prediction`: portion ID, original source ID, stratum,
  p95, standardized predictor, area, generating expected density, total
  count, and the 15 species-by-DBH counts.
- `stratum_effects`: the three magnitude effects and species scores.
- `true_intercepts`: the nine combined generating intercepts.

The polygons are shared inputs in `../plot-sampling/data/`. The article reuses
the forest prediction draws for stratum summaries and compares them with
the simulated stratum densities. It then intersects
`stratum_prediction_units.gpkg` with `stands.gpkg` to build an area table
using the portion IDs and compares stand estimates with area-weighted simulated
truth. The example table contains all 15 cells for one stand; the final figure
includes all 249 stands and all categories. The original population's
cell–stand table refers to whole cells and cannot be used directly with
these split-unit predictions.

Simulation objects, fits, and predictions are used in memory during article
rendering and are not saved as separate files.
