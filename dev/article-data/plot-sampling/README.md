# Shared plot-sampling data

These files support the plot-sampling articles. They are tracked in the
magComp GitHub repository under `dev/`, which is excluded from R package
builds. Quarto article sources live in `vignettes/articles/`, also excluded
from the built package.

## Reproduce the data

From this directory, run:

```sh
Rscript prepare_data.R
Rscript simulate_population.R
```

The scripts use tidyverse and sf. Preparation writes the forest, stand and
prediction-cell polygons, two predictor tables, and cell–stand overlap areas.
It also assigns stands to strata and saves separate plot and cell portions
for the stratum-effects example, as described below.
Simulation adds new total and species-by-DBH counts
to those tables. Each script overwrites its outputs; there are no dated copies.
Rerunning simulation alone replaces the counts using `set.seed(1)`.

To build the dataset article, run from the package root:

```sh
quarto render vignettes/articles/plot-sampling-data.qmd
```

Open `vignettes/articles/plot-sampling-data.html`, or use `quarto preview`
in place of `quarto render` for a preview that updates as the source changes.
Quarto, R, tidyverse, sf, knitr and rmarkdown are needed for rendering;
no package rebuild is required.

The article displays and runs chunks from `simulate_population.R`, keeping
the executable simulation and its explanation together. Render it after
running preparation. Its map reads the prepared cell polygons from `data/`.

## Article downloads

From `dev/article-data/`, run `Rscript build_downloads.R` after refreshing
the data or rebuilding the technical PDF. This writes the dataset ZIP archives and
the current technical PDF under `vignettes/articles/downloads/`. These are
website resources, excluded from R package builds.

`plot-sampling.zip` contains the seven shared files listed below. The data,
magnitude, composition, DBH, and species-by-DBH articles all link to that one archive.
`stratum-plot-sampling.zip` contains the prepared locations and polygons used
by the stratum article, which generates its effects and counts in memory.
The spatial, conditional DBH, and space-time examples have separate downloads,
each containing its generated population and optional generator. All archives
provide `data/`, matching the displayed relative paths. No fits or
prediction samples are included. Rebuilds overwrite the current downloads.

`_pkgdown.yml` puts the technical reference first on the article index and
groups the teaching articles in reading order, separating the shared-population
examples from the group, spatial, and space-time extensions. Download
resources are declared in each article's YAML so pkgdown copies them into
the site. From the package root, `Rscript dev/build_articles.R` refreshes
downloads and builds the article index and pages in `docs/articles/`, while
preserving the standalone HTML previews. Set `CI=true` for a local build
without pkgdown's optional online favicon generation.

## Inputs

The source geometry and covariates come from the Weisse Krug forest inventory
and UAV data used in the modeling paper. Data preparation uses
the supplied [source repository](https://github.com/Arkon74/masters-thesis)
at commit `baabf77af8f8df2b48d09c8d5f2f9e54b9da6728`.
All required inputs are copied here; running these scripts does not require
the manuscript repository. Source paths and checksums are in `inputs/provenance.csv`.

- `stands.gpkg`: the 249 original forest stand polygons in EPSG:25833.
- `plots.csv`: locations, forested plot areas, and p95 for all 151 measured
  plot locations. No observed tree counts are included.
- `prediction_cells.gpkg`: the 4,824 candidate forest-clipped square cells.
- `cell_predictors.csv`: cell coordinates, forested areas, p95, and the
  preparation's coverage flag. The script drops the 12 flagged cells.

p95 is the area-weighted empirical 95th percentile of the native 3 m
canopy-height raster within each plot or cell. The upstream preparation
calculates it separately over circular plots and square cells. Full prediction
cells have side length `sqrt(pi * 20^2)` m. Boundary cells and plots retain
their actual forested areas.

## Outputs and area convention

The shared files in `data/` are:

- `forest.gpkg`: one dissolved forest feature, preserving the original
  boundary, holes, and disconnected parts. It covers 513.066282 ha.
- `plot_locations.gpkg`: the 151 original plot centers, prepared as an sf
  point layer for mapping.
- `stands.gpkg`: the 249 original stand polygons, with `stand_id`,
  polygon `area_ha`, and `stratum_id` for the stratum-effects example.
- `plots.csv`: one row per actual plot location, identified by numeric `plot_id`.
- `prediction.csv`: one row per prediction cell, identified by numeric
  `plot_id`, using the numeric part of the source cell identifier.
- `prediction_cells.gpkg`: the 4,812 prepared prediction-cell polygons, with
  the same `plot_id` values as `prediction.csv`.
- `cell_stand_areas.csv`: one row per positive-area cell–stand intersection,
  with `prediction_id`, `stand_id`, and overlap `area_ha`. The prediction IDs
  match `plot_id` in `prediction.csv`.

Both location tables use `plot_id`, so the modeling articles read them without
renaming identifiers. IDs are unique within each table; observed plots and
prediction cells are separate sets of locations. Source files in `inputs/`
keep their original column names and identifiers.

Both tables have coordinates in meters, `area` in hectares and `p95` in meters.
Both also have `stand_id`, matching the stand polygons. This identifies
the stand containing each prediction location. For plots it identifies the
containing stand, or the nearest stand for the five centers just outside the forest
boundary (plots 49, 59, 104, 109, and 137). Their circular plots overlap the
forest. These labels do not change plot or cell areas or simulated counts.
Simulation adds the standardized predictor `x`, generating `expected_density`
in trees/ha, total `count`, and 15 species-by-DBH count columns. Species are
pine, oak, and birch; classes are 10–20, 20–30, 30–40, 40–50, and 50+ cm.
There are no latent spatial or stratum effects in this introductory population.

The generating magnitude model uses `beta_0 = log(350)`, `beta_1 = -0.2`
on standardized p95, and negative-binomial size 50. Expected densities range
from approximately 200 to 1,230 trees/ha across the prediction grid. The same
coefficients apply at plots and prediction cells; realized counts vary around
their negative-binomial means.

The prediction cells cover 512.770705 ha. Average cell densities with cell-area
weights, then multiply by the original forest polygon area if a forest total
is wanted. This extends the predicted average density to the uncovered
0.295577 ha (about 0.06% of the forest). It leaves density estimates unchanged
and treats that small area as having the same average density as the covered
forest. Simulation truth and model estimates must use this same convention.
Cell areas remain the exposures for count simulation and prediction.

All 151 plots are available for fitting; simulated prediction-cell counts
provide assessment observations without an additional holdout-plot dataset.

## Prepared inputs for the stratum example

`prepare_data.R` assigns the 249 stands to three strata of 83 stands each,
using `set.seed(1)`. It intersects plot circles and prediction cells with
the dissolved strata. This gives 206 observed portions from the 151 plots
and 6,575 prediction portions from the 4,812 cells. Forty-six plots and
1,636 cells cross stratum boundaries.

The additional files in `data/` are:

- `stratum_locations.csv`: the observed and prediction portions together,
  with `plot_id`, `source_id`, `stratum_id`, `area`, `p95`, standardized
  predictor `x`, and `dataset` (`plot` or `prediction`). It contains no
  simulated counts. Observed portions come first, followed by prediction
  portions, preserving the simulation order.
- `stratum_plot_units.gpkg`, `stratum_prediction_units.gpkg`: the matching
  portion polygons, with IDs, p95, and areas in hectares.

Each portion has one stratum and its own area. `plot_id` is unique within
each dataset; `source_id` identifies the original plot or prediction cell.
IDs are assigned after sorting by source ID and stratum. A multi-part
intersection within one stratum remains one unit. Each portion inherits its
parent plot or cell's p95 and uses the original 151 plots' standardization.
It is not assigned a newly calculated within-portion canopy-height percentile.

Prediction portions conserve the original cell areas and total covered
area. Reconstructed circular plot areas differ slightly from the original
tabulated areas because of the polygon approximation (at most 0.000058 ha).
The stratum example uses the calculated portion areas for both count
generation and fitting. In a field application, this model would require
counts by stratum portion for plots crossing boundaries.

The [stratum-effects article](../../../vignettes/articles/stratum-effects-plot-sampling.qmd)
reads these prepared inputs, generates its random effects and counts, and
uses that population in memory for fitting and comparison with truth.
The earlier modeling articles continue using `plots.csv` and `prediction.csv`.
Stratum formation, IDs, geometry, and the combination of observed and
prediction inputs are handled by `prepare_data.R`. The stratum article reads
the prepared locations and `plot_locations.gpkg`.

## Later stand summaries

A location's `stand_id` is useful for mapping and plot counts by stand.
For stand density or total estimates, use `cell_stand_areas.csv`: a prediction
cell can cross several stand boundaries, and its full area must not be assigned
to just the stand containing its location. Predict once for each cell, then
use the same density draw for every stand it overlaps. Multiplying that draw
by overlap `area_ha` gives the cell's contribution to that stand's total.
Sum these contributions within each draw and divide by the stand's covered
area to obtain stand density. The overlap areas sum to each prediction cell's
area and preserve the forest aggregation when summed across stands.

If a later example extends density estimates to each original stand polygon
area, scale the overlap weights within each stand to that polygon area. This
is the stand-level counterpart of the forest convention above; keep those
represented-area weights separate from the cell areas used for count
simulation and prediction. The data article shows how to construct the intersection table with
`st_intersection()` and `st_area()`. The magnitude article passes that table
to `standSummary()`, then uses `summary()` for stand means and intervals.
