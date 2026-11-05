# Harvard Forest inventory study

This example estimates forest and stratum stand tables from sampled plots in a mapped census. The population contains 20,258 living stems with DBH at least 12.7 cm. Five models add stratum and spatial effects in stages. Only sampled trees enter fitting; the census supplies truth for assessment.

## Start with one inventory

Use this directory as your working directory. Install magComp 0.0.7 or later, tidyverse, sf, terra, and tidyterra. The prepared first inventory contains 735 stems in 52 portions of 48 circular plots. The article reads `data/ladder/trees.csv` and `data/ladder/prediction.csv` and shows fitting, prediction, and aggregation. Their coordinates locate the centroids of the actual sampled and unsampled portions, in kilometers relative to the forest's southwest bounding-box corner.

To reproduce that inventory's five fits and summary tables from a terminal:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript repeat_ladder.R 1 1 1 3
Rscript export_ladder_example.R
```

The first four runner arguments are the first layout, last layout, number of concurrent layouts, and cores per fit. Each fit has three chains. `1 1 1 3` runs one inventory with its three chains in parallel. Set the fourth argument to 1 for sequential chains. An optional fifth argument sets a separate output directory. Parallel layouts and chains use forked processes on Linux and macOS; use one worker and one core on Windows.

The five models share the same counts and areas. H1 uses common intercepts. H2 adds a magnitude stratum random effect and a species stratum factor. H3 also adds a stratum factor in conditional DBH allocation. H4 adds a magnitude spatial effect; H5 also adds a species spatial factor. All composition process terms use one factor; spatial terms use NNGPs with 15 neighbors. No environmental covariates or DBH smoothing are used. Negative-binomial size is estimated in H1 and fixed at 1,000 in H2–H5.

## Cross-validation of the first inventory

Leave-one-plot-out cross-validation withholds each of the 48 physical plots in turn. Portions of a plot crossing stratum boundaries are withheld together. `magCompCV()` refits each model and generates posterior predictive counts for the omitted portions, using the same three-chain settings as above. Predictions use all 3,000 retained draws per chain.

```sh
for model in H1 H2 H3 H4 H5; do
  Rscript cross_validate_ladder.R "$model" 3
done
Rscript summarize_ladder_cv.R
```

The second argument specifies cores per fit. The commands above run models sequentially; they can also be run in separate terminals. Fits and prediction draws are saved in `output/cv/`. Summaries in `data/cv/` recombine portion counts within each draw and divide by the plot's sampled area. The table compares predicted and observed species-by-DBH densities, with equal weight for each of the 48 × 30 plot cells. `diagnostics.csv` records the parameter with the largest R-hat in each fit. This assesses prediction of withheld observations from the first inventory; the repeated-sampling experiment below assesses forest and stratum estimates against the census.

`portion_predictions.csv` retains the held-out density summaries for each stratum portion, using its own sampled area. `stratum_summary.csv` reports mean error, RMSPE, interval coverage, and interval width by model and stratum, with equal weight for each portion's 30 category cells. The article displays RMSPE and coverage alongside the number of sampled portions. These are prediction scores within strata, not estimates of stratum totals.

## Repeat the sampling

The census remains fixed while a random shift changes the systematic plot grid. The included results cover 1,000 layouts. To reproduce the experiment with four layouts in parallel and three cores per fit:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript repeat_ladder.R 1 1000 4 3
Rscript summarize_ladder_repetitions.R 1 1000
Rscript export_ladder_example.R
```

Use a smaller last-layout number to try a shorter experiment. `ladder_settings.R` uses three 50,000-iteration chains with 20,000 burn-in and thin 10. Predictions use 333 draws per chain. All 5,000 fits used these settings without extensions. Every fit uses chain seeds `1:3`; grid placement and prediction use the layout number as their seed. Every model uses the same inventory for a given layout.

Fits and predictions are saved in `output/repeated_ladder/`; the article reads compact CSV summaries in `data/ladder` and `data/repeated_ladder`. Completed fits are reused when a run is resumed. For a separate run, the runner accepts an output directory as its fifth argument; the summary script accepts fit and summary directories as its third and fourth arguments. Saving all fits for 1,000 layouts requires substantial disk space.

`data/repeated_ladder/run_counts.csv` records the fitting status. A status of `unresolved` indicates a diagnostic flag, not an execution failure. Sixteen fits exceeded the R-hat screening value of 1.05; their diagnostics remain in the download. Performance comparisons include all 1,000 layouts for all five models. `attempts.csv` records the fitting settings and diagnostics. Per-cell bias, relative bias, RMSE, interval coverage, and width are in `target_summary.csv`; `domain_summary.csv` averages cells within each area. Monte Carlo standard errors use inventories as the independent repetitions.

## Data and calculation

| File | Contents |
|---|---|
| `inputs/census.csv` | Original processed census records, including stem identifiers, coordinates, species, DBH, biomass, and status. |
| `inputs/forest.gpkg`, `inputs/strata.gpkg` | The 35-ha forest and four strata. |
| `inputs/satellite.tif` | The book's georeferenced satellite image, used as the map background. |
| `inputs/canopy_height.tif` | Source canopy-height raster; not used in these models. |
| `data/population.csv` | Eligible living stems assigned to strata. |
| `data/stratum_areas.csv` | Full polygon areas in hectares. |
| `data/plot_centers.gpkg`, `data/plot_units.gpkg` | Original plot centers and sampled stratum portions for layout 1. |
| `data/prediction_units.gpkg` | Unsampled prediction portions for layout 1. |
| `data/ladder/trees.csv`, `data/ladder/prediction.csv` | Model inputs with actual portion areas and spatial coordinates. |
| `data/area_summary.csv`, `data/truth.csv` | Sampled and unsampled areas, and census species-by-DBH counts and densities. |
| `data/ladder/` | First-inventory estimates, chain diagnostics, and settings. |
| `data/cv/` | Plot exclusions, held-out density summaries, predictive scores, and parameter diagnostics. |
| `data/repeated_ladder/` | Repeated-sampling performance and run accounting. |

Species groups are hemlock, red maple, red oak, white pine, and Other. DBH classes are [12.7,20), [20,30), [30,40), [40,50), [50,60), and [60,∞) cm. Each sampled portion has all 30 count cells, including zeros. `plot_id` identifies a circular plot; `unit_id` identifies its stratum portion.

Plots have radius 9 m, on a grid with spacing √(350,000/45) m and a uniform random shift. Centers may fall in the 9-m band outside the forest. Counts include only eligible stems inside the forest. Plots crossing strata are split and their actual areas recorded. Full prediction cells have area π9² m²; clipping to strata and removing sampled areas defines the unsampled remainder. Areas and counts are checked during preparation.

Each population-count draw adds the known sampled count to a joint posterior prediction over the unsampled remainder. Dividing by full stratum area gives stratum density. Forest density sums the stratum counts within each draw and divides by forest area. The calculation never replaces observed counts with simulated replicates.

| Script | Purpose |
|---|---|
| `prepare_data.R`, `prepare_sample.R` | Prepare the census and first plot sample. |
| `read_inputs.R`, `sample_functions.R` | Read spatial inputs, place each grid, count trees, and construct unsampled portions. |
| `ladder_settings.R` | MCMC lengths, thinning, prediction draws, and convergence checks. |
| `ladder_functions.R` | Five model formulas, fitting, and diagnostics. |
| `cross_validate_ladder.R`, `summarize_ladder_cv.R` | Leave-one-plot-out fitting and predictive assessment for the first inventory. |
| `population_summaries.R` | Joint prediction and observed-plus-predicted count aggregation. |
| `repeat_ladder.R` | Run layouts and save fits, attempts, and results. |
| `summarize_ladder_repetitions.R` | Calculate repeated-sampling metrics and Monte Carlo standard errors. |
| `export_ladder_example.R` | Export compact first-inventory results for the article. |

To regenerate the prepared census and map inputs, run `source("prepare_data.R")` and `source("prepare_sample.R")` in R. Coordinates in the GeoPackages are meters in WGS 84 / UTM zone 18N (EPSG:32618). Areas are hectares. Biomass is retained from the source but is not modeled here.

## Sources

Finley, A. O. and Doser, J. W. (2026). *Introduction to Forestry Data Analysis with R*. Chapman & Hall/CRC. [Book](https://ifdar-book.org/), [DOI](https://doi.org/10.1201/9781003222699). The buffered-frame approach is illustrated in [Section 17.7.2](https://www.finley-lab.com/ifdar/cluster#clusterHF).

Orwig, D., Foster, D., and Ellison, A. (2023). *Harvard Forest CTFS-ForestGEO Mapped Forest Plot since 2014*. Harvard Forest Data Archive: HF253 (v.6). Environmental Data Initiative. [DOI: 10.6073/pasta/818789a882a318c1d7f3fc43a2289e12](https://doi.org/10.6073/pasta/818789a882a318c1d7f3fc43a2289e12). This is the citation supplied by the [data archive](https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF253); the inputs here are the book's processed first-census data. We thank David Orwig, David Foster, and Aaron Ellison for making these data available. The archive releases its data under CC0 and requests acknowledgment of the creators.

The census, boundaries, and canopy-height raster came from the reader archive `IFDAR_book.zip` at book-repository commit `15cb88419a79855c032e0a78109a8a4ac1188076`. Archive SHA-256: `dc569b72b742adfdfba3c591e3cc3a0d1db40cabc791d30c0934a04a7fe4a1bf`. The satellite image is the book's `datasets/HF/HF_GES.tif`. `sources.csv` records the source paths and checksums.
