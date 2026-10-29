# Space-time plot-sampling example

`prepare_inputs.R` copies the shared forest geometry and standardized canopy
heights. Run it from this directory. `generate_population.R` then creates a
seed-1 population with 151 initial plots and 25 revisits in each of nine
additional years (376 measurements). Every visit has its own `plot_id`;
`site_id` identifies its physical location. Prediction contains the full
4,812-cell forest grid at years 0 and 9, with unique prediction measurement
IDs and `location_id` linking back to the cell polygons.

The two independent separable exponential fields are exact full-GP draws.
Separate spatial and temporal Cholesky factors avoid forming a full
site-by-year covariance matrix. Magnitude variance is 0.25, its spatial and
temporal decays are 3 and 0.35; the species factor has unit variance and decays
5 and 0.25. Spatial distances are kilometers and times are years. Count size
is 1000. Species loadings are 0.7 and -0.5; conditional DBH allocation uses
regression alone. The CSV files retain the true effects and parameters.

The article reads the supplied data, fits 15-neighbor NNGPs with three chains,
and predicts both dates jointly. Forest summaries group by year rather than
pooling dates. `build_downloads.R` in the parent directory refreshes the ZIP.
These development data and articles are excluded from the installed package.
