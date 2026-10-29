# Prepare the existing forest geometry and canopy heights for the time example.
# Run from dev/article-data/space-time-plot-sampling/.
library(tidyverse)

plots <- read_csv("../plot-sampling/data/plots.csv") |>
  select(plot_id, stand_id, easting, northing, area, p95, x) |>
  rename(location_id = plot_id) |>
  mutate(dataset = "plot")

prediction <- read_csv("../plot-sampling/data/prediction.csv") |>
  select(plot_id, stand_id, easting, northing, area, p95, x) |>
  rename(location_id = plot_id) |>
  mutate(dataset = "prediction")

locations <- bind_rows(plots, prediction) |>
  mutate(site_id = row_number(), easting_km = easting / 1000,
         northing_km = northing / 1000)

write_csv(locations, "data/locations.csv")
file.copy("../plot-sampling/data/forest.gpkg", "data/forest.gpkg", overwrite = TRUE)
file.copy("../plot-sampling/data/prediction_cells.gpkg", "data/prediction_cells.gpkg", overwrite = TRUE)
