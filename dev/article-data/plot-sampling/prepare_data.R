# Prepare the shared plot-sampling example.
# Run from dev/article-data/plot-sampling: Rscript prepare_data.R

library(tidyverse)
library(sf)

# Keep the original stand boundaries and dissolve them into one forest.
# Coordinates are in meters (EPSG:25833); convert square meters to hectares.
stands <- st_read("inputs/stands.gpkg", quiet = TRUE)
forest <- st_sf(forest_id = 1, geometry = st_union(stands))
forest_area_ha <- as.numeric(st_area(forest)) / 10000

# The 151 inventory plots supply locations, sampled areas and canopy height.
# Tree counts will be simulated separately.
plots <- read_csv("inputs/plots.csv") |>
  transmute(
    plot_id,
    easting = easting_m,
    northing = northing_m,
    area = area_ha,
    p95 = chm_p95_m
  )

# Label each plot by its location. For centers outside the forest boundary,
# use the nearest stand; their circular plots still overlap the forest.
plot_locations <- st_as_sf(plots, coords = c("easting", "northing"),
                           crs = st_crs(stands))
plots$stand_id <- stands$stand_id[st_nearest_feature(plot_locations, stands)]
plots <- plots |> relocate(stand_id, .after = plot_id)

# Drop the 12 cells with insufficient UAV coverage.
# Full cells have the area of a 20-m radius plot; boundary cells are smaller.
# Use the numeric part of each source cell ID as plot_id in prepared files.
prediction <- read_csv("inputs/cell_predictors.csv") |>
  filter(retained) |>
  transmute(
    plot_id = parse_number(prediction_id),
    easting = easting_m,
    northing = northing_m,
    area = cell_area_ha,
    p95 = chm_p95_m
  )

# Each prediction location lies in one stand. This ID labels the location;
# the cell itself can span several stands.
prediction <- prediction |>
  st_as_sf(coords = c("easting", "northing"), crs = st_crs(stands),
           remove = FALSE) |>
  st_join(stands |> select(stand_id)) |>
  st_drop_geometry() |>
  relocate(stand_id, .after = plot_id)

# Save the area of each cell within each stand for later stand summaries.
# The same cell's predicted density applies to each of its overlapping parts.
cells <- st_read("inputs/prediction_cells.gpkg", quiet = TRUE) |>
  mutate(plot_id = parse_number(prediction_id)) |>
  semi_join(prediction, by = "plot_id") |>
  select(plot_id)

## ---- cell-stand-areas ----
# A cell can contribute to several stands, with a separate area for each.
st_agr(cells) <- "constant"
st_agr(stands) <- "constant"

cell_stand_areas <- st_intersection(
  cells |> select(prediction_id = plot_id),
  stands |> select(stand_id)
)

cell_stand_areas$area_ha <- as.numeric(st_area(cell_stand_areas)) / 10000

cell_stand_areas <- cell_stand_areas |>
  st_drop_geometry() |>
  filter(area_ha > 0) |>
  arrange(prediction_id, stand_id)

write_csv(cell_stand_areas, "data/cell_stand_areas.csv")

## ---- save-spatial-data ----
# Use cell areas to average densities. Multiply the average by forest_area_ha
# to estimate a forest total, extending that average to the small UAV gaps.
st_write(forest, "data/forest.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(plot_locations, "data/plot_locations.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(cells, "data/prediction_cells.gpkg", delete_dsn = TRUE, quiet = TRUE)
write_csv(plots, "data/plots.csv")
write_csv(prediction, "data/prediction.csv")

## ---- define-strata ----
# Group stands for the stratum-effects example.
set.seed(1)

stands <- stands |>
  mutate(stratum_id = sample(rep(1:3, length.out = n())))

strata <- stands |>
  group_by(stratum_id) |>
  summarise()

st_write(stands, "data/stands.gpkg", delete_dsn = TRUE, quiet = TRUE)

## ---- split-stratum-units ----
# Give each plot or cell portion one stratum and its own area.
plot_circles <- plots |>
  select(plot_id, easting, northing, p95) |>
  st_as_sf(coords = c("easting", "northing"), crs = st_crs(stands)) |>
  st_buffer(dist = 20)

stratum_cells <- cells |>
  select(plot_id) |>
  left_join(prediction |> select(plot_id, p95), by = "plot_id")

# IDs and canopy heights are carried into the intersected portions.
st_agr(plot_circles) <- "constant"
st_agr(stratum_cells) <- "constant"
st_agr(strata) <- "constant"

plot_units <- st_intersection(plot_circles, strata |> select(stratum_id)) |>
  rename(source_id = plot_id)
plot_units$area <- as.numeric(st_area(plot_units)) / 10000

plot_units <- plot_units |>
  filter(area > 0) |>
  arrange(source_id, stratum_id) |>
  mutate(plot_id = row_number())

prediction_units <- st_intersection(stratum_cells, strata |> select(stratum_id)) |>
  rename(source_id = plot_id)
prediction_units$area <- as.numeric(st_area(prediction_units)) / 10000

prediction_units <- prediction_units |>
  filter(area > 0) |>
  arrange(source_id, stratum_id) |>
  mutate(plot_id = row_number())

## ---- save-stratum-inputs ----
# Use the original plots' standardization for every portion.
p95_mean <- mean(plots$p95)
p95_sd <- sd(plots$p95)

stratum_plots <- plot_units |>
  st_drop_geometry() |>
  mutate(x = as.numeric(scale(p95, center = p95_mean, scale = p95_sd))) |>
  select(plot_id, source_id, stratum_id, area, p95, x)

stratum_prediction <- prediction_units |>
  st_drop_geometry() |>
  mutate(x = as.numeric(scale(p95, center = p95_mean, scale = p95_sd))) |>
  select(plot_id, source_id, stratum_id, area, p95, x)

stratum_locations <- bind_rows(
  stratum_plots |> mutate(dataset = "plot"),
  stratum_prediction |> mutate(dataset = "prediction")
)

write_csv(stratum_locations, "data/stratum_locations.csv")
st_write(plot_units, "data/stratum_plot_units.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(prediction_units, "data/stratum_prediction_units.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
