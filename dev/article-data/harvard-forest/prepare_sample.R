# Run from this directory after prepare_data.R.
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
source("read_inputs.R")
source("sample_functions.R")

set.seed(1)
inventory <- draw_inventory(forest, strata, census, prediction_cells, spacing, radius)

write_csv(inventory$plots, "data/plots.csv")
write_csv(inventory$trees |> select(-stratum), "data/trees.csv")
write_csv(inventory$prediction |> select(-stratum), "data/prediction.csv")
write_csv(truth, "data/truth.csv")
write_csv(inventory$area_summary, "data/area_summary.csv")
st_write(inventory$centers, "data/plot_centers.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(inventory$plot_units, "data/plot_units.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(inventory$prediction_units, "data/prediction_units.gpkg", delete_dsn = TRUE, quiet = TRUE)

# Every census stem belongs to exactly one observed or unsampled portion.
observed_members <- lengths(st_within(census, inventory$plot_units))
unsampled_members <- lengths(st_within(census, inventory$prediction_units))
stopifnot(all(observed_members + unsampled_members == 1))

plot_extents <- st_buffer(inventory$centers, dist = radius)
imagery <- rast("inputs/satellite.tif")
sample_map <- ggplot() +
  geom_spatraster_rgb(data = imagery, r = 1, g = 2, b = 3) +
  geom_sf(data = strata, aes(fill = factor(stratum_id)), color = "white", alpha = 0.2) +
  geom_sf(data = plot_extents, fill = NA, color = "darkorange", linewidth = 0.5) +
  geom_sf(data = inventory$centers, color = "black", size = 0.7) +
  scale_fill_brewer(palette = "Pastel1", name = "Stratum") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  coord_sf(datum = NA) +
  theme_void()
ggsave("checks/sample_map.png", sample_map, width = 8, height = 5, dpi = 180, bg = "white")
print(inventory$area_summary)
cat(nrow(inventory$centers), "plot centers;", nrow(inventory$plots), "sampled portions;",
    sum(inventory$trees$count), "sampled stems;", nrow(inventory$prediction), "prediction portions.\n")
