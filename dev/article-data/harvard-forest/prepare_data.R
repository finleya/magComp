# Run from this directory: Rscript prepare_data.R
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)

dir.create("checks", showWarnings = FALSE)

forest <- st_read("inputs/forest.gpkg", quiet = TRUE)
strata <- st_read("inputs/strata.gpkg", quiet = TRUE) |>
  arrange(stratum_id)

# Living stems at least 12.7 cm DBH define the census population.
population <- read_csv("inputs/census.csv") |>
  filter(status == "A", dbh_cm >= 12.7) |>
  st_as_sf(coords = c("easting", "northing"), crs = 32618,
           remove = FALSE) |>
  st_join(strata |> select(stratum_id)) |>
  st_drop_geometry()

stratum_areas <- strata |>
  mutate(area_ha = as.numeric(st_area(strata)) / 10000) |>
  st_drop_geometry() |>
  select(stratum_id, area_ha)

species_density <- population |>
  count(stratum_id, species, name = "stems") |>
  complete(stratum_id = 1:4, species, fill = list(stems = 0L)) |>
  left_join(stratum_areas, by = "stratum_id") |>
  mutate(stems_ha = stems / area_ha)

write_csv(population, "data/population.csv")
write_csv(stratum_areas, "data/stratum_areas.csv")
write_csv(species_density, "data/species_density.csv")

# Plot locations are generated in prepare_sample.R.
imagery <- rast("inputs/satellite.tif")
forest_map <- ggplot() +
  geom_spatraster_rgb(data = imagery, r = 1, g = 2, b = 3) +
  geom_sf(data = strata, aes(fill = factor(stratum_id)), color = "white", alpha = 0.2) +
  scale_fill_brewer(palette = "Pastel1", name = "Stratum") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  coord_sf(datum = NA) +
  theme_void()
ggsave("checks/forest_map.png", forest_map, width = 7, height = 5, dpi = 160, bg = "white")
