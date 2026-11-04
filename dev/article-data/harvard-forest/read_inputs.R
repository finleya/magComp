# Shared inputs for the illustrated inventory and repeated sampling.
# Run from the folder containing inputs/ and data/.
forest <- st_read("inputs/forest.gpkg", quiet = TRUE)
strata <- st_read("inputs/strata.gpkg", quiet = TRUE) |>
  arrange(stratum_id) |>
  st_set_agr("constant")
population <- read_csv("data/population.csv", show_col_types = FALSE)

species_order <- c("Hemlock", "Red maple", "Red oak", "White pine", "Other")
dbh_classes <- c("12.7-20", "20-30", "30-40", "40-50", "50-60", "60+")
census <- population |>
  mutate(
    species = case_when(
      species == "Tsuga canadensis" ~ "Hemlock", species == "Acer rubrum" ~ "Red maple",
      species == "Quercus rubra" ~ "Red oak", species == "Pinus strobus" ~ "White pine",
      .default = "Other"),
    species = factor(species, levels = species_order),
    dbh_class = cut(dbh_cm, breaks = c(12.7, 20, 30, 40, 50, 60, Inf),
                    labels = dbh_classes, right = FALSE, ordered_result = TRUE)
  ) |>
  st_as_sf(coords = c("easting", "northing"), crs = st_crs(forest), remove = FALSE)

radius <- 9
spacing <- sqrt(as.numeric(st_area(forest)) / 45)
stratum_areas <- read_csv("data/stratum_areas.csv", show_col_types = FALSE)

# These cells and stratum boundaries stay fixed across inventories.
prediction_cells <- st_sf(geometry = st_make_grid(forest, cellsize = sqrt(pi * radius^2))) |>
  mutate(cell_id = row_number()) |>
  st_set_agr("constant") |>
  st_intersection(strata |> select(stratum_id)) |>
  st_set_agr("constant")

truth <- census |>
  st_drop_geometry() |>
  count(stratum_id, species, dbh_class, name = "count") |>
  complete(stratum_id = 1:4, species, dbh_class, fill = list(count = 0L)) |>
  left_join(stratum_areas, by = "stratum_id") |>
  mutate(density = count / area_ha)
