# A random shift over one grid cell gives equal sampling intensity across
# the forest and its outside buffer. Call set.seed() before this function.
systematic_centers <- function(forest, spacing, radius) {
  extent <- st_bbox(forest)
  shift <- runif(2, 0, spacing)
  expand_grid(
    easting = seq(extent["xmin"] - radius + shift[1], extent["xmax"] + radius, by = spacing),
    northing = seq(extent["ymin"] - radius + shift[2], extent["ymax"] + radius, by = spacing)
  ) |>
    st_as_sf(coords = c("easting", "northing"), crs = st_crs(forest), remove = FALSE) |>
    st_filter(st_buffer(forest, radius, nQuadSegs = 180)) |>
    mutate(plot_id = row_number())
}

# Count sampled stems and define unsampled prediction portions.
draw_inventory <- function(forest, strata, census, prediction_cells, spacing, radius) {
  centers <- systematic_centers(forest, spacing, radius)
  plot_units <- st_intersection(
    st_buffer(centers |> select(plot_id) |> st_set_agr("constant"), radius, nQuadSegs = 180),
    strata |> select(stratum_id)
  ) |>
    filter(as.numeric(st_area(geometry)) > 0) |>
    arrange(plot_id, stratum_id) |>
    mutate(unit_id = row_number(), area_ha = as.numeric(st_area(geometry)) / 10000)

  sampled <- st_join(census, plot_units |> select(unit_id, plot_id),
                     join = st_within, left = FALSE)
  plots <- plot_units |>
    st_drop_geometry() |>
    left_join(centers |> st_drop_geometry(), by = "plot_id") |>
    select(unit_id, plot_id, stratum_id, area_ha, easting, northing)

  # Complete all species-by-DBH cells, including empty sampled portions.
  trees <- sampled |>
    st_drop_geometry() |>
    count(unit_id, species, dbh_class, name = "count") |>
    complete(unit_id = plots$unit_id, species, dbh_class, fill = list(count = 0L)) |>
    left_join(plots, by = "unit_id") |>
    mutate(stratum = factor(stratum_id, levels = strata$stratum_id)) |>
    arrange(unit_id, species, dbh_class)

  prediction_units <- st_difference(prediction_cells, st_union(plot_units)) |>
    filter(as.numeric(st_area(geometry)) > 0) |>
    arrange(cell_id, stratum_id) |>
    mutate(unit_id = row_number(), area_ha = as.numeric(st_area(geometry)) / 10000)
  prediction <- prediction_units |>
    st_drop_geometry() |>
    select(unit_id, cell_id, stratum_id, area_ha) |>
    mutate(stratum = factor(stratum_id, levels = strata$stratum_id))

  # Include a row for any stratum with no sampled portion.
  area_summary <- strata |>
    mutate(area_ha = as.numeric(st_area(strata)) / 10000) |>
    st_drop_geometry() |>
    select(stratum_id, area_ha) |>
    left_join(plots |>
      group_by(stratum_id) |>
      summarise(sampled_area_ha = sum(area_ha), plot_portions = n(), .groups = "drop"),
      by = "stratum_id") |>
    mutate(sampled_area_ha = replace_na(sampled_area_ha, 0),
           plot_portions = replace_na(plot_portions, 0L)) |>
    left_join(prediction |>
      group_by(stratum_id) |>
      summarise(unsampled_area_ha = sum(area_ha), .groups = "drop"), by = "stratum_id") |>
    mutate(sampled_percent = 100 * sampled_area_ha / area_ha,
           no_sampled_portion = plot_portions == 0)

  # Check the sampling geometry, without looking at model estimates.
  stopifnot(
    max(abs(area_summary$sampled_area_ha + area_summary$unsampled_area_ha -
              area_summary$area_ha)) < 1e-7,
    n_distinct(paste(sampled$tree_id, sampled$stem_id)) == nrow(sampled),
    sum(trees$count) == nrow(sampled)
  )
  list(centers = centers, plot_units = plot_units, plots = plots, trees = trees,
       prediction_units = prediction_units, prediction = prediction,
       area_summary = area_summary)
}
