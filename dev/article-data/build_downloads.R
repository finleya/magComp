# Refresh article downloads after preparing the data or rebuilding magComp.pdf.
# Run from dev/article-data/: Rscript build_downloads.R

dir.create("../../vignettes/articles/downloads", showWarnings = FALSE)
downloads <- normalizePath("../../vignettes/articles/downloads")

# One shared archive for the introductory plot-sampling articles.
# Each archive extracts to data/, matching the paths shown in the articles.
zip::zipr(
  file.path(downloads, "plot-sampling.zip"),
  files = c("data/forest.gpkg", "data/plot_locations.gpkg", "data/stands.gpkg",
            "data/plots.csv", "data/prediction.csv", "data/prediction_cells.gpkg",
            "data/cell_stand_areas.csv"),
  root = "plot-sampling",
  mode = "mirror"
)

# The stratum article generates a new population from these prepared inputs.
zip::zipr(
  file.path(downloads, "stratum-plot-sampling.zip"),
  files = c("data/forest.gpkg", "data/plot_locations.gpkg", "data/stands.gpkg",
            "data/stratum_locations.csv", "data/stratum_prediction_units.gpkg"),
  root = "plot-sampling",
  mode = "mirror"
)

# The spatial example ships its full-GP population and optional generator.
zip::zipr(
  file.path(downloads, "spatial-plot-sampling.zip"),
  files = c("generate_population.R", "data/forest.gpkg", "data/plot_locations.gpkg",
            "data/stands.gpkg", "data/plots.csv", "data/prediction.csv",
            "data/true_parameters.csv", "data/prediction_cells.gpkg",
            "data/cell_stand_areas.csv"),
  root = "spatial-plot-sampling",
  mode = "mirror"
)

# The conditional DBH example has its own population and loading rows.
zip::zipr(
  file.path(downloads, "conditional-dbh-plot-sampling.zip"),
  files = c("generate_population.R", "data/forest.gpkg", "data/plot_locations.gpkg",
            "data/stands.gpkg", "data/plots.csv", "data/prediction.csv",
            "data/true_parameters.csv", "data/dbh_parameters.csv",
            "data/prediction_cells.gpkg", "data/cell_stand_areas.csv"),
  root = "conditional-dbh-plot-sampling",
  mode = "mirror"
)

# Repeated visits and joint predictions at two forest dates.
zip::zipr(
  file.path(downloads, "space-time-plot-sampling.zip"),
  files = c("generate_population.R", "data/locations.csv", "data/forest.gpkg",
            "data/plots.csv", "data/prediction.csv", "data/true_parameters.csv",
            "data/prediction_cells.gpkg"),
  root = "space-time-plot-sampling",
  mode = "mirror"
)

# Keep the website copy synchronized with the authoritative technical document.
file.copy("../modeling_and_software_details/magComp.pdf",
          file.path(downloads, "magComp.pdf"), overwrite = TRUE)
