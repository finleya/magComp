# Summarize leave-one-plot-out predictions from cross_validate_ladder.R.
Sys.setenv(TZ = "UTC")
library(tidyverse)
library(magComp)

trees <- read_csv("data/ladder/trees.csv", show_col_types = FALSE)
portions <- trees |> distinct(unit_id, plot_id, stratum_id, area_ha)
results <- portion_results <- diagnostics <- list()

for (model in paste0("H", 1:5)) {
  cat("Summarizing", model, "\n")
  cv <- readRDS(file.path("output/cv", paste0(model, ".rds")))

  results[[model]] <- map_dfr(names(cv$predictions), function(plot) {
    prediction <- cv$predictions[[plot]]
    units <- prediction$metadata$unit_ids
    categories <- dimnames(prediction$count)[[2]]
    area <- sum(portions$area_ha[match(units, portions$unit_id)])

    # Sum a plot's portions within each draw, then divide by its sampled area.
    # The remaining dimensions are category, draw, and chain.
    density <- apply(prediction$count, c(2, 3, 4), sum) / area
    observed <- colSums(cv$observed$count[units, categories, drop = FALSE]) / area
    tibble(model, plot_id = as.integer(plot), category = categories,
      area_ha = area, observed = unname(observed),
      mean = apply(density, 1, mean),
      lower = apply(density, 1, quantile, 0.025),
      upper = apply(density, 1, quantile, 0.975))
  })

  # Assess each portion within its own stratum, using its own sampled area.
  portion_results[[model]] <- map_dfr(names(cv$predictions), function(plot) {
    prediction <- cv$predictions[[plot]]
    units <- prediction$metadata$unit_ids
    categories <- dimnames(prediction$count)[[2]]
    map_dfr(seq_along(units), function(i) {
      portion <- portions[match(units[i], portions$unit_id), ]
      density <- prediction$count[i, , , ] / portion$area_ha
      observed <- cv$observed$count[units[i], categories] / portion$area_ha
      tibble(model, plot_id = as.integer(plot), unit_id = portion$unit_id,
        stratum_id = portion$stratum_id, category = categories,
        area_ha = portion$area_ha, observed = unname(observed),
        mean = apply(density, 1, mean),
        lower = apply(density, 1, quantile, 0.025),
        upper = apply(density, 1, quantile, 0.975))
    })
  })

  # Keep parameter diagnostics with the downloadable results.
  diagnostics[[model]] <- map_dfr(names(cv$fits), function(plot) {
    parameters <- summary(cv$fits[[plot]], thin = 3)
    parameters |>
      slice_max(rhat, n = 1, with_ties = FALSE) |>
      mutate(model = .env$model, plot_id = as.integer(plot), .before = 1)
  })
  rm(cv)
  gc()
}

predictions <- bind_rows(results) |>
  mutate(error = mean - observed,
    covered = (lower <= observed | near(lower, observed)) &
              (observed <= upper | near(observed, upper)),
    width = upper - lower)

# Each of the 48 plots has 30 species-by-DBH cells, including zero counts.
# Scores give equal weight to those 1,440 observed cells for each model.
assessment <- predictions |>
  group_by(model) |>
  summarise(plots = n_distinct(plot_id), cells = n(),
    mean_error = mean(error), rmspe = sqrt(mean(error^2)),
    coverage = 100 * mean(covered), mean_width = mean(width), .groups = "drop")

portion_predictions <- bind_rows(portion_results) |>
  mutate(error = mean - observed,
    covered = (lower <= observed | near(lower, observed)) &
              (observed <= upper | near(observed, upper)),
    width = upper - lower)

# Within a stratum, give equal weight to each portion's 30 category cells.
stratum_assessment <- portion_predictions |>
  group_by(stratum_id, model) |>
  summarise(portions = n_distinct(unit_id), cells = n(),
    mean_error = mean(error), rmspe = sqrt(mean(error^2)),
    coverage = 100 * mean(covered), mean_width = mean(width), .groups = "drop")

dir.create("data/cv", showWarnings = FALSE)
write_csv(portions, "data/cv/folds.csv")
write_csv(predictions, "data/cv/predictions.csv")
write_csv(assessment, "data/cv/summary.csv")
write_csv(portion_predictions, "data/cv/portion_predictions.csv")
write_csv(stratum_assessment, "data/cv/stratum_summary.csv")
write_csv(bind_rows(diagnostics), "data/cv/diagnostics.csv")
print(assessment)
print(stratum_assessment, n = 20)
