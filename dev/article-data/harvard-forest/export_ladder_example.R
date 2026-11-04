# Small tables for the article, using layout 1 from the repeated-study runner.
Sys.setenv(TZ = "UTC")
library(tidyverse)
library(magComp)
source("ladder_settings.R")
source("ladder_functions.R")

output_dir <- "output/repeated_ladder"
example_dir <- "data/ladder"
dir.create(example_dir, recursive = TRUE, showWarnings = FALSE)
inventory <- readRDS(file.path(output_dir, "fits/0001/inventory.rds"))
write_csv(inventory$trees, file.path(example_dir, "trees.csv"))
write_csv(inventory$prediction, file.path(example_dir, "prediction.csv"))

truth <- read_csv("data/truth.csv", show_col_types = FALSE)
forest_area <- sum(inventory$area_summary$area_ha)
cell_truth <- bind_rows(
  truth |> transmute(domain = paste("Stratum", stratum_id),
    category = paste(species, dbh_class, sep = ":"), truth = density),
  truth |> group_by(species, dbh_class) |>
    summarise(truth = sum(count) / forest_area, .groups = "drop") |>
    transmute(domain = "Forest", category = paste(species, dbh_class, sep = ":"), truth)
) |> mutate(target = "Species-by-DBH")
target_truth <- bind_rows(cell_truth,
  cell_truth |> group_by(domain) |> summarise(truth = sum(truth), .groups = "drop") |>
    mutate(category = "All", target = "Overall density"))

results <- map(paste0("H", 1:5), function(model) {
  readRDS(file.path(output_dir, "results", paste0("0001_", model, ".rds")))
})
stopifnot(all(map_lgl(results, ~ .x$status == "complete")))
comparison <- map_dfr(results, function(x) x$estimates |> mutate(model = x$model)) |>
  left_join(target_truth, by = c("domain", "category", "target")) |>
  mutate(error = mean - truth,
    covered = (lower <= truth | near(lower, truth)) &
              (truth <= upper | near(truth, upper)), width = upper - lower)
write_csv(comparison, file.path(example_dir, "comparison.csv"))
map_dfr(results, function(x) x$diagnostics |> mutate(model = x$model)) |>
  write_csv(file.path(example_dir, "diagnostics.csv"))
map_dfr(results, function(x) slice_tail(x$attempts, n = 1)) |>
  mutate(chains = n_chains, thin = n_thin, prediction_draws_per_chain = prediction_draws) |>
  write_csv(file.path(example_dir, "settings.csv"))
map_dfr(paste0("H", 1:5), function(model) {
  read_csv(file.path(output_dir, "fits/0001", paste0(model, "_parameters.csv")),
    show_col_types = FALSE) |> mutate(model = .env$model)
}) |> write_csv(file.path(example_dir, "parameters.csv"))

# Keep one draw for the worked observed-plus-predicted count calculation.
totals <- readRDS(file.path(output_dir, "fits/0001/H5_totals.rds"))
observed <- inventory$trees |> mutate(category = paste(species, dbh_class, sep = ":"))
known <- xtabs(count ~ stratum + category, data = observed)
tibble(observed = known["2", "Hemlock:12.7-20"],
  predicted_unsampled = totals["2", "Hemlock:12.7-20", 1, 1] - observed,
  population_total = observed + predicted_unsampled,
  density = population_total / inventory$area_summary$area_ha[inventory$area_summary$stratum_id == 2]) |>
  write_csv(file.path(example_dir, "one_draw.csv"))

# Plot-level fitted means, separate from census accuracy.
fitted_summaries <- list()
for (model in paste0("H", 1:5)) {
  fit <- readRDS(file.path(output_dir, "fits/0001", paste0(model, ".rds")))
  expected <- fitted(fit, target = "expected", type = "count", summary = "mean")
  means <- expected$count
  # Build the observed matrix from labeled input counts.
  input_counts <- inventory$trees |> mutate(category = paste(species, dbh_class, sep = ":"))
  observed <- xtabs(count ~ unit_id + category, data = input_counts)
  observed <- observed[rownames(means), colnames(means)]
  fitted_summaries[[model]] <- tibble(model, unit_id = rownames(means),
    observed = rowSums(observed), mean = rowSums(means))
  rm(fit, expected)
  gc()
}
bind_rows(fitted_summaries) |> write_csv(file.path(example_dir, "fitted_counts.csv"))
cat("Exported the first-layout summaries to", example_dir, "\n")
