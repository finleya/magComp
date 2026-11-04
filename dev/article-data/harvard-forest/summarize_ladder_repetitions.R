# Run after repeat_ladder.R, from the same directory.
# Optional arguments identify the first and last requested layouts.
# Third and fourth arguments set input and summary directories for a separate run.
Sys.setenv(TZ = "UTC")
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
first_layout <- if (length(args)) as.integer(args[1]) else 1L
last_layout <- if (length(args) > 1) as.integer(args[2]) else 1000L
output_dir <- if (length(args) > 2) args[3] else "output/repeated_ladder"
summary_dir <- if (length(args) > 3) args[4] else "data/repeated_ladder"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
results <- list.files(file.path(output_dir, "results"), pattern = "\\.rds$", full.names = TRUE) |>
  map(readRDS)
settings <- readRDS(file.path(output_dir, "settings.rds"))

status <- map_dfr(results, function(x) {
  tibble(repetition = x$repetition, model = x$model, status = x$status,
         message = x$message, elapsed_seconds = x$elapsed_seconds %||% NA_real_,
         cpu_seconds = x$cpu_seconds %||% NA_real_,
         peak_parent_mb = x$peak_parent_mb %||% NA_real_)
})
status <- expand_grid(repetition = first_layout:last_layout, model = paste0("H", 1:5)) |>
  left_join(status, by = c("repetition", "model")) |>
  mutate(status = replace_na(status, "pending"))
write_csv(status, file.path(summary_dir, "status.csv"))
run_counts <- status |> count(model, status, name = "layouts")
write_csv(run_counts, file.path(summary_dir, "run_counts.csv"))

# Include all completed fits, retaining diagnostic flags for review.
# Paired comparisons require estimates from all five models in each layout.
paired <- status |>
  group_by(repetition) |>
  summarise(complete = n() == 5 && all(status %in% c("complete", "unresolved")),
            .groups = "drop") |>
  filter(complete)
estimates <- map_dfr(results, function(x) {
  if (is.null(x$estimates)) return(tibble())
  x$estimates |> mutate(repetition = x$repetition, model = x$model, status = x$status)
})
write_csv(estimates, file.path(output_dir, "all_estimates.csv"))

truth <- read_csv("data/truth.csv", show_col_types = FALSE)
areas <- read_csv("data/stratum_areas.csv", show_col_types = FALSE)
cell_truth <- bind_rows(
  truth |> transmute(domain = paste("Stratum", stratum_id),
                    category = paste(species, dbh_class, sep = ":"), truth = density),
  truth |> group_by(species, dbh_class) |>
    summarise(truth = sum(count) / sum(areas$area_ha), .groups = "drop") |>
    transmute(domain = "Forest", category = paste(species, dbh_class, sep = ":"), truth)
) |> mutate(target = "Species-by-DBH")
target_truth <- bind_rows(cell_truth,
  cell_truth |> group_by(domain) |> summarise(truth = sum(truth), .groups = "drop") |>
    mutate(category = "All", target = "Overall density"))

assessment <- estimates |>
  semi_join(paired, by = "repetition") |>
  left_join(target_truth, by = c("domain", "category", "target")) |>
  mutate(error = mean - truth, squared_error = error^2,
         # Include equal endpoints despite rounding in count-to-density conversion.
         covered = (lower <= truth | near(lower, truth)) &
                   (truth <= upper | near(truth, upper)), width = upper - lower)
stopifnot(!anyNA(assessment$truth), all(assessment$lower <= assessment$upper))
saveRDS(assessment, file.path(output_dir, "assessment.rds"))

# Keep zero-count census cells in coverage and error calculations.
# Relative bias for an individual zero-density cell is undefined.
target_summary <- assessment |>
  group_by(model, target, domain, category) |>
  summarise(layouts = n(), truth = first(truth), bias = mean(error),
    bias_mcse = sd(error) / sqrt(n()), rmse = sqrt(mean(squared_error)),
    rmse_mcse = sd(squared_error) / (2 * rmse * sqrt(n())),
    coverage = 100 * mean(covered), coverage_mcse = 100 * sd(covered) / sqrt(n()),
    mean_width = mean(width), width_mcse = sd(width) / sqrt(n()), .groups = "drop") |>
  mutate(relative_bias = if_else(truth > 0, 100 * bias / truth, NA_real_))
write_csv(target_summary, file.path(summary_dir, "target_summary.csv"))

# Average cells within a layout before calculating Monte Carlo uncertainty.
# Cells from one inventory are dependent; layouts are the repetitions.
by_inventory <- assessment |>
  group_by(repetition, model, target, domain) |>
  summarise(n_targets = n(), truth = mean(truth), error = mean(error),
    squared_error = mean(squared_error), covered = mean(covered),
    width = mean(width), .groups = "drop")
domain_summary <- by_inventory |>
  group_by(model, target, domain) |>
  summarise(layouts = n(), n_targets = first(n_targets), truth = first(truth),
    bias = mean(error), bias_mcse = sd(error) / sqrt(n()),
    rmse = sqrt(mean(squared_error)),
    rmse_mcse = sd(squared_error) / (2 * rmse * sqrt(n())),
    coverage = 100 * mean(covered), coverage_mcse = 100 * sd(covered) / sqrt(n()),
    mean_width = mean(width), width_mcse = sd(width) / sqrt(n()), .groups = "drop") |>
  mutate(relative_bias = 100 * bias / truth, relative_bias_mcse = 100 * bias_mcse / truth)
write_csv(domain_summary, file.path(summary_dir, "domain_summary.csv"))

sampling <- map_dfr(results, function(x) {
  if (is.null(x$sampling)) return(tibble())
  x$sampling |> mutate(repetition = x$repetition)
}) |> distinct(repetition, stratum_id, .keep_all = TRUE)
sampling |> group_by(stratum_id) |>
  summarise(layouts = n(), min_portions = min(plot_portions),
    median_portions = median(plot_portions), max_portions = max(plot_portions),
    unsampled = sum(no_sampled_portion), .groups = "drop") |>
  write_csv(file.path(summary_dir, "sampling_summary.csv"))

diagnostics <- map_dfr(results, function(x) {
  if (is.null(x$diagnostics)) return(tibble())
  x$diagnostics |> mutate(repetition = x$repetition, model = x$model)
})
write_csv(diagnostics, file.path(summary_dir, "diagnostics.csv"))
diagnostics |> filter(!is.na(max_rhat)) |> group_by(model, scope) |>
  summarise(fits = n(), median_max_rhat = median(max_rhat),
    p95_max_rhat = quantile(max_rhat, .95), max_rhat = max(max_rhat), .groups = "drop") |>
  write_csv(file.path(summary_dir, "diagnostic_summary.csv"))
attempts <- map_dfr(results, function(x) {
  if (is.null(x$attempts)) return(tibble())
  x$attempts |> mutate(repetition = x$repetition)
})
write_csv(attempts, file.path(summary_dir, "attempts.csv"))
settings$models |>
  mutate(chains = settings$n_chains, thin = settings$n_thin,
    prediction_draws_per_chain = settings$prediction_draws,
    requested_layouts = last_layout - first_layout + 1,
    complete_pairs = nrow(paired), package_version = settings$package_version) |>
  write_csv(file.path(summary_dir, "settings.csv"))
print(run_counts)
print(domain_summary |> select(model, target, domain, bias, relative_bias, rmse, coverage, mean_width), n = 50)
