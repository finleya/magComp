# Functions used by the five-model census experiment.
source("population_summaries.R")

make_ladder_inventory <- function(repetition) {
  set.seed(repetition)
  inventory <- draw_inventory(forest, strata, census, prediction_cells, spacing, radius)
  origin <- st_bbox(forest)[c("xmin", "ymin")]
  observed_xy <- st_coordinates(st_centroid(st_geometry(inventory$plot_units)))
  prediction_xy <- st_coordinates(st_centroid(st_geometry(inventory$prediction_units)))
  observed_coordinates <- tibble(unit_id = inventory$plot_units$unit_id,
    easting_km = (observed_xy[, 1] - origin[1]) / 1000,
    northing_km = (observed_xy[, 2] - origin[2]) / 1000)
  inventory$trees <- inventory$trees |>
    left_join(observed_coordinates, by = "unit_id")
  inventory$prediction <- inventory$prediction |>
    mutate(easting_km = (prediction_xy[, 1] - origin[1]) / 1000,
           northing_km = (prediction_xy[, 2] - origin[2]) / 1000)
  stopifnot(!anyDuplicated(observed_coordinates[c("easting_km", "northing_km")]),
            !anyDuplicated(inventory$prediction[c("easting_km", "northing_km")]))
  inventory[c("plots", "trees", "prediction", "area_summary")]
}

ladder_formulas <- function(model) {
  magnitude <- list(
    H1 = ~ 1,
    H2 = ~ 1 + iid(stratum_id, name = "stratum"),
    H3 = ~ 1 + iid(stratum_id, name = "stratum"),
    H4 = ~ 1 + iid(stratum_id, name = "stratum") +
      nngp(easting_km, northing_km, m = 15, name = "space"),
    H5 = ~ 1 + iid(stratum_id, name = "stratum") +
      nngp(easting_km, northing_km, m = 15, name = "space")
  )
  species_models <- list(
    H1 = species ~ 1,
    H2 = species ~ 1 + factor_iid(stratum_id, K = 1, anchors = "Hemlock", name = "stratum"),
    H3 = species ~ 1 + factor_iid(stratum_id, K = 1, anchors = "Hemlock", name = "stratum"),
    H4 = species ~ 1 + factor_iid(stratum_id, K = 1, anchors = "Hemlock", name = "stratum"),
    H5 = species ~ 1 + factor_iid(stratum_id, K = 1, anchors = "Hemlock", name = "stratum") +
      factor_nngp(easting_km, northing_km, K = 1, m = 15, anchors = "Hemlock", name = "space")
  )
  dbh <- if (model %in% c("H1", "H2")) dbh_class | species ~ 1 else
    dbh_class | species ~ 1 + factor_iid(stratum_id, K = 1, name = "stratum")
  list(magnitude = magnitude[[model]],
       composition = magCompComposition(species_models[[model]], dbh))
}

fit_ladder_model <- function(trees, model, control) {
  dat <- magCompData(trees, unit = unit_id, categories = c(species, dbh_class),
                     tree_count = count, design = fixedArea(area_ha))
  formulas <- ladder_formulas(model)
  magComp(dat, magnitude = formulas$magnitude,
           composition = formulas$composition, control = control)
}

# Classical unsplit R-hat, matching the package. A common constant has no
# diagnostic; different constants in different chains signal disagreement.
chain_rhat <- function(x) {
  n <- nrow(x)
  within <- mean(apply(x, 2, var))
  between <- n * var(colMeans(x))
  if (within == 0) return(if (between == 0) NA_real_ else Inf)
  sqrt(max(1, ((n - 1) / n * within + between / n) / within))
}
array_rhats <- function(x) {
  dims <- dim(x)
  x <- array(x, c(prod(head(dims, -2)), tail(dims, 2)))
  apply(x, 1, chain_rhat)
}
largest_rhat <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

check_fit_chains <- function(fit, model) {
  retained <- (fit$control$n_samples - fit$control$n_burn) / fit$control$n_thin
  thin <- max(1, floor(retained / diagnostic_draws))
  parameters <- summary(fit, thin = thin)
  # Variances, NB size, spatial decay, and spatial loadings also govern
  # uncertainty away from the observed plots.
  covariance <- parameters |>
    filter(parameter %in% c("size", "variance", "phi") |
             (parameter == "loadings" & term == "space"))
  combined_max <- NA_real_
  if (model != "H1") {
    combined <- samples(fit, type = "combined", thin = thin)
    combined_max <- largest_rhat(vapply(combined, chain_rhat, numeric(1)))
  }
  fitted_values <- predict(fit, fit$data$unit_data, target = "expected",
                           type = "density", summary = "samples",
                           n_draws = min(diagnostic_draws, retained))
  stopifnot(all(is.finite(fitted_values$density)), all(fitted_values$density >= 0))
  fitted_max <- largest_rhat(array_rhats(fitted_values$density))
  diagnostics <- tibble(scope = c("Parameters", "Covariance and size", "Combined intercepts", "Fitted densities"),
    max_rhat = c(largest_rhat(parameters$rhat), largest_rhat(covariance$rhat),
                 combined_max, fitted_max))
  good <- all(diagnostics$max_rhat[diagnostics$scope != "Parameters"] <= rhat_limit,
              na.rm = TRUE)
  list(parameters = parameters, diagnostics = diagnostics, good = good)
}

check_total_chains <- function(totals) {
  # Include every cell, forest marginals, and overall stratum/forest counts.
  r <- c(array_rhats(totals), array_rhats(apply(totals, c(2, 3, 4), sum)),
         array_rhats(apply(totals, c(1, 3, 4), sum)),
         chain_rhat(apply(totals, c(3, 4), sum)))
  largest_rhat(r)
}

run_ladder_model <- function(repetition, model, inventory, output_dir, chain_cores) {
  started <- proc.time()
  folder <- file.path(output_dir, "fits", sprintf("%04d", repetition))
  choices <- bind_rows(model_settings |> filter(.data$model == .env$model) |>
                        select(n_samples, n_burn), extension_settings)
  attempts <- list()
  for (attempt in seq_len(nrow(choices))) {
    cat("Starting", repetition, model, "attempt", attempt, "at", format(Sys.time()), "\n")
    control <- magCompControl(
      n_samples = choices$n_samples[attempt], n_burn = choices$n_burn[attempt],
      n_thin = n_thin, n_chains = n_chains, n_cores = chain_cores,
      chain_seeds = chain_seeds,
      nb_size = fixed_size, verbose = FALSE
    )
    fit <- fit_ladder_model(inventory$trees, model, control)
    fit_file <- file.path(folder, paste0(model, ".rds"))
    saveRDS(fit, paste0(fit_file, ".tmp"))
    file.rename(paste0(fit_file, ".tmp"), fit_file)
    checked <- check_fit_chains(fit, model)
    total_rhat <- NA_real_
    totals <- NULL
    if (checked$good || attempt == nrow(choices)) {
      set.seed(repetition)
      totals <- predict_totals(fit, inventory, prediction_draws)
      total_rhat <- check_total_chains(totals)
    }
    good <- checked$good && !is.na(total_rhat) && total_rhat <= rhat_limit
    attempts[[attempt]] <- bind_cols(tibble(attempt, model,
      n_samples = choices$n_samples[attempt], n_burn = choices$n_burn[attempt],
      accepted = good, predictive_rhat = total_rhat),
      checked$diagnostics |> pivot_wider(names_from = scope, values_from = max_rhat))
    write_csv(bind_rows(attempts), file.path(folder, paste0(model, "_attempts.csv")))
    cat("Finished", repetition, model, "attempt", attempt, "accepted:", good,
        "at", format(Sys.time()), "\n")
    if (good || attempt == nrow(choices)) break
    file.rename(file.path(folder, paste0(model, ".rds")),
                file.path(folder, paste0(model, "_attempt", attempt, ".rds")))
    rm(fit, totals)
    gc()
  }
  saveRDS(totals, file.path(folder, paste0(model, "_totals.rds")))
  write_csv(checked$parameters, file.path(folder, paste0(model, "_parameters.csv")))
  elapsed <- proc.time() - started
  list(repetition = repetition, model = model,
       status = if (good) "complete" else "unresolved", message = "",
       estimates = summarize_inventory(totals, inventory$area_summary),
       diagnostics = bind_rows(checked$diagnostics,
         tibble(scope = "Population predictions", max_rhat = total_rhat)),
       attempts = bind_rows(attempts),
       elapsed_seconds = unname(elapsed["elapsed"]),
       cpu_seconds = sum(elapsed[c("user.self", "sys.self", "user.child", "sys.child")]))
}
