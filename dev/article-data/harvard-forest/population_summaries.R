# Predict the entire unsampled remainder in one call, preserving dependence.
predict_totals <- function(fit, inventory, n_draws) {
  predictions <- predict(fit, inventory$prediction, target = "observation",
                         type = "density", summary = "samples", n_draws = n_draws)
  unsampled <- standSummary(predictions, data = inventory$prediction, unit = unit_id,
                            group = stratum_id, represented_area = area_ha,
                            type = "total", summary = "samples")
  observed <- inventory$trees |>
    mutate(category = paste(species, dbh_class, sep = ":"))
  known <- xtabs(count ~ stratum + category, data = observed)
  ids <- dimnames(unsampled$total)
  totals <- sweep(unsampled$total, c(1, 2), known[ids[[1]], ids[[2]]], "+")
  # With these count exposures and represented areas, totals are integers.
  stopifnot(all(is.finite(totals)), all(totals >= 0),
            max(abs(totals - round(totals))) < 1e-7)
  totals
}

# Rows are categories; other dimensions are retained draw and chain.
summarize_draws <- function(draws) {
  tibble(category = dimnames(draws)[[1]], mean = apply(draws, 1, mean),
         lower = apply(draws, 1, quantile, 0.025),
         upper = apply(draws, 1, quantile, 0.975))
}

summarize_inventory <- function(totals, area_summary) {
  areas <- setNames(area_summary$area_ha, area_summary$stratum_id)
  forest <- apply(totals, c(2, 3, 4), sum) / sum(areas)
  cells <- bind_rows(
    summarize_draws(forest) |> mutate(domain = "Forest"),
    map_dfr(dimnames(totals)[[1]], function(h) {
      summarize_draws(totals[h, , , ] / areas[h]) |>
        mutate(domain = paste("Stratum", h))
    })
  ) |> mutate(target = "Species-by-DBH")
  overall <- map_dfr(c("Forest", dimnames(totals)[[1]]), function(h) {
    draws <- if (h == "Forest") apply(forest, c(2, 3), sum) else
      apply(totals[h, , , ], c(2, 3), sum) / areas[h]
    tibble(domain = if (h == "Forest") h else paste("Stratum", h),
           category = "All", mean = mean(draws),
           lower = quantile(draws, 0.025), upper = quantile(draws, 0.975),
           target = "Overall density")
  })
  bind_rows(cells, overall)
}
