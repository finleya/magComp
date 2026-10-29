test_that("prediction and aggregation include fitted latent terms", {
  rows <- expand.grid(plot = paste0("p", 1:8), cls = c("a", "b"))
  rows$count <- rep(1:8, 2)
  rows$area <- rep(seq(0.1, 0.8, 0.1), 2)
  rows$stand <- rep(rep(paste0("s", 1:4), each = 2), 2)
  rows$x <- rep(seq(-1, 1, length.out = 8), 2)
  dat <- magCompData(rows, unit = plot, category = cls,
                    tree_count = count, design = fixedArea(area))
  fit <- magComp(dat, magnitude = ~ x + iid(stand),
    composition = ~ x + factor_iid(stand, K = 1),
    control = magCompControl(n_samples = 10, n_burn = 5, n_chains = 2,
      chain_seeds = c(101, 102), verbose = FALSE))
  expected <- fitted(fit, summary = "samples")$count
  expect_equal(predict(fit, dat$unit_data, summary = "samples")$count, expected)
  magnitude <- mc_regression_magnitude_samples(fit)
  group <- fit$magnitude$process_terms[[1L]]$group
  effects <- fit$backend$magnitude_process_effect_samples$iid_1
  magnitude <- magnitude * exp(effects[group,,,drop = FALSE])
  expect_equal(fitted(fit, type = "magnitude", summary = "samples")$magnitude,
               magnitude, ignore_attr = TRUE)
  expect_equal(apply(expected, c(1L, 3L, 4L), sum),
               sweep(magnitude, 1L, dat$area, "*"), ignore_attr = TRUE)
  stand <- standSummary(predict(fit, dat$unit_data, summary = "samples"),
    data = dat$unit_data, unit = plot, group = stand, represented_area = area,
    type = "total", summary = "samples")
  expect_equal(stand$total,
    standSummary(fitted(fit, summary = "samples"), data = dat$unit_data,
      unit = plot, group = stand, represented_area = area, type = "total", summary = "samples")$total)

  expect_error(fitted(fit, newdata = dat$unit_data)$count, "unused fitted argument")
  expect_error(predict(fit, dat$unit_data, typo = TRUE)$count, "unused predict argument")

  expect_error(standSummary(fitted(fit, summary = "samples"), data = dat$unit_data, unit = plot, represented_area = area, typo = TRUE), "unused standSummary argument")
})
