mc_block_coef <- function(fit, block = names(fit$data$allocation_blocks)[1L])
{
  fit$backend$composition_block_coefficient_samples[[block]]
}

mc_block_loading <- function(fit, block = names(fit$data$allocation_blocks)[1L])
{
  fit$backend$composition_block_loading_samples[[block]]
}

mc_block_factor <- function(fit, block = names(fit$data$allocation_blocks)[1L])
{
  fit$backend$composition_block_score_samples[[block]]
}

mc_block_factor_phi <- function(fit,
                                block = names(fit$data$allocation_blocks)[1L])
{
  fit$backend$composition_block_covariance_samples[[block]]
}

mc_loading_entry <- function(x, row, column)
{
  if(length(dim(x)) == 2L)
    return(x[row,column])
  if(length(dim(x)) == 3L)
    return(x[row,column,])
  x[row,column,,]
}

expect_loading_constraints <- function(x, anchors, factor_names = anchors)
{
  for(k in seq_along(anchors)){
    diagonal <- as.numeric(mc_loading_entry(x, anchors[k], factor_names[k]))
    expect_true(all(is.finite(diagonal) & diagonal > 0))
    if(k > 1L){
      for(h in seq_len(k - 1L)){
        upper <- as.numeric(mc_loading_entry(x, anchors[h], factor_names[k]))
        expect_equal(upper, rep(0, length(upper)), tolerance = 0)
      }
    }
  }
}

test_that("factor starting transformation preserves the low-rank surface", {
  set.seed(67)

  loading <- matrix(rnorm(14), nrow = 7, ncol = 2)
  loading[2,] <- c(0.8, -0.4)
  loading[5,] <- c(0.3, 1.1)
  rownames(loading) <- paste0("species", seq_len(nrow(loading)))
  colnames(loading) <- c("factor1", "factor2")

  factor_score <- matrix(rnorm(40), nrow = 20, ncol = 2)
  colnames(factor_score) <- colnames(loading)

  fitted_surface <- factor_score %*% t(loading)
  constrained <- magComp:::magcomp_constrain_factor_start(
    loading,
    factor_score,
    anchor_index = c(2L, 5L),
    preserve_surface = TRUE
  )

  expect_gt(constrained$loading["species2", "factor1"], 0)
  expect_equal(constrained$loading["species2", "factor2"], 0)
  expect_gt(constrained$loading["species5", "factor2"], 0)
  expect_equal(
    constrained$factor %*% t(constrained$loading),
    fitted_surface,
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("constraint enforcement preserves valid lower-triangular starts", {
  loading <- matrix(
    c(0.45, 0.4, -0.2, 0, 1.7, 0.3),
    nrow = 3,
    ncol = 2,
    dimnames = list(paste0("species", 1:3), c("factor1", "factor2"))
  )
  factor_score <- matrix(1:8, nrow = 4, ncol = 2)

  constrained <- magComp:::magcomp_constrain_factor_start(
    loading,
    factor_score,
    anchor_index = c(1L, 2L)
  )

  expect_equal(constrained$loading["species2", "factor1"], 0.4)
  expect_equal(constrained$loading["species1", "factor2"], 0)
  expect_equal(constrained$loading, loading)
  expect_equal(constrained$factor, factor_score)
})

mc_mag_process <- function(fit, term = names(fit$backend$magnitude_process_effect_samples)[1L])
{
  fit$backend$magnitude_process_effect_samples[[term]]
}

mc_mag_sigma <- function(fit, term = names(fit$backend$magnitude_process_variance_samples)[1L])
{
  fit$backend$magnitude_process_variance_samples[[term]]
}

mc_mag_phi <- function(fit, term = names(fit$backend$magnitude_process_covariance_samples)[1L])
{
  fit$backend$magnitude_process_covariance_samples[[term]]
}

test_that("magComp returns an inspectable model specification", {
  trees <- data.frame(
    plot = rep(letters[1:3], each = 3),
    cls = rep(c("s", "m", "l"), 3),
    x = rep(1:3, each = 3),
    y = rep(1:3, each = 3),
    area = rep(0.1, 9),
    p95 = rep(c(12, 15, 18), each = 3)
  )
  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, y),
    design = fixedArea(area)
  )

  fit <- magComp(
    data = dat,
    family = fixedAreaStems(),
    magnitude = ~ p95 + nngp(x, y, m = 5),
    composition = ~ p95 + factor_nngp(x, y, K = 1),
    control = magCompControl(run_sampler = FALSE)
  )

  expect_s3_class(fit, "magComp")
  expect_false(fit$fitted)
  expect_equal(length(fit$magnitude$process_terms), 1L)
  expect_equal(length(fit$composition$process_terms), 1L)
  expect_equal(fit$magnitude$process_terms[[1L]]$term_type, "nngp")
  expect_equal(fit$magnitude$process_terms[[1L]]$m, 5L)
  expect_equal(length(fit$magnitude$process_terms[[1L]]$order), 3L)
})

test_that("negative-binomial size controls validate", {
  expect_identical(magCompPriors()$magnitude$size, "auto")
  expect_null(magCompPriors(magnitude = list(size = NULL))$magnitude$size)
  expect_equal(magCompPriors(magnitude = list(size = c(2, 0.2)))$magnitude$size,
               c(shape = 2, rate = 0.2))
  expect_error(magCompPriors(magnitude = list(size = c(2, -0.2))), "size")
  expect_error(magCompPriors(magnitude = list(size = "automatic")), "size")
  expect_error(magCompControl(nb_size_tuning = 0), "nb_size_tuning")
})

test_that("magComp parses block-specific composition formulas", {
  trees <- data.frame(
    plot = c("p1", "p1", "p1", "p2", "p2", "p3", "p3", "p3"),
    species_group = c("hw", "hw", "sw", "hw", "sw", "sw", "sw", "hw"),
    dbh_class = c("small", "large", "small", "large",
                  "large", "small", "large", "small"),
    z = c(1, 1, 1, 2, 2, 3, 3, 3),
    x = c(0, 0, 0, 1, 1, 2, 2, 2),
    y = c(0, 0, 0, 1, 1, 0, 0, 0),
    area = 0.1
  )
  dat <- magCompData(
    trees,
    unit = plot,
    categories = c(species_group, dbh_class),
    coords = c(x, y),
    design = fixedArea(area)
  )
  comp <- magCompComposition(
    species_group ~ z,
    dbh_class | species_group ~ z + factor_gp(x, y, K = 1)
  )
  fit <- magComp(
    dat,
    magnitude = ~ z,
    composition = comp,
    control = magCompControl(run_sampler = FALSE)
  )

  expect_s3_class(comp, "magCompComposition")
  expect_equal(names(comp), c("species_group", "dbh_class|species_group"))
  expect_equal(fit$composition$allocation_block_names, names(comp))
  expect_true(fit$composition$block_specific)
  expect_equal(nrow(fit$composition$allocation_components[[1L]]$X),
               dat$n_unit)
  expect_equal(nrow(fit$composition$allocation_components[[2L]]$X),
               dat$n_unit * length(dat$category_levels$species_group))
  expect_equal(fit$composition$allocation_components[[2L]]$process_terms[[1L]]$term_type,
               "factor_gp")
  expect_equal(
    fit$composition$allocation_components[[2L]]$process_terms[[1L]]$n_group,
    dat$n_unit
  )

  expect_error(
    magComp(
      dat,
      composition = magCompComposition(species_group ~ z),
      control = magCompControl(run_sampler = FALSE)
    ),
    "missing composition formula"
  )
  expect_error(
    magComp(
      dat,
      composition = magCompComposition(dbh_class ~ z,
                                       species_group | dbh_class ~ z),
      control = magCompControl(run_sampler = FALSE)
    ),
    "missing composition formula"
  )
  fit_sampler <- magComp(
    dat,
    composition = comp,
    control = magCompControl(n_samples = 3,
                             n_burn = 1,
                             verbose = FALSE),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )
  expect_true(fit_sampler$fitted)

  fit_anchor <- magComp(
    dat,
    composition = magCompComposition(
      species_group ~ 1,
      dbh_class | species_group ~
        factor_iid(plot, K = 2, anchors = c("hw:large", "sw:large"))
    ),
    control = magCompControl(run_sampler = FALSE)
  )
  anchor_component <-
    fit_anchor$composition$allocation_components[["dbh_class|species_group"]]
  anchor_term <- anchor_component$process_terms[[1L]]
  expect_equal(anchor_term$anchor_index, c(1L, 2L))
  expect_equal(anchor_term$anchors, c("hw:large", "sw:large"))
})

test_that("block-specific smoothing respects ordered category metadata", {
  trees <- data.frame(
    plot = rep(paste0("p", 1:6), each = 6),
    species_group = rep(rep(c("hw", "sw"), each = 3), times = 6),
    dbh_class = rep(c("small", "medium", "large"), times = 12),
    z = rep(seq_len(6), each = 6),
    area = 0.1
  )
  trees$dbh_class <- factor(trees$dbh_class, levels = c("large", "medium", "small"))
  dat <- magCompData(
    trees,
    unit = plot,
    categories = c(species_group, dbh_class),
    ordered = c(species_group = FALSE, dbh_class = TRUE),
    design = fixedArea(area)
  )
  comp <- magCompComposition(
    species_group ~ z,
    dbh_class | species_group ~ z
  )

  expect_error(
    magComp(
      dat,
      magnitude = ~ z,
      composition = comp,
      smooth = magCompSmooth(intercept = "rw1"),
      control = magCompControl(run_sampler = FALSE)
    ),
    "nominal allocation block"
  )

  fit <- magComp(
    dat,
    magnitude = ~ z,
    composition = comp,
    smooth = magCompSmooth(dbh_class = list(intercept = "rw1", slope = "rw2")),
    control = magCompControl(run_sampler = FALSE)
  )
  expect_equal(fit$data$allocation_blocks$species_group$smooth$intercept, "none")
  expect_equal(fit$data$allocation_blocks[["dbh_class|species_group"]]$smooth$intercept,
               "rw1")
  expect_equal(fit$data$allocation_blocks[["dbh_class|species_group"]]$smooth$slope,
               "rw2")

  fit_inherit <- magComp(
    dat,
    magnitude = ~ z,
    composition = comp,
    smooth = magCompSmooth(intercept = "none", slope = "rw2",
                           species_group = list(intercept = "none",
                                                slope = "none"),
                           dbh_class = list(intercept = "rw1")),
    control = magCompControl(run_sampler = FALSE)
  )
  expect_equal(fit_inherit$data$allocation_blocks[["dbh_class|species_group"]]$smooth$intercept,
               "rw1")
  expect_equal(fit_inherit$data$allocation_blocks[["dbh_class|species_group"]]$smooth$slope,
               "rw2")
})

test_that("prior variance controls require positive scalars", {
  expect_error(magCompPriors(magnitude = list(coefficient_variance = c(1, 2)), composition = list(slope_variance = c(1, 2))), "coefficient_variance")
  expect_error(magCompPriors(composition = list(intercept_smoothing_variance = 0)), "intercept_smoothing_variance")
  expect_error(magCompPriors(magnitude = list(iid_variance = c(2, -1))), "iid_variance")
  expect_error(magCompPriors(magnitude = list(spatial_variance = c(1, NA))),
               "spatial_variance")
  expect_s3_class(magCompPriors(composition = list(intercept_smoothing_variance = 2)), "magCompPriors")

  smooth <- magCompSmooth()
  expect_equal(smooth$intercept, "none")
  expect_equal(smooth$slope, "none")
  expect_equal(smooth$loading, "none")
  expect_error(magCompSmooth(intercept = "rw1", ordered = FALSE),
               "ordered categories")
})

test_that("MCMC control accepts R-side chain settings", {
  ctl <- magCompControl(n_chains = 2, n_cores = 2,
                        chain_seeds = c(11, 12),
                        metropolis = list(batch_length = 10,
                                          target_accept = 0.40))

  expect_s3_class(ctl, "magCompControl")
  expect_equal(ctl$n_chains, 2L)
  expect_equal(ctl$n_cores, 2L)
  expect_equal(ctl$chain_seeds, c(11L, 12L))
  expect_equal(ctl$composition_link, "probit")
  expect_true(ctl$metropolis$adapt)
  expect_equal(ctl$metropolis$batch_length, 10L)
  expect_equal(ctl$metropolis$target_accept, 0.40)
  expect_equal(magCompControl(composition_link = "probit")$composition_link,
               "probit")
  expect_equal(magCompControl(composition_link = "logit")$composition_link,
               "logit")
  expect_error(magCompControl(n_chains = 0), "n_chains")
  expect_error(magCompControl(n_cores = 0), "n_cores")
  expect_error(magCompControl(n_chains = 2, chain_seeds = 11),
               "chain_seeds")
  expect_error(magCompControl(composition_link = "cloglog"))
  expect_error(magCompControl(metropolis = list(target_accept = 1)),
               "target_accept")
  expect_error(magCompControl(metropolis = list(foo = 1)),
               "unknown metropolis")
})

test_that("sampler fits intercept-only fixed-area stem counts", {
  # Small count table for fitting and output checks.
  n <- 12
  cls <- c("s", "m", "l", "xl")
  area <- rep(0.08, n)
  rows <- expand.grid(plot = paste0("p", seq_len(n)), cls = cls)
  rows$cls <- factor(rows$cls, levels = sort(cls))
  i <- as.integer(rows$plot)
  rows$count <- 1L + (i + 2L * as.integer(rows$cls)) %% 7L
  rows$area <- area[i]
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     ordered = TRUE, design = fixedArea(area))

  set.seed(41)
  fit <- magComp(
    dat,
    control = magCompControl(
      n_samples = 40,
      n_burn = 10,
      n_thin = 1,
      nb_size = 1000,
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )

  expect_true(fit$fitted)
  expect_equal(fit$control$composition_link, "probit")
  expect_equal(as.character(fit$backend$composition_link), "probit")
  expect_true(fit$backend$profile_pg[["probit_composition_draws"]] > 0)
  expect_equal(fit$backend$profile_pg[["composition_draws"]], 0)
  expect_equal(dim(fit$backend$magnitude_coefficient_samples), c(1L, 30L))
  expect_equal(dim(mc_block_coef(fit)), c(length(cls) - 1L, 1L, 30L))
  expect_equal(dim(fit$backend$proportion_samples), c(n, length(cls), 30L))
  expect_equal(dim(fitted(fit, type = "count")$count), c(n, length(cls)))
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)
  expect_true(all(is.finite(fit$backend$magnitude_coefficient_samples)))
  expect_true(all(is.finite(mc_block_coef(fit))))
  expect_true(all(fitted(fit, type = "count")$count > 0))

  tab <- summary(fit)
  expect_s3_class(tab, "magCompSummary")
  expect_true(all(is.na(tab$rhat)))
  beta <- subset(tab, component == "magnitude" & parameter == "coefficients")
  beta_draws <- as.numeric(fit$backend$magnitude_coefficient_samples)
  expect_equal(beta$mean, mean(beta_draws))
  expect_equal(beta$sd, sd(beta_draws))
  expect_equal(as.numeric(beta[, c("q2_5", "q50", "q97_5")]),
               as.numeric(quantile(beta_draws, c(0.025, 0.5, 0.975))))
  expect_error(summary(fit, probs = numeric()), "probabilities")

  spatial_diag <- spatialDiagnostics(fit)
  expect_s3_class(spatial_diag, "magCompSpatialDiagnostics")
  expect_equal(nrow(spatial_diag$summary), 0L)
})

test_that("prediction and CV support point-sampling basal area", {
  set.seed(411)

  trees <- data.frame(
    point = rep(paste0("p", 1:6), each = 4),
    species = rep(c("a", "b", "a", "c"), 6),
    cls = factor(rep(c("small", "small", "large", "large"), 6),
                 levels = c("small", "large"), ordered = TRUE),
    baf = 10,
    tree_count = rep(c(2L, 1L, 1L, 2L), 6),
    z = rep(seq(-1, 1, length.out = 6), each = 4)
  )
  dat <- magCompData(
    trees,
    unit = point,
    categories = c(species, cls),
    tree_count = tree_count,
    design = pointSampling(baf),
    ordered = c(species = FALSE, cls = TRUE)
  )

  fit <- magComp(
    dat,
    family = pointBasalArea(),
    magnitude = ~ 1,
    composition = magCompComposition(
      species ~ 1,
      cls | species ~ 1
    ),
    control = magCompControl(
      n_samples = 6,
      n_burn = 2,
      n_thin = 1,
      n_chains = 1,
      nb_size = 1000,
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )

  expect_true(fit$fitted)
  expect_equal(unname(fit$data$sampler_area), rep(0.1, dat$n_unit))
  expect_true(all(is.finite(fitted(fit, type = "density")$density)))
  expect_true(all(is.finite(fitted(fit, type = "magnitude")$magnitude)))

  newdata <- trees[!duplicated(trees$point), c("point", "baf", "z")]
  pred_density <- predict(fit, newdata, type = "density",
                          summary = "samples")$density
  pred_count <- predict(fit, newdata, type = "count",
                        summary = "mean")$count

  expect_equal(dim(pred_density), c(dat$n_unit, ncol(dat$counts), 4L))
  expect_equal(dim(pred_count), c(dat$n_unit, ncol(dat$counts)))
  expect_true(all(is.finite(pred_density)))
  expect_true(all(is.finite(pred_count)))

  fold_id <- rep(1:2, length.out = dat$n_unit)
  names(fold_id) <- dat$unit_levels
  set.seed(91)
  cv <- magCompCV(
    dat,
    family = pointBasalArea(),
    magnitude = ~ 1,
    composition = magCompComposition(
      species ~ 1,
      cls | species ~ 1
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none"),
    control = magCompControl(
      n_samples = 4,
      n_burn = 2,
      n_thin = 1,
      n_chains = 1,
      nb_size = 1000,
      verbose = FALSE
    ),
    n_fold = 2,
    fold_id = fold_id,
    type = "density"
  )

  expect_s3_class(cv, "magCompCV")
  expect_equal(sum(vapply(cv$predictions, function(x) nrow(x$density), integer(1L))),
               dat$n_unit)
  expect_true(all(is.finite(unlist(lapply(cv$predictions, `[[`, "density")))))
  expect_true(all(is.finite(unlist(cv$totals))))
  expect_equal(cv$observed$density, dat$counts / 0.1)
})

test_that("sampler fits fixed-effect formulas", {
  # Small count table for fitting and output checks.
  n <- 12
  cls <- c("s", "m", "l", "xl")
  area <- rep(0.08, n)
  z <- seq(-1, 1, length.out = n)
  stand <- rep(c("A", "B", "C"), length.out = n)
  rows <- expand.grid(plot = paste0("p", seq_len(n)), cls = cls)
  rows$cls <- factor(rows$cls, levels = sort(cls))
  i <- as.integer(rows$plot)
  rows$count <- 1L + (i + 2L * as.integer(rows$cls)) %% 7L
  rows$area <- area[i]
  rows$z <- z[i]
  rows$stand <- stand[i]
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     ordered = TRUE, design = fixedArea(area))

  set.seed(42)
  fit <- magComp(
    dat,
    magnitude = ~ z,
    composition = ~ z,
    control = magCompControl(
      n_samples = 70,
      n_burn = 20,
      n_thin = 1,
      nb_size = 1000,
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )

  expect_equal(dim(fit$backend$magnitude_coefficient_samples), c(2L, 50L))
  expect_equal(dim(mc_block_coef(fit)),
               c(length(cls) - 1L, 2L, 50L))
  expect_equal(dim(fit$backend$magnitude_samples), c(n, 50L))
  expect_true(sd(as.numeric(fitted(fit, type = "magnitude")$magnitude)) > 0)
  expect_true(all(is.finite(mc_block_coef(fit))))
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)

  fit2 <- magComp(
    dat,
    magnitude = ~ z,
    composition = ~ z,
    control = magCompControl(
      n_samples = 30,
      n_burn = 10,
      n_thin = 1,
      n_chains = 2,
      chain_seeds = c(101, 202),
      nb_size = 1000,
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )

  expect_equal(fit2$backend$n_chains, 2L)
  expect_equal(dim(fit2$backend$magnitude_coefficient_samples), c(2L, 20L, 2L))
  expect_equal(dim(mc_block_coef(fit2)),
               c(length(cls) - 1L, 2L, 20L, 2L))
  expect_equal(dim(fit2$backend$proportion_samples),
               c(n, length(cls), 20L, 2L))
  expect_equal(dim(fitted(fit2, type = "count", summary = "samples")$count),
               c(n, length(cls), 20L, 2L))
  expect_equal(dim(fitted(fit2, type = "count", summary = "quantile")$count),
               c(n, length(cls), 3L))
  expect_equal(dim(fitted(fit2, type = "density", summary = "quantile",
                          probs = c(0.1, 0.9))$density),
               c(n, length(cls), 2L))
  expect_equal(dim(fitted(fit2, type = "proportion", summary = "samples")$proportion),
               c(n, length(cls), 20L, 2L))
  expect_equal(fitted(fit2, type = "proportion")$proportion,
               mc_regression_proportion_mean(fit2))
  expect_equal(dim(fitted(fit2, type = "magnitude", summary = "quantile")$magnitude),
               c(n, 3L))
  expect_equal(colnames(fitted(fit2, type = "magnitude",
                               summary = "quantile",
                               probs = c(0.1, 0.9))$magnitude),
               c("10%", "90%"))

  set.seed(71)
  y_rep <- fitted(fit2, target = "observation", type = "count", summary = "samples")$count
  expect_equal(dim(y_rep), c(n, length(cls), 20L, 2L))
  set.seed(71)
  ppc_unit <- ppcSummary(fit2, type = "unit_total")
  set.seed(71)
  ppc_category <- ppcSummary(fit2, type = "category_total")
  set.seed(71)
  ppc_cell <- ppcSummary(fit2, type = "unit_category",
                         probs = c(0.1, 0.9))
  set.seed(71)
  ppc_group <- ppcSummary(fit2, type = "group_total",
                          group = stand)
  stand_vec <- stand
  set.seed(71)
  ppc_group_category <- ppcSummary(fit2, type = "group_category",
                                   group = stand_vec)
  expect_s3_class(ppc_unit, "magCompPPCSummary")
  expect_equal(nrow(ppc_unit), n)
  expect_equal(nrow(ppc_category), length(cls))
  expect_equal(nrow(ppc_cell), n * length(cls))
  expect_equal(nrow(ppc_group), length(unique(stand)))
  expect_equal(nrow(ppc_group_category), length(unique(stand)) * length(cls))
  expect_true(all(c("group", "p_lower", "p_upper", "p_value") %in%
                    names(ppc_group)))
  expect_true(all(c("p_lower", "p_upper", "p_value") %in% names(ppc_unit)))
  expect_true(all(ppc_unit$n_draws == 40L))
  expect_equal(sum(ppc_group$observed), sum(dat$counts))
  expect_equal(sum(ppc_group_category$observed), sum(dat$counts))
  expect_equal(ppc_group_category$group,
               rep(unique(stand), times = length(cls)))
  expect_equal(ppc_group_category$category,
               rep(colnames(fit2$data$counts), each = length(unique(stand))))
  expect_error(ppcSummary(fit2, type = "group_total"), "group")

  tab <- summary(fit2)
  expect_s3_class(tab, "magCompSummary")
  expect_true("ess" %in% names(tab))
  expect_true(any(is.finite(tab$rhat)))
  expect_true(all(is.na(tab$rhat) | is.finite(tab$rhat)))
  expect_true(all(is.na(tab$rhat) | tab$rhat >= 1))
  expect_true(all(is.na(tab$ess) | tab$ess > 0))
  expect_equal(unique(tab$n_chains), 2L)

  trace_file <- tempfile(fileext = ".pdf")
  devices_before <- grDevices::dev.list()
  expect_equal(
    normalizePath(trace_file, mustWork = FALSE),
    plot(fit2, parameters = "(Intercept)", file = trace_file,
         n_col = 1, n_row = 2)
  )
  expect_identical(grDevices::dev.list(), devices_before)
  expect_true(file.exists(trace_file))
  expect_gt(file.info(trace_file)$size, 0)
  expect_error(plot(fit2, parameters = "unknown_parameter"), "no monitored")
  unlink(trace_file)
})

test_that("sampler fits the one-category model with smoothed coefficients", {
  # Small count table for fitting and output checks.
  n <- 12
  cls <- c("s", "m", "l", "xl")
  area <- rep(0.08, n)
  z <- seq(-1, 1, length.out = n)
  rows <- expand.grid(plot = paste0("p", seq_len(n)), cls = cls)
  rows$cls <- factor(rows$cls, levels = sort(cls))
  i <- as.integer(rows$plot)
  rows$count <- 1L + (i + 2L * as.integer(rows$cls)) %% 7L
  rows$area <- area[i]
  rows$z <- z[i]
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     ordered = TRUE, design = fixedArea(area))

  fit <- magComp(
    dat,
    magnitude = ~ z,
    composition = ~ z,
    control = magCompControl(
      n_samples = 18,
      n_burn = 6,
      n_thin = 1,
      n_chains = 2,
      chain_seeds = c(303, 404),
      nb_size = 1000,
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "rw1", slope = "rw2", loading = "none")
  )

  expect_true(fit$fitted)
  expect_equal(fit$backend$n_chains, 2L)
  expect_equal(fit$backend$nngp_graph_count, 0L)
  expect_equal(dim(fit$backend$magnitude_coefficient_samples), c(2L, 12L, 2L))
  expect_equal(dim(mc_block_coef(fit)),
               c(length(cls) - 1L, 2L, 12L, 2L))
  expect_equal(dim(fit$backend$proportion_samples),
               c(n, length(cls), 12L, 2L))
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)
  expect_true(all(is.finite(fitted(fit, type = "count")$count)))

  set.seed(43)
  y_rep <- fitted(fit, target = "observation", type = "count", summary = "samples")$count
  expect_equal(dim(y_rep), c(n, length(cls), 12L, 2L))

  newdata <- data.frame(plot = c("new1", "new2"),
                        stand = c("A", "B"),
                        area = c(0.08, 0.10),
                        z = c(-0.5, 0.5))
  pred_p <- predict(fit, newdata, type = "proportion")$proportion
  pred_count_samples <- predict(fit, newdata, type = "count",
                                summary = "samples", area = "area")$count
  pred_mag_q <- predict(fit, newdata, type = "magnitude",
                        summary = "quantile")$magnitude
  expect_equal(dim(pred_p), c(2L, length(cls)))
  expect_equal(rownames(pred_p), newdata$plot)
  expect_equal(unname(rowSums(pred_p)), rep(1, 2), tolerance = 1e-10)
  expect_equal(dim(pred_count_samples), c(2L, length(cls), 12L, 2L))
  expect_equal(dim(pred_mag_q), c(2L, 3L))
  expect_true(all(is.finite(pred_count_samples)))
  stand_density <- standSummary(fitted(fit, summary = "samples"),
    data = data.frame(id = fit$data$unit_levels, stand = rep(c("A", "B"), length.out = n), ha = 1),
    unit = id, group = stand, represented_area = ha,
    type = "density", summary = "samples")$density
  prediction <- predict(fit, newdata, summary = "samples")
  pred_stand_p <- standSummary(prediction, data = newdata, unit = plot, group = stand,
    type = "proportion", represented_area = area)$proportion
  expect_equal(dim(stand_density), c(2L, length(cls), 12L, 2L))
  expect_equal(dim(pred_stand_p), c(2L, length(cls)))
  expect_equal(unname(rowSums(pred_stand_p)), rep(1, 2), tolerance = 1e-10)
  expect_error(
    standSummary(prediction, data = newdata, unit = plot, group = stand,
               represented_area = area, area = "area"),
    "unused standSummary argument"
  )

})

test_that("sampler supports one-block composition factor_gp terms", {
  set.seed(437)

  n <- 20
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  x <- runif(n)
  ycoord <- runif(n)
  z <- as.numeric(scale(runif(n, -1, 1)))
  w <- 0.5 * sin(2 * pi * x)
  alpha <- c(-0.15, -0.70)
  lambda <- c(0.6, 1.0)
  y <- matrix(0L, n, length(cls))
  for(i in seq_len(n)){
    total_i <- max(1L, stats::rpois(1, area[i] * exp(log(140))))
    pi_i <- stats::plogis(alpha + lambda * z[i] * w[i])
    p_i <- c(pi_i[1L], (1 - pi_i[1L]) * pi_i[2L],
             prod(1 - pi_i))
    y[i,] <- as.vector(stats::rmultinom(1, total_i, p_i))
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, y[i,]),
      area = area[i],
      x = x[i],
      ycoord = ycoord[i],
      z = z[i]
    ))
  }

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit <- magComp(dat, composition = ~z + z:factor_gp(x, ycoord, K = 1, anchors = "m", factor_names = "mid"),
    priors = magCompPriors(composition = list(blocks = list(cls = list(terms = list(factor_gp_1 = list(covariance = list(phi = uniform(0.2,
        8)))))))), control = magCompControl(n_samples = 13, n_burn = 6, n_thin = 1, n_chains = 2,
        chain_seeds = c(1515, 1616), metropolis = list(batch_length = 3), nb_size = 1000, verbose = FALSE),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none"), tuning = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_gp_1 = list(covariance = list(phi = 0.15))))),
        names(dat$allocation_blocks)[1L]))))

  term <- fit$composition$process_terms[[1L]]
  expect_equal(term$term_type, "factor_gp")
  expect_equal(term$covariate_name, "z")
  expect_equal(dim(mc_block_loading(fit)),
               c(length(cls) - 1L, 1L, 7L, 2L))
  expect_equal(dim(mc_block_factor(fit)), c(n, 1L, 7L, 2L))
  expect_equal(dim(mc_block_factor_phi(fit)), c(1L, 7L, 2L))
  expect_loading_constraints(mc_block_loading(fit), "m", "mid")
  expect_equal(names(fit$backend$composition_block_covariance_samples),
               "cls")
  expect_equal(dim(fit$backend$composition_block_covariance_samples$cls),
               c(1L, 7L, 2L))
  expect_true(all(is.finite(mc_block_factor(fit))))
  expect_true(all(mc_block_factor_phi(fit) > term$covariance_lower[1,1]))
  expect_true(all(mc_block_factor_phi(fit) < term$covariance_upper[1,1]))
  expect_gt(mean(abs(fitted(fit, type = "proportion")$proportion -
                     mc_regression_proportion_mean(fit))), 0)
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)

  tab <- summary(fit)
  expect_true(any(tab$component == "composition" & tab$block == "cls" & tab$term == "factor_gp_1" & tab$parameter == "phi"))
})

test_that("sampler supports one-block composition factor_nngp terms", {
  set.seed(439)

  n <- 22
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  x <- runif(n)
  ycoord <- runif(n)
  alpha <- c(-0.20, -0.70)
  lambda <- c(0.5, 1.0)
  w <- 0.45 * sin(2 * pi * x)
  y <- matrix(0L, n, length(cls))
  for(i in seq_len(n)){
    total_i <- max(1L, stats::rpois(1, area[i] * 135))
    pi_i <- stats::plogis(alpha + lambda * w[i])
    p_i <- c(pi_i[1L], (1 - pi_i[1L]) * pi_i[2L],
             prod(1 - pi_i))
    y[i,] <- as.vector(stats::rmultinom(1, total_i, p_i))
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, y[i,]),
      area = area[i],
      x = x[i],
      ycoord = ycoord[i]
    ))
  }

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit <- magComp(dat, magnitude = ~nngp(x, ycoord, m = 4), composition = ~factor_nngp(x, ycoord, K = 1, m = 4,
    anchors = "m", factor_names = "mid"), control = magCompControl(n_samples = 13, n_burn = 6, n_thin = 1,
    n_chains = 2, chain_seeds = c(1919, 2020), metropolis = list(batch_length = 3), nb_size = 1000,
    verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none"),
    starting = list(magnitude = list(terms = list(nngp_1 = list(covariance = list(phi = fixed(2)))))),
    tuning = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_nngp_1 = list(covariance = list(phi = 0.15))))),
        names(dat$allocation_blocks)[1L]))))

  term <- fit$composition$process_terms[[1L]]
  expect_equal(term$term_type, "factor_nngp")
  expect_equal(dim(term$neighbor_index), c(n, 4L))
  expect_equal(fit$backend$nngp_graph_count, 1L)
  expect_loading_constraints(mc_block_loading(fit), "m", "mid")
  expect_equal(dim(mc_block_factor(fit)), c(n, 1L, 7L, 2L))
  expect_equal(dim(mc_block_factor_phi(fit)), c(1L, 7L, 2L))
  expect_equal(names(fit$backend$composition_block_covariance_samples),
               "cls")
  expect_true(all(is.finite(mc_block_factor(fit))))
  expect_true(all(mc_block_factor_phi(fit) > term$covariance_lower[1,1]))
  expect_true(all(mc_block_factor_phi(fit) < term$covariance_upper[1,1]))
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)

  tab <- summary(fit)
  expect_true(any(tab$component == "composition" & tab$block == "cls" & tab$term == "factor_nngp_1" &
                    tab$parameter == "phi"))
})

test_that("sampler fits the two-category conditional model", {
  set.seed(44)

  n <- 24
  spp <- c("hw", "sw")
  cls <- c("small", "medium", "large")
  area <- rep(0.1, n)
  z <- as.numeric(scale(runif(n, -1, 1)))
  gamma <- c(log(120), 0.25)
  alpha_a <- c(-0.10)
  beta_a <- c(0.40)
  alpha_b <- c(-0.20, -0.70)
  beta_b <- c(-0.25, 0.30)

  y <- array(0L, dim = c(n, length(spp), length(cls)),
             dimnames = list(NULL, spp, cls))
  for(i in seq_len(n)){
    total_i <- max(1L, rpois(1, area[i] *
                               exp(gamma[1] + gamma[2] * z[i])))
    pi_a <- plogis(alpha_a + beta_a * z[i])
    p_a <- c(pi_a, 1 - pi_a)
    y_a <- as.vector(rmultinom(1, total_i, p_a))
    for(a in seq_along(spp)){
      if(y_a[a] == 0)
        next
      pi_b <- plogis(alpha_b + beta_b * z[i])
      p_b <- numeric(length(cls))
      remaining <- 1
      for(j in seq_along(alpha_b)){
        p_b[j] <- remaining * pi_b[j]
        remaining <- remaining * (1 - pi_b[j])
      }
      p_b[length(cls)] <- remaining
      y[i,a,] <- as.vector(rmultinom(1, y_a[a], p_b))
    }
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    for(a in seq_along(spp)){
      for(j in seq_along(cls)){
        if(y[i,a,j] == 0)
          next
        trees <- rbind(trees, data.frame(
          plot = paste0("p", i),
          species_group = rep(spp[a], y[i,a,j]),
          dbh_class = rep(cls[j], y[i,a,j]),
          area = area[i],
          z = z[i]
        ))
      }
    }
  }

  dat <- magCompData(
    trees,
    unit = plot,
    categories = c(species_group, dbh_class),
    design = fixedArea(area)
  )

  fit <- magComp(
    dat,
    magnitude = ~ z,
    composition = magCompComposition(
      species_group ~ z,
      dbh_class | species_group ~ z
    ),
    control = magCompControl(
      n_samples = 15,
      n_burn = 5,
      n_thin = 1,
      n_chains = 2,
      chain_seeds = c(505, 606),
      nb_size = 1000,
      verbose = FALSE
    )
  )

  expect_true(fit$fitted)
  expect_equal(dim(fit$backend$proportion_samples),
               c(n, length(spp) * length(cls), 10L, 2L))
  expect_equal(dim(fit$backend$pi_samples),
               c(n, (length(spp) - 1L) +
                   length(spp) * (length(cls) - 1L), 10L, 2L))
  expect_equal(names(fit$backend$composition_block_coefficient_samples),
               names(dat$allocation_blocks))
  expect_equal(dim(fit$backend$composition_block_coefficient_samples$species_group),
               c(length(spp) - 1L, 2L, 10L, 2L))
  expect_equal(dim(fit$backend$composition_block_coefficient_samples[["dbh_class|species_group"]]),
               c(length(spp) * (length(cls) - 1L), 2L, 10L, 2L))
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)
  expect_true(all(is.finite(fitted(fit, type = "count")$count)))
  set.seed(44)
  expect_equal(dim(fitted(fit, target = "observation", type = "count", summary = "samples")$count),
               c(n, length(spp) * length(cls), 10L, 2L))

  newdata <- data.frame(plot = c("hold1", "hold2"),
                        stand = c("s1", "s1"),
                        area = c(0.1, 0.1),
                        z = c(-0.3, 0.4))
  pred_p <- predict(fit, newdata, type = "proportion")$proportion
  pred_count_q <- predict(fit, newdata, type = "count",
                          summary = "quantile", area = "area")$count
  expect_equal(dim(pred_p), c(2L, length(spp) * length(cls)))
  expect_equal(colnames(pred_p), colnames(dat$counts))
  expect_equal(unname(rowSums(pred_p)), rep(1, 2), tolerance = 1e-10)
  expect_equal(dim(pred_count_q),
               c(2L, length(spp) * length(cls), 3L))
  expect_true(all(is.finite(pred_count_q)))
  pred_stand_count <- standSummary(predict(fit, newdata, summary = "samples"),
    data = newdata, unit = plot, group = stand, represented_area = area,
    type = "total", summary = "quantile")$total
  expect_equal(dim(pred_stand_count),
               c(1L, length(spp) * length(cls), 3L))

  tab <- summary(fit)
  expect_true(any(grepl("dbh_class\\|species_group", tab$block)))
})

test_that("sampler supports conditional composition factor_gp terms", {
  set.seed(438)

  n <- 18
  spp <- c("hw", "sw")
  cls <- c("small", "medium", "large")
  area <- rep(0.1, n)
  x <- runif(n)
  ycoord <- runif(n)
  z <- as.numeric(scale(runif(n, -1, 1)))
  alpha_a <- c(-0.2)
  alpha_b <- c(-0.25, -0.75, -0.05, -0.55)
  lambda_b <- c(0.4, -0.2, 0.3, -0.1)
  w <- 0.4 * cos(2 * pi * x)

  y <- array(0L, dim = c(n, length(spp), length(cls)))
  for(i in seq_len(n)){
    total_i <- max(1L, stats::rpois(1, area[i] * 130))
    pi_a <- stats::plogis(alpha_a)
    p_a <- c(pi_a, 1 - pi_a)
    y_a <- as.vector(stats::rmultinom(1, total_i, p_a))
    for(a in seq_along(spp)){
      if(y_a[a] == 0L)
        next
      eta_b <- alpha_b[((a - 1L) * 2L + 1L):((a - 1L) * 2L + 2L)] +
        lambda_b[((a - 1L) * 2L + 1L):((a - 1L) * 2L + 2L)] * w[i]
      pi_b <- stats::plogis(eta_b)
      p_b <- c(pi_b[1L], (1 - pi_b[1L]) * pi_b[2L],
               prod(1 - pi_b))
      y[i,a,] <- as.vector(stats::rmultinom(1, y_a[a], p_b))
    }
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    for(a in seq_along(spp)){
      for(j in seq_along(cls)){
        if(y[i,a,j] == 0)
          next
        trees <- rbind(trees, data.frame(
          plot = paste0("p", i),
          species_group = rep(spp[a], y[i,a,j]),
          dbh_class = rep(cls[j], y[i,a,j]),
          area = area[i],
          x = x[i],
          ycoord = ycoord[i],
          z = z[i]
        ))
      }
    }
  }

  dat <- magCompData(
    trees,
    unit = plot,
    categories = c(species_group, dbh_class),
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit <- magComp(dat, magnitude = ~1, composition = magCompComposition(species_group ~ z, dbh_class | species_group ~
    z + factor_gp(x, ycoord, K = 1, anchors = "hw:medium", factor_names = "shape")), control = magCompControl(n_samples = 12,
    n_burn = 6, n_thin = 1, n_chains = 2, chain_seeds = c(1717, 1818), metropolis = list(batch_length = 3),
    nb_size = 1000, verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none",
    loading = "none"), tuning = list(composition = list(blocks = list(`dbh_class|species_group` = list(terms = list(factor_gp_1 = list(covariance = list(phi = 0.15))))))))

  block_name <- "dbh_class|species_group"
  expect_equal(names(fit$backend$composition_block_covariance_samples),
               block_name)
  expect_equal(dim(fit$backend$composition_block_loading_samples[[block_name]]),
               c(length(spp) * (length(cls) - 1L), 1L, 6L, 2L))
  expect_loading_constraints(
    fit$backend$composition_block_loading_samples[[block_name]],
    "hw:medium", "shape"
  )
  expect_equal(dim(fit$backend$composition_block_score_samples[[block_name]]),
               c(n, 1L, 6L, 2L))
  expect_equal(dim(fit$backend$composition_block_covariance_samples[[block_name]]),
               c(1L, 6L, 2L))
  expect_equal(dim(fit$backend$composition_block_covariance_acceptance),
               c(1L, 2L))
  expect_equal(rownames(fit$backend$composition_block_covariance_acceptance),
               block_name)
  expect_true(all(is.finite(
    fit$backend$composition_block_covariance_acceptance[block_name,]
  )))
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)

  tab <- summary(fit)
  expect_true(any(tab$component == "composition" & tab$block == block_name & tab$term == "factor_gp_1" &
                    tab$parameter == "phi"))

  spatial_diag <- spatialDiagnostics(fit)
  expect_true(any(spatial_diag$summary$block == block_name & spatial_diag$summary$term == "factor_gp_1" &
                    spatial_diag$summary$quantity == "phi"))
  expect_equal(nrow(spatial_diag$acceptance), 2L)

  newdata <- data.frame(plot = paste0("new", seq_len(3L)),
                        x = runif(3L),
                        ycoord = runif(3L),
                        z = runif(3L, -1, 1),
                        area = area[1L])
  set.seed(321)
  pred_full <- predict(fit, newdata, type = "proportion",
                       summary = "samples")$proportion
  expect_equal(dim(pred_full), c(3L, length(spp) * length(cls), 6L, 2L))
  expect_equal(as.numeric(apply(pred_full, c(1L, 3L, 4L), sum)),
               rep(1, 3L * 6L * 2L),
               tolerance = 1e-10)

  start <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(names(lapply(start$composition$blocks, "[[", "coefficients")), c("species_group",
                                                      block_name))
  expect_equal(dim(start$composition$blocks[[block_name]]$coefficients),
               c(length(spp) * (length(cls) - 1L), 2L))
  expect_equal(names(lapply(Filter(function(b) length(b$terms) > 0L, start$composition$blocks), function(b) lapply(b$terms, "[[", "loadings"))), block_name)
  expect_equal(dim(start$composition$blocks[[block_name]]$terms$factor_gp_1$loadings),
               c(length(spp) * (length(cls) - 1L), 1L))
  expect_loading_constraints(
    start$composition$blocks[[block_name]]$terms$factor_gp_1$loadings,
    "hw:medium", "shape"
  )
  expect_equal(dim(start$composition$blocks[[block_name]]$terms$factor_gp_1$scores),
               c(n, 1L))
  expect_equal(
    length(start$composition$blocks[[block_name]]$terms$factor_gp_1$covariance$phi),
    1L
  )

  fit_started <- magComp(dat, magnitude = ~1, composition = magCompComposition(species_group ~ z, dbh_class | species_group ~
    z + factor_gp(x, ycoord, K = 1, anchors = "hw:medium", factor_names = "shape")), starting = start,
    control = magCompControl(n_samples = 3, n_burn = 1, n_thin = 1, n_chains = 1, chain_seeds = 1919,
        nb_size = 1000, verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none",
        loading = "none"), tuning = list(composition = list(blocks = list(`dbh_class|species_group` = list(terms = list(factor_gp_1 = list(covariance = list(phi = 0.15))))))))
  expect_true(fit_started$fitted)
  expect_loading_constraints(
    fit_started$backend$composition_block_loading_samples[[block_name]],
    "hw:medium", "shape"
  )
})

test_that("sampler supports conditional composition factor_nngp terms", {
  set.seed(440)

  n <- 18
  spp <- c("hw", "sw")
  cls <- c("small", "medium", "large")
  area <- rep(0.1, n)
  x <- runif(n)
  ycoord <- runif(n)
  z <- as.numeric(scale(runif(n, -1, 1)))
  alpha_a <- c(-0.15)
  alpha_b <- c(-0.30, -0.80, -0.10, -0.55)
  lambda_b <- c(0.5, -0.2, 0.25, -0.15)
  w <- 0.35 * sin(2 * pi * x)

  y <- array(0L, dim = c(n, length(spp), length(cls)))
  for(i in seq_len(n)){
    total_i <- max(1L, stats::rpois(1, area[i] * 130))
    pi_a <- stats::plogis(alpha_a)
    p_a <- c(pi_a, 1 - pi_a)
    y_a <- as.vector(stats::rmultinom(1, total_i, p_a))
    for(a in seq_along(spp)){
      if(y_a[a] == 0L)
        next
      eta_b <- alpha_b[((a - 1L) * 2L + 1L):((a - 1L) * 2L + 2L)] +
        lambda_b[((a - 1L) * 2L + 1L):((a - 1L) * 2L + 2L)] * w[i]
      pi_b <- stats::plogis(eta_b)
      p_b <- c(pi_b[1L], (1 - pi_b[1L]) * pi_b[2L],
               prod(1 - pi_b))
      y[i,a,] <- as.vector(stats::rmultinom(1, y_a[a], p_b))
    }
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    for(a in seq_along(spp)){
      for(j in seq_along(cls)){
        if(y[i,a,j] == 0)
          next
        trees <- rbind(trees, data.frame(
          plot = paste0("p", i),
          species_group = rep(spp[a], y[i,a,j]),
          dbh_class = rep(cls[j], y[i,a,j]),
          area = area[i],
          x = x[i],
          ycoord = ycoord[i],
          z = z[i]
        ))
      }
    }
  }

  dat <- magCompData(
    trees,
    unit = plot,
    categories = c(species_group, dbh_class),
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit <- magComp(dat, magnitude = ~1, composition = magCompComposition(species_group ~ z, dbh_class | species_group ~
    z + factor_nngp(x, ycoord, K = 1, m = 4, anchors = "hw:medium", factor_names = "shape")), control = magCompControl(n_samples = 12,
    n_burn = 6, n_thin = 1, n_chains = 2, chain_seeds = c(2121, 2222), metropolis = list(batch_length = 3),
    nb_size = 1000, verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none",
    loading = "none"), tuning = list(composition = list(blocks = list(`dbh_class|species_group` = list(terms = list(factor_nngp_1 = list(covariance = list(phi = 0.15))))))))

  block_name <- "dbh_class|species_group"
  term <- fit$composition$allocation_components[[block_name]]$process_terms[[1L]]
  expect_equal(term$term_type, "factor_nngp")
  expect_equal(dim(term$neighbor_index), c(n, 4L))
  expect_equal(fit$backend$nngp_graph_count, 1L)
  expect_equal(dim(fit$backend$composition_block_loading_samples[[block_name]]),
               c(length(spp) * (length(cls) - 1L), 1L, 6L, 2L))
  expect_loading_constraints(
    fit$backend$composition_block_loading_samples[[block_name]],
    "hw:medium", "shape"
  )
  expect_equal(dim(fit$backend$composition_block_score_samples[[block_name]]),
               c(n, 1L, 6L, 2L))
  expect_equal(dim(fit$backend$composition_block_covariance_samples[[block_name]]),
               c(1L, 6L, 2L))
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)

  tab <- summary(fit)
  expect_true(any(tab$component == "composition" & tab$block == block_name & tab$term == "factor_nngp_1" &
                    tab$parameter == "phi"))
})

test_that("sampler supports conditional composition factor_iid terms", {
  set.seed(445)

  n <- 20
  spp <- c("hw", "sw")
  cls <- c("small", "medium", "large")
  area <- rep(0.1, n)
  y <- array(0L, dim = c(n, length(spp), length(cls)))
  for(i in seq_len(n)){
    total_i <- max(1L, stats::rpois(1, area[i] * 110))
    y_a <- as.vector(stats::rmultinom(1, total_i, c(0.55, 0.45)))
    for(a in seq_along(spp)){
      if(y_a[a] > 0){
        y[i,a,] <- as.vector(stats::rmultinom(1, y_a[a],
                                              c(0.50, 0.30, 0.20)))
      }
    }
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    for(a in seq_along(spp)){
      for(j in seq_along(cls)){
        if(y[i,a,j] == 0)
          next
        trees <- rbind(trees, data.frame(
          plot = paste0("p", i),
          species_group = rep(spp[a], y[i,a,j]),
          dbh_class = rep(cls[j], y[i,a,j]),
          area = area[i]
        ))
      }
    }
  }

  dat <- magCompData(
    trees,
    unit = plot,
    categories = c(species_group, dbh_class),
    design = fixedArea(area)
  )

  fit <- magComp(
    dat,
    magnitude = ~ 1,
    composition = magCompComposition(
      species_group ~ 1,
      dbh_class | species_group ~ factor_iid(plot, K = 1)
    ),
    control = magCompControl(
      n_samples = 12,
      n_burn = 4,
      n_thin = 1,
      n_chains = 2,
      chain_seeds = c(1111, 1212),
      nb_size = 1000,
      verbose = FALSE
    )
  )

  block_name <- "dbh_class|species_group"
  expect_equal(names(fit$backend$composition_block_loading_samples),
               block_name)
  expect_equal(dim(fit$backend$composition_block_loading_samples[[block_name]]),
               c(length(spp) * (length(cls) - 1L), 1L, 8L, 2L))
  term <- fit$composition$allocation_components[[block_name]]$process_terms[[1L]]
  expect_loading_constraints(
    fit$backend$composition_block_loading_samples[[block_name]],
    term$anchors,
    term$factor_names
  )
  expect_equal(dim(fit$backend$composition_block_score_samples[[block_name]]),
               c(n, 1L, 8L, 2L))
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)
  expect_true(all(is.finite(fitted(fit, type = "count")$count)))

  tab <- summary(fit)
  expect_true(any(grepl("dbh_class\\|species_group", tab$block)))
})

test_that("sampler supports magnitude iid random intercepts", {
  # Small count table for fitting and output checks.
  n <- 12
  cls <- c("s", "m", "l", "xl")
  area <- rep(0.08, n)
  group <- rep(paste0("g", 1:6), length.out = n)
  z <- seq(-1, 1, length.out = n)
  rows <- expand.grid(plot = paste0("p", seq_len(n)), cls = cls)
  rows$cls <- factor(rows$cls, levels = sort(cls))
  i <- as.integer(rows$plot)
  rows$count <- 1L + (i + 2L * as.integer(rows$cls)) %% 7L
  rows$area <- area[i]
  rows$group <- group[i]
  rows$z <- z[i]
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     ordered = TRUE, design = fixedArea(area))

  fit <- magComp(
    dat,
    magnitude = ~ z + iid(group),
    control = magCompControl(
      n_samples = 30,
      n_burn = 10,
      n_thin = 1,
      n_chains = 2,
      chain_seeds = c(505, 606),
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )

  expect_equal(names(fit$backend$magnitude_process_effect_samples), "iid_1")
  expect_equal(dim(mc_mag_process(fit, "iid_1")), c(6L, 20L, 2L))
  expect_equal(dim(fitted(fit, type = "magnitude", summary = "samples")$magnitude),
               c(n, 20L, 2L))
  expect_equal(dim(mc_regression_magnitude_samples(fit)),
               c(n, 20L, 2L))
  expect_equal(dim(fitted(fit, type = "count")$count),
               c(n, length(cls)))
  expect_gt(mean(abs(fitted(fit, type = "magnitude")$magnitude -
                     mc_regression_magnitude_mean(fit))), 0)

  start <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(length(start$magnitude$terms$iid_1$effects), 6L)
})

test_that("sampler supports magnitude iid varying slopes", {
  # Small count table for fitting and output checks.
  n <- 12
  cls <- c("s", "m", "l", "xl")
  area <- rep(0.08, n)
  group <- rep(paste0("g", 1:4), length.out = n)
  z <- seq(-1, 1, length.out = n)
  rows <- expand.grid(plot = paste0("p", seq_len(n)), cls = cls)
  rows$cls <- factor(rows$cls, levels = sort(cls))
  i <- as.integer(rows$plot)
  rows$count <- 1L + (i + 2L * as.integer(rows$cls)) %% 7L
  rows$area <- area[i]
  rows$group <- group[i]
  rows$z <- z[i]
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     ordered = TRUE, design = fixedArea(area))

  fit <- magComp(
    dat,
    magnitude = ~ z + z:iid(group),
    control = magCompControl(
      n_samples = 18,
      n_burn = 6,
      n_thin = 1,
      n_chains = 2,
      chain_seeds = c(515, 616),
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )

  term <- fit$magnitude$process_terms[[1L]]
  expect_equal(term$term_type, "iid")
  expect_equal(term$covariate_name, "z")
  expect_equal(names(fit$backend$magnitude_process_effect_samples), "iid_1")
  expect_equal(dim(mc_mag_process(fit, "iid_1")), c(4L, 12L, 2L))
  expect_equal(dim(fitted(fit, type = "magnitude", summary = "samples")$magnitude),
               c(n, 12L, 2L))
  expect_gt(mean(abs(fitted(fit, type = "magnitude")$magnitude -
                     mc_regression_magnitude_mean(fit))), 0)

  start <- magCompStarting(fit, use = "last", chain = 1)
  expect_equal(length(start$magnitude$terms$iid_1$effects), 4L)
})

test_that("magnitude latent defaults start from nonspatial residuals", {
  set.seed(465)

  n <- 24
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  group <- rep(paste0("g", 1:4), length.out = n)
  x <- runif(n)
  ycoord <- runif(n)
  z <- as.numeric(scale(runif(n, -1, 1)))
  counts <- matrix(0L, n, length(cls))
  for(i in seq_len(n)){
    total_i <- max(1L, rpois(1, area[i] *
                               exp(log(130) + 0.4 * z[i] +
                                     0.2 * (group[i] %in% c("g2", "g4")))))
    counts[i,] <- as.vector(rmultinom(1, total_i, c(0.45, 0.35, 0.20)))
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, counts[i,]),
      area = area[i],
      group = group[i],
      x = x[i],
      ycoord = ycoord[i],
      z = z[i]
    ))
  }

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit_iid <- magComp(
    dat,
    magnitude = ~ z + z:iid(group),
    control = magCompControl(run_sampler = FALSE)
  )
  expect_equal(length(fit_iid$starting$magnitude$coefficients), 2L)
  expect_equal(length(fit_iid$starting$magnitude$terms$iid_1$effects), 4L)
  expect_true(all(is.finite(fit_iid$starting$magnitude$terms$iid_1$effects)))

  supplied_start <- list(magnitude = list(coefficients = c(1, 2), terms = list(iid_1 = list(effects = stats::setNames(rep(0.25, 4),
                                                   paste0("g", 1:4))))))
  fit_supplied <- magComp(
    dat,
    magnitude = ~ z + z:iid(group),
    starting = supplied_start,
    control = magCompControl(run_sampler = FALSE)
  )
  expect_equal(unname(fit_supplied$starting$magnitude$coefficients), c(1, 2))
  expect_equal(unname(fit_supplied$starting$magnitude$terms$iid_1$effects),
               rep(0.25, 4))

  fit_nngp <- magComp(dat, magnitude = ~z + nngp(x, ycoord, m = 4), control = magCompControl(run_sampler = FALSE),
    starting = list(magnitude = list(terms = list(nngp_1 = list(covariance = list(phi = fixed(2)))))))
  expect_equal(length(fit_nngp$starting$magnitude$terms$nngp_1$effects), n)
  expect_true(is.finite(
    fit_nngp$starting$magnitude$terms$nngp_1$variance
  ))
  expect_gt(fit_nngp$starting$magnitude$terms$nngp_1$variance, 0)

  fit_gp <- magComp(dat, magnitude = ~z + gp(x, ycoord), control = magCompControl(run_sampler = FALSE), starting = list(magnitude = list(terms = list(gp_1 = list(covariance = list(phi = fixed(2)))))))
  expect_equal(length(fit_gp$starting$magnitude$terms$gp_1$effects), n)
  expect_true(is.finite(fit_gp$starting$magnitude$terms$gp_1$variance))
  expect_gt(fit_gp$starting$magnitude$terms$gp_1$variance, 0)
})

test_that("sampler supports magnitude GP varying slopes", {
  set.seed(466)

  n <- 18
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  x <- runif(n)
  ycoord <- runif(n)
  z <- as.numeric(scale(runif(n, -1, 1)))
  w <- 0.25 * sin(2 * pi * x)
  gamma <- c(log(125), 0.25)
  alpha <- c(-0.30, -0.80)

  counts <- matrix(0L, n, length(cls))
  for(i in seq_len(n)){
    total_i <- max(1L, stats::rpois(1, area[i] *
                                     exp(gamma[1L] + gamma[2L] * z[i] +
                                           z[i] * w[i])))
    pi_i <- stats::plogis(alpha)
    p_i <- c(pi_i[1L], (1 - pi_i[1L]) * pi_i[2L], prod(1 - pi_i))
    counts[i,] <- as.vector(stats::rmultinom(1, size = total_i, prob = p_i))
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, counts[i,]),
      area = area[i],
      x = x[i],
      ycoord = ycoord[i],
      z = z[i]
    ))
  }

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit <- magComp(dat, magnitude = ~z + z:gp(x, ycoord), priors = magCompPriors(magnitude = list(terms = list(gp_1 = list(covariance = list(phi = uniform(0.2,
    8)))))), control = magCompControl(n_samples = 13, n_burn = 6, n_thin = 1, n_chains = 2, chain_seeds = c(1313,
    1414), verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none"),
    tuning = list(magnitude = list(terms = list(gp_1 = list(covariance = list(phi = 0.15))))))

  term <- fit$magnitude$process_terms[[1L]]
  expect_equal(term$term_type, "gp")
  expect_equal(term$covariate_name, "z")
  expect_equal(dim(mc_mag_process(fit, "gp_1")), c(n, 7L, 2L))
  expect_equal(dim(mc_mag_sigma(fit, "gp_1")),
               c(1L, 7L, 2L))
  expect_equal(dim(mc_mag_phi(fit, "gp_1")),
               c(1L, 7L, 2L))
  expect_equal(fit$backend$nngp_graph_count, 0L)
  expect_true(all(is.finite(mc_mag_process(fit, "gp_1"))))
  expect_true(all(mc_mag_sigma(fit, "gp_1") > 0))
  expect_true(all(mc_mag_phi(fit, "gp_1") >= 0.2))
  expect_true(all(mc_mag_phi(fit, "gp_1") <= 8))
  expect_equal(length(fit$backend$magnitude_process_covariance_acceptance), 2L)
  expect_true(all(is.finite(
    fit$backend$magnitude_process_covariance_acceptance["gp_1",]
  )))
  expect_gt(mean(abs(fitted(fit, type = "magnitude")$magnitude -
                     mc_regression_magnitude_mean(fit))), 0)

  tab <- summary(fit)
  expect_true(any(tab$component == "magnitude" & tab$term == "gp_1" &
                    tab$parameter == "variance"))
  expect_true(any(tab$component == "magnitude" & tab$term == "gp_1" &
                    tab$parameter == "phi"))
  start <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(length(start$magnitude$terms$gp_1$effects), n)
  expect_equal(length(start$magnitude$terms$gp_1$variance), 1L)
  expect_equal(length(start$magnitude$terms$gp_1$covariance$phi), 1L)
})

test_that("sampler supports magnitude NNGP random intercepts", {
  set.seed(47)

  n <- 24
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  x <- runif(n)
  ycoord <- runif(n)
  z <- as.numeric(scale(runif(n, -1, 1)))
  w <- rnorm(n, 0, 0.35)
  gamma <- c(log(130), 0.30)
  alpha <- c(-0.25, -0.75)

  counts <- matrix(0L, n, length(cls))
  for(i in seq_len(n)){
    total_i <- max(1L, rpois(1, area[i] *
                               exp(gamma[1L] + gamma[2L] * z[i] + w[i])))
    pi_i <- plogis(alpha)
    p_i <- c(pi_i[1L], (1 - pi_i[1L]) * pi_i[2L],
             prod(1 - pi_i))
    counts[i,] <- as.vector(rmultinom(1, size = total_i, prob = p_i))
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, counts[i,]),
      area = area[i],
      x = x[i],
      ycoord = ycoord[i],
      z = z[i]
    ))
  }

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit <- magComp(dat, magnitude = ~z + nngp(x, ycoord, m = 5), priors = magCompPriors(magnitude = list(terms = list(nngp_1 = list(covariance = list(phi = uniform(0.2,
    8)))))), control = magCompControl(n_samples = 16, n_burn = 6, n_thin = 1, n_chains = 2, chain_seeds = c(707,
    808), metropolis = list(batch_length = 3), verbose = FALSE), smooth = magCompSmooth(intercept = "none",
    slope = "none", loading = "none"), tuning = list(magnitude = list(terms = list(nngp_1 = list(covariance = list(phi = 0.15))))))

  expect_equal(dim(mc_mag_process(fit, "nngp_1")), c(n, 10L, 2L))
  expect_equal(dim(mc_mag_sigma(fit, "nngp_1")),
               c(1L, 10L, 2L))
  expect_equal(dim(mc_mag_phi(fit, "nngp_1")),
               c(1L, 10L, 2L))
  expect_equal(fit$backend$nngp_graph_count, 1L)
  expect_true(all(is.finite(mc_mag_process(fit, "nngp_1"))))
  expect_true(all(mc_mag_sigma(fit, "nngp_1") > 0))
  expect_true(all(mc_mag_phi(fit, "nngp_1") >= 0.2))
  expect_true(all(mc_mag_phi(fit, "nngp_1") <= 8))
  mag_accept <- fit$backend$magnitude_process_covariance_acceptance["nngp_1",]
  expect_equal(length(mag_accept), 2L)
  expect_true(all(is.finite(mag_accept)))
  expect_true(all(mag_accept >= 0))
  expect_true(all(mag_accept <= 1))
  mag_tuning <- fit$backend$magnitude_process_covariance_tuning["nngp_1",]
  expect_equal(length(mag_tuning), 2L)
  expect_true(all(is.finite(mag_tuning)))
  expect_equal(length(fit$backend$magnitude_process_covariance_batch_acceptance),
               1L)
  expect_equal(length(
    fit$backend$magnitude_process_covariance_batch_acceptance$nngp_1[[1L]]
  ),
               5L)
  expect_equal(length(
    fit$backend$magnitude_process_covariance_tuning_history$nngp_1[[1L]]
  ),
               5L)
  expect_gt(mean(abs(fitted(fit, type = "magnitude")$magnitude -
                     mc_regression_magnitude_mean(fit))), 0)

  tab <- summary(fit)
  expect_true(any(tab$component == "magnitude" & tab$term == "nngp_1" &
                    tab$parameter == "variance"))
  expect_true(any(tab$component == "magnitude" & tab$term == "nngp_1" &
                    tab$parameter == "phi"))

  spatial_diag <- spatialDiagnostics(fit)
  expect_s3_class(spatial_diag, "magCompSpatialDiagnostics")
  expect_true(any(spatial_diag$summary$component == "magnitude" &
                    spatial_diag$summary$quantity == "variance"))
  expect_true(any(spatial_diag$summary$component == "magnitude" &
                    spatial_diag$summary$quantity == "phi"))
  expect_equal(nrow(spatial_diag$acceptance), 2L)
  expect_true(any(spatial_diag$correlations$component == "magnitude" &
                    spatial_diag$correlations$x == "phi" &
                    spatial_diag$correlations$y == "variance"))

  newdata <- data.frame(plot = paste0("new", seq_len(4L)),
                        x = runif(4L),
                        ycoord = runif(4L),
                        z = runif(4L, -1, 1),
                        area = area[1L])
  pred_fixed <- mc_regression_magnitude_samples(fit, newdata)
  set.seed(123)
  pred_full <- predict(fit, newdata, type = "magnitude",
                       summary = "samples")$magnitude
  expect_equal(dim(pred_full), c(4L, 10L, 2L))
  expect_true(all(is.finite(pred_full)))
  expect_gt(mean(abs(pred_full - pred_fixed)), 0)

  start <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(length(start$magnitude$terms$nngp_1$effects), n)
  expect_equal(length(start$magnitude$terms$nngp_1$variance), 1L)
  expect_equal(length(start$magnitude$terms$nngp_1$covariance$phi), 1L)
})

test_that("sampler supports magnitude dense GP random intercepts", {
  set.seed(49)

  n <- 18
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  x <- runif(n)
  ycoord <- runif(n)
  z <- as.numeric(scale(runif(n, -1, 1)))
  w <- rnorm(n, 0, 0.3)
  gamma <- c(log(130), 0.25)
  alpha <- c(-0.25, -0.75)

  counts <- matrix(0L, n, length(cls))
  for(i in seq_len(n)){
    total_i <- max(1L, rpois(1, area[i] *
                               exp(gamma[1L] + gamma[2L] * z[i] + w[i])))
    pi_i <- plogis(alpha)
    p_i <- c(pi_i[1L], (1 - pi_i[1L]) * pi_i[2L],
             prod(1 - pi_i))
    counts[i,] <- as.vector(rmultinom(1, size = total_i, prob = p_i))
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, counts[i,]),
      area = area[i],
      x = x[i],
      ycoord = ycoord[i],
      z = z[i]
    ))
  }

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit <- magComp(dat, magnitude = ~z + gp(x, ycoord), priors = magCompPriors(magnitude = list(terms = list(gp_1 = list(covariance = list(phi = uniform(0.2,
    8)))))), control = magCompControl(n_samples = 13, n_burn = 5, n_thin = 1, n_chains = 2, chain_seeds = c(709,
    810), verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none"),
    tuning = list(magnitude = list(terms = list(gp_1 = list(covariance = list(phi = 0.15))))))

  term <- fit$magnitude$process_terms[[1L]]
  expect_equal(term$term_type, "gp")
  expect_null(term$neighbor_index)
  expect_equal(dim(mc_mag_process(fit, "gp_1")), c(n, 8L, 2L))
  expect_equal(dim(mc_mag_sigma(fit, "gp_1")),
               c(1L, 8L, 2L))
  expect_equal(dim(mc_mag_phi(fit, "gp_1")),
               c(1L, 8L, 2L))
  expect_true(all(is.finite(mc_mag_process(fit, "gp_1"))))
  expect_true(all(mc_mag_sigma(fit, "gp_1") > 0))
  expect_true(all(mc_mag_phi(fit, "gp_1") >= 0.2))
  expect_true(all(mc_mag_phi(fit, "gp_1") <= 8))

  tab <- summary(fit)
  expect_true(any(tab$component == "magnitude" & tab$term == "gp_1" &
                    tab$parameter == "variance"))
  expect_true(any(tab$component == "magnitude" & tab$term == "gp_1" &
                    tab$parameter == "phi"))

  spatial_diag <- spatialDiagnostics(fit)
  expect_true(any(spatial_diag$summary$component == "magnitude" &
                    spatial_diag$summary$term == "gp_1" &
                    spatial_diag$summary$quantity == "variance"))

  start <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(length(start$magnitude$terms$gp_1$effects), n)
  expect_equal(length(start$magnitude$terms$gp_1$variance), 1L)
  expect_equal(length(start$magnitude$terms$gp_1$covariance$phi), 1L)
})

test_that("sampler supports colon spatially varying coefficients", {
  # Small count table for fitting and output checks.
  n <- 12
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  x <- rep(0:3, 3)
  ycoord <- rep(0:2, each = 4)
  z <- seq(-1, 1, length.out = n)
  rows <- expand.grid(plot = paste0("p", seq_len(n)), cls = cls)
  rows$cls <- factor(rows$cls, levels = sort(cls))
  i <- as.integer(rows$plot)
  rows$count <- 1L + (i + 2L * as.integer(rows$cls)) %% 7L
  rows$area <- area[i]
  rows$x <- x[i]
  rows$ycoord <- ycoord[i]
  rows$z <- z[i]
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     design = fixedArea(area))

  set.seed(50)
  fit_nngp <- magComp(dat, magnitude = ~z + z:nngp(x, ycoord, m = 5), control = magCompControl(n_samples = 13,
    n_burn = 5, n_thin = 1, verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none",
    loading = "none"), starting = list(magnitude = list(terms = list(nngp_1 = list(covariance = list(phi = fixed(2)))))))

  term <- fit_nngp$magnitude$process_terms[[1L]]
  expect_equal(term$term_type, "nngp")
  expect_equal(term$covariate_name, "z")
  expect_equal(dim(mc_mag_process(fit_nngp, "nngp_1")), c(n, 8L))
  expect_true(all(is.finite(mc_mag_process(fit_nngp, "nngp_1"))))
  expect_gt(mean(abs(fitted(fit_nngp, type = "magnitude")$magnitude -
                       mc_regression_magnitude_mean(fit_nngp))), 0)
  tab <- summary(fit_nngp)
  expect_true(any(tab$component == "magnitude" & tab$term == "nngp_1" &
                    tab$parameter == "variance"))
  start <- magCompStarting(fit_nngp, use = "last")
  expect_equal(length(start$magnitude$terms$nngp_1$effects), n)

  fit_gp <- magComp(dat, magnitude = ~z + z:gp(x, ycoord), control = magCompControl(n_samples = 10, n_burn = 4,
    n_thin = 1, verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none"),
    starting = list(magnitude = list(terms = list(gp_1 = list(covariance = list(phi = fixed(2)))))))

  term <- fit_gp$magnitude$process_terms[[1L]]
  expect_equal(term$term_type, "gp")
  expect_equal(term$covariate_name, "z")
  expect_null(term$neighbor_index)
  expect_equal(dim(mc_mag_process(fit_gp, "gp_1")), c(n, 6L))
  expect_true(all(is.finite(mc_mag_process(fit_gp, "gp_1"))))
  tab <- summary(fit_gp)
  expect_true(any(tab$component == "magnitude" & tab$term == "gp_1" &
                    tab$parameter == "variance"))
  start <- magCompStarting(fit_gp, use = "last")
  expect_equal(length(start$magnitude$terms$gp_1$effects), n)
})

test_that("sampler supports composition iid factor terms", {
  # Small count table for fitting and output checks.
  n <- 12
  cls <- c("s", "m", "l", "xl")
  area <- rep(0.08, n)
  group <- rep(paste0("g", 1:6), length.out = n)
  rows <- expand.grid(plot = paste0("p", seq_len(n)), cls = cls)
  rows$cls <- factor(rows$cls, levels = sort(cls))
  i <- as.integer(rows$plot)
  rows$count <- 1L + (i + 2L * as.integer(rows$cls)) %% 7L
  rows$area <- area[i]
  rows$group <- group[i]
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     ordered = TRUE, design = fixedArea(area))

  fit <- magComp(
    dat,
    composition = ~ factor_iid(group, K = 2, anchors = c("m", "s")),
    control = magCompControl(
      n_samples = 35,
      n_burn = 15,
      n_thin = 1,
      n_chains = 2,
      chain_seeds = c(303, 404),
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )

  expect_equal(dim(mc_block_loading(fit)),
               c(length(cls) - 1L, 2L, 20L, 2L))
  expect_equal(dim(mc_block_factor(fit)), c(6L, 2L, 20L, 2L))
  expect_equal(dim(fit$starting$composition$blocks$cls$terms$factor_iid_1$loadings),
               c(length(cls) - 1L, 2L))
  expect_equal(dim(fit$starting$composition$blocks$cls$terms$factor_iid_1$scores),
               c(6L, 2L))
  expect_equal(fit$composition$process_terms[[1L]]$anchors, c("m", "s"))
  expect_equal(fit$composition$process_terms[[1L]]$factor_names, c("m", "s"))
  expect_equal(colnames(mc_block_loading(fit)), c("m", "s"))
  expect_equal(colnames(mc_block_factor(fit)), c("m", "s"))
  start_loading <- fit$starting$composition$blocks$cls$terms$factor_iid_1$loadings
  expect_loading_constraints(start_loading, c("m", "s"))
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)
  fixed_only <- mc_regression_proportion_mean(fit)
  expect_equal(dim(fixed_only), c(n, length(cls)))
  expect_equal(unname(rowSums(fixed_only)), rep(1, n), tolerance = 1e-10)
  expect_gt(mean(abs(fitted(fit, type = "proportion")$proportion - fixed_only)), 0)

  intercept_diag <- samples(fit, type = "combined", group = ".all")
  expect_s3_class(intercept_diag, "magCompSamples")
  intercept_summary <- summary(intercept_diag)
  expect_true(any(intercept_summary$component == "composition"))
  expect_true(all(is.finite(intercept_summary$mean)))
  expect_true(all(is.na(intercept_summary$rhat) | intercept_summary$rhat >= 1))
  trace_file <- tempfile(fileext = ".pdf")
  devices_before <- grDevices::dev.list()
  plot(intercept_diag, file = trace_file, n_col = 1, n_row = 2)
  expect_identical(grDevices::dev.list(), devices_before)
  expect_true(file.exists(trace_file))
  unlink(trace_file)

  tab <- summary(fit)
  expect_true("loadings" %in% tab$parameter)
  expect_loading_constraints(mc_block_loading(fit), c("m", "s"))

  start <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(length(start$magnitude$coefficients), 1L)
  expect_equal(dim(start$composition$blocks$cls$coefficients),
               c(length(cls) - 1L, 1L))
  expect_equal(dim(start$composition$blocks$cls$terms$factor_iid_1$loadings),
               c(length(cls) - 1L, 2L))
  expect_equal(dim(start$composition$blocks$cls$terms$factor_iid_1$scores),
               c(6L, 2L))

  refit <- magComp(
    dat,
    composition = ~ factor_iid(group, K = 2, anchors = c("m", "s")),
    starting = start,
    control = magCompControl(
      n_samples = 2,
      n_burn = 0,
      n_thin = 1,
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )
  expect_true(refit$fitted)
})

test_that("sampler supports composition NNGP factor terms", {
  # Small count table for fitting and output checks.
  n <- 12
  cls <- c("s", "m", "l", "xl")
  area <- rep(0.08, n)
  x <- rep(0:3, 3)
  ycoord <- rep(0:2, each = 4)
  decay <- 0.4
  rows <- expand.grid(plot = paste0("p", seq_len(n)), cls = cls)
  rows$cls <- factor(rows$cls, levels = sort(cls))
  i <- as.integer(rows$plot)
  rows$count <- 1L + (i + 2L * as.integer(rows$cls)) %% 7L
  rows$area <- area[i]
  rows$x <- x[i]
  rows$y <- ycoord[i]
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     design = fixedArea(area))

  set.seed(46)
  fit <- magComp(dat, composition = ~factor_nngp(x, y, K = 1, m = 5, anchors = "m", factor_names = "mid"),
    control = magCompControl(n_samples = 20, n_burn = 8, n_thin = 1, verbose = FALSE), smooth = magCompSmooth(intercept = "none",
        slope = "none", loading = "none"), starting = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_nngp_1 = list(covariance = list(phi = fixed(decay)))))),
        names(dat$allocation_blocks)[1L]))))

  term <- fit$composition$process_terms[[1L]]
  expect_equal(term$term_type, "factor_nngp")
  expect_equal(term$factor_names, "mid")
  expect_equal(term$anchors, "m")
  expect_equal(term$m, 5L)
  expect_equal(as.numeric(term$covariance_values), decay)
  expect_equal(length(term$order), n)
  expect_equal(dim(term$neighbor_index), c(n, 5L))
  expect_equal(dim(mc_block_loading(fit)),
               c(length(cls) - 1L, 1L, 12L))
  expect_equal(dim(mc_block_factor(fit)), c(n, 1L, 12L))
  expect_equal(colnames(mc_block_loading(fit)), "mid")
  expect_equal(rownames(mc_block_factor(fit)), dat$unit_levels)
  expect_loading_constraints(mc_block_loading(fit), "m", "mid")
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)

  start <- magCompStarting(fit, use = "last")
  expect_equal(dim(start$composition$blocks$cls$terms$factor_nngp_1$scores),
               c(n, 1L))
})

test_that("spatial composition factors support factor-specific decay", {
  set.seed(124)

  n <- 16
  cls <- c("a", "b", "c", "d")
  area <- rep(1, n)
  x <- runif(n)
  ycoord <- runif(n)
  coords <- cbind(x, ycoord)
  phi_true <- c(1.2, 3.6)
  D <- as.matrix(stats::dist(coords))
  w <- matrix(NA_real_, nrow = n, ncol = 2L)
  for(k in seq_along(phi_true))
    w[,k] <- t(chol(exp(-phi_true[k] * D))) %*% rnorm(n)
  loading <- matrix(c(1.0, 0.25, -0.15,
                      0.0, 1.0, 0.35),
                    nrow = length(cls) - 1L)
  alpha <- c(-0.3, -0.6, -0.9)

  y <- matrix(0L, nrow = n, ncol = length(cls),
              dimnames = list(paste0("p", seq_len(n)), cls))
  for(i in seq_len(n)){
    eta <- alpha + as.vector(loading %*% w[i,])
    pi_i <- plogis(eta)
    p_i <- numeric(length(cls))
    remaining <- 1
    for(j in seq_along(eta)){
      p_i[j] <- remaining * pi_i[j]
      remaining <- remaining * (1 - pi_i[j])
    }
    p_i[length(cls)] <- remaining
    y[i,] <- as.vector(rmultinom(1, size = 80, prob = p_i))
  }

  trees <- do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(
      plot = rownames(y)[i],
      cls = cls,
      tree_count = as.integer(y[i,]),
      area = area[i],
      x = x[i],
      y = ycoord[i]
    )
  }))
  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    tree_count = tree_count,
    coords = c(x, y),
    design = fixedArea(area)
  )

  fit_gp <- magComp(dat, composition = ~factor_gp(x, y, K = 2, anchors = c("a", "b"), factor_names = c("f1",
    "f2")), priors = magCompPriors(composition = list(blocks = list(cls = list(terms = list(factor_gp_1 = list(covariance = list(phi = uniform(0.5,
    6)))))))), starting = list(composition = list(blocks = list(cls = list(terms = list(factor_gp_1 = list(covariance = list(phi = c(f1 = 1.5,
    f2 = 2.5)))))))), control = magCompControl(n_samples = 10, n_burn = 5, n_thin = 1, n_chains = 2,
    chain_seeds = c(501, 502), verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none",
    loading = "none"), tuning = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_gp_1 = list(covariance = list(phi = c(0.08,
    0.12)))))), names(dat$allocation_blocks)[1L]))))
  expect_equal(as.numeric(fit_gp$starting$composition$blocks$cls$terms$factor_gp_1$covariance$phi),
               c(1.5, 2.5))
  expect_equal(dim(mc_block_factor_phi(fit_gp)), c(2L, 5L, 2L))
  expect_equal(rownames(mc_block_factor_phi(fit_gp)), c("f1:phi", "f2:phi"))
  expect_equal(rownames(fit_gp$backend$composition_block_covariance_acceptance),
               c("cls:f1", "cls:f2"))
  expect_true(all(is.finite(mc_block_factor_phi(fit_gp))))
  expect_true(all(mc_block_factor_phi(fit_gp) >= 0.5))
  expect_true(all(mc_block_factor_phi(fit_gp) <= 6))

  summary <- summary(fit_gp)
  expect_true(any(summary$component == "composition" & summary$block == "cls" & summary$term == "factor_gp_1" &
                    summary$parameter == "phi" & summary$index == "f1"))
  expect_true(any(summary$component == "composition" & summary$block == "cls" & summary$term == "factor_gp_1" &
                    summary$parameter == "phi" & summary$index == "f2"))

  start <- magCompStarting(fit_gp, use = "last", chain = 2)
  expect_equal(length(start$composition$blocks$cls$terms$factor_gp_1$covariance$phi),
               2L)
  expect_equal(names(start$composition$blocks$cls$terms$factor_gp_1$covariance$phi),
               c("f1", "f2"))

  newdata <- data.frame(plot = c("new1", "new2"),
                        x = c(0.2, 0.8),
                        y = c(0.8, 0.2))
  set.seed(601)
  pred <- predict(fit_gp, newdata, type = "proportion")$proportion
  expect_equal(dim(pred), c(2L, length(cls)))
  expect_equal(unname(rowSums(pred)), rep(1, 2L), tolerance = 1e-10)

  set.seed(602)
  pred_samples <- predict(fit_gp, newdata, type = "proportion",
                          summary = "samples")$proportion
  fit_gp_phi2 <- fit_gp
  fit_gp_phi2$backend$composition_block_covariance_samples$cls[2L,,] <- 6
  set.seed(602)
  pred_samples_phi2 <- predict(fit_gp_phi2, newdata, type = "proportion",
                               summary = "samples")$proportion
  expect_gt(max(abs(pred_samples - pred_samples_phi2)), 1e-8)

  fit_nngp <- magComp(dat, composition = ~factor_nngp(x, y, K = 2, m = 4, anchors = c("a", "b"), factor_names = c("f1",
    "f2")), priors = magCompPriors(composition = list(blocks = list(cls = list(terms = list(factor_nngp_1 = list(covariance = list(phi = uniform(0.5,
    6)))))))), starting = list(composition = list(blocks = list(cls = list(terms = list(factor_nngp_1 = list(covariance = list(phi = c(f1 = 1.5,
    f2 = 2.5)))))))), control = magCompControl(n_samples = 10, n_burn = 5, n_thin = 1, n_chains = 2,
    chain_seeds = c(503, 504), verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none",
    loading = "none"), tuning = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_nngp_1 = list(covariance = list(phi = c(0.08,
    0.12)))))), names(dat$allocation_blocks)[1L]))))
  expect_equal(dim(mc_block_factor_phi(fit_nngp)), c(2L, 5L, 2L))
  expect_equal(rownames(mc_block_factor_phi(fit_nngp)), c("f1:phi", "f2:phi"))
  expect_equal(rownames(fit_nngp$backend$composition_block_covariance_acceptance),
               c("cls:f1", "cls:f2"))
  set.seed(603)
  pred <- predict(fit_nngp, newdata, type = "proportion")$proportion
  expect_equal(dim(pred), c(2L, length(cls)))
  expect_equal(unname(rowSums(pred)), rep(1, 2L), tolerance = 1e-10)

  set.seed(604)
  pred_samples <- predict(fit_nngp, newdata, type = "proportion",
                          summary = "samples")$proportion
  fit_nngp_phi2 <- fit_nngp
  fit_nngp_phi2$backend$composition_block_covariance_samples$cls[2L,,] <- 6
  set.seed(604)
  pred_samples_phi2 <- predict(fit_nngp_phi2, newdata, type = "proportion",
                               summary = "samples")$proportion
  expect_gt(max(abs(pred_samples - pred_samples_phi2)), 1e-8)
})

test_that("sampler keeps spatial factor coordinates valid when NNGP graphs are reused", {
  set.seed(631)

  n <- 121
  species <- paste0("sp", seq_len(12))
  dbh <- paste0("d", seq_len(8))
  plot <- paste0("p", seq_len(n))
  area <- rep(0.08, n)
  z <- rnorm(n)
  x <- rep(seq(0, 2, length.out = 11), times = 11) + rnorm(n, 0, 0.02)
  ycoord <- rep(seq(0, 2, length.out = 11), each = 11) + rnorm(n, 0, 0.02)

  cell_prob <- as.vector(outer(seq_along(species), seq_along(dbh),
                               function(a, b) exp(-0.18 * a - 0.28 * b)))
  cell_prob <- cell_prob / sum(cell_prob)
  counts <- matrix(0L, nrow = n, ncol = length(cell_prob))
  for(i in seq_len(n)){
    total_i <- max(1L, rpois(1, area[i] * exp(5.8 - 0.3 * z[i])))
    counts[i,] <- as.vector(rmultinom(1, total_i, cell_prob))
  }

  rows <- expand.grid(
    plot = plot,
    species = species,
    dbh = dbh,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  rows$area <- rep(area, times = length(species) * length(dbh))
  rows$z <- rep(z, times = length(species) * length(dbh))
  rows$x <- rep(x, times = length(species) * length(dbh))
  rows$ycoord <- rep(ycoord, times = length(species) * length(dbh))
  rows$tree_count <- as.integer(counts)
  rows$plot <- factor(rows$plot, levels = plot)
  rows$species <- factor(rows$species, levels = species)
  rows$dbh <- factor(rows$dbh, levels = dbh, ordered = TRUE)

  dat <- magCompData(
    rows,
    unit = plot,
    categories = c(species, dbh),
    coords = c(x, ycoord),
    tree_count = tree_count,
    design = fixedArea(area),
    ordered = c(species = FALSE, dbh = TRUE)
  )

  fit <- magComp(dat, magnitude = ~z + nngp(x, ycoord, m = 15), composition = magCompComposition(species ~
    z + factor_nngp(x, ycoord, K = 1, m = 15, anchors = "sp2"), dbh | species ~ z), control = magCompControl(n_samples = 2,
    n_burn = 1, n_thin = 1, n_chains = 1, n_cores = 1, verbose = FALSE), smooth = magCompSmooth(intercept = "none",
    slope = "none", loading = "none"), starting = list(magnitude = list(terms = list(nngp_1 = list(covariance = list(phi = fixed(3))))),
    composition = list(blocks = list(species = list(terms = list(factor_nngp_1 = list(covariance = list(phi = fixed(3)))))))))

  expect_true(fit$fitted)
  expect_equal(fit$backend$nngp_graph_count, 1L)
  expect_equal(dim(mc_mag_process(fit, "nngp_1")), c(n, 1L))
  expect_equal(dim(mc_block_factor(fit, "species")), c(n, 1L, 1L))
  expect_equal(dim(fit$starting$composition$blocks$species$terms$factor_nngp_1$loadings),
               c(length(species) - 1L, 1L))
  expect_equal(dim(fit$starting$composition$blocks$species$terms$factor_nngp_1$scores),
               c(n, 1L))
})

test_that("sampler updates factor_nngp decay", {
  set.seed(64)

  n <- 18
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  x <- runif(n)
  ycoord <- runif(n)
  alpha <- c(-0.2, -0.7)
  loading <- c(0.6, 1.0)
  w <- rnorm(n)

  counts <- matrix(0L, n, length(cls))
  for(i in seq_len(n)){
    total_i <- max(1L, rpois(1, area[i] * 120))
    eta <- alpha + loading * w[i]
    pi_i <- plogis(eta)
    p_i <- c(pi_i[1L], (1 - pi_i[1L]) * pi_i[2L],
             prod(1 - pi_i))
    counts[i,] <- as.vector(rmultinom(1, size = total_i, prob = p_i))
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, counts[i,]),
      area = area[i],
      x = x[i],
      ycoord = ycoord[i]
    ))
  }

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit <- magComp(dat, composition = ~factor_nngp(x, ycoord, K = 1, m = 4, anchors = "m", factor_names = "mid"),
    priors = magCompPriors(composition = list(blocks = list(cls = list(terms = list(factor_nngp_1 = list(covariance = list(phi = uniform(0.2,
        8)))))))), control = magCompControl(n_samples = 16, n_burn = 6, n_thin = 1, n_chains = 2,
        chain_seeds = c(641, 642), metropolis = list(batch_length = 3), verbose = FALSE), smooth = magCompSmooth(intercept = "none",
        slope = "none", loading = "none"), tuning = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_nngp_1 = list(covariance = list(phi = 0.15))))),
        names(dat$allocation_blocks)[1L]))))

  term <- fit$composition$process_terms[[1L]]
  expect_false(any(term$covariance_fixed))
  expect_equal(unname(c(term$covariance_lower, term$covariance_upper)), c(0.2, 8))
  expect_equal(dim(mc_block_factor_phi(fit)), c(1L, 10L, 2L))
  expect_equal(fit$backend$nngp_graph_count, 1L)
  expect_true(all(is.finite(mc_block_factor_phi(fit))))
  expect_true(all(mc_block_factor_phi(fit) > term$covariance_lower[1,1]))
  expect_true(all(mc_block_factor_phi(fit) < term$covariance_upper[1,1]))
  phi_accept <- fit$backend$composition_block_covariance_acceptance["cls",]
  expect_equal(length(phi_accept), 2L)
  expect_true(all(is.finite(phi_accept)))
  expect_true(all(phi_accept >= 0))
  expect_true(all(phi_accept <= 1))
  decay_tuning <- fit$backend$composition_block_covariance_tuning["cls",]
  expect_equal(length(decay_tuning), 2L)
  expect_true(all(is.finite(decay_tuning)))
  expect_equal(length(fit$backend$composition_block_covariance_batch_acceptance),
               1L)
  expect_equal(length(
    fit$backend$composition_block_covariance_batch_acceptance$cls[[1L]]
  ), 5L)
  expect_equal(length(
    fit$backend$composition_block_covariance_tuning_history$cls[[1L]]
  ), 5L)

  tab <- summary(fit)
  expect_true(any(tab$component == "composition" & tab$block == "cls" & tab$term == "factor_nngp_1" &
                    tab$parameter == "phi"))

  spatial_diag <- spatialDiagnostics(fit)
  expect_s3_class(spatial_diag, "magCompSpatialDiagnostics")
  expect_true(any(spatial_diag$summary$component == "composition" &
                    spatial_diag$summary$quantity == "phi"))
  expect_true(any(spatial_diag$summary$component == "composition" &
                    spatial_diag$summary$quantity == "loading_norm" & spatial_diag$summary$factor == "mid"))
  expect_equal(nrow(spatial_diag$acceptance), 2L)
  expect_true(any(spatial_diag$correlations$component == "composition" &
                    spatial_diag$correlations$x == "phi" &
                    spatial_diag$correlations$y == "loading_norm" & spatial_diag$correlations$factor == "mid"))

  newdata <- data.frame(plot = paste0("new", seq_len(4L)),
                        x = runif(4L),
                        ycoord = runif(4L),
                        area = area[1L],
                        stand = c("a", "a", "b", "b"))
  set.seed(456)
  pred_full <- predict(fit, newdata, type = "proportion",
                       summary = "samples")$proportion
  expect_equal(dim(pred_full), c(4L, length(cls), 10L, 2L))
  expect_equal(as.numeric(apply(pred_full, c(1L, 3L, 4L), sum)),
               rep(1, 4L * 10L * 2L),
               tolerance = 1e-10)

  set.seed(456)
  stand_density <- standSummary(predict(fit, newdata, summary = "samples"),
    data = newdata, unit = plot, group = stand, represented_area = area,
    type = "density", summary = "samples")$density
  expect_equal(dim(stand_density), c(2L, length(cls), 10L, 2L))

  start <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(length(start$composition$blocks$cls$terms$factor_nngp_1$covariance$phi),
               1L)
  expect_true(is.finite(
    start$composition$blocks$cls$terms$factor_nngp_1$covariance$phi
  ))
})

test_that("sampler supports composition dense GP factors", {
  # Small count table for fitting and output checks.
  n <- 12
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  x <- rep(0:3, 3)
  ycoord <- rep(0:2, each = 4)
  rows <- expand.grid(plot = paste0("p", seq_len(n)), cls = cls)
  rows$cls <- factor(rows$cls, levels = sort(cls))
  i <- as.integer(rows$plot)
  rows$count <- 1L + (i + 2L * as.integer(rows$cls)) %% 7L
  rows$area <- area[i]
  rows$x <- x[i]
  rows$ycoord <- ycoord[i]
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     design = fixedArea(area))

  fit <- magComp(dat, composition = ~factor_gp(x, ycoord, K = 1, anchors = "m", factor_names = "mid"), control = magCompControl(n_samples = 13,
    n_burn = 5, n_thin = 1, n_chains = 2, chain_seeds = c(651, 652), verbose = FALSE), smooth = magCompSmooth(intercept = "none",
    slope = "none", loading = "none"), tuning = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_gp_1 = list(covariance = list(phi = 0.15))))),
    names(dat$allocation_blocks)[1L]))))

  term <- fit$composition$process_terms[[1L]]
  expect_equal(term$term_type, "factor_gp")
  expect_equal(term$factor_names, "mid")
  expect_equal(term$anchors, "m")
  expect_null(term$neighbor_index)
  expect_equal(dim(mc_block_loading(fit)),
               c(length(cls) - 1L, 1L, 8L, 2L))
  expect_equal(dim(mc_block_factor(fit)), c(n, 1L, 8L, 2L))
  expect_equal(dim(mc_block_factor_phi(fit)), c(1L, 8L, 2L))
  expect_loading_constraints(mc_block_loading(fit), "m", "mid")
  expect_equal(unname(rowSums(fitted(fit, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)

  tab <- summary(fit)
  expect_true(any(tab$component == "composition" & tab$block == "cls" & tab$term == "factor_gp_1" & tab$parameter == "phi"))

  spatial_diag <- spatialDiagnostics(fit)
  expect_true(any(spatial_diag$summary$component == "composition" &
                    spatial_diag$summary$block == "cls" & spatial_diag$summary$term == "factor_gp_1" &
                    spatial_diag$summary$quantity == "phi"))

  start <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(dim(start$composition$blocks$cls$terms$factor_gp_1$scores),
               c(n, 1L))
  expect_equal(length(start$composition$blocks$cls$terms$factor_gp_1$covariance$phi),
               1L)
})

test_that("sampler supports composition factor varying coefficients", {
  set.seed(66)

  n <- 18
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  group <- rep(paste0("g", 1:3), length.out = n)
  x <- runif(n)
  ycoord <- runif(n)
  z <- runif(n, 0.5, 1.5)
  w <- c(g1 = 0.7, g2 = -0.4, g3 = 0.2)
  alpha <- c(-0.2, -0.7)
  loading <- c(0.6, 1.0)

  counts <- matrix(0L, n, length(cls))
  for(i in seq_len(n)){
    eta <- alpha + z[i] * loading * w[group[i]]
    pi_i <- plogis(eta)
    p_i <- c(pi_i[1L], (1 - pi_i[1L]) * pi_i[2L],
             prod(1 - pi_i))
    counts[i,] <- 1L + as.vector(rmultinom(1, size = 17, prob = p_i))
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, counts[i,]),
      area = area[i],
      group = group[i],
      x = x[i],
      ycoord = ycoord[i],
      z = z[i]
    ))
  }

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit_iid <- magComp(
    dat,
    composition = ~ z + z:factor_iid(group, K = 1, anchors = "m"),
    control = magCompControl(
      n_samples = 9,
      n_burn = 4,
      n_thin = 1,
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )
  term <- fit_iid$composition$process_terms[[1L]]
  expect_equal(term$term_type, "factor_iid")
  expect_equal(term$covariate_name, "z")
  expect_equal(length(term$covariate), n)
  expect_equal(dim(mc_block_factor(fit_iid)), c(3L, 1L, 5L))
  expect_loading_constraints(mc_block_loading(fit_iid), "m")
  expect_gt(mean(abs(fitted(fit_iid, type = "proportion")$proportion -
                       mc_regression_proportion_mean(fit_iid))), 0)

  fit_nngp <- magComp(dat, composition = ~z + z:factor_nngp(x, ycoord, K = 1, m = 4, anchors = "m"), control = magCompControl(n_samples = 9,
    n_burn = 4, n_thin = 1, verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none",
    loading = "none"), starting = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_nngp_1 = list(covariance = list(phi = fixed(2)))))),
    names(dat$allocation_blocks)[1L]))))
  term <- fit_nngp$composition$process_terms[[1L]]
  expect_equal(term$term_type, "factor_nngp")
  expect_equal(term$covariate_name, "z")
  expect_equal(dim(mc_block_factor(fit_nngp)), c(n, 1L, 5L))
  expect_loading_constraints(mc_block_loading(fit_nngp), "m")
  expect_equal(unname(rowSums(fitted(fit_nngp, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)

  fit_gp <- magComp(dat, composition = ~z + z:factor_gp(x, ycoord, K = 1, anchors = "m"), control = magCompControl(n_samples = 7,
    n_burn = 3, n_thin = 1, verbose = FALSE), smooth = magCompSmooth(intercept = "none", slope = "none",
    loading = "none"), starting = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_gp_1 = list(covariance = list(phi = fixed(2)))))),
    names(dat$allocation_blocks)[1L]))))
  term <- fit_gp$composition$process_terms[[1L]]
  expect_equal(term$term_type, "factor_gp")
  expect_equal(term$covariate_name, "z")
  expect_equal(dim(mc_block_factor(fit_gp)), c(n, 1L, 4L))
  expect_loading_constraints(mc_block_loading(fit_gp), "m")
  expect_equal(unname(rowSums(fitted(fit_gp, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)
})

test_that("sampler supports multiple composition factor terms per block", {
  set.seed(67)

  n <- 16
  cls <- c("s", "m", "l")
  area <- rep(0.08, n)
  group <- rep(paste0("g", 1:4), length.out = n)
  x <- runif(n)
  ycoord <- runif(n)
  z <- runif(n, 0.5, 1.5)
  alpha <- c(-0.1, -0.8)
  group_effect <- c(g1 = 0.6, g2 = -0.4, g3 = 0.2, g4 = -0.1)

  counts <- matrix(0L, n, length(cls))
  for(i in seq_len(n)){
    eta <- alpha + (0.4 + 0.3 * z[i]) * group_effect[group[i]]
    pi_i <- plogis(eta)
    p_i <- c(pi_i[1L], (1 - pi_i[1L]) * pi_i[2L],
             prod(1 - pi_i))
    counts[i,] <- 1L + as.vector(rmultinom(1, size = 16, prob = p_i))
  }

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, counts[i,]),
      area = area[i],
      group = group[i],
      x = x[i],
      ycoord = ycoord[i],
      z = z[i]
    ))
  }

  trees$cls <- factor(trees$cls, levels = cls)
  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )
  ctrl <- magCompControl(
    n_samples = 7,
    n_burn = 3,
    n_thin = 1,
    composition_link = "probit",
    verbose = FALSE
  )
  sm <- magCompSmooth(intercept = "none", slope = "none", loading = "none")

  mixed_formula <- ~ z +
    factor_nngp(x, ycoord, K = 1, m = 4, anchors = "m") +
    z:factor_iid(group, K = 1, anchors = "m") +
    factor_iid(group, K = 1, anchors = "m") +
    z:factor_nngp(x, ycoord, K = 1, m = 4, anchors = "m")

  fit_mixed_start <- magComp(
    dat,
    composition = mixed_formula,
    control = magCompControl(run_sampler = FALSE,
                             composition_link = "probit",
                             verbose = FALSE),
    smooth = sm
  )
  start_loading <- lapply(fit_mixed_start$starting$composition$blocks$cls$terms, "[[", "loadings")
  start_factor <- lapply(fit_mixed_start$starting$composition$blocks$cls$terms, "[[", "scores")
  expect_setequal(names(start_loading),
               c("factor_iid_1", "factor_iid_2",
                 "factor_nngp_1", "factor_nngp_2"))
  expect_equal(names(start_factor), names(start_loading))
  expect_equal(dim(fit_mixed_start$starting$composition$blocks$cls$coefficients),
               c(length(cls) - 1L, 2L))
  for(term_name in names(start_loading))
    expect_loading_constraints(start_loading[[term_name]], "m")
  expect_equal(dim(start_factor$factor_iid_1), c(4L, 1L))
  expect_equal(dim(start_factor$factor_iid_2), c(4L, 1L))
  expect_equal(dim(start_factor$factor_nngp_1), c(n, 1L))
  expect_equal(dim(start_factor$factor_nngp_2), c(n, 1L))
  expect_true(all(vapply(start_loading, function(x) all(is.finite(x)),
                         logical(1L))))
  expect_true(all(vapply(start_factor, function(x) all(is.finite(x)),
                         logical(1L))))

  true_start <- list(composition = list(blocks = list(cls = list(coefficients = matrix(c(-0.1, -0.8, 0.2, -0.1),
                   nrow = length(cls) - 1L,
                   dimnames = list(cls[-length(cls)],
                                   c("(Intercept)", "z"))), terms = list(factor_nngp_1 = list(loadings = matrix(c(0.2, 1), nrow = length(cls) - 1L,
                               dimnames = list(cls[-length(cls)], "m")), scores = matrix(seq(-0.3, 0.3, length.out = n),
                               nrow = n,
                               dimnames = list(paste0("p", seq_len(n)),
                                               "m"))), factor_iid_1 = list(loadings = matrix(c(-0.4, 1), nrow = length(cls) - 1L,
                              dimnames = list(cls[-length(cls)], "m")), scores = matrix(c(g1 = -0.2, g2 = 0.1, g3 = 0.2, g4 = -0.1),
                              ncol = 1L,
                              dimnames = list(paste0("g", 1:4), "m"))), factor_iid_2 = list(loadings = matrix(c(0.5, 1), nrow = length(cls) - 1L,
                              dimnames = list(cls[-length(cls)], "m")), scores = matrix(c(g1 = 0.3, g2 = -0.1, g3 = 0.0, g4 = -0.2),
                              ncol = 1L,
                              dimnames = list(paste0("g", 1:4), "m"))), factor_nngp_2 = list(loadings = matrix(c(0.3, 1), nrow = length(cls) - 1L,
                               dimnames = list(cls[-length(cls)], "m")), scores = matrix(seq(0.2, -0.2, length.out = n),
                               nrow = n,
                               dimnames = list(paste0("p", seq_len(n)),
                                               "m"))))))))
  fit_mixed_true_start <- magComp(
    dat,
    composition = mixed_formula,
    starting = true_start,
    control = magCompControl(run_sampler = FALSE,
                             composition_link = "probit",
                             verbose = FALSE),
    smooth = sm
  )
  expect_equal(sort(names(lapply(fit_mixed_true_start$starting$composition$blocks$cls$terms, "[[", "scores"))),
               c("factor_iid_1", "factor_iid_2",
                 "factor_nngp_1", "factor_nngp_2"))
  for(term_name in names(lapply(fit_mixed_true_start$starting$composition$blocks$cls$terms, "[[", "loadings"))){
    expect_loading_constraints(
      fit_mixed_true_start$starting$composition$blocks$cls$terms[[term_name]]$loadings,
      "m"
    )
  }

  fit_mixed <- magComp(
    dat,
    composition = mixed_formula,
    control = ctrl,
    smooth = sm
  )
  expect_equal(names(fit_mixed$backend$composition_block_score_samples),
               c("cls:factor_nngp_1", "cls:factor_iid_1",
                 "cls:factor_iid_2", "cls:factor_nngp_2"))
  expect_equal(names(fit_mixed$backend$composition_block_covariance_samples),
               c("cls:factor_nngp_1", "cls:factor_nngp_2"))
  for(sample_name in names(fit_mixed$backend$composition_block_loading_samples)){
    expect_loading_constraints(
      fit_mixed$backend$composition_block_loading_samples[[sample_name]], "m"
    )
  }
  newdata_mixed <- data.frame(plot = paste0("mixed", seq_len(3L)),
                              x = runif(3L),
                              ycoord = runif(3L),
                              z = runif(3L, 0.5, 1.5),
                              group = c("g1", "g2", "g3"),
                              area = area[1L])
  set.seed(68)
  pred_mixed <- predict(fit_mixed, newdata_mixed, type = "proportion",
                        summary = "samples")$proportion
  expect_equal(dim(pred_mixed), c(3L, length(cls), 4L))
  expect_equal(as.numeric(apply(pred_mixed, c(1L, 3L), sum)),
               rep(1, 3L * 4L),
               tolerance = 1e-10)

  fit_iid <- magComp(
    dat,
    composition = ~ z + factor_iid(group, K = 1, anchors = "m") +
      z:factor_iid(group, K = 1, anchors = "m"),
    control = ctrl,
    smooth = sm
  )
  expect_equal(names(fit_iid$backend$composition_block_score_samples),
               c("cls:factor_iid_1", "cls:factor_iid_2"))
  expect_equal(dim(fit_iid$backend$composition_block_score_samples[[1L]]),
               c(4L, 1L, 4L))
  expect_loading_constraints(
    fit_iid$backend$composition_block_loading_samples[["cls:factor_iid_1"]],
    "m"
  )
  expect_loading_constraints(
    fit_iid$backend$composition_block_loading_samples[["cls:factor_iid_2"]],
    "m"
  )
  expect_gt(mean(abs(fitted(fit_iid, type = "proportion")$proportion -
                       mc_regression_proportion_mean(fit_iid))), 0)

  start_iid <- magCompStarting(fit_iid, use = "last")
  expect_equal(names(lapply(start_iid$composition$blocks$cls$terms, "[[", "scores")),
               c("factor_iid_1", "factor_iid_2"))

  fit_nngp <- magComp(dat, composition = ~z + factor_nngp(x, ycoord, K = 1, m = 4, anchors = "m") + z:factor_nngp(x,
    ycoord, K = 1, m = 4, anchors = "m"), control = ctrl, smooth = sm, starting = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_nngp_1 = list(covariance = list(phi = fixed(2))),
    factor_nngp_2 = list(covariance = list(phi = fixed(2)))))), names(dat$allocation_blocks)[1L]))))
  expect_equal(names(fit_nngp$backend$composition_block_score_samples),
               c("cls:factor_nngp_1", "cls:factor_nngp_2"))
  expect_equal(names(fit_nngp$backend$composition_block_covariance_samples),
               c("cls:factor_nngp_1", "cls:factor_nngp_2"))
  expect_equal(fit_nngp$backend$nngp_graph_count, 1L)
  expect_loading_constraints(
    fit_nngp$backend$composition_block_loading_samples[["cls:factor_nngp_1"]],
    "m"
  )
  expect_loading_constraints(
    fit_nngp$backend$composition_block_loading_samples[["cls:factor_nngp_2"]],
    "m"
  )
  expect_equal(unname(rowSums(fitted(fit_nngp, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)

  newdata <- data.frame(plot = paste0("new", seq_len(3L)),
                        x = runif(3L),
                        ycoord = runif(3L),
                        z = runif(3L, 0.5, 1.5),
                        group = c("g1", "g2", "g3"),
                        area = area[1L])
  set.seed(67)
  pred <- predict(fit_nngp, newdata, type = "proportion",
                  summary = "samples")$proportion
  expect_equal(dim(pred), c(3L, length(cls), 4L))
  expect_equal(as.numeric(apply(pred, c(1L, 3L), sum)),
               rep(1, 3L * 4L),
               tolerance = 1e-10)

  fit_gp <- magComp(dat, composition = ~z + factor_gp(x, ycoord, K = 1, anchors = "m") + z:factor_gp(x, ycoord,
    K = 1, anchors = "m"), control = ctrl, smooth = sm, starting = list(composition = list(blocks = stats::setNames(list(list(terms = list(factor_gp_1 = list(covariance = list(phi = fixed(2))),
    factor_gp_2 = list(covariance = list(phi = fixed(2)))))), names(dat$allocation_blocks)[1L]))))
  expect_equal(names(fit_gp$backend$composition_block_score_samples),
               c("cls:factor_gp_1", "cls:factor_gp_2"))
  expect_equal(names(fit_gp$backend$composition_block_covariance_samples),
               c("cls:factor_gp_1", "cls:factor_gp_2"))
  expect_loading_constraints(
    fit_gp$backend$composition_block_loading_samples[["cls:factor_gp_1"]],
    "m"
  )
  expect_loading_constraints(
    fit_gp$backend$composition_block_loading_samples[["cls:factor_gp_2"]],
    "m"
  )
  expect_equal(unname(rowSums(fitted(fit_gp, type = "proportion")$proportion)),
               rep(1, n),
               tolerance = 1e-10)
})

test_that("sampler rejects unsupported fitted terms", {
  trees <- data.frame(
    plot = rep(letters[1:3], each = 3),
    cls = rep(c("s", "m", "l"), 3),
    area = rep(0.1, 9),
    p95 = rep(c(12, 15, 18), each = 3)
  )
  trees$cls <- factor(trees$cls, levels = c("l", "m", "s"))
  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    ordered = c(cls = TRUE),
    design = fixedArea(area)
  )

  expect_error(
    magComp(
      dat,
      composition = ~ iid(p95),
      control = magCompControl(verbose = FALSE),
      smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
    ),
    "factor_iid"
  )

  expect_error(
    magComp(
      dat,
      control = magCompControl(verbose = FALSE),
      smooth = magCompSmooth(loading = "rw1")
    ),
    "loading"
  )
})

test_that("native spatial prediction matches R reference helpers", {
  coords_obs <- cbind(x = c(0, 1, 0.5, 1.2),
                      ycoord = c(0, 0.2, 1, 1.1))
  coords_new <- cbind(x = c(0.25, 0.9),
                      ycoord = c(0.4, 0.8))
  w <- c(0.2, -0.4, 1.1, 0.5)
  decay <- 0.7
  variance <- 1.4

  r_gp <- getFromNamespace("magcomp_gp_predict_one", "magComp")
  r_nngp <- getFromNamespace("magcomp_nngp_predict_one", "magComp")

  set.seed(101)
  gp_r <- r_gp(w, coords_obs, coords_new, decay, variance)
  set.seed(101)
  gp_cpp <- .Call("magcomp_spatial_predict_cpp",
                  matrix(w, ncol = 1L),
                  coords_obs,
                  coords_new,
                  matrix(decay, 1L, 1L),
                  matrix(variance, 1L, 1L),
                  as.integer(0L),
                  "gp",
                  FALSE,
                  "exp", 1.0, PACKAGE = "magComp")[,1L]
  expect_equal(gp_cpp, gp_r, tolerance = 1e-10)

  set.seed(202)
  nngp_r <- r_nngp(w, coords_obs, coords_new, decay, variance, 3L)
  set.seed(202)
  nngp_cpp <- .Call("magcomp_spatial_predict_cpp",
                    matrix(w, ncol = 1L),
                    coords_obs,
                    coords_new,
                    matrix(decay, 1L, 1L),
                    matrix(variance, 1L, 1L),
                    as.integer(3L),
                    "nngp",
                    FALSE,
                    "exp", 1.0, PACKAGE = "magComp")[,1L]
  expect_equal(nngp_cpp, nngp_r, tolerance = 1e-10)

  # Reusing the neighbor geometry must still use each draw's covariance
  # parameters and the earlier new-location effects from that same draw.
  coords_new <- rbind(c(0.9, 0.8), c(0.25, 0.4), c(0.28, 0.44))
  for(is_factor in c(FALSE, TRUE)){
    K <- if(is_factor) 2L else 1L
    effect <- array(seq_len(4L * K * 3L * 2L) / 17,
                     c(4L, K, 3L, 2L))
    decay <- array(seq_len(K * 3L * 2L) / 7, c(K, 3L, 2L))
    variance <- array(if(is_factor) 1 else seq_len(6L) / 3,
                       c(1L, 3L, 2L))
    reference <- array(NA_real_, c(3L, K, 3L, 2L))
    set.seed(203)
    for(ch in 1:2) for(s in 1:3) for(k in seq_len(K))
      reference[,k,s,ch] <- r_nngp(effect[,k,s,ch], coords_obs, coords_new,
                                   decay[k,s,ch], variance[1L,s,ch], 2L)
    if(!is_factor){
      dim(effect) <- c(4L, 3L, 2L)
      dim(reference) <- c(3L, 3L, 2L)
    }
    set.seed(203)
    predicted <- .Call("magcomp_spatial_predict_cpp", effect, coords_obs,
                        coords_new, decay, variance, 2L, "nngp", is_factor,
                        "exp", 1.0, PACKAGE = "magComp")
    expect_equal(predicted, reference, tolerance = 1e-10)
  }
})

test_that("native prediction kernels match direct R calculations", {
  X <- cbind("(Intercept)" = c(1, 1, 1),
             z = c(-0.5, 0.2, 1.1))
  gamma <- array(c(0.3, -0.2,
                   0.5,  0.1,
                   0.4,  0.3,
                   0.1, -0.4,
                   0.2,  0.0,
                   0.6,  0.2),
                 dim = c(2, 3, 2))
  effect <- array(seq(-0.3, 0.4, length.out = 18),
                  dim = c(3, 3, 2))
  covariate <- c(1, -0.5, 0.25)

  mag_cpp <- .Call("magcomp_predict_magnitude_cpp",
                   X, gamma, effect, covariate,
                   PACKAGE = "magComp")
  mag_r <- array(NA_real_, dim = c(3, 3, 2))
  for(ch in 1:2){
    for(s in 1:3){
      mag_r[,s,ch] <- exp(drop(X %*% gamma[,s,ch]) +
                            covariate * effect[,s,ch])
    }
  }
  expect_equal(mag_cpp, mag_r, tolerance = 1e-12)

  coef <- array(seq(-0.8, 0.7, length.out = 24),
                dim = c(4, 2, 3))
  loading <- array(seq(0.1, 0.5, length.out = 12),
                   dim = c(4, 1, 3))
  factor <- array(c(0.2, -0.1,
                    0.0,  0.3,
                    0.4, -0.2),
                  dim = c(2, 1, 3))
  X_block <- cbind("(Intercept)" = c(1, 1),
                   z = c(-1, 0.5))
  parent_index <- as.integer(c(1, 2))

  prop_cpp <- .Call("magcomp_predict_block_proportion_cpp",
                    X_block, coef, loading, factor, c(1, 1),
                    parent_index, as.integer(3), "logit",
                    PACKAGE = "magComp")
  prop_r <- array(NA_real_, dim = c(2, 3, 3))
  for(s in 1:3){
    for(i in 1:2){
      remaining <- 1
      for(j in 1:2){
        row <- (parent_index[i] - 1L) * 2L + j
        eta <- sum(X_block[i,] * coef[row,,s]) +
          loading[row,1,s] * factor[i,1,s]
        pi_value <- plogis(eta)
        prop_r[i,j,s] <- remaining * pi_value
        remaining <- remaining * (1 - pi_value)
      }
      prop_r[i,3,s] <- remaining
    }
  }
  expect_equal(prop_cpp, prop_r, tolerance = 1e-12)

  parent <- array(c(0.6, 0.4,
                    0.4, 0.6,
                    0.7, 0.2,
                    0.3, 0.8),
                  dim = c(2, 2, 1, 2))
  child <- array(c(0.5, 0.2, 0.1, 0.3,
                   0.3, 0.5, 0.4, 0.2,
                   0.2, 0.3, 0.5, 0.5,
                   0.4, 0.6, 0.3, 0.2,
                   0.4, 0.2, 0.4, 0.5,
                   0.2, 0.2, 0.3, 0.3),
                 dim = c(4, 3, 1, 2))
  combined_cpp <- .Call("magcomp_combine_conditional_cpp",
                        parent, child, PACKAGE = "magComp")
  combined_r <- array(NA_real_, dim = c(2, 6, 1, 2))
  for(ch in 1:2){
    for(a in 1:2){
      for(b in 1:3){
        col <- (a - 1L) * 3L + b
        combined_r[,col,1,ch] <- parent[,a,1,ch] *
          child[(a - 1L) * 2L + (1:2),b,1,ch]
      }
    }
  }
  expect_equal(combined_cpp, combined_r, tolerance = 1e-12)

  magnitude <- array(c(10, 20,
                       30, 40),
                     dim = c(2, 1, 2))
  cat_cpp <- .Call("magcomp_predict_category_cpp",
                   magnitude, combined_cpp, c(0.08, 0.10),
                   as.integer(1L), PACKAGE = "magComp")
  cat_r <- combined_cpp
  for(ch in 1:2){
    for(s in seq_len(dim(combined_cpp)[3L])){
      cat_r[,,s,ch] <- combined_cpp[,,s,ch] *
        matrix(magnitude[,s,ch], nrow = dim(combined_cpp)[1L],
               ncol = dim(combined_cpp)[2L])
    }
  }
  expect_equal(cat_cpp, cat_r, tolerance = 1e-12)
  expect_error(
    .Call("magcomp_predict_category_cpp",
          magnitude, combined_cpp, numeric(0L), as.integer(99L),
          PACKAGE = "magComp"),
    "prediction type"
  )
})

test_that("conditional factor default anchors can use later parent levels", {
  y <- matrix(0L, nrow = 4L, ncol = 3L)
  anchor_index <- magComp:::magcomp_default_factor_anchor_index(
    y,
    K = 3L,
    parent_index = c(1L, 1L, 2L, 2L),
    coef_stick_count = 4L
  )

  expect_equal(length(anchor_index), 3L)
  expect_false(anyNA(anchor_index))
  expect_true(any(anchor_index > 2L))
})

test_that("sampler supports multiple additive magnitude process terms", {
  set.seed(44)
  n_unit <- 10L
  z_unit <- seq(-1, 1, length.out = n_unit)
  trees <- do.call(rbind, lapply(seq_len(n_unit), function(i) {
    n_tree <- 3L + (i %% 4L)
    data.frame(
      plot = paste0("p", i),
      group = paste0("g", 1L + (i %% 3L)),
      block = paste0("b", 1L + (i %% 2L)),
      cls = rep(c("s", "m", "l"), length.out = n_tree),
      z = rep(z_unit[i], n_tree),
      x = rep(i, n_tree),
      ycoord = rep((i %% 5L) + 0.1 * i, n_tree),
      area = rep(0.1, n_tree)
    )
  }))

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    design = fixedArea(area)
  )

  fit <- magComp(dat, family = fixedAreaStems(), magnitude = ~z + iid(group) + z:iid(block) + nngp(x, ycoord,
    m = 3) + z:nngp(x, ycoord, m = 3), composition = ~1, smooth = magCompSmooth(intercept = "none",
    slope = "none", loading = "none"), control = magCompControl(n_samples = 6, n_burn = 2, n_thin = 1,
    verbose = FALSE, n_report = 100), starting = list(magnitude = list(terms = list(nngp_1 = list(covariance = list(phi = fixed(2))),
    nngp_2 = list(covariance = list(phi = fixed(3)))))))

  expect_equal(names(fit$backend$magnitude_process_effect_samples),
               c("iid_1", "iid_2", "nngp_1", "nngp_2"))
  expect_equal(names(lapply(fit$starting$magnitude$terms, "[[", "effects")),
               c("iid_1", "iid_2", "nngp_1", "nngp_2"))
  expect_equal(names(lapply(fit$starting$magnitude$terms, "[[", "variance")),
               c("iid_1", "iid_2", "nngp_1", "nngp_2"))
  expect_true(all(vapply(lapply(fit$starting$magnitude$terms, "[[", "effects"), function(x)
    all(is.finite(x)), logical(1L))))
  expect_true(all(unlist(lapply(fit$starting$magnitude$terms, "[[", "variance")) > 0))
  expect_equal(dim(fit$backend$magnitude_process_effect_samples$iid_1), c(3L, 4L))
  expect_equal(dim(fit$backend$magnitude_process_effect_samples$iid_2), c(2L, 4L))
  expect_equal(dim(fit$backend$magnitude_process_effect_samples$nngp_1),
               c(n_unit, 4L))
  expect_equal(dim(fit$backend$magnitude_process_effect_samples$nngp_2),
               c(n_unit, 4L))
  expect_equal(names(fit$backend$magnitude_process_variance_samples),
               c("iid_1", "iid_2", "nngp_1", "nngp_2"))
  expect_equal(names(fit$backend$magnitude_process_covariance_samples),
               c("nngp_1", "nngp_2"))
  expect_equal(fit$backend$nngp_graph_count, 1L)
  expect_true(all(is.finite(fit$backend$magnitude_process_effect_samples$iid_1)))
  expect_true(all(is.finite(fit$backend$magnitude_process_effect_samples$nngp_2)))

  newdata <- data.frame(
    plot = paste0("new", 1:3),
    group = c("g1", "g2", "g3"),
    block = c("b1", "b2", "b1"),
    cls = "s",
    z = c(-0.8, 0.0, 0.8),
    x = c(1.5, 4.5, 8.5),
    ycoord = c(1.2, 3.3, 4.4),
    area = 0.1
  )
  set.seed(77)
  pred <- predict(fit, newdata, type = "magnitude", summary = "samples")$magnitude
  expect_equal(dim(pred), c(3L, 4L))
  expect_true(all(is.finite(pred)))

  start <- magCompStarting(fit, use = "last")
  expect_equal(names(lapply(start$magnitude$terms, "[[", "effects")),
               names(fit$backend$magnitude_process_effect_samples))
  expect_equal(names(lapply(start$magnitude$terms, "[[", "variance")),
               c("iid_1", "iid_2", "nngp_1", "nngp_2"))
  expect_true(all(vapply(Filter(function(t) !is.null(t$covariance), start$magnitude$terms), function(t) inherits(t$covariance$phi, "magCompFixed"), logical(1))))

  tab <- summary(fit)
  expect_true(any(tab$component == "magnitude" & tab$term == "nngp_1" &
                    tab$parameter == "variance"))
  expect_true(any(tab$component == "magnitude" & tab$term == "nngp_2" &
                    tab$parameter == "variance"))
})

test_that("sampler carries multiple term names through starts, summaries, and prediction", {
  set.seed(45)
  n_unit <- 12L
  cls <- c("s", "m", "l")
  unit <- paste0("p", seq_len(n_unit))
  group <- rep(c("g1", "g2", "g3"), length.out = n_unit)
  block <- rep(c("b1", "b2"), length.out = n_unit)
  z <- seq(-1, 1, length.out = n_unit)
  x <- rep(seq_len(4L), each = 3L) + stats::runif(n_unit, -0.05, 0.05)
  ycoord <- rep(seq_len(3L), times = 4L) + stats::runif(n_unit, -0.05, 0.05)
  area <- rep(0.1, n_unit)

  counts <- matrix(0L, nrow = n_unit, ncol = length(cls),
                   dimnames = list(unit, cls))
  for(i in seq_len(n_unit)){
    total <- 12L + (i %% 5L)
    eta <- c(-0.15 + 0.25 * z[i], -0.65 - 0.10 * z[i])
    pi_i <- stats::pnorm(eta)
    p_i <- c(pi_i[1L],
             (1 - pi_i[1L]) * pi_i[2L],
             prod(1 - pi_i))
    counts[i,] <- as.vector(stats::rmultinom(1L, total, p_i))
  }

  trees <- data.frame()
  for(i in seq_len(n_unit)){
    trees <- rbind(trees, data.frame(
      plot = unit[i],
      cls = rep(cls, counts[i,]),
      group = group[i],
      block = block[i],
      z = z[i],
      x = x[i],
      ycoord = ycoord[i],
      area = area[i]
    ))
  }

  trees$cls <- factor(trees$cls, levels = cls)
  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, ycoord),
    ordered = c(cls = TRUE),
    design = fixedArea(area)
  )
  pri <- magCompPriors(magnitude = list(coefficient_variance = 25, iid_variance = c(2, 0.7), spatial_variance = c(3, 2)), composition = list(slope_variance = 25, intercept_variance = 4, slope_smoothing_variance = 3, intercept_smoothing_variance = 2, loading_variance = 1.5))
  loading_row_names <- cls[-length(cls)]
  start <- list(magnitude = list(coefficients = c("(Intercept)" = log(140), z = 0.1), terms = list(iid_1 = list(effects = c(g1 = 0.02, g2 = -0.03, g3 = 0.01)), iid_2 = list(effects = c(b1 = 0.01, b2 = -0.01)), gp_1 = list(effects = rep(0, n_unit), variance = 0.35, covariance = list(phi = 1.1)), nngp_1 = list(effects = rep(0, n_unit), variance = 0.45, covariance = list(phi = 1.4)))), composition = list(blocks = list(cls = list(coefficients = matrix(c(-0.2, -0.7, 0.05, -0.03),
                   nrow = 2L,
                   dimnames = list(loading_row_names,
                                   c("(Intercept)", "z"))), terms = list(factor_iid_1 = list(loadings = matrix(c(0.10, 1), nrow = 2L,
                              dimnames = list(loading_row_names, "fiid")), scores = matrix(c(0.02, -0.01, 0.01), nrow = 3L,
                              dimnames = list(c("g1", "g2", "g3"), "fiid"))), factor_gp_1 = list(loadings = matrix(c(0.08, 1), nrow = 2L,
                             dimnames = list(loading_row_names, "fgp")), scores = matrix(rep(0, n_unit), nrow = n_unit,
                             dimnames = list(unit, "fgp")), covariance = list(phi = 1.2)), factor_nngp_1 = list(loadings = matrix(c(0.07, 1), nrow = 2L,
                               dimnames = list(loading_row_names, "fnngp")), scores = matrix(rep(0, n_unit), nrow = n_unit,
                               dimnames = list(unit, "fnngp")), covariance = list(phi = 1.3)))))))

  mag_formula <- ~ z + iid(group) + z:iid(block) +
    gp(x, ycoord) +
    z:nngp(x, ycoord, m = 4)
  comp_formula <- ~ z +
    factor_iid(group, K = 1, anchors = "m", factor_names = "fiid") +
    z:factor_gp(x, ycoord, K = 1, anchors = "m", factor_names = "fgp") +
    factor_nngp(x, ycoord, K = 1, m = 4, anchors = "m", factor_names = "fnngp")

  fit <- magComp(
    dat,
    family = fixedAreaStems(),
    magnitude = mag_formula,
    composition = comp_formula,
    priors = pri,
    starting = start,
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none"),
    control = magCompControl(n_samples = 6,
                             n_burn = 3,
                             n_thin = 1,
                             n_chains = 2,
                             chain_seeds = c(451, 452),
                             composition_link = "probit",
                             n_report = 100,
                             verbose = FALSE)
  )

  expect_true(fit$fitted)
  expect_equal(fit$control$composition_link, "probit")
  expect_equal(as.character(fit$backend$composition_link), "probit")
  expect_true(fit$backend$profile_pg[["probit_composition_draws"]] > 0)
  expect_equal(fit$backend$profile_pg[["composition_draws"]], 0)
  expect_equal(fit$priors$magnitude$coefficient_variance, 25)
  expect_equal(fit$priors$magnitude$spatial_variance, c(shape = 3, scale = 2))
  expect_equal(names(fit$backend$magnitude_process_effect_samples),
               c("iid_1", "iid_2", "gp_1", "nngp_1"))
  expect_equal(names(fit$backend$magnitude_process_variance_samples),
               c("iid_1", "iid_2", "gp_1", "nngp_1"))
  expect_equal(names(fit$backend$magnitude_process_covariance_samples),
               c("gp_1", "nngp_1"))
  expect_equal(names(fit$backend$composition_block_score_samples),
               c("cls:factor_iid_1", "cls:factor_gp_1",
                 "cls:factor_nngp_1"))
  expect_equal(names(fit$backend$composition_block_covariance_samples),
               c("cls:factor_gp_1", "cls:factor_nngp_1"))
  expect_equal(colnames(fit$backend$composition_block_loading_samples[["cls:factor_iid_1"]]),
               "fiid")
  expect_equal(colnames(fit$backend$composition_block_loading_samples[["cls:factor_gp_1"]]),
               "fgp")
  expect_equal(colnames(fit$backend$composition_block_loading_samples[["cls:factor_nngp_1"]]),
               "fnngp")
  expect_loading_constraints(
    fit$starting$composition$blocks$cls$terms$factor_iid_1$loadings, "m", "fiid"
  )
  expect_loading_constraints(
    fit$starting$composition$blocks$cls$terms$factor_gp_1$loadings, "m", "fgp"
  )
  expect_loading_constraints(
    fit$starting$composition$blocks$cls$terms$factor_nngp_1$loadings, "m", "fnngp"
  )
  expect_loading_constraints(
    fit$backend$composition_block_loading_samples[["cls:factor_iid_1"]],
    "m", "fiid"
  )
  expect_loading_constraints(
    fit$backend$composition_block_loading_samples[["cls:factor_gp_1"]],
    "m", "fgp"
  )
  expect_loading_constraints(
    fit$backend$composition_block_loading_samples[["cls:factor_nngp_1"]],
    "m", "fnngp"
  )
  expect_equal(dim(fit$backend$magnitude_coefficient_samples), c(2L, 3L, 2L))
  expect_equal(dim(fit$backend$magnitude_process_effect_samples$gp_1),
               c(n_unit, 3L, 2L))
  expect_equal(dim(fit$backend$composition_block_score_samples[["cls:factor_gp_1"]]),
               c(n_unit, 1L, 3L, 2L))

  summary <- summary(fit)
  expect_true(all(summary$n_chains == 2L))
  expect_true(any(summary$component == "magnitude" & summary$term == "gp_1" &
                    summary$parameter == "variance"))
  expect_true(any(summary$component == "magnitude" & summary$term == "nngp_1" &
                    summary$parameter == "phi"))
  expect_true(any(summary$component == "composition" & summary$block == "cls" & summary$term == "factor_gp_1" &
                    summary$parameter == "phi"))
  expect_true(any(summary$component == "composition" & summary$block == "cls" & summary$term == "factor_nngp_1" &
                    summary$parameter == "phi"))
  expect_true(any(summary$parameter == "loadings" &
                    summary$block == "cls" & summary$term == "factor_gp_1" & summary$index == "m:fgp"))

  last_start <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(names(lapply(last_start$magnitude$terms, "[[", "effects")),
               c("iid_1", "iid_2", "gp_1", "nngp_1"))
  expect_equal(names(lapply(last_start$composition$blocks$cls$terms, "[[", "loadings")),
               c("factor_iid_1", "factor_gp_1", "factor_nngp_1"))
  expect_equal(names(lapply(last_start$composition$blocks$cls$terms, "[[", "scores")),
               c("factor_iid_1", "factor_gp_1", "factor_nngp_1"))
  expect_equal(names(Filter(Negate(is.null), lapply(last_start$composition$blocks$cls$terms, "[[", "covariance"))),
               c("factor_gp_1", "factor_nngp_1"))

  newdata <- data.frame(
    plot = paste0("new", seq_len(4L)),
    group = c("g1", "g2", "g3", "g1"),
    block = c("b1", "b2", "b1", "b2"),
    z = c(-0.75, -0.1, 0.35, 0.9),
    x = c(1.2, 2.5, 3.2, 4.1),
    ycoord = c(1.1, 2.0, 2.8, 1.7),
    area = c(0.05, 0.08, 0.10, 0.12),
    stand = c("A", "A", "B", "B")
  )
  set.seed(46)
  pred_mag <- predict(fit, newdata, type = "magnitude",
                      summary = "samples")$magnitude
  set.seed(46)
  pred_prop <- predict(fit, newdata, type = "proportion",
                       summary = "samples")$proportion
  set.seed(46)
  pred_count <- predict(fit, newdata, type = "count",
                        summary = "samples",
                        area = "area")$count
  expect_equal(dim(pred_mag), c(4L, 3L, 2L))
  expect_equal(dim(pred_prop), c(4L, length(cls), 3L, 2L))
  expect_equal(dim(pred_count), c(4L, length(cls), 3L, 2L))
  expect_equal(as.numeric(apply(pred_prop, c(1L, 3L, 4L), sum)),
               rep(1, 4L * 3L * 2L),
               tolerance = 1e-10)

  fit_bad <- fit
  names(fit_bad$backend$composition_block_loading_samples)[1L] <- "bad_key"
  set.seed(46)
  expect_error(
    predict(fit_bad, newdata, type = "proportion",
            summary = "samples")$proportion,
    "posterior loading samples"
  )
  expect_true(all(is.finite(pred_count)))

  set.seed(46)
  stand_samples <- standSummary(predict(fit, newdata, summary = "samples"),
    data = newdata, unit = plot, group = stand, represented_area = area,
    type = "total", summary = "samples")$total
  expect_equal(dim(stand_samples), c(2L, length(cls), 3L, 2L))
  expect_true(all(is.finite(stand_samples)))
})

test_that("sampler estimates negative-binomial magnitude size", {
  set.seed(812)
  n <- 10
  cls <- c("small", "large")
  area <- rep(0.08, n)
  z <- seq(-1, 1, length.out = n)
  stand <- rep(c("a", "b"), length.out = n)
  total <- pmax(1L, stats::rnbinom(n, size = 4,
                                   mu = area * exp(log(160) + 0.3 * z)))
  prob <- plogis(-0.2 + 0.5 * z)
  y <- cbind(
    small = stats::rbinom(n, total, prob),
    large = total
  )
  y[,2L] <- y[,2L] - y[,1L]

  trees <- data.frame()
  for(i in seq_len(n)){
    trees <- rbind(trees, data.frame(
      plot = paste0("p", i),
      cls = rep(cls, y[i,]),
      area = area[i],
      z = z[i],
      stand = stand[i]
    ))
  }

  trees$cls <- factor(trees$cls, levels = sort(cls))
  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    ordered = c(cls = TRUE),
    design = fixedArea(area)
  )

  fit <- magComp(
    dat,
    magnitude = ~ z,
    composition = ~ z,
    priors = magCompPriors(magnitude = list(size = c(2, 0.2))),
    control = magCompControl(
      n_samples = 16,
      n_burn = 8,
      n_thin = 1,
      n_chains = 2,
      chain_seeds = c(701, 702),
      nb_size = 10,
      verbose = FALSE
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none")
  )

  expect_equal(dim(fit$backend$nb_size_samples), c(1L, 8L, 2L))
  expect_true(all(is.finite(fit$backend$nb_size_samples)))
  expect_true(all(fit$backend$nb_size_samples > 0))
  expect_equal(dim(fit$backend$nb_size_acceptance), c(1L, 2L))
  expect_true(all(is.finite(fit$backend$nb_size_acceptance)))
  expect_true(all(fit$backend$nb_size_acceptance >= 0))
  expect_true(all(fit$backend$nb_size_acceptance <= 1))
  expect_gt(length(unique(as.numeric(fit$backend$nb_size_samples))), 1L)
  expect_false(all(as.numeric(fit$backend$nb_size_samples) == 10))
  expect_equal(rownames(fit$backend$nb_size_samples), "size")

  summary <- summary(fit)
  expect_true(any(summary$component == "magnitude" &
                    summary$parameter == "size"))

  start <- magCompStarting(fit, use = "last", chain = 2)
  expect_true(is.finite(start$magnitude$size))
  expect_true(start$magnitude$size > 0)

  set.seed(91)
  y_rep <- fitted(fit, target = "observation", type = "count", summary = "samples")$count
  expect_equal(dim(y_rep), c(n, length(cls), 8L, 2L))
  expect_true(all(y_rep >= 0))

  fit_other_scalar_size <- fit
  fit_other_scalar_size$backend$nb_size <- 1000
  set.seed(91)
  expect_identical(fitted(fit_other_scalar_size, target = "observation", type = "count", summary = "samples")$count, y_rep)

  fit_bad_mean <- fit
  fit_bad_mean$backend$magnitude_samples[1L,1L,1L] <- Inf
  set.seed(91)
  expect_error(fitted(fit_bad_mean, target = "observation", type = "count", summary = "samples")$count,
               "predictive category means must be finite and nonnegative")

  fit_fixed_size <- fit
  fit_fixed_size$backend$nb_size_samples <- NULL
  fit_fixed_size$backend$nb_size <- 4
  predict_nb <- getFromNamespace("magcomp_posterior_predict_nb", "magComp")
  expected_fixed_rep <- local({
    set.seed(92)
    predict_nb(
      fitted(fit_fixed_size, type = "count", summary = "samples")$count,
      fit_fixed_size$backend$nb_size
    )
  })
  set.seed(92)
  expect_equal(fitted(fit_fixed_size, target = "observation", type = "count", summary = "samples")$count,
               expected_fixed_rep)

  expect_error(
    magComp(
      dat,
      magnitude = ~ z,
      composition = ~ z,
      starting = list(magnitude = list(size = 4)),
      priors = magCompPriors(magnitude = list(size = NULL)),
      control = magCompControl(run_sampler = FALSE, verbose = FALSE)
    ),
    "unknown setting in starting\\$magnitude: size"
  )

  expect_error(
    magComp(
      dat,
      magnitude = ~ z + iid(stand),
      priors = magCompPriors(magnitude = list(size = c(2, 0.2))),
      control = magCompControl(n_samples = 3,
                               n_burn = 1,
                               verbose = FALSE)
    ),
    "estimating nb_size"
  )
})
