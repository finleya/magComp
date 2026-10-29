test_that("each magnitude iid term samples its own conditional variance", {
  rows <- expand.grid(plot = paste0("p", 1:12), cls = c("a", "b"))
  rows$cls <- factor(rows$cls, levels = c("a", "b"))
  rows$stand <- rep(rep(paste0("s", 1:4), each = 3), 2)
  rows$block <- rep(rep(paste0("b", 1:3), 4), 2)
  rows$x <- rep(c(-1, 0.5, 1, 0.5, 1, -1, 1, -1, 0.5, -1, 1, 0.5), 2)
  rows$count <- c(4L, 8L, 12L, 16L, 6L, 20L, 3L, 14L, 25L, 9L, 30L, 5L,
                  2L, 4L, 8L, 7L, 3L, 8L, 2L, 6L, 10L, 3L, 12L, 2L)
  rows$area <- 0.1
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                    design = fixedArea(area))
  fit <- magComp(
    dat, magnitude = ~ x + iid(stand) + x:iid(block),
    priors = magCompPriors(magnitude = list(iid_variance = c(3, 0.8),
      terms = list(iid_2 = list(variance = c(5, 1.6))))),
    starting = list(magnitude = list(terms = list(iid_1 = list(variance = 0.2), iid_2 = list(variance = 2)))),
    control = magCompControl(n_samples = 900, n_burn = 300,
                             n_chains = 2, chain_seeds = c(1931, 1932),
                             nb_size = 20, verbose = FALSE)
  )
  variances <- fit$backend$magnitude_process_variance_samples
  expect_named(variances, c("iid_1", "iid_2"))
  expect_null(fit$backend$magnitude_process_covariance_samples)
  expect_equal(fit$starting$magnitude$terms$iid_1$variance, 0.2)
  expect_equal(fit$starting$magnitude$terms$iid_2$variance, 2)
  expect_false(identical(variances$iid_1, variances$iid_2))
  for(term in fit$magnitude$process_terms){
    sigma <- variances[[term$name]]
    expect_equal(dim(sigma), c(1L, 600L, 2L))
    expect_true(all(is.finite(sigma) & sigma > 0))
    expect_gt(sd(as.numeric(sigma)), 0)
    effects <- fit$backend$magnitude_process_effect_samples[[term$name]]
    prior <- term$variance_prior
    rate <- prior[2L] + apply(effects^2, c(2L, 3L), sum) / 2
    # If the saved variance is the correct Gibbs draw given these effects,
    # its conditional inverse-gamma CDF is uniform regardless of the data.
    pit <- pgamma(as.numeric(rate) / as.numeric(sigma),
                  shape = prior[1L] + term$n_group / 2)
    expect_lt(abs(mean(pit) - 0.5), 0.04)
    expect_lt(abs(mean(pit^2) - 1 / 3), 0.04)
  }
  tab <- summary(fit)
  expect_setequal(tab$term[tab$parameter == "variance"],
                  c("iid_1", "iid_2"))
  start <- magCompStarting(fit, use = "last", chain = 2)
  average <- magCompStarting(fit, use = "mean")
  for(name in names(variances)){
    expect_equal(as.numeric(start$magnitude$terms[[name]]$variance),
                  variances[[name]][1L,600L,2L])
    expect_equal(as.numeric(average$magnitude$terms[[name]]$variance),
                  mean(variances[[name]]))
  }
  restart <- magComp(
    dat, magnitude = ~ x + iid(stand) + x:iid(block), starting = start,
    priors = fit$priors,
    control = magCompControl(n_samples = 3, n_burn = 1, nb_size = 20,
                             chain_seeds = 1933, verbose = FALSE)
  )
  expect_equal(lapply(restart$starting$magnitude$terms, "[[", "variance"),
                lapply(start$magnitude$terms, "[[", "variance"))
  expect_equal(dim(restart$backend$magnitude_process_variance_samples$iid_1),
                c(1L, 2L))

  newdata <- data.frame(stand = "s1", block = "new", x = c(-1, 1))
  fixed <- mc_regression_magnitude_samples(fit, newdata)
  set.seed(1934)
  full <- predict(fit, newdata, type = "magnitude", summary = "samples")$magnitude
  stand_effect <- fit$backend$magnitude_process_effect_samples$iid_1[1L,,]
  slope1 <- -(log(full[1L,,]) - log(fixed[1L,,]) - stand_effect)
  slope2 <- log(full[2L,,]) - log(fixed[2L,,]) - stand_effect
  expect_equal(slope1, slope2, tolerance = 1e-12)
  changed <- fit
  changed$backend$magnitude_process_variance_samples$iid_2 <- 4 * variances$iid_2
  set.seed(1934)
  full_changed <- predict(changed, newdata, type = "magnitude", summary = "samples")$magnitude
  expect_equal(log(full_changed[2L,,]) - log(fixed[2L,,]) - stand_effect,
                2 * slope2, tolerance = 1e-12)

  newdata <- data.frame(stand = "new", block = "b1", x = 0)
  fixed <- mc_regression_magnitude_samples(fit, newdata)
  set.seed(1935)
  full <- predict(fit, newdata, type = "magnitude", summary = "samples")$magnitude
  changed <- fit
  changed$backend$magnitude_process_variance_samples$iid_1 <- 9 * variances$iid_1
  set.seed(1935)
  full_changed <- predict(changed, newdata, type = "magnitude", summary = "samples")$magnitude
  expect_equal(log(full_changed) - log(fixed), 3 * (log(full) - log(fixed)),
                tolerance = 1e-12)
})

test_that("iid variance priors and starts have automatic defaults and validation", {
  expect_equal(magCompPriors()$magnitude$iid_variance, c(shape = 2, scale = 1))
  for(prior in list(1, c(0, 1), c(1, -1), c(2, Inf), c(NA, 1)))
    expect_error(magCompPriors(magnitude = list(iid_variance = prior)), "iid_variance")
  rows <- data.frame(plot = rep(c("p1", "p2"), each = 2),
                     cls = rep(c("a", "b"), 2), count = c(2L, 3L, 4L, 5L),
                     area = 0.1, stand = rep(c("s1", "s2"), each = 2))
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                    design = fixedArea(area))
  control <- magCompControl(n_samples = 3, n_burn = 1, chain_seeds = 1936, verbose = FALSE)
  fit <- magComp(dat, magnitude = ~ iid(stand), control = control)
  expect_gt(fit$starting$magnitude$terms$iid_1$variance, 0)
  for(value in c(0, -1, Inf, NA_real_))
    expect_error(magComp(dat, magnitude = ~ iid(stand), control = control,
                          starting = list(magnitude = list(terms = list(iid_1 = list(variance = value))))),
                   "positive and finite|finite numeric")
})
