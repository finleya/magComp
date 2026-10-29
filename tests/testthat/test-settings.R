settings_data <- function(two = FALSE, zero = FALSE)
{
  rows <- if(two)
    expand.grid(plot = paste0("p", 1:8), species = c("a", "b", "c"), dbh = c("small", "large"))
  else expand.grid(plot = paste0("p", 1:8), cls = c("a", "b", "c"))
  i <- match(rows$plot, paste0("p", 1:8))
  rows$count <- if(zero) as.integer(i == 1 & rows$species == "a" & rows$dbh == "small") else as.integer(2 + i %% 4)
  rows$area <- 1
  rows$stand <- paste0("s", (i - 1L) %/% 2 + 1L)
  rows$x <- i / 8
  rows$y <- (i %% 3) / 3
  if(two)
    magCompData(rows, unit = plot, categories = c(species, dbh),
                ordered = c(species = FALSE, dbh = TRUE),
                tree_count = count, design = fixedArea(area))
  else
    magCompData(rows, unit = plot, category = cls, tree_count = count,
                design = fixedArea(area))
}

test_that("settings templates resolve named priors, starts, and estimated decays", {
  dat <- settings_data()
  priors <- magCompPriors(
    magnitude = list(coefficient_variance = 20,
      terms = list(stand = list(variance = c(3, 0.4)),
                   site = list(variance = c(4, 0.6), covariance = list(phi = uniform(0.2, 5))))),
    composition = list(slope_variance = 7, blocks = list(cls = list(
      intercept_variance = 2, terms = list(shape = list(
        loading_variance = 0.3, covariance = list(phi = uniform(0.5, 6))))))))
  mag <- ~ x + iid(stand, name = "stand") + gp(x, y, name = "site")
  comp <- ~ x + factor_nngp(x, y, K = 2, m = 3, name = "shape",
                            anchors = c("a", "b"), factor_names = c("f1", "f2"))
  model <- magComp(dat, magnitude = mag, composition = comp, priors = priors,
    control = magCompControl(run_sampler = FALSE))
  settings <- magCompSettings(model)
  expect_s3_class(settings, "magCompSettings")
  expect_equal(settings$terms$term, c("stand", "site", "shape"))
  expect_true(all(subset(settings$parameters, parameter == "phi")$status == "estimated"))
  expect_equal(settings$starting, magCompStarting(model))
  expect_output(print(settings), "anchor diagonal")
  expect_equal(model$magnitude$process_terms[[1L]]$variance_prior, c(shape = 3, scale = 0.4))
  expect_equal(model$magnitude$process_terms[[2L]]$variance_prior, c(shape = 4, scale = 0.6))
  expect_equal(model$config$composition_blocks$cls$coefficient_priors$intercept_variance, 2)
  expect_equal(model$config$composition_blocks$cls$coefficient_priors$slope_variance, 7)
  expect_equal(model$config$priors$coefficient_variance, 20)
  expect_equal(model$composition$process_terms[[1L]]$loading_variance, 0.3)

  starts <- settings$starting
  starts$magnitude$coefficients <- c(x = 0.1, "(Intercept)" = 2)
  starts$magnitude$terms$stand$variance <- 0.6
  starts$magnitude$terms$site$covariance$phi <- 1.2
  starts$composition$blocks$cls$terms$shape$covariance$phi <- c(f2 = 2.5, f1 = 1.5)
  control <- magCompControl(n_samples = 60, n_burn = 30, n_chains = 2,
                            chain_seeds = c(703, 704), nb_size = 20,
                            verbose = FALSE)
  fit <- magComp(dat, magnitude = mag, composition = comp, priors = settings$priors,
                  starting = starts, control = control)
  expect_equal(fit$starting$magnitude$coefficients, c("(Intercept)" = 2, x = 0.1))
  expect_equal(fit$starting$composition$blocks$cls$terms$shape$covariance$phi, c(f1 = 1.5, f2 = 2.5))
  mag_phi <- fit$backend$magnitude_process_covariance_samples$site
  comp_phi <- fit$backend$composition_block_covariance_samples$cls
  expect_gt(length(unique(as.numeric(mag_phi))), 1)
  expect_true(all(mag_phi > 0.2 & mag_phi < 5))
  for(k in 1:2) expect_gt(length(unique(as.numeric(comp_phi[k,,]))), 1)
  expect_true(all(comp_phi > 0.5 & comp_phi < 6))
  last <- magCompStarting(fit, use = "last", chain = 2)
  mean <- magCompStarting(fit, use = "mean")
  expect_equal(last$magnitude$terms$site$covariance$phi, unname(mag_phi[1L,30L,2L]))
  expect_equal(mean$magnitude$terms$site$covariance$phi, mean(mag_phi))
  control$n_samples <- 32L
  resumed <- magComp(dat, magnitude = mag, composition = comp, priors = settings$priors,
                      starting = last, control = control)
  expect_equal(resumed$starting, last, tolerance = 1e-12)
  set.seed(708)
  pred <- predict(resumed, data.frame(x = 0.2, y = 0.2, stand = "new"),
                    type = "magnitude", summary = "samples")$magnitude
  expect_true(all(is.finite(pred) & pred > 0))
})

test_that("decay defaults cover GP and NNGP and fixed decays stay fixed", {
  dat <- settings_data()
  for(kind in c("gp", "nngp")){
    mag <- as.formula(paste0("~ ", kind, "(x, y)"))
    comp <- as.formula(paste0("~ factor_", kind, "(x, y, K = 1)"))
    model <- magComp(dat, magnitude = mag, composition = comp,
    control = magCompControl(run_sampler = FALSE))
    expect_false(any(model$magnitude$process_terms[[1L]]$covariance_fixed))
    expect_false(any(model$composition$process_terms[[1L]]$covariance_fixed))
    starts <- magCompStarting(model)
    starts$magnitude$terms[[1L]]$covariance$phi <- fixed(2)
    starts$composition$blocks$cls$terms[[1L]]$covariance$phi <- fixed(3)
    fixed_fit <- magComp(dat, magnitude = mag, composition = comp, starting = starts,
      control = magCompControl(n_samples = 5, n_burn = 2, verbose = FALSE))
    settings <- magCompSettings(fixed_fit)
    expect_true(all(subset(settings$parameters, parameter == "phi")$status == "fixed"))
    expect_s3_class(settings$starting$magnitude$terms[[1L]]$covariance$phi, "magCompFixed")
    expect_true(all(fixed_fit$backend$magnitude_process_covariance_samples[[1L]] == 2))
    expect_true(all(fixed_fit$backend$composition_block_covariance_samples[[1L]] == 3))
    expect_s3_class(magCompStarting(fixed_fit, use = "last")$magnitude$terms[[1L]]$covariance$phi,
                    "magCompFixed")
  }
})

test_that("estimated decay starts are strictly inside default and custom bounds", {
  rows <- data.frame(plot = rep(1:2, 2), cls = rep(c("a", "b"), each = 2),
    count = 1, area = 1, x = rep(0:1, 2), y = 0)
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
    design = fixedArea(area))
  mag <- ~ gp(x, y, name = "space")
  comp <- ~ factor_nngp(x, y, K = 1, name = "space")
  control <- magCompControl(run_sampler = FALSE)
  automatic <- magComp(dat, magnitude = mag, composition = comp, control = control)
  for(term in c(automatic$magnitude$process_terms, automatic$composition$process_terms)){
    expect_gt(unname(term$covariance_values[1,1]), term$covariance_lower[1,1])
    expect_lt(unname(term$covariance_values[1,1]), term$covariance_upper[1,1])
  }

  priors <- magCompPriors(
    magnitude = list(terms = list(space = list(covariance = list(phi = uniform(1, 2))))),
    composition = list(blocks = list(cls = list(terms = list(space = list(covariance = list(phi = uniform(1, 2))))))))
  model <- magComp(dat, magnitude = mag, composition = comp, priors = priors, control = control)
  starts <- magCompStarting(model)
  expect_equal(starts$magnitude$terms$space$covariance$phi, 1.5)
  expect_equal(unname(starts$composition$blocks$cls$terms$space$covariance$phi), 1.5)
  for(boundary in c(1, 2)){
    bad <- starts
    bad$magnitude$terms$space$covariance$phi <- boundary
    expect_error(magComp(dat, magnitude = mag, composition = comp, priors = priors,
      starting = bad, control = control), "strictly inside.*magnitude")
    bad <- starts
    bad$composition$blocks$cls$terms$space$covariance$phi[] <- boundary
    expect_error(magComp(dat, magnitude = mag, composition = comp, priors = priors,
      starting = bad, control = control), "strictly inside.*composition")
  }
})

test_that("rounded decay proposals retain the current state in every component", {
  dat <- settings_data(two = TRUE)
  control <- magCompControl(n_samples = 4, n_burn = 0, verbose = FALSE,
    metropolis = list(adapt = FALSE))
  for(kind in c("gp", "nngp")){
    mag <- as.formula(paste0("~ ", kind, "(x, y)"))
    parent <- as.formula(paste0("species ~ factor_", kind,
      "(x, y, K = 2)"))
    child <- as.formula(paste0("dbh | species ~ factor_", kind,
      "(x, y, K = 1)"))
    composition <- magCompComposition(parent, child)
    set.seed(106)
    tune <- list(covariance = list(phi = 1e150))
    tuning <- list(magnitude = list(terms = setNames(list(tune), paste0(kind, "_1"))),
      composition = list(blocks = list(
        species = list(terms = setNames(list(tune), paste0("factor_", kind, "_1"))),
        "dbh|species" = list(terms = setNames(list(tune), paste0("factor_", kind, "_1"))))))
    fit <- magComp(dat, magnitude = mag, composition = composition, control = control, tuning = tuning)
    start <- fit$starting
    expect_true(all(fit$backend$magnitude_process_covariance_samples[[1]] ==
      start$magnitude$terms[[1]]$covariance$phi))
    for(block in c("species", "dbh|species")){
      draws <- fit$backend$composition_block_covariance_samples[[block]]
      initial <- start$composition$blocks[[block]]$terms[[1]]$covariance$phi
      expect_equal(as.numeric(draws), rep(as.numeric(initial), control$n_samples))
    }
  }
})

test_that("custom settings reject unused names and malformed arrays", {
  dat <- settings_data()
  model <- magComp(dat, magnitude = ~ iid(stand, name = "stand") + gp(x, y, name = "site"),
                    composition = ~ factor_iid(stand, K = 1, name = "shape"),
    control = magCompControl(run_sampler = FALSE))
  expect_equal(magCompPriors(magnitude = list(iid_variance = c(scale = 1, shape = 3)))$magnitude$iid_variance,
                c(shape = 3, scale = 1))
  expect_error(magCompPriors(magnitude = list(iid_variance = c(shape = 3, rate = 1))),
                "names must be")
  invalid <- magCompStarting(model)
  invalid$composition$blocks$cls$terms$shape$loadings[,] <- 0
  expect_error(magComp(dat, magnitude = model$magnitude$formula,
    composition = model$composition$formula, starting = invalid,
    control = magCompControl(run_sampler = FALSE)), "anchor constraints")
  expect_error(magCompPriors(magnitude = list(typo = 10)), "unknown setting")
  expect_error(magCompPriors(magnitude = list(terms = list(
    site = list(covariance = list(phi = uniform(2, 1)))))), "lower < upper")
  expect_error(magComp(dat, magnitude = ~ gp(x, y, typo = 2),
    control = magCompControl(run_sampler = FALSE)), "unknown setting")
  expect_error(magComp(dat, magnitude = ~ iid(stand, name = "same") +
                         x:iid(stand, name = "same"),
    control = magCompControl(run_sampler = FALSE)), "unique")
  expect_error(magComp(dat, priors = magCompPriors(magnitude = list(
    terms = list(typo = list(variance = c(2, 1))))),
    control = magCompControl(run_sampler = FALSE)), "unknown setting")
  expect_error(magComp(dat, priors = magCompPriors(composition = list(
    blocks = list(typo = list(intercept_variance = 1)))),
    control = magCompControl(run_sampler = FALSE)), "unknown setting")
  expect_error(magComp(dat, magnitude = ~ gp(x, y, name = "site"),
    priors = magCompPriors(magnitude = list(terms = list(site = list(
      covariance = list(phi = uniform(1, 4)))))),
    starting = list(magnitude = list(terms = list(site = list(covariance = list(phi = fixed(2)))))),
    control = magCompControl(run_sampler = FALSE)), "prior supplied for fixed")
  expect_error(magComp(dat, starting = list(typo = 0),
    control = magCompControl(run_sampler = FALSE)), "unknown setting")
  expect_error(magComp(dat, starting = list(magnitude = list(coefficients = c(1, 2))),
    control = magCompControl(run_sampler = FALSE)),
                "finite numeric")
  bad <- magCompStarting(model)
  bad$composition$blocks$cls$terms$shape$loadings <- matrix(1, 1, 2)
  expect_error(magComp(dat, magnitude = model$magnitude$formula,
    composition = model$composition$formula, starting = bad,
    control = magCompControl(run_sampler = FALSE)), "dimensions")
  bad <- magCompStarting(model)
  bad$magnitude$terms$site$covariance$phi <- 1e10
  expect_error(magComp(dat, magnitude = model$magnitude$formula,
    composition = model$composition$formula, starting = bad,
    control = magCompControl(run_sampler = FALSE)), "strictly inside")
  bad <- magCompStarting(model)
  names(bad$magnitude$coefficients) <- "typo"
  expect_error(magComp(dat, magnitude = model$magnitude$formula,
    composition = model$composition$formula, starting = bad,
    control = magCompControl(run_sampler = FALSE)), "names must match")

})

test_that("block and term priors reach the composition sampler", {
  dat <- settings_data(two = TRUE, zero = TRUE)
  priors <- magCompPriors(composition = list(intercept_variance = 2,
    blocks = list(species = list(terms = list(shape = list(loading_variance = 0.3))),
                  "dbh|species" = list(intercept_variance = 5,
                    terms = list(shape = list(loading_variance = 3))))))
  composition <- magCompComposition(
    species ~ factor_iid(stand, K = 1, name = "shape", anchors = "a"),
    dbh | species ~ factor_iid(stand, K = 1, name = "shape", anchors = "a:small"))
  # Sticks for species b and the child blocks of b/c have zero trials,
  # so their coefficient and free-loading draws follow the priors exactly.
  fit <- magComp(dat, composition = composition, priors = priors,
    control = magCompControl(n_samples = 810, n_burn = 10, nb_size = 20,
                              chain_seeds = 731, verbose = FALSE))
  a <- as.numeric(fit$backend$composition_block_coefficient_samples$species[2L,1L,])
  b <- as.numeric(fit$backend$composition_block_coefficient_samples[["dbh|species"]][2:3,1L,])
  expect_equal(mean(a^2), 2, tolerance = 0.12)
  expect_equal(mean(b^2), 5, tolerance = 0.12)
  la <- fit$backend$composition_block_loading_samples$species[2L,1L,]
  lb <- fit$backend$composition_block_loading_samples[["dbh|species"]][2L,1L,]
  # These are independent prior draws: check the exact chi-square law of
  # their squared norm, rather than a fixed relative Monte Carlo tolerance.
  for(value in list(list(draws = la, variance = 0.3), list(draws = lb, variance = 3))){
    pit <- pchisq(sum(value$draws^2) / value$variance, df = length(value$draws))
    expect_gt(pit, 1e-6)
    expect_lt(pit, 1 - 1e-6)
  }
  settings <- magCompSettings(fit)
  again <- magComp(dat, composition = composition, priors = settings$priors,
                    starting = settings$starting,
    control = magCompControl(run_sampler = FALSE))
  expect_equal(again$starting, settings$starting)
})

test_that("spatial variance overrides reach both GP and NNGP Gibbs updates", {
  dat <- settings_data()
  prior <- c(shape = 6, scale = 0.7)
  for(kind in c("gp", "nngp")){
    spatial <- if(kind == "gp") "gp(x, y, name = 'site')" else
      "nngp(x, y, m = 7, name = 'site')"
    formula <- as.formula(paste("~ iid(stand, name = 'stand') +", spatial))
    fit <- magComp(dat, magnitude = formula,
      starting = list(magnitude = list(terms = list(site = list(covariance = list(phi = fixed(2)))))),
      priors = magCompPriors(magnitude = list(terms = list(site = list(variance = prior)))),
      control = magCompControl(n_samples = 1100, n_burn = 100, nb_size = 20,
                                chain_seeds = 742, verbose = FALSE))
    term <- fit$magnitude$process_terms[[2L]]
    Q <- solve(exp(-2 * as.matrix(stats::dist(term$coords))))
    effects <- fit$backend$magnitude_process_effect_samples$site
    variance <- as.numeric(fit$backend$magnitude_process_variance_samples$site)
    rate <- prior[2L] + colSums(effects * (Q %*% effects)) / 2
    pit <- pgamma(rate / variance, shape = prior[1L] + nrow(effects) / 2)
    expect_lt(abs(mean(pit) - 0.5), 0.05)
    expect_lt(abs(mean(pit^2) - 1/3), 0.05)
  }
})

test_that("fixed and estimated spatial terms can share a block and its templates", {
  dat <- settings_data()
  magnitude <- ~ gp(x, y, name = "fixed") +
    nngp(x, y, m = 3, name = "estimated")
  composition <- ~ factor_gp(x, y, K = 1, name = "fixed") +
    factor_nngp(x, y, K = 1, m = 3, name = "estimated")
  model <- magComp(dat, magnitude = magnitude, composition = composition,
    starting = list(magnitude = list(terms = list(fixed = list(covariance = list(phi = fixed(2))))),
      composition = list(blocks = list(cls = list(terms = list(fixed = list(covariance = list(phi = fixed(3)))))))),
    control = magCompControl(run_sampler = FALSE))
  settings <- magCompSettings(model)
  fit <- magComp(dat, magnitude = magnitude, composition = composition,
    starting = settings$starting, priors = settings$priors,
    control = magCompControl(n_samples = 10, n_burn = 5, chain_seeds = 760,
                              verbose = FALSE))
  expect_true(all(fit$backend$magnitude_process_covariance_samples$fixed == 2))
  expect_true(all(fit$backend$composition_block_covariance_samples[["cls:fixed"]] == 3))
  starts <- magCompStarting(fit, use = "last")
  expect_s3_class(starts$magnitude$terms$fixed$covariance$phi, "magCompFixed")
  expect_s3_class(starts$composition$blocks$cls$terms$fixed$covariance$phi, "magCompFixed")
  expect_true(is.finite(starts$magnitude$terms$estimated$covariance$phi))
  expect_true(is.finite(starts$composition$blocks$cls$terms$estimated$covariance$phi))
  prepared <- magComp(dat, magnitude = magnitude, composition = composition,
                       starting = starts, priors = settings$priors,
    control = magCompControl(run_sampler = FALSE))
  expect_equal(prepared$starting, starts)
})


test_that("public samples match starts and preserve block, term, and factor identities", {
  dat <- settings_data(two = TRUE)
  comp <- magCompComposition(
    species ~ x + factor_iid(stand, K = 1, name = "shape"),
    dbh | species ~ x + factor_gp(x, y, K = 2, name = "shape",
                                  factor_names = c("f1", "f2")))
  fit <- magComp(dat, magnitude = ~x + iid(stand, name = "stand") + gp(x, y, name = "site"), composition = comp,
    smooth = magCompSmooth(dbh = list(intercept = "rw1", slope = "rw2")), control = magCompControl(n_samples = 7,
        n_burn = 3, n_chains = 2, chain_seeds = c(801, 802), verbose = FALSE), starting = list(magnitude = list(terms = list(site = list(covariance = list(phi = fixed(2)))))))
  samples <- samples(fit)
  last <- magCompStarting(fit, use = "last", chain = 2)
  expect_named(samples, c("magnitude", "composition"))
  expect_equal(dim(samples$magnitude$coefficients), c(2L, 4L, 2L))
  expect_named(samples$magnitude$terms, c("stand", "site"))
  expect_null(samples$magnitude$size)
  expect_null(samples$magnitude$terms$site$covariance$phi)
  expect_equal(last$magnitude$coefficients, samples$magnitude$coefficients[,4L,2L])
  expect_equal(last$magnitude$terms$stand$variance,
                unname(samples$magnitude$terms$stand$variance[1L,4L,2L]))
  for(b in names(samples$composition$blocks)){
    part <- samples$composition$blocks[[b]]
    expect_equal(last$composition$blocks[[b]]$coefficients, part$coefficients[,,4L,2L])
    expect_equal(last$composition$blocks[[b]]$terms$shape$loadings,
                  matrix(part$terms$shape$loadings[,,4L,2L],
                    nrow = dim(part$terms$shape$loadings)[1L],
                    dimnames = dimnames(part$terms$shape$loadings)[1:2]))
  }
  child <- samples$composition$blocks[["dbh|species"]]$terms$shape
  expect_equal(rownames(child$covariance$phi), c("f1", "f2"))
  expect_equal(last$composition$blocks[["dbh|species"]]$terms$shape$covariance$phi, child$covariance$phi[,4L,2L])
  tab <- summary(fit)
  expect_true(all(c("component", "block", "term", "parameter", "index", "symbol") %in% names(tab)))
  expect_equal(subset(tab, term == "stand" & parameter == "variance")$symbol, "sigma_delta^2")
  expect_true(all(startsWith(subset(tab, block == "species" & parameter == "loadings")$symbol, "gamma[")))
  expect_true(all(startsWith(subset(tab, block == "dbh|species" & parameter == "loadings")$symbol, "lambda[")))
  settings <- magCompSettings(fit)
  expect_equal(subset(settings$parameters, term == "site" & parameter == "variance")$symbol, "sigma_v^2")
  expect_equal(subset(settings$parameters, parameter == "phi")$symbol, c("phi_v", "phi_w", "phi_w"))
  spatial <- spatialDiagnostics(fit)
  expect_setequal(subset(spatial$summary, block == "dbh|species" & quantity == "phi")$factor, c("f1", "f2"))
  expect_equal(nrow(subset(spatial$acceptance, block == "dbh|species")), 4L)
  expect_error(samples(magComp(dat, composition = comp,
    control = magCompControl(run_sampler = FALSE))), "fitted magComp")
})

test_that("size defaults follow magnitude terms and settings can be reused", {
  dat <- settings_data()
  control <- magCompControl(run_sampler = FALSE, verbose = FALSE)
  for(mag in list(~ x, ~ iid(stand), ~ gp(x, y), ~ nngp(x, y, m = 3))){
    model <- magComp(dat, magnitude = mag, control = control)
    settings <- magCompSettings(model)
    estimated <- !length(model$magnitude$process_terms)
    expect_identical(subset(settings$parameters, parameter == "size")$status,
                     if(estimated) "estimated" else "fixed")
    expect_equal(model$priors$magnitude$size,
                 if(estimated) c(shape = 1, rate = 0.01) else NULL)
    again <- magComp(dat, magnitude = mag, priors = settings$priors,
                     starting = settings$starting, control = control)
    expect_equal(again$config$priors, model$config$priors)
    expect_equal(again$starting, model$starting)
  }
  fixed <- magComp(dat, priors = magCompPriors(magnitude = list(size = NULL)),
    control = magCompControl(nb_size = 20, run_sampler = FALSE))
  expect_null(fixed$config$priors$nb_size)
  expect_equal(subset(magCompSettings(fixed)$parameters, parameter == "size")$initial, "20")
  expect_error(magComp(dat, magnitude = ~ iid(stand),
    priors = magCompPriors(magnitude = list(size = c(1, 0.01))),
    control = control), "estimating nb_size")
})

test_that("sample extraction handles automatic size estimation and several factor terms", {
  fit <- magComp(settings_data(),
    composition = ~ factor_iid(stand, K = 1, name = "group") +
      factor_gp(x, y, K = 2, name = "space", factor_names = c("f1", "f2")),
    starting = list(magnitude = list(size = 10)),
    control = magCompControl(n_samples = 6, n_burn = 3, chain_seeds = 809,
                              verbose = FALSE))
  samples <- samples(fit)
  expect_equal(fit$priors$magnitude$size, c(shape = 1, rate = 0.01))
  expect_equal(fit$starting$magnitude$size, 10)
  expect_equal(dim(samples$magnitude$coefficients), c(1L, 3L))
  expect_equal(dim(samples$magnitude$size), c(1L, 3L))
  expect_equal(magCompStarting(fit, use = "mean")$magnitude$size, mean(samples$magnitude$size))
  expect_named(samples$composition$blocks$cls$terms, c("group", "space"))
  expect_equal(dim(samples$composition$blocks$cls$terms$group$loadings), c(2L, 1L, 3L))
  expect_equal(dim(samples$composition$blocks$cls$terms$space$covariance$phi), c(2L, 3L))
  expect_setequal(spatialDiagnostics(fit)$acceptance$factor, c("f1", "f2"))
  expect_identical(unique(spatialDiagnostics(fit)$acceptance$chain), "chain1")
  expect_equal(subset(summary(fit), parameter == "size")$symbol, "r_M")
})
