test_that("MCMC controls count total iterations and validate retained draws", {
  control <- magCompControl()
  expect_true(control$run_sampler)
  expect_equal(c(control$n_samples, control$n_burn, control$n_thin),
               c(2000L, 1000L, 1L))
  expect_output(print(control), "retained: +1000 per chain")
  expect_output(print(control), "total iterations per chain")
  expect_error(magCompControl(n_samples = 10, n_burn = 10), "exceed n_burn")
  expect_error(magCompControl(n_samples = 9, n_burn = 10), "exceed n_burn")
  expect_error(magCompControl(n_samples = 10, n_burn = 9, n_thin = 2),
               "at least one retained")
  for(name in c("n_samples", "n_burn", "n_thin")){
    for(value in list(NA_real_, Inf, 1.5, numeric(), c(1, 2), "10")){
      args <- list(n_samples = 10, n_burn = 2, n_thin = 1)
      args[[name]] <- value
      expect_error(do.call(magCompControl, args), paste0(name, " must be"))
    }
  }
})

test_that("sampling and post-burn thinning use the requested iteration schedule", {
  rows <- expand.grid(plot = paste0("p", 1:4), cls = c("a", "b"))
  rows$count <- c(1L, 2L, 3L, 4L, 2L, 1L, 4L, 3L)
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                    design = fixedArea(1))
  control <- magCompControl(n_samples = 14, n_burn = 3, n_chains = 2,
                            chain_seeds = c(851, 852), verbose = FALSE)
  full <- magComp(dat, control = control)
  expect_true(full$fitted)
  control$n_thin <- 4L
  thinned <- magComp(dat, control = control)
  expect_equal(thinned$backend$profile_stage_seconds[["iterations"]], 28)
  expect_equal(dim(thinned$backend$magnitude_coefficient_samples), c(1L, 2L, 2L))
  # With burn-in 3 and thinning 4, retain iterations 7 and 11; still run to 14.
  expect_equal(unname(thinned$backend$magnitude_coefficient_samples),
               unname(full$backend$magnitude_coefficient_samples[,c(4L, 8L),,drop = FALSE]))
  expect_equal(unname(thinned$backend$composition_coefficient_samples),
               unname(full$backend$composition_coefficient_samples[,,c(4L, 8L),,drop = FALSE]))
  pred <- fitted(thinned, summary = "samples")
  expect_equal(pred$metadata$n_samples, 2L)
  expect_equal(pred$metadata$n_draws, 4L)

  for(burn in c(0L, 2L)){
    one <- magComp(dat, control = magCompControl(n_samples = burn + 1L,
      n_burn = burn, verbose = FALSE))
    expect_equal(one$backend$profile_stage_seconds[["iterations"]], burn + 1L)
    expect_equal(ncol(one$backend$magnitude_coefficient_samples), 1L)
  }

  control$run_sampler <- FALSE
  prepared <- magComp(dat, control = control)
  expect_false(prepared$fitted)
  expect_null(prepared$backend)
  expect_s3_class(magCompSettings(prepared), "magCompSettings")
  config <- prepared$config
  config$control$n_samples <- config$control$n_burn
  expect_error(.Call("magcomp_sampler", config, PACKAGE = "magComp"), "invalid MCMC")
  config$control$n_samples <- config$control$n_burn + 1L
  expect_error(.Call("magcomp_sampler", config, PACKAGE = "magComp"), "at least one retained")
})
