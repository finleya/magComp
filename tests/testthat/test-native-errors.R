test_that("native failure cleans up sparse workspaces and records consumed RNG state", {
  rows <- expand.grid(plot = paste0("p", 1:6), cls = c("s", "l"))
  rows$x <- rep(1:6, 2)
  rows$y <- rows$x^2
  dat <- magCompData(rows, unit = plot, category = cls, design = fixedArea(1))
  model <- magComp(dat, magnitude = ~ nngp(x, y, m = 2),
    control = magCompControl(n_samples = 3, n_burn = 1, verbose = FALSE, run_sampler = FALSE))
  config <- model$config
  bad <- config
  bad$priors$coefficient_variance <- 0 # Fail after augmentation and sparse setup.
  set.seed(713)
  before <- .Random.seed
  expect_error(.Call("magcomp_sampler", bad, PACKAGE = "magComp"),
               "non-finite precision")
  expect_false(identical(before, .Random.seed))
  expect_no_error(.Call("magcomp_sampler", config, PACKAGE = "magComp"))

  expect_error(magCompControl(n_samples = .Machine$integer.max, n_thin = 2, run_sampler = FALSE),
               "integer indexing limit")
  config$control$n_samples <- .Machine$integer.max - 2L
  config$control$n_burn <- 0L
  expect_error(.Call("magcomp_sampler", config, PACKAGE = "magComp"),
               "array exceeds|matrix exceeds")
})

test_that("a terminated chain is reported before combining results", {
  check <- getFromNamespace("magcomp_check_backend_chain_errors", "magComp")
  expect_error(check(list(list(), NULL)), "chain 2 returned no result")
})

test_that("failed iid predictions release RNG state and allow subsequent calls", {
  for (factor in c(FALSE, TRUE)) {
    kernel <- if (factor) "magcomp_predict_iid_factor_cpp" else "magcomp_predict_iid_scalar_cpp"
    effect <- if (factor) array(0, c(1, 2, 1)) else matrix(0, 1, 1)
    args <- list(effect, 0L, 1L, 3L)
    if (!factor) args <- c(args, list(1.0))
    invoke <- function(args) do.call(.Call, c(list(kernel), args, list(PACKAGE = "magComp")))

    # Both invalid indices fail after the same normal draws as the valid call.
    set.seed(831)
    expected <- invoke(args)
    expected_seed <- .Random.seed
    for (observed in c(0L, 2L)) {
      bad <- args
      bad[[2]] <- observed
      bad[[3]] <- 0L
      set.seed(831)
      expect_error(invoke(bad), "invalid (new|observed) iid .* group index")
      expect_identical(.Random.seed, expected_seed)
    }
    set.seed(831)
    expect_identical(invoke(args), expected)
    expect_identical(.Random.seed, expected_seed)
  }
})

test_that("spatial prediction failure records earlier draws for GP and NNGP", {
  coords <- matrix(c(0, 1, 0, 1), 2, 2)
  new_coords <- matrix(c(0.5, 0.5), 1, 2)
  for (method in c("gp", "nngp")) {
    for (factor in c(FALSE, TRUE)) {
      one <- if (factor) array(0, c(2, 2, 1)) else matrix(0, 2, 1)
      two <- if (factor) array(0, c(2, 2, 2)) else matrix(0, 2, 2)
      invoke <- function(effect, variance) .Call("magcomp_spatial_predict_cpp",
        effect, coords, new_coords, 1.0, variance, 2L, method, factor,
        "exp", 1.0, PACKAGE = "magComp")
      set.seed(924)
      expected <- invoke(one, 1.0)
      expected_seed <- .Random.seed
      set.seed(924)
      expect_error(invoke(two, matrix(c(1, -1), 1, 2)),
                   "spatial prediction sigma_sq must be positive")
      expect_identical(.Random.seed, expected_seed)
      set.seed(924)
      expect_identical(invoke(one, 1.0), expected)
    }
  }
})

test_that("R API errors inside prediction return safely to R", {
  # STRING_ELT(integer) raises an R condition while native dimension vectors
  # are live; this exercises the R-to-C++ unwind bridge, not magcomp::error.
  coords <- matrix(c(0, 1), 2, 1)
  set.seed(925)
  before <- .Random.seed
  expect_error(.Call("magcomp_spatial_predict_cpp",
    matrix(0, 2, 1), coords, matrix(0.5, 1, 1), 1.0, 1.0, 1L, 1L, FALSE,
    "exp", 1.0, PACKAGE = "magComp"), "STRING_ELT")
  expect_identical(.Random.seed, before)
  expect_equal(dim(.Call("magcomp_predict_iid_scalar_cpp",
    matrix(0, 1, 1), 0L, 1L, 1L, 1.0, PACKAGE = "magComp")), c(1L, 1L))
})
