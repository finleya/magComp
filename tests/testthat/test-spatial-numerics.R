test_that("spatial prediction preserves small positive conditional variances", {
  observed <- matrix(c(0, 0), nrow = 1L)
  new <- matrix(c(1e-12, 0), nrow = 1L)
  correlation <- exp(-1e-12)
  variance <- 1 - correlation^2
  expect_gt(variance, 0)
  expect_lt(variance, 1e-10)
  set.seed(1)
  expected <- sqrt(variance) * rnorm(1L)

  for(method in c("gp", "nngp")){
    set.seed(1)
    actual <- .Call("magcomp_spatial_predict_cpp", matrix(0, 1L, 1L),
      observed, new, 1.0, 1.0, 1L, method, FALSE, "exp", 1.0, PACKAGE = "magComp")
    expect_equal(as.numeric(actual), expected, tolerance = 1e-12)
  }
})

test_that("spatial prediction reports singular matrices and zero variances", {
  # These distinct locations have correlations that round to one. Exact
  # coordinate ties are handled separately by the public prediction wrapper.
  observed <- cbind(c(0, 1e-20), c(0, 0))
  for(method in c("gp", "nngp")){
    expect_error(.Call("magcomp_spatial_predict_cpp", matrix(0, 2L, 1L),
      observed, matrix(c(1, 0), nrow = 1L), 1.0, 1.0, 2L, method, FALSE,
      "exp", 1.0, PACKAGE = "magComp"), "dpotrf failed")

    expect_error(.Call("magcomp_spatial_predict_cpp", matrix(0, 1L, 1L),
      observed[1L,,drop = FALSE], observed[2L,,drop = FALSE],
      1.0, 1.0, 1L, method, FALSE, "exp", 1.0, PACKAGE = "magComp"),
      if(method == "gp") "dpotrf failed" else "conditional variance")
  }
})

test_that("spatial fitting reports numerically singular initial correlations", {
  rows <- expand.grid(plot_id = 1:3, category = c("a", "b"))
  rows$count <- c(2, 3, 4, 4, 3, 2)
  rows$area <- 1
  rows$x <- rep(c(0, 1e-20, 2), 2)
  rows$y <- 0
  data <- magCompData(rows, unit = plot_id, category = category,
    tree_count = count, design = fixedArea(area), coords = c(x, y))
  control <- magCompControl(n_samples = 2, n_burn = 0, verbose = FALSE)

  for(method in c("gp", "nngp")){
    magnitude <- if(method == "gp") ~ gp(x, y, name = "space") else
      ~ nngp(x, y, m = 2, name = "space")
    composition <- if(method == "gp")
      ~ factor_gp(x, y, K = 1, name = "space") else
      ~ factor_nngp(x, y, K = 1, m = 2, name = "space")
    error <- if(method == "gp") "not positive definite" else "conditional variance"
    expect_error(magComp(data, magnitude = magnitude, composition = NULL,
      starting = list(magnitude = list(terms = list(space = list(covariance = list(phi = fixed(1)))))),
      control = control), error)
    expect_error(magComp(data, magnitude = NULL, composition = composition,
      starting = list(composition = list(blocks = list(category = list(terms = list(
        space = list(covariance = list(phi = fixed(1)))))))),
      control = control), error)
  }
})
