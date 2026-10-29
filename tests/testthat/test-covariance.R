covariance_test_data <- function(nested = FALSE)
{
  set.seed(421)
  sites <- data.frame(site = 1:8, x = runif(8), y = runif(8))
  visits <- merge(sites, data.frame(year = c(0, 1, 3)))
  visits$plot_id <- seq_len(nrow(visits))
  visits$z <- rnorm(nrow(visits))
  visits$area <- 0.1
  categories <- if(nested) expand.grid(A = c("a", "b", "c"), B = c("s", "m", "l")) else
    data.frame(A = c("a", "b", "c"))
  trees <- merge(visits, categories)
  trees$count <- rpois(nrow(trees), 5)
  if(nested) magCompData(trees, unit = plot_id, categories = c(A, B),
    tree_count = count, design = fixedArea(area)) else
    magCompData(trees, unit = plot_id, category = A, tree_count = count, design = fixedArea(area))
}

test_that("registry kernels agree with independent formulas and special cases", {
  registry <- magCompCovModels()
  expect_named(registry, c("exp", "matern", "sep_exp", "multi_res_sep_exp", "gneiting"))
  xy <- cbind(c(0, 0.2, 0.7, 1), c(0.3, 0, 0.5, 0.8))
  time <- c(0, 1, 3, 7)
  st <- cbind(xy, time)
  h <- as.matrix(dist(xy))
  u <- abs(outer(time, time, "-"))
  cor <- magComp:::magcomp_correlation
  expect_equal(cor(xy, "exp", 2), unname(exp(-2*h)))
  expect_equal(cor(xy, "matern", c(2, 0.5)), cor(xy, "exp", 2))
  nu <- 1.3
  reference <- (2*h)^nu * besselK(2*h, nu) / (2^(nu-1) * gamma(nu))
  diag(reference) <- 1
  expect_equal(cor(xy, "matern", c(2, nu)), unname(reference), tolerance = 1e-12)
  expect_equal(cor(st, "sep_exp", c(2, 0.3)), unname(exp(-2*h-0.3*u)))
  expect_equal(cor(st, "multi_res_sep_exp", c(0.6, 1, 0.3, 5, 0.7)),
    unname(0.6*exp(-h-0.3*u) + 0.4*exp(-5*h-0.7*u)))
  v <- 1 + 0.2*u^1.6
  expect_equal(cor(st, "gneiting", c(0.2, 2, 0.8, 0.7, 0.6, 0.3)),
    unname(v^(-1.3) * exp(-2*h^1.2 / v^(0.7*0.6))))
  expect_error(cor(st, "gneiting", c(0.2, 2, 0, 0.7, 0.6, 0.3)), "alpha")
})

test_that("Gaussian and all-neighbor predictions have the correct joint covariance", {
  # Three future locations, including two visits to the same site. Nonunit
  # variance catches double-scaling when symmetrizing a conditional matrix.
  xy <- cbind(c(0, 0.2, 0.7, 0.3, 0.3, 0.8), c(0.3, 0, 0.5, 0.8, 0.8, 0.5))
  st <- cbind(xy, c(0, 1, 2, 3, 4, 5))
  theta <- c(2, 0.3)
  R <- exp(-2*as.matrix(dist(xy)) - 0.3*abs(outer(st[,3], st[,3], "-")))
  observed <- c(0.4, -0.3, 0.8)
  mean <- R[4:6,1:3] %*% solve(R[1:3,1:3], observed)
  V <- 0.4*(R[4:6,4:6] - R[4:6,1:3] %*% solve(R[1:3,1:3], R[1:3,4:6]))
  for(method in c("gp", "nngp")){
    set.seed(97)
    draws <- .Call("magcomp_spatial_predict_cpp", matrix(rep(observed, 6000), 3),
      st[1:3,], st[4:6,], theta, 0.4, 5L, method, FALSE, "sep_exp", 1.0, PACKAGE = "magComp")
    expect_lt(max(abs(rowMeans(draws)-mean) / sqrt(diag(V)/6000)), 5)
    covariance_se <- sqrt((V^2 + outer(diag(V), diag(V)))/5999)
    expect_lt(max(abs(cov(t(draws))-V)/covariance_se), 5)
  }
})

test_that("GP and all-neighbor NNGP latent updates reproduce their covariance", {
  # A zero process covariate leaves the composition factors distributed as their
  # Gaussian priors. This checks the sampler's precision and sparse solve,
  # independently of the prediction kernel and without MCMC convergence noise.
  points <- data.frame(plot_id = 1:6, x = c(0, 0.2, 0.7, 0.3, 0.3, 0.8),
    y = c(0.3, 0, 0.5, 0.8, 0.8, 0.5), year = c(0, 1, 2, 3, 4, 5), area = 1)
  trees <- merge(points, data.frame(A = c("a", "b")))
  trees$count <- 2L
  trees$z <- 0
  data <- magCompData(trees, unit = plot_id, category = A, tree_count = count, design = fixedArea(area))
  h <- as.matrix(dist(points[c("x", "y")]))
  u <- abs(outer(points$year, points$year, "-"))
  for(method in c("gp", "nngp")){
    formula <- as.formula(paste0("~ z:factor_", method,
      "(x,y,year,cov_model='sep_exp',name='w'", if(method == "nngp") ",m=5" else "", ")"))
    fit <- magComp(data, magnitude = NULL, composition = formula,
      starting = list(composition = list(blocks = list(A = list(terms = list(w = list(
        covariance = list(phi = fixed(2), lambda = fixed(0.3)))))))),
      control = magCompControl(n_samples = 3000, n_burn = 0, chain_seeds = 41, verbose = FALSE))
    draws <- samples(fit, type = "effects")$composition$blocks$A$terms$w$scores[,1,]
    R <- exp(-2*h-0.3*u)
    se <- sqrt((R^2 + outer(diag(R), diag(R)))/2999)
    expect_lt(max(abs(cov(t(draws))-R)/se), 5)
  }
})

test_that("space-time priors, fixed parameters, labels, and restarts agree", {
  data <- covariance_test_data()
  prior <- magCompPriors(magnitude = list(terms = list(v = list(covariance = list(
    phi = uniform(0.2, 8), lambda = uniform(0.05, 2))))),
    composition = list(blocks = list(A = list(terms = list(w = list(covariance = list(
      phi = list(f1 = uniform(0.3, 5)), lambda = uniform(0.1, 3))))))))
  start <- list(magnitude = list(terms = list(v = list(covariance = list(phi = 2, lambda = 0.4)))),
    composition = list(blocks = list(A = list(terms = list(w = list(covariance = list(
      phi = list(f1 = 1, f2 = fixed(3)))))))))
  magnitude <- ~ z + nngp(x, y, year, cov_model = "sep_exp", name = "v", m = 7)
  composition <- ~ z + factor_nngp(x, y, year, cov_model = "sep_exp", name = "w", K = 2,
                                  factor_names = c("f1", "f2"), m = 7)
  control <- magCompControl(n_samples = 80, n_burn = 30, n_chains = 2, chain_seeds = 1:2, verbose = FALSE)
  fit <- magComp(data, magnitude = magnitude, composition = composition,
    priors = prior, starting = start, control = control)
  draws <- samples(fit)
  expect_named(draws$magnitude$terms$v$covariance, c("phi", "lambda"))
  expect_identical(rownames(draws$composition$blocks$A$terms$w$covariance$phi), "f1")
  expect_true(all(fit$backend$composition_block_covariance_samples$A[2,,] == 3))
  expect_true(all(draws$magnitude$terms$v$covariance$phi > 0.2 &
                  draws$magnitude$terms$v$covariance$phi < 8))
  expect_true(all(c("phi", "lambda") %in% summary(fit)$parameter))
  expect_equal(nrow(spatialDiagnostics(fit)$acceptance), 6)
  last <- magCompStarting(fit, use = "last", chain = 2)
  expect_equal(last$magnitude$terms$v$covariance$phi,
               unname(fit$backend$magnitude_process_covariance_samples$v[1,50,2]))
  expect_s3_class(last$composition$blocks$A$terms$w$covariance$phi$f2, "magCompFixed")
  prepared <- magComp(data, magnitude = magnitude, composition = composition,
    priors = magCompSettings(fit)$priors, starting = last,
    control = magCompControl(run_sampler = FALSE))
  expect_equal(magCompStarting(prepared), last)
  # Adaptation continues after burn-in; each covariance block has one history.
  expect_length(fit$backend$magnitude_process_covariance_batch_acceptance$v[[1]], 3)
})

test_that("revisits require time and exact sampling-unit duplicates remain errors", {
  data <- covariance_test_data()
  control <- magCompControl(run_sampler = FALSE)
  expect_error(magComp(data, magnitude = ~ gp(x, y), composition = NULL, control = control),
               "identical coordinates")
  expect_s3_class(magComp(data, magnitude = ~ gp(x, y, year, cov_model = "sep_exp"),
                          composition = NULL, control = control), "magComp")
  tied <- data$data
  tied$year[tied$plot_id == 9] <- 0
  tied$x[tied$plot_id == 9] <- tied$x[tied$plot_id == 1][1]
  tied$y[tied$plot_id == 9] <- tied$y[tied$plot_id == 1][1]
  tied <- magCompData(tied, unit = plot_id, category = A, tree_count = count, design = fixedArea(area))
  expect_error(magComp(tied, magnitude = ~ nngp(x, y, year, cov_model = "sep_exp"),
                       composition = NULL, control = control), "identical coordinates")
})

test_that("all space-time kernels fit conditional factors and predict new visits", {
  data <- covariance_test_data(nested = TRUE)
  for(model in c("sep_exp", "multi_res_sep_exp", "gneiting")){
    child <- as.formula(paste0("B | A ~ z + z:factor_nngp(x,y,year,K=2,m=8,cov_model='", model, "')"))
    fit <- magComp(data, magnitude = NULL, composition = magCompComposition(A ~ z, child),
      control = magCompControl(n_samples = 35, n_burn = 10, chain_seeds = 12, verbose = FALSE))
    new <- data$unit_data[1:4,]
    new$year <- 8
    new$plot_id <- 101:104
    prediction <- predict(fit, newdata = new, summary = "samples")
    expect_true(all(is.finite(prediction$proportion)))
    expect_equal(as.numeric(apply(prediction$proportion, c(1,3), sum)), rep(1, 4*25))
    expect_true(all(is.finite(summary(fit)$mean)))
    if(model == "gneiting"){
      parameters <- subset(summary(fit), term != "")$parameter
      expect_true(all(c("a", "c", "beta") %in% parameters))
      expect_false(any(c("alpha", "gamma", "delta") %in% parameters))
    }
  }
})

test_that("covariance settings respect support, fixed domains, and search units", {
  data <- covariance_test_data()
  control <- magCompControl(run_sampler = FALSE)
  form <- ~ nngp(x, y, year, cov_model = "sep_exp", name = "v", st_scale = 0.2)
  fixed_start <- list(magnitude = list(terms = list(v = list(covariance = list(lambda = fixed(0.3))))))
  fit <- magComp(data, magnitude = form, composition = NULL, starting = fixed_start, control = control)
  expect_equal(fit$magnitude$process_terms[[1]]$st_scale, 0.2)
  expect_equal(magCompSettings(fit)$tuning$magnitude$terms$v$covariance$phi, 0.1)
  expect_null(magCompSettings(fit)$tuning$magnitude$terms$v$covariance$lambda)
  bad <- fixed_start
  bad$magnitude$terms$v$covariance$lambda <- fixed(0)
  expect_error(magComp(data, magnitude = form, composition = NULL, starting = bad, control = control),
               "invalid covariance")
  expect_error(magComp(data, magnitude = form, composition = NULL, starting = fixed_start,
    tuning = list(magnitude = list(terms = list(v = list(covariance = list(lambda = 0.2)))))), "fixed")
  # Supplying a shape prior releases that Gneiting parameter from its default.
  released <- magComp(data, magnitude = ~ gp(x, y, year, cov_model = "gneiting", name = "v"),
    composition = NULL, priors = magCompPriors(magnitude = list(terms = list(v = list(
      covariance = list(alpha = uniform(0.2, 0.9)))))), control = control)
  expect_false(inherits(magCompStarting(released)$magnitude$terms$v$covariance$alpha, "magCompFixed"))
  expect_s3_class(magCompStarting(released)$magnitude$terms$v$covariance$delta, "magCompFixed")
})

test_that("cross-validation preserves space-time fixed parameters and tuning", {
  data <- covariance_test_data()
  form <- ~ z + nngp(x, y, year, cov_model = "sep_exp", name = "v", m = 5)
  start <- list(magnitude = list(terms = list(v = list(covariance = list(lambda = fixed(0.3))))))
  tuning <- list(magnitude = list(terms = list(v = list(covariance = list(phi = 0.25)))))
  cv <- magCompCV(data, magnitude = form, composition = NULL, starting = start, tuning = tuning,
    fold_id = setNames(rep(1:2, length.out = 24), data$unit_levels), keep_fits = TRUE,
    target = "observation", summary = "samples",
    control = magCompControl(n_samples = 12, n_burn = 4, chain_seeds = 9, verbose = FALSE))
  for(fit in cv$fits){
    expect_equal(fit$tuning$magnitude$terms$v$covariance$phi, 0.25)
    expect_s3_class(magCompStarting(fit)$magnitude$terms$v$covariance$lambda, "magCompFixed")
  }
  for(prediction in cv$predictions) expect_true(all(is.finite(prediction$count)))
})
