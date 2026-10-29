# Dense-GP truth and three-chain NNGP recovery for the additional families.
# Run from the package root, with one BLAS thread per R process.
lib <- Sys.getenv("MAGCOMP_DEV_LIB", "")
if(nzchar(lib)) .libPaths(c(lib, .libPaths()))
library(magComp)

output <- "dev/validation/outputs/covariance_family_recovery"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

run_family <- function(model) {
  set.seed(1)
  if(model == "matern") {
    locations <- data.frame(x = runif(300, 0, 4), y = runif(300, 0, 4), year = 0)
  } else {
    sites <- data.frame(x = runif(50, 0, 4), y = runif(50, 0, 4))
    locations <- merge(sites, data.frame(year = c(0, 1, 2, 4, 6, 9)))
  }
  locations$plot_id <- seq_len(nrow(locations))
  locations$area <- 0.12
  locations$z <- rnorm(nrow(locations))
  h <- as.matrix(dist(locations[c("x", "y")]))
  u <- abs(outer(locations$year, locations$year, "-"))

  # Independent covariance formulas, without calling the package evaluator.
  if(model == "matern") {
    R <- (1.5*h)^0.8 * besselK(1.5*h, 0.8) / (2^(-0.2) * gamma(0.8))
    diag(R) <- 1
    priors <- list(phi = uniform(0.2, 4), nu = uniform(0.05, 2))
    truth <- c(phi = 1.5, nu = 0.8)
  } else if(model == "multi_res_sep_exp") {
    R <- 0.6*exp(-0.8*h - 0.15*u) + 0.4*exp(-2.8*h - 0.6*u)
    priors <- list(alpha = uniform(0, 1), phi_1 = uniform(0.2, 1.5),
      lambda_1 = uniform(0.03, 1), phi_2 = uniform(1.5, 6), lambda_2 = uniform(0.03, 1))
    truth <- c(alpha = 0.6, phi_1 = 0.8, lambda_1 = 0.15, phi_2 = 2.8, lambda_2 = 0.6)
  } else {
    T <- 1 + 0.2*u^2
    R <- T^(-1) * exp(-1.5*h / T^(0.6*0.5))
    priors <- list(a = uniform(0.03, 1), c = uniform(0.2, 4), beta = uniform(0, 1))
    truth <- c(a = 0.2, c = 1.5, beta = 0.6)
  }
  locations$v <- drop(t(chol(R)) %*% rnorm(nrow(locations))) * sqrt(0.3)
  locations$count <- rnbinom(nrow(locations),
    mu = locations$area * exp(log(350) - 0.3*locations$z + locations$v), size = 1000)
  data <- magCompData(locations, unit = plot_id, tree_count = count, design = fixedArea(area))
  coordinates <- if(model == "matern") "x, y" else "x, y, year"
  formula <- as.formula(paste0("~ z + nngp(", coordinates,
    ", cov_model='", model, "', name='v', m=15)"))
  fit <- magComp(data, magnitude = formula, composition = NULL,
    priors = magCompPriors(magnitude = list(terms = list(v = list(covariance = priors)))),
    control = magCompControl(n_samples = 8000, n_burn = 3000, n_thin = 2,
      n_chains = 3, n_cores = 1, chain_seeds = 1:3, verbose = FALSE))
  parameters <- as.data.frame(summary(fit))
  parameters$truth <- c(log(350), -0.3, 0.3, truth)
  saveRDS(fit, file.path(output, paste0(model, "-fit.rds")))
  write.csv(parameters, file.path(output, paste0(model, "-parameters.csv")), row.names = FALSE)
  combined <- summary(fit, type = "combined", group = ".all")
  write.csv(combined, file.path(output, paste0(model, "-combined.csv")), row.names = FALSE)
  print(parameters)
  print(combined)
  parameters
}

# Families run independently; chains within each family run sequentially.
results <- parallel::mclapply(c("matern", "multi_res_sep_exp", "gneiting"),
                              run_family, mc.cores = 3)
stopifnot(!any(vapply(results, inherits, logical(1), "try-error")))
saveRDS(results, file.path(output, "results.rds"))
writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
