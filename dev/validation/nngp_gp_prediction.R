# Independent numerical validation, kept outside the routine test suite.
# Run from the package root; MAGCOMP_DEV_LIB can select a fresh installation.
dev_lib <- Sys.getenv("MAGCOMP_DEV_LIB", "")
if(nzchar(dev_lib)) .libPaths(c(dev_lib, .libPaths()))
library(magComp)
Sys.setenv(TZ = "UTC")
outdir <- file.path("dev", "validation", "outputs", "nngp_gp_prediction")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Independent joint Gaussian calculations. All new rows are sorted before use,
# so GP and NNGP use the same order for their independent normal innovations.
correlation <- function(a, b, phi) {
  d2 <- outer(a[, 1L], b[, 1L], "-")^2 +
    outer(a[, 2L], b[, 2L], "-")^2
  exp(-phi * sqrt(d2))
}
gp_distribution <- function(obs, new, phi, variance) {
  Coo <- correlation(obs, obs, phi)
  Cno <- correlation(new, obs, phi)
  H <- t(solve(Coo, t(Cno)))
  V <- variance * (correlation(new, new, phi) - H %*% t(Cno))
  V <- (V + t(V)) / 2
  list(H = H, V = V, L = t(chol(V)))
}
nngp_distribution <- function(obs, new, phi, variance, m) {
  n <- nrow(obs)
  q <- nrow(new)
  Bobs <- matrix(0, q, n)
  Bnew <- matrix(0, q, q)
  F <- numeric(q)
  coords <- rbind(obs, new)
  for(i in seq_len(q)) {
    eligible <- seq_len(n + i - 1L)
    d2 <- rowSums(sweep(coords[eligible, , drop = FALSE], 2L, new[i, ])^2)
    neighbors <- eligible[order(d2)[seq_len(min(m, length(eligible)))]]
    cc <- coords[neighbors, , drop = FALSE]
    r <- as.numeric(correlation(cc, new[i, , drop = FALSE], phi))
    b <- as.numeric(solve(correlation(cc, cc, phi), r))
    at_obs <- neighbors <= n
    Bobs[i, neighbors[at_obs]] <- b[at_obs]
    Bnew[i, neighbors[!at_obs] - n] <- b[!at_obs]
    F[i] <- variance * (1 - sum(b * r))
  }
  T <- solve(diag(q) - Bnew)
  stopifnot(all(is.finite(F)), all(F > 0))
  L <- sweep(T, 2L, sqrt(F), "*")
  list(H = T %*% Bobs, V = tcrossprod(L), L = L)
}
native <- function(W, obs, new, phi, variance, method, m, factor = FALSE) {
  .Call("magcomp_spatial_predict_cpp", W, obs, new, phi, variance,
        as.integer(m), method, factor, "exp", 1.0, PACKAGE = "magComp")
}
native_checks <- data.frame()
comparisons <- data.frame()
record_native <- function(case, method, actual, reference) {
  error <- max(abs(actual - reference))
  stopifnot(is.finite(error), error < 1e-8)
  native_checks <<- rbind(native_checks,
    data.frame(case, method, max_abs_error = error))
}
record_comparison <- function(case, method, g, a, W) {
  mu_gp <- g$H %*% W
  mu_approx <- a$H %*% W
  sd_gp <- sqrt(diag(g$V))
  sd_approx <- sqrt(diag(a$V))
  delta <- sweep(mu_approx - mu_gp, 1L, sd_gp, "/")
  comparisons <<- rbind(comparisons, data.frame(
    case, method,
    mean_rmse_in_gp_prediction_sd = sqrt(mean(delta^2)),
    mean_max_abs_in_gp_prediction_sd = max(abs(delta)),
    median_abs_point_sd_difference_pct = 100 * median(abs(sd_approx / sd_gp - 1)),
    max_abs_point_sd_difference_pct = 100 * max(abs(sd_approx / sd_gp - 1)),
    joint_covariance_relative_error_pct = 100 * norm(a$V - g$V, "F") / norm(g$V, "F"),
    spatial_average_sd_difference_pct = 100 * (sqrt(sum(a$V) / sum(g$V)) - 1),
    magnitude_mean_median_abs_difference_pct = if(grepl(":magnitude:", case))
      100 * median(abs(expm1(mu_approx - mu_gp + (diag(a$V) - diag(g$V)) / 2)))
      else NA_real_,
    max_abs_mean_difference = max(abs(mu_approx - mu_gp)),
    max_abs_covariance_difference = max(abs(a$V - g$V))
  ))
}

set.seed(7301)
nobs <- 60L
nnew <- 40L
nsamples <- 100L
sort_coords <- function(x) x[order(x[, 1L], x[, 2L]), , drop = FALSE]
layouts <- list(
  irregular = list(obs = matrix(runif(nobs * 2L), nobs, 2L),
                   new = sort_coords(matrix(runif(nnew * 2L), nnew, 2L))),
  gap_and_boundary = list(
    obs = cbind(c(runif(nobs / 2, 0, 0.3), runif(nobs / 2, 0.7, 1)), runif(nobs)),
    new = sort_coords(cbind(runif(nnew, -0.1, 1.1), runif(nnew, -0.1, 1.1))))
)
# phi = 1, 3, 10 covers long, intermediate, and short correlation ranges
# relative to this unit-square geometry. Variances span below/above one.
for(layout_name in names(layouts)) {
  obs <- layouts[[layout_name]]$obs
  new <- layouts[[layout_name]]$new
  for(phi in c(1, 3, 10)) for(variance in c(0.4, 1.7)) {
    label <- paste(layout_name, "magnitude", phi, variance, sep = ":")
    set.seed(7302)
    W <- sqrt(variance) * t(chol(correlation(obs, obs, phi))) %*%
      matrix(rnorm(nobs * nsamples), nobs, nsamples)
    g <- gp_distribution(obs, new, phi, variance)
    a15 <- nngp_distribution(obs, new, phi, variance, 15L)
    aall <- nngp_distribution(obs, new, phi, variance, nobs + nnew)
    set.seed(7303)
    Z <- matrix(rnorm(nnew * nsamples), nnew, nsamples)
    for(method in c("gp", "nngp15", "nngp_all")) {
      dist <- switch(method, gp = g, nngp15 = a15, nngp_all = aall)
      set.seed(7303)
      actual <- native(W, obs, new, phi, variance,
                       if(method == "gp") "gp" else "nngp",
                       if(method == "nngp15") 15L else 0L)
      record_native(label, method, actual, dist$H %*% W + dist$L %*% Z)
    }
    record_comparison(label, "nngp15", g, a15, W)
    record_comparison(label, "nngp_all", g, aall, W)
    cat("Checked", label, "\n")
    flush.console()
  }

  # The actual factor-shaped native call, with three unit-variance fields and
  # distinct decays. Check every factor, not only the scalar implementation.
  phis <- c(1, 3, 10)
  K <- length(phis)
  W <- array(NA_real_, c(nobs, K, nsamples))
  Z <- array(NA_real_, c(nnew, K, nsamples))
  set.seed(7304)
  for(k in seq_len(K))
    W[, k, ] <- t(chol(correlation(obs, obs, phis[k]))) %*%
      matrix(rnorm(nobs * nsamples), nobs, nsamples)
  set.seed(7305)
  Z[] <- rnorm(length(Z))
  for(method in c("gp", "nngp15", "nngp_all")) {
    set.seed(7305)
    actual <- native(W, obs, new, matrix(phis, K, nsamples), 1.0,
                     if(method == "gp") "gp" else "nngp",
                     if(method == "nngp15") 15L else 0L, TRUE)
    for(k in seq_len(K)) {
      g <- gp_distribution(obs, new, phis[k], 1.0)
      a <- if(method == "gp") g else
        nngp_distribution(obs, new, phis[k], 1.0,
                          if(method == "nngp15") 15L else nobs + nnew)
      label <- paste(layout_name, "factor", k, sep = ":")
      record_native(label, method, actual[, k, ], a$H %*% W[, k, ] + a$L %*% Z[, k, ])
      if(method != "gp") record_comparison(label, method, g, a, W[, k, ])
    }
  }
  cat("Checked", layout_name, "three composition factors\n")
  flush.console()
}
write.csv(native_checks, file.path(outdir, "native-reference-checks.csv"), row.names = FALSE)
write.csv(comparisons, file.path(outdir, "gp-comparison.csv"), row.names = FALSE)
cat("\nNative checks:", nrow(native_checks), "passed; maximum absolute error:",
    max(native_checks$max_abs_error), "\n")
cat("\n15-neighbor comparisons:\n")
print(comparisons[comparisons$method == "nngp15", 1:9], row.names = FALSE)
cat("\nAll-neighbor maximum mean/covariance differences:\n")
print(sapply(comparisons[comparisons$method == "nngp_all", 10:11], max))
stopifnot(all(comparisons$max_abs_mean_difference[comparisons$method == "nngp_all"] < 1e-8),
          all(comparisons$max_abs_covariance_difference[comparisons$method == "nngp_all"] < 1e-8))
capture.output(sessionInfo(), file = file.path(outdir, "session-info.txt"))
