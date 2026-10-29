# Independent posterior benchmark for a one-factor binomial mixed model.
# Integrate the iid scores by Gaussian quadrature, then integrate the
# intercept and positive loading on a grid. This does not use sampler kernels.
dev_lib <- Sys.getenv("MAGCOMP_DEV_LIB", "")
if(nzchar(dev_lib)) .libPaths(c(dev_lib, .libPaths()))
library(magComp)
out <- Sys.getenv("MAGCOMP_DIAG_OUTPUT",
                  "dev/validation/outputs/positive_diagonal_reference")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
gh <- function(n) {
  J <- matrix(0, n, n)
  J[cbind(1:(n-1),2:n)] <- sqrt(1:(n-1))
  e <- eigen(J + t(J), symmetric = TRUE)
  list(x = e$values, w = e$vectors[1,]^2)
}
grid_posterior <- function(y, link, nodes, step) {
  a <- seq(-2.5, 1.5, by = step)
  d <- seq(0, 3, by = step)
  grid <- expand.grid(alpha = a, diagonal = d)
  quad <- gh(nodes)
  eta <- outer(grid$diagonal, quad$x) + grid$alpha
  p <- if(link == 'probit') pnorm(eta) else plogis(eta)
  log_density <- dnorm(grid$alpha, 0, 3, log = TRUE) +
    dnorm(grid$diagonal, 0, 1, log = TRUE)
  frequencies <- table(y)
  for(s in names(frequencies)) {
    likelihood <- matrix(dbinom(as.integer(s), 20, p), nrow(grid)) %*% quad$w
    log_density <- log_density + frequencies[[s]] * log(pmax(likelihood, 1e-300))
  }
  density <- matrix(exp(log_density - max(log_density)), length(a), length(d))
  # Trapezoidal weights, including the half-normal's boundary at zero.
  density[c(1,length(a)),] <- density[c(1,length(a)),] / 2
  density[,c(1,length(d))] <- density[,c(1,length(d))] / 2
  density <- density / sum(density)
  marginal <- colSums(density)
  list(mean = sum(d * marginal), sd = sqrt(sum(d^2 * marginal) - sum(d*marginal)^2),
       alpha = sum(a * rowSums(density)), d = d, cdf = cumsum(marginal),
       edge_mass = sum(density[c(1,length(a)),]) + sum(density[,length(d)]))
}
results <- list()
for(link in c('probit', 'logit')) {
  set.seed(883)
  n <- 36L
  eta <- -0.5 + 0.75 * rnorm(n)
  y <- rbinom(n, 20, if(link == 'probit') pnorm(eta) else plogis(eta))
  rows <- expand.grid(plot = sprintf('p%02d', seq_len(n)), cls = c('a','b'))
  rows$count <- c(y, 20-y)
  rows$area <- 1
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     design = fixedArea(area))
  fit <- magComp(dat, composition = ~factor_iid(plot, K = 1, name = 'shape', anchors = 'a'),
    control = magCompControl(n_chains = 2, n_cores = 2, n_burn = 10000,
      n_samples = 70000, chain_seeds = c(444,555), composition_link = link,
      nb_size = 40, verbose = FALSE))
  samples <- samples(fit)$composition$blocks$cls
  diagonal <- samples$terms$shape$loadings[1,1,,]
  reference <- grid_posterior(y, link, 80, 0.02)
  refined <- grid_posterior(y, link, 160, 0.01)
  mcse <- sd(c(diagonal)) / sqrt(magComp:::magcomp_ess(diagonal))
  grid_cdf <- approx(refined$d + diff(refined$d)[1] / 2, refined$cdf,
                     xout = sort(c(diagonal)), rule = 2)$y
  empirical_cdf <- (seq_along(grid_cdf) - 0.5) / length(grid_cdf)
  result <- data.frame(link = link, mcmc = mean(diagonal), quadrature = refined$mean,
    numerical_difference = abs(reference$mean-refined$mean),
    mcse = mcse, standardized_error = (mean(diagonal)-refined$mean)/mcse,
    rhat = magComp:::magcomp_rhat(diagonal), ess = magComp:::magcomp_ess(diagonal),
    boundary_mass = refined$edge_mass,
    sd_ratio = sd(c(diagonal)) / refined$sd,
    cdf_difference = max(abs(grid_cdf - empirical_cdf)))
  print(result)
  stopifnot(result$numerical_difference < 1e-4, result$boundary_mass < 1e-6,
            abs(result$standardized_error) < 4, result$rhat < 1.01, result$ess > 1000,
            abs(result$sd_ratio - 1) < 0.04, result$cdf_difference < 0.03)
  saveRDS(list(result=result, reference=refined, draws=diagonal),
          file.path(out,paste0('quadrature-',link,'.rds')))
  results[[link]] <- result
}
write.csv(do.call(rbind, results), file.path(out,'quadrature.csv'), row.names=FALSE)
