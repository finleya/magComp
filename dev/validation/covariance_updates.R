# Check covariance updates against exactly known uniform posterior targets.
# Zero process covariates remove covariance information from the likelihood.
# Integrating the latent effects therefore leaves the independent uniform priors.
# This checks transformed block proposals, Jacobians, adaptation, and the
# magnitude, leading-allocation, and conditional-allocation updates together.
lib <- Sys.getenv("MAGCOMP_DEV_LIB", "")
if(nzchar(lib)) .libPaths(c(lib, .libPaths()))
library(magComp)

rows <- expand.grid(plot_id = 1:4, parent = c("a", "b", "c"), child = c("small", "large"))
rows$count <- 2L
rows$area <- 1
rows$x <- c(0, 1, 0, 1)[rows$plot_id]
rows$y <- c(0, 0, 1, 1)[rows$plot_id]
rows$year <- c(0, 1, 2, 4)[rows$plot_id]
rows$z <- 0
data <- magCompData(rows, unit = plot_id, categories = c(parent, child),
  tree_count = count, design = fixedArea(area))

supports <- list(
  exp = list(phi = uniform(0.2, 6)),
  matern = list(phi = uniform(0.2, 6), nu = uniform(0.1, 2)),
  sep_exp = list(phi = uniform(0.2, 6), lambda = uniform(0.05, 2)),
  multi_res_sep_exp = list(alpha = uniform(0.1, 0.9),
    phi_1 = uniform(0.2, 2), lambda_1 = uniform(0.05, 2),
    phi_2 = uniform(2, 6), lambda_2 = uniform(0.05, 2)),
  gneiting = list(a = uniform(0.1, 2), c = uniform(0.2, 6), beta = uniform(0.1, 0.9))
)
results <- list()
for(model in names(supports)) {
  covariance <- supports[[model]]
  term_prior <- list(terms = list(space = list(covariance = covariance)))
  priors <- magCompPriors(magnitude = term_prior,
    composition = list(blocks = list(parent = term_prior, `child|parent` = term_prior)))
  for(kind in c("gp", "nngp")) {
    coordinates <- if(model %in% c("exp", "matern")) "x, y" else "x, y, year"
    term <- paste0("(", coordinates, ", name='space', cov_model='", model, "'")
    magnitude <- as.formula(paste0("~ z:", kind, term, ")"))
    leading <- as.formula(paste0("parent ~ z:factor_", kind, term, ", K=2)"))
    conditional <- as.formula(paste0("child | parent ~ z:factor_", kind, term, ", K=1)"))
    for(link in if(model == "sep_exp") c("probit", "logit") else "probit") {
      set.seed(1)
      fit <- magComp(data, magnitude = magnitude,
        composition = magCompComposition(leading, conditional), priors = priors,
        control = magCompControl(n_samples = 22000, n_burn = 2000, n_thin = 1,
          composition_link = link, verbose = FALSE))
      draws <- samples(fit)
      for(parameter in names(covariance)) {
        values <- rbind(draws$magnitude$terms$space$covariance[[parameter]],
          draws$composition$blocks$parent$terms$space$covariance[[parameter]],
          draws$composition$blocks$`child|parent`$terms$space$covariance[[parameter]])
        bounds <- covariance[[parameter]]
        for(k in seq_len(nrow(values))) {
          x <- values[k,]
          ess <- magComp:::magcomp_ess(matrix(x, ncol = 1))
          mcse <- sd(x) / sqrt(ess)
          cdf_error <- max(abs(vapply(c(0.1, 0.5, 0.9), function(p)
            mean(x <= qunif(p, bounds[1], bounds[2])) - p, numeric(1))))
          result <- data.frame(model, kind, link, parameter,
            component = c("magnitude", "parent factor 1", "parent factor 2", "conditional factor")[k],
            mean = mean(x), expected_mean = mean(bounds), mcse, ess,
            sd_ratio = sd(x) / (diff(bounds) / sqrt(12)), cdf_error)
          print(result, row.names = FALSE)
          stopifnot(all(is.finite(x)), all(x > bounds[1] & x < bounds[2]),
            abs(result$mean - result$expected_mean) < 5 * mcse,
            ess > 200, abs(result$sd_ratio - 1) < 0.12, cdf_error < 0.08)
          results[[length(results) + 1L]] <- result
        }
      }
    }
  }
}
output <- "dev/validation/outputs/covariance_updates"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, results), file.path(output, "uniform-posterior.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
