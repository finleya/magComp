# Dense-GP truth, NNGP estimation, and prediction at held-out revisits.
# Run from the package root. Three chains use one BLAS thread each.
lib <- Sys.getenv("MAGCOMP_DEV_LIB", "")
if(nzchar(lib)) .libPaths(c(lib, .libPaths()))
library(magComp)
set.seed(1)
sites <- data.frame(site = 1:80, easting = runif(80, 0, 4), northing = runif(80, 0, 4))
locations <- merge(sites, data.frame(year = c(0, 1, 2, 4, 6, 9)))
locations$plot_id <- seq_len(nrow(locations))
locations$area <- 0.12
locations$x <- rnorm(nrow(locations))
h <- as.matrix(dist(locations[c("easting", "northing")]))
u <- abs(outer(locations$year, locations$year, "-"))
v <- drop(t(chol(exp(-1.5*h - 0.3*u))) %*% rnorm(nrow(locations))) * sqrt(0.3)
w <- drop(t(chol(exp(-2*h - 0.2*u))) %*% rnorm(nrow(locations)))
M <- exp(log(350) - 0.3*locations$x + v)
pi1 <- pnorm(-0.3 + 0.3*locations$x + 0.8*w)
pi2 <- pnorm(0.1 - 0.2*locations$x - 0.6*w)
p <- cbind(pi1, (1-pi1)*pi2, (1-pi1)*(1-pi2))
locations$count <- rnbinom(nrow(locations), mu = locations$area*M, size = 1000)
y <- t(vapply(seq_len(nrow(locations)), function(i)
  as.vector(rmultinom(1, locations$count[i], p[i,])), integer(3)))
rows <- locations[rep(seq_len(nrow(locations)), 3),]
rows$species <- rep(c("a", "b", "c"), each = nrow(locations))
rows$count <- as.vector(y)
train <- locations$plot_id[locations$year < 9 | locations$site <= 40]
data <- magCompData(rows[rows$plot_id %in% train,], unit = plot_id, category = species,
  tree_count = count, design = fixedArea(area), ordered = FALSE)
term_prior <- list(covariance = list(phi = uniform(0.2, 5), lambda = uniform(0.03, 1)))
fit <- magComp(data,
  magnitude = ~ x + nngp(easting, northing, year, cov_model = "sep_exp", name = "v", m = 15),
  composition = ~ x + factor_nngp(easting, northing, year, cov_model = "sep_exp",
                                 name = "w", m = 15, anchors = "a"),
  priors = magCompPriors(magnitude = list(terms = list(v = term_prior)),
    composition = list(blocks = list(species = list(terms = list(w = term_prior))))),
  control = magCompControl(n_samples = 10000, n_burn = 4000, n_thin = 3,
    n_chains = 3, n_cores = 3, chain_seeds = 1:3, verbose = FALSE))
output <- "dev/validation/outputs/space_time_recovery"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, file.path(output, "fit.rds"))
parameters <- as.data.frame(summary(fit))
parameters$truth <- c(log(350), -0.3, 0.3, 1.5, 0.3, -0.3, 0.3, 0.1, -0.2, 0.8, -0.6, 2, 0.2)
write.csv(parameters, file.path(output, "parameters.csv"), row.names = FALSE)
combined <- summary(fit, type = "combined", group = ".all")
write.csv(combined, file.path(output, "combined-intercepts.csv"), row.names = FALSE)
print(parameters)
print(combined)
new <- locations[!locations$plot_id %in% train,]
set.seed(1)
predictions <- predict(fit, newdata = new, target = "observation", summary = "samples")
intervals <- apply(predictions$count, c(1,2), quantile, c(0.025, 0.975))
truth <- y[!locations$plot_id %in% train,]
assessment <- data.frame(plots = nrow(new), cells = length(truth),
  coverage = mean(truth >= intervals[1,,] & truth <= intervals[2,,]))
write.csv(assessment, file.path(output, "prediction.csv"), row.names = FALSE)
print(assessment)
writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
