library(magComp)

set.seed(1001)

n_unit <- 80
classes <- c("dbh_10_15", "dbh_15_20", "dbh_20_25",
             "dbh_25_30", "dbh_30_plus")
area <- rep(0.05, n_unit)
p95 <- as.numeric(scale(runif(n_unit, 12, 30)))
unit_data <- data.frame(
  plot = paste0("p", seq_len(n_unit)),
  area = area,
  p95 = p95,
  x = runif(n_unit, 0, 1000),
  y = runif(n_unit, 0, 1000)
)

magnitude_coefficients <- c(log(170), 0.35)
composition_intercept <- c(-0.05, -0.45, -0.85, -1.25)
composition_slope <- c(-0.30, -0.05, 0.20, 0.45)

counts <- matrix(0L, nrow = n_unit, ncol = length(classes),
                 dimnames = list(unit_data$plot, classes))
for(i in seq_len(n_unit)){
  total_i <- max(1L, rpois(1, area[i] *
                             exp(magnitude_coefficients[1] + magnitude_coefficients[2] * p95[i])))
  pi_i <- plogis(composition_intercept + composition_slope * p95[i])
  p_i <- numeric(length(classes))
  remaining <- 1.0
  for(j in seq_along(composition_intercept)){
    p_i[j] <- remaining * pi_i[j]
    remaining <- remaining * (1.0 - pi_i[j])
  }
  p_i[length(classes)] <- remaining
  counts[i,] <- as.vector(rmultinom(1, size = total_i, prob = p_i))
}

trees <- do.call(rbind, lapply(seq_len(n_unit), function(i){
  data.frame(
    plot = unit_data$plot[i],
    cls = rep(classes, counts[i,]),
    area = unit_data$area[i],
    p95 = unit_data$p95[i],
    x = unit_data$x[i],
    y = unit_data$y[i]
  )
}))

trees$cls <- factor(trees$cls, levels = classes, ordered = TRUE)
dat <- magCompData(
  trees,
  unit = plot,
  category = cls,
  coords = c(x, y),
  design = fixedArea(area)
)

fit <- magComp(
  dat,
  magnitude = ~ p95,
  composition = ~ p95,
  smooth = magCompSmooth(intercept = "rw2", slope = "rw2", loading = "none"),
  control = magCompControl(
    n_samples = 400,
    n_burn = 200,
    n_thin = 1,
    n_chains = 2,
    chain_seeds = c(101, 202),
    verbose = FALSE
  )
)

print(summary(fit))

out_dir <- file.path(tempdir(), "magComp_fixed_area_stems")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

trace_file <- file.path(out_dir, "traceplots.pdf")
plot(fit, file = trace_file, n_col = 2, n_row = 3)

profile_file <- file.path(out_dir, "fitted_dbh_profiles.pdf")
observed <- dat$counts
fitted_mean <- fitted(fit, type = "count")$count
fitted_ci <- fitted(fit, type = "count", summary = "quantile")$count

grDevices::pdf(profile_file, width = 8, height = 8)
old_par <- graphics::par(no.readonly = TRUE)
graphics::par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
for(i in seq_len(min(9L, nrow(observed)))){
  ylim <- range(c(observed[i,], fitted_ci[i,,]), finite = TRUE)
  graphics::plot(seq_along(classes), observed[i,], type = "h",
                 lwd = 2, xaxt = "n", xlab = "DBH class",
                 ylab = "Trees", main = rownames(observed)[i],
                 ylim = ylim)
  graphics::axis(1, at = seq_along(classes), labels = classes, las = 2,
                 cex.axis = 0.7)
  graphics::segments(seq_along(classes), fitted_ci[i,,1L],
                     seq_along(classes), fitted_ci[i,,3L],
                     col = "steelblue", lwd = 2)
  graphics::points(seq_along(classes), fitted_mean[i,],
                   pch = 19, col = "steelblue")
}
graphics::par(old_par)
grDevices::dev.off()

cat("traceplots: ", trace_file, "\n", sep = "")
cat("fitted DBH profiles: ", profile_file, "\n", sep = "")
