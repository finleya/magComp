library(magComp)

# A small count table: two species and three DBH classes at each plot.
# The short runs here illustrate the interface, not convergence assessment.
set.seed(42)
rows <- expand.grid(plot = paste0("p", 1:20),
                    species = c("pine", "oak"), dbh = c("small", "medium", "large"))
rows$dbh <- ordered(rows$dbh, levels = c("small", "medium", "large"))
rows$x <- rep(seq(-1, 1, length.out = 20), 6)
rows$count <- rnbinom(nrow(rows), size = 5, mu = exp(1 + 0.4 * rows$x))
rows$area <- 0.05
control <- magCompControl(n_samples = 400, n_burn = 200, verbose = FALSE)

# Magnitude alone: category labels are unnecessary. Counts within a plot add up.
magnitude_data <- magCompData(rows, unit = plot, tree_count = count,
                              design = fixedArea(area))
magnitude_fit <- magComp(magnitude_data, magnitude = ~ x, composition = NULL,
                         control = control)
fitted(magnitude_fit)                 # Expected density and total count.
fitted(magnitude_fit, target = "observation")  # Negative-binomial replicates.

# Composition alone: no plot area is needed. Nest DBH allocation within species.
composition_data <- magCompData(rows, unit = plot, tree_count = count,
                                categories = c(species, dbh))
composition_fit <- magComp(composition_data, magnitude = NULL,
  composition = magCompComposition(species ~ x, dbh | species ~ x),
  control = control)
fitted(composition_fit)               # Species-by-DBH probabilities.
fitted(composition_fit, type = "count") # Expected counts given observed totals.
ppcSummary(composition_fit)           # Multinomial checks given observed totals.

# New locations: magnitude supplies totals; composition conditions on supplied totals.
newdata <- data.frame(plot = c("new1", "new2"), x = c(-0.5, 0.5),
                      area = 0.05, total = c(20L, 30L))
predict(magnitude_fit, newdata)
predict(composition_fit, newdata)      # Probabilities only by default.
conditional <- predict(composition_fit, newdata,
  type = c("proportion", "count", "density"), total_count = "total", area = "area",
  target = "observation", summary = "samples")
conditional$count                    # Unit x final category x retained draw.
conditional$metadata$count_conditioning

# Stand aggregation requires the full magnitude-composition model.
# Standalone fits use the same priors, starting values, and diagnostic functions.
magCompSettings(composition_fit)$parameters
summary(magnitude_fit)
