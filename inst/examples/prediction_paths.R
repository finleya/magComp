library(magComp)

# Serial fitting and subsequent simulations use the same external RNG stream.
set.seed(301)

# Small fixed-area example with two response categories.
rows <- expand.grid(plot = paste0("p", 1:8), cls = c("small", "large"))
rows$count <- c(2L, 4L, 1L, 0L, 5L, 3L, 2L, 1L,
                1L, 2L, 0L, 1L, 2L, 3L, 1L, 2L)
rows$area <- 0.05
rows$height <- rep(seq(10, 24, by = 2), 2)
dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                   design = fixedArea(area))
fit <- magComp(dat, magnitude = ~ height, composition = ~ height,
  control = magCompControl(n_samples = 80, n_burn = 40, n_chains = 2,
    nb_size = 10,
    verbose = FALSE))

# 1. Expected quantities at fitted units, using their original areas.
observed_expected <- fitted(fit, summary = "samples")

# 2. Replicate observations at those same units, from the fitted likelihood.
observed_replicates <- fitted(fit, target = "observation",
                              summary = "samples")

grid <- data.frame(plot = c("new1", "new2"), height = c(15, 20),
                    plot_area = c(0.1, 0.2))

# 3. Expected quantities at new units with explicit prediction exposures.
new_expected <- predict(fit, grid, area = "plot_area",
                         summary = "samples")

# 4. Simulated observations at new units, with the same model and exposures.
new_observations <- predict(fit, grid, area = "plot_area",
                             target = "observation", summary = "samples")
print(new_observations)
names(new_observations)
dim(new_observations$count)     # unit, category, retained sample, chain
dim(new_observations$magnitude) # unit, retained sample, chain

# Summarize these stored draws without simulating new observations.
summary(new_observations)

# A requested subset is still a list. Means and quantiles pool chains.
count_mean <- predict(fit, grid, area = "plot_area", type = "count")$count
count_median <- predict(fit, grid, area = "plot_area", type = "count",
                        summary = "quantile", probs = 0.5)$count

# Zero simulated totals have undefined proportions (NA). Proportion summaries
# condition on positive totals; metadata makes that conditioning explicit.
proportion_summary <- fitted(fit, type = "proportion", target = "observation")
proportion_summary$proportion
proportion_summary$metadata$zero_total_probability
proportion_summary$metadata$n_positive_draws

# Posterior predictive model-fit checks use the observed observation path.
ppcSummary(fit, type = "unit_total")

# Aggregate these same draws. Each grid point here represents 2 hectares,
# independently of the plot_area used to simulate its count.
areas <- data.frame(prediction_id = grid$plot, stand_id = "stand1", area_ha = 2)
stands <- standSummary(new_observations, data = areas, unit = prediction_id,
                       group = stand_id, represented_area = area_ha,
                       summary = "samples")
summary(stands)
stands$density
