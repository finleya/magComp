# Simulate the shared population for the plot-sampling articles.
# Run prepare_data.R first, then run this script from the same directory.
# The dataset article displays these code chunks directly.

## ---- read-data ----
library(tidyverse)
library(sf)

forest <- st_read("data/forest.gpkg", quiet = TRUE)
forest_area_ha <- as.numeric(st_area(forest)) / 10000

# Read only the prepared columns so rerunning this script replaces the counts.
plots <- read_csv("data/plots.csv") |>
  select(plot_id, stand_id, easting, northing, area, p95)

prediction <- read_csv("data/prediction.csv") |>
  select(plot_id, stand_id, easting, northing, area, p95)

## ---- scale-predictor ----
# Use the plots' mean and standard deviation at every location.
p95_mean <- mean(plots$p95)
p95_sd <- sd(plots$p95)

locations <- bind_rows(
  plots |> mutate(dataset = "plot"),
  prediction |> mutate(dataset = "prediction")
) |>
  mutate(x = as.numeric(scale(p95, center = p95_mean, scale = p95_sd)))

## ---- magnitude ----
beta_0 <- log(350)
beta_1 <- -0.2
size <- 50

set.seed(1)

locations <- locations |>
  mutate(
    expected_density = exp(beta_0 + beta_1 * x),
    count = rnbinom(n(), mu = area * expected_density, size = size)
  )

## ---- species ----
# Pine is allocated first, then oak, and birch receives the remainder.
species <- c("pine", "oak", "birch")
alpha_pine <- c(intercept = 0, slope = 0.5)
alpha_oak <- c(intercept = -0.5, slope = -0.5)

pi_pine <- pnorm(alpha_pine[1] + alpha_pine[2] * locations$x)
pi_oak <- pnorm(alpha_oak[1] + alpha_oak[2] * locations$x)

species_counts <- tibble(
  pine = rbinom(nrow(locations), locations$count, pi_pine)
)

species_counts$oak <- rbinom(
  nrow(locations), locations$count - species_counts$pine, pi_oak
)

species_counts$birch <- locations$count - rowSums(species_counts)

## ---- dbh ----
dbh_classes <- c("10-20", "20-30", "30-40", "40-50", "50+")

# Each species has four allocation steps; the fifth class gets the remainder.
alpha_dbh <- tibble(
  species = rep(species, each = 4),
  dbh_class = rep(dbh_classes[1:4], times = 3),
  intercept = c(-1, -0.5, 0, 0.5,
                -0.5, 0, 0.5, 1,
                 0, 0.5, 1, 1.5),
  slope = c(-0.4, -0.3, -0.2, -0.1,
            -0.3, -0.3, -0.2, -0.2,
            -0.2, -0.2, -0.1, -0.1)
)

categories <- expand_grid(species = species, dbh_class = dbh_classes) |>
  mutate(category = paste(species, dbh_class, sep = ":"))

counts <- matrix(0L, nrow = nrow(locations), ncol = nrow(categories),
                 dimnames = list(NULL, categories$category))

for (s in species) {
  coefficients <- alpha_dbh |> filter(species == s)
  remaining <- species_counts[[s]]

  for (b in 1:4) {
    probability <- pnorm(coefficients$intercept[b] +
                          coefficients$slope[b] * locations$x)
    category <- paste(s, dbh_classes[b], sep = ":")
    counts[, category] <- rbinom(nrow(locations), remaining, probability)
    remaining <- remaining - counts[, category]
  }

  counts[, paste(s, "50+", sep = ":")] <- remaining
}

## ---- save-population ----
locations <- bind_cols(locations, as_tibble(counts))

plots <- locations |>
  filter(dataset == "plot") |>
  select(-dataset)

prediction <- locations |>
  filter(dataset == "prediction") |>
  select(-dataset)

write_csv(plots, "data/plots.csv")
write_csv(prediction, "data/prediction.csv")

## ---- forest-truth ----
# Realized forest density comes from the simulated cell counts.
forest_density <- weighted.mean(prediction$count / prediction$area,
                                 prediction$area)
forest_total <- forest_density * forest_area_ha

forest_table <- prediction |>
  select(area, all_of(categories$category)) |>
  pivot_longer(-area, names_to = "category", values_to = "count") |>
  group_by(category) |>
  summarise(density = weighted.mean(count / area, area), .groups = "drop") |>
  left_join(categories, by = "category") |>
  mutate(total = density * forest_area_ha) |>
  select(species, dbh_class, density, total)
