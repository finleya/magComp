# Reproduce the spatial article's full-GP population with seed 1.
# Run from the folder containing data/. This replaces the supplied population.
# Generation uses dense covariance matrices; the article instead reads the
# supplied CSV files and fits a nearest-neighbor Gaussian process (NNGP).

## ---- spatial-inputs ----
library(tidyverse)

plots <- read_csv("data/plots.csv") |>
  select(plot_id, stand_id, easting, northing, area, p95, x)

prediction <- read_csv("data/prediction.csv") |>
  select(plot_id, stand_id, easting, northing, area, p95, x)

locations <- bind_rows(
  plots |> mutate(dataset = "plot"),
  prediction |> mutate(dataset = "prediction")
) |>
  mutate(easting_km = easting / 1000, northing_km = northing / 1000)

## ---- spatial-parameters ----
beta_0 <- log(350)
beta_1 <- -0.2
sigma_v <- 0.5
phi_v <- 3
phi_w <- 5
size <- 1000

alpha_pine <- c(intercept = 0, slope = 0.5)
alpha_oak <- c(intercept = -0.5, slope = -0.5)
lambda_pine <- 0.7
lambda_oak <- -0.5

## ---- spatial-fields ----
# Both fields include every observed and prediction location in one draw.
coordinates <- as.matrix(select(locations, easting_km, northing_km))
distance <- as.matrix(dist(coordinates))

set.seed(1)

locations$v <- drop(t(chol(sigma_v^2 * exp(-phi_v * distance))) %*%
                      rnorm(nrow(locations)))
locations$w <- drop(t(chol(exp(-phi_w * distance))) %*%
                      rnorm(nrow(locations)))

rm(distance)

## ---- spatial-counts ----
locations <- locations |>
  mutate(
    magnitude_intercept = beta_0 + v,
    pine_intercept = alpha_pine[1] + lambda_pine * w,
    oak_intercept = alpha_oak[1] + lambda_oak * w,
    expected_density = exp(magnitude_intercept + beta_1 * x),
    count = rnbinom(n(), mu = area * expected_density, size = size)
  )

species <- c("pine", "oak", "birch")
pi_pine <- pnorm(locations$pine_intercept + alpha_pine[2] * locations$x)
pi_oak <- pnorm(locations$oak_intercept + alpha_oak[2] * locations$x)

species_counts <- tibble(
  pine = rbinom(nrow(locations), locations$count, pi_pine)
)

species_counts$oak <- rbinom(
  nrow(locations), locations$count - species_counts$pine, pi_oak
)

species_counts$birch <- locations$count - rowSums(species_counts)

## ---- spatial-dbh ----
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

## ---- spatial-save ----
locations <- bind_cols(locations, as_tibble(counts))

plots <- locations |> filter(dataset == "plot") |> select(-dataset)
prediction <- locations |> filter(dataset == "prediction") |> select(-dataset)

write_csv(plots, "data/plots.csv")
write_csv(prediction, "data/prediction.csv")

true_parameters <- tribble(
  ~component, ~block, ~parameter, ~index, ~truth,
  "magnitude", "", "coefficients", "(Intercept)", beta_0,
  "magnitude", "", "coefficients", "x", beta_1,
  "magnitude", "", "variance", "", sigma_v^2,
  "magnitude", "", "phi", "", phi_v,
  "composition", "species", "coefficients", "pine:(Intercept)", alpha_pine[1],
  "composition", "species", "coefficients", "pine:x", alpha_pine[2],
  "composition", "species", "coefficients", "oak:(Intercept)", alpha_oak[1],
  "composition", "species", "coefficients", "oak:x", alpha_oak[2],
  "composition", "species", "loadings", "pine:pine", lambda_pine,
  "composition", "species", "loadings", "oak:pine", lambda_oak,
  "composition", "species", "phi", "pine", phi_w
)

write_csv(true_parameters, "data/true_parameters.csv")
