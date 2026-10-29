# Generate the population used in the stratum-effects article.
# Run from dev/article-data/plot-sampling/.
# The article displays and runs these chunks directly.

## ---- read-spatial-data ----
library(tidyverse)
library(sf)

stands <- st_read("data/stands.gpkg", quiet = TRUE)

forest <- st_read("data/forest.gpkg", quiet = TRUE)

plot_centers <- st_read("data/plot_locations.gpkg", quiet = TRUE)

locations <- read_csv("data/stratum_locations.csv")

## ---- stratum-effects ----
beta_0 <- log(350)
beta_1 <- -0.2
sigma_delta <- 0.5
size <- 1000

alpha_pine <- c(intercept = 0, slope = 0.5)
alpha_oak <- c(intercept = -0.5, slope = -0.5)
gamma_pine <- 0.7
gamma_oak <- -0.5

set.seed(1)

stratum_effects <- tibble(
  stratum_id = 1:3,
  delta = rnorm(3, sd = sigma_delta),
  u = rnorm(3)
)

true_intercepts <- stratum_effects |>
  transmute(
    stratum_id,
    Magnitude = beta_0 + delta,
    Pine = alpha_pine[1] + gamma_pine * u,
    Oak = alpha_oak[1] + gamma_oak * u
  ) |>
  pivot_longer(-stratum_id, names_to = "component", values_to = "truth")

locations <- locations |>
  left_join(stratum_effects, by = "stratum_id")

## ---- magnitude-counts ----
locations <- locations |>
  mutate(
    expected_density = exp(beta_0 + beta_1 * x + delta),
    count = rnbinom(n(), mu = area * expected_density, size = size)
  )

## ---- species-counts ----
species <- c("pine", "oak", "birch")

pi_pine <- pnorm(alpha_pine[1] + alpha_pine[2] * locations$x +
                   gamma_pine * locations$u)
pi_oak <- pnorm(alpha_oak[1] + alpha_oak[2] * locations$x +
                  gamma_oak * locations$u)

species_counts <- tibble(
  pine = rbinom(nrow(locations), locations$count, pi_pine)
)

species_counts$oak <- rbinom(
  nrow(locations), locations$count - species_counts$pine, pi_oak
)

species_counts$birch <- locations$count - rowSums(species_counts)

## ---- dbh-counts ----
dbh_classes <- c("10-20", "20-30", "30-40", "40-50", "50+")

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

locations <- bind_cols(locations, as_tibble(counts))

## ---- split-population ----
plots <- locations |>
  filter(dataset == "plot") |>
  select(-dataset, -delta, -u)

prediction <- locations |>
  filter(dataset == "prediction") |>
  select(-dataset, -delta, -u)
