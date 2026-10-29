# Generate the population used in the stratum-varying slopes article.
# Run from dev/article-data/plot-sampling/.
# The article displays and runs these chunks; nothing is saved to disk.

## ---- slope-data ----
library(tidyverse)

locations <- read_csv("data/stratum_locations.csv")

## ---- slope-effects ----
beta_0 <- log(350)
beta_1 <- -0.2
sigma_delta_0 <- 0.5
sigma_delta_1 <- 0.2
size <- 1000

alpha_pine <- c(intercept = 0, slope = 0.5)
alpha_oak <- c(intercept = -0.5, slope = -0.5)
gamma_pine <- c(intercept = 0.7, slope = 0.4)
gamma_oak <- c(intercept = -0.5, slope = -0.3)

set.seed(1)

stratum_effects <- tibble(
  stratum_id = 1:3,
  delta_0 = rnorm(3, sd = sigma_delta_0),
  delta_1 = rnorm(3, sd = sigma_delta_1),
  u_0 = rnorm(3),
  u_1 = rnorm(3)
)

true_coefficients <- stratum_effects |>
  transmute(
    stratum_id,
    Magnitude_intercept = beta_0 + delta_0,
    Magnitude_slope = beta_1 + delta_1,
    Pine_intercept = alpha_pine[1] + gamma_pine[1] * u_0,
    Pine_slope = alpha_pine[2] + gamma_pine[2] * u_1,
    Oak_intercept = alpha_oak[1] + gamma_oak[1] * u_0,
    Oak_slope = alpha_oak[2] + gamma_oak[2] * u_1
  ) |>
  pivot_longer(-stratum_id, names_to = c("component", "coefficient"),
               names_sep = "_", values_to = "truth") |>
  mutate(coefficient = recode(coefficient, intercept = "(Intercept)", slope = "x"))

locations <- locations |>
  left_join(stratum_effects, by = "stratum_id")

## ---- slope-totals ----
locations <- locations |>
  mutate(
    expected_density = exp(beta_0 + delta_0 + (beta_1 + delta_1) * x),
    count = rnbinom(n(), mu = area * expected_density, size = size)
  )

## ---- slope-species ----
species <- c("pine", "oak", "birch")

pi_pine <- pnorm(alpha_pine[1] + gamma_pine[1] * locations$u_0 +
                   (alpha_pine[2] + gamma_pine[2] * locations$u_1) * locations$x)

pi_oak <- pnorm(alpha_oak[1] + gamma_oak[1] * locations$u_0 +
                  (alpha_oak[2] + gamma_oak[2] * locations$u_1) * locations$x)

species_counts <- tibble(
  pine = rbinom(nrow(locations), locations$count, pi_pine)
)

species_counts$oak <- rbinom(
  nrow(locations), locations$count - species_counts$pine, pi_oak
)

species_counts$birch <- locations$count - rowSums(species_counts)

## ---- slope-dbh ----
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

## ---- slope-split ----
plots <- locations |>
  filter(dataset == "plot") |>
  select(-dataset, -delta_0, -delta_1, -u_0, -u_1)

prediction <- locations |>
  filter(dataset == "prediction") |>
  select(-dataset, -delta_0, -delta_1, -u_0, -u_1)
