# Reproduce the conditional DBH article's population with seed 1.
# Run from the folder containing data/. This replaces the supplied counts.
# Full-GP generation is optional; the article reads the supplied population.

## ---- dbh-spatial-inputs ----
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

## ---- dbh-spatial-parameters ----
beta_0 <- log(350)
beta_1 <- -0.2
size <- 50
phi_w <- 3

species <- c("pine", "oak", "birch")
dbh_classes <- c("10-20", "20-30", "30-40", "40-50", "50+")
alpha_pine <- c(intercept = 0, slope = 0.5)
alpha_oak <- c(intercept = -0.5, slope = -0.5)

# Each row describes one DBH allocation step within one species.
dbh_parameters <- tibble(
  species = rep(species, each = 4),
  dbh_class = rep(dbh_classes[1:4], times = 3),
  intercept = c(-1, -0.5, 0, 0.5,
                -0.5, 0, 0.5, 1,
                 0, 0.5, 1, 1.5),
  slope = c(-0.4, -0.3, -0.2, -0.1,
            -0.3, -0.3, -0.2, -0.2,
            -0.2, -0.2, -0.1, -0.1),
  loading = c(0.7, 0.5, 0.3, 0.1,
             -0.5, -0.4, -0.3, -0.2,
              0.4, 0.3, 0.2, 0.1)
)

## ---- dbh-spatial-field ----
# One unit-variance field shared by all species and DBH allocation steps.
coordinates <- as.matrix(select(locations, easting_km, northing_km))
distance <- as.matrix(dist(coordinates))

set.seed(1)

locations$w <- drop(t(chol(exp(-phi_w * distance))) %*%
                      rnorm(nrow(locations)))

rm(distance)

## ---- dbh-spatial-totals ----
# Magnitude and species allocation depend on canopy height alone.
locations <- locations |>
  mutate(
    expected_density = exp(beta_0 + beta_1 * x),
    count = rnbinom(n(), mu = area * expected_density, size = size)
  )

pi_pine <- pnorm(alpha_pine[1] + alpha_pine[2] * locations$x)
pi_oak <- pnorm(alpha_oak[1] + alpha_oak[2] * locations$x)

species_counts <- tibble(
  pine = rbinom(nrow(locations), locations$count, pi_pine)
)

species_counts$oak <- rbinom(
  nrow(locations), locations$count - species_counts$pine, pi_oak
)

species_counts$birch <- locations$count - rowSums(species_counts)

## ---- dbh-spatial-allocation ----
categories <- expand_grid(species = species, dbh_class = dbh_classes) |>
  mutate(category = paste(species, dbh_class, sep = ":"))

counts <- matrix(0L, nrow = nrow(locations), ncol = nrow(categories),
                 dimnames = list(NULL, categories$category))

for (s in species) {
  coefficients <- dbh_parameters |> filter(species == s)
  remaining <- species_counts[[s]]

  for (b in 1:4) {
    probability <- pnorm(coefficients$intercept[b] +
                          coefficients$slope[b] * locations$x +
                          coefficients$loading[b] * locations$w)
    category <- paste(s, dbh_classes[b], sep = ":")
    counts[, category] <- rbinom(nrow(locations), remaining, probability)
    remaining <- remaining - counts[, category]
  }

  counts[, paste(s, "50+", sep = ":")] <- remaining
}

## ---- dbh-spatial-save ----
locations <- bind_cols(locations, as_tibble(counts))

plots <- locations |> filter(dataset == "plot") |> select(-dataset)
prediction <- locations |> filter(dataset == "prediction") |> select(-dataset)

write_csv(plots, "data/plots.csv")
write_csv(prediction, "data/prediction.csv")
write_csv(dbh_parameters, "data/dbh_parameters.csv")

true_parameters <- tribble(
  ~component, ~block, ~parameter, ~index, ~truth,
  "magnitude", "", "coefficients", "(Intercept)", beta_0,
  "magnitude", "", "coefficients", "x", beta_1,
  "magnitude", "", "size", "", size,
  "composition", "species", "coefficients", "pine:(Intercept)", alpha_pine[1],
  "composition", "species", "coefficients", "pine:x", alpha_pine[2],
  "composition", "species", "coefficients", "oak:(Intercept)", alpha_oak[1],
  "composition", "species", "coefficients", "oak:x", alpha_oak[2],
  "composition", "dbh_class|species", "phi", "size", phi_w
)

write_csv(true_parameters, "data/true_parameters.csv")
