# Reproduce the space-time population. Run from the folder containing data/.
# Full-GP simulation uses separability to factor spatial and temporal matrices
# separately. The article reads the supplied population and fits NNGPs.

## ---- space-time-inputs ----
library(tidyverse)

locations <- read_csv("data/locations.csv")
years <- 0:9

set.seed(1)

# All 151 sites are observed initially. Revisit 25 sites in each later year.
plot_sites <- locations |> filter(dataset == "plot") |> pull(site_id)
visit_order <- sample(plot_sites)
visits <- bind_rows(
  tibble(site_id = plot_sites, year = 0),
  expand_grid(year = 1:9, visit = 1:25) |>
    mutate(site_id = visit_order[((year - 1) * 25 + visit - 1) %% 151 + 1]) |>
    select(site_id, year)
)

## ---- space-time-parameters ----
beta_0 <- log(350)
beta_1 <- -0.2
sigma_v <- 0.5
phi_v <- 3
lambda_v <- 0.35
phi_w <- 5
lambda_w <- 0.25
size <- 1000
alpha_pine <- c(intercept = 0, slope = 0.5)
alpha_oak <- c(intercept = -0.5, slope = -0.5)
loading_pine <- 0.7
loading_oak <- -0.5

## ---- space-time-fields ----
# One matrix draw contains every spatial location in every year.
# Rows index sites; columns index years. For a separable covariance,
# L_space %*% Z %*% t(L_time) is an exact full-GP draw.
distance <- as.matrix(dist(select(locations, easting_km, northing_km)))
time_lag <- abs(outer(years, years, "-"))

L_space <- t(chol(exp(-phi_v * distance)))
L_time <- t(chol(exp(-lambda_v * time_lag)))
Z <- matrix(rnorm(nrow(locations) * length(years)), nrow(locations))
v <- sigma_v * L_space %*% Z %*% t(L_time)

L_space <- t(chol(exp(-phi_w * distance)))
L_time <- t(chol(exp(-lambda_w * time_lag)))
Z <- matrix(rnorm(nrow(locations) * length(years)), nrow(locations))
w <- L_space %*% Z %*% t(L_time)
rm(distance, L_space, L_time, Z)

# Keep the sampled visits and the forest grid at the first and last years.
new_sites <- locations |> filter(dataset == "prediction") |> pull(site_id)
support <- bind_rows(visits, expand_grid(site_id = new_sites, year = c(0, 9)))
locations <- support |>
  left_join(locations, by = "site_id") |>
  mutate(plot_id = row_number(),
         v = v[cbind(site_id, year + 1)],
         w = w[cbind(site_id, year + 1)])

## ---- space-time-counts ----
locations <- locations |>
  mutate(
    magnitude_intercept = beta_0 + v,
    pine_intercept = alpha_pine[1] + loading_pine * w,
    oak_intercept = alpha_oak[1] + loading_oak * w,
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

## ---- space-time-dbh ----
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

## ---- space-time-save ----
locations <- bind_cols(locations, as_tibble(counts))
write_csv(filter(locations, dataset == "plot") |> select(-dataset), "data/plots.csv")
write_csv(filter(locations, dataset == "prediction") |> select(-dataset), "data/prediction.csv")

true_parameters <- tribble(
  ~component, ~block, ~parameter, ~index, ~truth,
  "magnitude", "", "coefficients", "(Intercept)", beta_0,
  "magnitude", "", "coefficients", "x", beta_1,
  "magnitude", "", "variance", "", sigma_v^2,
  "magnitude", "", "phi", "", phi_v,
  "magnitude", "", "lambda", "", lambda_v,
  "composition", "species", "coefficients", "pine:(Intercept)", alpha_pine[1],
  "composition", "species", "coefficients", "pine:x", alpha_pine[2],
  "composition", "species", "coefficients", "oak:(Intercept)", alpha_oak[1],
  "composition", "species", "coefficients", "oak:x", alpha_oak[2],
  "composition", "species", "loadings", "pine:pine", loading_pine,
  "composition", "species", "loadings", "oak:pine", loading_oak,
  "composition", "species", "phi", "pine", phi_w,
  "composition", "species", "lambda", "pine", lambda_w
)
write_csv(true_parameters, "data/true_parameters.csv")
