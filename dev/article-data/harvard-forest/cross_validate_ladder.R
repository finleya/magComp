# Run from this directory, for example: Rscript cross_validate_ladder.R H5 3
# Arguments: model (H1-H5) and cores for its three chains.
Sys.setenv(TZ = "UTC", OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library(tidyverse)
library(magComp)
source("ladder_settings.R")
source("ladder_functions.R")

args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
chain_cores <- if (length(args) >= 2) as.integer(args[2]) else 3L

trees <- read_csv("data/ladder/trees.csv", show_col_types = FALSE) |>
  mutate(species = factor(species,
    levels = c("Hemlock", "Red maple", "Red oak", "White pine", "Other")),
    dbh_class = ordered(dbh_class,
      levels = c("12.7-20", "20-30", "30-40", "40-50", "50-60", "60+")))
dat <- magCompData(trees, unit = unit_id, categories = c(species, dbh_class),
                   tree_count = count, design = fixedArea(area_ha))

# Each fold is one physical plot, including all its stratum portions.
portions <- trees |> distinct(unit_id, plot_id, stratum_id, area_ha)
fold_id <- setNames(portions$plot_id, portions$unit_id)
formulas <- ladder_formulas(model)
settings <- model_settings |> filter(.data$model == .env$model)
control <- magCompControl(
  n_samples = settings$n_samples, n_burn = settings$n_burn,
  n_thin = n_thin, n_chains = n_chains, n_cores = chain_cores,
  chain_seeds = chain_seeds, nb_size = fixed_size,
  n_report = settings$n_samples, verbose = TRUE
)

dir.create("output/cv", recursive = TRUE, showWarnings = FALSE)
cat("Starting", model, "leave-one-plot-out CV at", format(Sys.time()), "\n")
set.seed(1)
cv <- magCompCV(
  dat, magnitude = formulas$magnitude, composition = formulas$composition,
  control = control, fold_id = fold_id,
  target = "observation", type = "count", summary = "samples",
  keep_fits = TRUE
)
saveRDS(cv, file.path("output/cv", paste0(model, ".rds")))
capture.output(sessionInfo(), file = file.path("output/cv", paste0(model, "_sessionInfo.txt")))
cat("Finished", model, "at", format(Sys.time()), "\n")
