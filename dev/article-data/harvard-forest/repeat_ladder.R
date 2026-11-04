# Run from this directory: Rscript repeat_ladder.R 1 1000 4 1
# Arguments: first layout, last layout, concurrent layouts, cores per fit.
# An optional fifth argument sets the output directory for a separate run.
# Each fit has three chains. Keep workers * chain_cores within available CPUs.
Sys.setenv(TZ = "UTC", OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library(tidyverse)
library(sf)
library(magComp)
source("read_inputs.R")
source("sample_functions.R")
source("ladder_settings.R")
source("ladder_functions.R")

args <- commandArgs(trailingOnly = TRUE)
first <- as.integer(args[1])
last <- as.integer(args[2])
workers <- if (length(args) >= 3) as.integer(args[3]) else 4L
chain_cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
output_dir <- if (length(args) >= 5) args[5] else "output/repeated_ladder"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "results"), showWarnings = FALSE)
dir.create(file.path(output_dir, "fits"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

settings <- list(models = model_settings, extensions = extension_settings,
  n_chains = n_chains, chain_seeds = chain_seeds,
  n_thin = n_thin, prediction_draws = prediction_draws,
  diagnostic_draws = diagnostic_draws, rhat_limit = rhat_limit, fixed_size = fixed_size,
  package_version = as.character(packageVersion("magComp")),
  package_library = find.package("magComp"),
  package_dll_md5 = unname(tools::md5sum(list.files(
    system.file("libs", package = "magComp"), pattern = "\\.so$", full.names = TRUE))))
settings_file <- file.path(output_dir, "settings.rds")
if (file.exists(settings_file)) {
  stopifnot(identical(settings, readRDS(settings_file)))
} else {
  saveRDS(settings, settings_file)
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
}
cat("Layouts", first, "to", last, "with", workers, "workers and", chain_cores,
    "cores per fit; started", format(Sys.time()), "\n")

parallel::mclapply(first:last, function(repetition) {
  folder <- file.path(output_dir, "fits", sprintf("%04d", repetition))
  dir.create(folder, showWarnings = FALSE)
  log_file <- file.path(output_dir, "logs", sprintf("%04d.log", repetition))
  sink(log_file, append = TRUE)
  on.exit(sink())
  inventory_file <- file.path(folder, "inventory.rds")
  inventory <- if (file.exists(inventory_file)) readRDS(inventory_file) else {
    x <- make_ladder_inventory(repetition)
    saveRDS(x, inventory_file)
    x
  }
  for (model in model_settings$model) {
    result_file <- file.path(output_dir, "results", sprintf("%04d_%s.rds", repetition, model))
    if (file.exists(result_file)) next
    result <- tryCatch(run_ladder_model(repetition, model, inventory, output_dir, chain_cores),
      error = function(e) list(repetition = repetition, model = model, status = "failed",
                               message = conditionMessage(e)))
    result$sampling <- inventory$area_summary
    # An interrupted write never masquerades as a finished model.
    saveRDS(result, paste0(result_file, ".tmp"))
    file.rename(paste0(result_file, ".tmp"), result_file)
    cat("Saved", repetition, model, result$status, "at", format(Sys.time()), "\n")
    gc()
  }
  NULL
}, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE) |>
  invisible()
cat("Finished batch at", format(Sys.time()), "\n")
