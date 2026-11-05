# Refresh the downloadable data after summarizing the repeated study.
# Run from dev/article-data/harvard-forest in the source repository.
library(zip)

data_files <- list.files("data", full.names = TRUE)
data_files <- data_files[basename(data_files) != "repeated"]
files <- c("inputs", data_files, "prepare_data.R", "prepare_sample.R",
           "read_inputs.R", "sample_functions.R", "population_summaries.R",
           "ladder_settings.R", "ladder_functions.R", "repeat_ladder.R",
           "summarize_ladder_repetitions.R", "export_ladder_example.R",
           "cross_validate_ladder.R", "summarize_ladder_cv.R",
           "README.md", "sources.csv")

# The archive opens into one folder and contains no saved fits or run logs.
# zipr resolves the destination relative to root, here dev/article-data.
zipr("../../vignettes/articles/downloads/harvard-forest.zip",
     paste0("harvard-forest/", files), root = "..", mode = "mirror")
