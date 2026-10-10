# Builds the package's demo datasets from the CSVs in data-raw/
# (CODE_REVIEW.md R20.52). Run from the package root:
#   Rscript data-raw/create_datasets.R
#
# Each name is one dataset. Until v0.29.0 data/ held a .rda and a different
# .csv under the same name: data() loaded the .rda (volumes on a 0-7 scale,
# negative day-0 volumes, NaN), while the dashboard demo was the CSV. The CSVs
# here are the same files as the dashboard's demos
# (mouseExperimentDashboard/inst/sample_data/).
#
# master_synthetic_data has its own generator, create_master_dataset.R.
# "synthetic_data" and "necrotic_synthetic_data" were removed in v0.29.0:
# master_synthetic_data covers what they showed, necrosis included.

read_demo <- function(file) {
  utils::read.csv(file.path("data-raw", file), check.names = FALSE,
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}

combo_treatment_synthetic_data <- read_demo("combo_treatment_synthetic_data.csv")
dose_levels_synthetic_data     <- read_demo("dose_levels_synthetic_data.csv")
weight_synthetic_data          <- read_demo("weight_synthetic_data.csv")

usethis::use_data(combo_treatment_synthetic_data, dose_levels_synthetic_data,
                  weight_synthetic_data, overwrite = TRUE)
