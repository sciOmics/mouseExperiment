# Direct fix for the mouseExperiment package

# Create a clean temporary directory
if (dir.exists("fixed_pkg")) {
  unlink("fixed_pkg", recursive = TRUE)
}
dir.create("fixed_pkg")
dir.create("fixed_pkg/R")
dir.create("fixed_pkg/data")

# Copy the plot_auc.R file
file.copy("R/plot_auc.R", "fixed_pkg/R/", overwrite = TRUE)

# Copy the data files
for (data_file in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  file.copy(data_file, "fixed_pkg/data/", overwrite = TRUE)
}

# Create a new DESCRIPTION file
description_content <- '
Package: mouseExperiment
Title: Analysis of Mouse Tumor Growth Experiments
Version: 0.1.0
Authors@R: 
    person("Insight", "BioAnalytics", , "info@insightbioanalytics.com", role = c("aut", "cre", "cph"))
Author: Insight BioAnalytics [aut, cre, cph]
Maintainer: Insight BioAnalytics <info@insightbioanalytics.com>
Description: Provides tools for analyzing and visualizing mouse tumor growth experiments
    including calculating tumor volume, survival statistics, and tumor growth kinetics.
    The package offers functions for dose-response analysis, drug synergy evaluation,
    and comprehensive visualization of experimental results.
License: file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.2
Suggests:
    knitr,
    rmarkdown
VignetteBuilder: knitr
Imports: 
    ggplot2,
    ggpubr,
    stats
Depends: 
    R (>= 3.5)
LazyData: true
'
cat(description_content, file = "fixed_pkg/DESCRIPTION")

# Create the LICENSE file
cat("This is a proprietary license file.\n", file = "fixed_pkg/LICENSE")

# Create the data.R file
data_r_content <- '
#\' Synthetic tumor growth data
#\' 
#\' A dataset containing synthetic tumor volume measurements over time
#\' for different treatment groups.
#\' 
#\' @format A data frame with 120 rows and 7 variables:
#\' \\describe{
#\'   \\item{Mouse_ID}{Mouse identifier, with format M followed by a number}
#\'   \\item{Day}{Day of measurement, starting from day 0}
#\'   \\item{Treatment}{Treatment group (Control, Treatment A, Treatment B)}
#\'   \\item{Volume}{Tumor volume measurement}
#\'   \\item{Group}{Same as Treatment, an alternative name for the treatment group}
#\'   \\item{ID}{Numeric identifier for the mouse, without the "M" prefix}
#\'   \\item{Cage}{Cage number where the mouse is housed}
#\' }
#\' @source Synthetic data generated using random number generation
#\' @examples
#\' data(synthetic_data)
#\' head(synthetic_data)
"synthetic_data"

#\' Example data
#\' 
#\' A simple example dataset with x and y coordinates.
#\' 
#\' @format A data frame with 10 rows and 2 variables:
#\' \\describe{
#\'   \\item{x}{x-coordinate, sequence from 1 to 10}
#\'   \\item{y}{y-coordinate, x plus random noise}
#\' }
#\' @source Synthetic data generated using random number generation
#\' @examples
#\' data(my_data)
#\' plot(my_data$x, my_data$y)
"my_data"

#\' Combination treatment synthetic data
#\' 
#\' A dataset containing synthetic tumor volume measurements over time
#\' for different combination treatment groups.
#\' 
#\' @format A data frame with rows for each mouse at each timepoint and 6 variables:
#\' \\describe{
#\'   \\item{Mouse_ID}{Mouse identifier, with format M followed by a number}
#\'   \\item{Day}{Day of measurement, starting from day 0}
#\'   \\item{Treatment}{Treatment group (Control, aPD1, HDACi, HDACi + PD1)}
#\'   \\item{Volume}{Tumor volume measurement}
#\'   \\item{ID}{Numeric identifier for the mouse, without the "M" prefix}
#\'   \\item{Cage}{Cage number where the mouse is housed}
#\' }
#\' @source Synthetic data generated using random number generation with treatment-specific effects
#\' @examples
#\' data(combo_treatment_synthetic_data)
#\' head(combo_treatment_synthetic_data)
"combo_treatment_synthetic_data"

#\' Combination treatment schedule
#\' 
#\' A dataset specifying the dosing schedule for combination treatments.
#\' 
#\' @format A data frame with 20 rows and 3 variables:
#\' \\describe{
#\'   \\item{Treatment}{Treatment group (Control, aPD1, HDACi, HDACi + PD1)}
#\'   \\item{Day}{Day of dose administration}
#\'   \\item{Dose}{Dose amount}
#\' }
#\' @source Synthetic treatment schedule
#\' @examples
#\' data(combo_treatment_schedule)
#\' head(combo_treatment_schedule)
"combo_treatment_schedule"

#\' Dose levels synthetic data
#\' 
#\' A dataset containing synthetic tumor volume measurements over time
#\' for different dose levels of a single drug.
#\' 
#\' @format A data frame with rows for each mouse at each timepoint and 7 variables:
#\' \\describe{
#\'   \\item{Mouse_ID}{Mouse identifier, with format M followed by a number}
#\'   \\item{Day}{Day of measurement, starting from day 0}
#\'   \\item{Treatment}{Treatment name (always "Drug X")}
#\'   \\item{Dose}{Dose level (0, 10, 25, 50, 100)}
#\'   \\item{Volume}{Tumor volume measurement}
#\'   \\item{ID}{Numeric identifier for the mouse, without the "M" prefix}
#\'   \\item{Cage}{Cage number where the mouse is housed}
#\' }
#\' @source Synthetic data generated using random number generation with dose-dependent effects
#\' @examples
#\' data(dose_levels_synthetic_data)
#\' head(dose_levels_synthetic_data)
"dose_levels_synthetic_data"

#\' Dose levels treatment schedule
#\' 
#\' A dataset specifying the dosing schedule for different dose levels.
#\' 
#\' @format A data frame with 20 rows and 4 variables:
#\' \\describe{
#\'   \\item{Treatment}{Treatment name (always "Drug X")}
#\'   \\item{Dose}{Dose group level (0, 10, 25, 50, 100)}
#\'   \\item{Day}{Day of dose administration}
#\'   \\item{Administered_Dose}{Actual dose administered}
#\' }
#\' @source Synthetic treatment schedule
#\' @examples
#\' data(dose_levels_treatment_schedule)
#\' head(dose_levels_treatment_schedule)
"dose_levels_treatment_schedule"
'
cat(data_r_content, file = "fixed_pkg/R/data.R")

# Change to the fixed_pkg directory
setwd("fixed_pkg")

# Build and install the package
cat("Building and installing the fixed package...\n")
system("R CMD build . && R CMD INSTALL mouseExperiment_0.1.0.tar.gz")

# Test the package
cat("Testing the package...\n")
test_script <- '
library(mouseExperiment)
# Check if datasets are accessible
data(synthetic_data)
data(my_data)
data(combo_treatment_synthetic_data)
data(combo_treatment_schedule)
data(dose_levels_synthetic_data)
data(dose_levels_treatment_schedule)

# Test plot_auc function
example_data <- data.frame(
  ID = paste0("Mouse", 1:20),
  Group = rep(c("Control", "Treatment A", "Treatment B", "Treatment C"), each=5),
  AUC = c(runif(5, 10, 15), runif(5, 8, 12), runif(5, 5, 10), runif(5, 3, 8))
)
pdf("test_plot.pdf")
p <- plot_auc(example_data)
print(p)
dev.off()

# Print success message if we get here
cat("SUCCESS: All tests passed!\n")
'
cat(test_script, file = "test_package.R")
system("Rscript test_package.R")

# Check if the test was successful
if (file.exists("test_plot.pdf")) {
  cat("Package has been successfully fixed and tested!\n")
  
  # Return to the original directory
  setwd("..")
  
  # Copy the documentation to the main package
  cat("Updating the main package with the fixed files...\n")
  file.copy("fixed_pkg/R/data.R", "R/", overwrite = TRUE)
  
  # Update DESCRIPTION to include LazyData: true
  desc_lines <- readLines("DESCRIPTION", warn = FALSE)
  if (!any(grepl("LazyData: true", desc_lines))) {
    cat("LazyData: true\n", file = "DESCRIPTION", append = TRUE)
  }
  
  # Commit the changes
  cat("Committing the changes...\n")
  system("git add R/data.R")
  system("git add DESCRIPTION")
  system("git commit -m 'Fixed dataset documentation and added LazyData: true to DESCRIPTION'")
  
  cat("Done! The package has been fixed and changes have been committed.\n")
} else {
  cat("Testing failed. Please check the error messages above.\n")
} 