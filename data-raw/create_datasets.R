# This script regenerates all dataset files for the package

# Create directory if it doesn't exist
if (!dir.exists("data-raw")) {
  dir.create("data-raw")
}

# Load required packages
library(dplyr)

# Generate synthetic data for basic tumor growth
set.seed(123)

# Create synthetic_data
synthetic_data <- data.frame(
  Mouse_ID = rep(paste0("M", 1:12), each = 10),
  Day = rep(seq(0, 27, by = 3), 12),
  Treatment = rep(c("Control", "Treatment A", "Treatment B"), each = 40),
  Volume = runif(120, 0, 500)
)

# Add Group, ID and Cage columns
synthetic_data$Group <- synthetic_data$Treatment
synthetic_data$ID <- gsub("M", "", synthetic_data$Mouse_ID)
synthetic_data$Cage <- as.integer(as.numeric(synthetic_data$ID) %% 4) + 1

# Save the dataset
usethis::use_data(synthetic_data, overwrite = TRUE)

# Create my_data (simple example data)
my_data <- data.frame(
  x = 1:10,
  y = 1:10 + rnorm(10)
)

# Save the dataset
usethis::use_data(my_data, overwrite = TRUE)

# Create combo_treatment_synthetic_data
set.seed(456)
treatments <- c("Control", "aPD1", "HDACi", "HDACi + PD1")
n_per_group <- 8
n_timepoints <- 14
days <- seq(0, 26, by = 2)

# Function to generate growth trajectory with treatment effect
generate_trajectory <- function(treatment, mouse_id, cage) {
  base_growth <- 0.1

  effect <- switch(treatment,
                   "Control" = 0,
                   "aPD1" = -0.01,
                   "HDACi" = -0.05,
                   "HDACi + PD1" = -0.15)

  mouse_factor <- rnorm(1, mean = 1, sd = 0.2)

  volumes <- numeric(length(days))
  volumes[1] <- 0.1 + rnorm(1, mean = 0, sd = 0.05)

  for (i in 2:length(days)) {
    growth_rate <- (base_growth + effect) * mouse_factor
    noise <- rnorm(1, mean = 0, sd = 0.05 * volumes[i - 1])
    growth <- volumes[i - 1] * exp(growth_rate * (days[i] - days[i - 1]))
    volumes[i] <- max(0, growth) + noise
  }

  data.frame(
    Mouse_ID = rep(mouse_id, length(days)),
    Day = days,
    Treatment = rep(treatment, length(days)),
    Volume = volumes,
    ID = rep(substr(mouse_id, 2, 3), length(days)),
    Cage = rep(cage, length(days))
  )
}

# Generate data for all mice
combo_treatment_synthetic_data <- NULL
mouse_counter <- 1

for (treat in treatments) {
  for (i in 1:n_per_group) {
    mouse_id <- sprintf("M%02d", mouse_counter)
    cage <- ceiling(mouse_counter / 8)
    mouse_data <- generate_trajectory(treat, mouse_id, cage)
    combo_treatment_synthetic_data <- rbind(
      combo_treatment_synthetic_data, mouse_data
    )
    mouse_counter <- mouse_counter + 1
  }
}

# Save the dataset
usethis::use_data(combo_treatment_synthetic_data, overwrite = TRUE)


# Save the dataset

# Create dose_levels_synthetic_data
set.seed(789)
treatment <- "Drug X"
doses <- c(0, 10, 25, 50, 100)
n_per_dose <- 8
n_timepoints <- 14
days <- seq(0, 26, by = 2)

# Function to generate growth trajectory with dose-dependent effect
generate_dose_trajectory <- function(dose, mouse_id, cage) {
  base_growth <- 0.1

  max_effect <- -0.15
  effect <- max_effect * (dose / max(doses))

  mouse_factor <- rnorm(1, mean = 1, sd = 0.2)

  volumes <- numeric(length(days))
  volumes[1] <- 0.1 + rnorm(1, mean = 0, sd = 0.05)

  for (i in 2:length(days)) {
    growth_rate <- (base_growth + effect) * mouse_factor
    noise <- rnorm(1, mean = 0, sd = 0.05 * volumes[i - 1])
    growth <- volumes[i - 1] * exp(growth_rate * (days[i] - days[i - 1]))
    volumes[i] <- max(0, growth) + noise
  }

  data.frame(
    Mouse_ID = rep(mouse_id, length(days)),
    Day = days,
    Treatment = rep(treatment, length(days)),
    Dose = rep(dose, length(days)),
    Volume = volumes,
    ID = rep(substr(mouse_id, 2, 3), length(days)),
    Cage = rep(cage, length(days))
  )
}

# Generate data for all mice
dose_levels_synthetic_data <- NULL
mouse_counter <- 1

for (dose in doses) {
  for (i in 1:n_per_dose) {
    mouse_id <- sprintf("M%02d", mouse_counter)
    cage <- ceiling(mouse_counter / 8)
    mouse_data <- generate_dose_trajectory(dose, mouse_id, cage)
    dose_levels_synthetic_data <- rbind(
      dose_levels_synthetic_data, mouse_data
    )
    mouse_counter <- mouse_counter + 1
  }
}

# Save the dataset
usethis::use_data(dose_levels_synthetic_data, overwrite = TRUE)


# Save the dataset

cat("All datasets have been recreated and saved to the data/ directory.\n")
