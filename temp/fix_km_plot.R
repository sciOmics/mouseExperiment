# Test script to verify and fix the km_plot issue in survival_statistics
# The issue is that the KM plot shows too many events in each group

library(mouseExperiment)

# Load sample data 
data(combo_treatment_synthetic_data)

# Check the structure of the data
cat("Structure of dataset:\n")
str(combo_treatment_synthetic_data[, c("ID", "Day", "Treatment", "Cage")])

# The dataset doesn't have a Survival_Censor column, so we'll create synthetic survival data
# Set a seed for reproducibility
set.seed(123)

# Get unique IDs
unique_ids <- unique(combo_treatment_synthetic_data$ID)
cat("\nNumber of unique subjects: ", length(unique_ids), "\n")

# Create synthetic survival data with approximately 50% of subjects having events
survival_data <- data.frame(
  ID = unique_ids,
  Survival_Time = round(runif(length(unique_ids), 10, 50)),
  Survival_Censor = sample(c(0, 1), length(unique_ids), replace = TRUE, prob = c(0.5, 0.5)),
  stringsAsFactors = FALSE
)

# Merge survival data with original data
combo_data <- merge(combo_treatment_synthetic_data, survival_data, by = "ID")

# Count events per treatment group correctly (each subject should only be counted once)
treatment_groups <- unique(combo_data$Treatment)
cat("\nEvents per treatment group (counting each subject only once):\n")
for (group in treatment_groups) {
  # Get data for this treatment
  group_data <- combo_data[combo_data$Treatment == group, ]
  # Get unique IDs for this treatment
  group_ids <- unique(group_data$ID)
  # Count events
  events <- sum(unique(group_data[, c("ID", "Survival_Censor")])$Survival_Censor)
  
  cat(group, ": ", events, "/", length(group_ids), " subjects\n", sep="")
}

# Run the survival_statistics function
cat("\nRunning survival_statistics function...\n")
results <- survival_statistics(
  df = combo_data,
  time_column = "Survival_Time",
  censor_column = "Survival_Censor",
  treatment_column = "Treatment",
  id_column = "ID",
  reference_group = "Control"
)

# Check the generated plots
cat("\nChecking the Events/Total counts in the results:\n")
print(results$results[, c("Group", "Events", "Total")])

# Create a fixed version of the create_km_plot function
create_km_plot_fixed <- function(df, time_column, censor_column, treatment_column, id_column = "ID") {
  # Check for required packages
  if (!requireNamespace("survminer", quietly = TRUE) || 
      !requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    message("Required packages (survminer, survival, ggplot2) not available.")
    return(NULL)
  }
  
  # Make a completely new dataframe with one row per subject
  # First, get unique subjects
  subjects <- unique(df[[id_column]])
  
  # Create a dataframe to hold subject-level data
  subject_data <- data.frame(
    id = character(length(subjects)),
    time = numeric(length(subjects)),
    status = numeric(length(subjects)),
    group = character(length(subjects)),
    stringsAsFactors = FALSE
  )
  
  # For each subject, get their last observation and event status
  for (i in seq_along(subjects)) {
    id <- subjects[i]
    subject_rows <- df[df[[id_column]] == id, ]
    
    # Sort by time to get the last observation
    subject_rows <- subject_rows[order(subject_rows[[time_column]], decreasing = TRUE), ]
    
    # Check if subject had an event (if any row has an event, consider it an event)
    had_event <- any(subject_rows[[censor_column]] == 1)
    
    # Add to subject_data
    subject_data$id[i] <- id
    subject_data$time[i] <- subject_rows[[time_column]][1] # Last observation
    subject_data$status[i] <- ifelse(had_event, 1, 0)
    subject_data$group[i] <- subject_rows[[treatment_column]][1]
  }
  
  # Convert group to factor
  subject_data$group <- factor(subject_data$group)
  
  tryCatch({
    # Fit the survival model with explicit column names
    fit <- survival::survfit(survival::Surv(time, status) ~ group, data = subject_data)
    
    # Create a plot using ggsurvplot
    base_plot <- tryCatch({
      survminer::ggsurvplot(
        fit = fit,
        data = subject_data,
        risk.table = TRUE,
        conf.int = TRUE,
        pval = TRUE
      )
    }, error = function(e) {
      message("Error in creating plot with ggsurvplot: ", e$message)
      
      # Fall back to creating a very simple plot
      fit_summary <- summary(fit)
      plot_data <- data.frame(
        time = fit_summary$time,
        surv = fit_summary$surv,
        group = rep(names(fit$strata), fit$strata)
      )
      
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = surv, color = group)) +
        ggplot2::geom_step() +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Time", y = "Survival Probability", 
                     title = "Kaplan-Meier Survival Curve")
      
      return(list(plot = p))
    })
    
    return(base_plot)
    
  }, error = function(e) {
    message("Error in survival fit or plotting: ", e$message)
    message("Data dimensions: ", nrow(subject_data), " x ", ncol(subject_data))
    return(NULL)
  })
}

# Create a fixed version of the KM plot
cat("\nCreating fixed KM plot...\n")
fixed_km_plot <- create_km_plot_fixed(
  df = combo_data,
  time_column = "Survival_Time",
  censor_column = "Survival_Censor",
  treatment_column = "Treatment",
  id_column = "ID"
)

# Save the fixed KM plot
if (!is.null(fixed_km_plot)) {
  pdf("temp/fixed_km_plot.pdf", width = 10, height = 8)
  print(fixed_km_plot)
  dev.off()
  cat("Fixed KM plot saved to temp/fixed_km_plot.pdf\n")
}

# Fix the survival_statistics function by modifying it to use the fixed create_km_plot function
cat("\nTo fix the issue in the survival_statistics function:\n")
cat("1. Replace the create_km_plot function with the fixed version that handles subject-level data\n")
cat("2. The issue is that the current implementation doesn't account for repeated measurements for the same subject\n")
cat("3. The fix ensures each subject is only counted once in the KM plot\n") 