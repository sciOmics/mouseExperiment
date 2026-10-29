# Column names used inside ggplot2 aes() and dplyr verbs. Declared so that
# R CMD check does not report them as undefined globals (CODE_REVIEW.md R20.50).
utils::globalVariables(c(
  "Bliss_Expected_TGI", "Group", "Mouse_ID", "TGI", "TGI_Combo", "TGI_Drug_A",
  "TGI_Drug_B", "Time", "Time_Point", "Treatment", "Type", "ci95_lower",
  "ci95_upper", "estimate", "growth_rate", "growth_rate.CI_lower",
  "growth_rate.CI_upper", "growth_rate.Mean", "growth_rate.SD",
  "growth_rate.SEM", "mean_volume", "sd_volume", "sem_volume", "term", "x",
  "xend", "y"
))

# The same, in the dashboard's plots.
utils::globalVariables(c(
  "Alpha", "Animal", "Ann_Label", "Arm", "CI_Lower", "CI_Upper", "Contrast", "Day", "Effect",
  "Dose", "Effect_Size", "Estimate", "Fitted", "HR", "Label", "Lower", "Mean_Volume",
  "N_Per_Group", "Power", "Power_Lower", "Power_Upper", "Residuals", "Survival", "Upper",
  "Volume", "ctr", "grp", "hi", "label_text", "lo", "mean_tgi", "plot_group",
  "plot_val"
))
