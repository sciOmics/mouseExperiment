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
