
#' Synthetic tumor growth data
#' 
#' A dataset containing synthetic tumor volume measurements over time
#' for different treatment groups.
#' 
#' @format A data frame with 120 rows and 7 variables:
#' \describe{
#'   \item{Mouse_ID}{Mouse identifier, with format M followed by a number}
#'   \item{Day}{Day of measurement, starting from day 0}
#'   \item{Treatment}{Treatment group (Control, Treatment A, Treatment B)}
#'   \item{Volume}{Tumor volume measurement}
#'   \item{Group}{Same as Treatment, an alternative name for the treatment group}
#'   \item{ID}{Numeric identifier for the mouse, without the "M" prefix}
#'   \item{Cage}{Cage number where the mouse is housed}
#' }
#' @source Synthetic data generated using random number generation
#' @examples
#' data(synthetic_data)
#' head(synthetic_data)
"synthetic_data"

#' Example data
#' 
#' A simple example dataset with x and y coordinates.
#' 
#' @format A data frame with 10 rows and 2 variables:
#' \describe{
#'   \item{x}{x-coordinate, sequence from 1 to 10}
#'   \item{y}{y-coordinate, x plus random noise}
#' }
#' @source Synthetic data generated using random number generation
#' @examples
#' data(my_data)
#' plot(my_data$x, my_data$y)
"my_data"

#' Combination treatment synthetic data
#' 
#' A dataset containing synthetic tumor volume measurements over time
#' for different combination treatment groups.
#' 
#' @format A data frame with rows for each mouse at each timepoint and 6 variables:
#' \describe{
#'   \item{Mouse_ID}{Mouse identifier, with format M followed by a number}
#'   \item{Day}{Day of measurement, starting from day 0}
#'   \item{Treatment}{Treatment group (Control, aPD1, HDACi, HDACi + PD1)}
#'   \item{Volume}{Tumor volume measurement}
#'   \item{ID}{Numeric identifier for the mouse, without the "M" prefix}
#'   \item{Cage}{Cage number where the mouse is housed}
#' }
#' @source Synthetic data generated using random number generation with treatment-specific effects
#' @examples
#' data(combo_treatment_synthetic_data)
#' head(combo_treatment_synthetic_data)
"combo_treatment_synthetic_data"

#' Combination treatment schedule
#' 
#' A dataset specifying the dosing schedule for combination treatments.
#' 
#' @format A data frame with 20 rows and 3 variables:
#' \describe{
#'   \item{Treatment}{Treatment group (Control, aPD1, HDACi, HDACi + PD1)}
#'   \item{Day}{Day of dose administration}
#'   \item{Dose}{Dose amount}
#' }
#' @source Synthetic treatment schedule
#' @examples
#' data(combo_treatment_schedule)
#' head(combo_treatment_schedule)
"combo_treatment_schedule"

#' Dose levels synthetic data
#' 
#' A dataset containing synthetic tumor volume measurements over time
#' for different dose levels of a single drug.
#' 
#' @format A data frame with rows for each mouse at each timepoint and 7 variables:
#' \describe{
#'   \item{Mouse_ID}{Mouse identifier, with format M followed by a number}
#'   \item{Day}{Day of measurement, starting from day 0}
#'   \item{Treatment}{Treatment name (always "Drug X")}
#'   \item{Dose}{Dose level (0, 10, 25, 50, 100)}
#'   \item{Volume}{Tumor volume measurement}
#'   \item{ID}{Numeric identifier for the mouse, without the "M" prefix}
#'   \item{Cage}{Cage number where the mouse is housed}
#' }
#' @source Synthetic data generated using random number generation with dose-dependent effects
#' @examples
#' data(dose_levels_synthetic_data)
#' head(dose_levels_synthetic_data)
"dose_levels_synthetic_data"

#' Master synthetic dataset — all analysis scenarios
#'
#' A single long-format data frame designed to exercise every analysis
#' function in the package.  Six treatment groups provide known significant
#' and non-significant contrasts so that tests can validate direction of
#' effect without prescribing exact p-values.
#'
#' @format A data frame with up to 672 rows (48 mice × 14 time points; mice
#'   reaching the volume endpoint have fewer rows) and 13 variables:
#' \describe{
#'   \item{Mouse_ID}{Mouse identifier, format "M001"–"M048"}
#'   \item{ID}{Integer mouse number 1–48 (numeric suffix of Mouse_ID)}
#'   \item{Day}{Study day (0, 2, 4, …, 26)}
#'   \item{Length}{Tumor caliper long-axis measurement (mm)}
#'   \item{Width}{Tumor caliper short-axis measurement (mm)}
#'   \item{Height}{Tumor caliper third-axis measurement (mm)}
#'   \item{Volume}{Ellipsoid tumor volume from Length/Width: π·L·W²/6 (mm³)}
#'   \item{Weight}{Body weight (g)}
#'   \item{Treatment}{Treatment group: "Vehicle", "Drug_A Low", "Drug_A Mid",
#'     "Drug_A High", "Drug_B", "Drug_A Mid + Drug_B"}
#'   \item{Dose}{Drug_A dose component (mg/kg): 0, 5, 15, or 30;
#'     Drug_B fixed at 10 mg/kg (see Treatment column)}
#'   \item{Cage}{Cage identifier "C01"–"C12"; 2 cages per group, 4 mice each}
#'   \item{Necrotic}{0/1; 1 once a mouse develops necrosis (persistent);
#'     occurs only in "Drug_A High" and "Drug_A Mid + Drug_B" groups}
#'   \item{Survival_Censor}{0 (alive or study-end censored) or 1 (endpoint
#'     reached — volume > 2500 mm³); equal to 1 only on the last row of
#'     a mouse that reached the endpoint}
#' }
#'
#' @section Built-in statistical properties (set.seed(2026)):
#' \describe{
#'   \item{Tumor growth LME4}{Drug_A High vs Vehicle: clearly significant
#'     (large TGI ~85 %); Drug_B vs Vehicle: borderline; Drug_A Low vs
#'     Vehicle: non-significant (small TGI ~15 %)}
#'   \item{Survival}{Vehicle has the highest event rate (~80 %); Drug_A High
#'     and combination groups have near-zero event rates}
#'   \item{Body weight}{Drug_A High causes ~10 % weight loss (significant vs
#'     Vehicle); Drug_B causes ~3 % loss (non-significant)}
#'   \item{Drug synergy}{The combination "Drug_A Mid + Drug_B" exceeds Bliss
#'     Independence (excess ~+8 %) and has Loewe CI < 1 (synergistic);
#'     evaluate at \code{eval_time_point = 10} where all mice are present}
#'   \item{Dose-response}{Drug_A at doses 0/5/15/30 mg/kg exhibits a
#'     Hill-shaped dose-response; EC50 approximately 10–15 mg/kg}
#'   \item{Necrosis}{Approximately 20–35 % of mice in high-dose groups
#'     develop persistent necrosis from day 12 onward}
#' }
#'
#' @section Recommended usage per analysis:
#' \describe{
#'   \item{Tumor growth}{Pass full dataset; all six groups are available}
#'   \item{Survival}{\code{dplyr::slice_max(Day)} to one row per mouse before
#'     passing to \code{\link{survival_statistics}}}
#'   \item{Body weight / toxicity}{Use Weight column with
#'     \code{\link{analyze_body_weight}}}
#'   \item{Drug synergy}{\code{control_name = "Vehicle"},
#'     \code{drug_a_name = "Drug_A Mid"},
#'     \code{drug_b_name = "Drug_B"},
#'     \code{combo_name = "Drug_A Mid + Drug_B"},
#'     \code{eval_time_point = 10}}
#'   \item{Dose-response}{Filter to
#'     \code{Treatment \%in\% c("Vehicle","Drug_A Low","Drug_A Mid","Drug_A High")}}
#'   \item{Volume from dimensions}{Re-derive via
#'     \code{\link{calculate_volume}(Length, Width)}}
#' }
#'
#' @source Synthetic data generated with \code{data-raw/create_master_dataset.R}
#'   using \code{set.seed(2026)}.
#' @seealso \code{\link{combo_treatment_synthetic_data}},
#'   \code{\link{dose_levels_synthetic_data}}
#' @examples
#' data(master_synthetic_data)
#' head(master_synthetic_data)
#' table(master_synthetic_data$Treatment)
#'
#' # Survival subset (one row per mouse)
#' surv_df <- master_synthetic_data |>
#'   dplyr::group_by(Mouse_ID) |>
#'   dplyr::slice_max(Day, n = 1, with_ties = FALSE) |>
#'   dplyr::ungroup()
#'
#' # Dose-response subset
#' dr_df <- subset(master_synthetic_data,
#'                 Treatment %in% c("Vehicle","Drug_A Low","Drug_A Mid","Drug_A High"))
"master_synthetic_data"

#' Dose levels treatment schedule
#' 
#' A dataset specifying the dosing schedule for different dose levels.
#' 
#' @format A data frame with 20 rows and 4 variables:
#' \describe{
#'   \item{Treatment}{Treatment name (always "Drug X")}
#'   \item{Dose}{Dose group level (0, 10, 25, 50, 100)}
#'   \item{Day}{Day of dose administration}
#'   \item{Administered_Dose}{Actual dose administered}
#' }
#' @source Synthetic treatment schedule
#' @examples
#' data(dose_levels_treatment_schedule)
#' head(dose_levels_treatment_schedule)
"dose_levels_treatment_schedule"
