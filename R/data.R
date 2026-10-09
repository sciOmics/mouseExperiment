# Package datasets (CODE_REVIEW.md R20.52). Built by data-raw/create_datasets.R
# from the CSVs in data-raw/, which are also the dashboard's demo files;
# master_synthetic_data has its own generator, data-raw/create_master_dataset.R.

#' Combination treatment synthetic data
#'
#' Caliper measurements for a two-drug combination study, laid out as a
#' laboratory records them: calendar dates rather than study days, ear-punch
#' codes rather than unique animal numbers, and length and width rather than
#' volume. It is the dashboard's "Combo Treatment" demo.
#'
#' @format A data frame with 458 rows and 7 variables:
#' \describe{
#'   \item{Date}{Measurement date, "MM/DD/YYYY" (24 March to 25 April 2025).
#'     Convert with [calculate_dates()], start date "03/24/2025".}
#'   \item{Cage}{Cage number, 1 to 8; two cages per arm.}
#'   \item{ID}{Ear-punch code ("L", "R", "LR", "N"). It repeats in every
#'     cage, so an animal is identified by treatment, cage and ID together.}
#'   \item{Treatment}{"Control", "aPD1", "HDACi" or "HDACi + PD1"; 8 animals each.}
#'   \item{Length, Width}{Caliper measurements (mm). 0 means no palpable
#'     tumour yet. Volume comes from [calculate_volume()].}
#'   \item{Survival_Censor}{1 on the last row of an animal removed at the
#'     tumour endpoint, otherwise 0.}
#' }
#' @source Synthetic; `data-raw/combo_treatment_synthetic_data.csv`.
#' @seealso [master_synthetic_data], which has every column the package uses.
#' @examples
#' d <- calculate_dates(combo_treatment_synthetic_data, start_date = "03/24/2025")
#' d <- calculate_volume(d)
#' head(d)
"combo_treatment_synthetic_data"

#' Dose levels synthetic data
#'
#' Caliper measurements for one drug at three doses and a vehicle control. It
#' is the dashboard's "Dose Levels" demo.
#'
#' @format A data frame with 172 rows and 8 variables:
#' \describe{
#'   \item{Date}{Measurement date, "DD-Mon" without a year (24 March to
#'     18 April). Convert with [calculate_dates()], start date "24-Mar", year
#'     2025.}
#'   \item{Cage}{Cage number, 1 to 4; one cage per arm.}
#'   \item{ID}{Ear-punch code, repeated in every cage.}
#'   \item{Treatment}{"Control" or "HDACi".}
#'   \item{Dose}{0 for Control; 10, 50 or 100 for HDACi.}
#'   \item{Length, Width}{Caliper measurements (mm); 0 means no palpable
#'     tumour yet.}
#'   \item{Survival_Censor}{1 on the last row of an animal removed at the
#'     tumour endpoint, otherwise 0.}
#' }
#' @source Synthetic; `data-raw/dose_levels_synthetic_data.csv`.
#' @examples
#' d <- calculate_dates(dose_levels_synthetic_data, start_date = "24-Mar", year = 2025)
#' d <- calculate_volume(d)
#' table(d$Treatment, d$Dose)
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
#'     Independence (excess ~+8 %), i.e. synergistic;
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

#' Body-weight synthetic data
#'
#' Tumour volumes with body weights, for the toxicity functions
#' ([analyze_body_weight()], [therapeutic_window_metric()]). The column names
#' contain spaces, as exported spreadsheets often do. It is the dashboard's
#' "Weight" demo.
#'
#' @format A data frame with 100 rows and 5 variables:
#' \describe{
#'   \item{Date}{Measurement date, "Month D" without a year (1 to 15 March).
#'     Convert with [calculate_dates()], start date "March 1", year 2025.}
#'   \item{Ear Tag}{Animal number 1 to 5, repeated in every arm, so an animal
#'     is identified by treatment and ear tag together.}
#'   \item{Treatment}{"DMSO", "Drug A", "Drug B" or "Drug A + B".}
#'   \item{Tumor Volume}{Tumour volume (mm³), already calculated.}
#'   \item{Mouse Weight}{Body weight (g).}
#' }
#' @source Synthetic; `data-raw/weight_synthetic_data.csv`.
#' @examples
#' d <- calculate_dates(weight_synthetic_data, start_date = "March 1", year = 2025)
#' head(d)
"weight_synthetic_data"
