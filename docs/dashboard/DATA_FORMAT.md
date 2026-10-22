# Data format

The data the dashboard accepts, the column-detection logic, and the validation rules.

The dashboard reads CSV files in **long format**: one row per `(animal, timepoint)` observation. This document is the canonical reference for what columns the dashboard recognizes, how it auto-detects them, and what it does when columns are missing or malformed.

If you're a user uploading data: read sections 1-3.
If you're adding a new analysis module: read sections 4-5.

---

## 1. File format

| Extension | Notes |
|---|---|
| `.csv` | UTF-8 encoded. Comma separator. First row is the header. Quoted strings handled per `utils::read.csv` defaults |
| `.rda`, `.RData` | Not supported (removed in v26.09.001). Export the data frame with `write.csv()` first |
| `.xlsx`, `.xls` | Not supported. Convert to CSV first |
| `.tsv` | Not supported. Convert to CSV first |

Files up to 30 MB are accepted (`run_dashboard()` sets `shiny.maxRequestSize`). Comma-, semicolon- and tab-separated files are read, with a decimal comma when the delimiter is a semicolon or tab, in UTF-8, UTF-16 or Latin-1. A row whose number of fields differs from the header is reported with its line number; quote any value that contains the delimiter.

### Wide format

The dashboard does **not** accept wide format (one column per timepoint). If your spreadsheet has columns `Day_0`, `Day_7`, `Day_14`, …, pivot to long format before uploading:

```r
library(tidyr)
long_df <- wide_df %>%
  pivot_longer(starts_with("Day_"),
               names_to = "Day",
               names_prefix = "Day_",
               values_to = "Volume") %>%
  mutate(Day = as.numeric(Day))
```

---

## 2. Required vs optional columns

The Data Upload tab maps columns interactively, so the column names in your file don't have to match a fixed set — they just need to *exist*. The auto-detection logic (section 4) makes a best guess; you can override every guess in the dropdowns.

### Required for every analysis

| Concept | What it means | Default name candidates |
|---|---|---|
| **ID** | Unique identifier per animal | `ID`, `Mouse_ID`, `Animal_ID`, `Subject`, `mouse tag`, `ear tag`, `tag` |
| **Treatment** | Treatment-group label | `Treatment`, `Group`, `Treatment_Group` |
| **Time** (one of) | Either Day numeric or Date | `Day`, `Time`, `Timepoint`, `Day_Number` (numeric) — OR — `Date`, `Measurement_Date` |
| **Outcome** | Tumor volume **or** length + width | `Volume`, `Tumor_Volume`, `Tumor Volume`, `TV` — OR — `Length`/`Width`-like columns |

### Optional (used by specific analyses)

| Concept | When you need it | Default name candidates |
|---|---|---|
| **Cage** | Mixed-effects models with cage as a random intercept | `Cage`, `Cage_ID` |
| **Dose** | Dose-response analysis | `Dose`, `Dose_Level`, `Concentration`, `Conc` |
| **Weight** | Toxicity analysis (body weight) | `Weight`, `Mass`, `Body_Weight`, `Animal_Weight`, `Animal Mass` |
| **Survival event** | Survival analysis | `Survival_Censor`, `Event`, `Censor`, `Survival_Event` |
| **Survival time** | Survival analysis with a separate event column | `Survival_Time`, `Time_To_Event`, `Days_To_Event` |
| **Necrotic** | Necrosis handling in TG (exclude / covariate / none) | `Necrotic` (binary 0/1 or TRUE/FALSE) |
| **Sex** | Toxicity covariate | `Sex`, `Gender` |

### Example minimal valid input

For a basic tumor-growth analysis:

| ID | Day | Volume | Treatment | Cage |
|---|---|---|---|---|
| Mouse_01 | 0 | 45 | Control | A1 |
| Mouse_01 | 7 | 63 | Control | A1 |
| Mouse_02 | 0 | 41 | DrugA | A2 |
| Mouse_02 | 7 | 45 | DrugA | A2 |

The same mouse appears once per timepoint (long format). Treatment and Cage repeat for each row of a given mouse (constant per animal).

---

## 3. How the dashboard handles common data shapes

### Date columns → study day conversion

When the time comes from a Date column ("Calculate from Date column" on Data Upload), `study_days()` in `R/helpers_data.R` reads the dates. Study day is the number of days since the **reference date**, which you enter as YYYY-MM-DD. `mouseExperiment::calculate_dates()` then computes Day from unambiguous YYYY-MM-DD values.

The whole column is read, not just the first row (R20.D1):

| Layout | Examples | How it is read |
|---|---|---|
| YYYY-MM-DD | `2025-03-15`, `2025/03/15`, `2025-03-15 09:30` | Unambiguous; a trailing time is ignored |
| Day and month as numbers | `15/03/2025`, `03/15/2025`, `15.03.2025`, `15-03-2025` | See below |
| Month name with year | `15-Mar-2025`, `March 15, 2025` | Unambiguous (English month names) |
| Month name, no year | `15-Mar`, `March 15` | Takes the year that puts the date nearest the reference date |

**Day and month as numbers.** A value whose first number is above 12 (`15/03/2025`) can only be DD/MM/YYYY, and one whose second number is above 12 (`03/15/2025`) can only be MM/DD/YYYY; either decides the order for the whole column. When no value decides it (every number is 12 or less), the **Date order** selector under the date column decides; left at "Detect from the data", the tab asks you to choose. A column with values that fit only DD/MM and others that fit only MM/DD is an error, and so is a declared order that the data contradict.

Under the reference date, Data Upload says how the dates were read ("Dates read as DD/MM/YYYY ('15/03/2025' can only be DD/MM/YYYY)"), and it warns when:
- values could not be read as dates (those rows get no Day and are left out);
- measurements fall before the reference date (negative days);
- year-less dates were placed in a different year (the study crosses New Year);
- the column mixes layouts.

### Volume calculation from length + width

When `Length` and `Width` columns are mapped but `Volume` isn't, the dashboard runs `mouseExperiment::calculate_volume()` with one of four formulas:

| Formula | Equation | When to use |
|---|---|---|
| Ellipsoid (default) | V = (L × W² × π) / 6 | Standard solid-tumor estimate |
| Modified Ellipsoid | V = (L × W²) / 2 | "Rounded down" ellipsoid; used in some legacy preclinical SOPs |
| Cylinder | V = (π × W² × L) / 4 | Tumors that are clearly elongated rather than spherical |
| Sphere | V = (π × W³) / 6 | When you only measured width and the tumor is roughly spherical |

Dimensions are auto-corrected so the longer measurement is always L (handles "I measured them in the wrong order" data). The pre-correction values are preserved in the original columns; only the derived Volume reflects the corrected geometry.

### Missing values

Per-row NA handling depends on the analysis:

| Analysis | NA handling |
|---|---|
| TG LME4 | NA rows dropped per-mouse; mice with < 3 observations are excluded with a warning |
| TG AUC | NA volumes are dropped. The AUC is integrated from each arm's fitted curve, so a missing measurement does not remove the animal |
| Survival | NA in `Time` or `Event` drops the mouse from the survival fit |
| Toxicity | NA in weight at baseline drops the mouse (no baseline ⇒ no % change) |
| Dose-Response | NA dose drops the row (no missing-dose imputation) |

NA in the Treatment column is always fatal — the mouse can't be assigned to a group.

### Treatment-name normalization

The dashboard does not normalize treatment names automatically. `"Control"`, `"control"`, and `"CTRL"` are three different groups. If you have mixed-case data, normalize before upload:

```r
df$Treatment <- tolower(df$Treatment)
```

The control group is chosen once, in the **Study Design** panel (below). Its default is the first treatment whose name contains `vehicle`, `control`, `IgG` or `DMSO` (case-insensitive), or else the first treatment alphabetically.

### Study design (set once on Data Upload)

The Study Design panel below the column mapping holds the settings every tab uses (R20-D.H). They replace the reference and cage selectors each tab used to have:

| Setting | What it does |
|---|---|
| **Tumour volume units** (mm³ / cm³) | The units of the volume column, or of the volume calculated from length and width (millimetres give mm³). The Toxicity tab's tumour-mass correction needs them and will not run while they are undeclared. The demos, and uploads whose volume is calculated, start at mm³. An upload with a volume column starts undeclared. |
| **Control group** | The reference arm for every comparison. With a dose column mapped, Tumour Growth, Survival and Synergy label arms "Treatment Dose", and the control maps to its own dose label (e.g. "Vehicle 0"). Each tab states the reference it used. If the control is excluded by the subset filter, the tab says which group stands in for it. |
| **Treatment was assigned to** (individual animals / whole cages) | Shown when a cage column is mapped. "Whole cages" makes Survival compare arms by a cage-level permutation log-rank and makes the tumour-growth randomisation test permute cages. With 2 cages per arm no p-value can fall below 1/3, and with 3 below 0.1. |

**Evaluable days.** An arm is evaluable on a day when at least 50 % of its enrolled animals, and at least 3, are still on study (measured on that day or later; maintainer decision R20-K). Synergy, the therapeutic window, dose-response and the tumour-growth AUC use only days on which every arm they compare is evaluable. By default they use the last such day, and their day selectors list only evaluable days. Each result names the day it used and lists the days and arms the rule excluded. The rule is fixed, not a setting; the panel states it.

**One animal, one key.** Everywhere the dashboard identifies an animal (survival records, the subset filter, the tumour-growth drill-down), it uses treatment + ID + cage. An ear tag reused in another arm or cage is therefore a different animal.

### Cage assignment

Cage should be constant per animal across all timepoints. If a single mouse appears in multiple cages across rows, the dashboard uses the most-frequent cage and emits a warning.

---

## 4. Column auto-detection logic

`detect_columns()` in `R/helpers.R` runs at upload time. For each concept, it tries a list of candidate names against the actual headers (case-insensitive). The exact lists:

| Concept | Patterns tried |
|---|---|
| ID | `id`, `mouse_id`, `animal_id`, `subject`, `mouse tag`, `ear tag`, `tag` |
| Day | `day`, `time`, `timepoint`, `day_number` |
| Volume | `volume`, `tumor_volume`, `tumor volume`, `tv` |
| Length | `length`, `tumor_length`, `tumor length` (`startsWith` fallback catches `"Tumor Length (mm)"`) |
| Width | `width`, `tumor_width`, `tumor width` |
| Treatment | `treatment`, `group`, `treatment_group` |
| Cage | `cage`, `cage_id` |
| Date | `date`, `measurement_date` |
| Dose | `dose`, `dose_level`, `concentration`, `conc` |
| Survival time | `survival_time`, `time_to_event`, `days_to_event` |
| Survival event | `survival_event`, `event`, `survival_censor`, `censor` |
| Weight | `weight`, `mass`, `body_weight`, `animal_weight`, `animal mass` |
| Sex | `sex`, `gender` |

Detection is deterministic: the first matching candidate (in the order listed) wins. If you have both `Day` and `Timepoint` columns, the dashboard picks `Day`.

The auto-detection result is always overridable in the Data Upload tab's column-mapping dropdowns — you don't have to rename your data to match.

---

## 5. Validation rules

`validate_uploaded_file()` and `check_required_cols()` in `R/helpers_data.R` enforce these at upload time. Validation failures appear in the Data Upload status panel; the analysis tabs stay disabled until validation passes.

### File-level

| Check | Behavior |
|---|---|
| File extension | Must be `.csv`. Other extensions are rejected with a message asking for a CSV file |
| File parses | `read.csv` must succeed. Parse errors surface verbatim in the error alert |
| Non-empty | Must have ≥ 1 row of data after the header |
| Headers present | First row must contain non-NA strings |

### Column-level (per analysis)

| Analysis | Required columns mapped |
|---|---|
| TG, all model types | ID, Day (or Date + time_method = "computed"), Volume (or Length + Width), Treatment |
| Survival | ID, Treatment, Time (or Date), Survival event column OR a volume threshold to derive events |
| Toxicity | ID, Day, Treatment, Weight |
| Dose-Response | ID, Day, Treatment, Volume, Dose |
| Synergy | ID, Day, Treatment, Volume; treatment groups must include the named drug A, drug B, combo, and control |
| Power Analysis | None (purely simulation-based) |

### Cell-level

| Column | Validation |
|---|---|
| ID | Non-NA, non-empty string |
| Day | Non-negative numeric (or convertible from Date) |
| Volume | Numeric ≥ 0 (negatives rejected; 0 allowed for pre-treatment baseline) |
| Weight | Numeric > 0 |
| Treatment | Non-NA, non-empty string |
| Survival event | Either `0`/`1` integer OR `TRUE`/`FALSE` logical; other values rejected |
| Date | Parseable as a date in one of the four supported formats |
| Necrotic | Coerced from `0`/`1`, `TRUE`/`FALSE`, `Yes`/`No`, or `Y`/`N` |

### Cross-row consistency

| Check | Behavior |
|---|---|
| Same ID, same Treatment across timepoints | Required. Violation rejects the upload with the row numbers |
| Same ID, same Cage across timepoints | Warning (most-frequent wins) |
| Same ID, multiple identical timepoints | Warning (duplicate measurement; the dashboard averages by default) |

---

## 6. Subset filter (post-upload)

After upload, the **Subset Filter** panel on the Data Upload tab lets you exclude animals, treatments, or cages without re-uploading. Animals are listed as "ID (arm, cage C)" and excluded by treatment + ID + cage, so excluding ear tag 3 in one arm leaves animal 3 in the other arms. The exclusion lists are applied at the start of every analysis via `apply_subset_filter(df, rv)` (in `R/helpers_data.R`).

Exclusion lists are session-local: they are not part of the F.7 cookie-keyed registry. The settings preset (format 2, v26.09.003) saves them together with the column mappings, methods and study design. Load a preset after uploading the dataset it describes.

This was added in v26.06.002 as CODE_REVIEW.md F.5.

---

## 7. Example datasets

The dashboard ships sample datasets in `inst/sample_data/` for testing without uploading your own data:

| Demo name | File | Description |
|---|---|---|
| Combo Treatment | `combo_treatment_synthetic_data.csv` | Two-drug combination with tumor + weight + survival columns |
| Dose Levels | `dose_levels_synthetic_data.csv` | Multiple doses of a single drug |
| Weight | `weight_synthetic_data.csv` | Body-weight trajectories |
| All Features | `all_features_synthetic_data.csv` | Every optional column: six arms (vehicle, a three-dose series, a second drug and their combination), length/width/height, weight, dose, cage, necrosis flag, survival event and removal reason |

Load them from the Data Upload tab's "Try Demo Data" list, or read directly in R via `system.file("sample_data/combo_treatment_synthetic_data.csv", package = "mouseExperiment")`.

---

## See also

- `R/helpers_data.R::detect_columns()`, `detect_date_format()`, `study_days()`, `apply_subset_filter()` — the actual implementations
- `R/mod_data_upload.R` — the Data tab's upload + mapping UI
- [`ARCHITECTURE.md`](ARCHITECTURE.md) — how the data flows from Data Upload to each analysis module
- Backend `mouseExperiment::calculate_volume()`, `mouseExperiment::calculate_dates()` — the actual conversion functions
