# Data format

The data the dashboard accepts, the column-detection logic, and the validation rules.

The dashboard reads CSV and RDA files in **long format**: one row per `(animal, timepoint)` observation. This document is the canonical reference for what columns the dashboard recognizes, how it auto-detects them, and what it does when columns are missing or malformed.

If you're a user uploading data: read sections 1-3.
If you're adding a new analysis module: read sections 4-5.

---

## 1. File format

| Extension | Notes |
|---|---|
| `.csv` | UTF-8 encoded. Comma separator. First row is the header. Quoted strings handled per `utils::read.csv` defaults |
| `.rda` | A serialized R data frame. The file must contain exactly one object; if multiple, the first data frame is used |
| `.xlsx`, `.xls` | Not supported. Convert to CSV first |
| `.tsv` | Not supported. Convert to CSV first |

Max file size is configured at the Shiny server level (`shiny.maxRequestSize`). The default is 5 MB; the production VPS sets it to 100 MB.

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

When the time column contains dates (not numeric days), the dashboard converts to study day via `mouseExperiment::calculate_dates()`. Two methods:

| `time_method` | Behavior |
|---|---|
| `"direct"` | The Day column is already numeric and is used as-is |
| `"computed"` | The Date column is parsed; study day = days since the reference date (earliest date in the dataset, unless overridden) |

Date detection (`detect_date_format()` in `R/helpers.R`) tries multiple formats:

| Pattern | Example | Notes |
|---|---|---|
| ISO 8601 | `2025-03-15` | Default, unambiguous |
| US (`m/d/y`) | `3/15/2025` | Detected when day fields ≤ 12 in all rows OR a 4-digit year is present at position 3 |
| EU (`d/m/y`) | `15/3/2025` | Detected when day fields exceed 12 anywhere |
| Excel serial | `45366` | Numeric date that fits the Excel epoch (1899-12-30) |

If detection is ambiguous (all rows have day ≤ 12), the dashboard prompts you to confirm format in the column-mapping dropdowns.

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
| TG AUC | NA rows interpolated within the trapezoidal integration if `extrapolation_points > 0` |
| Survival | NA in `Time` or `Event` drops the mouse from the survival fit |
| Toxicity | NA in weight at baseline drops the mouse (no baseline ⇒ no % change) |
| Dose-Response | NA dose drops the row (no missing-dose imputation) |

NA in the Treatment column is always fatal — the mouse can't be assigned to a group.

### Treatment-name normalization

The dashboard does not normalize treatment names automatically. `"Control"`, `"control"`, and `"CTRL"` are three different groups. If you have mixed-case data, normalize before upload:

```r
df$Treatment <- tolower(df$Treatment)
```

Some control-group detection is forgiving: when no `reference_group` is specified, the dashboard auto-picks one if any of these patterns matches (case-insensitive): `vehicle`, `control`, `IgG`, `DMSO`. If none match, the first treatment alphabetically becomes the reference.

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

`validate_uploaded_file()` and `check_required_cols()` in `R/helpers.R` enforce these at upload time. Validation failures appear in the Data Upload status panel; the analysis tabs stay disabled until validation passes.

### File-level

| Check | Behavior |
|---|---|
| File extension | Must be `.csv` or `.rda`. Other extensions rejected with the supported-list message |
| File parses | `read.csv` / `load()` must succeed. Parse errors surface verbatim in the error alert |
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

After upload, the **Subset Filter** panel on the Data Upload tab lets you exclude IDs, treatments, or cages without re-uploading. The exclusion lists are applied at the start of every analysis via `apply_subset_filter(df, rv)` (in `R/helpers.R`).

State is preserved by the Save Session / Restore Session mechanism (`R/mod_log_session.R`) but is **not** part of the F.7 cookie-keyed registry — exclusion lists are session-local.

This was added in v26.06.002 as CODE_REVIEW.md F.5.

---

## 7. Example datasets

The dashboard ships sample datasets in `inst/sample_data/` for testing without uploading your own data:

| File | Description |
|---|---|
| `combo_treatment_synthetic_data.csv` | Two-drug combination with tumor + weight + survival columns |
| `dose_levels_synthetic_data.csv` | Multiple doses of a single drug |
| `weight_synthetic_data.csv` | Body-weight trajectories |
| `necrotic_synthetic_data.csv` | Includes necrosis indicator |

Load them via the Data Upload tab's "Load Sample Data" dropdown, or read directly in R via `system.file("sample_data/combo_treatment_synthetic_data.csv", package = "mouseExperimentDashboard")`.

---

## See also

- `R/helpers.R::detect_columns()`, `detect_date_format()`, `calculate_volume()` — the actual implementations
- `R/mod_data_upload.R` — the Data tab's upload + mapping UI
- [`ARCHITECTURE.md`](ARCHITECTURE.md) — how the data flows from Data Upload to each analysis module
- Backend `mouseExperiment::calculate_volume()`, `mouseExperiment::calculate_dates()` — the actual conversion functions
