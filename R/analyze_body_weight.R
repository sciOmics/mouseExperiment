#' Analyze Body Weight via Mixed Modeling
#'
#' Fits a linear mixed-effects model to longitudinal body weight data,
#' optionally adjusting for estimated tumor weight, sex, initial mass,
#' and tumor volume as covariates.
#'
#' @param df Data frame with longitudinal data.
#' @param weight_column Name of the body weight column (grams).
#' @param time_column Name of the time/day column.
#' @param treatment_column Name of the treatment group column.
#' @param id_column Name of the mouse/subject ID column.
#' @param volume_column Name of the tumor volume column (mm³). NULL to skip tumor adjustment.
#' @param sex_column Name of the sex column. NULL to omit.
#' @param cage_column Name of the cage column. NULL to omit.
#' @param adjust_tumor_weight Logical; subtract estimated tumor weight from body mass.
#' @param tumor_density Density in g/cm³ for tumor weight estimation (default 1.0).
#' @param covariates Character vector of optional covariates: "volume", "sex", "initial_mass".
#' @param estimation Character; "REML" (default) or "ML".
#' @param reference_group Name of the control/reference group. NULL auto-selects.
#' @return A list with components: model, fixed_effects, random_effects, emmeans_table,
#'   model_info, weight_data, summary_text.
#' @export
analyze_body_weight <- function(df,
                                weight_column    = "Weight",
                                time_column      = "Day",
                                treatment_column = "Treatment",
                                id_column        = "ID",
                                volume_column    = NULL,
                                sex_column       = NULL,
                                cage_column      = NULL,
                                adjust_tumor_weight = TRUE,
                                tumor_density    = 1.0,
                                covariates       = c("volume"),
                                estimation       = c("REML", "ML"),
                                reference_group  = NULL) {

  estimation <- match.arg(estimation)

  # --- Validate required columns ---
  required <- c(weight_column, time_column, treatment_column, id_column)
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # --- Prepare working data ---
  wd <- data.frame(
    ID        = as.factor(df[[id_column]]),
    Treatment = as.factor(df[[treatment_column]]),
    Day       = as.numeric(df[[time_column]]),
    Weight    = as.numeric(df[[weight_column]]),
    stringsAsFactors = FALSE
  )

  # Tumor volume (for adjustment and/or covariate)
  has_volume <- !is.null(volume_column) && volume_column %in% names(df)
  if (has_volume) {
    wd$Volume <- as.numeric(df[[volume_column]])
  }

  # Tumor weight adjustment
  if (adjust_tumor_weight && has_volume) {
    # Volume in mm³ → cm³ (÷1000) → grams (× density)
    wd$Tumor_Weight <- wd$Volume / 1000 * tumor_density
    wd$Net_Weight   <- wd$Weight - wd$Tumor_Weight
    response_col    <- "Net_Weight"
  } else {
    wd$Net_Weight   <- wd$Weight
    response_col    <- "Net_Weight"
  }

  # Sex covariate
  has_sex <- !is.null(sex_column) && sex_column %in% names(df)
  if (has_sex) {
    wd$Sex <- as.factor(df[[sex_column]])
  }

  # Cage
  has_cage <- !is.null(cage_column) && cage_column %in% names(df)
  if (has_cage) {
    wd$Cage <- as.factor(df[[cage_column]])
  }

  # Initial mass per mouse (baseline weight at earliest timepoint)
  wd <- wd[order(wd$ID, wd$Day), ]
  baseline <- stats::aggregate(Net_Weight ~ ID, data = wd,
                               FUN = function(x) x[1])
  names(baseline)[2] <- "Initial_Mass"
  wd <- merge(wd, baseline, by = "ID", all.x = TRUE)

  # Remove rows with missing response
  wd <- wd[!is.na(wd$Net_Weight) & !is.na(wd$Day), ]
  if (nrow(wd) == 0) stop("No valid weight observations after cleaning.")

  # Set reference group
  if (!is.null(reference_group) && reference_group %in% levels(wd$Treatment)) {
    wd$Treatment <- stats::relevel(wd$Treatment, ref = reference_group)
  }

  # --- Build formula ---
  fixed_terms <- "Treatment * Day"
  if ("volume" %in% covariates && has_volume) {
    fixed_terms <- paste(fixed_terms, "+ Volume")
  }
  if ("sex" %in% covariates && has_sex) {
    fixed_terms <- paste(fixed_terms, "+ Sex")
  }
  if ("initial_mass" %in% covariates) {
    fixed_terms <- paste(fixed_terms, "+ Initial_Mass")
  }

  # --- Fit model with auto-fallback ---
  model <- NULL
  model_simplified <- FALSE

  # Try random slope + intercept first
  formula_full <- stats::as.formula(
    paste(response_col, "~", fixed_terms, "+ (1 + Day | ID)")
  )
  tryCatch({
    model <- lme4::lmer(formula_full, data = wd, REML = (estimation == "REML"))
    # Check for singular fit
    if (lme4::isSingular(model)) {
      model <- NULL
    }
  }, error = function(e) {
    model <<- NULL
  }, warning = function(w) {
    if (grepl("singular|converge", w$message, ignore.case = TRUE)) {
      model <<- NULL
    }
  })

  # Fallback to random intercept only
  if (is.null(model)) {
    formula_simple <- stats::as.formula(
      paste(response_col, "~", fixed_terms, "+ (1 | ID)")
    )
    tryCatch({
      model <- lme4::lmer(formula_simple, data = wd, REML = (estimation == "REML"))
      model_simplified <- TRUE
    }, error = function(e) {
      stop("Mixed model failed to converge: ", conditionMessage(e))
    })
  }

  # --- Extract results ---
  model_summary <- summary(model)
  fe <- as.data.frame(stats::coef(model_summary))
  fe$Term <- rownames(fe)
  fe <- fe[, c("Term", setdiff(names(fe), "Term"))]
  rownames(fe) <- NULL

  re <- as.data.frame(lme4::VarCorr(model))

  # Marginal means per treatment
  emm <- tryCatch({
    em <- emmeans::emmeans(model, ~ Treatment)
    as.data.frame(em)
  }, error = function(e) NULL)

  # --- Summary text ---
  lines <- c(
    "=== BODY WEIGHT MIXED MODEL ===",
    "",
    sprintf("Response: %s", if (adjust_tumor_weight && has_volume) "Net Weight (body - tumor)" else "Body Weight"),
    sprintf("Estimation: %s", estimation),
    sprintf("Random effects: %s", if (model_simplified) "(1 | ID) [simplified]" else "(1 + Day | ID)"),
    sprintf("Fixed effects: %s", fixed_terms),
    sprintf("Observations: %d  |  Subjects: %d  |  Groups: %d",
            nrow(wd), length(unique(wd$ID)), length(levels(wd$Treatment))),
    ""
  )

  if (model_simplified) {
    lines <- c(lines,
      "NOTE: Model simplified due to sample size. Random slope removed; using random intercept only.",
      "")
  }

  summary_text <- paste(lines, collapse = "\n")

  list(
    model          = model,
    fixed_effects  = fe,
    random_effects = re,
    emmeans_table  = emm,
    model_info     = list(
      estimation       = estimation,
      response         = response_col,
      fixed_formula    = fixed_terms,
      model_simplified = model_simplified,
      adjust_tumor     = adjust_tumor_weight && has_volume,
      tumor_density    = tumor_density,
      n_obs            = nrow(wd),
      n_subjects       = length(unique(wd$ID)),
      n_groups         = length(levels(wd$Treatment))
    ),
    weight_data    = wd,
    summary_text   = summary_text
  )
}
