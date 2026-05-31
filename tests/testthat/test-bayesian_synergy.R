# =============================================================================
# Tests for bayesian_synergy()
#
# Uses make_synergy_additive() fixture: Control/DrugA/DrugB/DrugA+DrugB at
# day 21 only; combo designed for Bliss-neutral additivity.
# Ground truth: FE_A ≈ 0.50, FE_B ≈ 0.40, Bliss expected ≈ 0.70,
# combo ≈ 0.70 → Bliss excess ≈ 0, Loewe CI ≈ 1.29 (additive).
#
# make_synergy_synergistic() is used for the P(synergy) smoke tests.
# 2 chains × 600 iterations keeps each run < 90 s on CI.
# All tests are skipped when brms is not installed.
# =============================================================================

skip_bayes_syn <- function() skip_if_not_installed("brms")

local({
  .cached_additive <- NULL

  get_additive <- function() {
    if (is.null(.cached_additive)) {
      .cached_additive <<- suppressWarnings(suppressMessages(
        bayesian_synergy(
          make_synergy_additive(),
          volume_column    = "Volume",
          time_column      = "Day",
          treatment_column = "Treatment",
          id_column        = "ID",
          drug_a_name      = "DrugA",
          drug_b_name      = "DrugB",
          combo_name       = "DrugA+DrugB",
          control_name     = "Control",
          transform        = "log",
          prior_strength   = "skeptical",
          n_chains         = 2L,
          n_iter           = 600L,
          seed             = 42L,
          return_model     = TRUE,
          plots            = TRUE,
          verbose          = FALSE
        )
      ))
    }
    .cached_additive
  }

  # ── Return structure ────────────────────────────────────────────────────────
  test_that("bayesian_synergy: returns list with expected fields", {
    skip_bayes_syn()
    res      <- get_additive()
    required <- c(
      "model_type_used", "model", "transform_used",
      "tgi_summary", "bliss_summary", "loewe_summary",
      "synergy_table", "posterior_summary", "mcmc_diagnostics",
      "summary", "synergy_plot", "posterior_dist_plot"
    )
    expect_true(
      all(required %in% names(res)),
      info = paste(
        "Missing:",
        paste(setdiff(required, names(res)), collapse = ", ")
      )
    )
  })

  test_that("bayesian_synergy: model_type_used is 'bayes_synergy'", {
    skip_bayes_syn()
    expect_equal(get_additive()$model_type_used, "bayes_synergy")
  })

  test_that("bayesian_synergy: transform_used is 'log'", {
    skip_bayes_syn()
    expect_equal(get_additive()$transform_used, "log")
  })

  test_that("bayesian_synergy: model is brmsfit when return_model = TRUE", {
    skip_bayes_syn()
    expect_s3_class(get_additive()$model, "brmsfit")
  })

  # ── tgi_summary ────────────────────────────────────────────────────────────
  test_that("bayesian_synergy: tgi_summary has one row per group", {
    skip_bayes_syn()
    tgi <- get_additive()$tgi_summary
    expect_s3_class(tgi, "data.frame")
    expect_equal(nrow(tgi), 4L)
    expect_true(
      all(c("Control", "DrugA", "DrugB", "DrugA+DrugB") %in% tgi$Group)
    )
  })

  test_that("bayesian_synergy: tgi_summary has required columns", {
    skip_bayes_syn()
    tgi <- get_additive()$tgi_summary
    expect_true(
      all(c("Group", "TGI_Median", "TGI_Lower", "TGI_Upper") %in%
            colnames(tgi))
    )
  })

  test_that("bayesian_synergy: Control TGI_Median is near zero", {
    skip_bayes_syn()
    ctrl <- get_additive()$tgi_summary
    ctrl_tgi <- ctrl$TGI_Median[ctrl$Group == "Control"]
    expect_lt(abs(ctrl_tgi), 0.05)
  })

  test_that("bayesian_synergy: treated groups have positive TGI", {
    skip_bayes_syn()
    tgi <- get_additive()$tgi_summary
    trt <- tgi[tgi$Group != "Control", ]
    expect_true(all(trt$TGI_Median > 0))
  })

  # ── bliss_summary ──────────────────────────────────────────────────────────
  test_that("bayesian_synergy: bliss_summary has required fields", {
    skip_bayes_syn()
    bs <- get_additive()$bliss_summary
    required <- c(
      "Excess_Median", "Excess_Lower", "Excess_Upper",
      "P_Synergy", "Expected_FE_Median", "Observed_FE_Median"
    )
    expect_true(all(required %in% names(bs)))
  })

  test_that("bayesian_synergy: Bliss P_Synergy is in [0, 1]", {
    skip_bayes_syn()
    p <- get_additive()$bliss_summary$P_Synergy
    expect_gte(p, 0)
    expect_lte(p, 1)
  })

  test_that("bayesian_synergy: Bliss CrI brackets the median", {
    skip_bayes_syn()
    bs <- get_additive()$bliss_summary
    expect_lt(bs$Excess_Lower, bs$Excess_Median)
    expect_gt(bs$Excess_Upper, bs$Excess_Median)
  })

  # ── loewe_summary ──────────────────────────────────────────────────────────
  test_that("bayesian_synergy: loewe_summary has required fields", {
    skip_bayes_syn()
    ls <- get_additive()$loewe_summary
    required <- c(
      "CI_Median", "CI_Lower", "CI_Upper", "P_Synergy", "Interpretation"
    )
    expect_true(all(required %in% names(ls)))
  })

  test_that("bayesian_synergy: Loewe CI_Median is positive and finite", {
    skip_bayes_syn()
    ci <- get_additive()$loewe_summary$CI_Median
    expect_true(is.finite(ci) && ci > 0)
  })

  test_that("bayesian_synergy: Loewe P_Synergy is in [0, 1]", {
    skip_bayes_syn()
    p <- get_additive()$loewe_summary$P_Synergy
    expect_gte(p, 0)
    expect_lte(p, 1)
  })

  test_that("bayesian_synergy: Loewe Interpretation is a non-empty string", {
    skip_bayes_syn()
    interp <- get_additive()$loewe_summary$Interpretation
    expect_type(interp, "character")
    expect_gt(nchar(interp), 0L)
  })

  # ── synergy_table ──────────────────────────────────────────────────────────
  test_that("bayesian_synergy: synergy_table has 6 rows", {
    skip_bayes_syn()
    st <- get_additive()$synergy_table
    expect_s3_class(st, "data.frame")
    expect_equal(nrow(st), 6L)
  })

  test_that("bayesian_synergy: synergy_table includes expected rows", {
    skip_bayes_syn()
    st <- get_additive()$synergy_table
    expect_true(all(c("Bliss Expected", "Loewe Expected") %in% st$Group))
  })

  test_that("bayesian_synergy: synergy_table has Type column", {
    skip_bayes_syn()
    st <- get_additive()$synergy_table
    expect_true("Type" %in% colnames(st))
    expect_true(all(st$Type %in% c("Observed", "Expected")))
  })

  # ── MCMC diagnostics ───────────────────────────────────────────────────────
  test_that("bayesian_synergy: convergence per current Stan recommendations", {
    skip_bayes_syn()
    diag <- get_additive()$mcmc_diagnostics
    expect_s3_class(diag, "data.frame")
    expect_true("Rhat" %in% colnames(diag))
    expect_true("Converged" %in% colnames(diag))
    expect_type(diag$Converged, "logical")
    expect_true(all(diag$Rhat <= 1.01, na.rm = TRUE),
                info = "Rhat above 1.01 (Vehtari 2021 threshold)")
    expect_true(all(diag$Bulk_ESS >= 400, na.rm = TRUE),
                info = "Bulk_ESS below 400 (Stan recommendation)")
    expect_true(all(diag$Tail_ESS >= 400, na.rm = TRUE),
                info = "Tail_ESS below 400 (Stan recommendation)")
  })

  # ── summary metadata ───────────────────────────────────────────────────────
  test_that("bayesian_synergy: summary metadata is correct", {
    skip_bayes_syn()
    s <- get_additive()$summary
    expect_equal(s$data_description$control_group, "Control")
    expect_equal(s$data_description$drug_a, "DrugA")
    expect_equal(s$data_description$drug_b, "DrugB")
    expect_equal(s$data_description$combo, "DrugA+DrugB")
    expect_equal(s$model_specification$prior_strength, "skeptical")
  })

  # ── plots ──────────────────────────────────────────────────────────────────
  test_that("bayesian_synergy: synergy_plot is ggplot (plots = TRUE)", {
    skip_bayes_syn()
    expect_s3_class(get_additive()$synergy_plot, "gg")
  })

  test_that("bayesian_synergy: posterior_dist_plot is ggplot (plots = TRUE)", {
    skip_bayes_syn()
    expect_s3_class(get_additive()$posterior_dist_plot, "gg")
  })

  # ── Input validation ───────────────────────────────────────────────────────
  test_that("bayesian_synergy: errors on missing column", {
    skip_bayes_syn()
    df <- make_synergy_additive()
    df$Volume <- NULL
    expect_error(
      bayesian_synergy(df, drug_a_name = "DrugA", drug_b_name = "DrugB",
                       combo_name = "DrugA+DrugB"),
      "Missing required columns"
    )
  })

  test_that("bayesian_synergy: errors on missing group", {
    skip_bayes_syn()
    expect_error(
      bayesian_synergy(
        make_synergy_additive(),
        drug_a_name  = "DrugA",
        drug_b_name  = "DrugB",
        combo_name   = "NoSuchGroup",
        n_chains = 2L, n_iter = 200L, plots = FALSE, verbose = FALSE
      ),
      "not found"
    )
  })
})

# ── return_model = FALSE ────────────────────────────────────────────────────

test_that("bayesian_synergy: model is NULL when return_model = FALSE", {
  skip_bayes_syn()
  res <- suppressWarnings(suppressMessages(
    bayesian_synergy(
      make_synergy_additive(),
      drug_a_name  = "DrugA",
      drug_b_name  = "DrugB",
      combo_name   = "DrugA+DrugB",
      transform    = "log",
      n_chains     = 2L, n_iter = 600L, seed = 1L,
      return_model = FALSE, plots = FALSE, verbose = FALSE
    )
  ))
  expect_null(res$model)
})

# ── plots = FALSE ───────────────────────────────────────────────────────────

test_that("bayesian_synergy: plots are NULL when plots = FALSE", {
  skip_bayes_syn()
  res <- suppressWarnings(suppressMessages(
    bayesian_synergy(
      make_synergy_additive(),
      drug_a_name  = "DrugA",
      drug_b_name  = "DrugB",
      combo_name   = "DrugA+DrugB",
      transform    = "log",
      n_chains     = 2L, n_iter = 600L, seed = 2L,
      return_model = FALSE, plots = FALSE, verbose = FALSE
    )
  ))
  expect_null(res$synergy_plot)
  expect_null(res$posterior_dist_plot)
})
