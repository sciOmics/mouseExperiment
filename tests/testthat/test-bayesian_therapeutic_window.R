# =============================================================================
# Tests for bayesian_therapeutic_window()
#
# Fixtures:
#   make_tg_for_twm()  — Control log_slope=0.12, TreatmentA log_slope=0.03
#                        TreatmentA has positive TGI (>0)
#   make_bw_simple()   — Control stable, TreatmentA loses −0.08 g/day
#
# Both models fitted with 2 chains × 600 iterations and cached in a local()
# so heavy MCMC runs happen only once per test session.
# =============================================================================

# brms, bayesplot, gamm4 and mgcv are hard Imports as of v0.10.0, so this
# helper can no longer skip. Retained as a no-op because the call sites are
# numerous, and because a skip here is exactly what let bayesian_synergy()
# stay broken for five releases (CODE_REVIEW.md R3-L).
skip_bayes_twm <- function() invisible(NULL)

# Fixture: TG data where TreatmentA inhibits tumors (TGI > 0)
make_tg_for_twm <- function() {
  set.seed(7)
  days <- c(0, 7, 14, 21)
  make_mouse <- function(id, cage, treatment, log_slope) {
    noise <- rnorm(4, 0, 0.05)
    vol   <- exp(log(200) + log_slope * days + noise)
    data.frame(
      ID = id, Cage = cage, Treatment = treatment,
      Day = days, Volume = round(vol, 2),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    make_mouse("C01", "C1", "Control",    0.12),
    make_mouse("C02", "C1", "Control",    0.12),
    make_mouse("C03", "C2", "Control",    0.12),
    make_mouse("C04", "C2", "Control",    0.12),
    make_mouse("C05", "C3", "Control",    0.12),
    make_mouse("T01", "T1", "TreatmentA", 0.03),
    make_mouse("T02", "T1", "TreatmentA", 0.03),
    make_mouse("T03", "T2", "TreatmentA", 0.03),
    make_mouse("T04", "T2", "TreatmentA", 0.03),
    make_mouse("T05", "T3", "TreatmentA", 0.03)
  )
}

local({
  .tg_result  <- NULL
  .bw_result  <- NULL
  .cached_twm <- NULL

  get_tg_result <- function() {
    if (is.null(.tg_result)) {
      .tg_result <<- suppressWarnings(suppressMessages(
        bayesian_tumor_growth(
          make_tg_for_twm(),
          treatment_column = "Treatment",
          time_column      = "Day",
          id_column        = "ID",
          reference_group  = "Control",
          transform        = "log",
          prior_strength   = "skeptical",
          n_chains         = 2L,
          n_iter           = 600L,
          seed             = 42L,
          return_model     = TRUE,
          plots            = FALSE,
          verbose          = FALSE
        )
      ))
    }
    .tg_result
  }

  get_bw_result <- function() {
    if (is.null(.bw_result)) {
      .bw_result <<- suppressWarnings(suppressMessages(
        bayesian_body_weight(
          make_bw_simple(),
          weight_column    = "Weight",
          time_column      = "Day",
          treatment_column = "Treatment",
          id_column        = "ID",
          reference_group  = "Control",
          prior_strength   = "skeptical",
          n_chains         = 2L,
          n_iter           = 600L,
          seed             = 42L,
          return_model     = TRUE,
          plots            = FALSE,
          verbose          = FALSE
        )
      ))
    }
    .bw_result
  }

  get_twm <- function() {
    if (is.null(.cached_twm)) {
      .cached_twm <<- suppressWarnings(suppressMessages(
        bayesian_therapeutic_window(
          tg_result = get_tg_result(),
          bw_result = get_bw_result(),
          plots     = TRUE,
          verbose   = FALSE
        )
      ))
    }
    .cached_twm
  }

  # ── Return structure ────────────────────────────────────────────────────────
  test_that("bayesian_therapeutic_window: expected fields returned", {
    skip_bayes_twm()
    res      <- get_twm()
    required <- c(
      "model_type_used", "twm_table", "tgi_summary",
      "wl_summary", "summary", "twm_plot", "tgi_wl_plot"
    )
    expect_true(
      all(required %in% names(res)),
      info = paste(
        "Missing:",
        paste(setdiff(required, names(res)), collapse = ", ")
      )
    )
  })

  test_that("bayesian_therapeutic_window: model_type_used is 'bayes_twm'", {
    skip_bayes_twm()
    expect_equal(get_twm()$model_type_used, "bayes_twm")
  })

  # ── twm_table ──────────────────────────────────────────────────────────────
  test_that("bayesian_therapeutic_window: twm_table has required columns", {
    skip_bayes_twm()
    tt           <- get_twm()$twm_table
    required_cols <- c(
      "Group", "TWM_Median", "TWM_Lower", "TWM_Upper",
      "TGI_Median", "WL_Pct_Median", "N_Draws"
    )
    expect_s3_class(tt, "data.frame")
    expect_true(
      all(required_cols %in% colnames(tt)),
      info = paste(
        "Missing:",
        paste(setdiff(required_cols, colnames(tt)), collapse = ", ")
      )
    )
  })

  test_that("bayesian_therapeutic_window: twm_table excludes reference", {
    skip_bayes_twm()
    tt <- get_twm()$twm_table
    expect_false("Control" %in% tt$Group)
    expect_true("TreatmentA" %in% tt$Group)
  })

  test_that("bayesian_therapeutic_window: twm_table has one row per treated", {
    skip_bayes_twm()
    expect_equal(nrow(get_twm()$twm_table), 1L)
  })

  test_that("bayesian_therapeutic_window: TWM_Median is finite", {
    skip_bayes_twm()
    expect_true(is.finite(get_twm()$twm_table$TWM_Median))
  })

  test_that("bayesian_therapeutic_window: TGI_Median > 0 for inhibited group", {
    skip_bayes_twm()
    expect_gt(get_twm()$twm_table$TGI_Median, 0)
  })

  test_that("bayesian_therapeutic_window: TWM CrI brackets the median", {
    skip_bayes_twm()
    tt <- get_twm()$twm_table
    expect_lt(tt$TWM_Lower, tt$TWM_Median)
    expect_gt(tt$TWM_Upper, tt$TWM_Median)
  })

  test_that("bayesian_therapeutic_window: N_Draws is a positive integer", {
    skip_bayes_twm()
    nd <- get_twm()$twm_table$N_Draws
    expect_true(is.numeric(nd) && nd > 0L && nd == floor(nd))
  })

  # ── tgi_summary ────────────────────────────────────────────────────────────
  test_that("bayesian_therapeutic_window: tgi_summary one row per TG group", {
    skip_bayes_twm()
    tgi <- get_twm()$tgi_summary
    expect_s3_class(tgi, "data.frame")
    expect_equal(nrow(tgi), 2L)
    expect_true(all(c("Control", "TreatmentA") %in% tgi$Group))
  })

  test_that("bayesian_therapeutic_window: tgi_summary has correct columns", {
    skip_bayes_twm()
    tgi <- get_twm()$tgi_summary
    expect_true(
      all(c("Group", "TGI_Median", "TGI_Lower", "TGI_Upper") %in%
            colnames(tgi))
    )
  })

  test_that("bayesian_therapeutic_window: Control TGI_Median near zero", {
    skip_bayes_twm()
    ctrl_tgi <- get_twm()$tgi_summary$TGI_Median[
      get_twm()$tgi_summary$Group == "Control"
    ]
    expect_lt(abs(ctrl_tgi), 0.05)
  })

  # ── wl_summary ────────────────────────────────────────────────────────────
  test_that("bayesian_therapeutic_window: wl_summary one row per BW group", {
    skip_bayes_twm()
    wl <- get_twm()$wl_summary
    expect_s3_class(wl, "data.frame")
    expect_equal(nrow(wl), 2L)
    expect_true(all(c("Control", "TreatmentA") %in% wl$Group))
  })

  test_that("bayesian_therapeutic_window: wl_summary has correct columns", {
    skip_bayes_twm()
    wl <- get_twm()$wl_summary
    expect_true(
      all(c("Group", "WL_Median", "WL_Lower", "WL_Upper") %in% colnames(wl))
    )
  })

  test_that("bayesian_therapeutic_window: TreatmentA WL_Median is negative", {
    skip_bayes_twm()
    wl     <- get_twm()$wl_summary
    trt_wl <- wl$WL_Median[wl$Group == "TreatmentA"]
    expect_lt(trt_wl, 0)
  })

  # ── summary metadata ───────────────────────────────────────────────────────
  test_that("bayesian_therapeutic_window: summary metadata is correct", {
    skip_bayes_twm()
    s <- get_twm()$summary
    expect_equal(s$parameters$reference_group, "Control")
    expect_equal(s$parameters$noise_floor, 1.0)
    expect_equal(s$parameters$common_groups, "TreatmentA")
    expect_equal(s$tg_model$endpoint_day, 21)
    expect_equal(s$bw_model$first_day, 0)
    expect_equal(s$bw_model$last_day, 28)
  })

  # ── plots ──────────────────────────────────────────────────────────────────
  test_that("bayesian_therapeutic_window: twm_plot is ggplot (plots=TRUE)", {
    skip_bayes_twm()
    expect_s3_class(get_twm()$twm_plot, "gg")
  })

  test_that("bayesian_therapeutic_window: tgi_wl_plot is ggplot (plots=TRUE)", {
    skip_bayes_twm()
    expect_s3_class(get_twm()$tgi_wl_plot, "gg")
  })

  # ── Input validation ───────────────────────────────────────────────────────
  test_that("bayesian_therapeutic_window: errors when TG model is NULL", {
    skip_bayes_twm()
    bad_tg       <- get_tg_result()
    bad_tg$model <- NULL
    expect_error(
      bayesian_therapeutic_window(bad_tg, get_bw_result()),
      "brmsfit"
    )
  })

  test_that("bayesian_therapeutic_window: errors on invalid reference group", {
    skip_bayes_twm()
    expect_error(
      bayesian_therapeutic_window(
        get_tg_result(), get_bw_result(),
        reference_group = "NoSuchGroup"
      ),
      "not found"
    )
  })
})

# ── plots = FALSE ───────────────────────────────────────────────────────────

test_that("bayesian_therapeutic_window: plots NULL when plots = FALSE", {
  skip_bayes_twm()
  tg_r <- suppressWarnings(suppressMessages(
    bayesian_tumor_growth(
      make_tg_for_twm(),
      reference_group = "Control", transform = "log",
      prior_strength  = "skeptical",
      n_chains = 2L, n_iter = 600L, seed = 10L,
      return_model = TRUE, plots = FALSE, verbose = FALSE
    )
  ))
  bw_r <- suppressWarnings(suppressMessages(
    bayesian_body_weight(
      make_bw_simple(),
      reference_group = "Control", prior_strength = "skeptical",
      n_chains = 2L, n_iter = 600L, seed = 10L,
      return_model = TRUE, plots = FALSE, verbose = FALSE
    )
  ))
  res <- suppressWarnings(suppressMessages(
    bayesian_therapeutic_window(tg_r, bw_r, plots = FALSE, verbose = FALSE)
  ))
  expect_null(res$twm_plot)
  expect_null(res$tgi_wl_plot)
})
