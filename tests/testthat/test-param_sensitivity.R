# ============================================================================
# Parameter sensitivity tests (CODE_REVIEW.md Round 2 K.4)
#
# Every silent-ignore bug the prior reviews caught (Round 1 1.1
# handle_cage_effects, Round 2 A.3 auc_method, Round 2 J.11
# analyze_body_weight cage_column, …) would have been prevented by a
# small assertion that the argument actually changes the output. This
# file installs that regression net for the parameters most at risk:
# any future "this arg is documented but silently dropped" bug fails
# loudly here.
#
# Tests are intentionally low-stakes: they assert that the *output*
# differs between two parameter values, not that a specific numeric
# value matches. That avoids brittleness while still catching the
# entire class of silent-ignore bugs.
# ============================================================================

test_that("tumor_growth_statistics: random_effects_specification actually changes the model", {
  df <- make_tg_simple()
  res_int <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df, model_type = "lme4",
      random_effects_specification = "intercept_only",
      reference_group = "Control", plots = FALSE, verbose = FALSE
    )
  ))
  res_slope <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df, model_type = "lme4",
      random_effects_specification = "slope",
      reference_group = "Control", plots = FALSE, verbose = FALSE
    )
  ))

  # Variance-component structures must differ between the two specs.
  vc_int   <- tryCatch(as.data.frame(lme4::VarCorr(res_int$model)),   error = function(e) NULL)
  vc_slope <- tryCatch(as.data.frame(lme4::VarCorr(res_slope$model)), error = function(e) NULL)
  expect_false(is.null(vc_int))
  expect_false(is.null(vc_slope))
  expect_false(identical(nrow(vc_int), nrow(vc_slope)),
               info = "intercept_only vs slope should produce different VarCorr row counts")
})

test_that("tumor_growth_statistics: transform argument actually changes pairwise estimates", {
  df <- make_tg_simple()
  res_log  <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, model_type = "lme4", transform = "log",
                            reference_group = "Control", plots = FALSE)
  ))
  res_none <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, model_type = "lme4", transform = "none",
                            reference_group = "Control", plots = FALSE)
  ))

  pw_log  <- as.data.frame(res_log$pairwise_comparisons)
  pw_none <- as.data.frame(res_none$pairwise_comparisons)
  est_col <- intersect(c("estimate", "Estimate"), colnames(pw_log))[1]
  expect_false(is.na(est_col))
  est_log  <- as.numeric(pw_log[[est_col]])
  est_none <- as.numeric(pw_none[[est_col]])
  expect_false(isTRUE(all.equal(est_log, est_none, tolerance = 1e-3)),
               info = "log vs none should produce materially different pairwise estimates")
})

test_that("tumor_growth_statistics: auc_bootstrap_n actually populates boot CI columns", {
  df <- make_tg_simple()
  res_no  <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, model_type = "auc",
                            reference_group = "Control",
                            auc_bootstrap_n = 0L,
                            plots = FALSE)
  ))
  res_yes <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, model_type = "auc",
                            reference_group = "Control",
                            auc_bootstrap_n = 199L,
                            auc_bootstrap_seed = 1L,
                            plots = FALSE)
  ))

  pw_no  <- as.data.frame(res_no$pairwise_comparisons)
  pw_yes <- as.data.frame(res_yes$pairwise_comparisons)
  expect_true("boot_ci_lower" %in% colnames(pw_no))
  expect_true("boot_ci_lower" %in% colnames(pw_yes))
  expect_true(all(is.na(pw_no$boot_ci_lower)),
              info = "auc_bootstrap_n = 0 must leave boot_ci_lower NA")
  expect_true(any(!is.na(pw_yes$boot_ci_lower)),
              info = "auc_bootstrap_n > 0 must populate boot_ci_lower")
})

test_that("survival_statistics: ph_test is returned on the standard-Cox path", {
  set.seed(202606)
  n <- 24L
  df <- data.frame(
    Day              = pmax(round(stats::rexp(n, rate = 1 / 14)), 1L),
    Survival_Censor  = stats::rbinom(n, 1L, prob = 0.7),
    Treatment        = rep(c("Control", "Treatment"), each = n / 2L),
    ID               = paste0("M", seq_len(n)),
    Cage             = rep(c("C1", "C2", "C3", "C4"), length.out = n)
  )
  res <- suppressMessages(survival_statistics(df, verbose = FALSE))
  # ph_test only present when method_used == "cox"; if Firth/logrank
  # fallbacks fire (small fixture instability) the field is NULL.
  if (identical(res$method_used, "cox")) {
    expect_false(is.null(res$ph_test),
                 info = "survival_statistics cox path must return ph_test (CODE_REVIEW.md A.1)")
    expect_true(inherits(res$ph_test, "cox.zph"))
  } else {
    # CODE_REVIEW.md K.2 -- this used to skip(). A skip here is indistinguishable
    # from success, so a regression that stopped the Cox path from running (e.g.
    # by mis-detecting separation, as R3.27 found) would have gone unnoticed.
    # Assert the contract instead: whichever path ran, ph_test must be present
    # exactly when the Cox path was used.
    expect_true(res$method_used %in% c("coxphf", "logrank"),
                info = "non-Cox path must identify itself")
    expect_null(res$ph_test)
  }
})

test_that("analyze_body_weight: cage_column actually enters the model formula", {
  df <- make_bw_simple()
  res_no_cage <- suppressWarnings(suppressMessages(
    analyze_body_weight(df, weight_column = "Weight",
                        time_column = "Day",
                        treatment_column = "Treatment",
                        id_column = "ID",
                        adjust_tumor_weight = FALSE,
                        covariates = character(0))
  ))
  res_cage <- suppressWarnings(suppressMessages(
    analyze_body_weight(df, weight_column = "Weight",
                        time_column = "Day",
                        treatment_column = "Treatment",
                        id_column = "ID",
                        cage_column = "Cage",
                        adjust_tumor_weight = FALSE,
                        covariates = character(0))
  ))
  expect_false(isTRUE(res_no_cage$model_info$cage_in_model))
  expect_true(isTRUE(res_cage$model_info$cage_in_model),
              info = "analyze_body_weight must include Cage random effect when cage_column supplied (J.11)")

  # And the variance-component structure must differ.
  vc_no <- tryCatch(as.data.frame(lme4::VarCorr(res_no_cage$model)), error = function(e) NULL)
  vc_yes <- tryCatch(as.data.frame(lme4::VarCorr(res_cage$model)),   error = function(e) NULL)
  if (!is.null(vc_no) && !is.null(vc_yes)) {
    expect_false(identical(vc_no$grp, vc_yes$grp),
                 info = "VarCorr grouping must differ between with/without cage_column")
  }
})

test_that("therapeutic_window_metric: TGI < 0 is clamped to 0 (J.13)", {
  # Synthetic dataset with treatment that *accelerates* growth → negative TGI.
  set.seed(202607)
  ndays <- 5
  ids_per_group <- 6
  groups <- c("Control", "BadDrug")
  df <- do.call(rbind, lapply(groups, function(g) {
    rate <- if (g == "Control") 0.06 else 0.10  # BadDrug grows FASTER
    do.call(rbind, lapply(seq_len(ids_per_group), function(i) {
      data.frame(
        ID        = paste0(substr(g, 1, 1), i),
        Treatment = g,
        Day       = seq(0, 20, length.out = ndays),
        Volume    = round(100 * exp(rate * seq(0, 20, length.out = ndays))),
        Weight    = round(rnorm(ndays, 20, 0.5), 2),
        stringsAsFactors = FALSE
      )
    }))
  }))
  res <- suppressWarnings(suppressMessages(
    therapeutic_window_metric(df, reference_group = "Control")
  ))
  twm_tab <- res$twm_table
  bad_row <- twm_tab[twm_tab$Treatment == "BadDrug", , drop = FALSE]
  expect_equal(nrow(bad_row), 1L)
  expect_true(bad_row$TGI < 0,
              info = "BadDrug should have negative TGI on this fixture")
  expect_equal(bad_row$TWM, 0,
               info = "TWM must clamp TGI < 0 to 0, not return abs(TGI) (CODE_REVIEW.md J.13)")
})
