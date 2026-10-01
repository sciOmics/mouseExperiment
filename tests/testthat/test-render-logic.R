# =============================================================================
# Render-logic tests — CODE_REVIEW.md I.2
#
# I.2 has been the standing gap for several rounds. It was reported closed twice
# and is not: `test-module_servers.R` drives every module through testServer() but
# asserts only `inherits(session, "session_proxy")` and that `session$ns()` returns
# a non-empty string. That verifies testServer built a mock session. It does not
# touch the module.
#
# `test-backend-contract.R` closed the other half — the *field names* the renderers
# read. What neither covers is the render logic between them: the branches that
# decide whether a panel says "excludes zero" or "spans zero", whether a floor
# warning appears, whether an attrition figure is stated or suppressed. Round 3 and
# Round 8 added roughly a dozen such branches and all of them were verified by
# reading the source.
#
# These tests set rv$<X>_results to a real backend result and read the rendered
# output back. A branch that stops firing fails here.
# =============================================================================

skip_if_not_installed("shiny", minimum_version = "1.5.0")
skip_if_not_installed("mouseExperiment")

.pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperimentDashboard namespace not loadable")

# renderUI output comes back as a list carrying $html; renderPrint as a character
# scalar. Flatten either to plain text so assertions read the same way.
render_text <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}

rv_for <- function(...) {
  extra <- list(...)
  base <- list(
    uploaded_data = NULL, data_valid = TRUE,
    settings = list(col_id = "ID", col_treatment = "Treatment",
                    col_day = "Day", col_volume = "Volume", col_cage = "Cage")
  )
  do.call(shiny::reactiveValues, c(base, extra))
}

# ---- Power module: the dropout statement (K.15 / D.6) -----------------------

test_that("I.2: the Summary tab states the no-dropout assumption", {
  # The whole point of D.6 is that a complete-data power figure must not read as
  # unconditional. If this branch stops firing the number silently becomes a
  # claim it cannot support.
  fn <- get("power_server", envir = .pkg)
  res <- list(
    scenario_table = data.frame(N_Per_Group = c(5L, 10L), Power = c(0.45, 0.78),
                                Alpha = 0.05),
    power_curve_data = data.frame(N_Per_Group = c(5L, 10L), Power = c(0.45, 0.78),
                                  Alpha = 0.05),
    sensitivity_table = NULL, mode = "apriori_lmm", ap_submode = "find_power",
    lmm_params = list(dropout_limit = Inf),
    lmm_attrition = data.frame(N_Per_Group = c(5L, 10L), Alpha = 0.05,
                               Dropout_Limit = Inf, Pct_Obs_Lost = 0,
                               Pct_Mice_Truncated = 0)
  )
  shiny::testServer(fn, args = list(rv = rv_for(ap_results = res)), {
    txt <- render_text(output$pw_ap_summary_ui)
    expect_true(nzchar(txt))
    expect_match(txt, "none modelled")
    expect_match(txt, "upper bound")
  })
})

test_that("I.2: the Summary tab reports realised attrition when dropout is on", {
  fn <- get("power_server", envir = .pkg)
  res <- list(
    scenario_table = data.frame(N_Per_Group = 8L, Power = 0.71, Alpha = 0.05),
    power_curve_data = data.frame(N_Per_Group = 8L, Power = 0.71, Alpha = 0.05),
    sensitivity_table = NULL, mode = "apriori_lmm", ap_submode = "find_power",
    lmm_params = list(dropout_limit = 700),
    lmm_attrition = data.frame(N_Per_Group = 8L, Alpha = 0.05,
                               Dropout_Limit = 700, Pct_Obs_Lost = 30.2,
                               Pct_Mice_Truncated = 61.5)
  )
  shiny::testServer(fn, args = list(rv = rv_for(ap_results = res)), {
    txt <- render_text(output$pw_ap_summary_ui)
    expect_match(txt, "Dropout modelled")
    expect_match(txt, "30\\.2")
    expect_false(grepl("none modelled", txt, fixed = TRUE),
                 label = "must not claim complete data when dropout ran")
  })
})

test_that("I.2: an older result without attrition still renders", {
  # Backwards compatibility: a resumed session saved before K.15 has no
  # lmm_attrition. The branch must yield NULL, not error the whole panel out.
  fn <- get("power_server", envir = .pkg)
  res <- list(
    scenario_table = data.frame(N_Per_Group = 8L, Power = 0.71, Alpha = 0.05),
    power_curve_data = data.frame(N_Per_Group = 8L, Power = 0.71, Alpha = 0.05),
    sensitivity_table = NULL, mode = "apriori_lmm", ap_submode = "find_power"
  )
  shiny::testServer(fn, args = list(rv = rv_for(ap_results = res)), {
    txt <- render_text(output$pw_ap_summary_ui)
    expect_true(nzchar(txt))
    expect_match(txt, "Power-Analysis Summary")
  })
})

# ---- Tumor growth: comparison provenance and permutation panel --------------

test_that("I.2: the permutation panel distinguishes 'not run' from a result", {
  # tg_perm_test_ui has two branches and Round 3 added both. A permutation test
  # that silently reports nothing is worse than one that says it did not run.
  fn <- get("tumor_growth_server", envir = .pkg)

  not_run <- list(error = FALSE, model_type_used = "lme4",
                  permutation_test = NULL)
  shiny::testServer(fn, args = list(rv = rv_for(tg_results = not_run)), {
    txt <- render_text(output$tg_perm_test_ui)
    # Either an explicit "not run" statement or nothing at all; what it must not
    # do is render a p-value it does not have.
    expect_false(grepl("p +=", txt))
  })

  ran <- list(
    error = FALSE, model_type_used = "lme4",
    permutation_test = list(
      p_value = 0.002, unit = "mouse", n_units = 16L, n_perm_converged = 999L,
      exhaustive = FALSE, lrt_p_asymptotic = 1e-9, min_attainable_p = 0.001,
      n_distinct_statistics = 12870, method = "Randomisation test (mouse-level)")
  )
  shiny::testServer(fn, args = list(rv = rv_for(tg_results = ran)), {
    txt <- render_text(output$tg_perm_test_ui)
    expect_true(nzchar(txt))
    expect_match(txt, "mouse")
  })
})

# ---- Toxicity: the tolerability note (R20-K) --------------------------------

tox_res <- function(over_a) {
  # The toxicity module nests the therapeutic-window result under $twm. The
  # TWM ratio and its ranking note were replaced by two axes in v26.09.005.
  list(twm_reference = "Control", twm = list(
    error = FALSE,
    window_table = data.frame(
      Treatment = c("Control", "DrugA", "DrugB"), N_Animals = c(8L, 8L, 8L),
      TGI = c(0, 55, 40), TGI_Lower = c(0, 40, 22), TGI_Upper = c(0, 68, 55),
      Worst_Loss = c(3, 9, 24), Worst_Loss_Lower = c(2, 6, 21),
      Worst_Loss_Upper = c(4, 12, 27), N_Over_Threshold = c(0L, over_a, 6L),
      Tolerability = c("Tolerated", "Tolerated", "Not tolerated"),
      stringsAsFactors = FALSE),
    tolerability_threshold = 20,
    tolerability_rule = "Tolerated: the upper 95% bound ... below 20%.",
    interval_note = "Weight-loss intervals: bootstrap of animals within arm.",
    attrition = data.frame(Treatment = c("Control", "DrugA", "DrugB"),
                           N_Enrolled = c(8L, 8L, 8L), N_At_Endpoint = c(6L, 6L, 5L),
                           Pct_Lost = c(25, 25, 37.5)),
    endpoint_day = 28, endpoint_method = "model"))
}

test_that("R20-K: the tolerability note states the rule and flags a tolerated arm with animals over", {
  # A "Tolerated" arm describes its average animal; an animal over the
  # threshold in such an arm must be called out, not hidden behind the flag.
  fn <- get("toxicity_server", envir = .pkg)
  shiny::testServer(fn, args = list(rv = rv_for(tox_results = tox_res(0L))), {
    txt <- render_text(output$tox_twm_tolerability_note)
    expect_match(txt, "Tolerated: the upper 95% bound")
    expect_false(grepl("reached the threshold on their own", txt))
  })
  shiny::testServer(fn, args = list(rv = rv_for(tox_results = tox_res(2L))), {
    txt <- render_text(output$tox_twm_tolerability_note)
    expect_match(txt, "DrugA \\(2 of 8\\)")
    expect_match(txt, "reached the threshold on their own")
  })
})
