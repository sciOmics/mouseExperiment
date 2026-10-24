# =============================================================================
# helpers_format.R — CODE_REVIEW.md D.1 / D.9
#
# Presentation is not cosmetic when the number lives on a transformed scale:
# an axis must say which scale its values are on.
#
# This is the consumer side of R3.D6: that fix made `transform_used` report what
# the backend did rather than what was requested. These tests pin what the
# dashboard then does with it.
# =============================================================================

# ---- scale_suffix ----------------------------------------------------------

test_that("D.9: scale_suffix labels the three transforms the backend fits on", {
  # back_transform(), its converting partner, had no caller and was removed in
  # 0.42.0 (R22.11).
  for (tf in c("log", "sqrt", "none")) {
    expect_type(scale_suffix(tf), "character")
  }
  expect_identical(scale_suffix("log"),  " (log scale)")
  expect_identical(scale_suffix("sqrt"), " (sqrt scale)")
  expect_identical(scale_suffix("none"), "")
  # Unknown name: no label.
  expect_identical(scale_suffix("cuberoot"), "")
})

# ---- p-value formatting ----------------------------------------------------

test_that("D.9: tiny p-values are not rounded to zero", {
  # A p-value displayed as "0" is a claim no test can support. The whole reason
  # format_pval exists is to avoid it.
  expect_match(format_pval(1e-12), "e-12")
  expect_false(format_pval(1e-12) == "0")
  expect_match(format_pval(0), "^< ")
  expect_match(format_pval(1e-300), "^< ")
})

test_that("D.9: ordinary p-values format readably and NA survives", {
  expect_identical(format_pval(0.05), "0.05")
  expect_identical(format_pval(0.0499), "0.0499")
  expect_true(is.na(format_pval(NA_real_)))
})

test_that("D.9: format_pval is vectorised", {
  out <- format_pval(c(0.5, 1e-20, NA))
  expect_length(out, 3L)
  expect_true(is.na(out[3]))
  expect_match(out[2], "^< ")
})

# ---- table display ---------------------------------------------------------

test_that("D.9: p-value columns are detected by name and formatted as such", {
  # The pattern match is what decides whether a column is rounded to 4dp (and so
  # displayed as 0) or shown in scientific notation. Column naming across the
  # backend is inconsistent enough that this is worth pinning.
  df <- data.frame(Estimate = 1.23456789, p.value = 1e-14, p_adj = 1e-9,
                   padj = 1e-8, `Pr(>F)` = 1e-11, check.names = FALSE)
  out <- format_df_for_display(df)
  expect_identical(out$Estimate, round(1.23456789, 4))
  for (col in c("p.value", "p_adj", "padj", "Pr(>F)")) {
    expect_true(is.character(out[[col]]), info = col)
    expect_false(out[[col]] == "0", info = col)
  }
})

test_that("D.9: array- and table-class columns are flattened", {
  # survival_statistics() returns table-class columns; round() and datatable()
  # both mishandle them. This branch exists for that and had no test.
  df <- data.frame(Treatment = c("A", "B"), stringsAsFactors = FALSE)
  df$Events <- table(c("A", "B"))
  out <- format_df_for_display(df)
  expect_null(dim(out$Events))
  expect_false(inherits(out$Events, "table"))
})

test_that("D.9: non-numeric columns pass through untouched", {
  df <- data.frame(Contrast = "A - B", Note = NA_character_,
                   stringsAsFactors = FALSE)
  out <- format_df_for_display(df)
  expect_identical(out$Contrast, "A - B")
})

# ---- provenance ------------------------------------------------------------

test_that("D.9: comparisons_title reports the adjustment actually applied", {
  # R3.D1/R3.D2: the header used to assert "Dunnett-Adjusted" regardless of what
  # had been computed. It must now read the result.
  expect_match(comparisons_title(list(p_adjust_method_used = "tukey"), "All pairs"),
               "Tukey HSD", fixed = TRUE)
  expect_match(comparisons_title(list(p_adjust_method_used = "none"), "All pairs"),
               "unadjusted", fixed = TRUE)
  expect_match(comparisons_title(list(p_adjust_method_used = "fdr"), "All pairs"),
               "Benjamini-Hochberg", fixed = TRUE)
  # Dunnett must name the exact multivariate-t method, not the "dunnettx"
  # approximation the backend deliberately moved away from.
  expect_match(comparisons_title(list(p_adjust_method_used = "dunnett"), "vs control"),
               "multivariate-t", fixed = TRUE)
})

test_that("D.9: an unknown adjustment name is shown verbatim, not dropped", {
  # Better to display something unrecognised than to imply no adjustment.
  ttl <- comparisons_title(list(p_adjust_method_used = "sidak"), "All pairs")
  expect_match(ttl, "sidak", fixed = TRUE)
})

test_that("D.9: a result with no adjustment recorded gets no adjustment claim", {
  ttl <- comparisons_title(list(), "All pairs")
  expect_match(ttl, "All pairs", fixed = TRUE)
  expect_false(grepl("Tukey|Dunnett|unadjusted", ttl))
})

test_that("D.9: the adjustment scope is appended when the backend reports one", {
  res <- list(p_adjust_method_used = "holm",
              pairwise_comparisons = data.frame(
                contrast = "A - B", Adjust_Scope = "within day",
                stringsAsFactors = FALSE))
  expect_match(comparisons_title(res, "All pairs"), "within day", fixed = TRUE)
})

test_that("D.9: resolve_interval_cols honours a per-table override", {
  # Survival's tables use CI_Lower/CI_Upper while the headline result declares
  # Lower_CL/Upper_CL. The override branch exists for that and was untested.
  df <- data.frame(Group = "A", CI_Lower = 1, CI_Upper = 2)
  res <- list(meta = list(
    interval_type = "confidence",
    interval_columns = c(lower = "Lower_CL", upper = "Upper_CL"),
    interval_columns_override = c(lower = "CI_Lower", upper = "CI_Upper")))
  ic <- resolve_interval_cols(df, res)
  expect_identical(ic$lower, "CI_Lower")
  expect_identical(ic$upper, "CI_Upper")
})

test_that("D.9: a result declaring no intervals falls through to sniffing", {
  df <- data.frame(Group = "A", Lower_CL = 1, Upper_CL = 2)
  ic <- resolve_interval_cols(df, list(meta = list(interval_type = "none")))
  expect_identical(ic$lower, "Lower_CL")
})
