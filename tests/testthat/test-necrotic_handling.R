# =============================================================================
# Tests for necrotic tumor handling
#
# Tests tgs_handle_necrosis() directly, then integration via
# tumor_growth_statistics() to confirm parameter flow and backwards
# compatibility.
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: small two-group dataset with a Necrotic column
# ---------------------------------------------------------------------------
make_necrotic_df <- function() {
  set.seed(42)
  days <- c(0, 7, 14, 21)
  make_mouse <- function(id, trt, slope, necrotic_days = integer(0)) {
    vol   <- exp(log(200) + slope * days + rnorm(4, 0, 0.05))
    nec   <- ifelse(days %in% necrotic_days, "Y", "N")
    data.frame(ID = id, Cage = "C1", Treatment = trt,
               Day = days, Volume = round(vol, 2),
               Necrotic = nec, stringsAsFactors = FALSE)
  }
  rbind(
    make_mouse("C01", "Control",    0.04),
    make_mouse("C02", "Control",    0.04),
    make_mouse("C03", "Control",    0.04, necrotic_days = 21L),   # 1 necrotic obs
    make_mouse("T01", "TreatmentA", 0.14),
    make_mouse("T02", "TreatmentA", 0.14, necrotic_days = c(14L, 21L)),
    make_mouse("T03", "TreatmentA", 0.14)
  )
}

# ---------------------------------------------------------------------------
# tgs_handle_necrosis() — unit tests
# ---------------------------------------------------------------------------
test_that("tgs_handle_necrosis: parses Y/N correctly (exclude mode)", {
  df  <- make_necrotic_df()
  res <- tgs_handle_necrosis(df, "Necrotic", "exclude",
                             "Treatment", "ID", "Day")
  # 3 necrotic rows should be excluded
  expect_equal(nrow(res$data), nrow(df) - 3L)
  # necrotic_cov_flag should only be 0 in the returned data
  expect_true(all(res$data$necrotic_cov_flag == 0L))
})

test_that("tgs_handle_necrosis: none mode keeps all rows", {
  df  <- make_necrotic_df()
  res <- tgs_handle_necrosis(df, "Necrotic", "none",
                             "Treatment", "ID", "Day")
  expect_equal(nrow(res$data), nrow(df))
})

test_that("tgs_handle_necrosis: covariate mode keeps all rows and adds flag", {
  df  <- make_necrotic_df()
  res <- tgs_handle_necrosis(df, "Necrotic", "covariate",
                             "Treatment", "ID", "Day")
  expect_equal(nrow(res$data), nrow(df))
  expect_true("necrotic_cov_flag" %in% colnames(res$data))
  expect_equal(sum(res$data$necrotic_cov_flag), 3L)  # 3 necrotic observations
})

test_that("tgs_handle_necrosis: necrosis_summary has correct structure", {
  df  <- make_necrotic_df()
  res <- tgs_handle_necrosis(df, "Necrotic", "exclude",
                             "Treatment", "ID", "Day")
  sm  <- res$necrosis_summary
  expect_s3_class(sm, "data.frame")
  required_cols <- c("Treatment", "N_Observations", "N_Necrotic",
                     "Pct_Necrotic", "N_Mice_Total", "N_Mice_With_Necrosis",
                     "First_Necrotic_Day")
  expect_true(all(required_cols %in% colnames(sm)),
              info = paste("Missing:", paste(setdiff(required_cols, colnames(sm)), collapse = ", ")))
  expect_equal(nrow(sm), 2L)  # two treatment groups
})

test_that("tgs_handle_necrosis: N_Necrotic counts are correct per group", {
  df  <- make_necrotic_df()
  res <- tgs_handle_necrosis(df, "Necrotic", "exclude",
                             "Treatment", "ID", "Day")
  sm  <- res$necrosis_summary
  ctrl_row <- sm[sm$Treatment == "Control", ]
  tx_row   <- sm[sm$Treatment == "TreatmentA", ]
  expect_equal(ctrl_row$N_Necrotic, 1L)
  expect_equal(tx_row$N_Necrotic,   2L)
})

test_that("tgs_handle_necrosis: parses alternative necrotic encodings", {
  df <- make_necrotic_df()
  # Test numeric 1/0
  df2 <- df
  df2$Necrotic <- ifelse(df$Necrotic == "Y", "1", "0")
  res <- tgs_handle_necrosis(df2, "Necrotic", "exclude",
                              "Treatment", "ID", "Day")
  expect_equal(nrow(res$data), nrow(df) - 3L)
  # Test TRUE/FALSE strings
  df3 <- df
  df3$Necrotic <- ifelse(df$Necrotic == "Y", "TRUE", "FALSE")
  res3 <- tgs_handle_necrosis(df3, "Necrotic", "exclude",
                               "Treatment", "ID", "Day")
  expect_equal(nrow(res3$data), nrow(df) - 3L)
})

# ---------------------------------------------------------------------------
# tumor_growth_statistics() — integration with necrotic_column
# ---------------------------------------------------------------------------
test_that("tumor_growth_statistics: necrosis_summary is NULL when necrotic_column=NULL", {
  df  <- make_necrotic_df()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      reference_group  = "Control",
      plots            = FALSE,
      verbose          = FALSE
    )
  ))
  expect_null(res$necrosis_summary)
})

test_that("tumor_growth_statistics: necrosis_summary is non-NULL when necrotic_column supplied", {
  df  <- make_necrotic_df()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column       = "Day",
      volume_column     = "Volume",
      treatment_column  = "Treatment",
      id_column         = "ID",
      reference_group   = "Control",
      necrotic_column   = "Necrotic",
      necrotic_handling = "exclude",
      plots             = FALSE,
      verbose           = FALSE
    )
  ))
  expect_false(is.null(res$necrosis_summary))
  expect_s3_class(res$necrosis_summary, "data.frame")
})

test_that("tumor_growth_statistics: exclude mode fits model on reduced dataset", {
  df  <- make_necrotic_df()
  res_full <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      reference_group  = "Control",
      plots            = FALSE, verbose = FALSE
    )
  ))
  res_excl <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column       = "Day",
      volume_column     = "Volume",
      treatment_column  = "Treatment",
      id_column         = "ID",
      reference_group   = "Control",
      necrotic_column   = "Necrotic",
      necrotic_handling = "exclude",
      plots             = FALSE, verbose = FALSE
    )
  ))
  # Excluding 3 rows should yield a slightly different model (AIC will differ)
  expect_false(identical(res_full$anova, res_excl$anova))
})

test_that("tumor_growth_statistics: covariate mode includes necrotic_cov_flag in model frame", {
  df  <- make_necrotic_df()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column       = "Day",
      volume_column     = "Volume",
      treatment_column  = "Treatment",
      id_column         = "ID",
      reference_group   = "Control",
      necrotic_column   = "Necrotic",
      necrotic_handling = "covariate",
      plots             = FALSE, verbose = FALSE
    )
  ))
  # necrotic_cov_flag should appear in the fitted model's model frame
  mf <- tryCatch(stats::model.frame(res$model), error = function(e) NULL)
  expect_true(!is.null(mf) && "necrotic_cov_flag" %in% colnames(mf))
})
