# =============================================================================
# Tests for calculate_volume()
# Verifies the formula implementations against hand-computed ground truth.
# =============================================================================

test_that("ellipsoid formula: V = (L * W^2 * pi) / 6 when L > W", {
  # Hand computation: L=10, W=5  →  V = (10 * 25 * pi) / 6 ≈ 130.900
  df  <- data.frame(Length = 10, Width = 5)
  out <- calculate_volume(df)
  expected <- (10 * 5^2 * pi) / 6

  expect_true("Volume" %in% colnames(out))
  expect_equal(out$Volume, expected, tolerance = 1e-6)
})

test_that("ellipsoid formula: L and W are swapped when W > L", {
  # Function ensures longer dimension = L, shorter = W regardless of column labels.
  # Length=5, Width=10  →  internally L=10, W=5  →  same V as above.
  df_normal  <- data.frame(Length = 10, Width = 5)
  df_swapped <- data.frame(Length = 5,  Width = 10)

  expect_equal(
    calculate_volume(df_normal)$Volume,
    calculate_volume(df_swapped)$Volume,
    tolerance = 1e-10
  )
})

test_that("modified_ellipsoid formula: V = (L * W^2) / 2", {
  # L=10, W=5  →  V = (10 * 25) / 2 = 125
  df  <- data.frame(Length = 10, Width = 5)
  out <- calculate_volume(df, formula = "modified_ellipsoid")

  expect_equal(out$Volume, 125, tolerance = 1e-10)
})

test_that("sphere formula: V = (pi * W^3) / 6", {
  # W=6  →  V = (pi * 216) / 6 = 36*pi ≈ 113.097
  df  <- data.frame(Length = 6, Width = 6)
  out <- calculate_volume(df, formula = "sphere")
  expected <- (pi * 6^3) / 6

  expect_equal(out$Volume, expected, tolerance = 1e-6)
})

test_that("cylinder formula: V = (pi * W^2 * L) / 4", {
  # L=10, W=4  →  V = (pi * 16 * 10) / 4 = 40*pi ≈ 125.664
  df  <- data.frame(Length = 10, Width = 4)
  out <- calculate_volume(df, formula = "cylinder")
  expected <- (pi * 4^2 * 10) / 4

  expect_equal(out$Volume, expected, tolerance = 1e-6)
})

test_that("box formula requires height column; falls back to width when absent", {
  # When height_column not supplied, height = width.
  # L=10, W=5, H≈W=5  →  V = 10*5*5 = 250
  df  <- data.frame(Length = 10, Width = 5)
  out <- calculate_volume(df, formula = "box")

  expect_equal(out$Volume, 10 * 5 * 5, tolerance = 1e-10)
})

test_that("box formula with explicit height column", {
  df  <- data.frame(Length = 10, Width = 5, Height = 3)
  out <- calculate_volume(df, formula = "box", height_column = "Height")

  expect_equal(out$Volume, 10 * 5 * 3, tolerance = 1e-10)
})

test_that("Volume column is appended; original columns are preserved", {
  df  <- data.frame(ID = "M1", Length = 8, Width = 4)
  out <- calculate_volume(df)

  expect_true(all(c("ID", "Length", "Width", "Volume") %in% colnames(out)))
  expect_equal(nrow(out), 1L)
})

test_that("vectorised over multiple rows", {
  df <- data.frame(
    Length = c(10, 8, 6),
    Width  = c(5,  4, 3)
  )
  out      <- calculate_volume(df)
  expected <- (c(10, 8, 6) * c(5, 4, 3)^2 * pi) / 6

  expect_equal(out$Volume, expected, tolerance = 1e-6)
})

test_that("missing Length column raises an error", {
  df <- data.frame(Width = 5)
  expect_error(calculate_volume(df), "Length")
})

test_that("missing Width column raises an error", {
  df <- data.frame(Length = 10)
  expect_error(calculate_volume(df), "Width")
})
