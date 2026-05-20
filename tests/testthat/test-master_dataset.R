# =============================================================================
# test-master_dataset.R
# Directional property tests for master_synthetic_data (set.seed(2026))
#
# Tests validate the known built-in properties of the dataset rather than
# exact numeric values, ensuring the dataset remains suitable for demonstrating
# and testing all package analysis functions.
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

data(master_synthetic_data, package = "mouseExperiment")
dat <- master_synthetic_data

# Convenience: last observation per mouse (used for survival / volume ordering)
last_obs <- dat |>
  group_by(Mouse_ID) |>
  slice_max(Day, n = 1L, with_ties = FALSE) |>
  ungroup()

# ── Structure ─────────────────────────────────────────────────────────────────

test_that("master_synthetic_data has expected structure", {
  expect_s3_class(dat, "data.frame")
  expect_true(nrow(dat) > 0L)
  expect_setequal(
    colnames(dat),
    c("Mouse_ID","ID","Day","Length","Width","Height",
      "Volume","Weight","Treatment","Dose","Cage",
      "Necrotic","Survival_Censor")
  )
})

test_that("48 unique mice, 6 treatment groups, 12 cages", {
  expect_equal(n_distinct(dat$Mouse_ID), 48L)
  expect_equal(n_distinct(dat$Treatment), 6L)
  expect_equal(n_distinct(dat$Cage), 12L)
})

test_that("study days are 0–26 in steps of 2", {
  expect_equal(sort(unique(dat$Day)), seq(0L, 26L, by = 2L))
})

test_that("each group has 8 mice and 2 cages with 4 mice each", {
  per_group <- dat |>
    group_by(Treatment) |>
    summarise(
      n_mice = n_distinct(Mouse_ID),
      n_cages = n_distinct(Cage),
      .groups = "drop"
    )
  expect_true(all(per_group$n_mice == 8L))
  expect_true(all(per_group$n_cages == 2L))
})

test_that("Volume, Length, Width, Height are all positive", {
  expect_true(all(dat$Volume > 0))
  expect_true(all(dat$Length > 0))
  expect_true(all(dat$Width  > 0))
  expect_true(all(dat$Height > 0))
})

test_that("Necrotic and Survival_Censor are 0/1 integers", {
  expect_true(all(dat$Necrotic %in% c(0L, 1L)))
  expect_true(all(dat$Survival_Censor %in% c(0L, 1L)))
})

test_that("Survival_Censor = 1 only on last observation per mouse", {
  censor_rows <- dat |> filter(Survival_Censor == 1L)
  # Each censored mouse appears exactly once in censor_rows
  expect_equal(
    nrow(censor_rows),
    n_distinct(censor_rows$Mouse_ID)
  )
  # All censored rows are the last (max Day) for that mouse
  check <- censor_rows |>
    left_join(
      dat |> group_by(Mouse_ID) |> summarise(max_day = max(Day), .groups = "drop"),
      by = "Mouse_ID"
    )
  expect_true(all(check$Day == check$max_day))
})

# ── Tumor growth ordering ─────────────────────────────────────────────────────

test_that("mean final tumor volume decreases with treatment efficacy", {
  mean_vol <- last_obs |>
    group_by(Treatment) |>
    summarise(mv = mean(Volume), .groups = "drop")

  get_mv <- function(tx) mean_vol$mv[mean_vol$Treatment == tx]

  # Vehicle > Drug_A Low  (treatment reduces growth)
  expect_gt(get_mv("Vehicle"), get_mv("Drug_A Low"))
  # Drug_A Low > Drug_A Mid (dose-dependent reduction)
  expect_gt(get_mv("Drug_A Low"), get_mv("Drug_A Mid"))
  # Drug_A Mid > Drug_A High (dose-dependent reduction)
  expect_gt(get_mv("Drug_A Mid"), get_mv("Drug_A High"))
  # Drug_A Mid + Drug_B is smallest (tumor regression)
  expect_lt(get_mv("Drug_A Mid + Drug_B"), get_mv("Drug_A High"))
})

test_that("TGI at day 10 follows expected ordering", {
  day10 <- dat |>
    filter(Day == 10L) |>
    group_by(Treatment) |>
    summarise(mean_v = mean(Volume), .groups = "drop")

  v_ctrl <- day10$mean_v[day10$Treatment == "Vehicle"]
  tgi <- function(tx) 1 - day10$mean_v[day10$Treatment == tx] / v_ctrl

  expect_gt(tgi("Drug_A High"),         tgi("Drug_A Mid"))
  expect_gt(tgi("Drug_A Mid"),          tgi("Drug_A Low"))
  expect_gt(tgi("Drug_A Low"),          0)
  expect_gt(tgi("Drug_A Mid + Drug_B"), tgi("Drug_A High"))
})

# ── Survival ordering ─────────────────────────────────────────────────────────

test_that("survival event rates follow expected ordering", {
  events <- dat |>
    group_by(Treatment) |>
    summarise(
      n_mice   = n_distinct(Mouse_ID),
      n_events = sum(Survival_Censor),
      rate     = n_events / n_mice,
      .groups  = "drop"
    )

  get_rate <- function(tx) events$rate[events$Treatment == tx]

  # Vehicle: all mice hit endpoint
  expect_equal(get_rate("Vehicle"), 1.0)
  # Drug_A High and Combo: no events
  expect_equal(get_rate("Drug_A High"),         0.0)
  expect_equal(get_rate("Drug_A Mid + Drug_B"), 0.0)
  # Treated groups have lower event rates than Vehicle
  expect_lt(get_rate("Drug_A Low"), get_rate("Vehicle"))
  expect_lt(get_rate("Drug_B"),     get_rate("Vehicle"))
  expect_lt(get_rate("Drug_A Mid"), get_rate("Vehicle"))
})

# ── Body weight ordering ──────────────────────────────────────────────────────

test_that("mean final body weight decreases with treatment toxicity", {
  mean_bw <- last_obs |>
    group_by(Treatment) |>
    summarise(mbw = mean(Weight), .groups = "drop")

  get_bw <- function(tx) mean_bw$mbw[mean_bw$Treatment == tx]

  # Vehicle has the highest weight (positive bw_slope)
  expect_gt(get_bw("Vehicle"),             get_bw("Drug_A High"))
  expect_gt(get_bw("Drug_A High"),         get_bw("Drug_A Mid + Drug_B"))
  # Drug_B minimal loss, less than Drug_A Mid
  expect_gt(get_bw("Drug_B"),              get_bw("Drug_A Mid"))
})

# ── Drug synergy at day 10 ────────────────────────────────────────────────────

test_that("combination 'Drug_A Mid + Drug_B' shows positive Bliss excess at day 10", {
  day10 <- dat |>
    filter(Day == 10L) |>
    group_by(Treatment) |>
    summarise(mean_v = mean(Volume), .groups = "drop")

  get_v <- function(tx) day10$mean_v[day10$Treatment == tx]
  v_ctrl <- get_v("Vehicle")

  tgi_a    <- 1 - get_v("Drug_A Mid") / v_ctrl
  tgi_b    <- 1 - get_v("Drug_B")     / v_ctrl
  tgi_comb <- 1 - get_v("Drug_A Mid + Drug_B") / v_ctrl
  bliss    <- 1 - (1 - tgi_a) * (1 - tgi_b)

  # Each monotherapy has positive TGI
  expect_gt(tgi_a,    0)
  expect_gt(tgi_b,    0)
  # Combo exceeds both monotherapies
  expect_gt(tgi_comb, tgi_a)
  expect_gt(tgi_comb, tgi_b)
  # Positive Bliss excess (synergy)
  expect_gt(tgi_comb, bliss)
})

# ── Dose-response (Drug_A series) ─────────────────────────────────────────────

test_that("Drug_A dose-response is monotonically decreasing at day 10", {
  day10 <- dat |>
    filter(
      Day == 10L,
      Treatment %in% c("Vehicle", "Drug_A Low", "Drug_A Mid", "Drug_A High")
    ) |>
    group_by(Treatment, Dose) |>
    summarise(mean_v = mean(Volume), .groups = "drop") |>
    arrange(Dose)

  vols <- day10$mean_v
  expect_true(all(diff(vols) < 0),
    info = "Mean volumes should decrease with each dose step")
})

test_that("dose column correctly encodes Drug_A series doses", {
  doses <- dat |>
    filter(Treatment %in% c("Vehicle","Drug_A Low","Drug_A Mid","Drug_A High")) |>
    distinct(Treatment, Dose) |>
    arrange(Dose)

  expect_equal(sort(unique(doses$Dose)), c(0L, 5L, 15L, 30L))
})

# ── Necrosis ─────────────────────────────────────────────────────────────────

test_that("necrosis occurs only in Drug_A High and Drug_A Mid + Drug_B", {
  necrotic_groups <- dat |>
    filter(Necrotic == 1L) |>
    distinct(Treatment) |>
    pull(Treatment)

  expect_setequal(
    necrotic_groups,
    c("Drug_A High", "Drug_A Mid + Drug_B")
  )
})

test_that("necrosis is persistent: once onset, all subsequent observations are necrotic", {
  has_necrosis <- dat |>
    group_by(Mouse_ID) |>
    filter(any(Necrotic == 1L)) |>
    ungroup()

  # For each affected mouse: no non-necrotic observation after the first necrotic one
  all_persistent <- has_necrosis |>
    group_by(Mouse_ID) |>
    summarise(
      first_nec_day = min(Day[Necrotic == 1L]),
      has_non_nec_after = any(Necrotic == 0L & Day >= first_nec_day[1]),
      .groups = "drop"
    ) |>
    pull(has_non_nec_after)

  expect_true(!any(all_persistent))
})

test_that("necrotic observations start from day 12 onward (onset >= 12)", {
  first_nec_days <- dat |>
    filter(Necrotic == 1L) |>
    group_by(Mouse_ID) |>
    summarise(first_nec = min(Day), .groups = "drop") |>
    pull(first_nec)

  expect_true(all(first_nec_days >= 12L))
})

# ── Caliper dimensions consistency ───────────────────────────────────────────

test_that("Length >= Width for all observations", {
  expect_true(all(dat$Length >= dat$Width))
})

test_that("Volume is approximately consistent with caliper formula pi*L*W^2/6", {
  dat_check <- dat |>
    mutate(vol_calc = pi * Length * Width^2 / 6)

  # Measured Volume should be within 5% of the caliper-derived formula
  # (small discrepancies arise from rounding caliper values to 2 dp)
  rel_diff <- abs(dat_check$Volume - dat_check$vol_calc) / dat_check$vol_calc
  expect_true(mean(rel_diff) < 0.05)
})
