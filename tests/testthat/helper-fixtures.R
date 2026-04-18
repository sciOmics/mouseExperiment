# =============================================================================
# Test fixtures for mouseExperiment package
#
# All generators use set.seed() internally so that test runs produce identical
# data regardless of the calling order.  Helpers are automatically sourced by
# testthat before any test-*.R file is run.
# =============================================================================

# -----------------------------------------------------------------------------
# Tumor-growth: two-group linear exponential growth
#
# Ground-truth design:
#   - 2 groups × 6 mice × 4 time-points (days 0, 7, 14, 21)
#   - log-scale slopes: Control = 0.04/day, TreatmentA = 0.14/day
#   - Starting volume ≈ 200 mm³ (log-normal with SD_log = 0.05)
#   - SNR ≈ (0.14-0.04)*21 / 0.05 = 42 → LME interaction p << 0.001
# -----------------------------------------------------------------------------
make_tg_simple <- function() {
  set.seed(42)
  days <- c(0, 7, 14, 21)

  make_mouse <- function(id, cage, treatment, log_slope) {
    noise <- rnorm(4, mean = 0, sd = 0.05)
    vol   <- exp(log(200) + log_slope * days + noise)
    data.frame(
      ID        = id,
      Cage      = cage,
      Treatment = treatment,
      Day       = days,
      Volume    = round(vol, 2),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    make_mouse("C01", "C1", "Control",    0.04),
    make_mouse("C02", "C1", "Control",    0.04),
    make_mouse("C03", "C2", "Control",    0.04),
    make_mouse("C04", "C2", "Control",    0.04),
    make_mouse("C05", "C3", "Control",    0.04),
    make_mouse("C06", "C3", "Control",    0.04),
    make_mouse("T01", "T1", "TreatmentA", 0.14),
    make_mouse("T02", "T1", "TreatmentA", 0.14),
    make_mouse("T03", "T2", "TreatmentA", 0.14),
    make_mouse("T04", "T2", "TreatmentA", 0.14),
    make_mouse("T05", "T3", "TreatmentA", 0.14),
    make_mouse("T06", "T3", "TreatmentA", 0.14)
  )
}

# -----------------------------------------------------------------------------
# Survival: two groups with completely separate outcomes
#
# Ground-truth design:
#   - AggressiveTx: 7 mice that all die (Censor=1) at days 10, 15, 20, 25, 30, 35, 40
#     KM estimator: S(20)=4/7≈0.571, S(25)=3/7≈0.429
#     → KM median = 25 (first t where S(t) < 0.5, unambiguous – no tie at 0.5)
#   - Control: 7 mice all censored (Censor=0) at day 60
#     KM median = NA (never reaches 50%)
#   - Log-rank p << 0.0001
# -----------------------------------------------------------------------------
make_surv_simple <- function() {
  rbind(
    data.frame(
      ID               = paste0("A0", 1:7),
      Cage             = "CageA",
      Treatment        = "AggressiveTx",
      Day              = c(10, 15, 20, 25, 30, 35, 40),
      Survival_Censor  = 1L,
      stringsAsFactors = FALSE
    ),
    data.frame(
      ID               = paste0("C0", 1:7),
      Cage             = "CageC",
      Treatment        = "Control",
      Day              = 60L,
      Survival_Censor  = 0L,
      stringsAsFactors = FALSE
    )
  )
}

# -----------------------------------------------------------------------------
# AUC: exact trapezoid ground truth
#
# Two mice with perfectly controlled growth curves:
#   LowGrowth  (M1): Volume = 100, 200, 300, 400 at days 0, 7, 14, 21
#     → AUC = 7*(100+200)/2 + 7*(200+300)/2 + 7*(300+400)/2 = 5 250
#   HighGrowth (M2): Volume = 200, 400, 600, 800 at days 0, 7, 14, 21
#     → AUC = 7*(200+400)/2 + 7*(400+600)/2 + 7*(600+800)/2 = 10 500
# Expected AUC ratio: HighGrowth / LowGrowth = 2.0 exactly
# -----------------------------------------------------------------------------
make_auc_exact <- function() {
  data.frame(
    ID        = c(rep("M1", 4), rep("M2", 4)),
    Cage      = c(rep("C1", 4), rep("C2", 4)),
    Treatment = c(rep("LowGrowth", 4), rep("HighGrowth", 4)),
    Day       = rep(c(0L, 7L, 14L, 21L), 2),
    Volume    = c(100, 200, 300, 400,
                  200, 400, 600, 800),
    stringsAsFactors = FALSE
  )
}

# Exact expected individual AUC values (used in tests)
AUC_LOW_GROWTH  <- 5250   # 7*(150 + 250 + 350)
AUC_HIGH_GROWTH <- 10500  # 7*(300 + 500 + 700)

# -----------------------------------------------------------------------------
# Dose-response: clear monotonic inhibitory dose-response
#
# Design: 1 treatment (DrugX) + Control, 4 dose levels, 4 mice/group
#   dose = 0   → mean volume ~500  (Control)
#   dose = 1   → mean volume ~380
#   dose = 5   → mean volume ~220
#   dose = 25  → mean volume ~90
# Jonckheere trend test and linear model both expected significant (p < 0.01)
# -----------------------------------------------------------------------------
make_dose_response <- function() {
  set.seed(42)
  doses <- c(0, 1, 5, 25)
  means <- c(500, 380, 220, 90)
  n_per_dose <- 4

  do.call(rbind, lapply(seq_along(doses), function(i) {
    vols <- round(rnorm(n_per_dose, means[i], means[i] * 0.08), 1)
    data.frame(
      ID               = paste0("D", doses[i], "_M", seq_len(n_per_dose)),
      Cage             = paste0("Cage", i),
      Treatment        = ifelse(doses[i] == 0, "Control", "DrugX"),
      Dose             = doses[i],
      Day              = 21L,
      Volume           = pmax(vols, 5),
      stringsAsFactors = FALSE
    )
  }))
}

# -----------------------------------------------------------------------------
# Drug synergy: designed Bliss-independence data
#
# Design (all values at single time point, day 21):
#   Control   mean volume = 500  TGI = 0
#   DrugA     mean volume = 250  TGI = 0.50
#   DrugB     mean volume = 300  TGI = 0.40
#   Bliss expected FE = 0.50 + 0.40 - 0.20 = 0.70  → combo vol = 150
#   Loewe expected FE = min(0.50 + 0.40, 1.0) = 0.90  → combo vol = 50
#   Loewe CI = (FE_A + FE_B) / FE_combo = 0.90 / FE_combo
#
#   Bliss-neutral combo:  actual mean ≈ 150, FE=0.70, Loewe CI=1.29
#   Synergistic combo:    actual mean ≈ 25,  FE=0.95, Loewe CI=0.95 (CI < 1)
#   Antagonistic combo:   actual mean ≈ 325, FE=0.35, Loewe CI=2.57
# Uses 6 mice per group; 8% CV so group means are very close to target.
# -----------------------------------------------------------------------------
make_synergy_base <- function(combo_mean) {
  set.seed(42)
  n <- 6
  cv <- 0.08

  make_group <- function(name, mu) {
    data.frame(
      ID               = paste0(gsub("[^A-Za-z0-9]", "", name), "_M", seq_len(n)),
      Cage             = paste0(gsub("[^A-Za-z0-9]", "", name), "_C"),
      Treatment        = name,
      Day              = 21L,
      Volume           = round(rnorm(n, mu, mu * cv), 1),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    make_group("Control",    500),
    make_group("DrugA",      250),
    make_group("DrugB",      300),
    make_group("DrugA+DrugB", combo_mean)
  )
}

make_synergy_additive    <- function() make_synergy_base(150)   # Bliss-neutral
make_synergy_synergistic <- function() make_synergy_base(25)    # FE > Loewe expected → CI < 1
make_synergy_antagonist  <- function() make_synergy_base(325)   # Worse than additive

# Bliss ground-truth expected combo volume (used in assertions)
BLISS_EXPECTED_COMBO_VOL <- 150

# -----------------------------------------------------------------------------
# Four-group combo data (Control / HDACi / aPD1 / HDACi+aPD1)
# Mirrors the demo CSV structure; used for end-to-end pipeline tests.
# -----------------------------------------------------------------------------
make_combo_four_group <- function() {
  set.seed(7)
  groups <- c("Control", "HDACi", "aPD1", "HDACi+aPD1")
  slopes <- c(0.12, 0.07, 0.08, 0.03)   # log-scale growth rates
  n_per  <- 5
  days   <- c(0, 7, 14, 21, 28)

  do.call(rbind, lapply(seq_along(groups), function(g) {
    do.call(rbind, lapply(seq_len(n_per), function(m) {
      prefix <- substr(gsub("[^A-Za-z0-9]", "", groups[g]), 1, 3)
      cage   <- paste0(prefix, ceiling(m / 2))
      id     <- paste0(prefix, sprintf("%02d", m))
      noise  <- rnorm(length(days), 0, 0.06)
      vol    <- exp(log(200) + slopes[g] * days + noise)
      data.frame(
        ID               = id,
        Cage             = cage,
        Treatment        = groups[g],
        Day              = days,
        Volume           = round(vol, 2),
        Survival_Censor  = ifelse(days == max(days), 0L, NA_integer_),
        stringsAsFactors = FALSE
      )
    }))
  }))
}
