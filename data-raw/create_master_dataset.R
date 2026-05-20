# =============================================================================
# create_master_dataset.R
# Master synthetic dataset — all mouseExperiment analysis scenarios
# =============================================================================
#
# Study design
# ─────────────
#   6 treatment groups · 8 mice per group · 2 cages per group (4 mice/cage)
#   48 mice total · days 0–26 (every 2 days, up to 14 observations per mouse)
#   Mice reaching the volume endpoint (2500 mm³) are removed from study;
#   their last observation is marked Survival_Censor = 1.
#
# Treatment groups and growth parameters
# ────────────────────────────────────────
#   Group                       Dose  Growth rate  Expected event rate
#   Vehicle                        0   +0.220       ~100 % (all mice)
#   Drug_A Low                     5   +0.140       ~66 % (borderline)
#   Drug_A Mid                    15   +0.090       ~2 %  (nearly none)
#   Drug_A High                   30   +0.030        0 %  (no events)
#   Drug_B                        10   +0.155       ~88 % (most mice)
#   Drug_A Mid + Drug_B           15   -0.040        0 %  (tumor regression)
#
#   Survival ordering: Vehicle > Drug_B > Drug_A Low > Drug_A Mid
#                       >= Drug_A High = Combo
#
# Built-in statistical properties (set.seed(2026))
# ──────────────────────────────────────────────────
#   Tumor growth:   All treated groups grow slower than Vehicle (correct dir.)
#                   Drug_A High and Combo clearly significant vs Vehicle
#   Survival:       Vehicle worst, treated groups better (ordered by gr)
#   Body weight:    Drug_A High and Combo show significant weight loss
#                   Drug_B: minimal non-significant loss
#   Synergy:        Use control="Vehicle", drug_a="Drug_A Mid",
#                   drug_b="Drug_B", combo="Drug_A Mid + Drug_B".
#                   At day 10 (all mice alive):
#                     Drug_A Mid TGI ~73 %, Drug_B TGI ~48 %
#                     Bliss expected ~86 %, Combo actual ~93 % → positive synergy
#   Dose-response:  Filter to Drug_A series (Doses 0/5/15/30); volumes form
#                   a monotonically decreasing Hill-shaped dose-response curve
#   Necrosis:       Drug_A High (~25 %) and Combo (~50 %) develop persistent
#                   necrosis from day 12 onward (biologically: tumor stasis
#                   and regression are associated with central necrosis)
#
# Columns
# ────────
#   Mouse_ID, ID, Day, Length, Width, Height, Volume, Weight,
#   Treatment, Dose, Cage, Necrotic, Survival_Censor
#
# Usage guide
# ────────────
#   All analyses:        master_synthetic_data as-is
#   Survival:            filter to last row per mouse first
#   Synergy:             control="Vehicle", drug_a="Drug_A Mid",
#                        drug_b="Drug_B", combo="Drug_A Mid + Drug_B",
#                        eval_time_point=10
#   Dose-response:       filter Treatment %in%
#                          c("Vehicle","Drug_A Low","Drug_A Mid","Drug_A High")
#   Bayesian TWM:        Drug_A High (TG) + Drug_B (BW); or Drug_A Mid + Drug_B
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

set.seed(2026)

days           <- seq(0L, 26L, by = 2L)
n_per_group    <- 8L
mice_per_cage  <- 4L
surv_threshold <- 2500
v0_mean        <- 80
v0_sd          <- 12
bw0_mean       <- 22.0
bw0_sd         <- 1.2
re_gr_sd       <- 0.020
obs_log_sd     <- 0.07
bw_re_sd       <- 0.50
bw_obs_sd      <- 0.25

groups <- data.frame(
  treatment    = c(
    "Vehicle",
    "Drug_A Low",   "Drug_A Mid",  "Drug_A High",
    "Drug_B",
    "Drug_A Mid + Drug_B"
  ),
  dose         = c(0L, 5L, 15L, 30L, 10L, 15L),
  growth_rate  = c(0.220, 0.140, 0.090, 0.030, 0.155, -0.040),
  bw_slope     = c(0.040, -0.010, -0.050, -0.100, -0.030, -0.150),
  necrosis_day = c(NA_real_, NA_real_, NA_real_, 12, NA_real_, 12),
  necrosis_max = c(0, 0, 0, 0.25, 0, 0.50),
  stringsAsFactors = FALSE
)

all_mouse_rows <- vector("list", nrow(groups) * n_per_group)
mouse_counter  <- 1L

for (g in seq_len(nrow(groups))) {

  grp       <- groups[g, ]
  cage_base <- (g - 1L) * 2L + 1L

  for (m in seq_len(n_per_group)) {

    cage_id  <- cage_base + ((m - 1L) %/% mice_per_cage)
    mouse_id <- sprintf("M%03d", mouse_counter)

    re_gr <- rnorm(1L, 0, re_gr_sd)
    v0    <- max(20, rnorm(1L, v0_mean, v0_sd))
    bw0   <- rnorm(1L, bw0_mean, bw0_sd)
    bw_re <- rnorm(1L, 0, bw_re_sd)

    # Mouse-level necrosis onset: once past onset day all observations
    # for that mouse are marked necrotic (persistent necrosis model).
    has_necrosis <- !is.na(grp$necrosis_day) &&
      (grp$necrosis_max > 0) &&
      (runif(1L) < grp$necrosis_max)
    necrosis_onset <- if (has_necrosis) {
      grp$necrosis_day + sample(0:6, 1L)
    } else {
      Inf
    }

    mouse_obs <- vector("list", length(days))
    alive     <- TRUE

    for (di in seq_along(days)) {
      day <- days[di]
      if (!alive) break

      # Log-normal tumor volume
      v <- v0 * exp(
        (grp$growth_rate + re_gr) * day + rnorm(1L, 0, obs_log_sd)
      )
      v <- max(1, v)

      # Caliper dimensions: ellipsoid V = pi*L*W^2/6, L/W = 1.4
      # Invert: W = (6V / (1.4*pi))^(1/3); add per-axis measurement noise.
      w_base <- (6 * v / (1.4 * pi))^(1 / 3)
      w      <- round(w_base * exp(rnorm(1L, 0, 0.030)), 2)
      l      <- round(1.4 * w_base * exp(rnorm(1L, 0, 0.025)), 2)
      h      <- round(0.85 * w_base * exp(rnorm(1L, 0, 0.030)), 2)
      v_meas <- round(pi * l * w^2 / 6, 1)

      bw <- round(
        max(14, bw0 + grp$bw_slope * day + bw_re + rnorm(1L, 0, bw_obs_sd)),
        2
      )

      necrotic <- as.integer(day >= necrosis_onset)
      censor   <- as.integer(v > surv_threshold)
      if (censor == 1L) alive <- FALSE

      mouse_obs[[di]] <- data.frame(
        Mouse_ID        = mouse_id,
        ID              = mouse_counter,
        Day             = day,
        Length          = l,
        Width           = w,
        Height          = h,
        Volume          = v_meas,
        Weight          = bw,
        Treatment       = grp$treatment,
        Dose            = grp$dose,
        Cage            = sprintf("C%02d", cage_id),
        Necrotic        = necrotic,
        Survival_Censor = censor,
        stringsAsFactors = FALSE
      )
    }

    all_mouse_rows[[mouse_counter]] <- bind_rows(mouse_obs)
    mouse_counter <- mouse_counter + 1L
  }
}

master_synthetic_data <- bind_rows(all_mouse_rows)

# ── Quick sanity summary ──────────────────────────────────────────────────────

cat("\n=== master_synthetic_data summary ===\n")
cat("Total rows:", nrow(master_synthetic_data), "\n")
cat("Total mice:", length(unique(master_synthetic_data$Mouse_ID)), "\n\n")

cat("Survival events by group:\n")
events <- master_synthetic_data |>
  group_by(Treatment) |>
  summarise(
    n_mice    = n_distinct(Mouse_ID),
    n_events  = sum(Survival_Censor),
    event_pct = round(100 * n_events / n_mice, 1),
    .groups   = "drop"
  ) |>
  arrange(desc(event_pct))
print(as.data.frame(events))

cat("\nNecrotic observations by group:\n")
nec <- master_synthetic_data |>
  group_by(Treatment) |>
  summarise(
    n_necrotic = sum(Necrotic),
    pct        = round(100 * n_necrotic / n(), 1),
    .groups    = "drop"
  ) |>
  filter(n_necrotic > 0)
print(as.data.frame(nec))

cat("\nMean final Volume by group (last observation per mouse):\n")
final_vol <- master_synthetic_data |>
  group_by(Mouse_ID, Treatment) |>
  slice_max(Day, n = 1L, with_ties = FALSE) |>
  group_by(Treatment) |>
  summarise(
    mean_vol = round(mean(Volume), 0),
    sd_vol   = round(sd(Volume), 0),
    .groups  = "drop"
  )
print(as.data.frame(final_vol))

cat("\nMean final body weight by group (last observation per mouse):\n")
final_bw <- master_synthetic_data |>
  group_by(Mouse_ID, Treatment) |>
  slice_max(Day, n = 1L, with_ties = FALSE) |>
  group_by(Treatment) |>
  summarise(
    mean_bw = round(mean(Weight), 2),
    sd_bw   = round(sd(Weight), 2),
    .groups = "drop"
  )
print(as.data.frame(final_bw))

# Verify synergy properties at day 10 (all mice alive at this point)
cat("\nSynergy check at day 10 (all mice present):\n")
day10 <- master_synthetic_data |>
  filter(Day == 10) |>
  group_by(Treatment) |>
  summarise(mean_v = round(mean(Volume), 1), .groups = "drop")
v_ctrl <- day10$mean_v[day10$Treatment == "Vehicle"]
day10  <- day10 |>
  mutate(tgi = round(1 - mean_v / v_ctrl, 3))
print(as.data.frame(day10))
tgi_a    <- day10$tgi[day10$Treatment == "Drug_A Mid"]
tgi_b    <- day10$tgi[day10$Treatment == "Drug_B"]
tgi_comb <- day10$tgi[day10$Treatment == "Drug_A Mid + Drug_B"]
bliss    <- round(1 - (1 - tgi_a) * (1 - tgi_b), 3)
cat(sprintf(
  "  Drug_A Mid TGI=%.1f%%  Drug_B TGI=%.1f%%  Combo TGI=%.1f%%\n",
  100 * tgi_a, 100 * tgi_b, 100 * tgi_comb
))
cat(sprintf(
  "  Bliss expected=%.1f%%  Excess=%.1f%%  Synergistic: %s\n",
  100 * bliss, 100 * (tgi_comb - bliss), tgi_comb > bliss
))

usethis::use_data(master_synthetic_data, overwrite = TRUE)
cat("\nSaved: data/master_synthetic_data.rda\n")
