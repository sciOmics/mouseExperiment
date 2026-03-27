#!/usr/bin/env Rscript
# validate_all.R  ── comprehensive validation of dashboard model results
# against independent manual R calculations
# Run from mouseExperiment project directory

suppressWarnings(suppressMessages({
  library(mouseExperiment)
  library(lme4)
  library(emmeans)
  library(car)
  library(survival)
  library(dplyr)
}))

DATA <- "/home/david/Projects/mouseExperimentDashboard/inst/sample_data"
sep <- paste0(strrep("─", 70), "\n")

results <- list()   # collect for summary

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 1 – weight_synthetic_data (simple, no cages, pre-computed volume)
# ══════════════════════════════════════════════════════════════════════════════
cat(sep, "DATASET 1: weight_synthetic_data\n", sep)

d_w <- read.csv(file.path(DATA, "weight_synthetic_data.csv"))
d_w <- mouseExperiment::calculate_dates(d_w,
  date_column  = "Date",
  date_format  = "month_day",
  start_date   = "March 1"
)
cat("Days:", sort(unique(d_w$Day)), "\n")
cat("Column names after calculate_dates:", paste(colnames(d_w), collapse=", "), "\n")

# Rename to standard names expected by tumor_growth_statistics
d_w$Volume <- d_w$Tumor.Volume
d_w$ID     <- d_w$Ear.Tag
d_w$Cage   <- "1"   # no cage structure in this dataset

# ── Package: LME4 ─────────────────────────────────────────────────────────────
cat("\n--- Package LME4 ---\n")
pkg_lme4_w <- tryCatch(
  tumor_growth_statistics(
    df               = d_w,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    transform        = "log",
    model_type       = "lme4",
    reference_group  = "DMSO",
    plots            = FALSE,
    verbose          = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_lme4_w)) {
  cat("ANOVA (Type III):\n")
  print(pkg_lme4_w$anova)
  cat("\nPairwise comparisons (from pairwise_comparisons slot):\n")
  if (!is.null(pkg_lme4_w$pairwise_comparisons)) {
    pc <- as.data.frame(pkg_lme4_w$pairwise_comparisons)
    print(pc[, intersect(c("contrast","estimate","SE","df","t.ratio","p.value"), colnames(pc))])
  } else {
    cat("  NULL\n")
  }
  cat("\nTreatment effects (marginal means):\n")
  print(as.data.frame(pkg_lme4_w$treatment_effects))
}

# ── Manual: LME4 ──────────────────────────────────────────────────────────────
cat("\n--- Manual LME4 ---\n")
d_w$logVol <- log(d_w$Volume + 1e-6)   # guard zero
d_w$Treatment <- factor(d_w$Treatment)
d_w$Treatment <- relevel(d_w$Treatment, ref = "DMSO")

manual_lme4_w <- tryCatch(
  lme4::lmer(logVol ~ Treatment * Day + (1 | ID), data = d_w, REML = FALSE),
  error = function(e) { cat("lmer ERROR:", e$message, "\n"); NULL }
)

if (!is.null(manual_lme4_w)) {
  cat("Type III ANOVA (car::Anova):\n")
  print(car::Anova(manual_lme4_w, type = 3))
  
  emm_w <- emmeans::emmeans(manual_lme4_w, ~ Treatment)
  cat("\nEmmeans (averaged over Day):\n")
  print(summary(emm_w))
  
  cat("\nAll pairwise (Tukey):\n")
  pw_manual <- summary(pairs(emm_w, adjust = "tukey"), infer = c(TRUE, TRUE))
  print(as.data.frame(pw_manual)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
  
  cat("\nvs DMSO (Dunnett):\n")
  dun_manual <- summary(emmeans::contrast(emm_w, method = "trt.vs.ctrl",
                                          ref = which(levels(emm_w)[[1]] == "DMSO"),
                                          adjust = "dunnett"),
                        infer = c(TRUE, TRUE))
  print(as.data.frame(dun_manual)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
}

results$w_lme4_pkg_anova   <- if (!is.null(pkg_lme4_w)) as.data.frame(pkg_lme4_w$anova) else NULL
results$w_lme4_manual_anova <- if (!is.null(manual_lme4_w)) as.data.frame(car::Anova(manual_lme4_w, type=3)) else NULL

# ── Package: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Package AUC ---\n")
pkg_auc_w <- tryCatch(
  tumor_growth_statistics(
    df               = d_w,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    model_type       = "auc",
    reference_group  = "DMSO",
    plots            = FALSE,
    verbose          = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_auc_w)) {
  cat("AUC ANOVA:\n")
  print(pkg_auc_w$anova)
  cat("\nAUC treatment effects (Mean_AUC per group):\n")
  print(as.data.frame(pkg_auc_w$treatment_effects))
  cat("\nAUC pairwise (posthoc$pairwise):\n")
  if (!is.null(pkg_auc_w$posthoc$pairwise)) {
    print(pkg_auc_w$posthoc$pairwise)
  } else {
    cat("  NULL\n")
  }
}

# ── Manual: AUC ───────────────────────────────────────────────────────────────
cat("\n--- Manual AUC ---\n")
# Trapezoidal AUC per animal
manual_auc_w <- d_w %>%
  arrange(ID, Day) %>%
  group_by(ID, Treatment) %>%
  summarize(
    AUC = {
      tt <- Day; vv <- Volume
      # trapezoidal: sum of (dt * avg_v)
      n <- length(tt)
      if (n < 2) 0 else sum(diff(tt) * (vv[-n] + vv[-1]) / 2)
    },
    .groups = "drop"
  )

cat("Per-animal AUC values:\n")
print(as.data.frame(manual_auc_w))

cat("\nGroup mean AUC:\n")
manual_auc_grp <- manual_auc_w %>%
  group_by(Treatment) %>%
  summarize(Mean_AUC = mean(AUC), SD = sd(AUC), N = n(), .groups = "drop")
print(as.data.frame(manual_auc_grp))

cat("\nOne-way ANOVA on AUC:\n")
manual_aov_w <- aov(AUC ~ Treatment, data = manual_auc_w)
print(summary(manual_aov_w))

cat("\nPairwise Welch t-tests (Bonferroni-adjusted):\n")
pwt <- pairwise.t.test(manual_auc_w$AUC, manual_auc_w$Treatment,
                       pool.sd = FALSE, p.adjust.method = "bonferroni")
print(pwt)

results$w_auc_pkg  <- pkg_auc_w
results$w_auc_manual_grp <- manual_auc_grp

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 2 – combo_treatment_synthetic_data
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n", sep, "DATASET 2: combo_treatment_synthetic_data\n", sep)

d_c <- read.csv(file.path(DATA, "combo_treatment_synthetic_data.csv"))
cat("Treatments:", paste(sort(unique(d_c$Treatment)), collapse=", "), "\n")
cat("Unique dates:", paste(sort(unique(d_c$Date)), collapse=", "), "\n")
cat("Mice per treatment:\n")
print(table(d_c$Treatment))

# Calculate volume and dates
d_c <- mouseExperiment::calculate_volume(d_c,
  length_column = "Length", width_column = "Width")
d_c <- mouseExperiment::calculate_dates(d_c,
  date_column  = "Date",
  date_format  = "mdy",
  start_date   = "03/24/2025"
)
cat("Days after calculate_dates:", sort(unique(d_c$Day)), "\n")
cat("Volume sample (head):\n")
print(head(d_c[!is.na(d_c$Volume), c("Date","Day","ID","Treatment","Cage","Volume")], 6))

# Remove Day-0 zero volumes for LME4 (standard practice; keep zero for AUC)
d_c_nozero <- d_c[d_c$Volume > 0, ]

# ── Package: LME4 ─────────────────────────────────────────────────────────────
cat("\n--- Package LME4 ---\n")
pkg_lme4_c <- tryCatch(
  tumor_growth_statistics(
    df               = d_c_nozero,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    transform        = "log",
    model_type       = "lme4",
    reference_group  = "Control",
    plots            = FALSE,
    verbose          = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_lme4_c)) {
  cat("ANOVA:\n")
  print(pkg_lme4_c$anova)
  cat("\nPairwise vs Control (from pairwise_comparisons):\n")
  if (!is.null(pkg_lme4_c$pairwise_comparisons)) {
    print(as.data.frame(pkg_lme4_c$pairwise_comparisons)[,
      intersect(c("contrast","estimate","SE","df","t.ratio","p.value"),
                colnames(as.data.frame(pkg_lme4_c$pairwise_comparisons)))])
  } else cat("  NULL\n")
  cat("\nTreatment effects:\n")
  print(as.data.frame(pkg_lme4_c$treatment_effects))
} 

# ── Manual: LME4 ──────────────────────────────────────────────────────────────
cat("\n--- Manual LME4 ---\n")
d_c_nozero$logVol    <- log(d_c_nozero$Volume)
d_c_nozero$Treatment <- relevel(factor(d_c_nozero$Treatment), ref = "Control")
d_c_nozero$ID        <- factor(d_c_nozero$ID)

m_c <- tryCatch(
  lme4::lmer(logVol ~ Treatment * Day + (1 | ID), data = d_c_nozero, REML = FALSE),
  error = function(e) { cat("lmer ERROR:", e$message, "\n"); NULL }
)

if (!is.null(m_c)) {
  cat("Type III ANOVA:\n")
  print(car::Anova(m_c, type = 3))
  emm_c <- emmeans::emmeans(m_c, ~ Treatment)
  cat("\nvs Control (Dunnett):\n")
  dun_c <- summary(emmeans::contrast(emm_c, method = "trt.vs.ctrl",
                                     ref = which(levels(emm_c)[[1]] == "Control"),
                                     adjust = "dunnett"),
                   infer = c(TRUE, TRUE))
  print(as.data.frame(dun_c)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
  cat("\nMarginal means:\n")
  print(summary(emm_c))
}

# ── Package: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Package AUC ---\n")
pkg_auc_c <- tryCatch(
  tumor_growth_statistics(
    df               = d_c,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    model_type       = "auc",
    reference_group  = "Control",
    plots            = FALSE,
    verbose          = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_auc_c)) {
  cat("AUC treatment effects:\n")
  print(as.data.frame(pkg_auc_c$treatment_effects))
  cat("\nAUC pairwise:\n")
  if (!is.null(pkg_auc_c$posthoc$pairwise)) print(pkg_auc_c$posthoc$pairwise) else cat("  NULL\n")
}

# ── Manual: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Manual AUC ---\n")
manual_auc_c <- d_c %>%
  arrange(ID, Day) %>%
  group_by(ID, Treatment, Cage) %>%
  summarize(AUC = {
    tt <- Day; vv <- Volume; n <- length(tt)
    if (n < 2) NA_real_ else sum(diff(tt) * (vv[-n] + vv[-1]) / 2)
  }, .groups = "drop") %>%
  filter(!is.na(AUC))

cat("Group mean AUC:\n")
print(as.data.frame(manual_auc_c %>%
  group_by(Treatment) %>%
  summarize(Mean_AUC=round(mean(AUC),2), SD=round(sd(AUC),2), N=n(), .groups="drop")))

# ── Package: Survival ─────────────────────────────────────────────────────────
cat("\n--- Package Survival ---\n")
# Need last-day censoring info per animal
d_c_surv <- d_c %>%
  group_by(ID, Treatment, Cage) %>%
  summarize(
    MaxDay = max(Day),
    Event  = max(Survival_Censor),
    .groups = "drop"
  )
cat("Events per group:\n")
print(table(d_c_surv$Treatment, d_c_surv$Event))

pkg_surv_c <- tryCatch(
  survival_statistics(
    df               = d_c,
    time_column      = "Day",
    censor_column    = "Survival_Censor",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    reference_group  = "Control"
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_surv_c)) {
  cat("\nMethod used:", pkg_surv_c$method_used, "\n")
  cat("Survival results:\n")
  print(as.data.frame(pkg_surv_c$results))
}

# ── Manual: Survival ─────────────────────────────────────────────────────────
cat("\n--- Manual Survival (log-rank) ---\n")
# Use last observation per animal
d_c_surv_manual <- d_c %>%
  group_by(ID, Treatment) %>%
  summarize(
    time  = max(Day),
    event = max(Survival_Censor),
    .groups = "drop"
  )
manual_surv_c <- tryCatch({
  s <- Surv(d_c_surv_manual$time, d_c_surv_manual$event)
  sf <- survfit(s ~ Treatment, data = d_c_surv_manual)
  cat("Kaplan-Meier median survival:\n")
  print(summary(sf)$table[, c("records","events","median","0.95LCL","0.95UCL")])
  lr <- survival::survdiff(s ~ Treatment, data = d_c_surv_manual)
  cat("\nLog-rank test:\n")
  print(lr)
  cat("Log-rank p =", 1 - pchisq(lr$chisq, df = length(unique(d_c_surv_manual$Treatment))-1), "\n")
  lr
}, error = function(e) { cat("ERROR:", e$message, "\n"); NULL })

results$c_lme4_pkg <- pkg_lme4_c
results$c_surv_pkg <- pkg_surv_c

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 3 – dose_levels_synthetic_data
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n", sep, "DATASET 3: dose_levels_synthetic_data\n", sep)

d_d <- read.csv(file.path(DATA, "dose_levels_synthetic_data.csv"))
cat("Treatments:", paste(sort(unique(d_d$Treatment)), collapse=", "), "\n")
cat("Doses:     ", paste(sort(unique(d_d$Dose)), collapse=", "), "\n")
cat("Dates:     ", paste(sort(unique(d_d$Date)), collapse=", "), "\n")
cat("Mice per treatment×dose:\n")
print(table(d_d$Treatment, d_d$Dose))

d_d <- mouseExperiment::calculate_volume(d_d,
  length_column = "Length", width_column = "Width")
d_d <- mouseExperiment::calculate_dates(d_d,
  date_column  = "Date",
  date_format  = "dmy",
  start_date   = "24-Mar"
)
cat("Days:", sort(unique(d_d$Day)), "\n")
cat("Volume sample:\n")
print(head(d_d[d_d$Volume > 0, c("Date","Day","ID","Treatment","Dose","Cage","Volume")], 6))

d_d_nozero <- d_d[d_d$Volume > 0, ]

# ── Package: LME4 ─────────────────────────────────────────────────────────────
cat("\n--- Package LME4 ---\n")
pkg_lme4_d <- tryCatch(
  tumor_growth_statistics(
    df               = d_d_nozero,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    dose_column      = "Dose",
    transform        = "log",
    model_type       = "lme4",
    reference_group  = "Control",
    plots            = FALSE,
    verbose          = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_lme4_d)) {
  cat("ANOVA:\n")
  print(pkg_lme4_d$anova)
  cat("\nPairwise vs Control:\n")
  if (!is.null(pkg_lme4_d$pairwise_comparisons)) {
    print(as.data.frame(pkg_lme4_d$pairwise_comparisons)[,
      intersect(c("contrast","estimate","SE","df","t.ratio","p.value"),
                colnames(as.data.frame(pkg_lme4_d$pairwise_comparisons)))])
  } else cat("  NULL\n")
  cat("\nTreatment effects:\n")
  print(as.data.frame(pkg_lme4_d$treatment_effects))
}

# ── Manual: LME4 ──────────────────────────────────────────────────────────────
cat("\n--- Manual LME4 ---\n")
d_d_nozero$logVol    <- log(d_d_nozero$Volume)
d_d_nozero$Treatment <- relevel(factor(d_d_nozero$Treatment), ref = "Control")
d_d_nozero$ID        <- factor(d_d_nozero$ID)

m_d <- tryCatch(
  lme4::lmer(logVol ~ Treatment * Day + (1 | ID), data = d_d_nozero, REML = FALSE),
  error = function(e) { cat("lmer ERROR:", e$message, "\n"); NULL }
)

if (!is.null(m_d)) {
  cat("Type III ANOVA:\n")
  print(car::Anova(m_d, type = 3))
  emm_d <- emmeans::emmeans(m_d, ~ Treatment)
  cat("\nvs Control (Dunnett):\n")
  dun_d <- summary(
    emmeans::contrast(emm_d, method = "trt.vs.ctrl",
                      ref = which(levels(emm_d)[[1]] == "Control"),
                      adjust = "dunnett"),
    infer = c(TRUE, TRUE)
  )
  print(as.data.frame(dun_d)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
}

# ── Package: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Package AUC ---\n")
pkg_auc_d <- tryCatch(
  tumor_growth_statistics(
    df               = d_d,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    dose_column      = "Dose",
    model_type       = "auc",
    reference_group  = "Control",
    plots            = FALSE,
    verbose          = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_auc_d)) {
  cat("AUC treatment effects:\n")
  print(as.data.frame(pkg_auc_d$treatment_effects))
  cat("\nAUC pairwise:\n")
  if (!is.null(pkg_auc_d$posthoc$pairwise)) print(pkg_auc_d$posthoc$pairwise) else cat("  NULL\n")
}

# ── Manual: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Manual AUC ---\n")
manual_auc_d <- d_d %>%
  arrange(ID, Day) %>%
  group_by(ID, Treatment, Dose) %>%
  summarize(AUC = {
    tt <- Day; vv <- Volume; n <- length(tt)
    if (n < 2) NA_real_ else sum(diff(tt) * (vv[-n] + vv[-1]) / 2)
  }, .groups = "drop") %>%
  filter(!is.na(AUC))

cat("Group mean AUC:\n")
print(as.data.frame(manual_auc_d %>%
  group_by(Treatment, Dose) %>%
  summarize(Mean_AUC=round(mean(AUC),2), SD=round(sd(AUC),2), N=n(), .groups="drop")))

# ── Package: Survival ─────────────────────────────────────────────────────────
cat("\n--- Package Survival ---\n")
pkg_surv_d <- tryCatch(
  survival_statistics(
    df               = d_d,
    time_column      = "Day",
    censor_column    = "Survival_Censor",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    reference_group  = "Control"
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_surv_d)) {
  cat("Method:", pkg_surv_d$method_used, "\n")
  print(as.data.frame(pkg_surv_d$results))
}

# ── Manual: Survival ─────────────────────────────────────────────────────────
cat("\n--- Manual Survival ---\n")
d_d_surv <- d_d %>%
  group_by(ID, Treatment) %>%
  summarize(time=max(Day), event=max(Survival_Censor), .groups="drop")
tryCatch({
  s2 <- Surv(d_d_surv$time, d_d_surv$event)
  sf2 <- survfit(s2 ~ Treatment, data=d_d_surv)
  cat("KM medians:\n")
  print(summary(sf2)$table[, c("records","events","median","0.95LCL","0.95UCL")])
  lr2 <- survdiff(s2 ~ Treatment, data=d_d_surv)
  cat("Log-rank p =", round(1-pchisq(lr2$chisq, df=length(unique(d_d_surv$Treatment))-1), 4), "\n")
}, error=function(e) cat("ERROR:", e$message, "\n"))

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 4 – necrotic_synthetic_data
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n", sep, "DATASET 4: necrotic_synthetic_data\n", sep)

d_n <- read.csv(file.path(DATA, "necrotic_synthetic_data.csv"))
cat("Columns:", paste(colnames(d_n), collapse=", "), "\n")
cat("Treatments:", paste(sort(unique(d_n$Treatment)), collapse=", "), "\n")
cat("Days:", sort(unique(d_n$Day)), "\n")
cat("Necrotic:", paste(sort(unique(d_n$Necrotic)), collapse=", "), "\n")
cat("Mice per treatment:\n")
print(table(d_n$Treatment, d_n$Day))

# Compute volume and standardise column names
d_n$Volume    <- d_n$Tumor.Length..mm. * (d_n$Tumor.Width..mm.)^2 / 2
d_n$ID        <- d_n$Mouse.Tag
d_n$Cage      <- "1"

# Necrotic is a binary flag — treated as covariate or filter
# The package probably filters or notes it; let's compare both
cat("\nNecrotic observations:\n")
print(d_n[d_n$Necrotic == "Y", c("Day","ID","Treatment","Volume")])

d_n_nonec  <- d_n[d_n$Necrotic != "Y", ]
d_n_nozero <- d_n_nonec[d_n_nonec$Volume > 0, ]

# ── Package: LME4 (all data including necrotic) ───────────────────────────────
cat("\n--- Package LME4 (raw) ---\n")
pkg_lme4_n <- tryCatch(
  tumor_growth_statistics(
    df               = d_n,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    transform        = "log",
    model_type       = "lme4",
    reference_group  = "None",
    plots            = FALSE,
    verbose          = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_lme4_n)) {
  cat("ANOVA:\n"); print(pkg_lme4_n$anova)
  cat("\nPairwise vs None:\n")
  if (!is.null(pkg_lme4_n$pairwise_comparisons)) {
    print(as.data.frame(pkg_lme4_n$pairwise_comparisons)[,
      intersect(c("contrast","estimate","SE","df","t.ratio","p.value"),
                colnames(as.data.frame(pkg_lme4_n$pairwise_comparisons)))])
  } else cat("  NULL\n")
}

# ── Package: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Package AUC ---\n")
pkg_auc_n <- tryCatch(
  tumor_growth_statistics(
    df               = d_n,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    model_type       = "auc",
    reference_group  = "None",
    plots            = FALSE,
    verbose          = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_auc_n)) {
  cat("AUC treatment effects:\n")
  print(as.data.frame(pkg_auc_n$treatment_effects))
  cat("\nAUC pairwise:\n")
  if (!is.null(pkg_auc_n$posthoc$pairwise)) print(pkg_auc_n$posthoc$pairwise) else cat("  NULL\n")
}

# ── Manual: LME4 ──────────────────────────────────────────────────────────────
cat("\n--- Manual LME4 ---\n")
d_n2 <- d_n
d_n2$logVol    <- log(d_n2$Volume + 1e-6)
d_n2$Treatment <- relevel(factor(d_n2$Treatment), ref = "None")
d_n2$ID        <- factor(d_n2$ID)
m_n <- tryCatch(
  lme4::lmer(logVol ~ Treatment * Day + (1 | ID), data = d_n2, REML = FALSE),
  error = function(e) { cat("lmer ERROR:", e$message, "\n"); NULL }
)
if (!is.null(m_n)) {
  cat("Type III ANOVA:\n"); print(car::Anova(m_n, type=3))
  emm_n <- emmeans::emmeans(m_n, ~ Treatment)
  cat("\nAll pairwise (Tukey):\n")
  print(as.data.frame(summary(pairs(emm_n, adjust="tukey"), infer=c(TRUE,TRUE)))[,
             c("contrast","estimate","SE","df","t.ratio","p.value")])
}

# ── Manual: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Manual AUC ---\n")
manual_auc_n <- d_n %>%
  arrange(ID, Day) %>%
  group_by(ID, Treatment) %>%
  summarize(AUC = {
    tt <- Day; vv <- Volume; n <- length(tt)
    if (n < 2) NA_real_ else sum(diff(tt) * (vv[-n] + vv[-1]) / 2)
  }, .groups = "drop") %>%
  filter(!is.na(AUC))
cat("Per-animal AUC:\n"); print(as.data.frame(manual_auc_n))
cat("Group mean AUC:\n")
print(as.data.frame(manual_auc_n %>%
  group_by(Treatment) %>%
  summarize(Mean_AUC=round(mean(AUC),3), SD=round(sd(AUC),3), N=n(), .groups="drop")))
cat("\nOne-way ANOVA on AUC:\n")
print(summary(aov(AUC ~ Treatment, data=manual_auc_n)))
cat("\nPairwise Welch t-tests (Bonferroni):\n")
print(pairwise.t.test(manual_auc_n$AUC, manual_auc_n$Treatment,
                      pool.sd=FALSE, p.adjust.method="bonferroni"))

cat("\n\n====== DONE ======\n")
