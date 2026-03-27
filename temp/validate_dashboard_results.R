#!/usr/bin/env Rscript
# validate_dashboard_results.R
# Replicates dashboard processing exactly, then compares against independent
# manual R calculations for all four demo datasets.

suppressWarnings(suppressMessages({
  library(mouseExperiment)
  library(lme4)
  library(emmeans)
  library(car)
  library(survival)
  library(dplyr)
}))

DATA <- "/home/david/Projects/mouseExperimentDashboard/inst/sample_data"
sep  <- paste0(strrep("\u2550", 70), "\n")

# helper: compare two scalars and flag discrepancies
cmp <- function(label, a, b, tol = 1e-4) {
  if (is.na(a) || is.na(b)) {
    cat(sprintf("  %-45s  pkg=%-12s  manual=%-12s  [NA]\n",
                label, a, b))
  } else if (abs(a - b) / (abs(b) + 1e-10) > tol) {
    cat(sprintf("  %-45s  pkg=%-12.4g  manual=%-12.4g  *** MISMATCH ***\n",
                label, a, b))
  } else {
    cat(sprintf("  %-45s  pkg=%-12.4g  manual=%-12.4g  OK\n",
                label, a, b))
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 1 – weight_synthetic_data.csv
# Dashboard ref_date = "2025-03-01"; date format "March 1" → append year 2025
# Volume column is "Tumor.Volume" (pre-measured, no calculation)
# ══════════════════════════════════════════════════════════════════════════════
cat(sep, "DATASET 1: weight_synthetic_data\n", sep)

d_w_raw <- read.csv(file.path(DATA, "weight_synthetic_data.csv"))
# Replicate dashboard date handling: append year to "March X" dates
d_w_raw$Date_orig <- d_w_raw$Date
d_w_raw$Date      <- paste(d_w_raw$Date, "2025")   # "March 1 2025"
d_w <- mouseExperiment::calculate_dates(d_w_raw,
  date_column  = "Date",
  date_format  = "%B %d %Y",
  start_date   = "March 01 2025"
)
cat("Days:", sort(unique(d_w$Day)), "\n")

# Dashboard uses Tumor.Volume directly (not Length/Width); no cage column
d_w$Volume <- d_w$Tumor.Volume
d_w$ID     <- d_w$Ear.Tag
d_w$Cage   <- "1"    # dummy – no cage in this dataset

cat("Volume sample (first 3 rows):\n")
print(head(d_w[, c("Date_orig","Day","ID","Treatment","Volume")], 3))

d_w_nozero <- d_w[d_w$Volume > 0, ]

# ── Package: LME4 ─────────────────────────────────────────────────────────────
cat("\n--- Package LME4 ---\n")
pkg_w_lme4 <- tryCatch(
  tumor_growth_statistics(
    df               = d_w_nozero,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    transform        = "log",
    model_type       = "lme4",
    reference_group  = "DMSO",
    plots            = FALSE, verbose = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_w_lme4)) {
  cat("ANOVA:\n");                print(pkg_w_lme4$anova)
  cat("\nTreatment effects:\n");  print(as.data.frame(pkg_w_lme4$treatment_effects))
  pcs_w <- tryCatch(as.data.frame(pkg_w_lme4$pairwise_comparisons), error=function(e) NULL)
  if (!is.null(pcs_w)) {
    cat("\nPairwise vs DMSO:\n")
    print(pcs_w[, intersect(c("contrast","estimate","SE","df","t.ratio","p.value"), colnames(pcs_w))])
  }
}

# ── Manual: LME4 ──────────────────────────────────────────────────────────────
cat("\n--- Manual LME4 ---\n")
d_w2 <- d_w_nozero
d_w2$logVol    <- log(d_w2$Volume)
d_w2$Treatment <- relevel(factor(d_w2$Treatment), ref = "DMSO")
d_w2$ID        <- factor(d_w2$ID)

m_w <- tryCatch(
  lme4::lmer(logVol ~ Treatment * Day + (1 | ID), data = d_w2, REML = FALSE),
  error = function(e) { cat("lmer ERROR:", e$message, "\n"); NULL }
)

if (!is.null(m_w)) {
  cat("Type III ANOVA:\n");  print(car::Anova(m_w, type = 3))
  emm_w <- emmeans::emmeans(m_w, ~ Treatment)
  cat("\nMarginal means:\n"); print(summary(emm_w))
  dun_w <- summary(
    emmeans::contrast(emm_w, method = "trt.vs.ctrl",
                      ref = which(levels(emm_w)[[1]] == "DMSO"),
                      adjust = "dunnett"),
    infer = c(TRUE, TRUE)
  )
  cat("\nvs DMSO Dunnett:\n")
  print(as.data.frame(dun_w)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
  pw_w <- summary(pairs(emm_w, adjust = "tukey"), infer = c(TRUE, TRUE))
  cat("\nAll pairwise Tukey:\n")
  print(as.data.frame(pw_w)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
}

# ── Package: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Package AUC ---\n")
pkg_w_auc <- tryCatch(
  tumor_growth_statistics(
    df               = d_w,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    model_type       = "auc",
    reference_group  = "DMSO",
    plots            = FALSE, verbose = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_w_auc)) {
  cat("AUC ANOVA:\n");              print(pkg_w_auc$anova)
  cat("\nAUC treatment effects:\n"); print(as.data.frame(pkg_w_auc$treatment_effects))
  cat("\nAUC pairwise (posthoc$pairwise):\n")
  if (!is.null(pkg_w_auc$posthoc$pairwise)) print(pkg_w_auc$posthoc$pairwise)
}

# ── Manual: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Manual AUC ---\n")
manual_auc_w <- d_w %>% arrange(ID, Day) %>%
  group_by(ID, Treatment) %>%
  summarize(AUC = {
    n <- length(Day); if(n<2) NA_real_
    else sum(diff(Day) * (Volume[-n] + Volume[-1]) / 2)
  }, .groups="drop")
cat("Per-animal AUC:\n"); print(as.data.frame(manual_auc_w))
manual_auc_grp_w <- manual_auc_w %>%
  group_by(Treatment) %>%
  summarize(Mean_AUC=mean(AUC,na.rm=TRUE), SD=sd(AUC,na.rm=TRUE), N=n(), .groups="drop")
cat("\nGroup AUC means:\n"); print(as.data.frame(manual_auc_grp_w))
cat("\nOne-way ANOVA on AUC:\n")
print(summary(aov(AUC ~ Treatment, data = manual_auc_w)))
cat("\nWelch pairwise t-tests (Bonferroni):\n")
print(pairwise.t.test(manual_auc_w$AUC, manual_auc_w$Treatment,
                      pool.sd=FALSE, p.adjust.method="bonferroni"))

# ── Comparison: AUC treatment effects ──────────────────────────────────────
if (!is.null(pkg_w_auc)) {
  cat("\n--- AUC Group Means Comparison ---\n")
  pkg_te <- as.data.frame(pkg_w_auc$treatment_effects)
  for (grp in manual_auc_grp_w$Treatment) {
    pkg_val <- pkg_te$Mean_AUC[pkg_te$Treatment == grp]
    man_val <- manual_auc_grp_w$Mean_AUC[manual_auc_grp_w$Treatment == grp]
    if (length(pkg_val) && length(man_val))
      cmp(paste("Mean_AUC", grp), pkg_val, man_val)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 2 – combo_treatment_synthetic_data.csv
# Dashboard ref_date = "2025-03-24"; date format "%m/%d/%Y"
# Volume = (π/6) × L × W²  (ellipsoid, default)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n", sep, "DATASET 2: combo_treatment_synthetic_data\n", sep)

d_c_raw <- read.csv(file.path(DATA, "combo_treatment_synthetic_data.csv"))
cat("Raw treatments:", paste(sort(unique(d_c_raw$Treatment)), collapse=", "), "\n")
cat("Raw dates (first 3):", paste(head(unique(d_c_raw$Date), 3), collapse=", "), "\n")

d_c <- mouseExperiment::calculate_volume(d_c_raw,
  length_column = "Length", width_column = "Width", formula = "ellipsoid")
d_c <- mouseExperiment::calculate_dates(d_c,
  date_column  = "Date",
  date_format  = "%m/%d/%Y",
  start_date   = "03/24/2025"
)
cat("Days:", sort(unique(d_c$Day)), "\n")
cat("Volume sample:\n")
print(head(d_c[d_c$Volume > 0, c("Date","Day","ID","Treatment","Cage","Length","Width","Volume")], 4))

# Quick manual check: V = (π/6)*L*W²
samp <- head(d_c[d_c$Volume > 0, ], 1)
cat(sprintf("Volume check: (π/6)*%.2f*%.2f² = %.4f  (stored: %.4f)\n",
  samp$Length, samp$Width, (pi/6)*samp$Length*samp$Width^2, samp$Volume))

d_c_nozero <- d_c[d_c$Volume > 0, ]

# ── Package: LME4 ─────────────────────────────────────────────────────────────
cat("\n--- Package LME4 ---\n")
pkg_c_lme4 <- tryCatch(
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
    plots            = FALSE, verbose = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)

if (!is.null(pkg_c_lme4)) {
  cat("ANOVA:\n"); print(pkg_c_lme4$anova)
  cat("\nTreatment effects:\n"); print(as.data.frame(pkg_c_lme4$treatment_effects))
  pcs_c <- tryCatch(as.data.frame(pkg_c_lme4$pairwise_comparisons), error=function(e) NULL)
  if (!is.null(pcs_c)) {
    cat("\nPairwise vs Control:\n")
    print(pcs_c[, intersect(c("contrast","estimate","SE","df","t.ratio","p.value"), colnames(pcs_c))])
  }
}

# ── Manual: LME4 ──────────────────────────────────────────────────────────────
cat("\n--- Manual LME4 ---\n")
d_c2 <- d_c_nozero
d_c2$logVol    <- log(d_c2$Volume)
d_c2$Treatment <- relevel(factor(d_c2$Treatment), ref = "Control")
d_c2$ID        <- factor(d_c2$ID)
m_c <- tryCatch(
  lme4::lmer(logVol ~ Treatment * Day + (1 | ID), data = d_c2, REML = FALSE),
  error = function(e) { cat("lmer ERROR:", e$message, "\n"); NULL }
)
if (!is.null(m_c)) {
  cat("Type III ANOVA:\n"); print(car::Anova(m_c, type = 3))
  emm_c <- emmeans::emmeans(m_c, ~ Treatment)
  cat("\nMarginal means:\n"); print(summary(emm_c))
  dun_c <- summary(
    emmeans::contrast(emm_c, method="trt.vs.ctrl",
                      ref=which(levels(emm_c)[[1]]=="Control"), adjust="dunnett"),
    infer=c(TRUE,TRUE))
  cat("\nvs Control (Dunnett):\n")
  print(as.data.frame(dun_c)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
  pw_c <- summary(pairs(emm_c, adjust="tukey"), infer=c(TRUE,TRUE))
  cat("\nAll pairwise (Tukey):\n")
  print(as.data.frame(pw_c)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
}

# ── Package: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Package AUC ---\n")
pkg_c_auc <- tryCatch(
  tumor_growth_statistics(
    df               = d_c,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    model_type       = "auc",
    reference_group  = "Control",
    plots            = FALSE, verbose = FALSE
  ),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)
if (!is.null(pkg_c_auc)) {
  cat("AUC treatment effects:\n"); print(as.data.frame(pkg_c_auc$treatment_effects))
  cat("\nAUC pairwise:\n")
  if (!is.null(pkg_c_auc$posthoc$pairwise)) print(pkg_c_auc$posthoc$pairwise)
}

# ── Manual: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Manual AUC ---\n")
manual_auc_c <- d_c %>% arrange(ID, Day) %>%
  group_by(ID, Treatment) %>%
  summarize(AUC={n<-length(Day);if(n<2)NA_real_ else sum(diff(Day)*(Volume[-n]+Volume[-1])/2)},
            .groups="drop") %>% filter(!is.na(AUC))
cat("Group AUC means:\n")
manual_auc_grp_c <- manual_auc_c %>%
  group_by(Treatment) %>%
  summarize(Mean_AUC=round(mean(AUC),2),SD=round(sd(AUC),2),N=n(),.groups="drop")
print(as.data.frame(manual_auc_grp_c))
cat("\nOne-way ANOVA:\n"); print(summary(aov(AUC ~ Treatment, data=manual_auc_c)))
cat("\nWelch pairwise t-tests (Bonferroni):\n")
print(pairwise.t.test(manual_auc_c$AUC, manual_auc_c$Treatment,
                      pool.sd=FALSE, p.adjust.method="bonferroni"))

# ── Package: Survival ─────────────────────────────────────────────────────────
cat("\n--- Package Survival ---\n")
pkg_c_surv <- tryCatch(
  survival_statistics(df=d_c, time_column="Day", censor_column="Survival_Censor",
    treatment_column="Treatment", cage_column="Cage", id_column="ID",
    reference_group="Control"),
  error=function(e){cat("ERROR:",e$message,"\n");NULL})
if (!is.null(pkg_c_surv)) {
  cat("Method:", pkg_c_surv$method_used, "\n")
  print(as.data.frame(pkg_c_surv$results))
}

# ── Manual: Survival ─────────────────────────────────────────────────────────
cat("\n--- Manual Survival ---\n")
d_c_sv <- d_c %>% group_by(ID, Treatment) %>%
  summarize(time=max(Day), event=max(Survival_Censor), .groups="drop")
tryCatch({
  s_c <- Surv(d_c_sv$time, d_c_sv$event)
  sf_c <- survfit(s_c ~ Treatment, data=d_c_sv)
  cat("KM medians:\n")
  print(summary(sf_c)$table[, c("records","events","median","0.95LCL","0.95UCL")])
  lr_c <- survdiff(s_c ~ Treatment, data=d_c_sv)
  cat("Log-rank p =", round(1-pchisq(lr_c$chisq, df=length(unique(d_c_sv$Treatment))-1),4),"\n")
}, error=function(e) cat("ERROR:", e$message, "\n"))

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 3 – dose_levels_synthetic_data.csv
# Dashboard ref_date = "2025-03-24"; date format "%d-%b" + year 2025
# Volume = (π/6)*L*W²
# Dashboard creates TreatmentGroup = paste(Treatment, Dose)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n", sep, "DATASET 3: dose_levels_synthetic_data\n", sep)

d_d_raw <- read.csv(file.path(DATA, "dose_levels_synthetic_data.csv"))
cat("Treatments:", paste(sort(unique(d_d_raw$Treatment)), collapse=", "), "\n")
cat("Doses:     ", paste(sort(unique(d_d_raw$Dose)), collapse=", "), "\n")
cat("Date sample:", head(unique(d_d_raw$Date), 3), "\n")

d_d <- mouseExperiment::calculate_volume(d_d_raw,
  length_column="Length", width_column="Width", formula="ellipsoid")
d_d <- mouseExperiment::calculate_dates(d_d,
  date_column="Date", date_format="%d-%b", year=2025,
  start_date="24-Mar")
cat("Days:", sort(unique(d_d$Day)), "\n")

# Dashboard TreatmentGroup combining
d_d$TreatmentGroup <- paste(d_d$Treatment, d_d$Dose)
cat("TreatmentGroup values:", paste(sort(unique(d_d$TreatmentGroup)), collapse=", "), "\n")
cat("Volume sample:\n")
print(head(d_d[d_d$Volume>0, c("Date","Day","ID","TreatmentGroup","Cage","Volume")], 4))

d_d_nozero <- d_d[d_d$Volume > 0, ]

# ── Package: LME4 ─────────────────────────────────────────────────────────────
cat("\n--- Package LME4 (TreatmentGroup) ---\n")
ref_grp_d <- sort(unique(d_d$TreatmentGroup))[1]
cat("Reference group:", ref_grp_d, "\n")
pkg_d_lme4 <- tryCatch(
  tumor_growth_statistics(
    df               = d_d_nozero,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "TreatmentGroup",
    cage_column      = "Cage",
    id_column        = "ID",
    transform        = "log",
    model_type       = "lme4",
    reference_group  = ref_grp_d,
    plots            = FALSE, verbose = FALSE
  ),
  error=function(e){cat("ERROR:",e$message,"\n");NULL})
if (!is.null(pkg_d_lme4)) {
  cat("ANOVA:\n"); print(pkg_d_lme4$anova)
  cat("\nTreatment effects:\n"); print(as.data.frame(pkg_d_lme4$treatment_effects))
  pcs_d <- tryCatch(as.data.frame(pkg_d_lme4$pairwise_comparisons), error=function(e)NULL)
  if (!is.null(pcs_d)) {
    cat("\nPairwise vs", ref_grp_d, ":\n")
    print(pcs_d[, intersect(c("contrast","estimate","SE","df","t.ratio","p.value"), colnames(pcs_d))])
  }
}

# ── Manual: LME4 ──────────────────────────────────────────────────────────────
cat("\n--- Manual LME4 ---\n")
d_d2 <- d_d_nozero
d_d2$logVol         <- log(d_d2$Volume)
d_d2$TreatmentGroup <- relevel(factor(d_d2$TreatmentGroup), ref=ref_grp_d)
d_d2$ID             <- factor(d_d2$ID)
m_d <- tryCatch(
  lme4::lmer(logVol ~ TreatmentGroup * Day + (1 | ID), data=d_d2, REML=FALSE),
  error=function(e){cat("lmer ERROR:",e$message,"\n");NULL})
if (!is.null(m_d)) {
  cat("Type III ANOVA:\n"); print(car::Anova(m_d, type=3))
  emm_d <- emmeans::emmeans(m_d, ~ TreatmentGroup)
  cat("\nMarginal means:\n"); print(summary(emm_d))
  dun_d <- summary(
    emmeans::contrast(emm_d, method="trt.vs.ctrl",
                      ref=which(levels(emm_d)[[1]]==ref_grp_d), adjust="dunnett"),
    infer=c(TRUE,TRUE))
  cat("\nvs Reference (Dunnett):\n")
  print(as.data.frame(dun_d)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
}

# ── Package: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Package AUC ---\n")
pkg_d_auc <- tryCatch(
  tumor_growth_statistics(
    df               = d_d,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "TreatmentGroup",
    cage_column      = "Cage",
    id_column        = "ID",
    model_type       = "auc",
    reference_group  = ref_grp_d,
    plots            = FALSE, verbose = FALSE
  ),
  error=function(e){cat("ERROR:",e$message,"\n");NULL})
if (!is.null(pkg_d_auc)) {
  cat("AUC treatment effects:\n"); print(as.data.frame(pkg_d_auc$treatment_effects))
  cat("\nAUC pairwise:\n")
  if (!is.null(pkg_d_auc$posthoc$pairwise)) print(pkg_d_auc$posthoc$pairwise)
}

# ── Manual: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Manual AUC ---\n")
manual_auc_d <- d_d %>% arrange(ID, Day) %>%
  group_by(ID, TreatmentGroup) %>%
  summarize(AUC={n<-length(Day);if(n<2)NA_real_ else sum(diff(Day)*(Volume[-n]+Volume[-1])/2)},
            .groups="drop") %>% filter(!is.na(AUC))
cat("Group AUC means:\n")
print(as.data.frame(manual_auc_d %>%
  group_by(TreatmentGroup) %>%
  summarize(Mean_AUC=round(mean(AUC),2),SD=round(sd(AUC),2),N=n(),.groups="drop")))
cat("\nOne-way ANOVA:\n"); print(summary(aov(AUC ~ TreatmentGroup, data=manual_auc_d)))
cat("\nWelch pairwise t-tests (Bonferroni):\n")
print(pairwise.t.test(manual_auc_d$AUC, manual_auc_d$TreatmentGroup,
                      pool.sd=FALSE, p.adjust.method="bonferroni"))

# ── Package: Survival ─────────────────────────────────────────────────────────
cat("\n--- Package Survival ---\n")
# Survival uses Treatment column directly (not TreatmentGroup)
pkg_d_surv <- tryCatch(
  survival_statistics(df=d_d, time_column="Day", censor_column="Survival_Censor",
    treatment_column="Treatment", cage_column="Cage", id_column="ID",
    reference_group="Control"),
  error=function(e){cat("ERROR:",e$message,"\n");NULL})
if (!is.null(pkg_d_surv)) {
  cat("Method:", pkg_d_surv$method_used, "\n")
  print(as.data.frame(pkg_d_surv$results))
}

# ── Manual: Survival ─────────────────────────────────────────────────────────
cat("\n--- Manual Survival ---\n")
d_d_sv <- d_d %>% group_by(ID, Treatment) %>%
  summarize(time=max(Day), event=max(Survival_Censor), .groups="drop")
tryCatch({
  s_d <- Surv(d_d_sv$time, d_d_sv$event)
  sf_d <- survfit(s_d ~ Treatment, data=d_d_sv)
  cat("KM medians:\n")
  print(summary(sf_d)$table[, c("records","events","median","0.95LCL","0.95UCL")])
  lr_d <- survdiff(s_d ~ Treatment, data=d_d_sv)
  cat("Log-rank p =", round(1-pchisq(lr_d$chisq, df=length(unique(d_d_sv$Treatment))-1),4), "\n")
}, error=function(e) cat("ERROR:", e$message, "\n"))

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 4 – necrotic_synthetic_data.csv
# Dashboard ref_date = "" (empty); dataset has Day column already (numeric)
# Volume = (π/6)*L*W²  (dashboard calculates from Tumor.Length..mm./Width cols)
# The "Necrotic" flag is data to note, not a censor column per se
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n", sep, "DATASET 4: necrotic_synthetic_data\n", sep)

d_n_raw <- read.csv(file.path(DATA, "necrotic_synthetic_data.csv"))
cat("Columns:", paste(colnames(d_n_raw), collapse=", "), "\n")
cat("Treatments:", paste(sort(unique(d_n_raw$Treatment)), collapse=", "), "\n")
cat("Days:", sort(unique(d_n_raw$Day)), "\n")
cat("Necrotic Y/N:\n"); print(table(d_n_raw$Treatment, d_n_raw$Necrotic))

# Dashboard: Day column already exists; calculate volume from Length+Width cols
d_n <- mouseExperiment::calculate_volume(d_n_raw,
  length_column="Tumor.Length..mm.", width_column="Tumor.Width..mm.",
  formula="ellipsoid")
# Rename to standard names
colnames(d_n)[colnames(d_n)=="Mouse.Tag"] <- "ID"
d_n$Cage <- "1"

cat("Volume sample:\n")
print(head(d_n[d_n$Volume > 0, c("Day","ID","Treatment","Necrotic","Volume")], 5))

# Quick formula check
samp_n <- head(d_n[d_n$Volume > 0, ], 1)
cat(sprintf("Volume check: (π/6)*%.2f*%.2f² = %.4f (stored: %.4f)\n",
  samp_n$Tumor.Length..mm., samp_n$Tumor.Width..mm.,
  (pi/6)*samp_n$Tumor.Length..mm.*samp_n$Tumor.Width..mm.^2, samp_n$Volume))

d_n_nozero <- d_n[d_n$Volume > 0, ]

# ── Package: LME4 ─────────────────────────────────────────────────────────────
cat("\n--- Package LME4 ---\n")
pkg_n_lme4 <- tryCatch(
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
    plots            = FALSE, verbose = FALSE
  ),
  error=function(e){cat("ERROR:",e$message,"\n");NULL})
if (!is.null(pkg_n_lme4)) {
  cat("ANOVA:\n"); print(pkg_n_lme4$anova)
  cat("\nTreatment effects:\n"); print(as.data.frame(pkg_n_lme4$treatment_effects))
  pcs_n <- tryCatch(as.data.frame(pkg_n_lme4$pairwise_comparisons), error=function(e)NULL)
  if (!is.null(pcs_n)) {
    cat("\nPairwise vs None:\n")
    print(pcs_n[, intersect(c("contrast","estimate","SE","df","t.ratio","p.value"), colnames(pcs_n))])
  }
}

# ── Manual LME4 ───────────────────────────────────────────────────────────────
cat("\n--- Manual LME4 ---\n")
d_n2 <- d_n
d_n2$logVol    <- log(d_n2$Volume + 1e-6)
d_n2$Treatment <- relevel(factor(d_n2$Treatment), ref="None")
d_n2$ID        <- factor(d_n2$ID)
m_n <- tryCatch(
  lme4::lmer(logVol ~ Treatment * Day + (1 | ID), data=d_n2, REML=FALSE),
  error=function(e){cat("lmer ERROR:",e$message,"\n");NULL})
if (!is.null(m_n)) {
  cat("Type III ANOVA:\n"); print(car::Anova(m_n, type=3))
  emm_n <- emmeans::emmeans(m_n, ~ Treatment)
  cat("\nMarginal means:\n"); print(summary(emm_n))
  pw_n <- summary(pairs(emm_n, adjust="tukey"), infer=c(TRUE,TRUE))
  cat("\nAll pairwise (Tukey):\n")
  print(as.data.frame(pw_n)[, c("contrast","estimate","SE","df","t.ratio","p.value")])
}

# ── Package: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Package AUC ---\n")
pkg_n_auc <- tryCatch(
  tumor_growth_statistics(
    df               = d_n,
    time_column      = "Day",
    volume_column    = "Volume",
    treatment_column = "Treatment",
    cage_column      = "Cage",
    id_column        = "ID",
    model_type       = "auc",
    reference_group  = "None",
    plots            = FALSE, verbose = FALSE
  ),
  error=function(e){cat("ERROR:",e$message,"\n");NULL})
if (!is.null(pkg_n_auc)) {
  cat("AUC treatment effects:\n"); print(as.data.frame(pkg_n_auc$treatment_effects))
  cat("\nAUC pairwise:\n")
  if (!is.null(pkg_n_auc$posthoc$pairwise)) print(pkg_n_auc$posthoc$pairwise)
}

# ── Manual: AUC ──────────────────────────────────────────────────────────────
cat("\n--- Manual AUC ---\n")
manual_auc_n <- d_n %>% arrange(ID, Day) %>%
  group_by(ID, Treatment) %>%
  summarize(AUC={n<-length(Day);if(n<2)NA_real_ else sum(diff(Day)*(Volume[-n]+Volume[-1])/2)},
            .groups="drop") %>% filter(!is.na(AUC))
cat("Per-animal AUC:\n"); print(as.data.frame(manual_auc_n))
manual_auc_grp_n <- manual_auc_n %>%
  group_by(Treatment) %>%
  summarize(Mean_AUC=round(mean(AUC),3),SD=round(sd(AUC),3),N=n(),.groups="drop")
cat("\nGroup AUC means:\n"); print(as.data.frame(manual_auc_grp_n))
cat("\nOne-way ANOVA:\n"); print(summary(aov(AUC ~ Treatment, data=manual_auc_n)))
cat("\nWelch pairwise t-tests (Bonferroni):\n")
print(pairwise.t.test(manual_auc_n$AUC, manual_auc_n$Treatment,
                      pool.sd=FALSE, p.adjust.method="bonferroni"))

cat("\n\n====== ALL ANALYSES COMPLETE ======\n")
