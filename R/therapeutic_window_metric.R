#' Therapeutic window: efficacy and tolerability, side by side
#'
#' Reports two quantities per arm, each with a 95% interval: tumour growth
#' inhibition (TGI) at the evaluation day, and the worst body-weight loss,
#' measured in each animal from its own baseline and averaged over the arm's
#' animals. A tolerability flag compares the weight-loss interval with a
#' declared threshold.
#'
#' Before v0.27.0 the function reported a single ratio, TWM = TGI / weight
#' loss, and ranked the arms by it. The maintainer replaced it with this
#' two-axis summary (CODE_REVIEW.md R20-K). The ratio needed an arbitrary
#' floor for arms that lose no weight, and that floor sat below the scales'
#' own noise (R20.70). Its parts came from different estimators, so it could
#' not have a matching interval (R20.2). And one number cannot show both
#' efficacy and toxicity.
#'
#' @param df Data frame with longitudinal data.
#' @param weight_column Name of the body weight column.
#' @param volume_column Name of the tumor volume column.
#' @param time_column Name of the time/day column.
#' @param treatment_column Name of the treatment group column.
#' @param id_column Name of the mouse/subject ID column.
#' @param cage_column Name of the cage column. NULL to omit. Part of the
#'   composite mouse key so reused IDs across cages are not collapsed.
#' @param adjust_tumor_weight Logical; subtract estimated tumor weight.
#' @param tumor_density Density in g/cm³ (default 1.0).
#' @param volume_units Units of the volume column, \code{"mm3"} or \code{"cm3"}.
#'   Required when \code{adjust_tumor_weight = TRUE} and a volume column is used,
#'   because volume is converted to mass and subtracted from body weight: the
#'   wrong unit scales that correction by 1000. Units were inferred from the data
#'   before v0.25.0, which read small-tumour mm3 studies as cm3 (CODE_REVIEW.md
#'   R20.83); the data are now only checked against the declared unit, with a
#'   warning on disagreement and an error when the implied tumour mass exceeds
#'   half the body weight.
#' @param reference_group Name of the control/reference group.
#' @param tolerability_threshold Weight loss, in percent of baseline weight,
#'   that the tolerability flag tests against. Default 20. A value of 1 or less
#'   is an error, because it looks like a fraction (0.2 for 20 percent).
#' @param endpoint_day Day at which efficacy is evaluated. \code{NULL}
#'   (default): the last day on which every arm is evaluable (see
#'   \code{\link{evaluable_days}}); a requested day must be evaluable. Before
#'   v0.26.0 the default was the last observed day (CODE_REVIEW.md R20.1).
#' @param endpoint_method How the endpoint volume per arm is obtained
#'   (CODE_REVIEW.md R3.5 / G.3):
#'   \code{"model"} (default) takes each arm's geometric mean at
#'   \code{endpoint_day} from a mixed model of log volume fitted to every
#'   observation (natural spline in time per arm, per-animal random slopes),
#'   the same TGI as \code{\link{analyze_drug_synergy}}, so animals euthanised
#'   earlier still contribute;
#'   \code{"last_obs"} uses each animal's own last observation at or before the
#'   endpoint day; \code{"survivors"} reproduces the pre-0.8.0 raw mean among
#'   animals observed at the endpoint day, which conditions on survival and
#'   biases TGI downward — it warns when animals were lost.
#' @param n_boot Integer >= 0. Draws behind the 95% intervals. Under the model
#'   estimand the TGI intervals come from draws of the endpoint model's fixed
#'   effects; under the per-animal estimands, from a bootstrap of animals
#'   (control included). Weight-loss intervals always come from a bootstrap of
#'   animals within arm (CODE_REVIEW.md R3.6 / R3.7, R20.2). Default 2000; 0
#'   skips the intervals, and with them the tolerability flag.
#' @param boot_seed Optional integer seed for reproducible resampling.
#' @return A list:
#'   \describe{
#'     \item{window_table}{One row per arm, the reference arm first:
#'       \code{Treatment}, \code{N_Animals} (animals with weights),
#'       \code{TGI}, \code{TGI_Lower}, \code{TGI_Upper},
#'       \code{Worst_Loss} (the arm's mean of each animal's worst percent
#'       weight loss), \code{Worst_Loss_Lower}, \code{Worst_Loss_Upper},
#'       \code{N_Over_Threshold} (animals whose worst loss reached the
#'       threshold) and \code{Tolerability}.}
#'     \item{weight_loss_data}{One row per animal: \code{MouseKey},
#'       \code{Treatment}, \code{ID}, \code{Cage}, \code{Baseline_Weight},
#'       \code{Nadir_Weight}, \code{Nadir_Day}, \code{Pct_Loss} and
#'       \code{Over_Threshold}.}
#'     \item{tgi_data}{The endpoint means behind the TGI.}
#'     \item{tolerability_threshold, tolerability_rule}{The threshold, and
#'       the flag's rule in words.}
#'     \item{endpoint_day, endpoint_method, evaluability, endpoint_model,
#'       attrition, n_at_endpoint}{The evaluation day, the estimand, the
#'       evaluable-day record (the rule, the days used and the days and arms
#'       it excluded), the endpoint model, and the animals on study per arm.}
#'     \item{interval_note}{Where the intervals come from.}
#'   }
#'
#' @section Reading the table:
#' The two axes are read separately; the table is not a ranking.
#' \code{Tolerability} is "Tolerated" when the upper bound of the arm's mean
#' worst loss is below the threshold, "Not tolerated" when the lower bound is
#' at or above it, and "Unclear" when the interval contains it; it is
#' \code{NA} when there is no interval (\code{n_boot = 0}, or fewer than two
#' animals). The flag describes the arm's average animal:
#' \code{N_Over_Threshold} counts the animals that crossed the threshold on
#' their own, which can happen in a "Tolerated" arm. Worst loss is taken over
#' each animal's whole record, not only up to the evaluation day, because
#' toxicity after that day still counts.
#' @export
therapeutic_window_metric <- function(df,
                                      weight_column    = "Weight",
                                      volume_column    = "Volume",
                                      time_column      = "Day",
                                      treatment_column = "Treatment",
                                      id_column        = "ID",
                                      cage_column      = NULL,
                                      adjust_tumor_weight = TRUE,
                                      tumor_density    = 1.0,
                                      volume_units     = NULL,
                                      reference_group  = NULL,
                                      tolerability_threshold = 20,
                                      endpoint_day     = NULL,
                                      endpoint_method  = c("model", "last_obs", "survivors"),
                                      n_boot           = 2000L,
                                      boot_seed        = NULL) {

  endpoint_method <- match.arg(endpoint_method)
  if (!is.numeric(tolerability_threshold) || length(tolerability_threshold) != 1L ||
      !is.finite(tolerability_threshold) || tolerability_threshold >= 100) {
    stop("tolerability_threshold must be one number, in percent of baseline ",
         "weight (default 20).", call. = FALSE)
  }
  if (tolerability_threshold <= 1) {
    stop("tolerability_threshold is in percent of baseline weight: use 20 for ",
         "20 percent, not 0.2.", call. = FALSE)
  }

  # --- Validate ---
  required <- c(weight_column, volume_column, time_column, treatment_column, id_column)
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # --- Build working data ---
  cage_col_present <- !is.null(cage_column) && cage_column %in% names(df)
  cage_vec <- if (cage_col_present) as.character(df[[cage_column]]) else "1"
  wd <- data.frame(
    ID        = as.character(df[[id_column]]),
    Treatment = as.character(df[[treatment_column]]),
    Cage      = cage_vec,
    Day       = as.numeric(df[[time_column]]),
    Weight    = as.numeric(df[[weight_column]]),
    Volume    = as.numeric(df[[volume_column]]),
    stringsAsFactors = FALSE
  )
  # Composite mouse key prevents collapsing same-numeric-ID mice across cages
  wd$MouseKey <- make_mouse_key(wd$Treatment, wd$ID, wd$Cage)

  if (adjust_tumor_weight) {
    # CODE_REVIEW.md R3.30 — resolve units explicitly rather than assuming mm³.
    # R20.22: volume filled in on weighing days without a calliper reading.
    wd$Volume <- me_fill_volume(wd$MouseKey, wd$Day, wd$Volume)
    volume_units <- resolve_volume_units(wd$Volume, volume_units)
    tumor_mass <- volume_to_mass(wd$Volume, tumor_density, volume_units)
    check_tumor_mass_plausible(tumor_mass, wd$Weight, volume_units)
    wd$Weight <- wd$Weight - tumor_mass
  }

  # R20.22: weight rows are never filtered on volume -- a weighing without a
  # same-day calliper reading still counts toward the worst weight loss.
  wd <- wd[!is.na(wd$Weight) & !is.na(wd$Day), ]
  wd <- wd[order(wd$MouseKey, wd$Day), ]

  groups <- unique(wd$Treatment)
  if (is.null(reference_group)) {
    # Pick first alphabetically or common control names
    ctrl_patterns <- c("control", "vehicle", "dmso", "pbs", "saline", "placebo")
    ref_match <- groups[tolower(groups) %in% ctrl_patterns]
    reference_group <- if (length(ref_match) > 0) ref_match[1] else groups[1]
  }

  # --- Efficacy axis: TGI per arm ---
  # CODE_REVIEW.md R3.5 / G.3 -- each arm's geometric mean at the endpoint day
  # from the endpoint model fitted to every observation of every animal, the
  # same TGI as analyze_drug_synergy() and dose_response_statistics() (R20-K).
  # It uses every volume measurement, with or without a weight on that day.
  # R20.1: the default day is the last one on which every arm is evaluable,
  # and a requested day must be evaluable.
  vd <- data.frame(MouseKey  = make_mouse_key(as.character(df[[treatment_column]]),
                                              as.character(df[[id_column]]),
                                              cage_vec),
                   Treatment = as.character(df[[treatment_column]]),
                   Day       = as.numeric(df[[time_column]]),
                   Volume    = as.numeric(df[[volume_column]]),
                   stringsAsFactors = FALSE)
  vd <- vd[is.finite(vd$Day) & is.finite(vd$Volume), , drop = FALSE]
  ep <- endpoint_volumes(
    vd, id_column = "MouseKey", treatment_column = "Treatment",
    time_column = "Day", volume_column = "Volume",
    endpoint_day = endpoint_day, endpoint_method = endpoint_method,
    arms = sort(unique(vd$Treatment))
  )
  max_day  <- ep$endpoint_day
  tgi_data <- endpoint_tgi(ep$group_means, reference_group)
  model_based <- identical(ep$method, "model")

  # --- Tolerability axis: each animal's worst weight loss ---
  # R15.2: each animal's baseline is its own first weighing, not the global
  # first day, which dropped late-enrolled animals from the denominator and
  # understated the loss.
  baseline <- me_per_mouse_baseline(wd, c("MouseKey", "Treatment"), "Weight")
  nadir <- do.call(rbind, lapply(split(wd, wd$MouseKey), function(s) {
    i <- which.min(s$Weight)
    data.frame(MouseKey = s$MouseKey[1], ID = s$ID[1], Cage = s$Cage[1],
               Nadir_Weight = s$Weight[i], Nadir_Day = s$Day[i],
               stringsAsFactors = FALSE)
  }))
  mouse_wl <- merge(baseline, nadir, by = "MouseKey")
  mouse_wl$Pct_Loss <- pmax(
    (mouse_wl$Baseline_Weight - mouse_wl$Nadir_Weight) /
      mouse_wl$Baseline_Weight * 100, 0)   # weight gain is no loss
  mouse_wl$Over_Threshold <- mouse_wl$Pct_Loss >= tolerability_threshold
  mouse_wl <- mouse_wl[, c("MouseKey", "Treatment", "ID", "Cage",
                           "Baseline_Weight", "Nadir_Weight", "Nadir_Day",
                           "Pct_Loss", "Over_Threshold")]

  arms <- sort(unique(c(tgi_data$Treatment, mouse_wl$Treatment)))
  arms <- c(intersect(reference_group, arms), setdiff(arms, reference_group))
  wl_by <- split(mouse_wl, factor(mouse_wl$Treatment, levels = arms))
  win <- data.frame(
    Treatment        = arms,
    N_Animals        = vapply(wl_by, nrow, integer(1L)),
    TGI              = tgi_data$TGI[match(arms, tgi_data$Treatment)],
    Worst_Loss       = vapply(wl_by, function(g)
      if (nrow(g)) mean(g$Pct_Loss) else NA_real_, numeric(1L)),
    N_Over_Threshold = vapply(wl_by, function(g) sum(g$Over_Threshold),
                              integer(1L)),
    stringsAsFactors = FALSE)
  rownames(win) <- NULL

  # --- Intervals (R20.2: each describes its own point estimate) ---
  n_boot <- as.integer(n_boot)
  tgi_ci <- if (n_boot >= 2L) {
    if (model_based) {
      twm_model_tgi_ci(ep$model, tgi_data$Treatment, reference_group, max_day,
                       n_draws = n_boot, seed = boot_seed)
    } else {
      twm_animal_tgi_ci(ep$per_mouse, reference_group, n_boot = n_boot,
                        seed = boot_seed)
    }
  }
  wl_ci <- if (n_boot >= 2L) twm_wl_bootstrap(mouse_wl, n_boot = n_boot,
                                              seed = boot_seed)
  pick <- function(ci, col) {
    if (is.null(ci)) return(rep(NA_real_, length(arms)))
    ci[[col]][match(arms, ci$Treatment)]
  }
  win$TGI_Lower        <- pick(tgi_ci, "TGI_Lower")
  win$TGI_Upper        <- pick(tgi_ci, "TGI_Upper")
  win$Worst_Loss_Lower <- pick(wl_ci, "WL_Lower")
  win$Worst_Loss_Upper <- pick(wl_ci, "WL_Upper")
  lo <- win$Worst_Loss_Lower
  hi <- win$Worst_Loss_Upper
  win$Tolerability <- ifelse(!is.finite(lo) | !is.finite(hi), NA_character_,
                      ifelse(hi < tolerability_threshold, "Tolerated",
                      ifelse(lo >= tolerability_threshold, "Not tolerated",
                             "Unclear")))
  win <- win[, c("Treatment", "N_Animals", "TGI", "TGI_Lower", "TGI_Upper",
                 "Worst_Loss", "Worst_Loss_Lower", "Worst_Loss_Upper",
                 "N_Over_Threshold", "Tolerability")]

  thr <- format(tolerability_threshold)
  list(
    window_table     = win,
    weight_loss_data = mouse_wl,
    tgi_data         = tgi_data,
    tolerability_threshold = tolerability_threshold,
    tolerability_rule = paste0(
      "Tolerated: the upper 95% bound of the arm's mean worst weight loss is ",
      "below ", thr, "%. Not tolerated: its lower bound is at or above ", thr,
      "%. Unclear: the interval contains ", thr, "%."),
    # Per-group n at the endpoint day -- the number that makes survivor
    # attrition visible instead of implicit (see R3.5).
    n_at_endpoint    = ep$attrition,
    attrition        = ep$attrition,
    endpoint_day     = max_day,
    endpoint_method  = ep$method,
    # The evaluable-day record (R20-K) and how the TGI was modelled.
    evaluability     = ep$evaluability,
    endpoint_model   = me_endpoint_model_info(ep$model),
    interval_note    = if (n_boot < 2L) {
      "No intervals (n_boot = 0), so no tolerability flag."
    } else if (model_based) {
      paste("TGI intervals: draws from the endpoint model's fixed effects.",
            "Weight-loss intervals: bootstrap of animals within arm.")
    } else {
      "TGI and weight-loss intervals: bootstrap of animals within arm."
    }
  )
}

#' TGI intervals from the endpoint model
#' @noRd
#' @keywords internal
twm_model_tgi_ci <- function(em, arms, reference_group, t, n_draws = 2000L,
                             seed = NULL) {
  if (is.null(em) || n_draws < 2L || !reference_group %in% arms) return(NULL)
  lm_ <- me_draw_logmeans(em, me_beta_draws(em, n_draws, seed), arms, t)
  ref <- lm_[, reference_group]
  do.call(rbind, lapply(arms, function(a) {
    if (identical(a, reference_group)) {
      return(data.frame(Treatment = a, TGI_Lower = 0, TGI_Upper = 0))
    }
    q <- stats::quantile(100 * (1 - exp(lm_[, a] - ref)), c(0.025, 0.975),
                         names = FALSE)
    data.frame(Treatment = a, TGI_Lower = q[1], TGI_Upper = q[2])
  }))
}

#' Bootstrap intervals for each arm's mean weight loss
#' @noRd
#' @keywords internal
twm_wl_bootstrap <- function(mouse_wl, n_boot = 2000L, seed = NULL) {
  if (n_boot < 2L) return(NULL)
  by <- split(as.numeric(mouse_wl$Pct_Loss), mouse_wl$Treatment)
  by <- lapply(by, function(v) v[is.finite(v)])
  me_with_seed(seed, do.call(rbind, lapply(names(by), function(g) {
    v <- by[[g]]
    if (length(v) < 2L) {
      return(data.frame(Treatment = g, WL_Lower = NA_real_, WL_Upper = NA_real_))
    }
    m <- vapply(seq_len(n_boot), function(i) mean(sample(v, length(v), TRUE)),
                numeric(1L))
    q <- stats::quantile(m, c(0.025, 0.975), names = FALSE)
    data.frame(Treatment = g, WL_Lower = q[1], WL_Upper = q[2])
  })))
}

#' TGI intervals from a bootstrap of animals (per-animal estimands)
#'
#' @param per_mouse `endpoint_volumes()$per_mouse`: Treatment, Volume.
#' @param reference_group Control arm name.
#' @param n_boot,seed Resampling controls.
#' @return Data frame with Treatment, TGI_Lower, TGI_Upper, or NULL.
#' @noRd
#' @keywords internal
twm_animal_tgi_ci <- function(per_mouse, reference_group, n_boot = 2000L,
                              seed = NULL) {
  if (is.null(per_mouse) || n_boot < 2L) return(NULL)
  vol_by <- split(as.numeric(per_mouse$Volume), per_mouse$Treatment)
  vol_by <- lapply(vol_by, function(v) v[is.finite(v)])
  if (!reference_group %in% names(vol_by) ||
      length(vol_by[[reference_group]]) < 2L) return(NULL)
  rs <- function(v) mean(sample(v, length(v), replace = TRUE))
  me_with_seed(seed, do.call(rbind, lapply(names(vol_by), function(g) {
    if (identical(g, reference_group)) {
      return(data.frame(Treatment = g, TGI_Lower = 0, TGI_Upper = 0))
    }
    v <- vol_by[[g]]
    if (length(v) < 2L) {
      return(data.frame(Treatment = g, TGI_Lower = NA_real_, TGI_Upper = NA_real_))
    }
    d <- vapply(seq_len(n_boot), function(i) {
      ctrl <- rs(vol_by[[reference_group]])
      if (!is.finite(ctrl) || ctrl <= 0) NA_real_ else 100 * (1 - rs(v) / ctrl)
    }, numeric(1L))
    q <- stats::quantile(d, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    data.frame(Treatment = g, TGI_Lower = q[1], TGI_Upper = q[2])
  })))
}
