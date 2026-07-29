#' Therapeutic Window Metric (TWM)
#'
#' Computes TWM = max(TGI, 0) / max(MeanWeightLoss%, noise_floor) per treatment
#' group, so a well-tolerated arm scores its TGI and the metric stays continuous
#' at the floor.
#' The denominator uses the group mean of per-mouse maximum weight loss, not the
#' single worst-case mouse, to give a representative measure of typical toxicity.
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
#' @param volume_units Units of \code{volume_column}: "mm3", "cm3", or NULL
#'   (default) to auto-detect from the data and report the inference.
#' @param reference_group Name of the control/reference group.
#' @param noise_floor Minimum group-mean per-mouse maximum weight loss (percent)
#'   used as a floor on the ratio's denominator, so
#'   \code{TWM = max(TGI, 0) / max(weight_loss\%, noise_floor)}. This keeps TWM
#'   continuous for any \code{noise_floor} (CODE_REVIEW.md R3.18) and reduces to
#'   \code{TWM = max(TGI, 0)} for well-tolerated arms. Default 1.0 — an
#'   experimentally pragmatic threshold that avoids dividing by near-zero weight
#'   loss. No formal clinical basis; users with experiment-specific noise
#'   estimates (e.g. scale precision) should tune accordingly.
#' @param endpoint_day Day at which efficacy is evaluated. \code{NULL}
#'   (default) uses the maximum observed day.
#' @param endpoint_method How the endpoint volume per arm is obtained
#'   (CODE_REVIEW.md R3.5 / G.3):
#'   \code{"model"} (default) takes each arm's geometric mean at
#'   \code{endpoint_day} from a log-scale mixed model fitted to every
#'   observation, so animals euthanised earlier still contribute;
#'   \code{"last_obs"} uses each animal's own last observation at or before the
#'   endpoint day; \code{"survivors"} reproduces the pre-0.8.0 raw mean among
#'   animals observed at the endpoint day, which conditions on survival and
#'   biases TGI downward — it warns when animals were lost.
#' @param n_boot Integer >= 0. Mouse-level bootstrap resamples used to attach
#'   95\% percentile intervals to TGI, mean weight loss, and TWM. Mice are
#'   resampled within group including the control arm, so the interval
#'   propagates the TGI denominator's uncertainty (CODE_REVIEW.md R3.6 / R3.7).
#'   Default 2000; set 0 to skip.
#' @param boot_seed Optional integer seed for reproducible resampling.
#' @return A list with: \code{twm_table} (point estimates plus bootstrap
#'   bounds), \code{twm_ci}, \code{tgi_data}, \code{weight_loss_data},
#'   \code{n_at_endpoint} (animals contributing at the endpoint day, per arm),
#'   and \code{endpoint_day}.
#'
#' @section Interpreting the ranking:
#' \code{twm_table} is sorted by TWM, which reads as a ranking. Check
#' \code{TWM_Lower} / \code{TWM_Upper} before acting on the order: with typical
#' group sizes the intervals overlap heavily and the ordering of adjacent arms is
#' often not resolvable. \code{n_at_endpoint} shows how many animals in each arm
#' actually reached the endpoint day.
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
                                      noise_floor      = 1.0,
                                      endpoint_day     = NULL,
                                      endpoint_method  = c("model", "last_obs", "survivors"),
                                      n_boot           = 2000L,
                                      boot_seed        = NULL) {

  endpoint_method <- match.arg(endpoint_method)

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
  wd$MouseKey <- make_mouse_key(wd$ID, wd$Treatment, wd$Cage)

  if (adjust_tumor_weight) {
    # CODE_REVIEW.md R3.30 — resolve units explicitly rather than assuming mm³.
    volume_units <- resolve_volume_units(wd$Volume, volume_units)
    tumor_mass <- volume_to_mass(wd$Volume, tumor_density, volume_units)
    check_tumor_mass_plausible(tumor_mass, wd$Weight, volume_units)
    wd$Weight <- wd$Weight - tumor_mass
  }

  wd <- wd[!is.na(wd$Weight) & !is.na(wd$Day) & !is.na(wd$Volume), ]
  wd <- wd[order(wd$ID, wd$Day), ]

  groups <- unique(wd$Treatment)
  if (is.null(reference_group)) {
    # Pick first alphabetically or common control names
    ctrl_patterns <- c("control", "vehicle", "dmso", "pbs", "saline", "placebo")
    ref_match <- groups[tolower(groups) %in% ctrl_patterns]
    reference_group <- if (length(ref_match) > 0) ref_match[1] else groups[1]
  }

  # --- TGI per treatment group ---
  # CODE_REVIEW.md R3.5 / G.3 — this used to take the raw mean among animals
  # still observed at the global last day, which conditions on survival and
  # biases TGI downward (hardest against the control arm, which loses animals
  # first). The default now reads each arm's geometric mean at the endpoint day
  # off a log-scale LMM fitted to every observation from every animal.
  ep <- endpoint_volumes(
    wd, id_column = "MouseKey", treatment_column = "Treatment",
    time_column = "Day", volume_column = "Volume",
    endpoint_day = endpoint_day, endpoint_method = endpoint_method
  )
  max_day  <- ep$endpoint_day
  tgi_data <- endpoint_tgi(ep$group_means, reference_group)

  # Per-mouse endpoint volumes for the bootstrap. The model path has no
  # per-mouse draw, so resample the last-observation values, which use every
  # animal too.
  final <- if (!is.null(ep$per_mouse)) {
    ep$per_mouse
  } else {
    endpoint_volumes(wd, id_column = "MouseKey", treatment_column = "Treatment",
                     time_column = "Day", volume_column = "Volume",
                     endpoint_day = endpoint_day,
                     endpoint_method = "last_obs")$per_mouse
  }

  # --- Max % weight loss per group ---
  # Per mouse: baseline weight, nadir weight, max % loss
  # Filter to the earliest study day before aggregating so x[1] is ordered.
  # Aggregate by MouseKey (ID|||Treatment|||Cage) so reused IDs across cages
  # don't collapse — same fix class as Round 1 1.8 for weight_corrected_tgi.
  min_day <- min(wd$Day, na.rm = TRUE)
  baseline <- stats::aggregate(Weight ~ MouseKey + Treatment,
                               data = wd[wd$Day == min_day, ],
                               FUN = mean, na.rm = TRUE)
  names(baseline)[3] <- "Baseline_Weight"

  nadir <- stats::aggregate(Weight ~ MouseKey + Treatment, data = wd,
                            FUN = min, na.rm = TRUE)
  names(nadir)[3] <- "Nadir_Weight"

  mouse_wl <- merge(baseline, nadir, by = c("MouseKey", "Treatment"))
  mouse_wl$Pct_Loss <- (mouse_wl$Baseline_Weight - mouse_wl$Nadir_Weight) /
                        mouse_wl$Baseline_Weight * 100
  mouse_wl$Pct_Loss <- pmax(mouse_wl$Pct_Loss, 0)  # clamp negative (weight gain)

  # Group mean of per-mouse max weight loss: more representative of typical
  # toxicity than the single worst-case mouse (previously used max here).
  group_wl <- stats::aggregate(Pct_Loss ~ Treatment, data = mouse_wl,
                               FUN = mean, na.rm = TRUE)
  names(group_wl)[2] <- "Mean_Pct_Weight_Loss"

  # --- TWM ---
  twm <- merge(tgi_data, group_wl, by = "Treatment")
  # Clamp negative TGI to 0 — a treatment that *accelerates* tumour growth
  # has no efficacy benefit, so TWM should not be positive for it. Using
  # abs(TGI) here previously would have made a treatment that enhanced
  # disease AND caused weight loss appear safe (TWM > 0); now it scores 0.
  tgi_pos <- pmax(twm$TGI, 0)
  # CODE_REVIEW.md R3.18 — the previous two-branch form returned TGI (percentage
  # points) below the floor and TGI/WL% (a dimensionless ratio) above it, then
  # sorted both into one ranking. The two branches agree at the boundary only
  # when noise_floor == 1.0, because dividing by 1 is the identity — the
  # continuity was an accident of the default value, and any other setting
  # introduced a jump discontinuity that reordered treatments at the threshold.
  # Flooring the denominator is continuous for every noise_floor and reduces to
  # the old behaviour exactly at 1.0.
  twm$TWM <- tgi_pos / pmax(twm$Mean_Pct_Weight_Loss, noise_floor)
  twm$Safety_Note <- ifelse(
    twm$Mean_Pct_Weight_Loss <= noise_floor,
    "Negligible weight loss",
    ""
  )
  twm <- twm[order(-twm$TWM), ]

  # CODE_REVIEW.md R3.6 / R3.7 / G.6 — TWM is a ratio of two group means that
  # was previously reported as a bare point estimate and then *sorted*, so the
  # table read as a ranking of treatments with no indication of how much of the
  # ordering was noise. Resample mice within group (control included, since it
  # is the TGI denominator) and recompute TGI, mean weight loss, and TWM from
  # scratch on each resample.
  twm_ci <- if (n_boot > 0L) {
    twm_bootstrap(final, mouse_wl, reference_group, noise_floor,
                  n_boot = as.integer(n_boot), seed = boot_seed)
  } else NULL

  if (!is.null(twm_ci)) {
    twm <- merge(twm, twm_ci, by = "Treatment", all.x = TRUE)
    twm <- twm[order(-twm$TWM), ]
  }

  # Per-group n at the endpoint day — the number that makes survivor attrition
  # visible instead of implicit (see R3.5).
  n_at_endpoint <- ep$attrition

  list(
    twm_table        = twm,
    twm_ci           = twm_ci,
    tgi_data         = tgi_data,
    weight_loss_data = mouse_wl,
    n_at_endpoint    = n_at_endpoint,
    attrition        = ep$attrition,
    endpoint_day     = max_day,
    endpoint_method  = ep$method
  )
}

#' Mouse-level bootstrap CIs for TGI, mean weight loss, and TWM
#'
#' @param final Endpoint-day rows (one per surviving mouse) with `Treatment`,
#'   `Volume`.
#' @param mouse_wl Per-mouse weight-loss table with `Treatment`, `Pct_Loss`.
#' @param reference_group Control arm name.
#' @param noise_floor Denominator floor, as in the point estimate.
#' @param n_boot,seed Resampling controls.
#' @return Data frame with per-treatment percentile intervals, or NULL.
#' @noRd
#' @keywords internal
twm_bootstrap <- function(final, mouse_wl, reference_group, noise_floor,
                          n_boot = 2000L, seed = NULL) {
  groups <- unique(final$Treatment)
  if (!reference_group %in% groups || n_boot < 2L) return(NULL)

  vol_by  <- split(as.numeric(final$Volume),    final$Treatment)
  loss_by <- split(as.numeric(mouse_wl$Pct_Loss), mouse_wl$Treatment)
  vol_by  <- lapply(vol_by,  function(v) v[is.finite(v)])
  loss_by <- lapply(loss_by, function(v) v[is.finite(v)])
  if (any(vapply(vol_by, length, integer(1)) < 2L)) return(NULL)

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }

  rs <- function(v) mean(sample(v, length(v), replace = TRUE))
  draws <- lapply(seq_len(n_boot), function(i) {
    ctrl <- rs(vol_by[[reference_group]])
    if (!is.finite(ctrl) || ctrl <= 0) return(NULL)
    vapply(groups, function(g) {
      tgi <- if (identical(g, reference_group)) 0 else
        (1 - rs(vol_by[[g]]) / ctrl) * 100
      wl  <- if (length(loss_by[[g]]) >= 2L) rs(loss_by[[g]]) else NA_real_
      c(TGI = tgi, WL = wl, TWM = max(tgi, 0) / max(wl, noise_floor))
    }, numeric(3))
  })
  draws <- draws[!vapply(draws, is.null, logical(1))]
  if (length(draws) < 2L) return(NULL)

  do.call(rbind, lapply(seq_along(groups), function(gi) {
    q <- function(metric) {
      v <- vapply(draws, function(d) d[metric, gi], numeric(1))
      v <- v[is.finite(v)]
      if (length(v) < 2L) return(c(NA_real_, NA_real_))
      stats::quantile(v, c(0.025, 0.975), names = FALSE)
    }
    tgi_q <- q("TGI"); wl_q <- q("WL"); twm_q <- q("TWM")
    data.frame(
      Treatment    = groups[gi],
      TGI_Lower    = tgi_q[1], TGI_Upper    = tgi_q[2],
      WL_Lower     = wl_q[1],  WL_Upper     = wl_q[2],
      TWM_Lower    = twm_q[1], TWM_Upper    = twm_q[2],
      Boot_N       = length(draws),
      stringsAsFactors = FALSE
    )
  }))
}
