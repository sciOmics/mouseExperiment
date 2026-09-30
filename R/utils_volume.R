# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Volume-unit handling shared by the toxicity functions.
#
# CODE_REVIEW.md R3.30 — analyze_body_weight(), weight_loss_threshold() and
# therapeutic_window_metric() each hard-coded `Volume / 1000 * tumor_density`,
# which is correct only for mm3. On a cm3 upload the subtracted tumour mass came
# out 1000x too small, so Net_Weight was effectively unadjusted body weight while
# still being labelled "Net Weight (body - tumor)". Because the entire point of
# the adjustment is to stop a large tumour masking treatment-induced weight loss,
# the silent no-op inverted the conclusion for exactly the animals it matters
# most for. Units are now explicit, auto-detectable, and range-checked.

#' Infer whether tumour volumes are in mm3 or cm3
#'
#' Mouse tumour volumes in mm3 are conventionally in the tens to low thousands;
#' the same tumours in cm3 are fractions of a unit. The two scales differ by
#' 1000, so a median-based split is unambiguous in practice.
#'
#' @param volume Numeric vector of tumour volumes.
#' @return `"mm3"`, `"cm3"`, or `NA_character_` when there are no usable values.
#' @noRd
#' @keywords internal
detect_volume_units <- function(volume) {
  v <- volume[is.finite(volume) & volume > 0]
  if (length(v) == 0L) return(NA_character_)
  # CODE_REVIEW.md R20.83 -- this was the median over every observation, which
  # read small-tumour mm3 studies (many early, non-palpable measurements) as
  # cm3: the dashboard's weight demo, 0.1-1,092.6 mm3 with a median of 6.21,
  # came out cm3 and every net weight went 1000x wrong. The large tumours late
  # in a study set the scale, so use the 90th percentile. This is now only a
  # cross-check on the units the caller declares (resolve_volume_units()).
  if (stats::quantile(v, 0.9, names = FALSE) < 20) "cm3" else "mm3"
}

#' Convert tumour volume to estimated tumour mass in grams
#'
#' @param volume Numeric vector of tumour volumes.
#' @param tumor_density Density in g/cm3 (default 1.0, i.e. soft-tissue density).
#' @param volume_units `"mm3"` or `"cm3"`.
#' @return Numeric vector of estimated tumour masses in grams.
#' @noRd
#' @keywords internal
volume_to_mass <- function(volume, tumor_density = 1.0,
                           volume_units = c("mm3", "cm3")) {
  volume_units <- match.arg(volume_units)
  # 1 cm3 = 1000 mm3; mass (g) = volume (cm3) * density (g/cm3)
  volume_cm3 <- switch(volume_units, mm3 = volume / 1000, cm3 = volume)
  volume_cm3 * tumor_density
}

#' Resolve the volume unit, auto-detecting when the caller did not specify
#'
#' Emits a message when the unit is inferred so the choice is never invisible,
#' and warns when an explicitly supplied unit disagrees with what the data look
#' like — that disagreement is the signature of the R3.30 failure mode.
#'
#' @param volume Numeric vector of tumour volumes.
#' @param volume_units `NULL` to auto-detect, otherwise `"mm3"` / `"cm3"`.
#' @return The resolved unit string.
#' @noRd
#' @keywords internal
resolve_volume_units <- function(volume, volume_units = NULL) {
  # CODE_REVIEW.md R20.83 -- units are declared, not inferred. The inference
  # sent the dashboard's default Toxicity run to a net weight of -1,070 g.
  if (is.null(volume_units)) {
    stop("volume_units is required when tumour mass is subtracted from body ",
         "weight: pass volume_units = \"mm3\" or \"cm3\" (or set ",
         "adjust_tumor_weight = FALSE).", call. = FALSE)
  }
  volume_units <- match.arg(volume_units, c("mm3", "cm3"))
  detected <- detect_volume_units(volume)
  if (!is.na(detected) && detected != volume_units) {
    warning("volume_units was given as '", volume_units,
            "' but the data look like '", detected,
            "' (90th percentile of volume ", signif(stats::quantile(
              volume[is.finite(volume) & volume > 0], 0.9, names = FALSE), 3),
            "). Tumour-mass adjustment will be wrong by a factor of 1000 if ",
            "the declared unit is incorrect.", call. = FALSE)
  }
  volume_units
}

#' Sanity-check estimated tumour mass against body weight
#'
#' Catches a unit error from either direction: a tumour mass that is a large
#' fraction of body weight means the volume was probably cm3 treated as mm3, and
#' a mass that is negligible for a study-endpoint tumour means the reverse.
#'
#' @param tumor_mass Numeric vector of estimated tumour masses (g).
#' @param body_weight Numeric vector of body weights (g).
#' @param volume_units The unit used, for the message text.
#' @noRd
#' @keywords internal
check_tumor_mass_plausible <- function(tumor_mass, body_weight, volume_units) {
  ok <- is.finite(tumor_mass) & is.finite(body_weight) & body_weight > 0
  if (!any(ok)) return(invisible(NULL))
  frac <- tumor_mass[ok] / body_weight[ok]

  # R20.83 -- this was a warning, which the dashboard discarded, so a
  # 1000x unit error ran to "Analysis completed successfully!". A tumour
  # heavier than half the mouse is not a result; stop.
  if (max(frac, na.rm = TRUE) > 0.5) {
    stop("Estimated tumour mass exceeds 50% of body weight for at least one ",
         "animal (max ", round(100 * max(frac, na.rm = TRUE)), "%) with ",
         "volume_units = '", volume_units, "'. The volume units are probably ",
         "wrong: cm3 data treated as mm3 (or the reverse) is off by 1000x.",
         call. = FALSE)
  }
  invisible(NULL)
}
