# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Cage / treatment design structure.
#
# CODE_REVIEW.md R3.17 / G.2 — cage inclusion was previously decided by
# `chisq.test(table(cage, treatment))$p.value < 0.05`. Three things were wrong
# with that:
#
#   1. Whether cage is nested within treatment is a *structural* property of
#      the design, decided by the randomisation, not a hypothesis to test. A
#      significance test answers a different question and its answer depends on
#      sample size.
#   2. The table counted observations, not mice, so the statistic was inflated
#      by the number of timepoints and rejected for arbitrarily mild imbalance.
#   3. Most consequentially, it conflated two structurally different designs.
#      The maintainer's stated norm — one treatment per cage, two or more cages
#      per arm — has cage *estimable as a random intercept*, but the chi-square
#      flagged it "collinear" and the model dropped cage entirely. Ignoring
#      cage-level clustering that the design supports estimating understates
#      every standard error: with m mice per cage the design effect is
#      1 + (m-1)*ICC, so 4 mice/cage at ICC 0.1 understates SEs by ~14 %, and
#      at ICC 0.2 by ~26 %.
#
# The three cases and their correct handling:
#
#   crossed            cages hold >1 treatment       cage estimable as fixed or random
#   nested_replicated  >=2 cages/arm, 1 arm per cage cage NOT estimable as fixed
#                                                    (aliased); estimable as (1|Cage)
#   nested_confounded  exactly 1 cage per arm        cage not estimable at all; the
#                                                    cage effect is inseparable from
#                                                    the treatment effect

#' Classify the cage x treatment design structure
#'
#' Operates on **mice**, not observations — a mouse measured ten times must not
#' count ten times toward the design structure.
#'
#' @param df Data frame.
#' @param cage_column,id_column,treatment_column Column names.
#' @return A list with `structure` (one of `"none"`, `"crossed"`,
#'   `"nested_replicated"`, `"nested_confounded"`), `n_cages`,
#'   `cages_per_treatment` (named integer vector), `mice_per_cage` (named
#'   integer vector), `min_cages_per_treatment`, and a human-readable
#'   `description`.
#' @noRd
#' @keywords internal
classify_cage_structure <- function(df, cage_column, id_column,
                                    treatment_column) {
  none <- list(structure = "none", n_cages = 0L,
               cages_per_treatment = integer(0), mice_per_cage = integer(0),
               min_cages_per_treatment = NA_integer_,
               description = "No cage information supplied.")

  if (is.null(cage_column) || !cage_column %in% names(df)) return(none)
  if (length(unique(df[[cage_column]])) <= 1L) {
    none$description <- "All animals share a single cage label."
    return(none)
  }

  # One row per mouse — this is the unit the design structure lives at.
  mice <- unique(df[, c(id_column, treatment_column, cage_column), drop = FALSE])
  # R20.25: as labels, so an unused factor level is not counted as an arm with
  # no cages (it classified the whole study as confounded and dropped (1|cage)).
  ct   <- table(as.character(mice[[cage_column]]), as.character(mice[[treatment_column]]))

  treatments_per_cage <- rowSums(ct > 0)
  cages_per_treatment <- colSums(ct > 0)
  mice_per_cage       <- rowSums(ct)

  structure <- if (any(treatments_per_cage > 1L)) {
    "crossed"
  } else if (all(cages_per_treatment >= 2L)) {
    "nested_replicated"
  } else if (any(cages_per_treatment >= 2L)) {
    # R20.25: one arm in a single cage (a satellite arm, say) used to make the
    # whole study "confounded", and cage was dropped for every arm: the null
    # false-positive rate went from 0.02-0.035 to 0.16-0.29. The replicated
    # arms still estimate the cage variance.
    "nested_partial"
  } else {
    "nested_confounded"
  }

  description <- switch(structure,
    crossed = paste0(
      sum(treatments_per_cage > 1L), " of ", nrow(ct),
      " cages hold more than one treatment; cage is crossed with treatment ",
      "and is estimable as either a fixed or a random effect."),
    nested_replicated = paste0(
      "Each cage holds exactly one treatment, with ",
      paste(range(cages_per_treatment), collapse = "-"),
      " cages per arm. Cage is nested within treatment: it cannot be a fixed ",
      "effect (perfectly aliased with treatment) but IS estimable as a random ",
      "intercept, which is the correct way to account for cage-level ",
      "clustering."),
    nested_partial = paste0(
      "Each cage holds one treatment. ",
      paste(names(cages_per_treatment)[cages_per_treatment < 2L], collapse = ", "),
      if (sum(cages_per_treatment < 2L) == 1L) " is" else " are",
      " housed in a single cage, so for ",
      if (sum(cages_per_treatment < 2L) == 1L) "that arm" else "those arms",
      " cage and treatment cannot be separated; the other arms have two or ",
      "more cages, from which the cage variance is estimated."),
    nested_confounded = paste0(
      "At least one treatment is housed in a single cage (",
      paste(names(cages_per_treatment)[cages_per_treatment < 2L],
            collapse = ", "),
      "). Cage and treatment are completely confounded for those arms \u2014 any ",
      "cage effect is inseparable from the treatment effect and no model can ",
      "disentangle them.")
  )

  list(
    structure               = structure,
    n_cages                 = nrow(ct),
    cages_per_treatment     = cages_per_treatment,
    mice_per_cage           = mice_per_cage,
    min_cages_per_treatment = min(cages_per_treatment),
    description             = description
  )
}

#' Decide how cage enters the model, given the design and the user's request
#'
#' `handle_cage_effects = "auto"` (the new default) picks the statistically
#' correct treatment for the detected structure. The explicit options are
#' honoured where they are estimable and rejected where they are not, rather
#' than silently producing a rank-deficient design.
#'
#' @param cage_info Output of `classify_cage_structure()`.
#' @param handle_cage_effects One of `"auto"`, `"include_if_not_collinear"`,
#'   `"always_include"`, `"never_include"`, `"as_random_effect"`.
#' @param verbose Print the reasoning.
#' @return A list with `fixed` and `random` logicals and a `reason` string.
#' @noRd
#' @keywords internal
resolve_cage_handling <- function(cage_info, handle_cage_effects,
                                  verbose = FALSE) {
  st <- cage_info$structure

  if (st == "none") {
    return(list(fixed = FALSE, random = FALSE,
                reason = "No cage information; cage omitted."))
  }

  out <- switch(handle_cage_effects,
    auto = switch(st,
      crossed = list(
        fixed = TRUE, random = FALSE,
        reason = "Cage is crossed with treatment; included as a fixed effect."),
      nested_replicated = list(
        fixed = FALSE, random = TRUE,
        reason = paste0(
          "Cage is nested within treatment with replication (min ",
          cage_info$min_cages_per_treatment,
          " cages per arm); included as a random intercept (1|cage). ",
          "A fixed effect is not estimable here \u2014 it is aliased with ",
          "treatment.")),
      nested_partial = list(
        fixed = FALSE, random = TRUE,
        reason = paste0(
          "Cage is nested within treatment and replicated in some arms; ",
          "included as a random intercept (1|cage). For the single-cage arms ",
          "the treatment estimate includes any cage effect.")),
      nested_confounded = list(
        fixed = FALSE, random = FALSE,
        reason = paste0(
          "Cage is completely confounded with treatment; omitted. ",
          "Any cage effect is absorbed into the treatment estimate."))
    ),

    # Retained for backward compatibility. The old semantics were "include as a
    # fixed effect unless a chi-square says cage and treatment are associated",
    # which dropped cage in exactly the case where a random intercept was both
    # estimable and needed.
    include_if_not_collinear = if (st == "crossed") {
      list(fixed = TRUE, random = FALSE,
           reason = "Cage is crossed with treatment; included as a fixed effect.")
    } else {
      list(fixed = FALSE, random = FALSE,
           reason = paste0(
             "Cage is nested within treatment, so a fixed effect is not ",
             "estimable and this option omits it. NOTE: with ",
             cage_info$min_cages_per_treatment,
             "+ cages per arm, handle_cage_effects = 'auto' would fit ",
             "(1|cage) instead and account for cage-level clustering."))
    },

    always_include = if (st == "crossed") {
      list(fixed = TRUE, random = FALSE,
           reason = "Cage included as a fixed effect (crossed design).")
    } else {
      stop("handle_cage_effects = 'always_include' asks for cage as a fixed ",
           "effect, but the design is '", st, "': ", cage_info$description,
           "\nA fixed cage effect is not estimable. Use 'as_random_effect' ",
           "(or the default 'auto') instead.", call. = FALSE)
    },

    never_include = list(fixed = FALSE, random = FALSE,
                         reason = "Cage omitted at the caller's request."),

    as_random_effect = if (st == "nested_confounded") {
      # R20.25: one single-cage arm no longer counts as "confounded" (that is
      # nested_partial, below), so this branch is reached only when EVERY arm
      # is in a single cage. There the cage variance is not identifiable apart
      # from treatment, and fitting (1|cage) fails outright, so the request
      # cannot be honoured; say so rather than drop it quietly.
      warning("handle_cage_effects = 'as_random_effect' cannot be honoured: ",
              "every arm is housed in a single cage, so a cage variance cannot ",
              "be estimated apart from treatment. Cage is omitted.", call. = FALSE)
      list(fixed = FALSE, random = FALSE,
           reason = paste0(
             "Cage requested as a random effect, but every arm is in a single ",
             "cage (", cage_info$description, "); omitted."))
    } else {
      list(fixed = FALSE, random = TRUE,
           reason = "Cage included as a random intercept (1|cage).")
    }
  )

  if (st == "nested_partial" && !identical(handle_cage_effects, "never_include")) {
    warning("Some arms are housed in a single cage: ", cage_info$description,
            call. = FALSE)
  }
  if (st == "nested_confounded" && !identical(handle_cage_effects, "never_include")) {
    warning("Cage and treatment are completely confounded: ",
            cage_info$description,
            " Treatment effects therefore include any cage effect. ",
            "Housing each arm in >= 2 cages would make the two separable.",
            call. = FALSE)
  }
  if (isTRUE(verbose)) message("Cage handling: ", out$reason)
  out
}

#' Cage-level intraclass correlation from a fitted mixed model
#'
#' ICC_cage = var(cage) / (var(cage) + var(animal at day t) + var(residual)).
#' Reported so a user can see how much the clustering mattered: the design
#' effect on the effective sample size is 1 + (mice_per_cage - 1) * ICC.
#'
#' CODE_REVIEW.md R20.68: with random slopes the animal variance depends on
#' the day, Var(b0 + b1 t) = s0^2 + 2 t s01 + t^2 s1^2. Adding the slope
#' variance (per day^2) to the total made the ICC depend on the time units
#' (0.513 in days, 0.395 in weeks). It is now evaluated at a stated day, the
#' mean study day the Treatment test uses, which no change of units moves.
#'
#' @param model A fitted `merMod`.
#' @param cage_column Name of the cage grouping factor.
#' @param time_column Name of the time variable in the random slopes.
#' @param at Day at which the animal variance is evaluated.
#' @return One-row data frame, or NULL when not computable.
#' @noRd
#' @keywords internal
cage_icc <- function(model, cage_column, time_column = NULL, at = 0) {
  if (!inherits(model, "merMod")) return(NULL)
  vcl <- tryCatch(lme4::VarCorr(model), error = function(e) NULL)
  if (is.null(vcl)) return(NULL)
  grp_of <- function(nm) sub("\\.[0-9]+$", "", nm)   # "ID.1" from (Day || ID)
  var_at <- function(S) {
    S <- as.matrix(S)
    terms <- rownames(S)
    z <- ifelse(terms == "(Intercept)", 1,
                ifelse(!is.null(time_column) & terms == time_column, at, 0))
    as.numeric(t(z) %*% S %*% z)
  }
  parts <- vapply(vcl, var_at, numeric(1))
  is_cage <- grp_of(names(vcl)) == cage_column
  if (sum(is_cage) != 1L) return(NULL)
  cage_var <- parts[is_cage]
  total <- sum(parts) + stats::sigma(model)^2
  if (!is.finite(total) || total <= 0) return(NULL)
  data.frame(
    Grouping        = cage_column,
    Variance        = unname(cage_var),
    Total_Variance  = total,
    ICC             = unname(cage_var / total),
    At_Day          = at,
    stringsAsFactors = FALSE
  )
}

# CODE_REVIEW.md R22.3 -- cages numbered within each arm (C1 and C2 in every
# arm) read as cages that each hold several arms: the "crossed" design, with
# cage a fixed effect estimated across arms. From the labels alone that cannot
# be told apart from cages that really are shared, so the study says which
# (label_cages(), or Study Design on the dashboard). These flag the pattern.

#' The cage labels the arms reuse, when every arm's labels come from one set
#'
#' At least two arms share a label, and one arm's labels include every other
#' arm's (Vehicle C1-C3 with C1-C2 elsewhere counts). Unique labels, and
#' cages that each hold a different pair of arms, do not.
#' @param cage,arm Each row's cage and arm.
#' @return The labels, in order; empty when the pattern is absent.
#' @noRd
reused_cage_labels <- function(cage, arm) {
  ok <- !is.na(cage) & !is.na(arm)
  sets <- lapply(split(as.character(cage[ok]), as.character(arm[ok])), unique)
  if (length(sets) < 2L) return(character(0))
  labels  <- unlist(sets, use.names = FALSE)
  biggest <- sets[[which.max(lengths(sets))]]
  if (!anyDuplicated(labels) || !all(labels %in% biggest)) return(character(0))
  biggest[order(nchar(biggest), biggest)]
}

#' The warning for reused cage labels, or NULL
#' @noRd
reused_cage_message <- function(df, cage_column, treatment_column) {
  if (!has_col(df, cage_column) || !has_col(df, treatment_column)) return(NULL)
  lab <- reused_cage_labels(df[[cage_column]], df[[treatment_column]])
  if (!length(lab)) return(NULL)
  paste0("Every arm's cage labels come from one set (",
         paste(utils::head(lab, 5L), collapse = ", "), if (length(lab) > 5L) ", ...",
         "), so each label is read as one cage holding several arms. If cages are ",
         "numbered within each arm, make the labels unique first: label_cages(), or ",
         "'Numbered within each arm' under Study Design on the dashboard.")
}

#' Warn when the arms reuse cage labels
#' @noRd
warn_reused_cage_labels <- function(df, cage_column, treatment_column) {
  msg <- reused_cage_message(df, cage_column, treatment_column)
  if (!is.null(msg)) warning(msg, call. = FALSE)
  invisible(msg)
}
