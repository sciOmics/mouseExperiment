# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Global permutation test of the treatment x time interaction.
#
# CODE_REVIEW.md H.3. This was deliberately held open through Round 3 because a
# permutation test is a *randomisation* test: it is only valid if what you permute
# matches how the experiment actually randomised. Getting that wrong produces
# something that looks rigorous and is calibrated for the wrong design, which is
# worse than not having it.
#
# The maintainer's answer (2026-07-29) is that mice are most often randomised to
# treatment groups, though it varies. So the unit is an explicit argument rather
# than an assumption, `perm_spec(unit = "mouse")` is the default, and the
# cage-level alternative is available and validated against the design.
#
# Why this test is worth having at all: the lme4 path's p-values rest on the
# Kenward-Roger / Satterthwaite denominator-df approximation, which is known to be
# anti-conservative at the group sizes these studies use, and on normality of the
# random effects. A randomisation test assumes neither. It tests exactly the null
# the randomisation created.

#' Permutation-test configuration
#'
#' Bundles the settings for [trajectory_permutation_test()], following the same
#' config-helper pattern as [tg_priors()] and [tg_mcmc()] so entry-point
#' signatures stay readable (CODE_REVIEW.md Round 2 D.2).
#'
#' @param unit The unit of randomisation — **what was actually randomised to
#'   treatment**, not what is convenient to permute.
#'   \describe{
#'     \item{`"mouse"`}{(default) Individual animals were assigned to treatment.
#'       Whole mouse trajectories are relabelled. This is the common case.}
#'     \item{`"cage"`}{Cages were assigned to treatment, so every animal in a cage
#'       shares its label and moves with it. Requires cage to be nested within
#'       treatment; a crossed design cannot be permuted at cage level and is
#'       rejected.}
#'   }
#'   Getting this wrong invalidates the test. If mice were individually
#'   randomised and then re-housed by group, the unit is `"mouse"` and cage is a
#'   post-randomisation variable.
#' @param n_perm Number of label permutations. Each one refits the mixed model, so
#'   cost scales linearly; 999 is a reasonable default and 199 is enough for a
#'   quick look.
#' @param seed Optional integer seed. The caller's RNG state is restored on exit.
#' @return An object of class `perm_spec`.
#' @seealso [trajectory_permutation_test()], [tg_priors()], [tg_mcmc()]
#' @examples
#' perm_spec()
#' perm_spec(unit = "cage", n_perm = 199, seed = 42)
#' @export
perm_spec <- function(unit = c("mouse", "cage"), n_perm = 999L, seed = NULL) {
  unit <- match.arg(unit)
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) {
    stop("`n_perm` must be a positive integer.", call. = FALSE)
  }
  structure(list(unit = unit, n_perm = n_perm, seed = seed),
            class = "perm_spec")
}

#' @export
print.perm_spec <- function(x, ...) {
  cat("<perm_spec>\n")
  cat("  randomisation unit:", x$unit, "\n")
  cat("  permutations      :", x$n_perm, "\n")
  cat("  seed              :", if (is.null(x$seed)) "(caller's RNG)" else x$seed, "\n")
  invisible(x)
}

#' Enumerate all distinct permutations of a label multiset
#'
#' For a design small enough to enumerate, exhaustive enumeration is both exact
#' and cheaper than sampling: 6 cages in 3 arms of 2 has only 90 distinct
#' assignments, fewer than a typical `n_perm`. It also guarantees the reported
#' p-value respects the design's resolution floor, which a Monte-Carlo estimate
#' does not — a sampled p-value can land *below* the exact minimum simply because
#' too few of the sampled relabellings happened to tie with the observed one.
#'
#' @param labels Character vector of group labels (with repeats).
#' @return A list of character vectors, each a distinct permutation.
#' @noRd
#' @keywords internal
enumerate_label_perms <- function(labels) {
  out <- list()
  rec <- function(remaining, acc) {
    if (!length(remaining)) {
      out[[length(out) + 1L]] <<- acc
      return(invisible(NULL))
    }
    for (lv in sort(unique(remaining))) {
      i <- match(lv, remaining)
      rec(remaining[-i], c(acc, lv))
    }
    invisible(NULL)
  }
  rec(labels, character(0))
  out
}

#' Randomisation test for a treatment x time interaction
#'
#' Tests the global null that treatment has no effect on the growth trajectory, by
#' re-randomising treatment labels the way the experiment did and recomputing the
#' interaction statistic each time.
#'
#' @section Why use this:
#' The `lme4` path reports p-values that depend on the Kenward-Roger /
#' Satterthwaite denominator-df approximation and on normally-distributed random
#' effects. Both are approximations, and the df approximation is known to be
#' anti-conservative at the group sizes typical of these studies (8-10 animals per
#' arm). This test relies on neither: under the null created by the actual
#' randomisation, every relabelling is equally likely, so the reference
#' distribution is exact up to Monte-Carlo error in `n_perm`.
#'
#' @section What is permuted, and why it matters:
#' A permutation test is only valid if the permuted labels are exchangeable under
#' the null, which means the permutation has to mirror the randomisation. With
#' `unit = "mouse"` the treatment labels are shuffled across animals, holding each
#' animal's whole trajectory together. With `unit = "cage"` whole cages are
#' relabelled, so co-housed animals stay together.
#'
#' The choice is consequential. If cages were the unit of assignment, the number of
#' exchangeable units is the number of *cages*, not animals — a 6-arm study with 2
#' cages per arm has 12 units, not 48, and the test's resolution reflects that. The
#' returned `min_attainable_p` makes this explicit rather than leaving it to be
#' discovered.
#'
#' @section Test statistic:
#' The likelihood-ratio statistic comparing the full model
#' (`response ~ treatment * time + (1 | id)`) against the reduced model without the
#' interaction (`response ~ treatment + time + (1 | id)`), both fitted by ML. Large
#' values indicate the trajectories diverge. The statistic is only used to rank
#' permutations, so its own distributional properties do not matter — that is the
#' point of the approach.
#'
#' @param df Long data frame, one row per observation.
#' @param time_column,volume_column,treatment_column,id_column Column names.
#' @param cage_column Cage column. Required when `spec$unit == "cage"`.
#' @param transform Applied to the response before modelling: `"log"` (default),
#'   `"sqrt"`, or `"none"`. Should match what the primary analysis used.
#' @param spec A [perm_spec()] object.
#' @param verbose Print progress.
#' @return A list with:
#'   \item{p_value}{Two-sided permutation p-value, `(1 + #{T* >= T}) / (1 + n_perm)`
#'     so it is never exactly zero.}
#'   \item{observed_statistic}{The LRT statistic on the real labels.}
#'   \item{perm_statistics}{The `n_perm` null-distribution values.}
#'   \item{unit, n_units, n_perm}{What was permuted and how many of them.}
#'   \item{exhaustive}{`TRUE` when every distinct assignment was enumerated, so
#'     the p-value is exact rather than a Monte-Carlo estimate. Enumeration kicks
#'     in automatically for small designs, where it is both exact and usually
#'     cheaper than sampling.}
#'   \item{min_attainable_p}{The smallest p-value the design can produce
#'     regardless of effect size. This is **not** `1 / n_distinct_assignments`:
#'     the statistic is invariant to relabelling the groups, so with `g`
#'     equal-sized groups `g!` assignments tie with the observed one and the floor
#'     is `g! / n_distinct_assignments`. With 6 cages in 3 arms of 2 that is
#'     `6/90 = 0.067`, so p < 0.05 is unattainable — worth knowing before
#'     interpreting a null result as evidence of no effect.}
#'   \item{n_distinct_assignments}{Multinomial count of distinct labellings.}
#'   \item{n_distinct_statistics}{Distinct statistic values, i.e.
#'     `n_distinct_assignments` divided by the label symmetry.}
#'   \item{lrt_p_asymptotic}{The chi-square LRT p-value, for comparison. A large
#'     discrepancy against `p_value` is a signal that the asymptotic approximation
#'     is not to be trusted on this dataset.}
#' @seealso [perm_spec()], [tumor_growth_statistics()]
#' @examples
#' \dontrun{
#' data(master_synthetic_data)
#' trajectory_permutation_test(
#'   master_synthetic_data,
#'   spec = perm_spec(unit = "mouse", n_perm = 199, seed = 1)
#' )
#' }
#' @export
trajectory_permutation_test <- function(df,
                                        time_column      = "Day",
                                        volume_column    = "Volume",
                                        treatment_column = "Treatment",
                                        id_column        = "ID",
                                        cage_column      = "Cage",
                                        transform        = c("log", "sqrt", "none"),
                                        spec             = perm_spec(),
                                        verbose          = FALSE) {
  transform <- match.arg(transform)
  if (!inherits(spec, "perm_spec")) {
    stop("`spec` must come from perm_spec(). ",
         "See ?perm_spec for the randomisation-unit choice, which determines ",
         "whether this test is valid for your design.", call. = FALSE)
  }

  req <- c(time_column, volume_column, treatment_column, id_column)
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols)) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  has_cage <- !is.null(cage_column) && cage_column %in% names(df)
  if (spec$unit == "cage" && !has_cage) {
    stop("spec$unit = 'cage' needs a cage column, but '", cage_column,
         "' is not in the data.", call. = FALSE)
  }

  d <- data.frame(
    .y    = as.numeric(df[[volume_column]]),
    .t    = as.numeric(df[[time_column]]),
    .trt  = as.character(df[[treatment_column]]),
    .id   = as.character(df[[id_column]]),
    .cage = if (has_cage) as.character(df[[cage_column]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  d <- d[is.finite(d$.y) & is.finite(d$.t), , drop = FALSE]
  if (nrow(d) < 6L) stop("Too few usable observations.", call. = FALSE)

  if (transform == "log") {
    pos <- d$.y[d$.y > 0]
    if (!length(pos)) {
      stop("No positive values; cannot apply the log transform.", call. = FALSE)
    }
    d$.y[d$.y <= 0] <- min(pos) / 2
    d$.y <- log(d$.y)
  } else if (transform == "sqrt") {
    d$.y <- sqrt(pmax(d$.y, 0))
  }

  # Composite animal key so mice sharing an ID across arms/cages stay distinct.
  d$.mouse <- if (has_cage) {
    make_mouse_key(d$.trt, d$.id, d$.cage)
  } else {
    make_mouse_key(d$.trt, d$.id)
  }

  # ── Exchangeable units and their fixed labels ───────────────────────────────
  if (spec$unit == "cage") {
    cage_info <- classify_cage_structure(df, cage_column, id_column,
                                         treatment_column)
    if (identical(cage_info$structure, "crossed")) {
      stop("spec$unit = 'cage' but the design is crossed: ",
           cage_info$description,
           "\nA cage holding more than one treatment has no single label to ",
           "permute. If animals were randomised individually and then re-housed, ",
           "the randomisation unit is the mouse — use perm_spec(unit = 'mouse').",
           call. = FALSE)
    }
    unit_of <- d$.cage
  } else {
    unit_of <- d$.mouse
  }

  unit_tbl <- unique(data.frame(unit = unit_of, trt = d$.trt,
                                stringsAsFactors = FALSE))
  dup <- unit_tbl$unit[duplicated(unit_tbl$unit)]
  if (length(dup)) {
    stop("Unit '", dup[1], "' carries more than one treatment label, so it ",
         "cannot be relabelled as a block. ", call. = FALSE)
  }
  unit_labels <- stats::setNames(unit_tbl$trt, unit_tbl$unit)
  n_units <- length(unit_labels)
  if (n_units < 3L) {
    stop("Only ", n_units, " exchangeable ", spec$unit,
         "(s) — too few to permute.", call. = FALSE)
  }

  # Resolution ceiling. Two steps, and the second is easy to get wrong.
  #
  # (1) Distinct labellings = multinomial coefficient over the group sizes.
  # (2) The interaction statistic depends only on the *partition* of units into
  #     groups, not on which label each group receives. So any relabelling that
  #     maps the multiset of group sizes onto itself yields an IDENTICAL
  #     statistic, and those permutations tie with the observed one. With g
  #     equal-sized groups that is g! ties, so the floor is g! / n_assign, not
  #     1 / n_assign.
  #
  # Getting this wrong understates the floor by a factor of g!. It was caught by
  # a null-calibration run: at cage level with 6 cages in 3 arms of 2 the test
  # never produced p < 0.05 in 60 null datasets, because the true floor is
  # 6/90 = 0.067 while the naive 1/90 = 0.011 suggested p < 0.05 was reachable.
  grp_sizes <- as.integer(table(unit_labels))
  log_n_assign <- lgamma(n_units + 1) - sum(lgamma(grp_sizes + 1))
  n_assign <- exp(log_n_assign)
  # Label permutations that preserve the partition: product of factorials of the
  # counts of equal group sizes.
  n_sym <- prod(factorial(as.integer(table(grp_sizes))))
  n_distinct_stats <- n_assign / n_sym
  min_p <- n_sym / n_assign

  if (min_p > 0.01) {
    warning(sprintf(paste0(
      "Only %.0f distinct label assignments across %d %s(s), and the statistic ",
      "is invariant to relabelling the %d groups, leaving %.0f distinct ",
      "statistic values. The smallest attainable p-value is therefore %.3f -- ",
      "a result at that floor reflects the design's resolution, not the ",
      "strength of evidence.%s"),
      n_assign, n_units, spec$unit, length(grp_sizes), n_distinct_stats, min_p,
      if (min_p > 0.05) paste0(" Note p < 0.05 is UNATTAINABLE with ",
                               n_units, " ", spec$unit, "s at these group sizes.")
      else ""), call. = FALSE)
  }

  # ── Test statistic: ML likelihood-ratio for the interaction ─────────────────
  ctrl <- lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                            check.nobs.vs.nRE = "ignore",
                            calc.derivs = FALSE)
  lrt_stat <- function(dat) {
    fit_full <- tryCatch(
      suppressMessages(suppressWarnings(lme4::lmer(
        .y ~ .trt * .t + (1 | .mouse), data = dat, REML = FALSE, control = ctrl))),
      error = function(e) NULL)
    fit_red <- tryCatch(
      suppressMessages(suppressWarnings(lme4::lmer(
        .y ~ .trt + .t + (1 | .mouse), data = dat, REML = FALSE, control = ctrl))),
      error = function(e) NULL)
    if (is.null(fit_full) || is.null(fit_red)) return(NA_real_)
    as.numeric(2 * (stats::logLik(fit_full) - stats::logLik(fit_red)))
  }

  observed <- lrt_stat(d)
  if (!is.finite(observed)) {
    stop("The full or reduced model could not be fitted on the observed data, ",
         "so there is no statistic to permute.", call. = FALSE)
  }
  df_lrt <- (length(grp_sizes) - 1L)      # interaction df for a linear time term
  lrt_p_asym <- stats::pchisq(observed, df = df_lrt, lower.tail = FALSE)

  # ── Permute ────────────────────────────────────────────────────────────────
  if (!is.null(spec$seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(spec$seed)
  }

  unit_names <- names(unit_labels)
  base_labels <- unname(unit_labels)

  # Enumerate exhaustively when the design is small enough. Besides being exact,
  # this is often *cheaper* than sampling (90 assignments < a 999-permutation
  # run) and it guarantees the p-value cannot fall below `min_p` — a Monte-Carlo
  # estimate can, purely because too few sampled relabellings happened to tie
  # with the observed one, which is how the floor bug above was found.
  ENUM_LIMIT <- 2000
  exhaustive <- n_assign <= max(ENUM_LIMIT, spec$n_perm)

  assignments <- if (exhaustive) {
    enumerate_label_perms(base_labels)
  } else {
    lapply(seq_len(spec$n_perm), function(b) sample(base_labels))
  }
  n_draw <- length(assignments)

  if (isTRUE(verbose)) {
    message(if (exhaustive) "Enumerating all " else "Sampling ", n_draw,
            if (exhaustive) " distinct assignments" else " permutations",
            " at the ", spec$unit, " level (", n_units, " units)...")
  }

  perm_stats <- rep(NA_real_, n_draw)
  for (b in seq_len(n_draw)) {
    shuffled <- stats::setNames(assignments[[b]], unit_names)
    dp <- d
    dp$.trt <- unname(shuffled[unit_of])
    perm_stats[b] <- lrt_stat(dp)
    if (isTRUE(verbose) && b %% 100L == 0L) message("  ", b, "/", n_draw)
  }

  ok <- is.finite(perm_stats)
  n_ok <- sum(ok)
  if (n_ok < 2L) {
    stop("Almost no permutation refits converged (", n_ok, " of ", n_draw,
         "), so no null distribution could be built.", call. = FALSE)
  }
  if (n_ok < n_draw) {
    warning(n_draw - n_ok, " of ", n_draw,
            " permutation refits failed to converge and were dropped.",
            call. = FALSE)
  }

  p_value <- if (exhaustive) {
    # Exact: the observed labelling is itself one of the enumerated assignments,
    # so it is counted and p >= min_p holds by construction.
    sum(perm_stats[ok] >= observed - 1e-9) / n_ok
  } else {
    # (1 + count) / (1 + n) so a sampled p-value is never exactly zero.
    (1 + sum(perm_stats[ok] >= observed - 1e-9)) / (1 + n_ok)
  }

  list(
    p_value                = p_value,
    observed_statistic     = observed,
    perm_statistics        = perm_stats[ok],
    unit                   = spec$unit,
    n_units                = n_units,
    group_sizes            = stats::setNames(grp_sizes, names(table(unit_labels))),
    n_perm                 = n_draw,
    n_perm_requested       = spec$n_perm,
    n_perm_converged       = n_ok,
    exhaustive             = exhaustive,
    n_distinct_assignments = n_assign,
    n_distinct_statistics  = n_distinct_stats,
    min_attainable_p       = min_p,
    lrt_p_asymptotic       = lrt_p_asym,
    transform_used         = transform,
    method = paste0(
      "Randomisation test of the treatment x time interaction; ML ",
      "likelihood-ratio statistic, ",
      if (exhaustive) paste0("exhaustive enumeration of all ", n_draw,
                             " distinct assignments")
      else paste0(n_draw, " sampled permutations"),
      " of treatment labels at the ", spec$unit, " level")
  )
}
