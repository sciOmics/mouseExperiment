# =============================================================================
# helpers_plot.R / helpers_ui.R — CODE_REVIEW.md D.1 / D.9
#
# Lower stakes than the data and format helpers — a wrong colour misleads nobody
# — with two exceptions that are not cosmetic:
#
#   * `remap_*_to_treatment()` re-keys aesthetics onto treatment names, so a
#     mismatch mislabels which arm a curve belongs to.
#   * `get_ordered_groups()` feeds the reference-group choice, which is a
#     statistical input.
#
# On the UI side, `translate_error()` decides what a user is told when an
# analysis fails. A translation that swallows the cause turns a diagnosable
# error into a shrug, and it had no test.
# =============================================================================

# ---- colour and shape vectors ----------------------------------------------

test_that("D.9: build_color_vector returns NULL on partial coverage", {
  # The comment in the source explains why this matters: ggplot2's
  # scale_colour_manual maps unlisted groups to grey, so a partially-filled
  # vector silently renders an arm invisible rather than failing. All-or-nothing
  # is the correct contract and it was untested.
  groups <- c("Control", "DrugA")
  full <- list(tg_Control = "#111111", tg_DrugA = "#222222")
  expect_identical(unname(build_color_vector(full, "tg", groups)),
                   c("#111111", "#222222"))

  partial <- list(tg_Control = "#111111")            # DrugA missing
  expect_null(build_color_vector(partial, "tg", groups))

  blank <- list(tg_Control = "#111111", tg_DrugA = "")  # present but empty
  expect_null(build_color_vector(blank, "tg", groups))
})

test_that("D.9: build_color_vector handles names needing make.names()", {
  # Input IDs are built with make.names(), so a group called "Drug A (10 mg)"
  # resolves to a mangled id. If the two ever disagree the lookup silently
  # misses and the whole vector collapses to NULL.
  groups <- c("Drug A (10 mg)")
  vals <- stats::setNames(list("#333333"), paste0("tg_", make.names(groups)))
  expect_identical(unname(build_color_vector(vals, "tg", groups)), "#333333")
})

test_that("D.9: build_color_vector is NULL-safe on empty groups", {
  expect_null(build_color_vector(list(), "tg", NULL))
  expect_null(build_color_vector(list(), "tg", character(0)))
})

test_that("D.9: remapping to treatment collapses dose levels, keeping the first", {
  # Used when a dose column is mapped: "DrugA - Dose: 10" and "DrugA - Dose: 30"
  # both belong to treatment "DrugA 10"/"DrugA 30". The de-duplication is what
  # stops two entries with the same name reaching scale_*_manual, which errors.
  cols <- c("DrugA - Dose: 10" = "#AAAAAA", "DrugA - Dose: 30" = "#BBBBBB")
  out <- remap_colors_to_treatment(cols)
  expect_identical(names(out), c("DrugA 10", "DrugA 30"))
  expect_identical(unname(out), c("#AAAAAA", "#BBBBBB"))

  # Genuine duplicates after remapping keep the first occurrence.
  dup <- c("X - Dose: 1" = "#111111", "X - Dose: 1" = "#222222")
  expect_length(remap_colors_to_treatment(dup), 1L)
  expect_identical(unname(remap_colors_to_treatment(dup)), "#111111")

  expect_null(remap_colors_to_treatment(NULL))
})

test_that("D.9: shape remapping behaves identically to colour remapping", {
  # They are separate functions with the same body; if one is ever changed the
  # other must be too, or shapes and colours key differently and the legend
  # disagrees with itself.
  x <- c("DrugA - Dose: 10" = 16, "DrugB - Dose: 5" = 17)
  out <- remap_shapes_to_treatment(x)
  expect_identical(names(out), c("DrugA 10", "DrugB 5"))
  expect_null(remap_shapes_to_treatment(NULL))
})

# ---- palettes --------------------------------------------------------------

test_that("D.9: every palette returns exactly n colours", {
  # Recycling matters: okabe_ito has 8 entries and wes_anderson 5, so a study
  # with 10 arms must still get 10 colours rather than a short vector that
  # silently drops arms.
  for (p in c("viridis", "okabe_ito", "set2", "wes_anderson", "unknown_name")) {
    for (n in c(1L, 5L, 10L)) {
      expect_length(palette_colors(p, n), n)
    }
  }
})

test_that("D.9: palettes return valid hex and degrade to zero cleanly", {
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", palette_colors("viridis", 4))))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", palette_colors("okabe_ito", 4))))
  for (p in c("viridis", "okabe_ito", "unknown_name")) {
    expect_length(palette_colors(p, 0L), 0L)
  }
  expect_length(default_group_colors(0L), 0L)
  expect_length(default_group_colors(3L), 3L)
})

# ---- group ordering --------------------------------------------------------

test_that("D.9: a user reordering is honoured only when it is a permutation", {
  groups <- c("Control", "DrugA", "DrugB")

  # Valid permutation: applied.
  expect_identical(
    get_ordered_groups(list(tg_order = c("DrugB", "Control", "DrugA")), "tg", groups),
    c("DrugB", "Control", "DrugA"))

  # Wrong length, or different membership: ignored rather than partially
  # applied. This is the guard that stops a stale ordering — left over from a
  # previous dataset — from dropping an arm out of the plot entirely.
  expect_identical(get_ordered_groups(list(tg_order = c("DrugB", "Control")),
                                      "tg", groups), groups)
  expect_identical(get_ordered_groups(list(tg_order = c("X", "Y", "Z")),
                                      "tg", groups), groups)
  expect_identical(get_ordered_groups(list(), "tg", groups), groups)
  expect_null(get_ordered_groups(list(), "tg", NULL))
})

# ---- error translation -----------------------------------------------------

test_that("D.9: known failure modes are translated into actionable advice", {
  expect_match(translate_error("object 'Volume' not found"), "column")
  expect_match(translate_error("boundary (singular) fit: see help('isSingular')"),
               "singular")
  expect_match(translate_error("subscript out of bounds"), "indexing")
  expect_match(translate_error("non-conformable arguments"), "dimension")
  expect_match(translate_error("NaN produced"), "missing or invalid")
})

test_that("D.9: an unrecognised error is passed through, not swallowed", {
  # The important half. A translator that returned a generic message for
  # anything it did not recognise would destroy the only diagnostic the user
  # has. Unknown messages must survive verbatim.
  msg <- "Stan model 'x' does not contain samples after warmup"
  expect_identical(translate_error(msg), msg)
})

test_that("D.9: the singular-fit branch wins over the missing-values branch", {
  # lme4 emits "boundary (singular) fit" and the message contains no NA text, but
  # ordering between overlapping patterns is worth pinning: the advice differs
  # (simplify random effects vs check for NAs) and only one is right.
  out <- translate_error("boundary (singular) fit")
  expect_match(out, "random effects")
  expect_false(grepl("missing or invalid", out))
})

# ---- status alerts and JS escaping -----------------------------------------

test_that("D.9: status_alert maps each type to its own bootstrap class", {
  for (ty in c("danger", "warning", "info", "success")) {
    html <- as.character(status_alert("msg", ty))
    expect_match(html, paste0("alert-", ty))
    expect_match(html, "msg")
  }
  # An unknown type still renders rather than erroring.
  expect_match(as.character(status_alert("msg", "nonsense")), "alert-nonsense")
})

test_that("D.9: js_string escapes values interpolated into runjs()", {
  # These values reach the browser inside a JS string via runjs(sprintf(...)).
  # Quotes and backslashes must not be able to terminate the literal early.
  expect_match(as.character(js_string("plain")), '^"plain"$')

  quoted <- as.character(js_string('say "hi"'))
  expect_match(quoted, '\\\\"', info = "embedded double quotes must be escaped")

  backslash <- as.character(js_string("a\\b"))
  expect_match(backslash, "\\\\\\\\", info = "backslashes must be escaped")

  # A closing-tag payload must not survive unescaped into the page.
  expect_false(grepl('</script>', as.character(js_string("</script>")),
                     fixed = TRUE))

  # NULL and empty become an empty JS string rather than the literal "NULL".
  expect_identical(as.character(js_string(NULL)), '""')
  expect_identical(as.character(js_string("")), '""')
})
