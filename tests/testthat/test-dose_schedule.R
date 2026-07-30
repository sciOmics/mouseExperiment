# =============================================================================
# Dosing-schedule annotation — CODE_REVIEW.md D.15
#
# Replaces plot_treatments(), which was removed along with the two schedule
# datasets. The removal is not a loss of capability: the schedule is now drawn
# into the growth/weight/survival plots themselves rather than as a second panel
# whose x-axis had to be kept aligned by hand.
#
# The constraint driving the design is worth restating, because it is what makes
# the text-field input the right call rather than a shortcut: dosing timing is
# not recoverable from measurement data. `Dose` is constant within an animal in
# every dataset this package ships — a level, not a schedule — and the dosing
# days are not on the measurement grid.
# =============================================================================

suppressMessages(library(ggplot2))

tg_df <- function(seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(c("Control", "DrugA"), function(tx) {
    d <- rep(seq(0, 21, 3), 4)
    data.frame(Day = d, Treatment = tx,
               Volume = exp(log(120) + 0.1 * d + stats::rnorm(length(d), 0, 0.1)),
               stringsAsFactors = FALSE)
  }))
}
base_plot <- function(df = tg_df()) {
  ggplot(df, aes(Day, Volume, colour = Treatment)) + geom_point() + theme_classic()
}
# Count the annotation rows actually drawn (one geom_point layer of shape 25 each).
mark_rows <- function(p) {
  g <- ggplot_build(p)
  sum(vapply(g$data, function(L)
    if (!is.null(L$shape) && any(L$shape == 25)) 1L else 0L, integer(1)))
}
mark_positions <- function(p) {
  g <- ggplot_build(p)
  ys <- c()
  for (L in g$data) if (!is.null(L$shape) && any(L$shape == 25)) ys <- c(ys, unique(L$y))
  tr <- g$layout$panel_scales_y[[1]]$trans
  if (is.null(tr)) tr <- g$layout$panel_scales_y[[1]]$transformation
  inv <- if (!is.null(tr$inverse)) tr$inverse else identity
  inv(ys)
}

# ---- parsing ---------------------------------------------------------------

test_that("D.15: the bare form applies one schedule to every group", {
  s <- parse_dose_schedule("1, 5, 9, 13", groups = c("Control", "DrugA"))
  expect_s3_class(s, "data.frame")
  expect_identical(attr(s, "scope"), "shared")
  expect_setequal(unique(s$Treatment), c("Control", "DrugA"))
  expect_setequal(unique(s$Day), c(1, 5, 9, 13))
  expect_identical(nrow(s), 8L)
})

test_that("D.15: the per-group form keeps arms distinct", {
  s <- parse_dose_schedule("DrugA: 1,5,9\nCombo: 1,5", groups = c("DrugA", "Combo"))
  expect_identical(attr(s, "scope"), "per_group")
  expect_setequal(s$Day[s$Treatment == "DrugA"], c(1, 5, 9))
  expect_setequal(s$Day[s$Treatment == "Combo"], c(1, 5))
})

test_that("D.15: separators are flexible and days are de-duplicated and sorted", {
  # Users will type whatever they type; the parser should not be the obstacle.
  for (txt in c("1,5,9", "1 5 9", "1, 5,  9", "9,5,1,5")) {
    s <- parse_dose_schedule(txt, groups = "A")
    expect_identical(s$Day, c(1, 5, 9), info = txt)
  }
  # Semicolons separate per-group entries as well as newlines.
  s <- parse_dose_schedule("A: 1,3; B: 2,4", groups = c("A", "B"))
  expect_identical(attr(s, "scope"), "per_group")
  expect_setequal(s$Day[s$Treatment == "B"], c(2, 4))
})

test_that("D.15: group names match case-insensitively but return the data's casing", {
  s <- parse_dose_schedule("druga: 1,5", groups = c("DrugA"))
  expect_identical(unique(s$Treatment), "DrugA")
})

test_that("D.15: empty input is not an error", {
  # The toggle can be on with the field still blank; that must render the plot
  # unchanged rather than fail.
  expect_null(parse_dose_schedule(""))
  expect_null(parse_dose_schedule("   "))
  expect_null(parse_dose_schedule(NULL))
  expect_null(parse_dose_schedule(character(0)))
})

test_that("D.15: unmatched group names warn rather than vanish", {
  # §K.2 class: a name that matches nothing must not fail silently, or a user
  # who mistypes an arm sees an annotation quietly missing that arm.
  expect_warning(
    s <- parse_dose_schedule("Nope: 1,5\nDrugA: 1,5", groups = c("DrugA")),
    "not found in the data")
  expect_setequal(unique(s$Treatment), "DrugA")
})

test_that("D.15: unreadable days warn rather than producing an empty annotation", {
  expect_warning(expect_null(parse_dose_schedule("every other day", groups = "A")),
                 "No dosing days")
  expect_warning(parse_dose_schedule("A: 1,5\nnonsense", groups = "A"),
                 "without a group label")
})

# ---- style selection -------------------------------------------------------

test_that("D.15: intermittent schedules get marks, dense ones get a window", {
  g <- "A"
  expect_identical(dose_schedule_style(parse_dose_schedule("1,5,9,13", groups = g)), "rug")
  expect_identical(dose_schedule_style(parse_dose_schedule("1,8,15", groups = g)), "rug")
  # 21 daily doses as 21 triangles is noise.
  expect_identical(dose_schedule_style(
    parse_dose_schedule(paste(0:20, collapse = ","), groups = g)), "band")
  # Daily but short: still continuous dosing, so still a window.
  expect_identical(dose_schedule_style(parse_dose_schedule("1,2,3,4", groups = g)), "band")
  expect_identical(dose_schedule_style(NULL), "rug")
})

# ---- annotation ------------------------------------------------------------

test_that("D.15: a NULL or empty schedule returns the plot untouched", {
  p <- base_plot()
  expect_identical(length(annotate_dose_schedule(p, NULL)$layers), length(p$layers))
  expect_identical(
    length(annotate_dose_schedule(p, data.frame(Treatment = character(0),
                                                Day = numeric(0)))$layers),
    length(p$layers))
})

test_that("D.15: marks are drawn below every data point on linear, log and sqrt scales", {
  # The whole design rests on the strip sitting outside the data region. Placing
  # it requires inverting the y-scale transform, which is the part most likely to
  # be silently wrong — a mark landing among the curves is the failure mode.
  df <- tg_df()
  s  <- parse_dose_schedule("1,5,9,13", groups = unique(df$Treatment))
  for (sc in list(NULL, scale_y_log10(), scale_y_sqrt())) {
    p <- base_plot(df)
    if (!is.null(sc)) p <- p + sc
    pos <- mark_positions(annotate_dose_schedule(p, s))
    expect_length(pos, 1L)
    expect_lt(max(pos), min(df$Volume))
  }
})

test_that("D.15: a shared schedule draws one row, not one per group", {
  # The bare form expands to every group. Drawing a row each would imply a
  # per-arm difference that does not exist.
  df <- tg_df()
  p <- annotate_dose_schedule(base_plot(df),
         parse_dose_schedule("1,5,9,13", groups = unique(df$Treatment)))
  expect_identical(mark_rows(p), 1L)
})

test_that("D.15: differing schedules stack one row per group", {
  df <- tg_df()
  p <- annotate_dose_schedule(base_plot(df),
         parse_dose_schedule("Control: 1,5\nDrugA: 1,5,9,13",
                             groups = unique(df$Treatment)))
  expect_identical(mark_rows(p), 2L)
  # Rows must not overlap, or the two arms are indistinguishable.
  expect_identical(length(unique(mark_positions(p))), 2L)
})

test_that("D.15: per-group schedules that coincide collapse to one row", {
  df <- tg_df()
  p <- annotate_dose_schedule(base_plot(df),
         parse_dose_schedule("Control: 1,5\nDrugA: 1,5", groups = unique(df$Treatment)))
  expect_identical(mark_rows(p), 1L)
})

test_that("D.15: supplied colours are used so marks match their curves", {
  df <- tg_df()
  p <- annotate_dose_schedule(
    base_plot(df),
    parse_dose_schedule("Control: 1,5\nDrugA: 1,5,9,13", groups = unique(df$Treatment)),
    colors = c(Control = "#123456", DrugA = "#654321"))
  g <- ggplot_build(p)
  fills <- unlist(lapply(g$data, function(L)
    if (!is.null(L$shape) && any(L$shape == 25)) unique(L$fill) else NULL))
  expect_setequal(fills, c("#123456", "#654321"))
})

test_that("D.15: an unknown group falls back to grey rather than erroring", {
  df <- tg_df()
  p <- annotate_dose_schedule(
    base_plot(df),
    parse_dose_schedule("Control: 1,5\nDrugA: 1,9", groups = unique(df$Treatment)),
    colors = c(Control = "#123456"))
  expect_s3_class(ggplot_build(p), "ggplot_built")
})

test_that("D.15: band style draws rectangles instead of marks", {
  df <- tg_df()
  s <- parse_dose_schedule(paste(0:20, collapse = ","), groups = unique(df$Treatment))
  p <- annotate_dose_schedule(base_plot(df), s)   # auto -> band
  expect_identical(mark_rows(p), 0L)
  g <- ggplot_build(p)
  rects <- sum(vapply(g$data, function(L)
    if (all(c("xmin", "xmax", "ymin", "ymax") %in% names(L))) 1L else 0L, integer(1)))
  expect_gt(rects, 0L)
})

test_that("D.15: style can be forced against the automatic choice", {
  df <- tg_df()
  s <- parse_dose_schedule(paste(0:20, collapse = ","), groups = unique(df$Treatment))
  expect_gt(mark_rows(annotate_dose_schedule(base_plot(df), s, style = "rug")), 0L)
  expect_identical(mark_rows(annotate_dose_schedule(base_plot(df), s, style = "band")), 0L)
})

test_that("D.15: a non-ggplot input is returned unchanged", {
  expect_identical(annotate_dose_schedule("not a plot",
                                          parse_dose_schedule("1,5", groups = "A")),
                   "not a plot")
})

test_that("D.15: the annotated plot still renders end to end", {
  # ggplot errors are deferred to draw time, so building is not proof.
  df <- tg_df()
  p <- annotate_dose_schedule(base_plot(df),
         parse_dose_schedule("Control: 1,5\nDrugA: 1,5,9,13",
                             groups = unique(df$Treatment)))
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  expect_silent(suppressMessages(
    ggplot2::ggsave(f, p, width = 6, height = 4, dpi = 72)))
  expect_true(file.exists(f) && file.size(f) > 0)
})
