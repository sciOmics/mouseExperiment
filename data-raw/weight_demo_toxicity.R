# The Weight demo's drug arms: their body weights (CODE_REVIEW.md R23.24).
# Run from the package root, then data-raw/create_datasets.R:
#   Rscript data-raw/weight_demo_toxicity.R
#
# No animal crossed 20 % weight loss: the combination's worst was -19.7 %, and
# Drug B lost a moderate 12.6 %. Now the combination is severely toxic and the
# single agents mild, so the threshold analysis, its Kaplan-Meier curve and
# the therapeutic window have something to show:
#   * Drug A + B: three of five animals pass -20 %. Two pass it on March 11
#     and leave the study there, as a humane endpoint requires (their March 15
#     rows are dropped); the third passes it on March 15. The other two lose
#     13 % and 16 %.
#   * Drug B: about 7.5 % by March 15. Drug A: about 3 %.
#   * DMSO is unchanged.
# Each weight is the animal's first-day weight times the planned change,
# rounded to 0.1 g, so the file is reproducible. The tumour volumes are not
# changed.

path <- "data-raw/weight_synthetic_data.csv"
# As text, so that unchanged values are written exactly as they were.
d <- utils::read.csv(path, check.names = FALSE, colClasses = "character")
dates <- c("March 1", "March 4", "March 8", "March 11", "March 15")
stopifnot(all(d$Date %in% dates))

# Percentage change from the first weighing, on each date.
plan <- list(
  "Drug A" = list(N = c(0, -0.4, -1.3, -2.4, -2.9), L = c(0, -0.6, -1.6, -2.6, -3.2),
                  R = c(0, -0.3, -1.2, -2.2, -2.7), LR = c(0, -0.7, -1.8, -2.8, -3.4),
                  LL = c(0, -0.5, -1.5, -2.5, -3.0)),
  "Drug B" = list(N = c(0, -1.4, -3.8, -5.8, -7.2), L = c(0, -1.7, -4.3, -6.4, -7.9),
                  R = c(0, -1.2, -3.6, -5.6, -7.0), LR = c(0, -1.9, -4.5, -6.6, -8.1),
                  LL = c(0, -1.5, -4.0, -6.0, -7.5)),
  "Drug A + B" = list(N = c(0, -6, -13, -21, NA), R = c(0, -7, -14, -22, NA),
                      L = c(0, -5, -11, -16, -22), LR = c(0, -4, -8, -11, -13),
                      LL = c(0, -5, -9, -13, -16)))

for (arm in names(plan)) for (tag in names(plan[[arm]])) {
  rows <- which(d$Treatment == arm & d[["Ear Tag"]] == tag)
  base <- as.numeric(d[["Mouse Weight"]][rows[d$Date[rows] == "March 1"]])
  pct <- plan[[arm]][[tag]][match(d$Date[rows], dates)]
  w <- round(base * (1 + pct / 100), 1)
  d[["Mouse Weight"]][rows] <- ifelse(is.na(w), NA_character_, sprintf("%.1f", w))
}
# Past -20 % the animal leaves the study: its later rows go.
d <- d[!is.na(d[["Mouse Weight"]]), , drop = FALSE]

utils::write.csv(d, path, row.names = FALSE, quote = FALSE)   # as the file was
file.copy(path, "inst/sample_data/weight_synthetic_data.csv", overwrite = TRUE)
cat(nrow(d), "rows written to", path, "and inst/sample_data/\n")
