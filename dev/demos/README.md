# Personal demo notebooks (not shipped)

These Quarto notebooks were moved out of `vignettes/` in v0.29.0
(CODE_REVIEW.md R20.49). They are not built with the package, and they do not run
as they stand:

- `mouseExperiment_dose_demo.qmd` installs the package from a local iCloud path,
  and reads `./data/*.csv` files that do not exist.
- `mouseExperiment_combo_demo.qmd` installs from a different GitHub organization
  (`SciOmicsLab`). It also refers to a local logo, and presents the Combination
  Index result that was withdrawn in R17.2.

The package vignette (`vignettes/mouseExperiment.Rmd`) covers the same workflow
on the shipped `master_synthetic_data`.
