# =============================================================================
# Shinytest2 apps
#
# The apps driven by test-shiny_app.R are app directories under
# tests/testthat/apps/, each with its own app.R that loads the package (see
# apps/load-package.R). They replaced the shinyApp-object builders that lived
# here: AppDriver serialises an app object into a fresh R process, which cannot
# see a package loaded with devtools::load_all(), and the builders still used
# the pre-moduleServer signatures (CODE_REVIEW.md R20.D11).
# =============================================================================
