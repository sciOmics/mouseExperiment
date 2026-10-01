# The whole dashboard, as runDashboard() serves it. A start-up error in its
# server function closes every session while the page itself still loads, which
# only a real browser session shows (CODE_REVIEW.md R20.D3).
source(file.path("..", "load-package.R"))
asNamespace("mouseExperimentDashboard")$dashboard_app()
