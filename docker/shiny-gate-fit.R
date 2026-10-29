# Image gate: a Bayesian tumour-growth fit at the dashboard's defaults, compiled
# and sampled inside production's memory budget for a fit (CODE_REVIEW.md
# R23-A).
#
# On 2026-10-25 production's compiler was OOM-killed: a newer brms and Stan
# made the dashboard's model need more memory to compile than the container
# had left beside the app (sciOmicsVPS CODE_REVIEW.md V3-T). docker/shiny-gate.R
# compiles a one-line model, which never came close, and the staging server
# that ran the suite skipped Bayesian fits. CI runs this script in a container
# whose memory limit is that budget, with no swap, so a model that would no
# longer compile on the server fails here first.
#
# The call is the one the dashboard's worker makes (tg_resume.R) at the
# dashboard's defaults: 4 chains of 1,000 warmup and 500 draws, rstan. It
# prints the container's peak memory (cgroup v2 memory.peak) and the time.

options(warn = 1)
t0 <- Sys.time()
step <- function(...) cat(sprintf("[fit %6.1fs] ", as.numeric(Sys.time() - t0, units = "secs")), ..., "\n", sep = "")

suppressPackageStartupMessages(library(mouseExperiment))
step("mouseExperiment ", as.character(utils::packageVersion("mouseExperiment")),
     "; R_MAKEVARS_USER: ", Sys.getenv("R_MAKEVARS_USER", "(unset)"))

fit <- suppressMessages(bayesian_tumor_growth(
  master_synthetic_data, reference_group = "Vehicle", cage_column = "Cage",
  n_chains = 4, n_warmup = 1000, n_iter = 500, backend = "rstan",
  return_model = FALSE, plots = FALSE, verbose = FALSE))
if (is.null(fit$posterior_summary) || !nrow(fit$posterior_summary))
  stop("the fit returned no posterior summary")
step("fitted; the model was ", fit$summary$methods$compiled_model)

peak <- tryCatch(as.numeric(readLines("/sys/fs/cgroup/memory.peak", warn = FALSE)) / 2^20,
                 error = function(e) NA_real_)
step("container peak memory: ", if (is.na(peak)) "not available" else sprintf("%.0f MB", peak))

cat("FIT_OK\n")
