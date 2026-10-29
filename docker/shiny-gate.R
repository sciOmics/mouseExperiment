# Gate for a candidate server image (sciOmicsVPS CODE_REVIEW.md V3.9, V3.10).
#
# CI (.github/workflows/image.yaml) and the server's rebuild-shiny.sh pipe this
# into `R --no-save` inside a throwaway container of the candidate image, as
# the image's own user and without a network. Any error exits non-zero and the
# candidate is never deployed. It lived in sciOmicsVPS until 0.49.2 (R23-A).
#
# The old gate ran library(mouseExperimentDashboard). That proves every hard
# Import is installed, but it loads only the namespaces the NAMESPACE imports,
# so brms, rstan, lme4, lmerTest, coxphf, clinfun, coin, car, future and
# promises were never loaded, and no Stan model was ever compiled.
#
# Since mouseExperiment 0.34.0 the dashboard is part of the package
# (run_dashboard()), and the image installs that one package (V3-S).
#
# Environment (optional): EXPECT_PKG_SHA -- the commit the image must have been
# installed from.

options(warn = 1)
t0 <- Sys.time()
step <- function(...) cat(sprintf("[gate %5.1fs] ", as.numeric(Sys.time() - t0, units = "secs")), ..., "\n", sep = "")

# 1. The package comes from the commit that was asked for --------------------
got <- utils::packageDescription("mouseExperiment")[["RemoteSha"]]
want <- Sys.getenv("EXPECT_PKG_SHA")
step("mouseExperiment ", as.character(utils::packageVersion("mouseExperiment")), " @ ", got)
if (nzchar(want) && !identical(got, want))
  stop("mouseExperiment was installed from ", got, ", not ", want)
ns <- asNamespace("mouseExperiment")
if (!exists("run_dashboard", envir = ns, inherits = FALSE))
  stop("this mouseExperiment has no dashboard (run_dashboard() arrived in 0.34.0)")

# 2. Every Import, and every package the dashboard needs, loads --------------
imports_of <- function(pkg) {
  d <- utils::packageDescription(pkg)
  f <- paste(c(d[["Depends"]], d[["Imports"]]), collapse = ",")
  x <- trimws(sub("\\(.*", "", strsplit(f, ",", fixed = TRUE)[[1]]))
  setdiff(x[nzchar(x)], "R")
}
# The dashboard's packages are Suggests; dashboard_packages() lists them.
pkgs <- unique(c("mouseExperiment", imports_of("mouseExperiment"),
                 get("dashboard_packages", envir = ns)()))
failed <- character(0)
for (p in pkgs) {
  ok <- tryCatch({ loadNamespace(p); TRUE }, error = function(e) {
    message("  cannot load ", p, ": ", conditionMessage(e))
    FALSE
  })
  if (!ok) failed <- c(failed, p)
}
if (length(failed)) stop("namespaces that do not load: ", paste(failed, collapse = ", "))
step("loaded ", length(pkgs), " namespaces")

# 3. rstan compiles and samples (R23.3: the dashboard's only backend; the
#    image has no CmdStan since 0.50.0) ----------------------------------------
code <- "parameters { real y; } model { y ~ normal(0, 1); }"
m <- rstan::stan_model(model_code = code, verbose = FALSE)
fit <- rstan::sampling(m, chains = 1, iter = 200, refresh = 0, seed = 1)
stopifnot(length(rstan::extract(fit)$y) == 100)
step("rstan compiled and sampled")

# 4. The dashboard's UI builds (its www/ check and UI construction used to run
#    only when the container booted) -----------------------------------------
app <- get("dashboard_app", envir = ns)()
stopifnot(inherits(app, "shiny.appobj"))
step("dashboard app object built")

cat("GATE_OK\n")
