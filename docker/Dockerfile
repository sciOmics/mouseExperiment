# syntax=docker/dockerfile:1
################################################################################
# docker/Dockerfile
# The dashboard's server image (mouseexperiment.com)
#
# Built from this repository at one commit, with the repository as the build
# context. Until 0.49.2 this file lived in sciOmicsVPS (docker/Dockerfile.shiny)
# and installed the package from GitHub; it moved here so that the image is
# versioned, and tested, with the code it runs (CODE_REVIEW.md R23-A).
#
# Two builders use it:
#   * CI (.github/workflows/image.yaml) builds it for every pull request and
#     every push to development, and gates the image: docker/shiny-gate.R, a
#     Bayesian fit inside production's memory budget (docker/shiny-gate-fit.R),
#     the test suite (docker/shiny-gate-tests.R) and an HTTP check.
#   * The server (sciOmicsVPS scripts/webhook-handlers/rebuild-shiny.sh) builds
#     it from a source tarball of the production commit, runs shiny-gate.R and
#     the HTTP check, and only then swaps it in.
#
# Manual build, from the repository root:
#   DOCKER_BUILDKIT=1 docker build --build-arg PKG_REF=$(git rev-parse HEAD) \
#     -f docker/Dockerfile -t mouseexperiment:test .
# PKG_REF must name the commit of the working tree: the image's R scripts tell
# users to install exactly that commit.
################################################################################

# sciOmicsVPS CODE_REVIEW.md V3.3/V3.11 -- R 4.5.2 matches development, and
# rocker pins its CRAN to a frozen Posit Package Manager binary snapshot
# (noble/2026-03-10). Do NOT override `repos`: on R 4.4.2 with the live CRAN,
# Deriv 4.3.5 (R >= 4.5) was skipped with only a warning, a dozen dependants
# failed, and the build still exited 0. The snapshot also makes every build
# resolve the same versions and installs binaries (minutes, not the 20-40 min
# source build).
FROM rocker/r-ver:4.5.2

# remotes otherwise turns dependency install errors into warnings and exits 0.
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=false

# System libraries for the binary packages (from Posit's sysreqs for the
# closure below, for Ubuntu 24.04) plus curl for the healthcheck.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libicu-dev \
    libuv1-dev \
    zlib1g-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    pandoc \
    cmake \
    make \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# The non-root runtime user exists before CmdStan is installed, so the
# toolchain can be handed to it (V3.10).
RUN useradd -m -u 1000 shiny

# Every CRAN dependency, in one cached layer. The list is the package's
# Imports (minus base packages) plus the dashboard's packages, which are
# Suggests (dashboard_packages()). An Import missed here is still installed
# below, from the same frozen snapshot.
# warn = 2: a package that fails to install fails the build.
RUN R -q -e " \
  options(warn = 2); \
  pkgs <- c( \
    'remotes', \
    'anytime', 'bayesplot', 'brms', 'car', 'clinfun', 'coin', 'coxphf', \
    'dplyr', 'drc', 'emmeans', 'ggplot2', 'ggpubr', 'lme4', 'lmerTest', \
    'posterior', 'rlang', 'rstan', 'survival', \
    'shiny', 'shinycssloaders', 'shinyFeedback', 'shinyjs', 'DT', 'tidyr', \
    'scales', 'survminer', 'colourpicker', 'sortable', 'htmlwidgets', \
    'logger', 'jsonlite', 'wesanderson', 'future', 'promises' \
  ); \
  install.packages(pkgs); \
  missing <- setdiff(pkgs, rownames(installed.packages())); \
  if (length(missing)) stop('not installed: ', paste(missing, collapse = ', ')) \
"

# CmdStan for the dashboard's 'cmdstanr' backend (V3.10). It used to be
# installed as root into /root/.cmdstan, which the shiny user can neither find
# nor read, so choosing cmdstanr failed with 'no CmdStan toolchain was found'.
# Now: a pinned cmdstanr release and CmdStan version, installed in /opt/cmdstan,
# found through $CMDSTAN, and owned by shiny (CmdStan writes its precompiled
# header and model binaries next to itself). The smoke compile below runs as
# shiny, so the image cannot build unless the shiny user can compile a model.
# CMDSTANR_NO_VER_CHECK: cmdstanr otherwise asks GitHub for the latest CmdStan
# whenever it is attached. CMDSTAN is unset while installing, because cmdstanr
# warns at load (an error under warn = 2) when it names a missing directory.
ARG CMDSTANR_REF=v0.9.0
ARG CMDSTAN_VERSION=2.36.0
ENV CMDSTAN=/opt/cmdstan/cmdstan-${CMDSTAN_VERSION} \
    CMDSTANR_NO_VER_CHECK=true
# install_cmdstan(dir =) requires the directory to exist already.
RUN mkdir -p /opt/cmdstan && CMDSTAN= R -q -e " \
  options(warn = 2); \
  remotes::install_github('stan-dev/cmdstanr@${CMDSTANR_REF}', upgrade = 'never'); \
  options(warn = 1); \
  cmdstanr::check_cmdstan_toolchain(fix = TRUE); \
  cmdstanr::install_cmdstan(dir = '/opt/cmdstan', version = '${CMDSTAN_VERSION}', \
                            cores = 2, quiet = TRUE) \
" && test -x "${CMDSTAN}/bin/stanc" && chown -R shiny:shiny /opt/cmdstan
USER shiny
RUN R -q -e " \
  options(warn = 1); \
  stopifnot(normalizePath(cmdstanr::cmdstan_path()) == normalizePath(Sys.getenv('CMDSTAN'))); \
  f <- cmdstanr::write_stan_file('parameters { real y; } model { y ~ normal(0, 1); }'); \
  invisible(cmdstanr::cmdstan_model(f)) \
"
USER root

# testthat for the gate's test suites (docker/shiny-gate-tests.R). Its own
# layer after CmdStan, so adding it left the slow layers above cached.
# Production images carry it too (a few MB) but never run the tests.
RUN R -q -e " \
  options(warn = 2); \
  pkgs <- c('testthat', 'withr'); \
  install.packages(pkgs); \
  missing <- setdiff(pkgs, rownames(installed.packages())); \
  if (length(missing)) stop('not installed: ', paste(missing, collapse = ', ')) \
"

# Stan models are compiled when a Bayesian fit runs, inside the container's
# memory limit (V3-T). With R's flags (-g ... -O2 ... -g), compiling the
# dashboard's tumour-growth model grew past 3.8 GB, and production's
# compiler was OOM-killed on 2026-10-25. Without debug symbols and at -O1 it
# peaks at 3.0 GB and compiles about 25 s faster; sampling is about 15 %
# slower. R_MAKEVARS_USER applies to compiles at run time only: every package
# is built above, with R's flags. Its own layer, so the layers above stay
# cached.
RUN printf 'CXX17FLAGS = %s\n' \
      "$(R CMD config CXX17FLAGS | sed -e 's/ -g / /g; s/^-g //; s/ -g$//; s/-O2/-O1/')" \
      > /usr/local/lib/R/etc/Makevars.runtime && \
    cat /usr/local/lib/R/etc/Makevars.runtime && \
    ! grep -qE '(^| )-g( |$)' /usr/local/lib/R/etc/Makevars.runtime
ENV R_MAKEVARS_USER=/usr/local/lib/R/etc/Makevars.runtime

# Everything above is shared and cached; the package starts here, so a new
# commit rebuilds only the layers below. .dockerignore keeps the repository's
# docs, review logs and these docker files out of the context, so editing
# them does not rebuild the package either.
# PKG_REF is required: the image's R scripts name the commit to install.
ARG PKG_REF

LABEL org.opencontainers.image.source="https://github.com/sciOmics/mouseExperiment"
LABEL org.opencontainers.image.revision="${PKG_REF}"
LABEL app="mouseExperiment"

# The package, from the build context. install_github() used to add the
# Remote* fields to the installed DESCRIPTION; they are written here instead,
# because the dashboard's R scripts read RemoteSha for their install line and
# package_lock() records the package as a GitHub install. The DESCRIPTION is
# given a final newline first, or the first field appended would join its
# last line. install_deps() installs any Import missing above, from the same
# snapshot. --install-tests keeps tests/ in the library, so the gate runs the
# suite of exactly this commit.
COPY . /tmp/mouseExperiment
RUN : "${PKG_REF:?pass --build-arg PKG_REF=<commit sha>}" && \
    cd /tmp/mouseExperiment && \
    sed -i -e '$a\' DESCRIPTION && \
    printf 'RemoteType: github\nRemoteHost: api.github.com\nRemoteUsername: sciOmics\nRemoteRepo: mouseExperiment\nRemoteRef: %s\nRemoteSha: %s\n' \
      "${PKG_REF}" "${PKG_REF}" >> DESCRIPTION && \
    R -q -e "options(warn = 2); remotes::install_deps('.', upgrade = 'never')" && \
    R CMD INSTALL --install-tests . && \
    R -q -e "stopifnot(identical(utils::packageDescription('mouseExperiment')[['RemoteSha']], '${PKG_REF}'))" && \
    cd / && rm -rf /tmp/mouseExperiment

# What this image actually contains, for the rebuild log and later debugging
# (an image label cannot hold computed output).
RUN mkdir -p /usr/local/share/mouseexperiment && \
    R -q -e " \
      sink('/usr/local/share/mouseexperiment/build-info.txt'); \
      cat('mouseExperiment', as.character(packageVersion('mouseExperiment')), \
          packageDescription('mouseExperiment')[['RemoteSha']], '\n'); \
      cat('CmdStan', as.character(cmdstanr::cmdstan_version()), '\n'); \
      print(sessionInfo()); \
      sink() \
    "

EXPOSE 3838

WORKDIR /srv
RUN chown -R shiny:shiny /srv
USER shiny

CMD ["R", "-e", "mouseExperiment::run_dashboard(host='0.0.0.0', port=3838)"]
