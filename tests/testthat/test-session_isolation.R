# Tests for CODE_REVIEW.md F.7 — session-isolation guarantees.
#
# These are the security-critical tests. If any of them ever fail, one
# user can see another's jobs / results. Every layer of the isolation
# story (registry keying, server-side cookie reading, sid unguessability,
# no-leak audit) has at least one assertion here.

skip_if_not_installed("shiny")
skip_if_not_installed("promises")

pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(pkg), "mouseExperiment namespace not loadable")

g <- function(nm) get(nm, envir = pkg)

# Convenience: build a promise that resolves to `value` immediately so
# tests don't need to spin up a future worker.
resolved_promise <- function(value) {
  promises::promise(function(resolve, reject) resolve(value))
}
rejected_promise <- function(err) {
  promises::promise(function(resolve, reject) reject(err))
}

# Force-resolve a promise's callbacks. promises schedules onFulfilled /
# onRejected on the event loop; for sync tests we flush it.
flush_promises <- function() {
  later::run_now(timeoutSecs = 1)
}

# ---------------------------------------------------------------------------
# Layer 1 — registry is keyed by sid; cross-sid lookup returns NULL
# ---------------------------------------------------------------------------

test_that("get_job returns NULL when job_id belongs to a different sid", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)

  sid_a <- "a" %>% strrep(64L)
  sid_b <- "b" %>% strrep(64L)

  job_a <- g("launch_job")(sid_a, "bayes_tg", resolved_promise("A's result"))
  job_b <- g("launch_job")(sid_b, "bayes_tg", resolved_promise("B's result"))

  expect_false(identical(job_a, job_b))
  # A reading A is fine.
  expect_equal(g("get_job")(sid_a, job_a)$kind, "bayes_tg")
  # B reading A returns NULL — this is THE isolation primitive.
  expect_null(g("get_job")(sid_b, job_a))
  # And A reading B returns NULL too.
  expect_null(g("get_job")(sid_a, job_b))
})

test_that("list_jobs returns only this sid's jobs", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)

  sid_a <- strrep("a", 64L)
  sid_b <- strrep("b", 64L)

  g("launch_job")(sid_a, "bayes_tg", resolved_promise(1))
  g("launch_job")(sid_a, "bayes_tg", resolved_promise(2))
  g("launch_job")(sid_b, "bayes_surv", resolved_promise(3))

  jobs_a <- g("list_jobs")(sid_a)
  jobs_b <- g("list_jobs")(sid_b)

  expect_length(jobs_a, 2L)
  expect_length(jobs_b, 1L)
  # No bleed-through: B's only job is the one kind A never launched.
  expect_true(all(vapply(jobs_a, function(j) j$kind, character(1L))
                  == "bayes_tg"))
  expect_equal(jobs_b[[1L]]$kind, "bayes_surv")
})

test_that("list_jobs returns empty list for an unknown sid", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)
  # Empty registry: any sid is unknown.
  expect_equal(g("list_jobs")(strrep("z", 64L)), list())
})

test_that("get_job rejects malformed sid / job_id inputs", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)
  expect_null(g("get_job")(NULL,   "anything"))
  expect_null(g("get_job")("",     "anything"))
  expect_null(g("get_job")(c("a", "b"), "anything"))
  expect_null(g("get_job")(NA_character_, "anything"))
})

# ---------------------------------------------------------------------------
# Layer 2 — server reads cookie from HTTP_COOKIE, not from input$
# ---------------------------------------------------------------------------

test_that("get_sid parses a well-formed cookie header", {
  good_sid <- strrep("a", 64L)
  fake_session <- list(request = list(
    HTTP_COOKIE = paste0("other=1; mE_session=", good_sid, "; theme=dark")))
  expect_equal(g("get_sid")(fake_session), good_sid)
})

test_that("get_sid returns NULL for missing header", {
  fake_session <- list(request = list(HTTP_COOKIE = NULL))
  expect_null(g("get_sid")(fake_session))

  fake_session2 <- list(request = list(HTTP_COOKIE = ""))
  expect_null(g("get_sid")(fake_session2))
})

test_that("get_sid rejects malformed cookie values", {
  # Too short.
  fake_session <- list(request = list(
    HTTP_COOKIE = "mE_session=abc123"))
  expect_null(g("get_sid")(fake_session))

  # Non-hex characters.
  fake_session2 <- list(request = list(
    HTTP_COOKIE = paste0("mE_session=", strrep("x", 64L))))
  expect_null(g("get_sid")(fake_session2))

  # Wrong cookie name.
  fake_session3 <- list(request = list(
    HTTP_COOKIE = paste0("session=", strrep("a", 64L))))
  expect_null(g("get_sid")(fake_session3))
})

test_that("get_sid is the only sid source — input$mE_session_id is never read", {
  # This is enforced by code structure (no exported function reads
  # input$mE_session_id), but assert it: grep for any read of that
  # input across R/ would catch a regression where someone wired the
  # JS path back in. Doc-comments / roxygen lines (lines starting with
  # whitespace + #) are stripped before the scan.
  r_files <- list.files("../../R", pattern = "\\.R$", full.names = TRUE)
  hits <- unlist(lapply(r_files, function(f) {
    src <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    # Drop comment-only lines (the docstring in session_identity.R
    # legitimately MENTIONS input$mE_session_id while explaining we
    # don't read it).
    code <- src[!grepl("^\\s*#", src)]
    grep("input\\$mE_session_id", code, value = TRUE)
  }))
  expect_length(hits, 0L)
})

# ---------------------------------------------------------------------------
# Layer 3 — sid never leaks into logs or banner output
# ---------------------------------------------------------------------------

test_that("no R source references rv$session_id in a log_* / message / stop call", {
  r_files <- list.files("../../R", pattern = "\\.R$", full.names = TRUE)
  bad_patterns <- c(
    "log_(info|warn|error|debug)\\([^)]*session_id",
    "message\\([^)]*session_id",
    "stop\\([^)]*session_id",
    "warning\\([^)]*session_id",
    "showNotification\\([^)]*session_id"
  )
  hits <- unlist(lapply(r_files, function(f) {
    src <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    unlist(lapply(bad_patterns, function(p) grep(p, src, value = TRUE, perl = TRUE)))
  }))
  expect_length(hits, 0L)
})

test_that("make_sid produces 64 hex chars (unguessable on its own)", {
  # 32 fresh sids — assert format and uniqueness.
  sids <- replicate(32L, g("make_sid")())
  expect_true(all(nchar(sids) == 64L))
  expect_true(all(grepl("^[a-f0-9]{64}$", sids)))
  expect_equal(length(unique(sids)), 32L)  # no collisions in 32 draws
})

# ---------------------------------------------------------------------------
# Layer 4 — job lifecycle (registry mechanics, not security per se but
# tested here so the security tests can rely on a working registry).
# ---------------------------------------------------------------------------

test_that("then() callbacks write completion status to the registry", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)
  sid <- strrep("c", 64L)
  job_id <- g("launch_job")(sid, "bayes_tg", resolved_promise(list(value = 42)))
  flush_promises()
  job <- g("get_job")(sid, job_id)
  expect_equal(job$status, "completed")
  expect_equal(job$result$value, 42)
  expect_false(is.na(job$completed_at))
})

test_that("then() callbacks write failure status for rejected promises", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)
  sid <- strrep("d", 64L)
  job_id <- g("launch_job")(sid, "bayes_surv",
    rejected_promise(simpleError("oops")))
  flush_promises()
  job <- g("get_job")(sid, job_id)
  expect_equal(job$status, "failed")
  expect_match(job$error, "oops")
})

test_that("cancel_job before resolution → result is discarded", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)
  sid <- strrep("e", 64L)
  # Build a manual promise we can resolve later.
  resolver <- NULL
  p <- promises::promise(function(resolve, reject) { resolver <<- resolve })
  job_id <- g("launch_job")(sid, "bayes_tg", p)
  expect_true(g("cancel_job")(sid, job_id))
  # Now resolve the underlying future. The onFulfilled callback should
  # detect status == "cancelled" and NOT overwrite the cancellation.
  resolver("late result")
  flush_promises()
  job <- g("get_job")(sid, job_id)
  expect_equal(job$status, "cancelled")
  expect_null(job$result)
})

test_that("cancel_job is idempotent on terminal-state jobs", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)
  sid <- strrep("f", 64L)
  job_id <- g("launch_job")(sid, "bayes_tg", resolved_promise("done"))
  flush_promises()
  expect_false(g("cancel_job")(sid, job_id))  # already completed
})

test_that("evict_stale_jobs drops old done jobs but keeps recent running ones", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)
  sid <- strrep("0", 64L)
  # Job 1: completed long ago.
  jid1 <- g("launch_job")(sid, "bayes_tg", resolved_promise(1))
  flush_promises()
  # Manually backdate via the package-private field setter so the deep
  # assignment goes through the registry's env-aware writer rather than
  # creating a detached copy of the bucket list.
  g(".registry_set")(sid, jid1, "completed_at", Sys.time() - 7200)
  # Job 2: running, started recently.
  resolver <- NULL
  p2 <- promises::promise(function(resolve, reject) { resolver <<- resolve })
  jid2 <- g("launch_job")(sid, "bayes_surv", p2)

  n <- g("evict_stale_jobs")(max_done_age = 1800, max_run_age = 14400)
  expect_equal(n, 1L)
  expect_null(g("get_job")(sid, jid1))
  expect_false(is.null(g("get_job")(sid, jid2)))
})

test_that("retire_job removes a job and cleans up empty sid bucket", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)
  sid <- strrep("1", 64L)
  jid <- g("launch_job")(sid, "bayes_tg", resolved_promise(1))
  expect_true(g("retire_job")(sid, jid))
  expect_null(g("get_job")(sid, jid))
  # Bucket emptied → entirely removed.
  expect_null(g(".job_registry")$sessions[[sid]])
})
