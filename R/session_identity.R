#' Session-identity cookie helpers (CODE_REVIEW.md F.7, layer 2).
#'
#' The server-authoritative reader of the `mE_session` cookie. Reads
#' the HTTP `Cookie` header directly from `session$request$HTTP_COOKIE`,
#' never trusting any JS-side input value. This is the anti-spoofing
#' boundary that prevents user B from claiming to be user A by sending
#' a forged websocket message setting `input$mE_session_id`.
#'
#' The cookie is set client-side by `inst/www/js/session-id.js` on
#' first connect, using `crypto.getRandomValues` (256 bits of CSPRNG
#' output) and attributes `SameSite=Strict; Secure; Path=/;
#' max-age=86400`. The server never sees the JS code or the client-side
#' generation step — it only reads whatever the browser presents in the
#' Cookie header.
#'
#' Cookie attributes:
#'   * SameSite=Strict — no cross-site request can submit the cookie
#'   * Secure          — HTTPS only; never transmits on plain HTTP
#'   * Path=/          — scoped to this app
#'   * max-age=86400   — 24-hour expiry
#'   * Not HttpOnly    — JS must read it to set it the first time. XSS
#'                       audit covers the residual risk.

# Regex matches "mE_session=<64 hex>" anywhere in a Cookie header. The
# header is a semicolon-separated list of name=value pairs; we don't
# need to split because the perl-regex lookbehind anchors on the name.
.SESSION_COOKIE_RE <- "(?<=mE_session=)[a-f0-9]{64}"

#' Read the session-identity cookie from a Shiny session.
#'
#' Returns the sid as a character scalar, or NULL when the cookie is
#' absent / malformed. NEVER reads from `input$...` — the cookie value
#' comes only from the HTTP Cookie header that the browser sent.
#'
#' @param session A Shiny session object (e.g. from `moduleServer()`).
#' @return Character scalar (64 hex chars) or NULL.
#' @keywords internal
get_sid <- function(session) {
  cookie_header <- tryCatch(session$request$HTTP_COOKIE, error = function(e) NULL)
  if (is.null(cookie_header) || !nzchar(cookie_header)) return(NULL)
  m <- regmatches(cookie_header,
                  regexpr(.SESSION_COOKIE_RE, cookie_header, perl = TRUE))
  if (length(m) == 1L && nchar(m) == 64L) m else NULL
}

#' Generate a session ID server-side.
#'
#' Used as a fallback when the JS path is unavailable (e.g. in unit
#' tests with `shiny::testServer()` where there's no browser to set a
#' cookie). 256 bits of entropy from R's PRNG — for tests this is
#' enough; in production the browser-side `crypto.getRandomValues`
#' path is the primary source.
#'
#' @return Character scalar (64 hex chars).
#' @keywords internal
make_sid <- function() {
  # Sample 32 bytes (= 64 hex chars) from R's PRNG. Not as strong as
  # the browser CSPRNG but adequate for fallback identity that doesn't
  # cross the network.
  bytes <- sample.int(255L, 32L, replace = TRUE)
  paste0(format(as.hexmode(bytes), width = 2L), collapse = "")
}

#' Tag for a missing/unset sid.
#'
#' Returned by get_sid_or_warn() so callers can distinguish "the user
#' hasn't established a sid yet" from "the user is sid X". Useful for
#' the banner observer which should silently no-op when the sid hasn't
#' arrived yet (the JS init runs after Shiny's first server tick).
#' @keywords internal
SID_NONE <- ""

#' Best-effort sid resolution suitable for top-of-server use.
#'
#' Reads the cookie (the only authoritative source). If absent, returns
#' SID_NONE so callers can pattern-match against it without testing for
#' NULL.
#'
#' @keywords internal
get_sid_or_none <- function(session) {
  sid <- get_sid(session)
  if (is.null(sid)) SID_NONE else sid
}

#' Registry key for a session that presented no identity cookie
#'
#' CODE_REVIEW.md R20.D3. Such sessions used to take a "legacy inline" path
#' whose promise operators were never imported, so every Bayesian fit failed
#' before it started and left the "fit running" modal open. They now use the job
#' registry like every other session, under this key.
#'
#' The key starts with "local-". `get_sid()` accepts only 64 hex characters
#' from the Cookie header, so no browser can ever present it: its jobs are
#' reachable only from the session that created them, whatever the quality of
#' the random suffix. Resuming after a disconnect still needs the cookie, so
#' these sessions cannot resume -- the same as before.
#'
#' @return Character scalar: "local-" followed by 64 hex characters.
#' @keywords internal
make_local_sid <- function() paste0("local-", make_sid())
