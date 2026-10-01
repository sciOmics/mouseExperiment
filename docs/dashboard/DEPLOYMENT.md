# Deployment

How the dashboard is deployed to production. This document covers the dashboard-side concerns; the VPS infrastructure (NPM, Cloudflare, Docker base image, webhooks, monitoring) is detailed in the companion `mouseExperimentVPS` repository.

The production URL is `https://mouseexperiment.com`; the staging URL is `https://staging.mouseexperiment.com` (Basic Auth required). Both are R Shiny applications served from long-lived Docker containers behind Nginx Proxy Manager (NPM) on a self-hosted VPS.

---

## Architecture at a glance

```
                           Cloudflare (DNS + caching for static assets)
                                          │
                                          ▼
                                  Public IP : 443
                                          │
                                          ▼
                          ┌─────────────────────────────┐
                          │  Nginx Proxy Manager (NPM)  │
                          │  HTTPS termination          │
                          │  WebSocket upgrade          │
                          │  Long timeouts              │
                          │  Basic Auth (staging only)  │
                          └────────────┬────────────────┘
                                       │
              ┌────────────────────────┼────────────────────────┐
              ▼                                                 ▼
   ┌─────────────────────┐                          ┌─────────────────────┐
   │ shiny-prod:3838     │                          │ shiny-staging:3838  │
   │ Docker container    │                          │ Docker container    │
   │ rocker/r-ver:4.5.2  │                          │ rocker/r-ver:4.5.2  │
   │ runDashboard()      │                          │ runDashboard()      │
   │ branch: production  │                          │ branch: staging     │
   └─────────────────────┘                          └─────────────────────┘
              ▲                                                 ▲
              │                                                 │
              │  GitHub webhook                                 │  GitHub webhook
              │  on push to production                          │  on push to staging
              │                                                 │
   ┌──────────┴─────────────────────────────────────────────────┴──────────┐
   │  webhook receiver (webhook container on deploy.sciomics.bio:9000)     │
   │  runs rebuild-shiny.sh, which:                                        │
   │    1. pins the pushed dashboard SHA and the package's branch head     │
   │    2. builds a candidate image at those two commits                   │
   │    3. gates it (deps load, Stan compiles, UI builds, HTTP 200)        │
   │    4. swaps it in, waits for health, rolls back if unhealthy          │
   └───────────────────────────────────────────────────────────────────────┘
```

Key design decisions inherited from the VPS setup:

- **No ShinyProxy** — each environment runs one shared Shiny process, not per-user containers. Lower overhead; the F.7 session-identity cookie still keeps users isolated.
- **Direct WebSocket proxy** — NPM routes WebSocket traffic straight to the Shiny port (3838).
- **Gated automated rebuilds** — a push to `production` or `staging` of this repository triggers a rebuild that must pass a gate before it replaces the running container. A push to the backend alone triggers nothing; it ships with the next dashboard push or a manual rebuild.

---

## NPM (Nginx Proxy Manager) settings

The dashboard depends on three NPM behaviors. Misconfiguring any of them produces the bugs we saw in CODE_REVIEW.md Round 2 A.3 / A.4.

### 1. WebSocket support

Shiny's reactive layer is a WebSocket. Without `Upgrade` / `Connection: upgrade` headers, the dashboard loads but every reactive input fails silently.

```nginx
proxy_set_header Upgrade $http_upgrade;
proxy_set_header Connection "upgrade";
proxy_cache_bypass $http_upgrade;
```

In the NPM UI: **Edit Proxy Host → Details → Websockets Support: ON**.

### 2. Long proxy timeouts

Bayesian fits take 3-12 minutes. NPM's default `proxy_read_timeout` is 60 seconds; anything over that triggers the Shiny disconnect overlay and the F.7 banner becomes the only way to recover the result.

```nginx
proxy_read_timeout 3600s;
proxy_send_timeout 3600s;
```

In the NPM UI: **Edit Proxy Host → Advanced → Custom Nginx Configuration**.

If you must change this, **don't drop below 900 seconds (15 minutes)** or the F.7 banner can't catch all disconnect scenarios.

### 3. Basic Auth for staging

Staging is behind an Access Control List (ACL) in NPM. Users with the staging URL must enter HTTP Basic Auth credentials before the dashboard loads. Configured via NPM's Access Lists.

Production is public.

See `mouseExperimentVPS/docs/npm-proxy-rules.md` for the canonical NPM settings.

---

## Docker image structure

One Dockerfile in the VPS repo, `docker/Dockerfile.shiny`, serves both environments:

```
FROM rocker/r-ver:4.5.2        R 4.5.2, CRAN pinned to rocker's frozen Posit binary snapshot

Shared cached layers           (identical for prod and staging)
  apt-get: system libraries
  install.packages(<every Import of both packages>)   options(warn = 2): a failure fails the build
  CmdStan 2.36.0 in /opt/cmdstan (cmdstanr v0.9.0), owned by the app user

Per-build layers               keyed by commit SHA (ARG PKG_REF / ARG DASH_REF)
  remotes::install_github('sciOmics/mouseExperiment@<sha>')
  remotes::install_github('sciOmics/mouseExperimentDashboard@<sha>')

CMD R -e "mouseExperimentDashboard::runDashboard(host='0.0.0.0', port=3838)"
```

The two `install_github` calls need `GITHUB_PAT`: a fine-grained token with *Contents:
Read-only* on both private repos and *Commit statuses: Read and write* on this one. It reaches
the build as a BuildKit secret and is never stored in a layer.

### Base image

`rocker/r-ver:4.5.2`, the R version used in development. Its CRAN is a frozen snapshot, so every
build, cold or cached, prod or staging, resolves the same package versions. **Do not override
`repos`:** against the live CRAN, a dependency that required a newer R was skipped with only a
warning and the build still succeeded (VPS CODE_REVIEW.md V3.3). Moving to newer packages means
moving to a newer `rocker/r-ver` tag, and running both test suites against it first.

### What an image contains

`/usr/local/share/mouseexperiment/build-info.txt` lists both packages' versions and commit SHAs,
the CmdStan version and `sessionInfo()`; the label `org.opencontainers.image.revision` holds
the dashboard commit.

---

## Deployment workflow

```
You push to staging
        │
        ▼
GitHub webhook → POST https://deploy.sciomics.bio/hooks/rebuild-shiny-staging
        │
        ▼
webhook container runs scripts/webhook-handlers/rebuild-shiny.sh staging <sha>
        │
        ├─ one rebuild at a time; a push superseded while waiting is skipped
        ├─ build shiny-mouseexperiment-staging:<dash12>-<pkg12>
        ├─ gate 1: docker/shiny-gate.R (every Import loads, rstan and cmdstanr compile
        │          and sample, the UI builds, both packages from the expected commits)
        ├─ gate 2: the candidate serves the app over HTTP
        ├─ tag :latest, recreate shiny-staging, wait for its healthcheck
        │          (unhealthy → back to the previous image)
        └─ GitHub commit status on this commit: success / failure / error
```

If the gate fails, the running container is not touched. Logs: `/var/log/shiny-rebuild/` and
`docker logs webhook` on the VPS.

### Manual rebuild

```bash
ssh -p 2222 <user>@<vps>          # key-only; root login is disabled
docker exec webhook /hooks/rebuild-shiny.sh staging     # or prod
```

### Rollback

```bash
docker image ls shiny-mouseexperiment-prod      # :previous and the last 3 <dash12>-<pkg12> builds
docker tag shiny-mouseexperiment-prod:<tag> shiny-mouseexperiment-prod:latest
cd /opt/vps && docker compose --profile shiny up -d --no-deps --no-build --force-recreate --wait shiny-prod
```

For more, see `mouseExperimentVPS/docs/shiny-deployment.md`.

---

## Cloudflare

DNS is served by Cloudflare. Cloudflare also caches static assets (CSS, JS, images) at the edge — important for the cache-busting `?v=<version>` query string that the dashboard appends to its CSS/JS includes (CODE_REVIEW.md Round 2 A.4).

**Important caveats:**

- **Don't enable Cloudflare's "Rocket Loader" or aggressive JavaScript transformations** — Shiny's reactive bindings break when JS is rewritten.
- **Don't enable Cloudflare's "Auto Minify" for HTML** — Shiny's UI is server-rendered with whitespace that matters.
- **Cache rules:** static `?v=<n>` URLs cache aggressively (1 year). Bare URLs (`/`) cache for 0 seconds. The dashboard's `app.R` appends the version automatically so static assets are immutable per-deploy.

---

## Environment configuration

The dashboard reads no configuration from the environment beyond two variables that only
decide whether to print the proxy-timeout reminder at startup (`SHINY_PORT`,
`SHINY_SERVER_VERSION`; the reminder also prints whenever the app binds a non-local address,
as it does in the container). Everything else is fixed in code:

| Setting | Value | Where |
|---|---|---|
| Host / port | `0.0.0.0:3838` in the container | the image's `CMD` |
| Upload limit | Shiny's default, 5 MB | not set |
| Bayesian fits | one background R worker, one fit at a time, queued | `use_background_worker()` |
| Log format | JSON lines on stdout (`logger`) | `R/logging.R` |
| Asset cache-busting | `?v=<DESCRIPTION Version>` | `dashboard_app()` |

The container's memory limit (production 3.75 GiB, staging 3.25 GiB, no swap) is set in the VPS
`docker-compose.yml`. The background worker marks itself as the kernel's preferred OOM victim,
so an over-limit model compile fails that fit instead of killing the session-serving process.

---

## What can go wrong

| Symptom | Likely cause | Fix |
|---|---|---|
| Disconnect overlay during a Bayesian fit | NPM timeout < fit duration | Raise `proxy_read_timeout` to ≥ 3600s |
| Sidebar inputs do nothing | WebSocket not upgrading | NPM "Websockets Support: ON" |
| CSS / JS not updating after deploy | Browser cached old asset | Force-refresh; check `?v=<version>` is appended |
| `Failed to read shiny-prod from container` after deploy | Container started but Shiny didn't bind | Check the container's logs; usually a missing R dep |
| 502 Bad Gateway from NPM | Container is down or unreachable | `docker ps`, `docker logs shiny-prod` |
| Slow first request after deploy | Shiny process still warming up | Normal — first request triggers package load + Stan model compile cache |
| F.7 banner doesn't show on reconnect | Cookie not set (no HTTPS / dev mode) | Use HTTPS in production. Without the cookie (blocked, or a non-localhost host over HTTP) the session uses a local registry key: fits run but cannot be resumed |
| Bayesian fit fails immediately | brms / Stan toolchain broken | Cannot ship any more: the gate compiles and samples with both backends. On a running container, check `docker logs`; an OOM-killed compile shows as a failed fit |

For deeper troubleshooting: `mouseExperimentVPS/docs/troubleshooting.md`.

---

## Monitoring

Three sources of operational visibility:

1. **Structured logs (stdout):** the dashboard logs every analysis click, every download, every session start/end as a JSON line. Docker keeps them (`docker logs shiny-prod`); there is no central log aggregator.
2. **Uptime checks:** `https://status.sciomics.bio` (Uptime Kuma) checks production and staging every 5 minutes (monitors are added by hand; see `mouseExperimentVPS/docs/monitoring-setup.md`).
3. **Umami analytics:** `https://analytics.mouseexperiment.com` records page views and the custom events emitted by `inst/www/js/umami-events.js` (analysis-start, analysis-complete, download).

If the dashboard isn't responding, the first place to look is the Docker container's stdout:

```bash
ssh -p 2222 <user>@<vps>
docker logs --tail 200 -f shiny-prod
```

For each session there's a `session_start` and `session_end` JSON event with the session token (the Shiny one, not the F.7 sid) so you can correlate "user reported X" with the corresponding logs.

---

## Backup and disaster recovery

The dashboard itself has **no persistent state** — every artifact (uploaded CSV, analysis result, settings preset JSON) is either ephemeral (in-memory) or lives client-side. Container loss doesn't lose user data; it loses any in-flight F.7 jobs.

What is backed up (to Backblaze B2):
- A weekly snapshot of the services' configuration under `/opt/services` (NPM hosts,
  certificates and access lists, Uptime Kuma, Homepage), with consistent SQLite copies
  (`mouseExperimentVPS/scripts/08-backup-server-image.sh`)
- Nightly SQL dumps of the Umami database

What is **not** backed up:
- The Shiny images: they are rebuilt from GitHub by the gated pipeline
- The F.7 `.job_registry` (intentional; in-memory only)
- User-uploaded data in the dashboard (intentional; ephemeral)

If the VPS dies: deploy a fresh VPS with `scripts/00-deploy-all.sh`, restore the snapshot and
dumps, then `scripts/12-deploy-shiny.sh` rebuilds both images (Scenario C in the runbook).

See `mouseExperimentVPS/docs/backup-restore.md` for the full procedure.

---

## Local development vs production

| Concern | Local | Production |
|---|---|---|
| Host | `127.0.0.1` | `0.0.0.0` (Docker) |
| Port | `3838` | `3838` |
| HTTPS | No (plain HTTP) | Yes (Let's Encrypt via NPM) |
| F.7 cookie `Secure` flag | Can't be set (the cookie won't transmit) | Set; cookie works normally |
| F.7 resume behavior | Usually full (browsers treat localhost as a secure context); without the cookie, a local registry key: fits run, no resume | Full resume support |
| Max upload | 5 MB (Shiny default) | Same |
| Bayesian fits | one background worker | Same |
| Opens a browser | yes, in an interactive session | no |
| Logging | stdout | stdout → `docker logs` |

**The plain-HTTP localhost fallback is intentional.** You can't develop F.7 resume against `http://localhost:3838` because the `Secure` cookie can't be set on plain HTTP. To test F.7 locally, either:

- Use a local HTTPS proxy (e.g., `caddy` or `mkcert` + `nginx`) in front of the dashboard
- Test against the staging environment

---

## Security touchpoints

A short hardening checklist; see [`SECURITY.md`](SECURITY.md) for the threat model.

- **HTTPS only:** the F.7 `mE_session` cookie has the `Secure` flag. Plain HTTP can't read it. NPM redirects 80 → 443.
- **HSTS:** set via NPM advanced config. `Strict-Transport-Security: max-age=63072000; includeSubDomains`.
- **CSP:** not yet set. Recommended: `default-src 'self'; script-src 'self' 'unsafe-inline'; connect-src 'self'; frame-ancestors 'none'`. Add via NPM advanced config to mitigate the residual XSS-to-cookie risk noted in [`SECURITY.md`](SECURITY.md).
- **Cloudflare WAF:** enabled. Blocks common XSS / SQLi patterns at the edge before they reach the VPS.
- **SSH:** key-only; root login disabled; CrowdSec on the VPS bans repeat-offender IPs. See `mouseExperimentVPS/docs/security.md`.

---

## See also

- `mouseExperimentVPS/docs/shiny-deployment.md` — full deployment workflow
- `mouseExperimentVPS/docs/npm-proxy-rules.md` — NPM configuration reference
- `mouseExperimentVPS/docs/architecture.md` — VPS architecture overview
- `mouseExperimentVPS/docs/troubleshooting.md` — troubleshooting playbook
- `mouseExperimentVPS/docs/backup-restore.md` — disaster recovery
- [`ARCHITECTURE.md`](ARCHITECTURE.md) — dashboard architecture
- [`SECURITY.md`](SECURITY.md) — threat model and F.7 isolation
- [`DATA_FORMAT.md`](DATA_FORMAT.md) — what data the dashboard accepts
