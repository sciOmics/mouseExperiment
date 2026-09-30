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
   │ rocker/r-ver:4.4    │                          │ rocker/r-ver:4.4    │
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
   │  triggers rebuild-shiny.sh which:                                     │
   │    1. git pull both repos                                             │
   │    2. docker build (rocker base + R deps + the two packages)          │
   │    3. docker compose up -d the relevant container                     │
   └───────────────────────────────────────────────────────────────────────┘
```

Key design decisions inherited from the VPS setup:

- **No ShinyProxy** — each environment runs one shared Shiny process, not per-user containers. Lower overhead; the F.7 session-identity cookie still keeps users isolated.
- **Direct WebSocket proxy** — NPM routes WebSocket traffic straight to the Shiny port (3838).
- **Automated rebuilds** — GitHub webhooks trigger image rebuilds on every push to `production` or `staging`.

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

The Dockerfile lives in the VPS repo (`docker/Dockerfile.shiny-prod` and `Dockerfile.shiny-staging`). Layered for cache efficiency:

```
Layer 1 (slow, cached):  OS packages
  apt-get install libcurl4-openssl-dev libssl-dev ...

Layer 2 (slow, cached):  CRAN packages
  R -e "install.packages(c('shiny', 'ggplot2', 'DT', 'brms', ...))"
  → ~15-20 minutes first time, cached afterwards

Layer 3 (fast, rebuilt per push):  Both packages from GitHub
  R -e "devtools::install_github('sciOmics/mouseExperiment',          ref='production')"
  R -e "devtools::install_github('sciOmics/mouseExperimentDashboard', ref='production')"

Layer 4 (entrypoint):
  CMD R -e "mouseExperimentDashboard::runDashboard(host='0.0.0.0', port=3838)"
```

The two GitHub `install_github` calls require a `GITHUB_PAT` with read scope on both private repos.

### Base image

`rocker/r-ver:4.4` — a maintained R image with the toolchain pre-installed. Pinned R version so a CRAN-side R version bump doesn't break our build silently.

### Image size

~2.5 GB. Most of the size is R + Stan + CRAN dependencies; the dashboard package itself is < 5 MB. Layer caching makes incremental rebuilds (after a code-only push) take 1-2 minutes.

---

## Deployment workflow

```
You push to staging
        │
        ▼
GitHub webhook → POST https://deploy.sciomics.bio/hooks/rebuild-shiny-staging
        │
        ▼
webhook container runs rebuild-shiny.sh
        │
        ├─ git fetch + checkout staging branch in both repos
        ├─ docker build --target staging
        ├─ docker compose up -d shiny-staging
        ├─ Wait for HTTP 200 from staging.mouseexperiment.com
        └─ Slack notification on success/failure
```

`rebuild-shiny.sh` is part of the VPS repo (`mouseExperimentVPS/scripts/12-deploy-shiny.sh`). Important: the script has a "verify HTTP 200" gate, added after CODE_REVIEW.md Round 2 B.1 (a CRAN miss had left the prod image broken because the build failed silently). If the rebuild's smoke test fails, the script keeps the previous image running and the deployment is reported as failed.

### Manual rebuild

```bash
ssh <user>@<server>
cd /home/admin/mouseExperimentVPS
bash scripts/12-deploy-shiny.sh staging       # or production
```

### Rollback

```bash
# List the last N image tags
docker images shiny-mouseexperiment-prod

# Restart with a specific tag
docker tag shiny-mouseexperiment-prod:<previous-sha> shiny-mouseexperiment-prod:latest
docker compose restart shiny-prod
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

The container reads environment variables for configurable behavior:

| Env var | Purpose | Default |
|---|---|---|
| `SHINY_PORT` | Port to listen on | `3838` |
| `SHINY_HOST` | Host to bind to | `0.0.0.0` (Docker) |
| `SHINY_LOG_LEVEL` | `logger` package log level | `INFO` |
| `SHINY_LOG_BUFFER_LIMIT` | In-memory log ring buffer size | `5000` |
| `SHINY_FUTURE_WORKERS` | Worker count for `future::plan(multisession)` | `availableCores() - 1` |
| `SHINY_MAX_UPLOAD_MB` | Max upload size (MB) | `5` (dev) / `100` (prod) |
| `MOUSE_EXP_DASHBOARD_VERSION` | Override version for cache-busting | DESCRIPTION's Version |

Set via `docker compose.yml` or the compose-override file on the VPS.

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
| F.7 banner doesn't show on reconnect | Cookie not set (no HTTPS / dev mode) | Use HTTPS in production; localhost dev falls back to legacy inline |
| Bayesian fit fails immediately | brms / Stan toolchain broken | Check Dockerfile's Layer 2; rebuild from scratch |

For deeper troubleshooting: `mouseExperimentVPS/docs/troubleshooting.md`.

---

## Monitoring

Three sources of operational visibility:

1. **Structured logs (stdout):** the dashboard logs every analysis click, every download, every session start/end as a JSON line. The VPS captures these via Docker logging driver and ships them to a central log aggregator.
2. **Uptime checks:** `https://status.sciomics.bio` (Uptime Kuma) pings the prod and staging URLs every 60 seconds.
3. **Umami analytics:** `https://analytics.mouseexperiment.com` records page views and the custom events emitted by `inst/www/js/umami-events.js` (analysis-start, analysis-complete, download).

If the dashboard isn't responding, the first place to look is the Docker container's stdout:

```bash
ssh <user>@<server>
docker logs --tail 200 -f shiny-prod
```

For each session there's a `session_start` and `session_end` JSON event with the session token (the Shiny one, not the F.7 sid) so you can correlate "user reported X" with the corresponding logs.

---

## Backup and disaster recovery

The dashboard itself has **no persistent state** — every artifact (uploaded CSV, analysis result, settings preset JSON) is either ephemeral (in-memory) or lives client-side. Container loss doesn't lose user data; it loses any in-flight F.7 jobs.

What is backed up:
- The VPS image itself (via Backblaze; see `mouseExperimentVPS/scripts/08-backup-server-image.sh`)
- The Docker images (rebuildable from GitHub)
- The NPM configuration database

What is **not** backed up:
- The F.7 `.job_registry` (intentional; in-memory only)
- User-uploaded data in the dashboard (intentional; ephemeral)

If the VPS dies:
1. Restore the VPS image from Backblaze (~20 minutes)
2. Re-run `scripts/00-deploy-all.sh` (~30 minutes; rebuilds images from GitHub)
3. Total downtime: ~1 hour

See `mouseExperimentVPS/docs/backup-restore.md` for the full procedure.

---

## Local development vs production

| Concern | Local | Production |
|---|---|---|
| Host | `127.0.0.1` | `0.0.0.0` (Docker) |
| Port | `3838` | `3838` |
| HTTPS | No (plain HTTP) | Yes (Let's Encrypt via NPM) |
| F.7 cookie `Secure` flag | Can't be set (the cookie won't transmit) | Set; cookie works normally |
| F.7 resume behavior | Falls back to legacy inline (no resume) | Full resume support |
| Max upload | 5 MB | 100 MB |
| Workers | `availableCores() - 1` | Same |
| Logging | stdout | stdout → Docker → aggregator |

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
