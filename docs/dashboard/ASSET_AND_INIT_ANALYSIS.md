# mouseExperimentDashboard: Asset Loading & Initialization Analysis

## Executive Summary

The dashboard uses an explicit asset registration pattern via `shiny::addResourcePath()` in the `runDashboard()` function, combined with `tags$head()` declarations in the UI. This ensures that CSS/JS assets are accessible regardless of the working directory or deployment environment. The key initialization points occur in `runDashboard()`, not in `global.R` or start-up hooks.

---

## 1. Asset Loading & Registration

### 1.1 Asset Registration in `runDashboard()` (app.R:16-36)

**Location:** [R/app.R](R/app.R#L16-L36)

```r
pkg_root <- system.file(package = "mouseExperimentDashboard")
www_dir  <- file.path(pkg_root, "www")
if (!dir.exists(www_dir)) www_dir <- file.path(getwd(), "www")
if (dir.exists(www_dir)) {
  shiny::addResourcePath("css",    file.path(www_dir, "css"))
  shiny::addResourcePath("images", file.path(www_dir, "images"))
  shiny::addResourcePath("js",     file.path(www_dir, "js"))
}
```

**How it works:**
1. **Primary path:** Uses `system.file(package=)` to locate the installed package root
   - Under `devtools::load_all()`: returns project root (e.g., `/Users/david/Projects/mouseExperimentDashboard/`)
   - Under production (installed package): returns installed location (e.g., `/Library/Frameworks/R.framework/...`)
2. **Fallback path:** If www/ is not found at `$pkg_root/www`, falls back to `$cwd/www`
3. **Registration:** Maps logical paths (`css`, `images`, `js`) to physical directories via `addResourcePath()`
   - These mappings persist for the entire app session
   - Assets are accessed as `css/custom.css`, `images/mouse_icon.png`, etc. in HTML

### 1.2 Asset Declaration in UI Header (app.R:46-53)

**Location:** [R/app.R](R/app.R#L46-L53)

```r
header = tagList(
  shinyjs::useShinyjs(),
  shinyFeedback::useShinyFeedback(),
  tags$head(
    tags$link(rel = "stylesheet", href = "css/custom.css"),
    tags$link(rel = "stylesheet", href = "css/dark.css"),
    tags$script(src = "js/navbar.js"),
    tags$script(src = "js/sidebar.js"),
    tags$script(src = "js/darkmode.js"),
    tags$script(src = "js/overlay.js")
  )
),
```

**Key points:**
- CSS files are loaded first (Bootstrap styling + custom overrides)
- JS files are loaded at the end (after DOM is ready for DOM manipulation)
- No inline scripts — all functionality is in external `.js` files
- Asset paths are relative (e.g., `js/navbar.js`), not absolute

### 1.3 Assets Directory Structure

```
www/
├── css/
│   ├── custom.css    (custom styling, CSS variables)
│   └── dark.css      (dark mode overrides)
├── images/
│   ├── mouse_icon.png
│   ├── mouse_outline1.png (background layer)
│   └── mouse_outline2.png (background layer)
└── js/
    ├── navbar.js     (tab reorganization, enable/disable logic)
    ├── sidebar.js    (collapsible sidebar panels)
    ├── darkmode.js   (dark mode toggle)
    └── overlay.js    (background images, diagnostics tab disabling)
```

---

## 2. The runDashboard() Function

### 2.1 Location & Signature

**File:** [R/app.R](R/app.R#L11-L200)
**Exported:** Yes (NAMESPACE line 3: `export(runDashboard)`)

```r
runDashboard <- function(host = "127.0.0.1", port = 3838) {
```

### 2.2 Initialization Steps (in order)

| Step | Code | Purpose |
|------|------|---------|
| 1 | Check `shiny` availability | Fail fast if Shiny not installed |
| 2 | Set `APP_VERSION` if not defined | Fallback version constant |
| 3 | **Asset registration** (lines 26-36) | Map logical paths to physical directories |
| 4 | Bind `ui` (navbarPage + modules) | Build UI structure |
| 5 | Bind `server` (line 96+) | Attach server-side logic |
| 6 | Create `shinyApp` object | Encapsulate UI + server |
| 7 | `shiny::runApp()` with `launch.browser=TRUE` | Start Shiny server |

### 2.3 Server-Side Initialization (lines 96-190)

**Location:** [R/app.R](R/app.R#L96-L190)

Key initialization tasks in `server` function:
```r
server <- function(input, output, session) {
  # 1. Setup logging (JSON to stdout + in-process buffer)
  setup_logging_with_buffer()
  log_session(session, "start")
  session$onSessionEnded(function() { log_session(session, "end") })

  # 2. Create reactive values store
  rv <- reactiveValues()

  # 3. Show session disclaimer modal (blocks interaction until dismissed)
  shiny::showModal(shiny::modalDialog(...))

  # 4. Initialize all module servers
  data_upload_server("data", rv)
  tumor_growth_server("tg", rv)
  survival_server("surv", rv)
  synergy_server("syn", rv)
  power_server("pw", rv)
  dose_response_server("dr", rv)
  analysis_log_server("log", rv)

  # 5. Setup reactive observers for tab state management
  observeEvent(input$main_nav, { ... })  # Guard redirects
  observe({ ... })                        # Enable/disable analysis tabs
  observeEvent(rv$proceed_to_analysis, { ... })  # Auto-navigate
}
```

### 2.4 Special Setup & Differences from Standard Shiny

| Aspect | runDashboard() | Standard shinyApp/runApp |
|--------|---|---|
| **Asset paths** | Explicitly registered via `addResourcePath()` | Auto-discovered from `www/` if in app dir |
| **Working directory** | Paths resolved from package root OR cwd | Must match app working directory |
| **Browser launch** | Always `launch.browser = TRUE` | Can be configured |
| **Host/port** | Parameterized with defaults | Must be specified in `runApp()` |
| **Session disclaimer** | Shown to every user (via `showModal()`) | Application-specific |
| **Logging** | Structured JSON logging to stdout | Built-in Shiny logging only |
| **Package check** | Validates `shiny` and `mouseExperiment` | No dependency checks |

---

## 3. Main App.R Initialization & UI Setup

### 3.1 Entry Points

**File:** [app_standalone.R](app_standalone.R)
```r
library(mouseExperimentDashboard)
mouseExperimentDashboard::runDashboard()
```

This is the simplest entry point — just loads the package and calls `runDashboard()`.

**Development entry point:**
```r
devtools::load_all()  # Loads package code into memory
runDashboard()        # Calls the function
```

### 3.2 UI Structure (navbarPage Layout)

**Location:** [R/app.R](R/app.R#L38-L91)

```r
ui <- navbarPage(
  title = tags$span(id = "nav-brand", "mouseExperiment", 
                    tags$img(src="images/mouse_icon.png", ...)),
  id       = "main_nav",
  position = "fixed-top",
  header = tagList(
    # shinyjs & shinyFeedback initialization
    # CSS + JS includes
  ),
  
  # Main tabs (left/center)
  tabPanel("Data",          dataUploadUI("data")),
  tabPanel("Tumor Growth",  tumorGrowthUI("tg")),
  tabPanel("Survival",      survivalUI("surv")),
  tabPanel("Toxicity",      ""),       # Disabled (placeholder)
  tabPanel("Dose Response", ""),       # Disabled (placeholder)
  tabPanel("Drug Synergy",  ""),       # Disabled (placeholder)
  tabPanel("Power Analysis",""),       # Disabled (placeholder)
  
  # Right-aligned tabs (moved by navbar.js)
  tabPanel("Log",           analysisLogUI("log")),
  tabPanel("Help",          helpUI()),
  tabPanel("About",         citationUI())
)
```

**Layout points:**
- **Fixed navbar at top** (`position = "fixed-top"`)
- **Tabs reorganized by JS** (navbar.js moves Log, Help, About to right side)
- **Active tab tracking** via `input$main_nav`
- **Tab enable/disable** via `window.enableNavTab()` / `window.disableNavTab()` (JavaScript)

### 3.3 Global Setup (global.R)

**File:** [global.R](global.R)

This file is sourced **once** when the app starts (before any sessions connect).

**Responsibilities:**
```r
# 1. Dependency check
if (!requireNamespace("mouseExperiment", quietly = TRUE)) {
  stop("mouseExperiment package is required...")
}

# 2. Load required libraries (because runApp() bypasses NAMESPACE)
library(shiny)
library(shinycssloaders)
library(shinyFeedback)
library(shinyjs)
library(DT)
library(dplyr)
library(ggplot2)
library(scales)
library(mouseExperiment)

# 3. Source helper modules
source("R/helpers.R")
source("R/logging.R")
source("R/log_session_modules.R")

# 4. Define application constants
APP_VERSION <- "0.1.0 (Development)"
MAX_FILE_SIZE <- 50 * 1024^2
SUPPORTED_FORMATS <- c("csv", "rda", "RData")
DASHBOARD_COLORS <- list(...)
SAMPLE_DATASETS <- list(...)
DECIMAL_PLACES <- 3
```

**Key difference from runDashboard():**
- `global.R` uses `library()` instead of `requireNamespace()` because the app runs via `shinyAppDir()` or `runApp()`, which does not activate the NAMESPACE
- `runDashboard()` uses `requireNamespace()` and does not rely on global.R (can run standalone)

---

## 4. JavaScript Files & Initialization

### 4.1 JavaScript Loading Sequence

| File | Load Order | Purpose | Initialization |
|------|---|---|---|
| navbar.js | 1st (script tag) | Tab reorganization | Runs on `$(document).ready()` |
| sidebar.js | 2nd | Collapsible sidebar | Runs on `$(document).ready()` |
| darkmode.js | 3rd | Dark mode toggle | Runs on `$(document).ready()` |
| overlay.js | 4th | Background images | Runs on `$(document).ready()` + `shiny:sessioninitialized` |

### 4.2 navbar.js - Tab Management

**File:** [www/js/navbar.js](www/js/navbar.js)

**Init tasks:**
1. Find Log, Help, About tabs in the left `navbar-nav`
2. Detach them from original list
3. Create new `#nav-right` right-aligned list and append to it
4. Mark main nav for CSS centering
5. **Cross-nav active-state sync:** Click handler clears `active` class from other nav list
   - Fixes double-click issue when navigating between main and right tabs
6. **Disabled tab click prevention:** Stops clicks on `.nav-tab-disabled` elements

**Key event handler:**
```javascript
$(document).on('click', '.navbar-center a[data-value]', function() {
  $('#nav-right li.active').removeClass('active');
});
```

This is crucial for proper tab switching when mixing left and right-aligned tabs.

### 4.3 sidebar.js - Collapsible Panels

**File:** [www/js/sidebar.js](www/js/sidebar.js)

**Init tasks:**
1. Listen for clicks on `.sidebar-toggle-btn`
2. Toggle `.sb-collapsed` class on nearest `.collapsible-sidebar-row`
3. Update button text and title based on state

**No date input handlers** — this file only manages sidebar visibility.

### 4.4 darkmode.js - Dark Mode Toggle

**File:** [www/js/darkmode.js](www/js/darkmode.js)

**Init tasks:**
1. Check `localStorage` for `meDashboard_darkMode` key
2. If set to `'1'`, add `dark-mode` class to `<body>` on page load
3. Click handler on `#nav-brand-icon` toggles `dark-mode` class and updates localStorage

**No date input handlers.**

### 4.5 overlay.js - Background Images & Diagnostics

**File:** [www/js/overlay.js](www/js/overlay.js)

**Init tasks:**
1. Create `#mouse-bg-layer` div with two background images (mouse outlines)
2. Append to `<body>`
3. Wait for `shiny:sessioninitialized` event
4. Find Diagnostics tab in `#tg_results_tabs` and disable it
   - Adds `.disabled` class
   - Sets `pointer-events: none` + `cursor: not-allowed`

**⚠️ Key observation:** This uses `$(document).one('shiny:sessioninitialized', ...)` — waits for Shiny to finish initializing before accessing tab structure. This means it's **not** safe to manipulate Shiny-generated tabs in `$(document).ready()` because Shiny hasn't rendered them yet.

---

## 5. Development vs. Production Differences

### 5.1 Asset Path Resolution

#### Development (devtools::load_all())

```r
# In runDashboard():
pkg_root <- system.file(package = "mouseExperimentDashboard")
# Returns: "/Users/david/Projects/mouseExperimentDashboard"
#
# Why: Under devtools::load_all(), the package is loaded from the source directory,
# so system.file() returns the source root, not an installed package path.

www_dir <- file.path(pkg_root, "www")
# Resolves to: "/Users/david/Projects/mouseExperimentDashboard/www"
```

#### Production (Installed Package)

```r
# system.file() returns the installed package directory
pkg_root <- system.file(package = "mouseExperimentDashboard")
# Returns: "/Library/Frameworks/R.framework/Versions/4.x/Resources/library/mouseExperimentDashboard"
#
# OR on Shiny Server:
# Returns: "/srv/shiny-server/.../mouseExperimentDashboard"

www_dir <- file.path(pkg_root, "www")
# Resolves to: "/Library/Frameworks/.../mouseExperimentDashboard/www"
```

**Fallback logic ensures portability:**
```r
if (!dir.exists(www_dir)) www_dir <- file.path(getwd(), "www")
# If package root www/ not found, fall back to current working directory
# This allows development in project folder without installation
```

### 5.2 Dependency Handling

| Scenario | Handling |
|---|---|
| `devtools::load_all()` | NAMESPACE directives are **inactive**; must use `library()` in global.R |
| `library(mouseExperimentDashboard)` + `runDashboard()` | NAMESPACE directives are **active**; package dependencies auto-loaded |
| Shiny Server app (app.R via shinyAppDir()) | NAMESPACE inactive; must use `library()` in global.R |

**Why?** Shiny app directories (including when served via Shiny Server) do NOT activate the NAMESPACE. They treat code as a script, not a loaded package.

### 5.3 Logging & Session Management

Both development and production use the same logging setup:
```r
# Identical in both environments:
setup_logging_with_buffer()      # Called in server function
log_session(session, "start")    # Log session start
session$onSessionEnded(...)      # Log session end
```

### 5.4 Disclaimer Modal

**Identical in both:**
```r
# Same disclaimer shown to every user in both development and production
shiny::showModal(shiny::modalDialog(...))
```

### 5.5 Tab Enable/Disable Logic

| Environment | Method | Reliability |
|---|---|---|
| Development | `window.enableNavTab()` / `window.disableNavTab()` | ✅ Works (JS loads early) |
| Production | `window.enableNavTab()` / `window.disableNavTab()` | ✅ Works (JS loads early) |

The reactive observer that calls these functions is safe:
```r
observe({
  # This runs on startup (no data loaded yet)
  # Disables analysis tabs initially
  shinyjs::runjs("window.disableNavTab('Tumor Growth');")
  shinyjs::runjs("window.disableNavTab('Survival');")
  
  # Enables tabs conditionally as data is loaded
  if (isTRUE(rv$data_valid)) {
    shinyjs::runjs("window.enableNavTab('Tumor Growth');")
  }
})
```

### 5.6 Working Directory Dependency

**Potential issue:** If `runApp()` is called from a different working directory than the app root.

```r
# BEFORE runDashboard() or runApp()
setwd("/Users/david/Projects/mouseExperimentDashboard")
runDashboard()  # ✅ Works — finds www/

# vs.

setwd("/Users/david/Projects")
runDashboard()  # ❌ Fails — cwd/www doesn't exist, but pkg_root/www exists ✅
# Actually works because addResourcePath() uses pkg_root as primary source
```

**Verdict:** `runDashboard()` is **working-directory independent** due to `system.file()` fallback.

---

## 6. Potential Issues & Edge Cases

### 6.1 Date Input Handlers

**Finding:** There are **NO date input event listeners in any JS file.**

This is consistent with the codebase design:
- [helpers.R](R/helpers.R#L38-L70): Contains `detect_date_format()` which handles date parsing **on the R side**
- No custom JS for date inputs means Shiny handles date inputs via Shinyjs reactivity

**If date inputs have issues**, the problem is likely:
1. Date format detection failing (check `detect_date_format()` logic)
2. R-side date parsing issue (not JS-side)
3. Missing validation feedback in UI

### 6.2 Tab Reorganization Timing

**Potential race condition in navbar.js:**
```javascript
// Runs on $(document).ready()
var $mainNav = $('.navbar-nav:first');
var logLi = $mainNav.find('a[data-value="Log"]').closest('li').detach();
```

If the navbar tabs aren't fully rendered by `$(document).ready()`, this could fail.

**Current protection:** Shiny renders the navbar synchronously before JS execution, so this is safe.

### 6.3 Diagnostics Tab Disabling Timing

**overlay.js uses:** `$(document).one('shiny:sessioninitialized', ...)`

This ensures the Shiny session is ready before accessing tabs. ✅ Correct approach.

### 6.4 Asset Loading Under Different Hostnames

If Shiny Server is accessed via:
- `localhost:3838` (development)
- `https://myserver.com/shiny/dashboard` (production with reverse proxy)

**Asset paths remain relative** (`css/custom.css`), so they work correctly in both cases. ✅ No hostname dependency.

### 6.5 CORS Issues (if applicable)

All assets are loaded from the same origin (Shiny server), so CORS is not an issue. ✅

---

## 7. Summary Table: Initialization Flow

| Stage | File | When | What Happens |
|---|---|---|---|
| **1. App Load** | app_standalone.R or global.R | Before any user connects | Libraries loaded, functions sourced |
| **2. runDashboard() called** | R/app.R | Start of session | Assets registered, version set |
| **3. UI Rendered** | R/app.R + R/ui_modules.R | Server listens | navbarPage structure built, CSS/JS tags added to head |
| **4. Browser requests /app** | Shiny Server | HTTP request | HTML sent with `<link>` and `<script>` tags |
| **5. Assets loaded** | www/css + www/js | Browser render phase | CSS applied, JS evaluated |
| **6. JS $(document).ready()** | navbar.js, sidebar.js, darkmode.js, overlay.js | DOM ready | Tab reorganization, click handlers attached |
| **7. Shiny JS loaded** | shinyjs + shinyFeedback | After JS evaluation | Establishes WebSocket connection |
| **8. shiny:sessioninitialized** | (Shiny event) | Session ready | `overlay.js` disables Diagnostics tab, `runjs()` calls safe |
| **9. Server handler runs** | R/app.R server function | Per-session startup | Logging setup, modules initialized, disclaimer shown |
| **10. User interaction** | Various | Ongoing | Reactive observers respond to clicks, form inputs, etc. |

---

## Key Takeaways

1. **Asset registration is explicit** and happens early in `runDashboard()` — not automatic
2. **Working directory independent** — uses `system.file()` with fallback
3. **JS initialization uses proper timing** — `$(document).ready()` for DOM access, `shiny:sessioninitialized` for Shiny access
4. **No date input JS handlers** — date parsing is R-side only
5. **Development and production use identical code paths** — differences are in how packages are loaded, not in app logic
6. **Tab reorganization sync is critical** — navbar.js clears active state between nav lists to prevent double-click issue
7. **Shiny event sequencing is respected** — overlay.js waits for `sessioninitialized` before manipulating Shiny-rendered elements

