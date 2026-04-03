# Troubleshooting mouseExperimentDashboard Production Issues

## Quick Diagnosis

Run the diagnostic script on your VPS:
```bash
bash /path/to/VPS/scripts/99-diagnose-shiny-dashboard.sh
```

This will identify:
- Where the package is installed
- Where the www/ folder is located
- Whether assets are loading correctly
- Package dependencies status

---

## Common Issues & Fixes

### Issue 1: CSS/JS Assets Not Loading (Formatting Broken, Icons Missing)

**Symptoms:**
- No mouse icon in navbar
- Tabs not aligned properly
- Dark mode toggle not working
- Colors look wrong

**Root Cause:**
The `runDashboard()` function in `R/app.R` looks for the `www/` folder using:
```r
pkg_root <- system.file(package = "mouseExperimentDashboard")
www_dir  <- file.path(pkg_root, "www")
if (!dir.exists(www_dir)) www_dir <- file.path(getwd(), "www")
```

If neither path has `www/`, assets won't load.

**Quick Fix:**
Verify the www/ folder exists and is in the right place:
```bash
# Check if package installation includes www/
R -e "cat(system.file(package='mouseExperimentDashboard'), '\n')" 
# Then check if www/ exists there:
ls -la <PKG_PATH>/www/

# If not found, create symlink in working directory:
ln -s /absolute/path/to/www ./www
```

**Fix for Shiny Server:**
Modify how the app is launched. Two options:

**Option A: Use runApp() with app directory (recommended)**
Create `/srv/shiny-server/mouseexperiment/app.R`:
```r
# app.R - Shiny Server entry point
library(shiny)
library(mouseExperimentDashboard)

# Launch the dashboard
mouseExperimentDashboard::runDashboard(host = "0.0.0.0", port = 3838)
```

**Option B: Manually register assets in shiny-server.conf**
```
location /mouseexperiment {
  app_dir /srv/shiny-server/mouseexperiment;
  # If assets are elsewhere, use symbolic link or copy
}
```

---

### Issue 2: Date Input Changes Don't Trigger Calculations

**Symptoms:**
- Changing date picker value doesn't update "Day" calculations
- Reference date update doesn't work
- Date format detection seems inactive

**Root Cause:**
The date input reactivity depends on:
1. JavaScript event handlers being attached properly
2. Shiny session fully initialized
3. The `observe()` blocks in the module servers executing

**Fix Steps:**

1. **Check browser console for JavaScript errors** (F12 → Console):
   - Look for `ReferenceError`, `SyntaxError`, or jQuery errors
   - Check that `shiny` object is defined

2. **Verify Shiny session is initializing correctly:**
   Add this to `R/app.R` temporarily for debugging:
   ```r
   # Add after setup_logging_with_buffer():
   log_info(sprintf("Shiny session initialized: %s", session$ns("")))
   ```

3. **Check if the dateInput is properly connected:**
   In `R/mod_data_upload.R`, the `reference_date` textInput should trigger:
   - If you're using a date picker elsewhere, convert to `dateInput()` 
   - Ensure the input ID matches between UI and server

4. **Force reactivity reset** (add to module server):
   ```r
   observe({
     # This creates a reactive dependency on reference_date
     log_debug(sprintf("Reference date: %s", input$reference_date))
   })
   ```

---

### Issue 3: App Works Locally but Not in Production

**Typical Causes:**
- Package not properly installed (using `devtools::load_all()` vs. production install)
- Working directory different in Shiny Server
- Asset paths hard-coded for development
- Permission issues on files/directories

**Solution:**
1. **Install package properly on VPS:**
   ```bash
   # Option 1: From GitHub
   R -e "remotes::install_github('sciOmics/mouseExperimentDashboard@your-branch')"
   
   # Option 2: From local source
   cd /path/to/mouseExperimentDashboard
   R CMD INSTALL .
   ```

2. **Verify installation:**
   ```r
   library(mouseExperimentDashboard)
   system.file(package = "mouseExperimentDashboard")
   # Should return a path like /usr/local/lib/R/site-library/mouseExperimentDashboard
   
   # Check www/ folder exists:
   list.files(system.file("www", package = "mouseExperimentDashboard"))
   ```

3. **Check permissions:**
   ```bash
   # Shiny Server usually runs as 'shiny' user
   sudo ls -la /usr/local/lib/R/site-library/mouseExperimentDashboard/www/
   ```

---

## Advanced: Shiny Server Configuration

### `/etc/shiny-server/shiny-server.conf` example:

```
# Run as shiny user
run_as shiny;

# Listening port
listen 3838;

# Allow up to 100 concurrent connections
max_connections 100;

# Application directory
location / {
  site_dir /srv/shiny-server;
  logs_dir /var/log/shiny-server;
  directory_index on;
}

# Specific app with higher resource limit
location /mouseexperiment {
  app_dir /srv/shiny-server/mouseexperiment;
  
  # Increase memory limit for R processes
  # (create /srv/shiny-server/mouseexperiment/.shiny_app config file)
}
```

**/srv/shiny-server/mouseexperiment/.shiny_app** (if using environment vars):
```
SHINY_TRACE_CONSOLE=TRUE
R_MAX_MEM_SIZE=2048M
```

---

## Debugging Steps

### Step 1: Check if package loads correctly
```bash
R --quiet <<'EOF'
library(mouseExperimentDashboard)
cat("✓ Package loaded\n")
cat("Package location:", system.file(package="mouseExperimentDashboard"), "\n")
cat("www/ exists:", dir.exists(system.file("www", package="mouseExperimentDashboard")), "\n")
EOF
```

### Step 2: Test asset resolution
```r
# In R console:
pkg_root <- system.file(package = "mouseExperimentDashboard")
www_dir <- file.path(pkg_root, "www")
print(paste0("www/ path: ", www_dir))
print(paste0("exists: ", dir.exists(www_dir)))
if (dir.exists(www_dir)) print(list.files(www_dir))
```

### Step 3: Check Shiny Server logs
```bash
# Real-time logs
sudo tail -f /var/log/shiny-server.log

# Application-specific errors
sudo journalctl -u shiny-server -f

# Shiny app's console output
# Check /var/log/shiny-server/<app-name>-<session>/ directories
```

### Step 4: Browser Developer Tools (F12)
1. **Console tab:** Look for errors (red X)
2. **Network tab:** Check for 404s on CSS/JS files
3. **Application tab → Cookies:** Verify session state persists
4. **Network → XHR/Fetch:** Watch for failed requests to R backend

### Step 5: Start app with verbose logging
In `/srv/shiny-server/mouseexperiment/app.R`:
```r
options(shiny.trace = TRUE)  # Log all shiny operations
options(shiny.sandbox.mode = FALSE)  # Better for debugging

mouseExperimentDashboard::runDashboard(host = "0.0.0.0", port = 3838)
```

---

## Testing from Command Line

Test if the app can start:
```bash
cd /srv/shiny-server/mouseexperiment
R <<'EOF'
library(shiny)
library(mouseExperimentDashboard)

# Simulate Shiny Server environment
Sys.setenv(SHINY_PORT=3838, SHINY_HOST="127.0.0.1")

# Try to run
tryCatch({
  mouseExperimentDashboard::runDashboard(host="127.0.0.1", port=3838, launch.browser=FALSE)
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
})
EOF
```

---

## Next Steps

1. **Run the diagnostic script** to gather baseline info
2. **Share the diagnostic output** for specific guidance
3. **Check browser console** (F12) for JavaScript errors  
4. **Verify installation** using the R commands above
5. **Check Shiny Server logs** for detailed error messages

---

## Contact / More Help

If issues persist after following these steps, gather:
- Output from `99-diagnose-shiny-dashboard.sh`
- Browser console errors (F12 → Console)
- Shiny Server logs (tail -f /var/log/shiny-server.log)
- Package installation path and www/ folder location
