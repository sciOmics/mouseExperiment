# Deploying mouseExperimentDashboard to Shiny Server (VPS)

This guide covers deploying mouseExperimentDashboard to your VPS using Shiny Server.

## Prerequisites

- Ubuntu 24.04 (or compatible Linux)
- Shiny Server installed
- Nginx Proxy Manager (or nginx) for reverse proxy
- mouseExperiment package installed
- Domain configured (e.g., `mouseexperiment.sciomics.bio`)

## Step-by-Step Deployment

### 1. Install Shiny Server

```bash
# Update package list
sudo apt update

# Install R and dependencies
sudo apt install -y r-base r-base-dev gdebi-core

# Download and install Shiny Server
cd /tmp
wget https://www.rstudio.org/download/latest/shiny-server-latest-amd64.deb
sudo gdebi shiny-server-latest-amd64.deb

# Verify installation
sudo systemctl status shiny-server
```

### 2. Install Package Dependencies

```bash
# Install required R packages
sudo R <<'EOF'
packages <- c(
  "shiny", "shinycssloaders", "shinyFeedback", "shinyjs", "DT",
  "dplyr", "ggplot2", "tidyr", "scales", "emmeans", "survminer",
  "colourpicker", "sortable", "htmlwidgets", "ggbeeswarm", "logger"
)

install.packages(packages, repos = "https://cloud.r-project.org/")

# Install mouseExperiment from GitHub
remotes::install_github("sciOmics/mouseExperiment@production")

# Install mouseExperimentDashboard from GitHub
remotes::install_github("sciOmics/mouseExperimentDashboard@production")

# Verify
library(mouseExperimentDashboard)
cat("✓ mouseExperimentDashboard loaded successfully\n")
EOF
```

### 3. Create Application Directory

```bash
# Create app directory
sudo mkdir -p /srv/shiny-server/mouseexperiment

# Set permissions
sudo chown -R shiny:shiny /srv/shiny-server/mouseexperiment
sudo chmod -R 755 /srv/shiny-server/mouseexperiment

# Create subdirectories
sudo mkdir -p /srv/shiny-server/mouseexperiment/{logs,cache}
sudo chown -R shiny:shiny /srv/shiny-server/mouseexperiment/{logs,cache}
```

### 4. Configure Shiny App

Create the app entry point at `/srv/shiny-server/mouseexperiment/app.R`:

```r
# app.R - Shiny Server entry point for mouseExperimentDashboard

library(shiny)
library(mouseExperimentDashboard)

# Ensure www/ folder is accessible
run_shiny_server_app <- function() {
  # Try multiple locations for www folder
  www_paths_to_try <- c(
    # From installed package
    file.path(system.file(package = "mouseExperimentDashboard"), "www"),
    # Current working directory
    file.path(getwd(), "www"),
    # Parent directory (in case of substructure)
    file.path(dirname(getwd()), "www"),
    # Absolute path
    "/srv/shiny-server/mouseexperiment/www"
  )
  
  www_dir <- NULL
  for (path in www_paths_to_try) {
    if (dir.exists(path)) {
      www_dir <- path
      cat(sprintf("[INFO] Using www/ from: %s\n", www_dir))
      break
    }
  }
  
  if (is.null(www_dir)) {
    stop(
      "Could not find www/ folder in any of these locations:\n",
      paste(sprintf("  - %s", www_paths_to_try), collapse = "\n"),
      "\n\nTo fix:\n",
      "1. Check package installation: system.file(package='mouseExperimentDashboard')\n",
      "2. Or create symlink in /srv/shiny-server/mouseexperiment/www"
    )
  }
  
  # Register assets
  shiny::addResourcePath("css",    file.path(www_dir, "css"))
  shiny::addResourcePath("images", file.path(www_dir, "images"))
  shiny::addResourcePath("js",     file.path(www_dir, "js"))
  
  # Launch dashboard
  shiny::runApp(
    mouseExperimentDashboard::runDashboard(host = "0.0.0.0", port = 3838),
    launch.browser = FALSE
  )
}

# Run
run_shiny_server_app()
```

Set correct permissions:
```bash
sudo chown shiny:shiny /srv/shiny-server/mouseexperiment/app.R
sudo chmod 644 /srv/shiny-server/mouseexperiment/app.R
```

### 5. Verify www/ Folder Accessibility

```bash
# Option A: Check if package installation includes www/
pkg_dir=$(R --quiet --vanilla -e "cat(system.file(package='mouseExperimentDashboard'))")
echo "Package location: $pkg_dir"
ls -la "$pkg_dir/www/" 2>/dev/null || echo "www/ not found in package"

# Option B: Create symlink if needed
if [[ ! -d "/srv/shiny-server/mouseexperiment/www" ]]; then
  if [[ -d "$pkg_dir/www" ]]; then
    sudo ln -s "$pkg_dir/www" /srv/shiny-server/mouseexperiment/www
    echo "Created symlink for www/ folder"
  fi
fi

# Verify permissions
sudo ls -la /srv/shiny-server/mouseexperiment/www/
```

### 6. Update Shiny Server Configuration

Edit `/etc/shiny-server/shiny-server.conf`:

```
run_as shiny;
listen 3838;
tcp_backlog 50;
max_connections 100;
max_client_connections 100;

location / {
  site_dir /srv/shiny-server;
  logs_dir /var/log/shiny-server;
  directory_index on;
}

location /mouseexperiment {
  app_dir /srv/shiny-server/mouseexperiment;
}
```

Apply changes:
```bash
sudo systemctl restart shiny-server
```

### 7. Configure Nginx Proxy

In Nginx Proxy Manager or nginx config, add:

```nginx
server {
  server_name mouseexperiment.sciomics.bio;
  listen 80;
  listen [::]:80;

  # Redirect to HTTPS (NPM handles this)
  return 301 https://$server_name$request_uri;
}

server {
  server_name mouseexperiment.sciomics.bio;
  listen 443 ssl http2;
  listen [::]:443 ssl http2;

  # SSL handled by NPM or certbot
  ssl_certificate /etc/letsencrypt/live/mouseexperiment.sciomics.bio/fullchain.pem;
  ssl_certificate_key /etc/letsencrypt/live/mouseexperiment.sciomics.bio/privkey.pem;

  location / {
    # Proxy to Shiny Server
    proxy_pass http://127.0.0.1:3838;
    proxy_http_version 1.1;
    
    # WebSocket support (for Shiny)
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection "upgrade";
    
    # Headers
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
    
    # Timeouts
    proxy_connect_timeout 600s;
    proxy_send_timeout 600s;
    proxy_read_timeout 600s;
  }
}
```

Or if using Nginx Proxy Manager (NPM):
1. Go to Proxy Hosts in NPM admin panel
2. Create new proxy:
   - Domain: `mouseexperiment.sciomics.bio`
   - Scheme: `http`
   - Forward Hostname: `127.0.0.1`
   - Forward Port: `3838`
   - Enable WebSocket Support ✓
   - Enable SSL Certificate ✓
   - Use Let's Encrypt ✓

### 8. Verify Deployment

```bash
# Check Shiny Server is running
sudo systemctl status shiny-server

# Monitor logs
sudo tail -f /var/log/shiny-server.log

# Check specific app logs
sudo ls -la /var/log/shiny-server/

# Test locally
curl -I http://127.0.0.1:3838

# Test via proxy
curl -I https://mouseexperiment.sciomics.bio
```

### 9. Monitor and Troubleshoot

**Real-time logs:**
```bash
sudo tail -f /var/log/shiny-server.log
```

**Check if assets load:**
```bash
# In browser, open Developer Tools (F12)
# Network tab → Filter by CSS/JS
# Look for 404 errors
```

**Restart Shiny Server if needed:**
```bash
sudo systemctl restart shiny-server
```

**Check package installation:**
```bash
R <<'EOF'
library(mouseExperimentDashboard)
cat("Package location:", system.file(package="mouseExperimentDashboard"), "\n")
cat("www/ exists:", dir.exists(system.file("www", package="mouseExperimentDashboard")), "\n")
EOF
```

## Troubleshooting

### Issue: "www/ folder not found"
```bash
# Force symlink
pkg_dir=$(R --quiet --vanilla -e "cat(system.file(package='mouseExperimentDashboard'))")
sudo ln -sf "$pkg_dir/www" /srv/shiny-server/mouseexperiment/www

# Restart
sudo systemctl restart shiny-server
```

### Issue: Formatting broken or CSS not loading
1. Check browser console (F12) for 404s
2. Verify www/ folder: `ls -la /srv/shiny-server/mouseexperiment/www/`
3. Check can be read by shiny user: `sudo su - shiny -s /bin/bash`

### Issue: Date input not working or calculations not triggering
1. Check browser console for JavaScript errors
2. Verify Shiny connection in Network tab (F12)
3. Check for timeouts: tail -f /var/log/shiny-server.log | grep -i error

### Issue: WebSocket connection failing
1. Ensure nginx proxy has WebSocket support enabled
2. Check headers include `Upgrade` and `Connection`
3. Verify proxy timeout settings are sufficient (600s recommended)

## Additional Resources

- [Shiny Server Documentation](https://docs.rstudio.com/shiny-server/)
- [mouseExperimentDashboard GitHub](https://github.com/sciOmics/mouseExperimentDashboard)
- [Nginx Proxy Manager](https://nginxproxymanager.com/)
