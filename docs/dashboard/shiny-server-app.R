# run.R - Alternative entry point for Shiny Server deployment
# Place in: /srv/shiny-server/mouseexperiment/app.R
#
# This file properly handles asset paths for Shiny Server deployment
# and provides better debugging information.

#' Alternative runDashboard wrapper for Shiny Server
#' 
#' This version:
#' - Handles asset paths robustly
#' - Provides clear error messages if www/ is missing
#' - Works correctly when launched by Shiny Server
#' - Logs asset resolution for debugging

run_shiny_server_app <- function() {
  library(shiny)
  library(mouseExperimentDashboard)
  
  # Get environment info
  is_shiny_server <- nzchar(Sys.getenv("SHINY_PORT"))
  current_wd <- getwd()
  
  # Ensure www/ folder is accessible
  www_paths_to_try <- c(
    # 1. Package installation
    file.path(system.file(package = "mouseExperimentDashboard"), "www"),
    # 2. Current working directory
    file.path(current_wd, "www"),
    # 3. Parent of current working directory (if in subdirectory)
    file.path(dirname(current_wd), "www"),
    # 4. Explicit /srv/shiny-server path
    "/srv/shiny-server/mouseexperiment/www"
  )
  
  # Find first existing www/ path
  www_dir <- NULL
  for (path in www_paths_to_try) {
    if (dir.exists(path)) {
      www_dir <- path
      cat(sprintf("[INFO] Found www/ at: %s\n", www_dir))
      break
    }
  }
  
  # Fail loudly if www/ not found
  if (is.null(www_dir)) {
    stop(
      "mouseExperimentDashboard www/ folder not found.\n",
      "Searched in:\n",
      paste(sprintf("  - %s", www_paths_to_try), collapse = "\n"),
      "\n\nTo fix:\n",
      "1. Reinstall package: devtools::install_local('../../mouseExperimentDashboard')\n",
      "2. OR create symlink: ln -s /path/to/www ./www\n",
      "3. Check package installation location: system.file(package='mouseExperimentDashboard')"
    )
  }
  
  # Verify required asset files exist
  required_files <- c(
    file.path(www_dir, "css", "custom.css"),
    file.path(www_dir, "js", "navbar.js"),
    file.path(www_dir, "images", "mouse_icon.png")
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    warning(sprintf(
      "[WARN] Some required files missing from www/:\n%s",
      paste(sprintf("  - %s", missing_files), collapse = "\n")
    ))
  }
  
  # Register resource paths
  shiny::addResourcePath("css",    file.path(www_dir, "css"))
  shiny::addResourcePath("images", file.path(www_dir, "images"))
  shiny::addResourcePath("js",     file.path(www_dir, "js"))
  
  cat(sprintf(
    "[INFO] Shiny Server Mode: %s | Working Dir: %s\n",
    is_shiny_server, current_wd
  ))
  
  # Call the main dashboard function
  mouseExperimentDashboard::runDashboard(
    host = "0.0.0.0",
    port = 3838
  )
}

# Launch
tryCatch({
  run_shiny_server_app()
}, error = function(e) {
  cat(sprintf("[ERROR] Failed to start dashboard:\n%s\n", e$message))
  stop(e)
})
