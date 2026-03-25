#' Run mouseExperimentDashboard
#'
#' Launch the interactive dashboard for analyzing mouse tumor growth experiments.
#'
#' @param host The host IP address to listen on. Default: "127.0.0.1" (localhost)
#' @param port The port number. Default: 3838
#'
#' @importFrom shiny shinyApp navbarPage tabPanel fluidPage tags observeEvent updateNavbarPage reactiveValues observe
#' @export
#'
#' @examples
#' \dontrun{
#' runDashboard()
#' runDashboard(port = 8080)
#' }
runDashboard <- function(host = "127.0.0.1", port = 3838) {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required but not installed")
  }

  if (!exists("APP_VERSION")) APP_VERSION <- "0.1.0 (Development)"

  # Mount static assets explicitly so they resolve regardless of working
  # directory. system.file(package=) returns the project root under
  # devtools::load_all(), so "www" is always one level below it.
  pkg_root <- system.file(package = "mouseExperimentDashboard")
  www_dir  <- file.path(pkg_root, "www")
  if (!dir.exists(www_dir)) www_dir <- file.path(getwd(), "www")
  if (dir.exists(www_dir)) {
    shiny::addResourcePath("css",    file.path(www_dir, "css"))
    shiny::addResourcePath("images", file.path(www_dir, "images"))
    shiny::addResourcePath("js",     file.path(www_dir, "js"))
  }

  # JS injected on page load:
  #  1. Moves Help + Citation to a right-aligned nav
  #  2. Adds "navbar-center" class to the main nav for CSS centering
  #  3. Exposes disableNavTab / enableNavTab helpers for server-side use
  nav_js <- HTML("
$(document).ready(function() {
  // Move Help and Citation to a right-aligned nav list
  var $mainNav = $('.navbar-nav:first');
  var helpLi     = $mainNav.find('a[data-value=\"Help\"]').closest('li').detach();
  var citationLi = $mainNav.find('a[data-value=\"About\"]').closest('li').detach();
  var $rightNav  = $('<ul class=\"nav navbar-nav navbar-right\" id=\"nav-right\"></ul>');
  $rightNav.append(helpLi).append(citationLi);
  // Bootstrap 5: no .navbar-collapse wrapper — append directly to .container-fluid
  var $container = $mainNav.closest('.container-fluid, .navbar-collapse');
  if ($container.length === 0) $container = $mainNav.parent();
  $container.append($rightNav);


  // Mark main nav for centering via CSS
  $mainNav.addClass('navbar-center');

  // Cross-nav active-state sync:
  // Bootstrap only clears 'active' within the same <ul>. Because Help/Citation
  // live in #nav-right (a separate <ul>), clicking a main-nav tab leaves them
  // visually active, causing a double-click requirement.
  // Solution: clear active from the *other* nav whenever a tab is clicked.
  $(document).on('click', '#nav-right a[data-value]', function() {
    $mainNav.find('li.active').removeClass('active');
  });
  $(document).on('click', '.navbar-center a[data-value]', function() {
    $('#nav-right li.active').removeClass('active');
  });
});

// Prevent clicks on disabled tabs
$(document).on('click', '.nav-tab-disabled > a, .nav-tab-disabled > .nav-link', function(e) {
  e.preventDefault(); e.stopPropagation(); return false;
});

function disableNavTab(val) {
  $('.navbar-nav li a[data-value=\"' + val + '\"]').closest('li').addClass('nav-tab-disabled');
  $('#nav-right li a[data-value=\"' + val + '\"]').closest('li').addClass('nav-tab-disabled');
}
function enableNavTab(val) {
  $('.navbar-nav li a[data-value=\"' + val + '\"]').closest('li').removeClass('nav-tab-disabled');
  $('#nav-right li a[data-value=\"' + val + '\"]').closest('li').removeClass('nav-tab-disabled');
}
window.disableNavTab = disableNavTab;
window.enableNavTab  = enableNavTab;
")

  ui <- navbarPage(
    title = tags$span(id = "nav-brand", "mouseExperiment",
               tags$img(src   = "images/mouse_icon.png",
                        id    = "nav-brand-icon",
                        title = "Toggle dark mode",
                        alt   = "Toggle dark mode")),
    id    = "main_nav",
    header = tagList(
      shinyjs::useShinyjs(),
      tags$head(
        tags$link(rel = "stylesheet", href = "css/custom.css"),
        tags$link(rel = "stylesheet", href = "css/dark.css"),
        tags$script(nav_js),
        tags$script(src = "js/sidebar.js"),
        tags$script(src = "js/darkmode.js"),
        # Inject mouse overlay as a direct child of <body> to avoid
        # z-index conflicts with the Bootstrap navbar.
        tags$script(HTML("
$(document).ready(function() {
  var layer = $('<div id=\"mouse-bg-layer\">' +
    '<img id=\"mouse-bg-img1\" class=\"mouse-bg-img\" src=\"images/mouse_outline1.png\" alt=\"\">' +
    '<img id=\"mouse-bg-img2\" class=\"mouse-bg-img\" src=\"images/mouse_outline2.png\" alt=\"\">' +
    '</div>');
  $('body').append(layer);
});
"))
      )
    ),

    tabPanel("Data",          dataUploadUI()),
    tabPanel("Tumor Growth",  tumorGrowthUI()),
    tabPanel("Survival",      survivalUI()),
    tabPanel("Dose Response", doseResponseUI()),
    tabPanel("Drug Synergy",  synergyUI()),
    tabPanel("Power Analysis",powerAnalysisUI()),
    tabPanel("Help",          helpUI()),
    tabPanel("About",          citationUI())
  )

  server <- function(input, output, session) {

    rv <- reactiveValues()

    # Startup disclaimer modal — shown once per session
    shiny::showModal(
      shiny::modalDialog(
        title = shiny::tagList(
          shiny::icon("triangle-exclamation"), " Disclaimer"
        ),
        shiny::p(
          shiny::strong("mouseExperiment"), " is provided \"as is\" for research purposes only,",
          " without warranty of any kind, express or implied. The authors and contributors",
          " accept no liability for errors, inaccuracies, or any consequences — including",
          " but not limited to incorrect results, data loss, or decisions made on the basis",
          " of outputs from this software.",
          " All results must be independently validated before use in publications, reports,",
          " regulatory submissions, clinical decisions, or any other application."
        ),
        shiny::p(
          "By clicking ", shiny::strong("I Understand"),
          ", you acknowledge that you have read this disclaimer and accept full responsibility",
          " for validating any results produced by this software."
        ),
        footer = shiny::modalButton("I Understand"),
        easyClose = FALSE
      )
    )

    data_upload_server(input, output, session, rv)
    tumor_growth_server(input, output, session, rv)
    survival_server(input, output, session, rv)
    synergy_server(input, output, session, rv)
    power_server(input, output, session, rv)
    dose_response_server(input, output, session, rv)

    # Analysis tabs that require data to be uploaded first
    analysis_tabs <- c("Tumor Growth", "Survival", "Dose Response",
                       "Drug Synergy", "Power Analysis")

    # Disable all analysis tabs on startup
    shinyjs::runjs(paste(
      sprintf("window.disableNavTab('%s');", analysis_tabs),
      collapse = "\n"
    ))

    # Enable / disable analysis tabs reactively
    observe({
      if (isTRUE(rv$data_valid)) {
        # Always available once data is loaded
        shinyjs::runjs("window.enableNavTab('Tumor Growth');")
        shinyjs::runjs("window.enableNavTab('Power Analysis');")
        shinyjs::runjs("window.enableNavTab('Drug Synergy');")

        # Survival: only when a censor/event column is mapped
        surv_mapped <- !is.null(input$col_survival_event) &&
                       nzchar(input$col_survival_event) &&
                       input$col_survival_event != "None"
        if (surv_mapped) {
          shinyjs::runjs("window.enableNavTab('Survival');")
        } else {
          shinyjs::runjs("window.disableNavTab('Survival');")
        }

        # Dose Response: only when a dose column is mapped
        dose_mapped <- (!is.null(input$col_dose) &&
                        nzchar(input$col_dose) &&
                        input$col_dose != "None") ||
                       (!is.null(rv$dose_col) && nzchar(rv$dose_col))
        if (dose_mapped) {
          shinyjs::runjs("window.enableNavTab('Dose Response');")
        } else {
          shinyjs::runjs("window.disableNavTab('Dose Response');")
        }

      } else {
        shinyjs::runjs(paste(
          sprintf("window.disableNavTab('%s');", analysis_tabs),
          collapse = "\n"
        ))
      }
    })

    # Navigate to Tumor Growth after data configuration
    observeEvent(input$proceed_analysis, {
      updateNavbarPage(session, "main_nav", selected = "Tumor Growth")
    })
  }

  app <- shinyApp(ui = ui, server = server)
  shiny::runApp(app, host = host, port = port, launch.browser = TRUE)
}
