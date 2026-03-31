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

  # JS files are loaded from www/js/ via tags$script(src = ...) in the UI header

  ui <- navbarPage(
    title = tags$span(id = "nav-brand", "mouseExperiment",
               tags$img(src   = "images/mouse_icon.png",
                        id    = "nav-brand-icon",
                        title = "Toggle dark mode",
                        alt   = "Toggle dark mode")),
    id       = "main_nav",
    position = "fixed-top",
    header = tagList(
      shinyjs::useShinyjs(),
      tags$head(
        tags$link(rel = "stylesheet", href = "css/custom.css"),
        tags$link(rel = "stylesheet", href = "css/dark.css"),
        tags$style(HTML("
          body { padding-top: 70px; }
          .navbar-nav li a[data-value='__toxicity_disabled__'],
          .navbar-nav li a[data-value='__dose_response_disabled__'],
          .navbar-nav li a[data-value='__drug_synergy_disabled__'],
          .navbar-nav li a[data-value='__power_analysis_disabled__'] {
            pointer-events: none;
            opacity: 0.45;
          }
        ")),
        tags$script(src = "js/navbar.js"),
        tags$script(src = "js/sidebar.js"),
        tags$script(src = "js/darkmode.js"),
        tags$script(src = "js/overlay.js")
      )
    ),

    tabPanel("Data",          dataUploadUI("data")),
    tabPanel("Tumor Growth",  tumorGrowthUI("tg")),
    tabPanel("Survival",      survivalUI("surv")),
    tabPanel(
      title = tags$span("Toxicity",
                        style = "color: #999; cursor: not-allowed;",
                        title = "Coming soon"),
      value = "__toxicity_disabled__",
      ""
    ),
    tabPanel(
      title = tags$span("Dose Response",
                        style = "color: #999; cursor: not-allowed;",
                        title = "Coming soon"),
      value = "__dose_response_disabled__",
      ""
    ),
    tabPanel(
      title = tags$span("Drug Synergy",
                        style = "color: #999; cursor: not-allowed;",
                        title = "Coming soon"),
      value = "__drug_synergy_disabled__",
      ""
    ),
    tabPanel(
      title = tags$span("Power Analysis",
                        style = "color: #999; cursor: not-allowed;",
                        title = "Coming soon"),
      value = "__power_analysis_disabled__",
      ""
    ),
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

    data_upload_server("data", rv)
    tumor_growth_server("tg", rv)
    survival_server("surv", rv)
    synergy_server("syn", rv)
    power_server("pw", rv)
    dose_response_server("dr", rv)

    # Guard: redirect away from placeholder tabs if reached via keyboard
    observeEvent(input$main_nav, {
      disabled_vals <- c("__toxicity_disabled__", "__dose_response_disabled__",
                         "__drug_synergy_disabled__", "__power_analysis_disabled__")
      if (input$main_nav %in% disabled_vals)
        updateNavbarPage(session, "main_nav", selected = "Survival")
    }, ignoreInit = TRUE)

    # Analysis tabs that require data to be uploaded first
    analysis_tabs <- c("Tumor Growth", "Survival")

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

        # Survival: only when a censor/event column is mapped
        surv_mapped <- !is.null(rv$settings$col_survival_event) &&
                       nzchar(rv$settings$col_survival_event) &&
                       rv$settings$col_survival_event != "None"
        if (surv_mapped) {
          shinyjs::runjs("window.enableNavTab('Survival');")
        } else {
          shinyjs::runjs("window.disableNavTab('Survival');")
        }

      } else {
        shinyjs::runjs(paste(
          sprintf("window.disableNavTab('%s');", analysis_tabs),
          collapse = "\n"
        ))
      }
    })

    # Navigate to Tumor Growth after data configuration
    observeEvent(rv$proceed_to_analysis, {
      if (isTRUE(rv$proceed_to_analysis))
        updateNavbarPage(session, "main_nav", selected = "Tumor Growth")
      rv$proceed_to_analysis <- FALSE
    })
  }

  app <- shinyApp(ui = ui, server = server)
  shiny::runApp(app, host = host, port = port, launch.browser = TRUE)
}
