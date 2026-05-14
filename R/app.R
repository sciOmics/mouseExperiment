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

  # Mount static assets from inst/www/ (bundled with package)
  # system.file() returns the installation path during normal use, and the
  # source directory root during devtools::load_all(), so inst/www resolves correctly in both cases.
  www_dir <- system.file("www", package = "mouseExperimentDashboard")
  
  if (!nzchar(www_dir) || !dir.exists(www_dir)) {
    stop(
      "Static assets (www/) folder not found in package installation.\n",
      "Expected location: ", system.file(package = "mouseExperimentDashboard"), "/www\n",
      "To fix: Reinstall the package to ensure inst/www/ is included:\n",
      "  devtools::install() or R CMD INSTALL ."
    )
  }
  
  shiny::addResourcePath("css",    file.path(www_dir, "css"))
  shiny::addResourcePath("images", file.path(www_dir, "images"))
  shiny::addResourcePath("js",     file.path(www_dir, "js"))

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
      shinyFeedback::useShinyFeedback(),
      tags$head(
        tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
        tags$meta(property = "og:type",        content = "website"),
        tags$meta(property = "og:url",         content = "https://mouseexperiment.com/"),
        tags$meta(property = "og:title",       content = "mouseExperiment Dashboard"),
        tags$meta(property = "og:description", content = "Interactive survival and phenotype analysis for mouse experiment data. Easy data upload, curve fitting, and statistical analysis for preclinical research."),
        tags$meta(property = "og:image",       content = "https://mouseexperiment.com/images/og_preview.png"),
        tags$meta(property = "og:image:width",  content = "1200"),
        tags$meta(property = "og:image:height", content = "630"),
        tags$meta(name = "twitter:card",        content = "summary_large_image"),
        tags$meta(name = "twitter:title",       content = "mouseExperiment Dashboard"),
        tags$meta(name = "twitter:description", content = "Interactive survival and phenotype analysis for mouse experiment data."),
        tags$meta(name = "twitter:image",       content = "https://mouseexperiment.com/images/og_preview.png"),
        tags$link(rel = "stylesheet", href = "css/custom.css"),
        tags$link(rel = "stylesheet", href = "css/dark.css"),
        tags$script(src = "js/navbar.js"),
        tags$script(src = "js/sidebar.js"),
        tags$script(src = "js/darkmode.js"),
        tags$script(src = "js/overlay.js"),
        tags$script(HTML(
          "$(document).ready(function(){",
          "  if(window.innerWidth>768)return;",
          "  var o=$('<div style=\"position:fixed;top:0;left:0;right:0;bottom:0;z-index:99999;",
          "    background:#f5f5f5;color:#333;display:flex;flex-direction:column;",
          "    align-items:center;justify-content:center;text-align:center;",
          "    padding:2rem 1.5rem\">",
          "    <div style=\"font-size:4rem;margin-bottom:1rem\">&#x1F4BB;</div>",
          "    <h2 style=\"font-size:1.4rem;margin-bottom:.75rem\">Desktop Recommended</h2>",
          "    <p style=\"font-size:1rem;line-height:1.5;max-width:360px;",
          "      margin:0 auto 1.5rem;color:#555\">This dashboard requires a larger screen ",
          "      for the best experience. Please visit on a desktop or laptop computer.</p>",
          "    <button id=\"mobile-dismiss\" style=\"padding:.5rem 1.5rem;border:1px solid #999;",
          "      border-radius:4px;background:transparent;color:#555;font-size:.85rem;",
          "      cursor:pointer\">Continue anyway</button></div>');",
          "  $('body').append(o);",
          "  o.on('click','#mobile-dismiss',function(){o.remove();});",
          "});"
        )),
        tags$script(
          async = NA,
          src = "https://analytics.headpudding.com/script.js",
          `data-website-id` = "2bab3ca9-ec31-4ff2-9492-bd261dbedd48"
        ),
        tags$script(src = "js/umami-events.js")
      )
    ),

    tabPanel("Data",          dataUploadUI("data")),
    tabPanel("Power Analysis",  powerAnalysisUI("pw")),
    tabPanel("Tumor Growth",  tumorGrowthUI("tg")),
    tabPanel("Survival",           survivalUI("surv")),
    tabPanel("Toxicity",            toxicityUI("tox")),
    tabPanel("Dose Response",       doseResponseUI("dr")),
    tabPanel("Drug Synergy",        synergyUI("syn")),
    tabPanel("Log",           analysisLogUI("log")),
    tabPanel("Help",          helpUI()),
    tabPanel("About",          citationUI())
  )

  server <- function(input, output, session) {

    # Initialise structured logging (JSON to stdout + in-process buffer)
    setup_logging_with_buffer()
    log_session(session, "start")
    session$onSessionEnded(function() {
      log_session(session, "end")
    })

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
    toxicity_server("tox", rv)
    synergy_server("syn", rv)
    power_server("pw", rv)
    dose_response_server("dr", rv)
    analysis_log_server("log", rv)

    # Guard: redirect away from placeholder tabs if reached via keyboard
    observeEvent(input$main_nav, {
      log_navigation(session, input$main_nav)
      shinyjs::runjs(sprintf("window.trackTabNavigation('%s');", input$main_nav))
    }, ignoreInit = TRUE)

    # Analysis tabs that require data to be uploaded first
    analysis_tabs <- c("Tumor Growth", "Survival", "Toxicity",
                       "Dose Response", "Drug Synergy")

    # Disable all analysis tabs on startup
    shinyjs::runjs(paste(
      sprintf("window.disableNavTab('%s');", analysis_tabs),
      collapse = "\n"
    ))
    # Power Analysis is always accessible (a priori modes need no data)
    shinyjs::runjs("window.enableNavTab('Power Analysis');")

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

        # Toxicity: only when a weight column is mapped
        wt_mapped <- !is.null(rv$settings$col_weight) &&
                     nzchar(rv$settings$col_weight) &&
                     rv$settings$col_weight != "None"
        if (wt_mapped) {
          shinyjs::runjs("window.enableNavTab('Toxicity');")
        } else {
          shinyjs::runjs("window.disableNavTab('Toxicity');")
        }

        # Drug Synergy: available once data is loaded
        shinyjs::runjs("window.enableNavTab('Drug Synergy');")


        # Dose Response: only when a dose column is mapped
        dose_mapped <- !is.null(rv$settings$col_dose) &&
                       nzchar(rv$settings$col_dose) &&
                       rv$settings$col_dose != "None"
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
    observeEvent(rv$proceed_to_analysis, {
      if (isTRUE(rv$proceed_to_analysis))
        updateNavbarPage(session, "main_nav", selected = "Tumor Growth")
      rv$proceed_to_analysis <- FALSE
    })
  }

  # Enable background async execution for Bayesian models (1 worker to avoid
  # overloading the server; each Bayesian fit runs Stan in parallel internally)
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::multisession, workers = 1L)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  app <- shinyApp(ui = ui, server = server)
  shiny::runApp(app, host = host, port = port, launch.browser = TRUE)
}
