# =============================================================================
# Shinytest2 app driver helper
#
# Builds a minimal Shiny app using the dashboard's internal UI/server
# modules — without the jQuery navbar manipulation so the test environment
# stays simple.  This helper is sourced automatically by testthat before
# running test-shiny_app.R.
# =============================================================================

# Build a testable version of the Data-upload module only.
# Returns a shiny::shinyApp object (does not call runApp).
build_test_app <- function() {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyjs")

  pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
  if (is.null(pkg)) return(NULL)

  get_fn <- function(nm) tryCatch(get(nm, envir = pkg), error = function(e) NULL)

  dataUploadUI_fn       <- get_fn("dataUploadUI")
  data_upload_server_fn <- get_fn("data_upload_server")

  if (is.null(dataUploadUI_fn) || is.null(data_upload_server_fn)) return(NULL)

  ui <- shiny::fluidPage(
    shinyjs::useShinyjs(),
    dataUploadUI_fn()
  )

  server <- function(input, output, session) {
    rv <- shiny::reactiveValues(
      uploaded_data = NULL,
      data_valid    = FALSE,
      detected_cols = NULL,
      schedule_data = NULL
    )
    data_upload_server_fn(input, output, session, rv)

    # Expose rv$data_valid for AppDriver$get_value()
    output$test_data_valid <- shiny::renderText(as.character(isTRUE(rv$data_valid)))
    output$test_nrow       <- shiny::renderText(
      if (!is.null(rv$uploaded_data)) as.character(nrow(rv$uploaded_data)) else "0"
    )
  }

  shiny::shinyApp(ui = ui, server = server)
}

#' Build a fuller test app that wires the Data Upload + Tumor Growth modules
#'
#' Used by the I.5 workflow integration test. Loads combo data, fits LME4,
#' asserts results render. Frequentist LME4 on the small combo fixture
#' finishes in < 5s so this is CI-viable.
#'
#' Returns NULL if module functions can't be loaded (graceful skip).
#' @keywords internal
build_workflow_test_app <- function() {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyjs")

  pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
  if (is.null(pkg)) return(NULL)

  get_fn <- function(nm) tryCatch(get(nm, envir = pkg), error = function(e) NULL)

  dataUploadUI_fn       <- get_fn("dataUploadUI")
  data_upload_server_fn <- get_fn("data_upload_server")
  tumorGrowthUI_fn      <- get_fn("tumorGrowthUI")
  tumor_growth_server_fn <- get_fn("tumor_growth_server")

  if (is.null(dataUploadUI_fn) || is.null(data_upload_server_fn) ||
      is.null(tumorGrowthUI_fn) || is.null(tumor_growth_server_fn)) {
    return(NULL)
  }

  ui <- shiny::fluidPage(
    shinyjs::useShinyjs(),
    shiny::tabsetPanel(
      id = "main_tabs",
      shiny::tabPanel("Data",    dataUploadUI_fn()),
      shiny::tabPanel("TG",      tumorGrowthUI_fn(id = "tg"))
    )
  )

  server <- function(input, output, session) {
    rv <- shiny::reactiveValues(
      uploaded_data    = NULL,
      data_valid       = FALSE,
      detected_cols    = NULL,
      schedule_data    = NULL,
      processed_data   = NULL,
      time_col         = NULL,
      vol_col          = NULL,
      tx_col           = NULL,
      id_col           = NULL,
      cage_col         = NULL,
      dose_col         = NULL,
      tg_results       = NULL,
      tg_running       = FALSE,
      tg_col_names     = NULL,
      tg_plot_groups   = NULL
    )
    data_upload_server_fn(input, output, session, rv)
    tumor_growth_server_fn(id = "tg", rv = rv)

    # Expose rv state for AppDriver
    output$test_data_valid <- shiny::renderText(as.character(isTRUE(rv$data_valid)))
    output$test_tg_has_results <- shiny::renderText(
      as.character(!is.null(rv$tg_results))
    )
  }

  shiny::shinyApp(ui = ui, server = server)
}
