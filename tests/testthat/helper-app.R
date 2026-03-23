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
