# The Data Upload module on its own, with two test hooks.
source(file.path("..", "load-package.R"))
ns <- asNamespace("mouseExperiment")

shiny::shinyApp(
  ui = shiny::fluidPage(
    shinyjs::useShinyjs(),
    ns$dataUploadUI("data"),
    shiny::textOutput("test_data_valid"),
    shiny::textOutput("test_nrow")
  ),
  server = function(input, output, session) {
    rv <- shiny::reactiveValues(uploaded_data = NULL, data_valid = FALSE,
                                detected_cols = NULL, schedule_data = NULL)
    ns$data_upload_server("data", rv)
    output$test_data_valid <- shiny::renderText(as.character(isTRUE(rv$data_valid)))
    output$test_nrow <- shiny::renderText(
      if (!is.null(rv$uploaded_data)) as.character(nrow(rv$uploaded_data)) else "0")
  }
)
