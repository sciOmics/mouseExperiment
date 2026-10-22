# Data Upload + Tumor Growth, for the upload -> fit workflow test (Round 2 I.5).
source(file.path("..", "load-package.R"))
ns <- asNamespace("mouseExperimentDashboard")

shiny::shinyApp(
  ui = shiny::fluidPage(
    shinyjs::useShinyjs(),
    shiny::tabsetPanel(
      id = "main_tabs",
      shiny::tabPanel("Data", ns$dataUploadUI("data")),
      shiny::tabPanel("TG",   ns$tumorGrowthUI("tg"))
    ),
    shiny::textOutput("test_data_valid"),
    shiny::textOutput("test_tg_has_results")
  ),
  server = function(input, output, session) {
    rv <- shiny::reactiveValues(
      uploaded_data = NULL, data_valid = FALSE, detected_cols = NULL,
      schedule_data = NULL, processed_data = NULL, time_col = NULL,
      vol_col = NULL, tx_col = NULL, id_col = NULL, cage_col = NULL,
      dose_col = NULL, tg_results = NULL, tg_running = FALSE,
      tg_col_names = NULL, tg_plot_groups = NULL)
    ns$data_upload_server("data", rv)
    ns$tumor_growth_server("tg", rv)
    output$test_data_valid <- shiny::renderText(as.character(isTRUE(rv$data_valid)))
    # A failed fit also leaves a result (an error list), so require success.
    output$test_tg_has_results <- shiny::renderText(as.character(
      !is.null(rv$tg_results) && !isTRUE(rv$tg_results$error)))
  }
)
