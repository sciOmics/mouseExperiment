#' Run mouseExperimentDashboard
#'
#' Launch the interactive dashboard for analyzing mouse tumor growth experiments.
#'
#' @param host The host IP address to listen on. Default: "127.0.0.1" (localhost)
#' @param port The port number. Default: 3838
#'
#' @importFrom shiny shinyApp navbarPage tabPanel fluidPage tags observeEvent updateNavbarPage reactiveValues
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch the dashboard
#' runDashboard()
#' 
#' # Or with custom host/port
#' runDashboard(port = 8080)
#' }
runDashboard <- function(host = "127.0.0.1", port = 3838) {
  
  # Ensure shiny is loaded
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required but not installed")
  }
  
  # Ensure constant is available (from global.R)
  if (!exists("APP_VERSION")) {
    APP_VERSION <- "0.1.0 (Development)"
  }
  
  # Note: UI and server functions (homeUI, dataUploadUI, etc.) are already
  # loaded as part of the package namespace, so we don't need to source them
  
  # Define UI
  ui <- navbarPage(
    title = "mouseExperiment Dashboard",
    id = "main_nav",
    
    # Initialize shinyjs for dynamic UI manipulation
    shinyjs::useShinyjs(),
    
    # Home tab
    tabPanel("Home",
      homeUI()
    ),
    
    # Data Upload tab
    tabPanel("Data Upload",
      dataUploadUI()
    ),
    
    # Tumor Growth Analysis tab
    tabPanel("Tumor Growth",
      tumorGrowthUI()
    ),
    
    # Help tab
    tabPanel("Help",
      helpUI()
    )
  )
  
  # Define server logic
  server <- function(input, output, session) {
    
    # Reactive values container
    rv <- reactiveValues()
    
    # Initialize data upload server
    data_upload_server(input, output, session, rv)
    
    # Initialize tumor growth server
    tumor_growth_server(input, output, session, rv)
    
    # Navigation handlers (for buttons on Home tab)
    observeEvent(input$goto_upload, {
      updateNavbarPage(session, "main_nav", selected = "Data Upload")
    })
    
    observeEvent(input$goto_analysis, {
      updateNavbarPage(session, "main_nav", selected = "Tumor Growth")
    })
    
    # Navigation handler for proceed_analysis button on Data Upload tab
    observeEvent(input$proceed_analysis, {
      cat("DEBUG: Proceed to Analysis button clicked\n")
      updateNavbarPage(session, "main_nav", selected = "Tumor Growth")
    })
  }
  
  # Create and run Shiny app
  app <- shinyApp(ui = ui, server = server)
  
  # Run the app
  shiny::runApp(app, host = host, port = port, launch.browser = TRUE)
}
