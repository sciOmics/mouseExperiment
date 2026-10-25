# Tumor Growth's "Response & Delay" tab (CODE_REVIEW.md R22-E, F1/F2): each
# animal's mRECIST response with a best-response waterfall, and the tumour
# growth delay to a volume. Both are computed from the data the displayed run
# analysed, with the run's columns, so they follow every run; their own
# settings recompute them without refitting the model.

#' The Response & Delay tab
#' @param ns Namespace function.
#' @keywords internal
tg_endpoints_tab <- function(ns) {
  tabPanel("Response & Delay",
    tab_info(
      strong("Per-animal response and growth delay, "),
      "from the data this run analysed. Neither uses the model.",
      tags$ul(
        tags$li(strong("mRECIST response"), " (Gao et al. 2015): from each animal's % change ",
                "in volume from its first measurement, the best response is the smallest ",
                "change, and the best average response the smallest running mean of the ",
                "changes, both from the day set below. mCR: best < \u221295 % and average < ",
                "\u221240 %; mPR: best < \u221250 % and average < \u221220 %; mSD: best < 35 % and ",
                "average < 30 %; otherwise mPD. The response rate is mCR + mPR among ",
                "evaluable animals."),
        tags$li(strong("Growth delay"), " (T \u2212 C): the Kaplan-Meier median time for an ",
                "arm's tumours to reach the volume below, minus the reference arm's. Each ",
                "animal's time is interpolated on the log scale between measurements; an ",
                "animal that never reached it is censored at its last measurement. ",
                "\"Not reached\" means fewer than half the arm's animals reached it. The ",
                "interval is a bootstrap over animals.")
      )
    ),
    fluidRow(
      column(4, numericInput(ns("tg_resp_min_day"), "Responses count from day (after baseline)",
                             value = 10, min = 0, step = 1)),
      column(4, uiOutput(ns("tg_delay_volume_ui"))),
      column(4, numericInput(ns("tg_delay_n_boot"), "Growth delay interval resamples",
                             value = 1000, min = 0, max = 10000, step = 100))
    ),
    h4("mRECIST response by arm"),
    DT::DTOutput(ns("tg_response_summary")),
    br(),
    fluidRow(
      column(4, radioButtons(ns("tg_waterfall_fill"), "Colour the bars by",
                             choices = c("Arm" = "treatment", "Response" = "category"),
                             inline = TRUE)),
      column(4, numericInput(ns("tg_waterfall_cap"), "Cap the axis at (%)", value = 200,
                             min = 50, step = 50))
    ),
    plotOutput(ns("tg_response_waterfall"), height = "420px"),
    h5("Each animal"),
    DT::DTOutput(ns("tg_response_animals")),
    tags$hr(),
    h4("Growth delay"),
    uiOutput(ns("tg_delay_note")),
    DT::DTOutput(ns("tg_delay_table"))
  )
}

#' Server side of the Response & Delay tab
#'
#' @param input,output,session Module arguments.
#' @param rv Shared reactiveValues.
#' @param colors A reactive of the run's arm colours, or NULL.
#' @return list(responses, delay, waterfall): reactives. `responses()` and
#'   `delay()` return list(result, record), the record for the R script.
#' @keywords internal
add_tg_endpoints <- function(input, output, session, rv, colors = NULL) {
  # The displayed run's data and columns, as recorded (R20-Z).
  run_args <- reactive({
    res <- rv$tg_results
    req(res, !isTRUE(res$error), res$call_record, res$call_record$data)
    a <- res$call_record$args
    list(df = res$call_record$data,
         cols = a[intersect(names(a), c("time_column", "volume_column", "treatment_column",
                                        "id_column", "cage_column"))],
         reference = a$reference_group)
  })

  responses <- reactive({
    ra <- run_args()
    md <- input$tg_resp_min_day
    md <- if (is.numeric(md) && length(md) == 1L && is.finite(md) && md >= 0) md else 10
    args <- c(list(df = ra$df), ra$cols, list(min_day = md))
    res <- tryCatch(do.call(response_categories, args), error = function(e) e)
    list(result = res, record = record_call("response_categories", args))
  })

  output$tg_delay_volume_ui <- renderUI({
    u <- if (identical(rv$settings$volume_units, "cm3")) "cm\u00b3" else "mm\u00b3"
    numericInput(session$ns("tg_delay_volume"), paste0("Growth delay: time to volume (", u, ")"),
                 value = isolate(input$tg_delay_volume) %||%
                   if (identical(rv$settings$volume_units, "cm3")) 1 else 1000,
                 min = 0)
  })

  delay <- reactive({
    ra <- run_args()
    v <- input$tg_delay_volume
    req(is.numeric(v), length(v) == 1L, is.finite(v), v > 0)
    nb <- input$tg_delay_n_boot
    nb <- if (is.numeric(nb) && is.finite(nb) && nb >= 0) as.integer(nb) else 1000L
    args <- c(list(df = ra$df, volume = v), ra$cols,
              list(reference_group = ra$reference, n_boot = nb, seed = 20260730L))
    res <- tryCatch(do.call(growth_delay, args), error = function(e) e)
    list(result = res, record = record_call("growth_delay", args))
  })

  failed <- function(x) inherits(x, "error")

  output$tg_response_summary <- DT::renderDT({
    r <- responses()$result
    if (failed(r)) return(dt_table(data.frame(Message = conditionMessage(r))))
    s <- r$summary
    names(s)[names(s) == "Response_Rate"] <- "Response rate (%)"
    DT::datatable(format_df_for_display(s), rownames = FALSE,
                  options = list(paging = FALSE, searching = FALSE, scrollX = TRUE))
  })

  output$tg_response_animals <- DT::renderDT({
    r <- responses()$result
    if (failed(r)) return(NULL)
    pa <- r$per_animal
    for (col in c("Baseline_Volume", "Best_Response", "Best_Avg_Response")) pa[[col]] <- round(pa[[col]], 1)
    DT::datatable(pa, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  waterfall <- reactive({
    r <- responses()$result
    if (failed(r) || !any(is.finite(r$per_animal$Best_Response))) return(NULL)
    cap <- input$tg_waterfall_cap
    cap <- if (is.numeric(cap) && is.finite(cap) && cap > 0) cap else 200
    fill <- input$tg_waterfall_fill %||% "treatment"
    cols <- if (fill == "treatment" && is.function(colors)) tryCatch(colors(), error = function(e) NULL)
    plot_response_waterfall(r, fill = fill, cap = cap,
                            colors = if (length(cols)) cols)
  })

  # F9: the waterfall as a call on the script's `responses`.
  waterfall_record <- reactive({
    r <- responses()$result
    if (failed(r) || !any(is.finite(r$per_animal$Best_Response))) return(NULL)
    cap <- input$tg_waterfall_cap
    cols <- if ((input$tg_waterfall_fill %||% "treatment") == "treatment" && is.function(colors))
      tryCatch(colors(), error = function(e) NULL)
    record_call("plot_response_waterfall",
                list(responses = NULL, fill = input$tg_waterfall_fill %||% "treatment",
                     colors = if (length(cols)) cols,
                     cap = if (is.numeric(cap) && is.finite(cap) && cap > 0) cap else 200),
                data_arg = NULL, exprs = list(responses = quote(responses)))
  })

  output$tg_response_waterfall <- renderPlot({
    p <- waterfall()
    req(!is.null(p))
    p
  })

  output$tg_delay_note <- renderUI({
    d <- delay()$result
    if (failed(d)) return(status_alert(conditionMessage(d), "warning"))
    tags$p(class = "text-muted", style = "font-size:12px;",
           "Time to ", strong(format(d$volume)), " ", if (identical(rv$settings$volume_units, "cm3"))
             "cm\u00b3" else "mm\u00b3", " against ", strong(d$reference_group), ". ",
           if (d$n_boot >= 2L) sprintf("Intervals: 95%%, from %d bootstrap resamples.", d$n_boot)
           else "No intervals (0 resamples).")
  })

  output$tg_delay_table <- DT::renderDT({
    d <- delay()$result
    if (failed(d)) return(NULL)
    DT::datatable(delay_display_table(d$delay_table, d$reference_group), rownames = FALSE,
                  options = list(paging = FALSE, searching = FALSE, scrollX = TRUE))
  })

  list(responses = responses, delay = delay, waterfall = waterfall,
       waterfall_record = waterfall_record)
}

#' The growth-delay table as shown: "not reached" for infinite medians
#' @param t `growth_delay()$delay_table`.
#' @param reference The reference arm.
#' @keywords internal
delay_display_table <- function(t, reference) {
  dash <- "\u2014"
  days <- function(x) ifelse(is.infinite(x), "not reached", ifelse(is.na(x), dash, sprintf("%.1f", x)))
  is_ref <- t$Treatment %in% reference
  out <- data.frame(
    t$Treatment, t$N, t$N_Reached, days(t$Median_Days),
    ifelse(is_ref, "reference", days(t$Delay_Days)),
    ifelse(is_ref | is.na(t$Delay_Lower) | is.infinite(t$Delay_Lower), dash,
           paste0("[", days(t$Delay_Lower), ", ", days(t$Delay_Upper), "]")),
    ifelse(is_ref, dash, ifelse(is.infinite(t$Delay_Pct), "not reached",
                                ifelse(is.na(t$Delay_Pct), dash, sprintf("%.0f", t$Delay_Pct)))),
    stringsAsFactors = FALSE)
  names(out) <- c("Treatment", "Animals", "Reached the volume", "Median days",
                  "Delay T \u2212 C (days)", "Delay 95% interval", "Delay (% of C)")
  out
}
