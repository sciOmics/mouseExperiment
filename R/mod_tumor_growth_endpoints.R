# Tumor Growth's "Growth Delay" tab (CODE_REVIEW.md R22-E, F2; R23.1, R23.2):
# the tumour growth delay to a volume, with its test and the specific growth
# delay, and each animal's time. It is computed from the data the displayed
# run analysed, with the run's columns, so it follows every run; its own
# settings recompute it without refitting the model. The mRECIST responses
# and their waterfall, which shared the tab, were removed (R23.1).

#' The Growth Delay tab
#' @param ns Namespace function.
#' @keywords internal
tg_endpoints_tab <- function(ns) {
  tabPanel("Growth Delay",
    tab_info(
      strong("Growth delay (T \u2212 C), "),
      "from the data this run analysed, without the model: the Kaplan-Meier median time for an ",
      "arm's tumours to reach the volume below, minus the reference arm's. Each animal's time is ",
      "interpolated on the log scale between measurements; an animal that never reached the volume ",
      "is censored at its last measurement. \"Not reached\" means fewer than half the arm's animals ",
      "reached it.",
      tags$ul(
        tags$li(strong("Delay 95% interval"), " \u2014 a bootstrap over animals within each arm."),
        tags$li(strong("Specific growth delay"), " \u2014 T \u2212 C divided by the reference arm's ",
                "volume doubling time: the delay in doublings, comparable across tumour models that ",
                "grow at different rates."),
        tags$li(strong("p (log-rank)"), " \u2014 each arm's times against the reference arm's, ",
                "corrected over the arms by the method shown. It treats animals as independent, so it ",
                "overstates the evidence when whole cages were randomised.")
      )
    ),
    fluidRow(
      column(4, uiOutput(ns("tg_delay_volume_ui"))),
      column(4, numericInput(ns("tg_delay_n_boot"), "Interval resamples",
                             value = 1000, min = 0, max = 10000, step = 100)),
      column(4, selectInput(ns("tg_delay_p_adjust"), "Multiple-comparison correction",
                            choices = c("Holm" = "holm", "Bonferroni" = "bonferroni",
                                        "FDR (Benjamini-Hochberg)" = "fdr", "None" = "none"),
                            selected = "holm"))
    ),
    uiOutput(ns("tg_delay_note")),
    DT::DTOutput(ns("tg_delay_table")),
    br(),
    h5(strong("Each animal")),
    DT::DTOutput(ns("tg_delay_animals"))
  )
}

#' Server side of the Growth Delay tab
#'
#' @param input,output,session Module arguments.
#' @param rv Shared reactiveValues.
#' @return list(delay): a reactive returning list(result, record), the record
#'   for the R script.
#' @keywords internal
add_tg_endpoints <- function(input, output, session, rv) {
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

  output$tg_delay_volume_ui <- renderUI({
    u <- if (identical(rv$settings$volume_units, "cm3")) "cm\u00b3" else "mm\u00b3"
    numericInput(session$ns("tg_delay_volume"), paste0("Time to volume (", u, ")"),
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
    pm <- input$tg_delay_p_adjust %||% "holm"
    args <- c(list(df = ra$df, volume = v), ra$cols,
              list(reference_group = ra$reference, n_boot = nb, seed = 20260730L,
                   p_adjust_method = pm))
    res <- tryCatch(do.call(growth_delay, args), error = function(e) e)
    list(result = res, record = record_call("growth_delay", args))
  })

  failed <- function(x) inherits(x, "error")

  output$tg_delay_note <- renderUI({
    d <- delay()$result
    if (failed(d)) return(status_alert(conditionMessage(d), "warning"))
    cage <- identical(rv$settings$randomisation_unit, "cage")
    tags$p(class = "text-muted", style = "font-size:12px;",
           "Time to ", strong(format(d$volume)), " ", if (identical(rv$settings$volume_units, "cm3"))
             "cm\u00b3" else "mm\u00b3", " against ", strong(d$reference_group), ". ",
           if (d$n_boot >= 2L) sprintf("Intervals: 95%%, from %d bootstrap resamples. ", d$n_boot)
           else "No intervals (0 resamples). ",
           if (is.finite(d$doubling_time_reference %||% NA))
             sprintf("The reference arm's doubling time: %.1f days. ", d$doubling_time_reference)
           else "The reference arm's doubling time could not be estimated, so there is no specific growth delay. ",
           if (cage) strong("Whole cages were randomised (Study Design): the log-rank test treats animals as independent, so its p-values are too small."))
  })

  output$tg_delay_table <- DT::renderDT({
    d <- delay()$result
    if (failed(d)) return(NULL)
    dt_table(delay_display_table(d$delay_table, d$reference_group), rownames = FALSE,
                  options = list(paging = FALSE, searching = FALSE, scrollX = TRUE))
  })

  output$tg_delay_animals <- DT::renderDT({
    d <- delay()$result
    if (failed(d)) return(NULL)
    dt_table(delay_animal_table(d$per_animal), rownames = FALSE,
             options = list(pageLength = 10, scrollX = TRUE))
  })

  list(delay = delay)
}

#' The growth-delay table as shown: "not reached" for infinite medians
#' @param t `growth_delay()$delay_table`.
#' @param reference The reference arm.
#' @keywords internal
delay_display_table <- function(t, reference) {
  dash <- "\u2014"
  days <- function(x) ifelse(is.infinite(x), "not reached", ifelse(is.na(x), dash, sprintf("%.1f", x)))
  num <- function(x, fmt) ifelse(is.infinite(x), "not reached", ifelse(is.na(x), dash, sprintf(fmt, x)))
  interval <- function(lo, hi, f) ifelse(is.na(lo) | is.infinite(lo), dash, paste0("[", f(lo), ", ", f(hi), "]"))
  pv <- function(x) ifelse(is.na(x), dash, format_pval(x))
  is_ref <- t$Treatment %in% reference
  # A result from before 0.52.0 has no test or specific growth delay.
  col <- function(nm) if (is.null(t[[nm]])) rep(NA_real_, nrow(t)) else t[[nm]]
  sgd <- col("SGD")
  out <- data.frame(
    t$Treatment, t$N, t$N_Reached, days(t$Median_Days),
    ifelse(is_ref, "reference", days(t$Delay_Days)),
    ifelse(is_ref, dash, interval(t$Delay_Lower, t$Delay_Upper, days)),
    ifelse(is_ref, dash, ifelse(is.infinite(t$Delay_Pct), "not reached",
                                ifelse(is.na(t$Delay_Pct), dash, sprintf("%.0f", t$Delay_Pct)))),
    ifelse(is_ref, dash, num(sgd, "%.2f")),
    ifelse(is_ref, dash, interval(col("SGD_Lower"), col("SGD_Upper"), function(x) num(x, "%.2f"))),
    ifelse(is_ref, dash, pv(col("P_Value"))),
    ifelse(is_ref, dash, pv(col("P_Adjusted"))),
    stringsAsFactors = FALSE)
  names(out) <- c("Treatment", "Animals", "Reached the volume", "Median days",
                  "Delay T \u2212 C (days)", "Delay 95% interval", "Delay (% of C)",
                  "Specific growth delay", "SGD 95% interval", "p (log-rank)", "p (adjusted)")
  out
}

#' Each animal's time to the volume, as shown (R23.2)
#' @param pa `growth_delay()$per_animal`.
#' @keywords internal
delay_animal_table <- function(pa) {
  out <- pa
  out$Time <- round(out$Time, 1)
  out$Reached <- ifelse(out$Reached == 1L, "reached", "censored (never reached)")
  names(out)[names(out) == "Time"] <- "Day (reached or last measured)"
  names(out)[names(out) == "Reached"] <- "Outcome"
  out[order(out$Treatment, out[["Day (reached or last measured)"]]), , drop = FALSE]
}
