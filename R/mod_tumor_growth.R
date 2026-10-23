
#' Tumor Growth Analysis Server Logic
#' @keywords internal
tumor_growth_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

  
  # D.7: `rv$tg_results <- NULL` removed here for the same reason as the other
  # five modules — it wiped the slot at module init, so testServer() could never
  # drive this module with a pre-populated result. `rv` is created empty in app.R
  # and the module mounts once, so it never cleared anything in production.
  # The two below are genuine initialisation of state this module owns.
  rv$tg_running    <- FALSE

  ordered_tg_grps <- reactive(
    get_ordered_groups(input, "tg", rv$tg_plot_groups)
  )

  # The reference group is the study's control group (Study Design panel on
  # Data Upload, R20-D.H), in this tab's labels: "Treatment Dose" when a dose
  # column is mapped, so the control maps to its own dose label.
  tg_ref <- reactive({
    req(rv$uploaded_data, rv$data_valid)
    study_reference(rv, analysis_group_labels(rv))
  })

  output$tg_reference_note <- renderUI({
    ref <- tryCatch(tg_ref(), error = function(e) NULL)
    reference_note(rv, ref)
  })

  output$tg_perm_unit_note <- renderUI({
    unit <- rv$settings$randomisation_unit %||% "mouse"
    tags$p(style = "font-size:12px;", "Unit permuted: ",
           strong(if (identical(unit, "cage")) "whole cages" else "individual animals"),
           " (set under Study Design on Data Upload).",
           help_tip(paste(
             "Valid only if this matches how treatment was assigned: if animals were",
             "randomised and then housed by group, the unit is animals. Cage-level",
             "tests have little power; with two cages per arm no effect can reach",
             "p < 0.05. Small designs are enumerated exactly.")))
  })

  # Run analysis button handler
  observeEvent(input$run_tg_analysis, {
    req(rv$uploaded_data, rv$data_valid)
    
    # The reference is the study's control group (R20-D.H).
    ref <- tg_ref()
    if (is.null(ref) || !nzchar(ref)) {
      output$tg_status_msg <- renderUI({
        status_alert(paste("Error: no reference group. Map the treatment column",
                           "and choose the control group under Study Design on",
                           "Data Upload."), "danger")
      })
      shinyjs::runjs("window.trackError('validation', 'Reference group not selected', 'tumor_growth');")
      return()
    }
    rv$tg_reference_used <- ref

    # Track analysis start (no study strings: R20.D12)
    shinyjs::runjs(sprintf(
      "window.trackAnalysisStart('tumor_growth', %s, 'study_control');",
      js_string(input$tg_model_type)
    ))
    
    rv$tg_running <- TRUE
    output$tg_status_msg <- renderUI({
      status_alert("Running analysis...", "info")
    })
    
    # Caliper columns, checked first so the message names the problem.
    if (rv$settings$volume_method == "calculate") {
      len_col <- rv$settings$col_length %||% ""
      wid_col <- rv$settings$col_width  %||% ""
      if (!nzchar(len_col) || !nzchar(wid_col)) {
        output$tg_status_msg <- renderUI({
          status_alert("Please select both Length and Width columns before running analysis.", "danger")
        })
        shinyjs::runjs("window.trackError('validation', 'Length/Width columns not selected', 'tumor_growth');")
        rv$tg_running <- FALSE
        return()
      }
      if (!len_col %in% colnames(rv$uploaded_data) || !wid_col %in% colnames(rv$uploaded_data)) {
        missing_col <- if (!len_col %in% colnames(rv$uploaded_data)) len_col else wid_col
        output$tg_status_msg <- renderUI({
          status_alert(paste0("Column not found in data: '", missing_col, "'"), "danger")
        })
        # R20.D12: no column names to analytics (and none pasted into JS).
        shinyjs::runjs("window.trackError('validation', 'Column not found', 'tumor_growth');")
        rv$tg_running <- FALSE
        return()
      }
    }

    # Prepare data: the subset filter (CODE_REVIEW.md Round 2 F.5) before any
    # column derivation, then study days, then volume, each step recorded for
    # the R script (prepare_upload(), R20-AC).
    # R20.D4: a date error used to escape this observer and leave
    # "Running analysis..." on screen. Report it and stop instead.
    pu <- tryCatch(prepare_upload(rv), error = function(e) {
      log_caught_error("tumor_growth", "date preprocessing", e)
      output$tg_status_msg <- renderUI(status_alert(conditionMessage(e), "danger"))
      rv$tg_running <- FALSE
      NULL
    })
    if (is.null(pu)) return()
    if (length(pu$notes)) {
      showNotification(paste(pu$notes, collapse = " "), type = "warning", duration = 10)
    }
    data     <- pu$data
    time_col <- pu$time_col
    vol_col  <- pu$vol_col
    
    # Store preprocessed data and resolved column names for reuse by other modules
    rv$processed_data  <- data
    rv$processed_steps <- pu$steps
    rv$time_col       <- time_col
    rv$vol_col        <- vol_col
    rv$tx_col         <- rv$settings$col_treatment
    rv$id_col         <- rv$settings$col_id
    rv$cage_col       <- if (!is.null(rv$settings$col_cage)     && rv$settings$col_cage != "")     rv$settings$col_cage     else NULL
    rv$dose_col       <- if (!is.null(rv$settings$col_dose)     && rv$settings$col_dose != "")     rv$settings$col_dose     else NULL
    rv$necrotic_col   <- if (!is.null(rv$settings$col_necrotic) && rv$settings$col_necrotic != "") rv$settings$col_necrotic else NULL

    # When a dose column is mapped, each treatment-dose combination is its own
    # arm (label_arms()). R20-AB: on this run's copy only. It used to replace
    # the shared rv$processed_data and rv$tx_col, so Toxicity and Dose
    # Response saw "TreatmentGroup" arms only after Tumor Growth had run.
    tg_tx_col <- rv$settings$col_treatment
    prep_steps <- pu$steps
    if (!is.null(rv$dose_col)) {
      la_args <- list(df = data, treatment_column = rv$settings$col_treatment,
                      dose_column = rv$dose_col, name = "TreatmentGroup")
      data <- do.call(label_arms, la_args)
      prep_steps <- c(prep_steps, list(record_call("label_arms", la_args, keep_data = FALSE)))
      tg_tx_col <- "TreatmentGroup"
    }

    # Prepare column names
    col_names <- list(
      time_column = time_col,
      volume_column = vol_col,
      id_column = rv$settings$col_id,
      treatment_column = tg_tx_col,
      cage_column = if (rv$settings$col_cage != "") rv$settings$col_cage else NULL
    )

    # Validate all required columns exist in data
    required <- c(Treatment = col_names$treatment_column,
                  ID        = col_names$id_column,
                  Time      = col_names$time_column,
                  Volume    = col_names$volume_column)
    missing_cols <- check_required_cols(data, required)
    if (length(missing_cols) > 0) {
      error_msg <- paste0("Required column(s) not found in data: ",
                   paste(names(missing_cols), "=", sQuote(missing_cols), collapse = ", "))
      output$tg_status_msg <- renderUI({
        status_alert(error_msg, "danger")
      })
      shinyjs::runjs("window.trackError('validation', 'Mapped columns missing', 'tumor_growth');")  # R20.D12
      rv$tg_running <- FALSE
      return()
    }
    
    # ── Bayesian async branch ──────────────────────────────────────────────
    if (isTRUE(input$tg_model_type == "bayes")) {
      pkg_ok <- requireNamespace("future",   quietly = TRUE) &&
                requireNamespace("promises", quietly = TRUE) &&
                requireNamespace("brms",     quietly = TRUE)
      if (!pkg_ok) {
        output$tg_status_msg <- renderUI(
          status_alert(paste(
            "Bayesian mode requires the 'brms', 'future', and 'promises' packages.",
            "Install them with: install.packages(c('brms', 'future', 'promises'))"
          ), "danger")
        )
        rv$tg_running <- FALSE
        return(NULL)
      }

      # Random effects: the Bayesian model offers correlated slopes or an
      # intercept only. "Uncorrelated" maps to correlated slopes (the prior
      # regularises the correlation); "none" to an intercept only.
      re_spec <- switch(input$tg_random_effects %||% "slope",
                        none = "intercept_only",
                        intercept_only = "intercept_only",
                        "slope")

      # Capture all reactive values before entering the future
      # Build manual prior strings from structured UI inputs if needed
      prior_str <- input$tg_prior_strength %||% "skeptical"
      bayes_manual_priors <- if (prior_str == "manual") {
        build_fe_prior_string <- function(dist, mu, sd, df) {
          if (isTRUE(dist == "student_t"))
            sprintf("student_t(%s, %s, %s)", df %||% 3, mu %||% 0, sd %||% 0.25)
          else
            sprintf("normal(%s, %s)", mu %||% 0, sd %||% 0.25)
        }
        build_scale_prior_string <- function(dist, param) {
          if (isTRUE(dist == "half_normal"))
            sprintf("normal(0, %s)", param %||% 0.5)
          else
            sprintf("exponential(%s)", param %||% 2)
        }
        list(
          prior_b         = build_fe_prior_string(
            input$tg_prior_b_dist, input$tg_prior_b_mu,
            input$tg_prior_b_sd,   input$tg_prior_b_df),
          prior_intercept = build_fe_prior_string(
            input$tg_prior_int_dist, input$tg_prior_int_mu,
            input$tg_prior_int_sd,   input$tg_prior_int_df),
          prior_sd        = build_scale_prior_string(
            input$tg_prior_sd_dist,    input$tg_prior_sd_param),
          prior_sigma     = build_scale_prior_string(
            input$tg_prior_sigma_dist, input$tg_prior_sigma_param)
        )
      } else NULL

      bayes_params <- list(
        data              = data,
        col_names         = col_names,
        reference_group   = ref,
        random_effects    = re_spec,
        prior_strength    = prior_str,
        manual_priors     = bayes_manual_priors,
        n_chains          = as.integer(input$tg_bayes_chains  %||% 4L),
        n_warmup          = as.integer(input$tg_bayes_warmup  %||% 1000L),
        n_iter            = as.integer(input$tg_bayes_iter    %||% 500L),
        backend           = input$tg_bayes_backend %||% "rstan",
        extrap_pts        = as.integer(input$tg_extrapolation_points %||% 0L),
        dose_col          = rv$dose_col,
        necrotic_col      = rv$necrotic_col,
        necrotic_handling = tg_necrotic_handling(),
        volume_units      = backend_volume_units(rv),
        plot_groups       = rv$tg_plot_groups,
        prep_steps        = prep_steps,       # R20-AC: for the R script
        source_file       = rv$source_file
      )

      # F.7 — route through the job registry so the fit survives tab
      # disconnect / NPM proxy timeout. Sessions without the identity
      # cookie use a local registry key (R20.D3) and cannot resume.
      launch_tg_bayes_job(rv$session_id, bayes_params,
                          input, output, session, rv)

      return(NULL)   # return immediately — Shiny knows this is non-blocking
    }

    show_analysis_modal("Fitting tumor growth model\u2026", "Results will appear automatically when complete.")
    on.exit(shiny::removeModal(), add = TRUE)

    # R20.D6: the arguments are built once. The fit shown now runs without
    # the randomisation test; when the test is requested, the same call with
    # it runs on the background worker (22 s at 999 permutations, 145 s at
    # 9,999, during which every session on the server used to freeze), and
    # its result joins this run's when it arrives.
    # Built as a call and evaluated, so the R script writes the same call
    # (R20-Z) instead of a deparsed perm_spec object.
    perm_call <- if (isTRUE(input$tg_run_perm)) {
      bquote(perm_spec(unit   = .(rv$settings$randomisation_unit %||% "mouse"),
                       n_perm = .(as.integer(input$tg_perm_n %||% 999L)),
                       seed   = 20260730L))
    } else NULL
    perm <- if (!is.null(perm_call)) eval(perm_call, asNamespace("mouseExperiment")) else NULL
    tg_args <- list(
        df = data,
        time_column = col_names$time_column,
        volume_column = col_names$volume_column,
        id_column = col_names$id_column,
        treatment_column = col_names$treatment_column,
        dose_column = rv$dose_col,
        cage_column = col_names$cage_column,
        model_type = input$tg_model_type,
        transform = input$tg_transform,
        random_effects_specification = input$tg_random_effects,
        handle_cage_effects = input$tg_cage_effect,
        reference_group = ref,  # the study's control group, validated above
        # CODE_REVIEW.md R3.D1 / R3.D3 — pass the user's comparison family and
        # correction to the backend, which now honours both, instead of
        # recomputing contrasts here with a hardcoded adjustment.
        comparison_family = if (identical(input$tg_comparison_type, "tukey")) {
          "all_pairs"
        } else "vs_reference",
        p_adjust_method = local({
          req_adj <- input$tg_p_adjust %||% "dunnett"
          fam <- if (identical(input$tg_comparison_type, "tukey")) {
            "all_pairs"
          } else "vs_reference"
          # Dunnett is only defined for vs-reference and Tukey only for
          # all-pairs; the backend rejects a mismatch, so substitute the closest
          # valid default rather than surfacing an error the user cannot act on.
          # R20.D39: neither is defined on the AUC path (independent t-tests),
          # which the dropdown no longer offers; this covers a stale input.
          if (identical(input$tg_model_type, "auc") &&
              req_adj %in% c("dunnett", "tukey")) "holm"
          else if (req_adj == "dunnett" && fam != "vs_reference") "tukey"
          else if (req_adj == "tukey" && fam != "all_pairs") "dunnett"
          else req_adj
        }),
        # extrapolation_points is deprecated in backend v0.8.0 — the model-based
        # endpoint means extrapolate inside the model and carry the uncertainty.
        # Plot-side extrapolation is unaffected.
        # Backend H.3 -- a perm_spec only when the user opted in and declared
        # the unit; it is added for the background run (R20.D6).
        necrotic_column   = rv$necrotic_col,
        necrotic_handling = tg_necrotic_handling(),
        volume_units      = backend_volume_units(rv),   # backend R20.32 zero floor
        return_model = TRUE,
        include_diagnostics = TRUE,
        plots = FALSE
      )
    result <- safe_analysis(module = "tumor_growth",
                            do.call(mouseExperiment::tumor_growth_statistics, tg_args))
    # R20-Z: the call exactly as made, for the R script in Download Results;
    # R20-AC: with the preparation before it and the uploaded file.
    if (!isTRUE(result$error)) {
      result$call_record <- record_call("tumor_growth_statistics", tg_args)
      result$prep_steps  <- prep_steps
      result$source_file <- rv$source_file
    }
    
    # Generate plots from the results
    if (!isTRUE(result$error)) {
      tryCatch({
        # Plot 1: Tumor growth curves from original data.
        # Enable extrapolation when extrapolation_points > 0.
        # If a survival/censor column is mapped, use it directly.
        # Otherwise synthesize one by flagging mice that disappear before
        # the study's final timepoint.
        if (!is.null(result$growth_rates)) {
          surv_col_mapped <- !is.null(rv$settings$col_survival_event) &&
                             rv$settings$col_survival_event != "" &&
                             rv$settings$col_survival_event != "None" &&
                             rv$settings$col_survival_event %in% colnames(data)
          extra_pts <- if (!is.null(input$tg_extrapolation_points) &&
                           !is.na(input$tg_extrapolation_points))
                         as.integer(input$tg_extrapolation_points) else 0L
          do_extrapolate <- extra_pts > 0L

          # Resolve which data frame + censor column to use
          plot_data  <- data
          surv_col   <- NULL
          if (do_extrapolate) {
            if (surv_col_mapped) {
              surv_col <- rv$settings$col_survival_event
            } else {
              syn <- synthesize_censor_column(
                data     = data,
                id_col   = col_names$id_column,
                time_col = col_names$time_column
              )
              if (syn$n_flagged > 0L) {
                plot_data <- syn$data
                surv_col  <- syn$col
              } else {
                # All mice survived to study end — nothing to extrapolate
                do_extrapolate <- FALSE
              }
            }
          }

          result$tumor_growth_plot <- mouseExperiment::plot_tumor_growth(
            df = plot_data,
            volume_column = col_names$volume_column,
            day_column = col_names$time_column,
            treatment_column = if (!is.null(rv$dose_col)) rv$settings$col_treatment else col_names$treatment_column,
            dose_column = rv$dose_col,
            ID_column = col_names$id_column,
            cage_column = col_names$cage_column,
            survival_column = surv_col,
            extrapolate_volumes = do_extrapolate,
            extrapolation_points = if (extra_pts > 0L) extra_pts else "all",
            group_summary_line = TRUE,
            necrotic_column = rv$necrotic_col
          ) + ggplot2::theme_classic()
        }
        
        # R20.D61: the growth-rate plot built here was never read; the Growth
        # Rates tab draws tg_rates_plot.
        
        # Plot 3: AUC / treatment effects plot
        # When model_type="auc" the result contains auc_analysis$individual with
        # per-subject AUC values — use the package's plot_auc() to match the demo.
        # For lme4 models fall back to a bar chart of adjusted means.
        # R20.4: on the AUC path the plotted quantity is the model-based AUC
        # per arm (treatment_effects, below); the per-animal trapezoids are
        # descriptive and stay in the AUC tab's table.
        model_auc <- isTRUE(result$model_type_used == "auc") &&
          "AUC" %in% names(result$treatment_effects %||% list())
        if (!model_auc && !is.null(result$auc_analysis) &&
            !is.null(result$auc_analysis$individual)) {
          tryCatch({
            extra_pts2 <- if (!is.null(input$tg_extrapolation_points) &&
                              !is.na(input$tg_extrapolation_points))
                            as.integer(input$tg_extrapolation_points) else 0L
            auc_ind <- result$auc_analysis$individual
            # Safety: if the Group column values don't match known treatment groups
            # (can happen when IDs with underscores confuse the package's string-split),
            # re-derive Group by joining on the original ID column.
            known_groups <- unique(as.character(data[[col_names$treatment_column]]))
            if (!all(auc_ind$Group %in% known_groups)) {
              id_tx <- unique(data[, c(col_names$id_column, col_names$treatment_column),
                                   drop = FALSE])
              merged <- merge(auc_ind, id_tx,
                              by.x = "ID", by.y = col_names$id_column,
                              all.x = TRUE)
              new_grp <- merged[[col_names$treatment_column]]
              if (!is.null(new_grp) && !all(is.na(new_grp))) {
                auc_ind$Group    <- as.character(new_grp)
                auc_ind$Treatment <- as.character(new_grp)
                result$auc_analysis$individual <- auc_ind
              }
            }
            result$auc_plot <- mouseExperiment::plot_auc(
              result$auc_analysis$individual,
              show_mean = TRUE,
              error_bar_type = "CI",
              extrapolated_column = if (extra_pts2 > 0L &&
                                       "Extrapolated" %in% colnames(result$auc_analysis$individual))
                                      "Extrapolated" else NULL
            ) + ggplot2::theme_classic()
          }, error = function(e) {
            log_caught_error("tumor_growth", "AUC plot generation", e)
            NULL
          })
        }
        if (is.null(result$auc_plot) && !is.null(result$treatment_effects)) {
          effects_df <- as.data.frame(result$treatment_effects)
          if (nrow(effects_df) > 0) {
            # treatment_effects uses "Group" (lme4) or "Treatment" (AUC) as the label column;
            # do NOT fall back to rownames() since rbind() resets them to "1","2","3"...
            label_col <- intersect(c("Group", "Treatment"), colnames(effects_df))[1]
            # Choose a meaningful numeric column: prefer adjusted mean / mean AUC over
            # count/DF/Reference columns
            preferred <- c("Adjusted_Mean", "AUC", "Mean_AUC", "emmean", "estimate", "mean")
            effect_col <- intersect(preferred, colnames(effects_df))[1]
            if (is.na(effect_col)) {
              # fallback: first numeric that isn't N / DF / Reference
              skip_cols <- c("N", "DF", "Reference")
              num_cols <- colnames(effects_df)[sapply(effects_df, is.numeric)]
              num_cols <- setdiff(num_cols, c(skip_cols, label_col))
              effect_col <- if (length(num_cols) > 0) num_cols[1] else NA
            }
            if (!is.na(label_col) && !is.na(effect_col)) {
              # CODE_REVIEW.md R4.D5 -- read the interval columns the backend
              # says it used. The old sniff list omitted "Lower_CL", which is
              # exactly what treatment_effects contains, so has_ci was always
              # FALSE here and these bars were SE-derived rather than the model's
              # confidence interval.
              ic <- resolve_interval_cols(effects_df, rv$tg_results)
              lower_col <- ic$lower
              upper_col <- ic$upper
              se_col    <- intersect(c("SE", "se", "SEM", "sem"), colnames(effects_df))[1]
              has_ci    <- !is.na(lower_col) && !is.na(upper_col)
              has_se    <- !is.na(se_col)

              plot_df <- data.frame(
                Treatment = as.character(effects_df[[label_col]]),
                Effect    = as.numeric(effects_df[[effect_col]]),
                stringsAsFactors = FALSE
              )
              if (has_ci) {
                plot_df$Lower <- as.numeric(effects_df[[lower_col]])
                plot_df$Upper <- as.numeric(effects_df[[upper_col]])
              } else if (has_se) {
                se_vals <- as.numeric(effects_df[[se_col]])
                df_col  <- intersect(c("df", "DF"), colnames(effects_df))[1]
                t_mult  <- if (!is.na(df_col)) qt(0.975, pmax(as.numeric(effects_df[[df_col]]), 1)) else qnorm(0.975)
                plot_df$Lower <- plot_df$Effect - t_mult * se_vals
                plot_df$Upper <- plot_df$Effect + t_mult * se_vals
              }

              # Build beeswarm using per-subject growth rates when available,
              # overlaid with the LME adjusted mean ± CI from treatment_effects.
              gr_df <- tryCatch(as.data.frame(result$growth_rates), error = function(e) {
                log_caught_error("tumor_growth", "growth_rates coerce to data.frame", e)
                NULL
              })
              grp_col_gr <- if (!is.null(gr_df)) intersect(c("Group","Treatment","TreatmentGroup"), colnames(gr_df))[1] else NA
              val_col_gr <- if (!is.null(gr_df)) intersect(c("growth_rate","GrowthRate"), colnames(gr_df))[1] else NA

              if (!model_auc && !is.null(gr_df) && !is.na(grp_col_gr) && !is.na(val_col_gr)) {
                gr_df$plot_group <- as.character(gr_df[[grp_col_gr]])
                gr_df$plot_val   <- as.numeric(gr_df[[val_col_gr]])
                grp_levels <- as.character(plot_df$Treatment)
                gr_df$plot_group <- factor(gr_df$plot_group, levels = grp_levels)

                ci_half_fn <- function(x) {
                  n <- sum(!is.na(x))
                  qt(0.975, df = max(n - 1L, 1L)) * sd(x, na.rm = TRUE) / sqrt(n)
                }

                p_fallback <- ggplot2::ggplot(gr_df,
                    ggplot2::aes(x = plot_group, y = plot_val, colour = plot_group)) +
                  # Layer 1: 95% CI error bars with horizontal end caps
                  ggplot2::stat_summary(
                    fun.min = function(x) mean(x, na.rm = TRUE) - ci_half_fn(x),
                    fun.max = function(x) mean(x, na.rm = TRUE) + ci_half_fn(x),
                    geom = "errorbar", colour = "black", width = 0.25, linewidth = 0.9
                  ) +
                  # Layer 2: Mean marker as horizontal bar (degenerate crossbar)
                  ggplot2::stat_summary(
                    fun = function(x) mean(x, na.rm = TRUE),
                    fun.min = function(x) mean(x, na.rm = TRUE),
                    fun.max = function(x) mean(x, na.rm = TRUE),
                    geom = "crossbar", colour = "black", width = 0.25, linewidth = 0.9
                  ) +
                  # Layer 3: Individual points with horizontal jitter only
                  ggplot2::geom_point(
                    position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 42),
                    size = 2.5, alpha = 0.85
                  ) +
                  ggplot2::labs(title = "Growth Rate Comparison",
                               x = "Treatment", y = "Growth Rate (log-volume / day)") +
                  ggplot2::theme_classic() +
                  ggplot2::theme(legend.position = "none")
              } else {
                # No individual data — point + errorbar (still not a bar chart)
                p_fallback <- ggplot2::ggplot(plot_df,
                                              ggplot2::aes(x = reorder(Treatment, Effect),
                                                           y = Effect, colour = Treatment)) +
                  ggplot2::geom_point(size = 4) +
                  ggplot2::labs(
                    title = if (model_auc) "Area under each arm's fitted curve" else
                      "Treatment Effects (Adjusted Means)",
                    subtitle = if (model_auc && !is.null(result$auc_window)) sprintf(
                      "Days %s to %s, the last day on which every arm is evaluable; 95%% intervals from model draws",
                      format(result$auc_window[["start"]]), format(result$auc_window[["end"]])),
                    x = "Treatment",
                    y = if (model_auc) "AUC (volume \u00d7 day)" else effect_col) +
                  ggplot2::theme_classic() +
                  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                                legend.position = "none")

                if ("Lower" %in% colnames(plot_df)) {
                  p_fallback <- p_fallback +
                    ggplot2::geom_errorbar(
                      ggplot2::aes(ymin = Lower, ymax = Upper),
                      width = 0.2, colour = "black", linewidth = 0.8
                    )
                }
              }
              result$auc_plot <- p_fallback
            }
          }
        }
        
        # Plot 4: QQ-plot + residuals-vs-fitted for diagnostics tab
        # Diagnostics are inherently on the model scale; the axis labels are
        # tagged so users don't read raw-volume meaning into them.
        diag_scale_tag <- scale_suffix(input$tg_transform %||% "none")
        # R20.D42: the Q-Q plot is the backend's (0.31.0, R20.84), on
        # standardised residuals. This copy drew y = x over raw residuals, so
        # every model looked non-normal (4 % of points near the line).
        # Residuals vs Fitted
        if (!is.null(result$model)) {
          tryCatch({
            fit_vals  <- stats::fitted(result$model)
            resid_vals <- stats::residuals(result$model)
            rf_df <- data.frame(Fitted = fit_vals, Residuals = resid_vals)
            result$diag_resid_fitted_plot <- ggplot2::ggplot(
              rf_df, ggplot2::aes(x = Fitted, y = Residuals)
            ) +
              ggplot2::geom_point(alpha = 0.5) +
              ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
              ggplot2::geom_smooth(method = "loess", se = FALSE,
                                  colour = "steelblue", linewidth = 0.8) +
              ggplot2::theme_classic() +
              ggplot2::labs(title = "Residuals vs Fitted",
                           x = paste0("Fitted Values", diag_scale_tag),
                           y = paste0("Residuals", diag_scale_tag))
          }, error = function(e) {
            log_caught_error("tumor_growth", "residuals-vs-fitted plot", e)
            NULL
          })
        }
        
      }, error = function(e) {
        log_caught_error("tumor_growth", "plot generation", e)
      })
    }
    
    # Generate all comparison types for comparisons tab
    if (!isTRUE(result$error) && !is.null(result$model)) {
      tryCatch({
        # Extract fixed effects terms and treatment comparisons from emmeans
        treatment_col_name <- col_names$treatment_column

        if (!is.null(result$posthoc) && !is.null(result$posthoc$pairwise) &&
            "comparison" %in% colnames(result$posthoc$pairwise)) {
          # --- AUC mode: the backend's model-based AUC comparisons (R20.4) ---
          # Each arm's AUC is the area under its fitted curve over the common
          # evaluable window; comparisons are ratios, with the difference in
          # volume x day alongside, intervals and Wald p-values from model draws.
          pw <- result$posthoc$pairwise
          # Reformat to emmeans-compatible column names for downstream display/plot code
          # R20.D43: the intervals are adjusted like p_adjusted (backend
          # 0.31.0, R20.85); the unadjusted 95 % ones sat beside adjusted
          # p-values, so an interval could clear zero next to "not significant".
          pw_df <- data.frame(
            contrast    = pw$comparison,
            estimate    = pw$estimate,
            lower.CL    = pw$CI_Adjusted_Lower %||% pw$ci_lower,
            upper.CL    = pw$CI_Adjusted_Upper %||% pw$ci_upper,
            ratio       = pw$ratio,
            ratio.lower = pw$Ratio_CI_Adjusted_Lower %||% pw$ratio_lower,
            ratio.upper = pw$Ratio_CI_Adjusted_Upper %||% pw$ratio_upper,
            z.ratio     = pw$z_value,
            p.value     = pw$p_adjusted,
            CI_Adjustment = pw$CI_Adjustment %||% NA_character_,
            stringsAsFactors = FALSE
          )
          result$all_pairwise_comparisons <- pw_df

          ref_group <- if (!is.null(ref) && nzchar(ref))
            ref else {
              lvls <- sort(unique(as.character(
                if (!is.null(result$treatment_effects$Treatment))
                  result$treatment_effects$Treatment
                else
                  data[[col_names$treatment_column]]
              )))
              if (length(lvls) > 0) lvls[1] else NULL
            }

          if (!is.null(ref_group)) {
            is_ref <- grepl(ref_group, pw_df$contrast, fixed = TRUE)
            ref_comp_df <- if (any(is_ref)) pw_df[is_ref, ] else pw_df
            
            # Standardize contrast ordering: ensure reference is always on the RIGHT side
            # If contrast is "Control - Treatment", swap to "Treatment - Control" and negate estimate
            for (i in seq_len(nrow(ref_comp_df))) {
              parts <- strsplit(ref_comp_df$contrast[i], " - ", fixed = TRUE)[[1]]
              if (length(parts) == 2) {
                left_part <- trimws(parts[1])
                right_part <- trimws(parts[2])
                
                # If reference is on the left, swap and negate
                if (identical(left_part, ref_group) && !identical(right_part, ref_group)) {
                  ref_comp_df$contrast[i] <- paste0(right_part, " - ", left_part)
                  ref_comp_df$estimate[i] <- -ref_comp_df$estimate[i]
                  # Swap CI bounds and negate
                  temp_lower <- ref_comp_df$lower.CL[i]
                  ref_comp_df$lower.CL[i] <- -ref_comp_df$upper.CL[i]
                  ref_comp_df$upper.CL[i] <- -temp_lower
                  if (!is.na(ref_comp_df$z.ratio[i])) {
                    ref_comp_df$z.ratio[i] <- -ref_comp_df$z.ratio[i]
                  }
                  r_lo <- ref_comp_df$ratio.lower[i]
                  ref_comp_df$ratio[i]       <- 1 / ref_comp_df$ratio[i]
                  ref_comp_df$ratio.lower[i] <- 1 / ref_comp_df$ratio.upper[i]
                  ref_comp_df$ratio.upper[i] <- 1 / r_lo
                }
              }
            }
            result$reference_comparisons_dunnett <- ref_comp_df
          } else {
            result$reference_comparisons_dunnett <- pw_df
          }

        } else {
          # CODE_REVIEW.md R3.D1 / R3.D2 / R3.D4 — the dashboard used to discard
          # result$pairwise_comparisons and rebuild contrasts here from
          # result$model. Three problems, all now fixed upstream:
          #
          #  * The backend applied no adjustment at all (backend R3.1) while this
          #    block applied Tukey/Dunnett, so the same contrast carried one
          #    p-value in the UI and a different one in an R session. As of
          #    backend v0.6.0 the package applies the user's chosen
          #    comparison_family + p_adjust_method itself, so there is exactly one
          #    implementation and the two agree by construction.
          #  * adjust = "dunnett" is not a value emmeans recognises; it silently
          #    downgraded to the "dunnettx" approximation while the help text
          #    claimed exact Dunnett. The backend now maps "dunnett" to the exact
          #    multivariate-t ("mvt").
          #  * emmeans cannot dispatch on a gamm4 object, so for model_type =
          #    "gam" this errored into log_caught_error() and the "All
          #    Comparisons" tab could never populate. Rendering the backend's
          #    table removes that failure mode entirely.
          pw_backend <- result$pairwise_comparisons
          if (!is.null(pw_backend)) {
            pw_df <- as.data.frame(
              if (inherits(pw_backend, "emmGrid")) {
                summary(pw_backend, infer = c(TRUE, TRUE))
              } else pw_backend
            )
            if ("p_value" %in% names(pw_df) && !"p.value" %in% names(pw_df)) {
              pw_df$p.value <- pw_df$p_value
            }
            fam <- result$comparison_family %||% "vs_reference"
            if (identical(fam, "all_pairs")) {
              result$all_pairwise_comparisons  <- pw_df
            } else {
              result$reference_comparisons_dunnett <- pw_df
            }
          }
        }

        # R20.D61: a second forest plot was built here on every run and never
        # shown; the displayed one is tg_forest_plot (mod_tumor_growth_plots.R).

      }, error = function(e) {
        log_caught_error("tumor_growth", "comparison generation", e)
      })
    }

    # Persist what the backend actually did, so renderDTs condition on the
    # outcome rather than the request (CODE_REVIEW.md R3.D6). These diverge: the
    # AUC path ignores `transform` and computes on raw volumes, so an AUC run
    # submitted with the default transform = "log" reports transform_used
    # "none". Recording the input would label that output "log" and mis-scale
    # any consumer not gated on model_type_used == "lme4".
    #
    # $meta is the backend's own declaration (B7.1/B7.2); the top-level field and
    # then the input are fallbacks for a backend older than that.
    result$transform_used  <- result$meta$transform_used %||%
                              result$transform_used      %||%
                              input$tg_transform         %||% "none"
    result$model_type_used <- result$meta$model_type_used %||%
                              result$model_type_used      %||%
                              input$tg_model_type         %||% "lme4"

    # R20.D46: what this run was asked for, for the report.
    result$run_settings <- list(
      model_type      = input$tg_model_type     %||% "",
      transform       = input$tg_transform      %||% "",
      random_effects  = input$tg_random_effects %||% "",
      reference_group = ref,
      cage_effect     = input$tg_cage_effect    %||% "",
      necrotic_handling = tg_necrotic_handling(),
      n_rows          = nrow(data))

    run_token <- paste0(format(Sys.time(), "%H%M%OS6"), sample.int(1e6, 1))
    result$run_token <- run_token
    # R20.D6: the randomisation test's call (run below), recorded for the
    # R script; the script includes it once the test's result has arrived.
    if (!is.null(perm) && !isTRUE(result$error)) {
      perm_args <- tg_args
      perm_args[c("permutation_test", "return_model", "include_diagnostics")] <-
        list(perm, FALSE, FALSE)
      result$perm_call_record <- record_call("tumor_growth_statistics", perm_args,
                                             exprs = list(permutation_test = perm_call))
    }
    rv$tg_results <- result
    updateTabsetPanel(session, "tg_results_tabs", selected = "Summary")
    rv$tg_running <- FALSE

    # R20.D6: the randomisation test, on the background worker.
    rv$tg_perm_status <- NULL
    if (!is.null(perm) && !isTRUE(result$error)) {
      launched <- launch_background_job(
        kind = "tg_perm",
        task = background_task(.tg_perm_only, perm_args),
        context = list(run_token = run_token),
        input = input, output = output, session = session, rv = rv,
        job_id_field = "tg_perm_job_id",
        finalize_success      = .finalize_perm_success,
        finalize_failure      = .finalize_perm_failure,
        finalize_cancellation = .finalize_perm_cancellation,
        finalize_lost         = .finalize_perm_lost)
      rv$tg_perm_status <- if (is.null(launched)) "busy" else "running"
    }

    # Track analysis completion with Umami
    if (isTRUE(result$error)) {
      shinyjs::runjs(sprintf("window.trackAnalysisComplete('tumor_growth', false, %s);",
                             js_string(result$message)))
    } else {
      shinyjs::runjs("window.trackAnalysisComplete('tumor_growth', true);")
    }

    # Store the data and column names needed for reactive plot re-rendering
    rv$tg_plot_data   <- data
    rv$tg_col_names   <- col_names
    rv$tg_plot_groups <- if (!is.null(rv$dose_col) && rv$dose_col %in% colnames(data) &&
                              !is.null(rv$settings$col_treatment) && rv$settings$col_treatment %in% colnames(data)) {
      # plot_tumor_growth builds Group as paste(treatment, dose, sep = " - Dose: ")
      # so colour-picker names must match that format exactly
      unique(paste(data[[rv$settings$col_treatment]], data[[rv$dose_col]], sep = " - Dose: "))
    } else if (!is.null(data) && !is.null(col_names$treatment_column) &&
                col_names$treatment_column %in% colnames(data)) {
      unique(as.character(data[[col_names$treatment_column]]))
    } else NULL
    
    # Update status
    if (isTRUE(result$error)) {
      output$tg_status_msg <- renderUI({
        status_alert(paste("Error:", result$message), "danger")
      })
    } else {
      # R20.D34: the backend's warnings are shown, not dropped.
      output$tg_status_msg <- renderUI(completion_status(list(result)))
    }
  })
  
  # Dynamic labels for manual scale prior parameter inputs
  output$tg_prior_sd_param_label <- renderUI({
    if (isTRUE(input$tg_prior_sd_dist == "half_normal")) "Scale \u03c3" else "Rate"
  })
  output$tg_prior_sigma_param_label <- renderUI({
    if (isTRUE(input$tg_prior_sigma_dist == "half_normal")) "Scale \u03c3" else "Rate"
  })

  # R20-D.V: how necrotic observations are handled. The selector used to appear
  # only once a run had set rv$necrotic_col, and until then the run fell back to
  # "none": the first analysis kept necrotic observations, then the selector
  # appeared showing "Exclude", and the next run excluded them. It now appears
  # as soon as the column is mapped, and the fallback is its displayed default.
  necrotic_mapped <- reactive({
    col <- rv$settings$col_necrotic
    is.character(col) && length(col) == 1L && nzchar(col) && !identical(col, "None")
  })
  tg_necrotic_handling <- function() {
    if (is.null(rv$necrotic_col)) "none" else input$tg_necrotic_handling %||% "exclude"
  }

  output$necrotic_handling_ui <- renderUI({
    if (!isTRUE(necrotic_mapped())) return(NULL)
    tagList(
      hr(),
      h5(strong("Necrotic Observations"), help_tip(paste(
        "Exclude removes flagged observations before fitting. Covariate adds a 0/1",
        "Necrotic term to the fixed effects. No adjustment ignores the flag and keeps",
        "every observation."))),
      selectInput(session$ns("tg_necrotic_handling"), NULL,
        choices = c(
          "Exclude from analysis"    = "exclude",
          "Include as covariate"     = "covariate",
          "No adjustment (raw data)" = "none"
        ),
        selected = "exclude"
      )
    )
  })

  # Summary output
  output$tg_summary <- renderText({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) {
      return(paste("Error:", rv$tg_results$message))
    }
    
    # Build a summary text from the summary object
    summary_list <- rv$tg_results$summary
    
    # Format summary information
    summary_lines <- c(
      "=== TUMOR GROWTH ANALYSIS SUMMARY ===",
      "",
      "MODEL SPECIFICATION:",
      paste("  Model Type:", summary_list$analysis_type),
      paste("  Transformation:", paste(summary_list$methods$volume_transformation)),
      zero_floor_line(rv$tg_results),
      "",
      "DATA DESCRIPTION:",
      paste("  Subjects:", summary_list$data_description$subjects),
      paste("  Treatment Groups:", summary_list$data_description$treatment_groups),
      paste("  Time Points:", summary_list$data_description$time_points),
      paste("  Reference Group:", summary_list$data_description$reference_group)
    )

    # Model specification details only available for LME4
    if (!is.null(summary_list$model_specification)) {
      summary_lines <- c(summary_lines,
        "",
        "MODEL DETAILS:",
        paste("  Fixed Effects:", summary_list$model_specification$fixed_effects),
        paste("  Random Effects:", summary_list$model_specification$random_effects),
        paste("  Tests:", rv$tg_results$anova_method %||% "Type III")
      )
      # R20.5 / R20-K: say when the requested random effects could not be fitted.
      re <- rv$tg_results$random_effects
      if (!is.null(re) && !identical(re$used, re$requested)) {
        summary_lines <- c(summary_lines, paste0(
          "  NOTE: ", re$requested, " could not be fitted (",
          paste(re$fallback, collapse = "; "), "); the model uses ", re$used, "."))
      }
    }

    # AUC (R20.4): the model and the window behind it.
    if (identical(rv$tg_results$model_type_used, "auc") && !is.null(rv$tg_results$auc_window)) {
      w <- rv$tg_results$auc_window
      summary_lines <- c(summary_lines,
        "",
        "AUC:",
        paste0("  Window: day ", format(w[["start"]]), " to day ", format(w[["end"]]),
               ", the last day on which every arm is evaluable"),
        paste("  Model:", rv$tg_results$endpoint_model$formula %||% ""),
        paste("  Rule:", rv$tg_results$evaluability$rule %||% ""))
    }

    # CODE_REVIEW.md R3.17 / G.2 — surface the detected cage design, where cage
    # ended up, and the cage ICC. The ICC is what tells a reader how much the
    # clustering mattered: the design effect on effective sample size is
    # 1 + (mice per cage - 1) x ICC, so an ICC of 0.2 with 4 mice per cage means
    # the effective n is about 1.6x smaller than the animal count suggests.
    cage_struct <- summary_list$model_specification$cage_structure
    if (!is.null(cage_struct)) {
      icc_df <- rv$tg_results$cage_analysis$icc
      icc_line <- if (!is.null(icc_df) && is.finite(icc_df$ICC[1])) {
        paste0("  Cage ICC: ", round(icc_df$ICC[1], 3),
               "  (variance ", signif(icc_df$Variance[1], 3), " of ",
               signif(icc_df$Total_Variance[1], 3), " total)")
      } else {
        "  Cage ICC: not estimated (cage is not in the model)"
      }
      summary_lines <- c(summary_lines,
        "",
        "CAGE STRUCTURE:",
        paste("  Detected design:", cage_struct),
        paste("  Cage entered as:",
              summary_list$model_specification$cage_placement %||% "(unknown)"),
        icc_line,
        paste("  Reason:", summary_list$model_specification$cage_reason %||% "")
      )
    }

    # Multiplicity provenance — state what was applied rather than leaving the
    # reader to assume a default (CODE_REVIEW.md R3.D1 / R3.D7).
    if (!is.null(rv$tg_results$comparison_family)) {
      summary_lines <- c(summary_lines,
        "",
        "COMPARISONS:",
        paste("  Family:", rv$tg_results$comparison_family),
        paste("  Correction:", rv$tg_results$p_adjust_method_used %||% "(none)")
      )
    }

    summary_lines <- c(summary_lines,
      "",
      "NOTES:",
      paste("  ", paste(summary_list$notes, collapse="\n  "))
    )

    elapsed_str <- format_elapsed(rv$tg_results$elapsed_seconds)
    if (!is.null(elapsed_str)) {
      summary_lines <- c(summary_lines, "",
        paste("COMPUTATION TIME:", elapsed_str))
    }

    paste(summary_lines, collapse = "\n")
  })

  # Update plot selector choices based on model type (AUC plot only in AUC mode;

  # forest plot only in LME4 mode)
  observeEvent(input$tg_model_type, {
    req(input$tg_model_type)
    shinyjs::runjs(sprintf("window.trackEvent('model_type_selected', {analysis_type: 'tumor_growth', method: %s});",
                           js_string(input$tg_model_type)))  # R20.D22
  }, ignoreInit = TRUE)

  observeEvent(input$tg_model_type, {
    base_choices <- c(
      "Tumor Growth Curves"       = "growth_curves",
      "Growth Rates Plot"         = "growth_rate",
      "Waterfall Plot"            = "waterfall"
    )
    if (isTRUE(input$tg_model_type == "auc")) {
      base_choices <- append(
        base_choices,
        c("Area Under the Curve Plot" = "auc"),
        after = 2L  # insert after "Growth Rates Plot"
      )
    } else {
      base_choices <- c(base_choices, "Effect Size Forest Plot" = "effect_forest")
    }
    current <- isolate(input$tg_plot_selector) %||% "growth_curves"
    selected <- if (current %in% base_choices) current else "growth_curves"
    updateSelectInput(session, "tg_plot_selector",
      choices  = base_choices,
      selected = selected
    )
  }, ignoreNULL = FALSE)

  # R20.D39: offer only corrections the selected model supports. Holm replaces a
  # Dunnett/Tukey selection on the way into AUC, and Dunnett comes back on the
  # way out if the user has not changed that Holm since.
  p_adjust_auto_holm <- reactiveVal(FALSE)
  observeEvent(input$tg_model_type, {
    mt      <- input$tg_model_type %||% "lme4"
    choices <- tg_p_adjust_choices(mt)
    current <- isolate(input$tg_p_adjust) %||% "dunnett"
    selected <- if (!current %in% choices) {
      p_adjust_auto_holm(TRUE)
      "holm"
    } else if (mt != "auc" && isTRUE(p_adjust_auto_holm()) && current == "holm") {
      p_adjust_auto_holm(FALSE)
      "dunnett"
    } else current
    updateSelectInput(session, "tg_p_adjust", choices = choices, selected = selected)
  }, ignoreInit = TRUE)

  # Dynamically insert / remove AUC tabs based on model type
  auc_tabs_inserted <- reactiveVal(FALSE)

  observeEvent(input$tg_model_type, {
    tabs_id <- "tg_results_tabs"
    if (isTRUE(input$tg_model_type == "auc") && !isTRUE(auc_tabs_inserted())) {
      insertTab(inputId = tabs_id, target = "Growth Rates", position = "after",
        tabPanel("AUC Individual Data",
          tab_info(
            strong("Per-animal area under the growth curve (AUC), for description only. "), "Each animal's trapezoidal AUC over its own measurements. An animal removed early has less area because it was measured for less time, so these values are not compared between arms; the comparisons use the model-based AUC over a common window.",
            tags$ul(
              tags$li(strong("ID"), " \u2014 individual animal identifier."),
              tags$li(strong("Treatment"), " \u2014 treatment group."),
              tags$li(strong("Cage"), " \u2014 cage identifier."),
              tags$li(strong("Group"), " \u2014 group assignment."),
              tags$li(strong("AUC"), " \u2014 area under the volume-time curve (mm\u00b3 \u00d7 days)."),
              tags$li(strong("Last_Day"), " \u2014 final measurement timepoint for this animal."),
              tags$li(strong("First_Day"), " \u2014 initial measurement timepoint for this animal."),
              tags$li(strong("Extrapolated"), " \u2014 logical flag indicating if AUC was extended beyond observed data."),
              tags$li(strong("NumPoints"), " \u2014 number of observed volume measurements used to compute AUC.")
            )
          ),
          shinycssloaders::withSpinner(
            DT::DTOutput(session$ns("tg_auc_individual_table"))
          ),
          downloadButton(session$ns("download_tg_auc_individual"), "Download AUC Individual Data")
        )
      )
      insertTab(inputId = tabs_id, target = "AUC Individual Data", position = "after",
        tabPanel("AUC Summary",
          tab_info(
            strong("Group summaries of the per-animal AUCs, for description only"), " (mean, standard deviation, sample size and standard error). They inherit the per-animal AUC's dependence on follow-up time; the model-based AUC by arm is in the Treatment Effects tab.",
            tags$ul(
              tags$li(strong("Treatment"), " \u2014 treatment group name."),
              tags$li(strong("Mean_AUC"), " \u2014 group mean area under curve (mm\u00b3 \u00d7 days)."),
              tags$li(strong("SD_AUC"), " \u2014 standard deviation of AUC within the group."),
              tags$li(strong("N"), " \u2014 number of animals in the group with valid AUC data."),
              tags$li(strong("SEM_AUC"), " \u2014 standard error of the mean AUC (SD / \u221aN).")
            ),
            em("Use these to check data and sample sizes, not to compare arms.")
          ),
          shinycssloaders::withSpinner(
            DT::DTOutput(session$ns("tg_auc_summary_table"))
          ),
          downloadButton(session$ns("download_tg_auc_summary"), "Download AUC Summary")
        )
      )
      auc_tabs_inserted(TRUE)
    } else if (!isTRUE(input$tg_model_type == "auc") && isTRUE(auc_tabs_inserted())) {
      removeTab(inputId = tabs_id, target = "AUC Individual Data")
      removeTab(inputId = tabs_id, target = "AUC Summary")
      auc_tabs_inserted(FALSE)
    }
  }, ignoreNULL = FALSE)

  # Dynamically insert / remove Bayesian tabs based on model type
  bayes_tabs_inserted <- reactiveVal(FALSE)

  observeEvent(input$tg_model_type, {
    tabs_id <- "tg_results_tabs"
    if (isTRUE(input$tg_model_type == "bayes") && !isTRUE(bayes_tabs_inserted())) {
      appendTab(inputId = tabs_id,
        tabPanel("Posterior Check",
          tab_info(
            strong("Posterior predictive check (PPC)"), " \u2014 the model generates new datasets from the posterior and overlays their densities on the observed data. Good model fit means the simulated densities closely envelope the observed density.",
            tags$ul(
              tags$li(strong("y"), " \u2014 observed (transformed) tumor volume values."),
              tags$li(strong("y_rep"), " \u2014 50 simulated datasets drawn from the posterior predictive distribution.")
            ),
            em("If the observed density falls well outside the simulated envelope, consider a different transformation or random-effects specification.")
          ),
          shinycssloaders::withSpinner(
            plotOutput(session$ns("tg_pp_check_plot"), height = "450px")
          ),
          downloadButton(session$ns("download_tg_pp_check"), "Download Plot"),
          # R20.D59: rendered, but never placed anywhere.
          tags$hr(),
          h4("Posterior summary"),
          p(class = "text-muted", "Posterior medians and 95% credible intervals of the model's parameters."),
          DT::DTOutput(session$ns("tg_posterior_summary"))
        ),
        select = FALSE
      )
      bayes_tabs_inserted(TRUE)
    } else if (!isTRUE(input$tg_model_type == "bayes") && isTRUE(bayes_tabs_inserted())) {
      removeTab(inputId = tabs_id, target = "Posterior Check")
      bayes_tabs_inserted(FALSE)
    }
  }, ignoreNULL = FALSE)

  # Dynamically insert / remove Necrosis results tab when necrosis_summary is available
  necrosis_tab_inserted <- reactiveVal(FALSE)

  observe({
    has_necrotic <- !is.null(rv$tg_results) && !isTRUE(rv$tg_results$error) &&
                   !is.null(rv$tg_results$necrosis_summary)
    if (isTRUE(has_necrotic) && !isTRUE(necrosis_tab_inserted())) {
      insertTab(inputId = "tg_results_tabs",
        target = "ANOVA", position = "before",
        tabPanel("Necrosis",
          tab_info(
            strong("Per-group necrosis summary"), " \u2014 counts and percentages of necrotic
            observations, number of affected mice, and the first day necrosis was observed.
            The summary always reflects the full uploaded dataset, regardless of which
            handling mode was selected.",
            tags$ul(
              tags$li(strong("N_Observations"), " \u2014 total measurement rows for this group."),
              tags$li(strong("N_Necrotic"), " \u2014 number flagged as necrotic."),
              tags$li(strong("Pct_Necrotic"), " \u2014 percentage of observations that were necrotic."),
              tags$li(strong("N_Mice_Total"), " \u2014 total animals in this group."),
              tags$li(strong("N_Mice_With_Necrosis"), " \u2014 animals with at least one necrotic observation."),
              tags$li(strong("First_Necrotic_Day"), " \u2014 earliest day a necrotic observation was recorded.")
            )
          ),
          shinycssloaders::withSpinner(
            DT::DTOutput(session$ns("tg_necrosis_table"))
          ),
          downloadButton(session$ns("download_tg_necrosis"), "Download Necrosis Summary")
        ),
        select = FALSE
      )
      necrosis_tab_inserted(TRUE)
    } else if (!isTRUE(has_necrotic) && isTRUE(necrosis_tab_inserted())) {
      removeTab(inputId = "tg_results_tabs", target = "Necrosis")
      necrosis_tab_inserted(FALSE)
    }
  })

  # Track plot type switching
  observeEvent(input$tg_plot_selector, {
    req(input$tg_plot_selector)
    shinyjs::runjs(sprintf("window.trackEvent('plot_type_selected', {analysis_type: 'tumor_growth', plot_type: %s});",
                           js_string(input$tg_plot_selector)))  # R20.D22
  }, ignoreInit = TRUE)

  # ANOVA table
  output$tg_download_plot_btn <- renderUI({
    btn_id <- switch(input$tg_plot_selector %||% "growth_curves",
      growth_curves  = "download_tg_plot_growth",
      growth_rate    = "download_tg_plot_rates",
      auc            = "download_tg_plot_auc",
      waterfall      = "download_tg_plot_waterfall",
      effect_forest  = "download_tg_plot_effect_forest",
      "download_tg_plot_growth"
    )
    downloadButton(session$ns(btn_id), "Download")
  })

  output$tg_summary_tab_info <- renderUI({
    model_type <- input$tg_model_type %||% "lme4"
    if (model_type == "auc") {
      tab_info(
        strong("AUC analysis summary: "), "Reports the analysis type (AUC), number of subjects and treatment groups, reference group, the model behind the AUC, and its window: from the first study day to the last day on which every arm is evaluable."
      )
    } else {
      tab_info(
        strong("Model summary: "), "Reports the fitted model \u2014 model type (LME4 / AUC / Bayesian), data transformation applied, random effects structure (subject and cage), number of subjects and treatment groups, reference group, and any convergence warnings or singular-fit notes."
      )
    }
  })

  output$tg_comparisons_tab_info <- renderUI({
    # R20.D44: describe the run on screen, not the sidebar's current values.
    res <- if (!isTRUE(rv$tg_results$error)) rv$tg_results else NULL
    model_type <- res$model_type_used %||% input$tg_model_type %||% "lme4"
    if (model_type == "auc") {
      tab_info(
        strong("Comparisons of the model-based AUCs "), "(the area under each arm's fitted curve over the common evaluable window), with the correction selected in the sidebar.",
        tags$ul(
          tags$li(strong("contrast"), " \u2014 the two groups being compared (\u2018Group A \u2212 Group B\u2019)."),
          tags$li(strong("ratio"), " \u2014 AUC of the first group divided by the second's. Below 1 means less total tumour burden; ratio.lower / ratio.upper are its interval."),
          tags$li(strong("estimate"), " \u2014 the difference in AUC (volume \u00d7 days); lower.CL / upper.CL are its interval."),
          tags$li(strong("z.ratio"), " \u2014 Wald statistic for the log ratio."),
          tags$li(strong("p.value"), " \u2014 ", p_value_help(res), " It tests the ratio."),
          tags$li("The intervals are adjusted like the p-values where the correction has an interval (", em("CI_Adjustment"), " says which); the table's caption names the correction.")
        ),
        em("Intervals and p-values come from draws of the fitted model's coefficients.")
      )
    } else {
      tab_info(
        strong("Post-hoc pairwise contrasts "), "between treatment groups, computed from estimated marginal means (emmeans). Use the Analysis Parameters dropdown to switch between comparison types.",
        tags$p(
          strong("What is being compared:"), " these contrasts compare model-predicted ",
          strong("tumor volume"), " between groups \u2014 not growth rate. Each group\u2019s marginal mean is evaluated at the ",
          em("average study day"), ", so ", strong("estimate"), " answers: \u201cAt the typical timepoint, how much larger or smaller is one group\u2019s predicted volume relative to another?\u201d",
          " To test whether ", em("growth trajectories diverge"), " across groups, see the first (treatment \u00d7 time) row in the ANOVA tab. Per-animal growth rates are in the Growth Rates tab."
        ),
        tags$ul(
          tags$li(strong("contrast"), " \u2014 the two groups being compared (\u2018Treatment \u2212 Reference\u2019 for Dunnett, \u2018A \u2212 B\u2019 for Tukey)."),
          tags$li(strong("df"), " \u2014 Satterthwaite or Kenward\u2013Roger approximate degrees of freedom used for the t-test."),
          tags$li(strong("estimate"), " \u2014 mean difference in tumor volume between the two groups on the analysis scale (log or sqrt if a transformation was applied). Negative values mean the first group has lower volume."),
          tags$li(strong("lower.CL / upper.CL"), " \u2014 confidence interval for the mean difference on the analysis scale, adjusted like the p-values (", em("CI_Adjustment"), " says how; Holm and FDR have no interval of their own, so theirs are Bonferroni-level)."),
          tags$li(strong("SE"), " \u2014 standard error of the estimated difference."),
          tags$li(strong("Fold_Change"), " (if shown; log-transform models only) \u2014 exp(estimate): the multiplicative ratio of group geometric means on the original volume scale. Values below 1.0 indicate the treatment group has lower volume than reference. Fold_Change_Lower / Fold_Change_Upper give the back-transformed CI bounds."),
          tags$li(strong("t.ratio"), " \u2014 t-statistic (estimate \u00f7 SE). Large absolute values indicate a larger difference relative to its uncertainty."),
          tags$li(strong("p.value"), " \u2014 ", p_value_help(res)),
          tags$li(strong("Lower_CrI / Upper_CrI"), " (Bayesian path only) \u2014 2.5\u202f% and 97.5\u202f% bounds of the posterior distribution of the contrast on the model scale."),
          tags$li(strong("P_direction"), " (Bayesian path only) \u2014 posterior probability that the contrast is in its dominant direction, > 0.5 by construction. P_direction = 0.97 means 97\u202f% of posterior mass agrees on the sign of the effect; equivalent to a one-sided posterior probability and replaces the frequentist \u2018does the CrI exclude zero?\u2019 dichotomy.")
        )
      )
    }
  })

  output$tg_treatment_effects_tab_info <- renderUI({
    model_type <- input$tg_model_type %||% "lme4"
    if (model_type == "auc") {
      tab_info(
        strong("Model-based AUC by arm: "), "the area under each arm's fitted growth curve from the first study day to the last day on which every arm is evaluable (at least 50 % of its enrolled animals, and at least 3, on study).",
        tags$ul(
          tags$li(strong("Treatment"), " \u2014 treatment group name."),
          tags$li(strong("AUC"), " \u2014 area under the arm's fitted geometric-mean curve (volume \u00d7 days). Higher values mean more total tumour burden over the window."),
          tags$li(strong("Lower_CL / Upper_CL"), " \u2014 95% interval, from draws of the model."),
          tags$li(strong("N"), " \u2014 animals enrolled in the arm."),
          tags$li(strong("Reference"), " \u2014 TRUE for the reference group.")
        ),
        p("Animals that died or were removed early count through all their measurements: the ",
          "model keeps each one's growth rate in its arm's curve, so early deaths of large tumours ",
          "do not pull an arm's AUC down. See the Info tab for how far this holds."),
        em("For comparisons between arms see the Comparisons tab; per-animal trapezoidal AUCs, for description, are in the AUC Individual Data tab.")
      )
    } else {
      tab_info(
        strong("Estimated marginal means "), "(EMMs, also called least-squares means) from the mixed-effects model, averaged over the time course. These are model-adjusted group means that control for random animal effects and cage structure.",
        tags$ul(
          tags$li(strong("Group"), " \u2014 treatment group name."),
          tags$li(strong("Adj Mean"), " \u2014 model-estimated marginal mean tumor volume for the group. If a transformation was applied, the header will read \u2018Adj Mean (log scale)\u2019 or \u2018Adj Mean (sqrt scale)\u2019."),
          tags$li(strong("SE"), " \u2014 standard error of the estimated marginal mean."),
          tags$li(strong("DF"), " \u2014 Kenward\u2013Roger approximate degrees of freedom for the confidence interval."),
          tags$li(strong("Lower CI / Upper CI"), " \u2014 95% confidence interval bounds on the analysis scale."),
          tags$li(strong("Note"), " \u2014 flags the reference group row.")
        ),
        em("For pairwise statistical comparisons between groups use the Comparisons tab.")
      )
    }
  })

  output$tg_anova_tab_info <- renderUI({
    model_type <- input$tg_model_type %||% "lme4"
    if (model_type == "auc") {
      tab_info(
        strong("Wald test "), "that every arm has the same model-based AUC, on the log scale.",
        tags$ul(
          tags$li(strong("Chisq"), " \u2014 the Wald chi-square statistic."),
          tags$li(strong("Df"), " \u2014 number of arms minus 1."),
          tags$li(strong("p_value"), " \u2014 p-value; below 0.05 indicates that the arms' AUCs differ.")
        ),
        "For pairwise differences see the Comparisons tab."
      )
    } else {
      # Backend 0.30.0 (R20.31): growth rate first, and Treatment at the mean
      # study day. The old Type III "Treatment" row tested day 0, so its
      # p-value depended on where time was counted from.
      tab_info(
        strong("F-tests "), "from the linear mixed-effects model, with Satterthwaite degrees of freedom. Each row answers the question in its ", em("Tests"), " column:",
        tags$ul(
          tags$li(strong("First row, treatment \u00d7 time"), " (e.g. ", em("Treatment:Day"), ") \u2014 do the arms' growth rates differ? This is the primary test: a small p-value means the arms' tumours grow at different rates on the analysis scale."),
          tags$li(strong("Second row, treatment at day X"), " \u2014 do the arms differ in volume at the mean study day, the day the adjusted means use? It does not depend on where time is counted from."),
          tags$li(strong("Third row, time (average over arms)"), " \u2014 is there growth on average? Almost always significant; it confirms the model sees growth."),
          tags$li("Cage or necrosis rows, when present, are adjustment terms.")
        ),
        tags$ul(
          tags$li(strong("F value"), " \u2014 F statistic; larger values are stronger evidence against no difference."),
          tags$li(strong("NumDF / DenDF"), " \u2014 numerator degrees of freedom (number of constraints tested) and Satterthwaite denominator degrees of freedom."),
          tags$li(strong("Pr(>F)"), " \u2014 p-value for the question in that row.")
        ),
        "With per-animal random slopes (the default), the growth-rate test is judged against how much animals' growth rates vary; with a random intercept only it is not, and it finds differences that are not there."
      )
    }
  })

  # R20.D30 / R20.D61: the Welch variance panel was removed; the backend's
  # model-based AUC (0.26.0) returns no variance test, so it never rendered.


  output$tg_perm_test_ui <- renderUI({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    pt <- rv$tg_results$permutation_test
    # R20.D6: the test runs on the background worker after the fit.
    st <- rv$tg_perm_status
    if (is.null(pt) && !is.null(st)) {
      return(tagList(
        tags$hr(),
        h4("Randomisation test"),
        switch(st,
          running = tagList(
            p(class = "text-muted", icon("spinner", class = "fa-spin"),
              " Running on the background worker; the result appears here when ",
              "it finishes. The rest of this analysis, and other users, are not held up."),
            actionButton(session$ns("tg_cancel_perm"), "Cancel the randomisation test",
                         class = "btn-default btn-sm")),
          busy = p(class = "text-muted",
            "Not run: another analysis from this browser was already running. ",
            "Run the analysis again when it has finished."),
          cancelled = p(class = "text-muted", "Cancelled."),
          p(class = "text-danger", st))
      ))
    }
    if (is.null(pt)) {
      return(tagList(
        tags$hr(),
        h4("Randomisation test"),
        p(class = "text-muted",
          "Not run. Enable it in the sidebar to test the treatment \u00d7 time ",
          "interaction without relying on the denominator-df approximation or on ",
          "normally-distributed random effects. You must state the unit that was ",
          "randomised \u2014 animals or cages \u2014 because the test is only valid if it ",
          "mirrors how the study assigned treatment.")
      ))
    }

    floor_note <- if (is.finite(pt$min_attainable_p) && pt$min_attainable_p > 0.05) {
      status_alert(sprintf(paste0(
        "This design cannot produce p < 0.05 at the %s level: %d %ss in these group ",
        "sizes give only %.0f distinct statistic values, so the smallest attainable ",
        "p-value is %.3f. A non-significant result here is a statement about the ",
        "design's resolution, not about the treatment."),
        pt$unit, pt$n_units, pt$unit, pt$n_distinct_statistics, pt$min_attainable_p),
        "warning")
    } else NULL

    tagList(
      tags$hr(),
      h4("Randomisation test of the treatment \u00d7 time interaction"),
      p(class = "text-muted", pt$method),
      status_alert(sprintf(
        "Permutation p = %.4f (%s). Asymptotic likelihood-ratio p = %.3g.",
        pt$p_value,
        if (isTRUE(pt$exhaustive)) "exact \u2014 every assignment enumerated"
        else paste0(pt$n_perm_converged, " permutations"),
        pt$lrt_p_asymptotic),
        if (pt$p_value < 0.05) "success" else "info"),
      if (is.finite(pt$lrt_p_asymptotic) && pt$lrt_p_asymptotic > 0 &&
          pt$p_value / pt$lrt_p_asymptotic > 100) {
        status_alert(paste0(
          "The asymptotic p-value is more than 100\u00d7 smaller than the ",
          "randomisation p-value. At these group sizes the chi-square ",
          "approximation to the mixed-model likelihood-ratio test is unreliable, ",
          "and the randomisation result is the one to report."),
          "warning")
      },
      floor_note,
      p(class = "text-muted",
        strong("Unit randomised: "), pt$unit,
        " \u2014 ", pt$n_units, " exchangeable units.")
    )
  })

  output$tg_anova_table <- DT::renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$model_type_used == "bayes")) {
      return(DT::datatable(data.frame(
        Message = paste(
          "ANOVA is not applicable for Bayesian models.",
          "See the Posterior Check tab for model fit validation",
          "and the Diagnostics tab for convergence assessment."
        )
      ), options = list(paging = FALSE, searching = FALSE)))
    }
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$anova)) {
      return(DT::datatable(data.frame(Message = "ANOVA table not available")))
    }

    anova_df <- as.data.frame(rv$tg_results$anova)
    anova_df <- format_df_for_display(anova_df)

    DT::datatable(anova_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Comparisons title (reactive based on selection)
  # Comparisons table
  output$tg_comparisons_table <- DT::renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) {
      return(DT::datatable(data.frame(Message = "Comparisons not available")))
    }
    
    # R20.D44: the family is the run's, not the dropdown's current value. A
    # dropdown changed after the run showed the other family's rows, or
    # "Comparisons not available for tukey", with no hint to re-run.
    fam <- rv$tg_results$comparison_family %||% "vs_reference"
    all_pairs <- identical(fam, "all_pairs")
    comp_df <- if (all_pairs) {
      rv$tg_results$all_pairwise_comparisons
    } else {
      rv$tg_results$reference_comparisons_dunnett %||% rv$tg_results$pairwise_comparisons
    }
    if (is.null(comp_df)) {
      return(DT::datatable(data.frame(Message = "Comparisons not available for this run")))
    }
    comp_df <- as.data.frame(comp_df)
    table_title <- comparisons_title(rv$tg_results,
                                     if (all_pairs) "All Pairwise" else "vs Reference Group")
    shown <- if (all_pairs) "tukey" else "reference"
    if (!identical(input$tg_comparison_type %||% shown, shown)) {
      table_title <- paste0(table_title, ". The Post-hoc Comparisons setting has ",
                            "changed since this run; run the analysis again to see ",
                            "those comparisons.")
    }

    # Add fold-change back-transform for LME4 log-scale models
    transform_used   <- rv$tg_results$transform_used  %||% "none"
    model_type_used  <- rv$tg_results$model_type_used %||% "lme4"
    if (model_type_used == "lme4" && transform_used == "log") {
      est_col <- intersect(c("estimate", "Estimate"), colnames(comp_df))[1]
      ic_c    <- resolve_interval_cols(comp_df, rv$tg_results)
      lo_col  <- ic_c$lower
      hi_col  <- ic_c$upper
      if (!is.na(est_col)) {
        comp_df$Fold_Change <- round(exp(as.numeric(comp_df[[est_col]])), 3)
        if (!is.na(lo_col)) comp_df$Fold_Change_Lower <- round(exp(as.numeric(comp_df[[lo_col]])), 3)
        if (!is.na(hi_col)) comp_df$Fold_Change_Upper <- round(exp(as.numeric(comp_df[[hi_col]])), 3)
      }
    }

    # Enforce column order. Bayesian path also surfaces Lower_CrI / Upper_CrI
    # (credible intervals) and P_direction (CODE_REVIEW.md Round 2 D —
    # mE v0.4.6 directional posterior probability; replaces the awkward
    # "does the 95% CrI exclude 0?" check with a quantitative posterior).
    preferred_order <- c("contrast", "Contrast",
                         "df",
                         "estimate", "Estimate",
                         "lower.CL", "lower_CL",
                         "upper.CL", "upper_CL",
                         "Lower_CrI", "Upper_CrI",
                         "SE",
                         "t.ratio", "z.ratio",
                         "p.value",
                         "P_direction",
                         "Fold_Change", "Fold_Change_Lower", "Fold_Change_Upper")
    ordered_cols <- intersect(preferred_order, colnames(comp_df))
    rest_cols    <- setdiff(colnames(comp_df), ordered_cols)
    comp_df      <- comp_df[, c(ordered_cols, rest_cols), drop = FALSE]

    comp_df <- format_df_for_display(comp_df)
    # R20.D44: comparisons_title() was computed and never shown.
    DT::datatable(comp_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE),
              caption = tags$caption(
                style = "caption-side: top; text-align: left; font-weight: bold;",
                table_title))
  })
  
  # Treatment effects table
  output$tg_effects_table <- DT::renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$treatment_effects)) {
      return(DT::datatable(data.frame(Message = "Treatment effects not available")))
    }

    effects_df <- as.data.frame(rv$tg_results$treatment_effects)
    transform_used   <- rv$tg_results$transform_used  %||% "none"
    model_type_used  <- rv$tg_results$model_type_used %||% "lme4"

    if (model_type_used == "lme4") {
      # Keep only the columns the user needs; drop any extras / _orig columns
      keep <- c("Group", "DF", "Adjusted_Mean", "Lower_CL", "Upper_CL", "SE", "Note")
      effects_df <- effects_df[, intersect(keep, colnames(effects_df)), drop = FALSE]

      # Enforce display order
      effects_df <- effects_df[, intersect(keep, colnames(effects_df)), drop = FALSE]

      effects_df <- format_df_for_display(effects_df)

      scale_note <- switch(transform_used, "log" = " (log scale)", "sqrt" = " (sqrt scale)", "")
      colnames(effects_df) <- gsub("^Adjusted_Mean$", paste0("Adj Mean", scale_note), colnames(effects_df))
      colnames(effects_df) <- gsub("^Lower_CL$",      paste0("Lower CI", scale_note), colnames(effects_df))
      colnames(effects_df) <- gsub("^Upper_CL$",      paste0("Upper CI", scale_note), colnames(effects_df))
    } else {
      # Non-lme4 models: show all available columns, formatted
      effects_df <- format_df_for_display(effects_df)
    }

    DT::datatable(effects_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Growth rates table
  output$tg_growth_rates_table <- DT::renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$growth_rates)) {
      return(DT::datatable(data.frame(Message = "Growth rates not available")))
    }
    
    rates_df <- as.data.frame(rv$tg_results$growth_rates)
    if ("growth_rate" %in% colnames(rates_df)) {
      # Back-transformed point estimates
      rates_df$daily_fold_change  <- exp(rates_df$growth_rate)
      rates_df$pct_change_per_day <- (exp(rates_df$growth_rate) - 1) * 100

      # Back-transformed 95% CIs (from per-animal lm confint)
      if ("growth_rate_lower" %in% colnames(rates_df)) {
        rates_df$fold_change_lower    <- exp(rates_df$growth_rate_lower)
        rates_df$fold_change_upper    <- exp(rates_df$growth_rate_upper)
        rates_df$pct_change_lower     <- (exp(rates_df$growth_rate_lower) - 1) * 100
        rates_df$pct_change_upper     <- (exp(rates_df$growth_rate_upper) - 1) * 100
      }

      # Enforce column order:
      # growth_rate, growth_rate_lower, growth_rate_upper,
      # daily_fold_change, fold_change_lower, fold_change_upper,
      # pct_change_per_day, pct_change_lower, pct_change_upper,
      # R_squared, then ID/Treatment/Cage
      preferred <- c("growth_rate", "growth_rate_lower", "growth_rate_upper",
                     "daily_fold_change", "fold_change_lower", "fold_change_upper",
                     "pct_change_per_day", "pct_change_lower", "pct_change_upper",
                     "R_squared")
      ordered_cols <- intersect(preferred, colnames(rates_df))
      rest_cols    <- setdiff(colnames(rates_df), ordered_cols)
      rates_df     <- rates_df[, c(rest_cols, ordered_cols), drop = FALSE]
    }
    rates_df <- format_df_for_display(rates_df)
    # Human-readable column headers
    rename_map <- c(
      growth_rate        = "Growth Rate",
      growth_rate_lower  = "GR 95% CI Lower Bound",
      growth_rate_upper  = "GR 95% CI Upper Bound",
      daily_fold_change  = "Fold Change",
      fold_change_lower  = "FC 95% CI Lower Bound",
      fold_change_upper  = "FC 95% CI Upper Bound",
      pct_change_per_day = "Percent Change",
      pct_change_lower   = "PC 95% CI Lower Bound",
      pct_change_upper   = "PC 95% CI Upper Bound",
      R_squared          = "R\u00b2"
    )
    new_names <- colnames(rates_df)
    for (old in names(rename_map)) new_names[new_names == old] <- rename_map[[old]]
    # CODE_REVIEW.md Round 2 F.4 — single-row selection enables the
    # per-mouse drill-down plot rendered just below the table.
    DT::datatable(rates_df,
      colnames = new_names,
      selection = list(mode = "single", target = "row"),
      options  = list(scrollX = TRUE, paging = FALSE, searching = FALSE, scrollY = "500px", scrollCollapse = TRUE)
    )
  })

  # Per-mouse drill-down plot (CODE_REVIEW.md Round 2 F.4)
  # Shows the selected mouse's volume trajectory overlaid against its
  # treatment-group mean trajectory and the fitted log-linear slope.
  output$tg_drilldown_plot <- renderPlot({
    req(rv$tg_results)
    req(rv$tg_plot_data, rv$tg_col_names)
    rates_df <- as.data.frame(rv$tg_results$growth_rates)
    if (is.null(rates_df) || nrow(rates_df) == 0L) {
      plot.new(); title(main = "No per-animal data available"); return()
    }
    sel <- input$tg_growth_rates_table_rows_selected
    if (length(sel) == 0L) {
      plot.new()
      title(main = "Select a row in the table above to drill down",
            sub  = "Each row = one mouse. Click to see that mouse's trajectory.",
            cex.main = 1.0, cex.sub = 0.8)
      return()
    }
    sel <- as.integer(sel[1])
    if (sel < 1L || sel > nrow(rates_df)) return()

    id_col   <- rv$tg_col_names$id_column        %||% "ID"
    tx_col   <- rv$tg_col_names$treatment_column %||% "Treatment"
    time_col <- rv$tg_col_names$time_column      %||% "Day"
    vol_col  <- rv$tg_col_names$volume_column    %||% "Volume"

    cage_col <- rv$tg_col_names$cage_column

    target_id <- as.character(rates_df$ID[sel])
    target_tx <- as.character(rates_df$Treatment[sel])
    raw_data  <- rv$tg_plot_data

    target_cage <- if ("Cage" %in% names(rates_df) && has_col(raw_data, cage_col))
                     as.character(rates_df$Cage[sel])
    mouse_rows <- drilldown_rows(raw_data, rates_df[sel, , drop = FALSE],
                                 id_col, tx_col, cage_col)
    group_rows <- raw_data[as.character(raw_data[[tx_col]]) == target_tx, , drop = FALSE]

    if (nrow(mouse_rows) == 0L) {
      plot.new(); title(main = paste("No raw data found for", target_id)); return()
    }

    group_mean <- stats::aggregate(
      stats::as.formula(paste(vol_col, "~", time_col)),
      data = group_rows, FUN = mean, na.rm = TRUE
    )
    names(group_mean)[2] <- "Mean_Volume"

    gr_rate <- rates_df$growth_rate[sel]

    p <- ggplot2::ggplot() +
      ggplot2::geom_line(data = group_mean,
        ggplot2::aes(x = .data[[time_col]], y = Mean_Volume),
        colour = "grey60", linewidth = 0.7, linetype = "dashed") +
      ggplot2::geom_point(data = mouse_rows,
        ggplot2::aes(x = .data[[time_col]], y = .data[[vol_col]]),
        size = 3, colour = "steelblue") +
      ggplot2::geom_line(data = mouse_rows,
        ggplot2::aes(x = .data[[time_col]], y = .data[[vol_col]]),
        colour = "steelblue", linewidth = 1) +
      ggplot2::labs(
        title = sprintf("Drill-down: %s", animal_label(target_tx, target_id, target_cage)),
        subtitle = sprintf(
          "Per-animal growth rate = %.4f log-vol/day (= %.2fx daily fold change). Grey dashed = group mean.",
          gr_rate, exp(gr_rate)),
        x = time_col, y = vol_col
      ) +
      ggplot2::theme_classic(base_size = 14)
    print(p)
  })

  # CODE_REVIEW.md Round 2 F.1 — standalone HTML report download
  output$download_tg_report <- downloadHandler(
    filename = function() {
      paste0("tumor_growth_report_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
    },
    content = function(file) {
      req(rv$tg_results)
      # R20.D46: the run's settings, not the sidebar's current values.
      build_html_report(
        result   = rv$tg_results,
        title    = "Tumor Growth Analysis Report",
        settings = tg_report_settings(rv$tg_results),
        file     = file
      )
      shinyjs::runjs("window.trackDownload('tumor_growth', 'report', 'html');")
    }
  )

  # AUC individual data table
  output$tg_auc_individual_table <- DT::renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$auc_analysis$individual)) {
      return(DT::datatable(data.frame(Message = "AUC individual data not available")))
    }
    
    auc_ind_df <- as.data.frame(rv$tg_results$auc_analysis$individual)
    auc_ind_df <- format_df_for_display(auc_ind_df)
    DT::datatable(auc_ind_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # AUC summary table
  output$tg_auc_summary_table <- DT::renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$auc_analysis$summary)) {
      return(DT::datatable(data.frame(Message = "AUC summary not available")))
    }
    
    auc_sum_df <- as.data.frame(rv$tg_results$auc_analysis$summary)
    auc_sum_df <- format_df_for_display(auc_sum_df)
    DT::datatable(auc_sum_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Data summary table
  output$tg_data_summary_table <- DT::renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$data_summary)) {
      return(DT::datatable(data.frame(Message = "Data summary not available")))
    }
    
    # R20.D45: the backend's summary is on the volume scale since 0.31.0
    # (R20.57); before, log-scale means were shown under "mm³" (4.31 for a
    # mean of 75.9). The headers carry the declared unit.
    data_sum_df <- data_summary_display(rv$tg_results$data_summary, volume_unit_label(rv))
    data_sum_df <- format_df_for_display(data_sum_df)
    DT::datatable(data_sum_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE),
              rownames = FALSE)
  })
  
  # Plot reactives + plot renderers + diagnostics live in
  # `mod_tumor_growth_plots.R` (CODE_REVIEW.md C.1/C.3). Returned list is
  # captured locally so download handlers can receive the same reactives.
  tg_plots <- add_tg_plots(input, output, session, rv, ordered_tg_grps)
  tg_growth_plot <- tg_plots$growth
  tg_rates_plot  <- tg_plots$rates
  tg_auc_plot    <- tg_plots$auc
  tg_wf_plot     <- tg_plots$waterfall
  tg_forest_plot <- tg_plots$forest

  # --- Bayesian-specific output renderers ---

  output$tg_posterior_summary <- DT::renderDT({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    ps <- rv$tg_results$posterior_summary
    req(!is.null(ps))
    DT::datatable(
      format_df_for_display(as.data.frame(ps)),
      options = list(scrollX = TRUE, paging = FALSE, searching = FALSE)
    )
  })

  output$tg_pp_check_plot <- renderPlot({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    # R20.D18: from the result; the fitted model is no longer kept.
    req(!is.null(rv$tg_results$pp_check_plot))
    rv$tg_results$pp_check_plot
  })

  output$tg_mcmc_diag_table <- DT::renderDT({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    diag <- rv$tg_results$mcmc_diagnostics
    req(!is.null(diag))
    DT::datatable(
      format_df_for_display(as.data.frame(diag)),
      options = list(scrollX = TRUE, paging = FALSE, searching = FALSE)
    )
  })

  # CODE_REVIEW.md Round 2 D — Bayesian diagnostics panel (NUTS banner +
  # Bayes R² + PPC coverage + LOO summary).
  output$tg_bayes_diagnostics <- renderUI({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    bayes_diagnostics_panel(rv$tg_results)
  })

  output$tg_bayes_influential <- DT::renderDT({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    inf_df <- bayes_influential_obs(rv$tg_results)
    if (nrow(inf_df) == 0L) {
      return(DT::datatable(
        data.frame(Note = sprintf("No observations exceed the Pareto-k threshold (%s).",
                                  format(rv$tg_results$loo_diagnostics$k_threshold %||% 0.7))),
        options = list(dom = "t", paging = FALSE, searching = FALSE), rownames = FALSE
      ))
    }
    DT::datatable(
      inf_df,
      options = list(scrollX = TRUE, paging = FALSE, searching = FALSE), rownames = FALSE
    )
  })

  output$tg_mcmc_trace_plot <- renderPlot({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    req(!is.null(rv$tg_results$posterior_draws))
    bayesplot::mcmc_trace(rv$tg_results$posterior_draws)
  })

  # CODE_REVIEW.md DIAGNOSTICS gaps (10), (11), (12) — Bayesian extras.
  output$tg_ppc_coverage_plot <- renderPlot({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    p <- build_ppc_coverage_plot(rv$tg_results)
    req(!is.null(p))
    p
  })

  output$tg_pareto_k_plot <- renderPlot({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    p <- build_pareto_k_plot(rv$tg_results)
    req(!is.null(p))
    p
  })

  output$tg_mcmc_rank_plot <- renderPlot({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!is.null(rv$tg_results$posterior_draws))
    p <- build_mcmc_rank_plot(rv$tg_results$posterior_draws)
    req(!is.null(p))
    p
  })

  # CODE_REVIEW.md DIAGNOSTICS gap (2) — render the previously-wasted
  # backend plots. They're populated only when bayesian_tumor_growth()
  # was called with plots = TRUE; the dashboard's F.7 path passes
  # plots = FALSE, so users see "not available" until we either flip the
  # default or extend the F.7 path to compute these on demand.
  output$tg_posterior_dist_plot <- renderPlot({
    req(rv$tg_results)
    p <- rv$tg_results$posterior_dist_plot
    if (is.null(p)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Posterior densities not available\n(fit with plots = TRUE to generate)")
      return(invisible(NULL))
    }
    p
  })

  output$tg_bayes_residuals_plot <- renderPlot({
    req(rv$tg_results)
    p <- rv$tg_results$residuals_plot
    if (is.null(p)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Bayesian residuals not available\n(fit with plots = TRUE to generate)")
      return(invisible(NULL))
    }
    p
  })

  # Necrosis table render + all download handlers live in
  # `mod_tumor_growth_downloads.R` (CODE_REVIEW.md C.1/C.3).
  add_tg_download_handlers(
    input, output, session, rv,
    tg_growth_plot, tg_rates_plot,
    tg_auc_plot, tg_forest_plot, tg_wf_plot
  )

  # F.7 — register the resume handler. The banner module dispatches to
  # this when the user clicks Resume on a bayes_tg job. Captures the
  # module's input / output / session / rv via closure.
  register_tg_resume(input, output, session, rv)
  # R20.D6: the randomisation test's Cancel button, and Resume from the banner.
  register_bayes_resume(
    kind = "tg_perm", tab_name = "Tumor Growth",
    modal_title = NULL, modal_detail = NULL,
    input = input, output = output, session = session, rv = rv,
    job_id_field = "tg_perm_job_id",
    finalize_success      = .finalize_perm_success,
    finalize_failure      = .finalize_perm_failure,
    finalize_cancellation = .finalize_perm_cancellation,
    finalize_lost         = .finalize_perm_lost,
    cancel_input          = "tg_cancel_perm")

  })
}