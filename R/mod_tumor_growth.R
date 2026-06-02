
#' Tumor Growth Analysis Server Logic
#' @keywords internal
tumor_growth_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
  
  rv$tg_results    <- NULL
  rv$tg_running    <- FALSE
  rv$bayes_models  <- list()

  ordered_tg_grps <- reactive(
    get_ordered_groups(input, "tg", rv$tg_plot_groups)
  )

  # Update reference group choices when data is available or column mapping changes.
  # When a dose column is also mapped, group choices are shown as "Treatment Dose"
  # combinations so each dose level is independently selectable as the reference.
  observe({
    if (!is.null(rv$uploaded_data) && rv$data_valid) {
      # Prefer the user-selected treatment column; fall back to auto-detected one
      treatment_col <- if (!is.null(rv$settings$col_treatment) && nzchar(rv$settings$col_treatment))
                         rv$settings$col_treatment
                       else
                         rv$detected_cols$treatment
      dose_col <- if (!is.null(rv$settings$col_dose) && nzchar(rv$settings$col_dose))
                    rv$settings$col_dose
                  else NULL

      if (!is.null(treatment_col) && nzchar(treatment_col)) {
        df <- rv$uploaded_data

        if (!is.null(dose_col) && dose_col %in% colnames(df)) {
          # Build unique treatment+dose combination labels
          treatment_values <- sort(unique(na.omit(
            paste(as.character(df[[treatment_col]]),
                  as.character(df[[dose_col]]), sep = " ")
          )))
        } else {
          treatment_values <- sort(unique(na.omit(as.character(df[[treatment_col]]))))
        }

        if (length(treatment_values) > 0) {
          new_choices <- c("-- Select a reference group --" = "",
                          setNames(treatment_values, treatment_values))
          default_ref <- pick_reference_group(treatment_values)
          updateSelectInput(session, "tg_reference", choices = new_choices,
                            selected = default_ref)
        }
      }
    }
  })
  
  # Run analysis button handler
  observeEvent(input$run_tg_analysis, {
    req(rv$uploaded_data, rv$data_valid)
    
    # Validate reference group selection
    if (input$tg_reference == "" || is.null(input$tg_reference)) {
      output$tg_status_msg <- renderUI({
        status_alert("Error: Please select a reference group before running analysis", "danger")
      })
      shinyjs::runjs("window.trackError('validation', 'Reference group not selected', 'tumor_growth');")
      return()
    }
    
    # Track analysis start and reference group selection with Umami
    shinyjs::runjs(sprintf(
      "window.trackAnalysisStart('tumor_growth', '%s', '%s'); window.trackReferenceGroupSelected('tumor_growth', '%s');",
      input$tg_model_type, input$tg_reference, input$tg_reference
    ))
    
    rv$tg_running <- TRUE
    output$tg_status_msg <- renderUI({
      status_alert("Running analysis...", "info")
    })
    
    # Prepare data
    # CODE_REVIEW.md Round 2 F.5: apply the global Subset Filter (set on the
    # Data Upload tab) before any column derivation, so excluded mice /
    # treatments / cages are gone from every downstream computation.
    data <- apply_subset_filter(rv$uploaded_data, rv)
    
    # Handle time column - either use direct day or calculate from date
    if (rv$settings$time_method == "from_date") {
      date_info <- detect_date_format(data, rv$settings$col_date_for_day,
                                      rv$settings$reference_date)
      data <- mouseExperiment::calculate_dates(
        date_info$data,
        start_date  = date_info$start_date,
        date_column = rv$settings$col_date_for_day,
        date_format = date_info$date_format
      )
      time_col <- "Day"  # calculate_dates produces a "Day" column
    } else {
      time_col <- rv$settings$col_day
    }
    
    # Apply volume calculation if needed
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
      if (!len_col %in% colnames(data) || !wid_col %in% colnames(data)) {
        missing_col <- if (!len_col %in% colnames(data)) len_col else wid_col
        output$tg_status_msg <- renderUI({
          status_alert(paste0("Column not found in data: '", missing_col, "'"), "danger")
        })
        shinyjs::runjs(sprintf("window.trackError('validation', 'Column not found: %s', 'tumor_growth');",
                               missing_col))
        rv$tg_running <- FALSE
        return()
      }
      ht_col_tg <- if (!is.null(rv$settings$col_height) && nzchar(rv$settings$col_height)) rv$settings$col_height else NULL
      data <- mouseExperiment::calculate_volume(
        data,
        length_column = len_col,
        width_column  = wid_col,
        height_column = ht_col_tg,
        formula       = rv$settings$volume_formula
      )
      vol_col <- "Volume"
    } else {
      vol_col <- rv$settings$col_volume
    }
    
    # Store preprocessed data and resolved column names for reuse by other modules
    rv$processed_data <- data
    rv$time_col       <- time_col
    rv$vol_col        <- vol_col
    rv$tx_col         <- rv$settings$col_treatment
    rv$id_col         <- rv$settings$col_id
    rv$cage_col       <- if (!is.null(rv$settings$col_cage)     && rv$settings$col_cage != "")     rv$settings$col_cage     else NULL
    rv$dose_col       <- if (!is.null(rv$settings$col_dose)     && rv$settings$col_dose != "")     rv$settings$col_dose     else NULL
    rv$necrotic_col   <- if (!is.null(rv$settings$col_necrotic) && rv$settings$col_necrotic != "") rv$settings$col_necrotic else NULL

    # When a dose column is mapped, combine treatment + dose into a single
    # group label so each dose level is treated as its own group by all analyses.
    if (!is.null(rv$dose_col)) {
      data[["TreatmentGroup"]] <- paste(data[[rv$settings$col_treatment]],
                                        data[[rv$dose_col]], sep = " ")
      rv$processed_data <- data
      rv$tx_col         <- "TreatmentGroup"
    }

    # Prepare column names
    col_names <- list(
      time_column = time_col,
      volume_column = vol_col,
      id_column = rv$settings$col_id,
      treatment_column = rv$tx_col,
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
      shinyjs::runjs(sprintf("window.trackError('validation', '%s', 'tumor_growth');",
                             gsub("'", "\\\\'", error_msg)))
      rv$tg_running <- FALSE
      return()
    }
    
    # ── Bayesian async branch (covers both bayes-LMM and bayes-GAM) ─────────
    if (isTRUE(input$tg_model_type == "bayes") ||
        isTRUE(input$tg_model_type == "bayes_gam")) {
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

      # Random effects: Bayesian does not support "none" — fall back to intercept_only
      re_spec <- if (isTRUE(input$tg_random_effects == "none")) "intercept_only" else
                   input$tg_random_effects %||% "intercept_only"

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

      bayes_model_type <- if (isTRUE(input$tg_model_type == "bayes_gam")) "gam" else "lmm"

      bayes_params <- list(
        data              = data,
        col_names         = col_names,
        reference_group   = input$tg_reference,
        model_type        = bayes_model_type,
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
        necrotic_handling = input$tg_necrotic_handling %||% "none",
        plot_groups       = rv$tg_plot_groups
      )

      show_analysis_modal(
        if (bayes_model_type == "gam") "Bayesian GAMM tumor growth fit running…"
        else "Bayesian tumor growth fit running…",
        if (bayes_model_type == "gam")
          "Fitting Bayesian GAMM with brms smoothers. This may take 3–12 minutes."
        else
          "Fitting Bayesian LMM with brms. This may take 3–12 minutes."
      )

      promises::future_promise({
        mouseExperiment::bayesian_tumor_growth(
          df                           = bayes_params$data,
          time_column                  = bayes_params$col_names$time_column,
          volume_column                = bayes_params$col_names$volume_column,
          id_column                    = bayes_params$col_names$id_column,
          treatment_column             = bayes_params$col_names$treatment_column,
          cage_column                  = bayes_params$col_names$cage_column,
          dose_column                  = bayes_params$dose_col,
          model_type                   = bayes_params$model_type,
          random_effects_specification = bayes_params$random_effects,
          reference_group              = bayes_params$reference_group,
          prior_strength               = bayes_params$prior_strength,
          prior_b                      = bayes_params$manual_priors$prior_b,
          prior_intercept              = bayes_params$manual_priors$prior_intercept,
          prior_sd                     = bayes_params$manual_priors$prior_sd,
          prior_sigma                  = bayes_params$manual_priors$prior_sigma,
          n_chains                     = bayes_params$n_chains,
          n_warmup                     = bayes_params$n_warmup,
          n_iter                       = bayes_params$n_iter,
          backend                      = bayes_params$backend,
          necrotic_column              = bayes_params$necrotic_col,
          necrotic_handling            = bayes_params$necrotic_handling,
          return_model                 = TRUE,
          plots                        = FALSE,
          verbose                      = FALSE
        )
      }) %...>% (function(result) {
        shiny::removeModal()
        result$model_type_used <- "bayes"
        if (!is.null(rv$bayes_models)) rv$bayes_models$tg <- result$model
        else rv$bayes_models <- list(tg = result$model)
        result$model <- NULL
        gc()
        # CODE_REVIEW.md Round 2 F.3 — preserve previous result for compare
        # tab when the model_type just changed (e.g. user ran LME4 then bayes).
        if (!is.null(rv$tg_results) &&
            !identical(rv$tg_results$model_type_used %||% "", result$model_type_used)) {
          rv$tg_results_prev <- isolate(rv$tg_results)
        }
        rv$tg_results  <- result
        updateTabsetPanel(session, "tg_results_tabs", selected = "Summary")
        rv$tg_running  <- FALSE
        rv$tg_plot_data <- bayes_params$data
        rv$tg_col_names <- bayes_params$col_names

        output$tg_status_msg <- renderUI(
          status_alert("Bayesian analysis completed successfully!", "success")
        )
        shinyjs::runjs("window.trackAnalysisComplete('tumor_growth', true);")
      }) %...!% (function(err) {
        shiny::removeModal()
        rv$tg_running <- FALSE
        err_msg <- conditionMessage(err)
        output$tg_status_msg <- renderUI(
          status_alert(paste("Bayesian analysis failed:", err_msg), "danger")
        )
        shinyjs::runjs(sprintf(
          "window.trackAnalysisComplete('tumor_growth', false, %s);",
          js_string(err_msg)
        ))
      })

      return(NULL)   # return immediately — Shiny knows this is non-blocking
    }

    show_analysis_modal("Fitting tumor growth model…", "Results will appear automatically when complete.")
    on.exit(shiny::removeModal(), add = TRUE)

    result <- safe_analysis(module = "tumor_growth", {
      res <- mouseExperiment::tumor_growth_statistics(
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
        reference_group = input$tg_reference,  # Already validated to not be empty
        extrapolation_points = if (!is.null(input$tg_extrapolation_points) &&
                                    !is.na(input$tg_extrapolation_points))
                                 as.integer(input$tg_extrapolation_points) else 0L,
        necrotic_column   = rv$necrotic_col,
        necrotic_handling = input$tg_necrotic_handling %||% "none",
        return_model = TRUE,
        include_diagnostics = TRUE,
        plots = FALSE
      )
      res
    })
    
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
        
        # Plot 2: Growth rate distribution
        if (!is.null(result$growth_rates) && nrow(result$growth_rates) > 0) {
          result$growth_rate_plot <- mouseExperiment::plot_growth_rate(
            growth_data = result$growth_rates,
            show_mean = TRUE,
            error_bar_type = "SEM"
          ) + ggplot2::theme_classic()
        }
        
        # Plot 3: AUC / treatment effects plot
        # When model_type="auc" the result contains auc_analysis$individual with
        # per-subject AUC values — use the package's plot_auc() to match the demo.
        # For lme4 models fall back to a bar chart of adjusted means.
        if (!is.null(result$auc_analysis) && !is.null(result$auc_analysis$individual)) {
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
            preferred <- c("Adjusted_Mean", "Mean_AUC", "emmean", "estimate", "mean")
            effect_col <- intersect(preferred, colnames(effects_df))[1]
            if (is.na(effect_col)) {
              # fallback: first numeric that isn't N / DF / Reference
              skip_cols <- c("N", "DF", "Reference")
              num_cols <- colnames(effects_df)[sapply(effects_df, is.numeric)]
              num_cols <- setdiff(num_cols, c(skip_cols, label_col))
              effect_col <- if (length(num_cols) > 0) num_cols[1] else NA
            }
            if (!is.na(label_col) && !is.na(effect_col)) {
              # Detect CI or SE columns for error bars
              lower_col <- intersect(c("lower.CL", "lower_CL", "Lower", "lower"), colnames(effects_df))[1]
              upper_col <- intersect(c("upper.CL", "upper_CL", "Upper", "upper"), colnames(effects_df))[1]
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

              if (!is.null(gr_df) && !is.na(grp_col_gr) && !is.na(val_col_gr)) {
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
                  ggplot2::labs(title = "Treatment Effects (Adjusted Means)",
                               x = "Treatment", y = effect_col) +
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
        if (!is.null(result$diagnostics) && !is.null(result$diagnostics$residuals)) {
          resid_vec <- result$diagnostics$residuals$residuals
          # QQ plot
          qq_data <- qqnorm(resid_vec, plot.it = FALSE)
          qq_df <- data.frame(theoretical = qq_data$x, sample = qq_data$y)
          qq_p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = sample)) +
            ggplot2::geom_point(alpha = 0.6) +
            ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
            ggplot2::theme_classic() +
            ggplot2::labs(title = "Q-Q Plot of Residuals",
                         x = "Theoretical Quantiles",
                         y = paste0("Sample Quantiles", diag_scale_tag))
          result$diag_qq_plot     <- qq_p
        } else if (!is.null(result$model)) {
          tryCatch({
            resid_vec <- stats::residuals(result$model)
            qq_data   <- qqnorm(resid_vec, plot.it = FALSE)
            qq_df     <- data.frame(theoretical = qq_data$x, sample = qq_data$y)
            qq_p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = sample)) +
              ggplot2::geom_point(alpha = 0.6) +
              ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
              ggplot2::theme_classic() +
              ggplot2::labs(title = "Q-Q Plot of Residuals",
                           x = "Theoretical Quantiles",
                           y = paste0("Sample Quantiles", diag_scale_tag))
            result$diag_qq_plot     <- qq_p
          }, error = function(e) {
            log_caught_error("tumor_growth", "QQ plot from model", e)
            NULL
          })
        }
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
          # --- AUC mode: use pre-computed Welch t-tests (avoids pooled-error bias) ---
          pw <- result$posthoc$pairwise
          # Reformat to emmeans-compatible column names for downstream display/plot code
          pw_df <- data.frame(
            contrast = pw$comparison,
            estimate = pw$estimate,
            SE       = abs(pw$estimate / pw$t_value),
            df       = pw$df,
            lower.CL = pw$ci_lower,
            upper.CL = pw$ci_upper,
            t.ratio  = pw$t_value,
            p.value  = pw$p_adjusted,   # Bonferroni-adjusted, all-pairs
            stringsAsFactors = FALSE
          )
          result$all_pairwise_comparisons <- pw_df

          ref_group <- if (!is.null(input$tg_reference) && nzchar(input$tg_reference))
            input$tg_reference else {
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
                  if (!is.na(ref_comp_df$t.ratio[i])) {
                    ref_comp_df$t.ratio[i] <- -ref_comp_df$t.ratio[i]
                  }
                }
              }
            }
            result$reference_comparisons_dunnett <- ref_comp_df
          } else {
            result$reference_comparisons_dunnett <- pw_df
          }

        } else {
          # --- LME4 / other modes: derive comparisons via emmeans ---
          emm <- emmeans::emmeans(result$model, as.formula(paste("~", treatment_col_name)))
          emm_levels <- levels(emm)[[1]]   # character vector of factor levels

          # Type 1: All pairwise comparisons with Tukey adjustment
          # infer=c(TRUE,TRUE) returns both CIs and test statistics (t.ratio, p.value)
          result$all_pairwise_comparisons <- as.data.frame(
            summary(pairs(emm, adjust = "tukey"), infer = c(TRUE, TRUE))
          )

          # Type 2: vs Reference Group with Dunnett correction
          ref_group <- if (!is.null(input$tg_reference) && nzchar(input$tg_reference))
            input$tg_reference else emm_levels[1]
          ref_idx <- which(emm_levels == ref_group)
          if (length(ref_idx) == 1) {
            result$reference_comparisons_dunnett <- as.data.frame(
              summary(emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_idx, adjust = "dunnett"),
                      infer = c(TRUE, TRUE))
            )
          }
        }

        # Plot 5: Effect size forest plot — built here so reference_comparisons_dunnett exists
        pc <- if (!is.null(result$reference_comparisons_dunnett)) result$reference_comparisons_dunnett
              else if (!is.null(result$all_pairwise_comparisons)) result$all_pairwise_comparisons
              else NULL
        if (!is.null(pc)) {
          tryCatch({
            pc_df <- as.data.frame(pc)
            label_col   <- intersect(c("contrast", "Contrast", "comparison", "pair"), colnames(pc_df))[1]
            est_col     <- intersect(c("estimate", "Estimate"), colnames(pc_df))[1]
            se_col      <- intersect(c("SE"), colnames(pc_df))[1]
            lower_col   <- intersect(c("lower.CL", "Lower", "lower_CL", "lower"), colnames(pc_df))[1]
            upper_col   <- intersect(c("upper.CL", "Upper", "upper_CL", "upper"), colnames(pc_df))[1]
            p_col       <- intersect(c("p.value", "p_value", "p.adj", "adj.p.value", "P_Value"), colnames(pc_df))[1]
            if (!is.na(label_col) && !is.na(est_col)) {
              if (is.na(lower_col) && !is.na(se_col)) {
                df_col <- intersect(c("df", "DF"), colnames(pc_df))[1]
                t_mult <- if (!is.na(df_col)) qt(0.975, pmax(as.numeric(pc_df[[df_col]]), 1)) else qnorm(0.975)
                pc_df$lower.CL <- pc_df[[est_col]] - t_mult * pc_df[[se_col]]
                pc_df$upper.CL <- pc_df[[est_col]] + t_mult * pc_df[[se_col]]
                lower_col <- "lower.CL"; upper_col <- "upper.CL"
              }
              p_label <- if (!is.na(p_col)) vapply(pc_df[[p_col]], format_pval, character(1)) else NULL
              using_dunnett <- !is.null(result$reference_comparisons_dunnett) &&
                               identical(pc, result$reference_comparisons_dunnett)
              ref_grp_ef <- if (!is.null(input$tg_reference) && nzchar(input$tg_reference))
                              input$tg_reference else NULL
              contrast_labels <- pc_df[[label_col]]
              if (isTRUE(using_dunnett) && !is.null(ref_grp_ef)) {
                suffix_pat <- paste0("\\s*[-\u2013]\\s*", gsub("([+()])", "\\\\\\1", ref_grp_ef), "$")
                clean_labels <- sub(suffix_pat, "", contrast_labels)
                if (all(clean_labels != contrast_labels)) contrast_labels <- clean_labels
              }
              forest_df <- data.frame(
                Contrast = contrast_labels,
                Estimate = pc_df[[est_col]],
                Lower    = if (!is.na(lower_col)) pc_df[[lower_col]] else pc_df[[est_col]],
                Upper    = if (!is.na(upper_col)) pc_df[[upper_col]] else pc_df[[est_col]],
                P_Label  = if (is.null(p_label)) "" else paste0("p=", p_label),
                stringsAsFactors = FALSE
              )
              if (!is.null(ref_grp_ef)) {
                forest_df <- rbind(forest_df, data.frame(
                  Contrast = paste0(ref_grp_ef, "\n(Reference)"),
                  Estimate = 0, Lower = NA_real_, Upper = NA_real_,
                  P_Label  = "", stringsAsFactors = FALSE
                ))
              }
              forest_df$is_ref <- seq_len(nrow(forest_df)) == nrow(forest_df) & !is.null(ref_grp_ef)
              all_vals <- c(forest_df$Estimate[!forest_df$is_ref],
                            forest_df$Lower[!forest_df$is_ref],
                            forest_df$Upper[!forest_df$is_ref])
              max_abs  <- ceiling(max(abs(all_vals), na.rm = TRUE))
              x_lim    <- c(-max_abs, max_abs * 1.20)
              x_breaks <- pretty(c(-max_abs, max_abs), n = 6)
              x_breaks <- x_breaks[x_breaks >= -max_abs & x_breaks <= max_abs]
              result$effect_forest_plot <- ggplot2::ggplot(
                forest_df,
                ggplot2::aes(x = Estimate, y = stats::reorder(Contrast, -Estimate))
              ) +
                ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.6) +
                ggplot2::geom_errorbarh(
                  data = ~ subset(., !is_ref),
                  ggplot2::aes(xmin = Lower, xmax = Upper),
                  height = 0.25, colour = "steelblue"
                ) +
                ggplot2::geom_point(data = ~ subset(., !is_ref), size = 3, colour = "steelblue") +
                ggplot2::geom_point(data = ~ subset(., is_ref),  size = 4, colour = "grey50", shape = 18) +
                ggplot2::geom_text(
                  data = ~ subset(., !is_ref & nzchar(P_Label) & !is.na(Upper)),
                  ggplot2::aes(label = P_Label, x = Upper),
                  hjust = -0.2, size = 3.2
                ) +
                ggplot2::scale_x_continuous(limits = x_lim, breaks = x_breaks) +
                ggplot2::labs(
                  title = "Effect Size Forest Plot",
                  x = if (isTRUE(using_dunnett)) "Mean Difference (treatment \u2212 reference)" else "Mean Difference",
                  y = NULL
                ) +
                ggplot2::theme_classic() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
            }
          }, error = function(e) {
            log_caught_error("tumor_growth", "forest plot generation", e)
          })
        }

      }, error = function(e) {
        log_caught_error("tumor_growth", "comparison generation", e)
      })
    }

    # Persist the settings used so renderDTs can condition on them correctly
    # (the backend result$methods list is not reliably populated for all paths)
    result$transform_used  <- input$tg_transform   %||% "none"
    result$model_type_used <- input$tg_model_type  %||% "lme4"

    # CODE_REVIEW.md Round 2 F.3 — preserve previous result for compare tab
    # when the model_type just changed.
    if (!is.null(rv$tg_results) &&
        !identical(rv$tg_results$model_type_used %||% "", result$model_type_used)) {
      rv$tg_results_prev <- isolate(rv$tg_results)
    }
    rv$tg_results <- result
    updateTabsetPanel(session, "tg_results_tabs", selected = "Summary")
    rv$tg_running <- FALSE

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
      output$tg_status_msg <- renderUI({
        status_alert("Analysis completed successfully!", "success")
      })
    }
  })
  
  # Dynamic labels for manual scale prior parameter inputs
  output$tg_prior_sd_param_label <- renderUI({
    if (isTRUE(input$tg_prior_sd_dist == "half_normal")) "Scale σ" else "Rate"
  })
  output$tg_prior_sigma_param_label <- renderUI({
    if (isTRUE(input$tg_prior_sigma_dist == "half_normal")) "Scale σ" else "Rate"
  })

  # Necrotic handling selector — visible only when a necrotic column is mapped
  output$necrotic_handling_ui <- renderUI({
    if (is.null(rv$necrotic_col)) return(NULL)
    tagList(
      hr(),
      h5(strong("Necrotic Observations")),
      selectInput(session$ns("tg_necrotic_handling"), NULL,
        choices = c(
          "Exclude from analysis"    = "exclude",
          "Include as covariate"     = "covariate",
          "No adjustment (raw data)" = "none"
        ),
        selected = "exclude"
      ),
      helpText(
        strong("Exclude"), " — removes flagged observations before fitting the model.",
        tags$br(),
        strong("Covariate"), " — adds a 0/1 Necrotic term to the fixed effects.",
        tags$br(),
        strong("No adjustment"), " — ignores the flag (keeps all observations)."
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
        paste("  Smoother basis k:", summary_list$model_specification$smoother_k %||% "(linear time — LME4)")
      )
    }

    summary_lines <- c(summary_lines,
      "",
      "NOTES:",
      paste("  ", paste(summary_list$notes, collapse="\n  "))
    )
    
    paste(summary_lines, collapse = "\n")
  })

  # Update plot selector choices based on model type (AUC plot only in AUC mode;

  # forest plot only in LME4 mode)
  observeEvent(input$tg_model_type, {
    req(input$tg_model_type)
    shinyjs::runjs(sprintf("window.trackEvent('model_type_selected', {analysis_type: 'tumor_growth', method: '%s'});",
                           input$tg_model_type))
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

  # Dynamically insert / remove AUC tabs based on model type
  auc_tabs_inserted <- reactiveVal(FALSE)

  observeEvent(input$tg_model_type, {
    tabs_id <- "tg_results_tabs"
    if (isTRUE(input$tg_model_type == "auc") && !isTRUE(auc_tabs_inserted())) {
      insertTab(inputId = tabs_id, target = "Growth Rates", position = "after",
        tabPanel("AUC Individual Data",
          tab_info(
            strong("Per-animal area under the growth curve (AUC) "), "computed for each subject across their observed measurement timepoints. Only displayed in AUC analysis mode.",
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
            DTOutput(session$ns("tg_auc_individual_table"))
          ),
          downloadButton(session$ns("download_tg_auc_individual"), "Download AUC Individual Data")
        )
      )
      insertTab(inputId = tabs_id, target = "AUC Individual Data", position = "after",
        tabPanel("AUC Summary",
          tab_info(
            strong("Group-level AUC summary statistics"), " showing mean area under the growth curve, standard deviation, sample size, and standard error for each treatment group. Only displayed in AUC analysis mode.",
            tags$ul(
              tags$li(strong("Treatment"), " \u2014 treatment group name."),
              tags$li(strong("Mean_AUC"), " \u2014 group mean area under curve (mm\u00b3 \u00d7 days)."),
              tags$li(strong("SD_AUC"), " \u2014 standard deviation of AUC within the group."),
              tags$li(strong("N"), " \u2014 number of animals in the group with valid AUC data."),
              tags$li(strong("SEM_AUC"), " \u2014 standard error of the mean AUC (SD / \u221aN).")
            ),
            em("Use these summaries for descriptive group comparisons and to verify sample sizes and variability before statistical testing.")
          ),
          shinycssloaders::withSpinner(
            DTOutput(session$ns("tg_auc_summary_table"))
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
            strong("Posterior predictive check (PPC)"), " — the model generates new datasets from the posterior and overlays their densities on the observed data. Good model fit means the simulated densities closely envelope the observed density.",
            tags$ul(
              tags$li(strong("y"), " — observed (transformed) tumor volume values."),
              tags$li(strong("y_rep"), " — 50 simulated datasets drawn from the posterior predictive distribution.")
            ),
            em("If the observed density falls well outside the simulated envelope, consider a different transformation or random-effects specification.")
          ),
          shinycssloaders::withSpinner(
            plotOutput(session$ns("tg_pp_check_plot"), height = "450px")
          ),
          downloadButton(session$ns("download_tg_pp_check"), "Download Plot")
        ),
        select = FALSE
      )
      appendTab(inputId = tabs_id,
        tabPanel("MCMC Diagnostics",
          tab_info(
            strong("MCMC convergence diagnostics"), " for all fixed-effect parameters. All chains should have converged (R̂ ≈ 1.00) and sufficient effective samples before interpreting results.",
            tags$ul(
              tags$li(strong("Parameter"), " — model parameter name."),
              tags$li(strong("Rhat"), " — potential scale-reduction factor; values > 1.01 indicate non-convergence."),
              tags$li(strong("Bulk_ESS"), " — bulk effective sample size; aim for > 400 per chain."),
              tags$li(strong("Tail_ESS"), " — tail effective sample size; important for credible interval accuracy."),
              tags$li(strong("Converged"), " — TRUE when R̂ ≤ 1.01.")
            ),
            em("If Rhat > 1.01, try increasing iterations or chains before reporting results.")
          ),
          # NUTS sampler banner + Bayes R² + PPC coverage + PSIS-LOO summary
          # (CODE_REVIEW.md Round 2 D-series, surfaces mouseExperiment v0.4.6
          # fields: nuts_diagnostics, bayes_R2, ppc_coverage, loo_diagnostics).
          uiOutput(session$ns("tg_bayes_diagnostics")),
          shinycssloaders::withSpinner(
            DTOutput(session$ns("tg_mcmc_diag_table"))
          ),
          br(),
          # Influential observations (Pareto-k > 0.7)
          h5(strong("Influential observations (Pareto-k > 0.7)")),
          em(class = "text-muted",
             "Observations the LOO approximation can't reliably estimate. Investigate these mice individually."),
          shinycssloaders::withSpinner(
            DTOutput(session$ns("tg_bayes_influential"))
          ),
          br(),
          shinycssloaders::withSpinner(
            plotOutput(session$ns("tg_mcmc_trace_plot"), height = "500px")
          )
        ),
        select = FALSE
      )
      bayes_tabs_inserted(TRUE)
    } else if (!isTRUE(input$tg_model_type == "bayes") && isTRUE(bayes_tabs_inserted())) {
      removeTab(inputId = tabs_id, target = "Posterior Check")
      removeTab(inputId = tabs_id, target = "MCMC Diagnostics")
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
            strong("Per-group necrosis summary"), " — counts and percentages of necrotic
            observations, number of affected mice, and the first day necrosis was observed.
            The summary always reflects the full uploaded dataset, regardless of which
            handling mode was selected.",
            tags$ul(
              tags$li(strong("N_Observations"), " — total measurement rows for this group."),
              tags$li(strong("N_Necrotic"), " — number flagged as necrotic."),
              tags$li(strong("Pct_Necrotic"), " — percentage of observations that were necrotic."),
              tags$li(strong("N_Mice_Total"), " — total animals in this group."),
              tags$li(strong("N_Mice_With_Necrosis"), " — animals with at least one necrotic observation."),
              tags$li(strong("First_Necrotic_Day"), " — earliest day a necrotic observation was recorded.")
            )
          ),
          shinycssloaders::withSpinner(
            DTOutput(session$ns("tg_necrosis_table"))
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
    shinyjs::runjs(sprintf("window.trackEvent('plot_type_selected', {analysis_type: 'tumor_growth', plot_type: '%s'});",
                           input$tg_plot_selector))
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
        strong("AUC analysis summary: "), "Reports the analysis type (AUC), number of subjects and treatment groups, reference group, and the per-animal AUC computation method (trapezoidal integration of the volume\u2013time curve)."
      )
    } else {
      tab_info(
        strong("Model summary: "), "Reports the fitted model — model type (LME4 / GAMM / AUC / Bayesian), data transformation applied, random effects structure (subject and cage), number of subjects and treatment groups, reference group, smoother basis dimension when GAMM is selected, and any convergence warnings or singular-fit notes."
      )
    }
  })

  output$tg_comparisons_tab_info <- renderUI({
    model_type <- input$tg_model_type %||% "lme4"
    if (model_type == "auc") {
      tab_info(
        strong("Pairwise Welch\u2019s t-tests "), "on per-animal AUC values between all treatment group pairs, with Bonferroni correction for multiple comparisons.",
        tags$ul(
          tags$li(strong("contrast"), " \u2014 the two groups being compared (\u2018Group A \u2212 Group B\u2019)."),
          tags$li(strong("estimate"), " \u2014 mean difference in AUC between the two groups (mm\u00b3 \u00d7 days). Negative values indicate the first group has lower total tumour burden."),
          tags$li(strong("SE"), " \u2014 standard error of the estimated AUC difference."),
          tags$li(strong("df"), " \u2014 Welch\u2013Satterthwaite approximate degrees of freedom."),
          tags$li(strong("lower.CL / upper.CL"), " \u2014 95% confidence interval bounds for the mean AUC difference."),
          tags$li(strong("t.ratio"), " \u2014 Welch t-statistic. Large absolute values indicate a larger difference relative to within-group variability."),
          tags$li(strong("p.value"), " \u2014 Bonferroni-corrected p-value across all pairwise comparisons. Values below 0.05 indicate a significant difference after correction for multiple testing.")
        )
      )
    } else {
      tab_info(
        strong("Post-hoc pairwise contrasts "), "between treatment groups, computed from estimated marginal means (emmeans). Use the Analysis Parameters dropdown to switch between comparison types.",
        tags$p(
          strong("What is being compared:"), " these contrasts compare model-predicted ",
          strong("tumor volume"), " between groups \u2014 not growth rate. Each group\u2019s marginal mean is evaluated at the ",
          em("average study day"), ", so ", strong("estimate"), " answers: \u201cAt the typical timepoint, how much larger or smaller is one group\u2019s predicted volume relative to another?\u201d",
          " To test whether ", em("growth trajectories diverge"), " across groups, see the ", strong("Treatment:Time"), " interaction row in the ANOVA tab. Per-animal growth rates are in the Growth Rates tab."
        ),
        tags$ul(
          tags$li(strong("contrast"), " \u2014 the two groups being compared (\u2018Treatment \u2212 Reference\u2019 for Dunnett, \u2018A \u2212 B\u2019 for Tukey)."),
          tags$li(strong("df"), " \u2014 Satterthwaite or Kenward\u2013Roger approximate degrees of freedom used for the t-test."),
          tags$li(strong("estimate"), " \u2014 mean difference in tumor volume between the two groups on the analysis scale (log or sqrt if a transformation was applied). Negative values mean the first group has lower volume."),
          tags$li(strong("lower.CL / upper.CL"), " (if shown) \u2014 95% confidence interval bounds for the mean difference on the analysis scale."),
          tags$li(strong("SE"), " \u2014 standard error of the estimated difference."),
          tags$li(strong("Fold_Change"), " (if shown; log-transform models only) \u2014 exp(estimate): the multiplicative ratio of group geometric means on the original volume scale. Values below 1.0 indicate the treatment group has lower volume than reference. Fold_Change_Lower / Fold_Change_Upper give the back-transformed CI bounds."),
          tags$li(strong("t.ratio"), " \u2014 t-statistic (estimate \u00f7 SE). Large absolute values indicate a larger difference relative to its uncertainty."),
          tags$li(strong("p.value"), " \u2014 adjusted p-value (Dunnett correction for \u2018vs Reference\u2019; Tukey HSD for \u2018All Comparisons\u2019). Values below 0.05 indicate a significant difference."),
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
        strong("Group-level AUC means "), "showing the raw mean area under the growth curve for each treatment group, its variability, and sample size.",
        tags$ul(
          tags$li(strong("Treatment"), " \u2014 treatment group name."),
          tags$li(strong("Mean_AUC"), " \u2014 group mean area under the volume\u2013time curve (mm\u00b3 \u00d7 days). Higher values indicate greater total tumour burden over the study."),
          tags$li(strong("SD"), " \u2014 standard deviation of per-animal AUC within the group."),
          tags$li(strong("N"), " \u2014 number of animals with valid AUC data."),
          tags$li(strong("Reference"), " \u2014 TRUE for the reference group, FALSE for all others.")
        ),
        em("These are descriptive group summaries. For statistical comparisons between groups see the Comparisons tab; for individual animal values see the AUC Individual Data tab.")
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
        strong("One-way ANOVA "), "on per-animal AUC (area under the growth curve) values. Tests whether mean AUC differs significantly across treatment groups.",
        tags$ul(
          tags$li(strong("Df"), " — degrees of freedom (Treatment row = number of groups minus 1; Residuals row = total animals minus number of groups)."),
          tags$li(strong("Sum Sq"), " — sum of squares: variance attributable to that source."),
          tags$li(strong("Mean Sq"), " — Sum Sq ÷ Df."),
          tags$li(strong("F value"), " — ratio of the Treatment Mean Sq to the Residuals Mean Sq. Larger values indicate more separation between group means relative to within-group variability."),
          tags$li(strong("Pr(>F)"), " — p-value. Values below 0.05 indicate a statistically significant difference in mean AUC across treatment groups.")
        ),
        "This test assumes independent observations and equal variances across groups. For pairwise differences see the Comparisons tab."
      )
    } else {
      tab_info(
        strong("Type III ANOVA "), "from the linear mixed-effects model. Tests whether each fixed effect (treatment group, time, and their interaction) significantly explains variance in tumor volume after accounting for all other terms.",
        tags$ul(
          tags$li(strong("Chisq"), " — likelihood-ratio chi-square statistic for the effect. Larger values indicate stronger evidence against the null hypothesis."),
          tags$li(strong("Df"), " — degrees of freedom for the test (number of parameters being tested)."),
          tags$li(strong("Pr(>Chisq)"), " — p-value. Values below 0.05 indicate a statistically significant effect at the 5% level.")
        ),
        "The ", em("Treatment:Time"), " interaction row is the most biologically informative: a significant interaction means treatment groups diverge in their growth trajectories over time, not just in overall volume."
      )
    }
  })

  output$tg_anova_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$model_type_used == "bayes")) {
      return(datatable(data.frame(
        Message = paste(
          "ANOVA is not applicable for Bayesian models.",
          "See the Posterior Check tab for model fit validation",
          "and the MCMC Diagnostics tab for convergence assessment."
        )
      ), options = list(paging = FALSE, searching = FALSE)))
    }
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$anova)) {
      return(datatable(data.frame(Message = "ANOVA table not available")))
    }

    anova_df <- as.data.frame(rv$tg_results$anova)
    anova_df <- format_df_for_display(anova_df)

    datatable(anova_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Comparisons title (reactive based on selection)
  output$tg_comparisons_title <- renderUI({
    comp_type <- input$tg_comparison_type
    
    title_text <- switch(comp_type,
      "reference" = h4("Pairwise Comparisons vs Reference Group"),
      "pairwise" = h4("All Pairwise Comparisons (Unadjusted)"),
      "tukey" = h4("All Pairwise Comparisons (Tukey HSD Adjusted)"),
      h4("Pairwise Comparisons")
    )
    
    tagList(
      title_text,
      p(class = "text-muted", style = "margin-bottom: 15px;",
        em("Change comparison type in the 'Post-hoc Comparisons' dropdown in Analysis Parameters"))
    )
  })
  
  # Comparisons table
  output$tg_comparisons_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) {
      return(datatable(data.frame(Message = "Comparisons not available")))
    }
    
    # Select comparison type based on user input
    comp_type <- input$tg_comparison_type
    
    comp_df <- NULL
    table_title <- NULL
    
    if (comp_type == "reference") {
      # vs Reference with Dunnett correction (preferred), fallback to unadjusted
      if (!is.null(rv$tg_results$reference_comparisons_dunnett)) {
        comp_df <- as.data.frame(rv$tg_results$reference_comparisons_dunnett)
        table_title <- "Comparisons: vs Reference Group (Dunnett-Adjusted)"
      } else if (!is.null(rv$tg_results$pairwise_comparisons)) {
        comp_df <- as.data.frame(rv$tg_results$pairwise_comparisons)
        table_title <- "Comparisons: vs Reference Group"
      }
    } else if (comp_type == "tukey") {
      # All pairwise with Tukey HSD adjustment
      if (!is.null(rv$tg_results$all_pairwise_comparisons)) {
        comp_df <- as.data.frame(rv$tg_results$all_pairwise_comparisons)
        table_title <- "Comparisons: All Pairwise (Tukey HSD Adjusted)"
      }
    }
    
    if (is.null(comp_df)) {
      return(datatable(data.frame(Message = paste("Comparisons not available for", comp_type))))
    }

    # Add fold-change back-transform for LME4 log-scale models
    transform_used   <- rv$tg_results$transform_used  %||% "none"
    model_type_used  <- rv$tg_results$model_type_used %||% "lme4"
    if (model_type_used == "lme4" && transform_used == "log") {
      est_col <- intersect(c("estimate", "Estimate"), colnames(comp_df))[1]
      lo_col  <- intersect(c("lower.CL", "lower_CL", "Lower"), colnames(comp_df))[1]
      hi_col  <- intersect(c("upper.CL", "upper_CL", "Upper"), colnames(comp_df))[1]
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
    datatable(comp_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Treatment effects table
  output$tg_effects_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$treatment_effects)) {
      return(datatable(data.frame(Message = "Treatment effects not available")))
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

    datatable(effects_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Growth rates table
  output$tg_growth_rates_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$growth_rates)) {
      return(datatable(data.frame(Message = "Growth rates not available")))
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
      R_squared          = "R²"
    )
    new_names <- colnames(rates_df)
    for (old in names(rename_map)) new_names[new_names == old] <- rename_map[[old]]
    # CODE_REVIEW.md Round 2 F.4 — single-row selection enables the
    # per-mouse drill-down plot rendered just below the table.
    datatable(rates_df,
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

    target_id <- as.character(rates_df$ID[sel])
    target_tx <- as.character(rates_df$Treatment[sel])
    raw_data  <- rv$tg_plot_data

    mouse_rows <- raw_data[as.character(raw_data[[id_col]]) == target_id &
                            as.character(raw_data[[tx_col]]) == target_tx, , drop = FALSE]
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
        title = sprintf("Drill-down: %s (%s)", target_id, target_tx),
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
      settings <- list(
        model_type      = input$tg_model_type     %||% "",
        transform       = input$tg_transform      %||% "",
        random_effects  = input$tg_random_effects %||% "",
        reference_group = input$tg_reference      %||% "",
        cage_effect     = input$tg_cage_effect    %||% "",
        n_rows          = if (!is.null(rv$tg_plot_data)) nrow(rv$tg_plot_data) else NA
      )
      build_html_report(
        result   = rv$tg_results,
        title    = "Tumor Growth Analysis Report",
        settings = settings,
        file     = file
      )
      shinyjs::runjs("window.trackDownload('tumor_growth', 'report', 'html');")
    }
  )

  # CODE_REVIEW.md Round 2 F.3 — side-by-side compare tab
  output$tg_compare_panel <- renderUI({
    if (is.null(rv$tg_results)) {
      return(div(class = "alert alert-info",
                 icon("info-circle"),
                 " Run a Tumor Growth analysis first; the most recent two runs ",
                 "with different model types will appear here for comparison."))
    }
    if (is.null(rv$tg_results_prev)) {
      return(div(class = "alert alert-info",
                 icon("info-circle"),
                 " Only one analysis result is stored so far. Run again with a ",
                 "different ", strong("Statistical Model"),
                 " (e.g. switch from LME4 to Bayesian, or AUC to GAMM) to ",
                 "populate this comparison view."))
    }

    cur  <- rv$tg_results
    prev <- rv$tg_results_prev

    tag_for <- function(res, label) {
      te <- res$treatment_effects
      if (is.null(te)) {
        df_tag <- div(class = "alert alert-warning",
                      "No treatment_effects available for ", strong(label))
      } else {
        df <- format_df_for_display(as.data.frame(te))
        df_tag <- DT::datatable(
          df,
          options = list(scrollX = TRUE, paging = FALSE, searching = FALSE,
                         dom = "t"),
          rownames = FALSE
        )
      }
      div(
        h5(strong(label),
           tags$small(style = "color:#6c757d; margin-left:8px;",
                      paste("model_type_used:",
                            res$model_type_used %||% "unknown"))),
        df_tag
      )
    }

    fluidRow(
      column(6, tag_for(prev,
                        sprintf("Previous (%s)", prev$model_type_used %||% "?"))),
      column(6, tag_for(cur,
                        sprintf("Current (%s)", cur$model_type_used %||% "?")))
    )
  })

  # AUC individual data table
  output$tg_auc_individual_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$auc_analysis$individual)) {
      return(datatable(data.frame(Message = "AUC individual data not available")))
    }
    
    auc_ind_df <- as.data.frame(rv$tg_results$auc_analysis$individual)
    auc_ind_df <- format_df_for_display(auc_ind_df)
    datatable(auc_ind_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # AUC summary table
  output$tg_auc_summary_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$auc_analysis$summary)) {
      return(datatable(data.frame(Message = "AUC summary not available")))
    }
    
    auc_sum_df <- as.data.frame(rv$tg_results$auc_analysis$summary)
    auc_sum_df <- format_df_for_display(auc_sum_df)
    datatable(auc_sum_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Data summary table
  output$tg_data_summary_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$data_summary)) {
      return(datatable(data.frame(Message = "Data summary not available")))
    }
    
    data_sum_df <- as.data.frame(rv$tg_results$data_summary)
    data_sum_df <- format_df_for_display(data_sum_df)
    datatable(data_sum_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
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

  output$tg_posterior_summary <- renderDT({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    ps <- rv$tg_results$posterior_summary
    req(!is.null(ps))
    datatable(
      format_df_for_display(as.data.frame(ps)),
      options = list(scrollX = TRUE, paging = FALSE, searching = FALSE)
    )
  })

  output$tg_pp_check_plot <- renderPlot({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    req(!is.null(rv$bayes_models$tg))
    brms::pp_check(rv$bayes_models$tg, type = "dens_overlay", ndraws = 50)
  })

  output$tg_mcmc_diag_table <- renderDT({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    diag <- rv$tg_results$mcmc_diagnostics
    req(!is.null(diag))
    datatable(
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

  output$tg_bayes_influential <- renderDT({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    inf_df <- bayes_influential_obs(rv$tg_results)
    if (nrow(inf_df) == 0L) {
      return(datatable(
        data.frame(Note = "No observations exceed the Pareto-k = 0.7 threshold."),
        options = list(dom = "t", paging = FALSE, searching = FALSE), rownames = FALSE
      ))
    }
    datatable(
      inf_df,
      options = list(scrollX = TRUE, paging = FALSE, searching = FALSE), rownames = FALSE
    )
  })

  output$tg_mcmc_trace_plot <- renderPlot({
    req(rv$tg_results)
    req(isTRUE(rv$tg_results$model_type_used == "bayes"))
    req(!isTRUE(rv$tg_results$error))
    req(!is.null(rv$bayes_models$tg))
    bayesplot::mcmc_trace(brms::as.array(rv$bayes_models$tg))
  })

  # Necrosis table render + all download handlers live in
  # `mod_tumor_growth_downloads.R` (CODE_REVIEW.md C.1/C.3).
  add_tg_download_handlers(
    input, output, session, rv,
    tg_growth_plot, tg_rates_plot,
    tg_auc_plot, tg_forest_plot, tg_wf_plot
  )

  })
}