
#' Tumor Growth Analysis Server Logic
#' @keywords internal
tumor_growth_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
  
  rv$tg_results <- NULL
  rv$tg_running <- FALSE
  
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
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error: Please select a reference group before running analysis")
      })
      return()
    }
    
    rv$tg_running <- TRUE
    output$tg_status_msg <- renderUI({
      div(class = "alert alert-info",
          icon("spinner", class = "fa-spin"), "Running analysis...")
    })
    
    # Prepare data
    data <- rv$uploaded_data
    
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
          div(class = "alert alert-danger",
              icon("exclamation-circle"),
              "Please select both Length and Width columns before running analysis.")
        })
        rv$tg_running <- FALSE
        return()
      }
      if (!len_col %in% colnames(data) || !wid_col %in% colnames(data)) {
        output$tg_status_msg <- renderUI({
          div(class = "alert alert-danger",
              icon("exclamation-circle"),
              paste0("Column not found in data: '",
                     if (!len_col %in% colnames(data)) len_col else wid_col, "'"))
        })
        rv$tg_running <- FALSE
        return()
      }
      data <- mouseExperiment::calculate_volume(
        data,
        length_column = len_col,
        width_column  = wid_col,
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
    rv$cage_col       <- if (!is.null(rv$settings$col_cage) && rv$settings$col_cage != "") rv$settings$col_cage else NULL
    rv$dose_col       <- if (!is.null(rv$settings$col_dose) && rv$settings$col_dose != "") rv$settings$col_dose else NULL

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
      output$tg_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"),
            paste0("Required column(s) not found in data: ",
                   paste(names(missing_cols), "=", sQuote(missing_cols), collapse = ", ")))
      })
      rv$tg_running <- FALSE
      return()
    }
    
    # Run analysis with progress indicator
    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = "Running tumor growth analysis...", value = 0.1)
    on.exit(progress$close(), add = TRUE)

    result <- safe_analysis(module = "tumor_growth", {
      progress$set(value = 0.3, detail = "Fitting model")
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
        polynomial_degree = 1L,
        random_effects_specification = input$tg_random_effects,
        handle_cage_effects = input$tg_cage_effect,
        reference_group = input$tg_reference,  # Already validated to not be empty
        extrapolation_points = if (!is.null(input$tg_extrapolation_points) &&
                                    !is.na(input$tg_extrapolation_points))
                                 as.integer(input$tg_extrapolation_points) else 0L,
        return_model = TRUE,
        include_diagnostics = TRUE,
        plots = FALSE
      )
      progress$set(value = 0.7, detail = "Generating plots")
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
            group_summary_line = TRUE
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
                         x = "Theoretical Quantiles", y = "Sample Quantiles")
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
                           x = "Theoretical Quantiles", y = "Sample Quantiles")
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
                           x = "Fitted Values", y = "Residuals")
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

        if (!is.null(result$posthoc) && !is.null(result$posthoc$pairwise)) {
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

    rv$tg_results <- result
    rv$tg_running <- FALSE

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
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error:", result$message)
      })
    } else {
      output$tg_status_msg <- renderUI({
        div(class = "alert alert-success",
            icon("check"), "Analysis completed successfully!")
      })
    }
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
        paste("  Polynomial Degree:", summary_list$model_specification$polynomial_degree %||% 1L)
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
        strong("LME4 model summary: "), "Reports the fitted linear mixed-effects model including model type (LME4), data transformation applied, polynomial degree, random effects structure (subject and cage), number of subjects and treatment groups, reference group, and any convergence warnings or singular-fit notes."
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
          tags$li(strong("p.value"), " \u2014 adjusted p-value (Dunnett correction for \u2018vs Reference\u2019; Tukey HSD for \u2018All Comparisons\u2019). Values below 0.05 indicate a significant difference.")
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

    # Enforce column order: contrast, df, estimate, lower.CL, upper.CL, SE,
    # t.ratio, p.value, Fold_Change, Fold_Change_Lower, Fold_Change_Upper, then any remaining cols
    preferred_order <- c("contrast", "Contrast",
                         "df",
                         "estimate", "Estimate",
                         "lower.CL", "lower_CL",
                         "upper.CL", "upper_CL",
                         "SE",
                         "t.ratio", "z.ratio",
                         "p.value",
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
    datatable(rates_df,
      colnames = new_names,
      options  = list(scrollX = TRUE, paging = FALSE, searching = FALSE, scrollY = "500px", scrollCollapse = TRUE)
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
  
  # Group order and color combined control for Plots tab
  output$tg_group_order_color_ui <- renderUI({
    req(rv$tg_plot_groups)
    group_order_and_color_ui("tg", rv$tg_plot_groups, ns = session$ns)
  })

  # Growth curves plot  — reactive: re-renders on colour/option changes
  tg_growth_plot <- reactive({
    tryCatch({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    req(rv$tg_plot_data, rv$tg_col_names)
    data      <- rv$tg_plot_data
    col_names <- rv$tg_col_names
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    colors    <- build_color_vector(input, "tg_color", ordered_grps)
    # For dose data, apply factor level order on the original treatment column;
    # plot_tumor_growth will build Group = 'Treatment - Dose: X' internally.
    # For non-dose data, apply order on the treatment column directly.
    if (isTRUE(length(ordered_grps) > 0)) {
      if (!is.null(rv$dose_col) && rv$dose_col %in% colnames(data) &&
          !is.null(rv$settings$col_treatment) && rv$settings$col_treatment %in% colnames(data)) {
        # Derive treatment-level order from the combined group order
        tx_order <- unique(sub(" - Dose: .*$", "", ordered_grps))
        data[[rv$settings$col_treatment]] <- factor(data[[rv$settings$col_treatment]], levels = tx_order)
      } else if (isTRUE(col_names$treatment_column %in% colnames(data))) {
        data[[col_names$treatment_column]] <- factor(data[[col_names$treatment_column]],
                                                     levels = ordered_grps)
      }
    }

    surv_col_mapped <- local({
      col <- rv$settings$col_survival_event
      is.character(col) && length(col) == 1L &&
        nzchar(col) && col != "None" && isTRUE(col %in% colnames(data))
    })
    extra_pts <- local({
      ep <- input$tg_extrapolation_points
      if (is.numeric(ep) && length(ep) == 1L && !is.na(ep)) as.integer(ep) else 0L
    })
    do_extrapolate <- isTRUE(surv_col_mapped) && extra_pts > 0L && isTRUE(input$tg_opt_show_extrap)

    p <- tryCatch({
      mouseExperiment::plot_tumor_growth(
        df                  = data,
        volume_column       = col_names$volume_column,
        day_column          = col_names$time_column,
        # When a dose column is present, pass the original treatment column so
        # plot_tumor_growth can build 'Group = Treatment - Dose: X' itself,
        # matching the names used in rv$tg_plot_groups and the colour pickers.
        treatment_column    = if (!is.null(rv$dose_col)) rv$settings$col_treatment else col_names$treatment_column,
        dose_column         = rv$dose_col,
        ID_column           = col_names$id_column,
        cage_column         = col_names$cage_column,
        survival_column     = if (isTRUE(surv_col_mapped)) rv$settings$col_survival_event else NULL,
        extrapolate_volumes = do_extrapolate,
        extrapolation_points = if (extra_pts > 0L) extra_pts else "all",
        group_summary_line  = isTRUE(input$tg_opt_group_line),
        point_size          = input$tg_opt_growth_pt_size %||% 2
      ) + ggplot2::theme_classic()
    }, error = function(e) {
      if (!is.null(rv$tg_results$tumor_growth_plot)) return(rv$tg_results$tumor_growth_plot)
      NULL
    })
    req(!is.null(p))
    shapes <- build_shape_vector(input, "tg_shape", ordered_grps)
    if (!is.null(colors)) {
      p <- p + ggplot2::scale_color_manual(values = colors)
      # scale_fill_manual keeps dual-fill shapes (e.g. pch 25, down-triangle)
      # rendered with the correct group colour interior
      p <- p + ggplot2::scale_fill_manual(values = colors)
    }
    if (!is.null(shapes)) {
      p <- p + ggplot2::scale_shape_manual(values = shapes)
    }
    # Apply text size, title/subtitle, and legend position
    text_size  <- input$tg_opt_text_size %||% 16
    legend_pos <- input$tg_opt_legend_pos %||% "right"
    p <- apply_plot_labels(p, input$tg_plot_title %||% "", input$tg_plot_subtitle %||% "",
                           input$tg_plot_title_align %||% "center", text_size,
                           input$tg_plot_xaxis_label %||% "", input$tg_plot_yaxis_label %||% "",
                           input$tg_plot_legend_label %||% "")
    p <- apply_legend_pos(p, legend_pos)
    # Apply y-axis scale
    yscale <- input$tg_opt_yscale %||% "linear"
    if (yscale == "log10") p <- p + ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.05)))
    # Apply x-axis start and ensure 5% right padding.
    # coord_cartesian is always applied so the explicit x_upper value
    # (data max + 5%) is honoured regardless of scale expand settings.
    x_start <- local({
      v <- input$tg_opt_xstart
      if (is.numeric(v) && length(v) == 1L && !is.na(v) && v > 0) v else 0
    })
    x_upper <- tryCatch({
      time_col <- rv$tg_col_names$time_column
      x_max <- max(rv$tg_plot_data[[time_col]], na.rm = TRUE)
      x_max + 0.05 * (x_max - x_start)
    }, error = function(e) NA)
    if (!is.na(x_upper))
      p <- suppressMessages(p + ggplot2::coord_cartesian(xlim = c(x_start, x_upper)))
    # Annotate hollow circles = extrapolated when extrapolation is active
    if (isTRUE(do_extrapolate))
      p <- p + ggplot2::labs(caption = "Hollow circles / dashed lines = extrapolated data points after mouse death")
    p
    }, error = function(e) {
      msg <- conditionMessage(e)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
          label = paste0("Plot error:\n", msg),
          hjust = 0.5, vjust = 0.5, size = 4, color = "red") +
        ggplot2::theme_void()
    })
  })

  output$tg_plot_growth <- renderPlot({
    tg_growth_plot()
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)
  
  tg_rates_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    req(rv$tg_results$growth_rates)
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    colors     <- remap_colors_to_treatment(build_color_vector(input, "tg_color", ordered_grps))
    shapes     <- remap_shapes_to_treatment(build_shape_vector(input, "tg_shape", ordered_grps))
    center_bar <- input$tg_opt_center_bar %||% "mean"
    err_type   <- input$tg_opt_error_bar  %||% "SEM"
    pt_size    <- input$tg_opt_point_size %||% 2.5
    jit_width  <- input$tg_opt_jitter_width %||% 0.2
    text_size  <- input$tg_opt_text_size %||% 12
    legend_pos <- input$tg_opt_legend_pos %||% "auto"
    rates      <- rv$tg_results$growth_rates
    grp_col_r  <- intersect(c("Treatment", "Group"), colnames(rates))[1]
    req(!is.na(grp_col_r))
    val_col_r  <- intersect(c("growth_rate", "GrowthRate"), colnames(rates))[1]
    req(!is.na(val_col_r))
    # Derive treatment-level order (convert " - Dose: " → " " to match TreatmentGroup column)
    gr_order <- if (length(ordered_grps) > 0) {
      raw <- unique(gsub(" - Dose: ", " ", ordered_grps))
      raw[raw %in% as.character(unique(rates[[grp_col_r]]))]
    } else NULL
    if (!is.null(gr_order) && length(gr_order) > 0)
      rates[[grp_col_r]] <- factor(rates[[grp_col_r]], levels = gr_order)

    # Helper: build summary stats
    make_gr_summary <- function(df, grp, val) {
      grps <- if (is.factor(df[[grp]])) levels(df[[grp]]) else unique(as.character(df[[grp]]))
      do.call(rbind, lapply(grps, function(g) {
        x <- df[[val]][as.character(df[[grp]]) == g]
        x <- x[!is.na(x)]
        n <- length(x)
        ctr <- switch(center_bar,
          mean = mean(x), median = median(x),
          geometric_mean = exp(mean(log(pmax(x, 1e-9)))), mean(x))
        half <- switch(err_type,
          SEM = sd(x) / sqrt(n), SD = sd(x),
          CI = qt(0.975, df = max(n - 1L, 1L)) * sd(x) / sqrt(n), 0)
        data.frame(grp = g, ctr = ctr, lo = ctr - half, hi = ctr + half,
                   stringsAsFactors = FALSE)
      }))
    }
    summ_r <- make_gr_summary(rates, grp_col_r, val_col_r)
    summ_r$grp <- factor(summ_r$grp, levels = if (is.factor(rates[[grp_col_r]])) levels(rates[[grp_col_r]]) else unique(as.character(rates[[grp_col_r]])))

    # Build plot from scratch — ensures individual points are drawn LAST (on top)
    p <- ggplot2::ggplot(rates,
           ggplot2::aes(x = .data[[grp_col_r]], y = .data[[val_col_r]],
                        colour = .data[[grp_col_r]], fill = .data[[grp_col_r]],
                        shape = .data[[grp_col_r]])) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Treatment Group", y = "Growth Rate (per day)",
                    title = "Growth Rates by Treatment Group") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5)

    # Error bars
    if (err_type != "none")
      p <- p + ggplot2::geom_errorbar(
        data = summ_r,
        ggplot2::aes(x = grp, ymin = lo, ymax = hi),
        inherit.aes = FALSE,
        colour = "black", width = 0.25, linewidth = 0.9)

    # Center bar
    if (center_bar != "none")
      p <- p + ggplot2::geom_crossbar(
        data = summ_r,
        ggplot2::aes(x = grp, y = ctr, ymin = ctr, ymax = ctr),
        inherit.aes = FALSE,
        colour = "black", width = 0.35, linewidth = 0.9)

    # Individual points — drawn LAST so they appear on top of bars
    p <- p + ggplot2::geom_point(
      position = ggplot2::position_jitter(width = jit_width, height = 0, seed = 42),
      size = pt_size, alpha = 0.85)

    if (!is.null(colors)) p <- p + ggplot2::scale_colour_manual(values = colors)
    if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
    if (!is.null(shapes)) p <- p + ggplot2::scale_shape_manual(values = shapes)
    p <- apply_plot_labels(p, input$tg_plot_title %||% "", input$tg_plot_subtitle %||% "",
                           input$tg_plot_title_align %||% "center", text_size,
                           input$tg_plot_xaxis_label %||% "", input$tg_plot_yaxis_label %||% "",
                           input$tg_plot_legend_label %||% "")
    p <- apply_legend_pos(p, legend_pos)
    p
  })

  output$tg_plot_rates <- renderPlot({
    tg_rates_plot()
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)
  
  # AUC/Treatment effects plot — reactive
  tg_auc_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    # Shared controls — remap colours/shapes to treatment-level names
    # so scale_*_manual values match the Treatment column, even when
    # rv$tg_plot_groups contains dose-composite names.
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    colors     <- remap_colors_to_treatment(build_color_vector(input, "tg_color", ordered_grps))
    shapes     <- remap_shapes_to_treatment(build_shape_vector(input, "tg_shape", ordered_grps))
    center_bar <- input$tg_opt_center_bar %||% "mean"
    err_type   <- input$tg_opt_error_bar  %||% "SEM"
    pt_size    <- input$tg_opt_point_size %||% 2.5
    jit_width  <- input$tg_opt_jitter_width %||% 0.1
    text_size  <- input$tg_opt_text_size %||% 16
    legend_pos <- input$tg_opt_legend_pos %||% "right"

    # Helper: build a summary df (per group: center, lower, upper)
    make_summary <- function(df, grp_col, val_col) {
      groups <- if (is.factor(df[[grp_col]])) levels(df[[grp_col]]) else unique(as.character(df[[grp_col]]))
      do.call(rbind, lapply(groups, function(g) {
        x <- df[[val_col]][as.character(df[[grp_col]]) == g]
        x <- x[!is.na(x)]
        n <- length(x)
        ctr <- switch(center_bar,
          mean           = mean(x),
          median         = median(x),
          geometric_mean = exp(mean(log(pmax(x, 1e-9)))),
                           mean(x))
        half <- switch(err_type,
          SEM = sd(x) / sqrt(n),
          SD  = sd(x),
          CI  = qt(0.975, df = max(n - 1L, 1L)) * sd(x) / sqrt(n),
                0)
        data.frame(grp = g, ctr = ctr, lo = ctr - half, hi = ctr + half,
                   stringsAsFactors = FALSE)
      }))
    }

    # Helper: y-axis limits — floor at 0 (unless data dips below), 10% top buffer
    y_limits <- function(point_vals, lo_vals = NULL, hi_vals = NULL) {
      all_vals <- c(point_vals, lo_vals, hi_vals)
      all_vals <- all_vals[!is.na(all_vals)]
      ymin <- min(0, min(all_vals))
      ymax <- max(all_vals) * 1.10
      c(ymin, ymax)
    }

    # Helper: apply colours + labels + legend
    finish_plot <- function(p) {
      if (!is.null(colors)) p <- p + ggplot2::scale_colour_manual(values = colors)
      if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
      if (!is.null(shapes)) p <- p + ggplot2::scale_shape_manual(values = shapes)
      p <- apply_plot_labels(p, input$tg_plot_title %||% "", input$tg_plot_subtitle %||% "",
                             input$tg_plot_title_align %||% "center", text_size = text_size,
                             xlab = input$tg_plot_xaxis_label %||% "",
                             ylab = input$tg_plot_yaxis_label %||% "",
                             legend_lab = input$tg_plot_legend_label %||% "")
      apply_legend_pos(p, legend_pos)
    }

    # ── Path A: AUC model ── per-animal AUC values available ──────────────
    if (!is.null(rv$tg_results$auc_analysis) &&
        !is.null(rv$tg_results$auc_analysis$individual)) {

      auc_ind   <- rv$tg_results$auc_analysis$individual
      grp_col_a <- intersect(c("Treatment", "Group"), colnames(auc_ind))[1]
      req(!is.na(grp_col_a))

      # Factor ordering
      if (length(ordered_grps) > 0) {
        auc_levels <- unique(gsub(" - Dose: ", " ", ordered_grps))
        auc_levels <- auc_levels[auc_levels %in% as.character(unique(auc_ind[[grp_col_a]]))]
        if (length(auc_levels) > 0)
          auc_ind[[grp_col_a]] <- factor(auc_ind[[grp_col_a]], levels = auc_levels)
      }

      summ <- make_summary(auc_ind, grp_col_a, "AUC")
      summ$grp <- factor(summ$grp, levels = if (is.factor(auc_ind[[grp_col_a]])) levels(auc_ind[[grp_col_a]]) else unique(as.character(auc_ind[[grp_col_a]])))

      p <- ggplot2::ggplot(auc_ind,
             ggplot2::aes(x = .data[[grp_col_a]], y = .data[["AUC"]],
                          colour = .data[[grp_col_a]], fill = .data[[grp_col_a]],
                          shape = .data[[grp_col_a]])) +
        ggplot2::theme_classic() +
        ggplot2::labs(x = "Treatment", y = "AUC", title = "Area Under the Curve") +
        ggplot2::theme(legend.title = ggplot2::element_blank())

      if (err_type != "none")
        p <- p + ggplot2::geom_errorbar(
          data = summ,
          ggplot2::aes(x = grp, ymin = lo, ymax = hi),
          inherit.aes = FALSE,
          colour = "black", width = 0.25, linewidth = 0.9)

      if (center_bar != "none")
        p <- p + ggplot2::geom_crossbar(
          data = summ,
          ggplot2::aes(x = grp, y = ctr, ymin = ctr, ymax = ctr),
          inherit.aes = FALSE,
          colour = "black", width = 0.25, linewidth = 0.9)

      p <- p + ggplot2::geom_point(
        position = ggplot2::position_jitter(width = jit_width, height = 0, seed = 42),
        size = pt_size, alpha = 0.85)

      ylims <- y_limits(auc_ind[["AUC"]], summ$lo, summ$hi)
      p <- p + ggplot2::coord_cartesian(ylim = ylims)
      return(finish_plot(p))
    }

    # ── Path B: LME4 model ── same visual style as Path A, using growth rates ──
    rates <- rv$tg_results$growth_rates
    if (!is.null(rates) && nrow(rates) > 0) {
      grp_col_b <- intersect(c("Treatment", "Group"), colnames(rates))[1]
      val_col_b <- intersect(c("growth_rate", "GrowthRate"), colnames(rates))[1]
      if (!is.na(grp_col_b) && !is.na(val_col_b)) {
        gr_order <- if (length(ordered_grps) > 0) {
          raw <- unique(gsub(" - Dose: ", " ", ordered_grps))
          raw[raw %in% as.character(unique(rates[[grp_col_b]]))]
        } else NULL
        if (!is.null(gr_order) && length(gr_order) > 0)
          rates[[grp_col_b]] <- factor(rates[[grp_col_b]], levels = gr_order)

        summ_b <- make_summary(rates, grp_col_b, val_col_b)
        summ_b$grp <- factor(summ_b$grp, levels = if (is.factor(rates[[grp_col_b]])) levels(rates[[grp_col_b]]) else unique(as.character(rates[[grp_col_b]])))

        p <- ggplot2::ggplot(rates,
               ggplot2::aes(x = .data[[grp_col_b]], y = .data[[val_col_b]],
                            colour = .data[[grp_col_b]], fill = .data[[grp_col_b]],
                            shape = .data[[grp_col_b]])) +
          ggplot2::theme_classic() +
          ggplot2::labs(x = "Treatment", y = "Growth Rate (log-volume / day)",
                        title = "Growth Rate Comparison") +
          ggplot2::theme(legend.title = ggplot2::element_blank())

        if (err_type != "none")
          p <- p + ggplot2::geom_errorbar(
            data = summ_b,
            ggplot2::aes(x = grp, ymin = lo, ymax = hi),
            inherit.aes = FALSE,
            colour = "black", width = 0.25, linewidth = 0.9)

        if (center_bar != "none")
          p <- p + ggplot2::geom_crossbar(
            data = summ_b,
            ggplot2::aes(x = grp, y = ctr, ymin = ctr, ymax = ctr),
            inherit.aes = FALSE,
            colour = "black", width = 0.25, linewidth = 0.9)

        p <- p + ggplot2::geom_point(
          position = ggplot2::position_jitter(width = jit_width, height = 0, seed = 42),
          size = pt_size, alpha = 0.85)

        ylims <- y_limits(rates[[val_col_b]], summ_b$lo, summ_b$hi)
        p <- p + ggplot2::coord_cartesian(ylim = ylims)
        return(finish_plot(p))
      }
    }

    # ── Path C: last resort — stored plot ─────────────────────────────────
    p <- rv$tg_results$auc_plot
    if (is.null(p)) return(NULL)
    finish_plot(p)
  })

  output$tg_plot_auc <- renderPlot({
    tg_auc_plot()
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)
  
  # ── Waterfall Plot ────────────────────────────────────────────────────────
  tg_wf_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    rates <- rv$tg_results$growth_rates
    req(!is.null(rates), nrow(rates) > 0)

    tryCatch({
      wf_data <- rates |>
        dplyr::filter(!is.na(.data$growth_rate)) |>
        dplyr::mutate(
          animal_label = if ("Cage" %in% colnames(rates))
                           paste0(.data$Cage, "/", .data$ID)
                         else
                           as.character(.data$ID),
          pct_per_day = (exp(.data$growth_rate) - 1) * 100
        ) |>
        dplyr::arrange(dplyr::desc(.data$pct_per_day)) |>
        dplyr::mutate(rank = dplyr::row_number())

      req(nrow(wf_data) > 0)

      trt_col      <- "Treatment"
      ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
      colors       <- remap_colors_to_treatment(build_color_vector(input, "tg_color", ordered_grps))
      text_size    <- input$tg_opt_text_size %||% 16

      p <- ggplot2::ggplot(
             wf_data,
             ggplot2::aes(x = .data$rank,
                          y = .data$pct_per_day,
                          fill = .data[[trt_col]])) +
        ggplot2::geom_col(width = 0.8, colour = "grey30", linewidth = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.6) +
        { if (!is.null(colors)) ggplot2::scale_fill_manual(values = colors, name = trt_col)
          else ggplot2::scale_fill_discrete(name = trt_col) } +
        ggplot2::scale_x_continuous(
          breaks = wf_data$rank,
          labels = wf_data$animal_label,
          expand = c(0.01, 0)) +
        ggplot2::labs(
          x = "Animal",
          y = "Tumor growth rate (% change per day)") +
        ggplot2::theme_classic(base_size = text_size) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

      p <- apply_legend_pos(p, input$tg_opt_legend_pos %||% "right")
      apply_plot_labels(p,
        input$tg_plot_title    %||% "",
        input$tg_plot_subtitle %||% "",
        input$tg_plot_title_align %||% "center",
        text_size,
        input$tg_plot_xaxis_label %||% "",
        input$tg_plot_yaxis_label %||% "",
        input$tg_plot_legend_label %||% "")

    }, error = function(e) {
      NULL
    })
  })

  output$tg_plot_waterfall <- renderPlot({
    p <- tg_wf_plot()
    req(!is.null(p))
    p
  }, width  = function() input$tg_plot_width  %||% 800,
     height = function() input$tg_plot_height %||% 550)

  # Effect size forest plot — rebuilt from raw comparison data each render
  # so group order changes take effect without re-running analysis
  tg_forest_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    # Pick best comparison source (Dunnett preferred)
    using_dunnett <- !is.null(rv$tg_results$reference_comparisons_dunnett)
    pc <- if (using_dunnett) rv$tg_results$reference_comparisons_dunnett
          else rv$tg_results$all_pairwise_comparisons
    if (is.null(pc)) return(NULL)

    text_size    <- input$tg_opt_text_size %||% 12
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    grps_clean   <- unique(gsub(" - Dose: ", " ", ordered_grps))
    colors       <- remap_colors_to_treatment(build_color_vector(input, "tg_color", ordered_grps))
    shapes       <- remap_shapes_to_treatment(build_shape_vector(input, "tg_shape", ordered_grps))

    pc_df     <- as.data.frame(pc)
    label_col <- intersect(c("contrast", "Contrast", "comparison", "pair"), colnames(pc_df))[1]
    est_col   <- intersect(c("estimate", "Estimate"), colnames(pc_df))[1]
    se_col    <- intersect(c("SE"), colnames(pc_df))[1]
    lower_col <- intersect(c("lower.CL", "Lower", "lower_CL", "lower"), colnames(pc_df))[1]
    upper_col <- intersect(c("upper.CL", "Upper", "upper_CL", "upper"), colnames(pc_df))[1]
    p_col     <- intersect(c("p.value", "p_value", "p.adj", "adj.p.value", "P_Value"), colnames(pc_df))[1]
    if (is.na(label_col) || is.na(est_col)) return(NULL)

    if (is.na(lower_col) && !is.na(se_col)) {
      df_col <- intersect(c("df", "DF"), colnames(pc_df))[1]
      t_mult <- if (!is.na(df_col)) qt(0.975, pmax(as.numeric(pc_df[[df_col]]), 1)) else qnorm(0.975)
      pc_df$lower.CL <- pc_df[[est_col]] - t_mult * pc_df[[se_col]]
      pc_df$upper.CL <- pc_df[[est_col]] + t_mult * pc_df[[se_col]]
      lower_col <- "lower.CL"; upper_col <- "upper.CL"
    }

    p_label <- if (!is.na(p_col)) vapply(pc_df[[p_col]], format_pval, character(1)) else NULL
    ref_grp <- if (!is.null(input$tg_reference) && nzchar(input$tg_reference)) input$tg_reference else NULL

    contrast_labels <- as.character(pc_df[[label_col]])
    if (isTRUE(using_dunnett) && !is.null(ref_grp)) {
      suffix_pat <- paste0("\\s*[-\u2013]\\s*", gsub("([+()])", "\\\\\\1", ref_grp), "$")
      clean_tmp  <- sub(suffix_pat, "", contrast_labels)
      if (all(clean_tmp != contrast_labels)) contrast_labels <- clean_tmp
    }

    forest_df <- data.frame(
      Contrast = contrast_labels,
      Estimate = as.numeric(pc_df[[est_col]]),
      Lower    = if (!is.na(lower_col)) as.numeric(pc_df[[lower_col]]) else as.numeric(pc_df[[est_col]]),
      Upper    = if (!is.na(upper_col)) as.numeric(pc_df[[upper_col]]) else as.numeric(pc_df[[est_col]]),
      P_Label  = if (is.null(p_label)) "" else paste0("p=", p_label),
      is_ref   = FALSE,
      stringsAsFactors = FALSE
    )
    if (!is.null(ref_grp)) {
      forest_df <- rbind(forest_df, data.frame(
        Contrast = paste0(ref_grp, "\n(Reference)"),
        Estimate = 0, Lower = NA_real_, Upper = NA_real_,
        P_Label = "", is_ref = TRUE,
        stringsAsFactors = FALSE
      ))
    }

    # Sort non-reference rows by user group order; reference always at bottom of plot
    non_ref <- forest_df[!forest_df$is_ref, , drop = FALSE]
    ref_row <- forest_df[ forest_df$is_ref, , drop = FALSE]
    if (length(grps_clean) > 0) {
      rank_by_grp <- match(non_ref$Contrast, grps_clean)
      rank_by_grp[is.na(rank_by_grp)] <- .Machine$integer.max
      non_ref <- non_ref[order(rank_by_grp), , drop = FALSE]
    } else {
      non_ref <- non_ref[order(-non_ref$Estimate), , drop = FALSE]
    }
    # Factor levels bottom-to-top: reference at bottom, first group-order item at top
    forest_df <- rbind(ref_row, non_ref)
    forest_df$Contrast <- factor(
      forest_df$Contrast,
      levels = c(as.character(ref_row$Contrast), rev(as.character(non_ref$Contrast)))
    )

    show_pval   <- isTRUE(input$tg_opt_forest_show_pval   %||% TRUE)
    show_effect <- isTRUE(input$tg_opt_forest_show_effect %||% FALSE)
    show_ci     <- isTRUE(input$tg_opt_forest_show_ci     %||% FALSE)
    any_ann     <- show_pval || show_effect || show_ci

    # Build per-row annotation label from whichever elements are enabled
    forest_df$Ann_Label <- vapply(seq_len(nrow(forest_df)), function(j) {
      if (forest_df$is_ref[j]) return("")
      parts <- character(0)
      if (show_effect)
        parts <- c(parts, paste0("est: ", signif(forest_df$Estimate[j], 3)))
      if (show_ci && !is.na(forest_df$Lower[j]))
        parts <- c(parts,
          paste0("[", signif(forest_df$Lower[j], 3), ", ", signif(forest_df$Upper[j], 3), "]"))
      if (show_pval && nzchar(forest_df$P_Label[j]))
        parts <- c(parts, forest_df$P_Label[j])
      paste(parts, collapse = "\n")
    }, character(1L))

    # Build per-contrast colour/shape maps (contrast labels correspond to group names)
    contrast_lvls <- setdiff(as.character(levels(forest_df$Contrast)),
                             as.character(ref_row$Contrast))
    forest_colours <- if (!is.null(colors)) {
      vals <- ifelse(contrast_lvls %in% names(colors), unname(colors[contrast_lvls]), "steelblue")
      setNames(vals, contrast_lvls)
    } else setNames(rep("steelblue", length(contrast_lvls)), contrast_lvls)
    forest_shapes <- if (!is.null(shapes)) {
      vals <- ifelse(contrast_lvls %in% names(shapes), unname(shapes[contrast_lvls]), 16L)
      setNames(as.integer(vals), contrast_lvls)
    } else setNames(rep(16L, length(contrast_lvls)), contrast_lvls)

    all_vals <- c(non_ref$Estimate, non_ref$Lower, non_ref$Upper)
    max_abs  <- ceiling(max(abs(all_vals), na.rm = TRUE))
    # Widen right margin when annotations are shown so text is not clipped
    x_lim    <- c(-max_abs, max_abs * if (any_ann) 1.55 else 1.20)
    x_breaks <- pretty(c(-max_abs, max_abs), n = 6)
    x_breaks <- x_breaks[x_breaks >= -max_abs & x_breaks <= max_abs]

    x_lab <- if (isTRUE(using_dunnett)) "Mean Difference (treatment \u2212 reference)" else "Mean Difference"
    p <- ggplot2::ggplot(forest_df, ggplot2::aes(x = Estimate, y = Contrast)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.6) +
      ggplot2::geom_errorbarh(
        data = ~ subset(., !is_ref),
        ggplot2::aes(xmin = Lower, xmax = Upper, colour = Contrast),
        height = 0.25) +
      ggplot2::geom_point(data = ~ subset(., !is_ref),
        ggplot2::aes(colour = Contrast, shape = Contrast),
        size = input$tg_opt_forest_pt_size %||% 3) +
      ggplot2::geom_point(data = ~ subset(., is_ref),
        size = input$tg_opt_forest_pt_size %||% 3, colour = "grey50", shape = 18) +
      ggplot2::geom_text(
        data = ~ subset(., !is_ref & nzchar(Ann_Label) & !is.na(Upper)),
        ggplot2::aes(label = Ann_Label, x = Upper),
        hjust = -0.15, vjust = 0.5,
        size = (input$tg_opt_forest_ann_size %||% 12) / ggplot2::.pt,
        lineheight = 0.9) +
      ggplot2::scale_colour_manual(values = forest_colours, guide = "none") +
      ggplot2::scale_shape_manual(values = forest_shapes, guide = "none") +
      ggplot2::scale_x_continuous(limits = x_lim, breaks = x_breaks) +
      ggplot2::labs(x = x_lab, y = NULL, title = "Effect Size Forest Plot") +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
    apply_plot_labels(p,
      input$tg_plot_title    %||% "",
      input$tg_plot_subtitle %||% "",
      input$tg_plot_title_align %||% "center",
      text_size,
      input$tg_plot_xaxis_label  %||% "",
      input$tg_plot_yaxis_label  %||% "",
      input$tg_plot_legend_label %||% "")
  })

  output$tg_plot_effect_forest <- renderPlot({
    p <- tg_forest_plot()
    if (!is.null(p)) {
      p
    } else {
      plot(1, 1, type="n", main="Effect size forest plot not available\n(requires pairwise comparison data)", xlab="", ylab="")
    }
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)
  
  # Diagnostics
  # Diagnostics tab: Q-Q plot
  output$tg_diag_qq <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    if (!is.null(rv$tg_results$diag_qq_plot)) {
      return(rv$tg_results$diag_qq_plot)
    }
    plot(1, 1, type = "n", main = "Q-Q plot not available", xlab = "", ylab = "")
  })

  # Diagnostics tab: Residuals vs Fitted
  output$tg_diag_resid_fitted <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    if (!is.null(rv$tg_results$diag_resid_fitted_plot)) {
      return(rv$tg_results$diag_resid_fitted_plot)
    }
    plot(1, 1, type = "n", main = "Residuals vs Fitted not available", xlab = "", ylab = "")
  })
  
  # Download handlers (placeholder implementations)
  output$download_tg_anova <- downloadHandler(
    filename = "anova_results.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$anova)) {
        write.csv(as.data.frame(rv$tg_results$anova), file, row.names = TRUE)
      }
    }
  )
  
  output$download_tg_comparisons <- downloadHandler(
    filename = "pairwise_comparisons.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$pairwise_comparisons)) {
        write.csv(as.data.frame(rv$tg_results$pairwise_comparisons), file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_effects <- downloadHandler(
    filename = "treatment_effects.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$treatment_effects)) {
        write.csv(as.data.frame(rv$tg_results$treatment_effects), file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_rates <- downloadHandler(
    filename = "growth_rates.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$growth_rates)) {
        df <- as.data.frame(rv$tg_results$growth_rates)
        if ("growth_rate" %in% colnames(df)) {
          df$daily_fold_change  <- exp(df$growth_rate)
          df$pct_change_per_day <- (exp(df$growth_rate) - 1) * 100
        }
        write.csv(df, file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_auc_individual <- downloadHandler(
    filename = "auc_individual_data.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$auc_analysis$individual)) {
        write.csv(as.data.frame(rv$tg_results$auc_analysis$individual), file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_auc_summary <- downloadHandler(
    filename = "auc_summary.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$auc_analysis$summary)) {
        write.csv(as.data.frame(rv$tg_results$auc_analysis$summary), file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_data_summary <- downloadHandler(
    filename = "data_summary.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$data_summary)) {
        write.csv(as.data.frame(rv$tg_results$data_summary), file, row.names = FALSE)
      }
    }
  )
  
  # Plot download handlers - save as PNG
  output$download_tg_plot_growth <- downloadHandler(
    filename = "tumor_growth_curves.png",
    content = function(file) {
      p <- tg_growth_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_rates <- downloadHandler(
    filename = "growth_rate_distribution.png",
    content = function(file) {
      p <- tg_rates_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_auc <- downloadHandler(
    filename = "treatment_effects.png",
    content = function(file) {
      p <- tg_auc_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_effect_forest <- downloadHandler(
    filename = "effect_size_forest_plot.png",
    content = function(file) {
      p <- tg_forest_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 900,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )

  output$download_tg_plot_waterfall <- downloadHandler(
    filename = "waterfall_plot.png",
    content = function(file) {
      p <- tg_wf_plot()
      if (!is.null(p)) {
        png(file, width  = input$tg_plot_width  %||% 800,
                  height = input$tg_plot_height %||% 550,
                  units  = "px", res = input$tg_plot_dpi %||% 100)
        tryCatch(print(p), finally = dev.off())
      }
    }
  )

  })
}