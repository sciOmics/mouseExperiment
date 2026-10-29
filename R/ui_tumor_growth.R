# The Tumor Growth tab's UI (split from ui_modules.R, CODE_REVIEW.md R22-K).

#' Tumor Growth Analysis UI
#' @keywords internal
tumorGrowthUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Tumor Growth Analysis"),
    
    p("Configure analysis parameters and run statistical modeling. ",
      "Ensure data is uploaded first from the Data Upload tab."),
    
    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Analysis Parameters"),

          # 1. Statistical Model
          # Explanations sit in help_tip() tooltips, so the controls and the
          # Run button stay in view (26.10.007); the Info tab has the detail.
          selectInput(ns("tg_model_type"),
                     tagList("Statistical Model", help_tip(paste(
                       "LME4: a linear mixed-effects model of each animal's growth over time.",
                       "AUC: compares total tumour burden, the area under each arm's fitted",
                       "curve, over a window that every arm covers; animals that die early",
                       "still count through all their measurements. Bayesian LMM: the LME4",
                       "model with posterior distributions and credible intervals; it runs in",
                       "the background and takes seconds to a few minutes."))),
                     choices = c("Linear Mixed-Effects (LME4)"           = "lme4",
                                "Area Under Curve (AUC)"                 = "auc",
                                "Bayesian LMM (brms)"                    = "bayes"),
                     selected = "lme4"),

          # 2. Reference Group: the study's control group (R20-D.H)
          uiOutput(ns("tg_reference_note")),

          # 3. Post-hoc Comparisons
          selectInput(ns("tg_comparison_type"),
                     tagList("Post-hoc Comparisons", help_tip(paste(
                       "vs Reference Group compares each arm with the control group.",
                       "All Comparisons compares every pair of arms."))),
                     choices = c("vs Reference Group" = "reference",
                                "All Comparisons" = "tukey"),
                     selected = "reference"),
          # CODE_REVIEW.md R3.D3 — p_adjust_method was documented by the backend
          # with four options and exposed by the UI with none, so users got
          # whichever correction each code path happened to hardcode. The
          # backend now honours the argument (backend R3.1), so surface it.
          # Backend H.3 -- opt-in because each permutation refits the model, and
          # because the unit has to be a deliberate statement about the design
          # rather than a default someone accepts without reading.
          checkboxInput(ns("tg_run_perm"),
                       tagList("Run randomisation test of treatment \u00d7 time", help_tip(paste(
                         "A permutation test of whether the arms' growth rates differ, which",
                         "needs neither the model's degrees-of-freedom approximation nor",
                         "normal random effects. Each permutation refits the model, so it",
                         "takes longer."))),
                       value = FALSE),
          conditionalPanel(ns = ns, condition = "input.tg_run_perm == true",
            uiOutput(ns("tg_perm_unit_note")),
            numericInput(ns("tg_perm_n"), "Permutations",
                        value = 999, min = 99, max = 9999, step = 100)
          ),

          selectInput(ns("tg_p_adjust"),
                     tagList("Multiple-comparison correction", help_tip(paste(
                       "Applied over exactly the comparisons shown. Dunnett (vs reference)",
                       "and Tukey (all pairs) use the correlation between comparisons, so",
                       "they are less conservative than Bonferroni or Holm. The AUC model",
                       "offers Holm, Bonferroni, FDR or none."))),
                     choices = tg_p_adjust_choices("lme4"),
                     selected = "dunnett"),

          # 4. Data Transformation (LME4 only)
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'lme4'",
            selectInput(ns("tg_transform"),
                       tagList("Data Transformation", help_tip(paste(
                         "Applied to tumour volume before modelling. Log is recommended for",
                         "most data: it stabilises the variance and makes exponential growth",
                         "linear."))),
                       choices = c("Log" = "log",
                                  "Square Root" = "sqrt",
                                  "None" = "none"),
                       selected = "log")
          ),

          # 5. Extrapolation Points moved to the Plots tab, beside "Show
          # extrapolated curves/points": it changes only the plots (R23.9).

          # 6. Cage Effect Handling (LME4 only)
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'lme4'",
            selectInput(ns("tg_cage_effect"),
                       tagList("Cage Effect Handling", help_tip(paste(
                         "How cage enters the model. Auto follows the design: a fixed effect",
                         "when cages mix treatments, a random intercept when each arm has two",
                         "or more cages of its own, and left out with a warning when an arm",
                         "has a single cage. The legacy option drops nested cages, which",
                         "understates standard errors."))),
                       choices = c("Auto (match the design)" = "auto",
                                  "Include if not collinear (legacy)" = "include_if_not_collinear",
                                  "Always include" = "always_include",
                                  "Never include" = "never_include",
                                  "As random effect" = "as_random_effect"),
                       selected = "auto")
          ),

          # 7. Random Effects (LME4 and Bayesian)
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'lme4' || input.tg_model_type == 'bayes'",
            # CODE_REVIEW.md R20.5 / R20-D.H: random slopes by default. With a
            # random intercept only, animals' differing growth rates land in
            # the residual and the Treatment x Day test rejected a true null in
            # 73 % of simulated studies.
            selectInput(ns("tg_random_effects"),
                       tagList("Random Effects", help_tip(paste(
                         "Random slopes give each animal its own baseline and growth rate.",
                         "Intercept only assumes every animal in an arm grows at the same",
                         "rate; when they do not, the treatment \u00d7 time test finds",
                         "differences that are not there. If a model cannot be fitted, the",
                         "analysis falls back to simpler random effects and says so. Bayesian",
                         "mode treats Uncorrelated as Random slopes and None as Intercept only."))),
                       choices = c("Random slopes (recommended)" = "slope",
                                   "Uncorrelated random slopes" = "slope_uncorrelated",
                                   "Intercept only (all animals in an arm grow at the same rate)" = "intercept_only",
                                   "None" = "none"),
                       selected = "slope")
          ),

          # 8. Bayesian-specific controls
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'bayes'",
            hr(),
            h5(strong("Bayesian Settings")),
            selectInput(ns("tg_prior_strength"),
              tagList("Prior Strength", help_tip(paste(
                "The prior on the treatment effects. Skeptical, N(0, 0.25), needs strong",
                "data to support a large effect; Weakly informative is N(0, 1);",
                "Informative N(0, 0.5); Diffuse N(0, 2.5). Manual sets each prior below."))),
              choices = c(
                "Skeptical (recommended)" = "skeptical",
                "Weakly informative"      = "weakly_informative",
                "Informative"             = "informative",
                "Diffuse"                 = "diffuse",
                "Manual specification"    = "manual"
              ),
              selected = "skeptical"
            ),
            # Manual prior specification — shown only when "Manual" is selected
            conditionalPanel(ns = ns,
              condition = "input.tg_prior_strength == 'manual'",
              tags$hr(style = "margin:8px 0"),
              tags$p(tags$strong("Fixed effects (b)"), style = "margin-bottom:4px; font-size:0.9em;"),
              fluidRow(
                column(6, selectInput(ns("tg_prior_b_dist"), NULL,
                  choices = c("normal" = "normal", "student_t" = "student_t"),
                  selected = "normal")),
                column(3, numericInput(ns("tg_prior_b_mu"),  "\u03bc",  value = 0,    step = 0.05)),
                column(3, numericInput(ns("tg_prior_b_sd"),  "\u03c3",  value = 0.25, step = 0.05, min = 0.01))
              ),
              conditionalPanel(ns = ns,
                condition = "input.tg_prior_b_dist == 'student_t'",
                numericInput(ns("tg_prior_b_df"), "df (student_t)", value = 3, min = 1, step = 1)
              ),
              tags$p(tags$strong("Intercept"), style = "margin:6px 0 4px; font-size:0.9em;"),
              fluidRow(
                column(6, selectInput(ns("tg_prior_int_dist"), NULL,
                  choices = c("normal" = "normal", "student_t" = "student_t"),
                  selected = "normal")),
                column(3, numericInput(ns("tg_prior_int_mu"), "\u03bc",  value = 0,    step = 0.1)),
                column(3, numericInput(ns("tg_prior_int_sd"), "\u03c3",  value = 0.625, step = 0.05, min = 0.01))
              ),
              conditionalPanel(ns = ns,
                condition = "input.tg_prior_int_dist == 'student_t'",
                numericInput(ns("tg_prior_int_df"), "df (student_t)", value = 3, min = 1, step = 1)
              ),
              tags$p(tags$strong("Random-effect SD (sd)"), style = "margin:6px 0 4px; font-size:0.9em;"),
              fluidRow(
                column(6, selectInput(ns("tg_prior_sd_dist"), NULL,
                  choices = c("exponential" = "exponential", "half-normal" = "half_normal"),
                  selected = "exponential")),
                column(6, numericInput(ns("tg_prior_sd_param"),
                  label = uiOutput(ns("tg_prior_sd_param_label")),
                  value = 2, min = 0.01, step = 0.1))
              ),
              tags$p(tags$strong("Residual SD (sigma)"), style = "margin:6px 0 4px; font-size:0.9em;"),
              fluidRow(
                column(6, selectInput(ns("tg_prior_sigma_dist"), NULL,
                  choices = c("exponential" = "exponential", "half-normal" = "half_normal"),
                  selected = "exponential")),
                column(6, numericInput(ns("tg_prior_sigma_param"),
                  label = uiOutput(ns("tg_prior_sigma_param_label")),
                  value = 2, min = 0.01, step = 0.1))
              )
            ),
            numericInput(ns("tg_bayes_chains"), "MCMC chains",
                         value = 4, min = 1, max = 8, step = 1),
            numericInput(ns("tg_bayes_warmup"), "Warmup per chain",
                         value = 1000, min = 200, max = 5000, step = 100),
            numericInput(ns("tg_bayes_iter"), "Post-warmup draws per chain",
                         value = 500, min = 200, max = 5000, step = 100),
            # R23.3: no Stan backend menu; the dashboard fits with rstan.
            # R23.10: the note on compiling and queueing is in Help.
          ),

          # Necrotic handling: rendered server-side, shown once a necrotic column
          # is mapped on Data Upload (R20-D.V)
          uiOutput(ns("necrotic_handling_ui")),

          br(),

          actionButton(ns("run_tg_analysis"),
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),

          uiOutput(ns("tg_status_msg"))
        )
      ),
      
      tags$div(class = "col-sm-9 main-panel-col",
        tags$button(
          class = "btn btn-default btn-xs sidebar-toggle-btn",
          type  = "button",
          title = "Toggle options panel",
          tags$i(class = "fa fa-chevron-left fa-fw"), " Collapse"
        ),
        tabsetPanel(
          id = ns("tg_results_tabs"),

          tabPanel("Info",
            fluidRow(
              column(10, offset = 1,
                br(),

                # ── Overview ────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Overview")),
                  div(class = "panel-body",
                    p("The Tumor Growth tab provides three statistical modeling approaches selectable from the ",
                      strong("Statistical Model"), " dropdown:"),
                    tags$ul(
                      tags$li(strong("Linear Mixed-Effects Model (LME4)"), " \u2014 fits a mixed-effects model to repeated tumor volume measurements. Accounts for inter-animal variability and unequal observation times. Recommended for most longitudinal datasets."),
                      tags$li(strong("Area Under Curve (AUC)"), " \u2014 compares each arm\u2019s total tumor burden as the area under its fitted growth curve, over the window in which every arm still has enough animals on study."),
                      tags$li(strong("Bayesian LMM (brms)"), " \u2014 the LME4 model fitted by MCMC, returning full posterior distributions instead of point estimates and p-values.")
                    ),
                    p(tags$small(class = "text-muted",
                      "R14.3: this tab does ", strong("not"), " report Tumor Growth Inhibition. ",
                      "TGI is an endpoint summary \u2014 a ratio of two group means on one day \u2014 ",
                      "while every model here estimates a growth ", em("trajectory"), ", and ",
                      "collapsing a trajectory to a single-day ratio discards most of what was ",
                      "measured. TGI is computed where it belongs to the question being asked: ",
                      "the ", strong("Toxicity"), " tab (as the efficacy half of the therapeutic ",
                      "window) and the ", strong("Drug Synergy"), " tab (as the fractional effect ",
                      "the Bliss model acts on). Both expose an endpoint-estimand ",
                      "selector, because a ratio of endpoint means is exactly where early ",
                      "removals bias the answer."
                    )),
                    div(class = "alert alert-info", style = "font-size:0.9em; margin-bottom:0;",
                      icon("circle-info"), " ",
                      strong("Which endpoint?"), " Animals removed before study end make \u201cthe final timepoint\u201d ",
                      "ambiguous, and taking it literally computes TGI among survivors only \u2014 which flatters ",
                      "the control arm, because the fastest-growing controls are exactly the ones that left. ",
                      "The Toxicity and Synergy tabs expose an ", strong("Endpoint estimand"), " selector for this; ",
                      "here TGI is derived from the fitted model, which uses every observation each animal ",
                      "contributed rather than only those present at the end."
                    )
                  )
                ),

                # ── LME4 Mode ───────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-info", "LME4 Mode"), " \u2014 Analyses & Outputs")),
                  div(class = "panel-body",

                    h5(strong("Model")),
                    p("A linear mixed-effects model is fitted to transformed tumor volumes:"),
                    tags$blockquote(tags$code("log(Volume) ~ Treatment \u00d7 Day + (Day | Animal)")),
                    p("The Treatment \u00d7 Day interaction captures whether treatment groups differ in their growth ", em("trajectories"), " over time, not just at a single endpoint."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Transformation"), tags$dd("Log (default) stabilises variance and linearises exponential growth. Square root suits moderate heteroscedasticity. Use 'None' only for already-linear data."),
                      tags$dt("Random Effects"), tags$dd("'Random slopes' (default) gives each animal its own baseline volume and its own growth rate. If that model cannot be fitted, the analysis falls back to uncorrelated slopes and then to an intercept only, and the Summary says so. 'Intercept only' assumes every animal in an arm grows at the same rate, which makes the treatment \u00d7 time test find differences that are not there. 'None' drops the animal-level term entirely and is appropriate only when each animal contributes a single observation."),
                      tags$dt("Cage Effect"),
                      tags$dd(strong("Auto (default)"), " classifies the housing design and applies the matching treatment, because the correct handling depends on the design rather than on preference:",
                              tags$ul(
                                tags$li(em("Crossed"), " \u2014 a cage holds more than one treatment. Cage is separable from treatment and enters as a fixed effect."),
                                tags$li(em("Nested, replicated"), " \u2014 one treatment per cage, several cages per arm. Cage nests within treatment as a random effect; ignoring it treats cagemates as independent and understates standard errors."),
                                tags$li(em("Nested, partly replicated"), " \u2014 some arms have several cages and others one. Cage stays in the model as a random effect, estimated from the replicated arms, and the Summary names the single-cage arms."),
                                tags$li(em("Nested, confounded"), " \u2014 one cage per treatment. Cage and treatment are the same variable and no model can separate them. The result reports the effect as inseparable rather than silently crediting it to treatment.")
                              ),
                              "The remaining options force a specific handling. ", em("Include if not collinear (legacy)"),
                              " is the pre-v0.12 behaviour, kept so earlier analyses can be reproduced.")
                    ),
                    effects_help("animal"),   # R23.21

                    hr(),
                    h5(strong("ANOVA Tab")),
                    p("F-tests with Satterthwaite degrees of freedom, reported as ",
                      code("NumDF"), " / ", code("DenDF"), " / ", code("F value"), " / ", code("Pr(>F)"),
                      ", in this order: whether the arms' growth rates differ (treatment \u00d7 time, the primary ",
                      "analysis target); whether the arms differ at the mean study day; and whether there is growth ",
                      "on average. A ", code("Tests"), " column states each question."),
                    p(tags$small(class = "text-muted",
                      "The treatment row compares the arms at the mean study day, not at day 0, so it does not ",
                      "change when time is counted from implantation instead of randomisation. With few animals ",
                      "per arm, the randomisation test below makes fewer assumptions."
                    )),

                    hr(),
                    h5(strong("Comparisons Tab")),
                    p("Pairwise group comparisons at the mean observation day using estimated marginal means (EMMs). ",
                      strong("Two independent settings"), " control this, and conflating them is the usual mistake:"),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Post-hoc Comparisons"),
                      tags$dd("Which comparisons are made \u2014 ", em("vs Reference Group"), " (each treatment against the reference only) or ", em("All Comparisons"), " (every pair). This defines the family."),
                      tags$dt("Correction"),
                      tags$dd("How that family is corrected: Dunnett (exact, multivariate-t), Tukey HSD, Bonferroni, Holm, Benjamini-Hochberg FDR, or None. Dunnett is the default and is the right partner for a vs-reference family \u2014 it accounts for the correlation between comparisons that share a control arm, which Bonferroni ignores and therefore over-corrects.")
                    ),
                    p(tags$small(class = "text-muted",
                      "The correction is applied over exactly the set of comparisons requested, and the table header ",
                      "states which correction was actually used rather than asserting one. Comparisons are ",
                      "marginalised at the mean study day; for differences in trajectory across time, see the ",
                      "Treatment Effects tab."
                    )),

                    hr(),
                    h5(strong("Randomisation Test")),
                    p("An optional exact test of the Treatment \u00d7 Day interaction that avoids both the ",
                      "denominator-degrees-of-freedom approximation and the assumption of normally distributed ",
                      "random effects. Treatment labels are permuted and the model refitted; the p-value is the ",
                      "proportion of permutations giving a statistic at least as extreme as the observed one."),
                    div(class = "alert alert-warning", style = "font-size:0.9em;",
                      icon("triangle-exclamation"), " ",
                      strong("You must state the unit that was randomised."), " If animals were assigned to arms ",
                      "individually, permute animals. If whole cages were assigned, permute cages \u2014 permuting ",
                      "animals then breaks the exchangeability the test relies on and yields a p-value that is too ",
                      "small. This cannot be inferred from the data; it is a fact about how the experiment was run."),
                    p(tags$small(class = "text-muted",
                      "The panel reports the number of exchangeable units, how many permutations converged, whether ",
                      "the permutation set was enumerated exhaustively, and the smallest p-value the design can ",
                      "attain. With few units that floor can exceed 0.05, in which case no amount of signal can ",
                      "produce a significant result \u2014 and the asymptotic p-value beside it should not be trusted."
                    )),

                    hr(),
                    h5(strong("Treatment Effects Tab")),
                    p("Estimated marginal means (EMMs) of transformed volume per treatment group, at the mean observation day: the day the treatment row of the ANOVA compares the arms at. How the arms diverge over time is in the growth-rate test (ANOVA) and the Growth Rates tab."),

                    hr(),
                    h5(strong("Growth Rates Tab")),
                    p("Per-animal exponential growth rates estimated by log-linear regression of each animal\u2019s volume trajectory independently. Reported as growth rate (log-scale slope), fold change, and % change per day with 95\u202f% CIs and R\u00b2. These are descriptive; for group-level inference use the ANOVA and Comparisons tabs."),

                    hr(),
                    h5(strong("Diagnostics Tab")),
                    p("Four residual plots, each answering a different question about whether the model's ",
                      "assumptions hold \u2014 which is what makes the p-values above mean anything:"),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Q\u2013Q plot"),
                      tags$dd("Are residuals normally distributed? Systematic curvature at the tails suggests the transformation choice should be revisited."),
                      tags$dt("Residuals vs fitted"),
                      tags$dd("Is the mean structure right? Curvature means log-linear growth does not hold over the whole study \u2014 read growth-rate comparisons with care."),
                      tags$dt("Scale\u2013location"),
                      tags$dd("Is the variance constant on the modelling scale? A rising trend means heteroscedasticity survived the transformation."),
                      tags$dt("Random-effects Q\u2013Q"),
                      tags$dd("Are the per-animal effects themselves normal? This is a separate assumption from residual normality and the one most often overlooked. If it fails, prefer the randomisation test, which assumes neither.")
                    ),
                    p(tags$small(class = "text-muted",
                      "Influence measures (Cook's distance, DFBETAS) are also available: they identify individual ",
                      "animals whose removal would materially change the conclusion. A result resting on one mouse ",
                      "is worth knowing about before publication rather than after."
                    ))
                  )
                ),

                # ── AUC Mode ─────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-warning", "AUC Mode"), " \u2014 Analyses & Outputs")),
                  div(class = "panel-body",
                    p("Each arm's AUC is the area under its fitted growth curve, the geometric-mean curve of a ",
                      "mixed model of log volume with a smooth (natural spline) time course per arm and ",
                      "per-animal random slopes. It is integrated over one window common to every arm: from the ",
                      "first study day to the last day on which every arm is evaluable (at least 50 % of its ",
                      "enrolled animals, and at least 3, still on study). AUC is in volume \u00d7 days."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("AUC by arm"), tags$dd("Each arm's AUC with a 95 % interval, from draws of the model."),
                      tags$dt("Comparisons"),
                      tags$dd("The ratio of two arms' AUCs and their difference, with 95 % intervals and a Wald ",
                              "p-value from the same draws, corrected over the family selected in the sidebar. ",
                              "A ratio below 1 means less tumor burden than the comparator."),
                      tags$dt("ANOVA"),
                      tags$dd("A Wald test that every arm has the same AUC (on the log scale)."),
                      tags$dt("Per-animal AUC"),
                      tags$dd("Each animal's own trapezoidal AUC, for description only.")
                    ),
                    div(class = "alert alert-info", style = "font-size:0.9em;",
                      icon("circle-info"), " ",
                      strong("Animals that die or are removed early "), "still count. Every measurement an ",
                      "animal had goes into the model, which estimates its growth rate from them and keeps it in ",
                      "the arm's curve as if the animal had continued. So the AUC is not pulled down when the ",
                      "largest tumours leave the study, as an average of the animals still on study would be. ",
                      "In simulations with 10 animals per arm and about a third of the treated arm lost, removal at ",
                      "the volume limit moved the AUC by \u22120.2 % (the average of the survivors: \u221213 %). ",
                      "Deaths that the earlier measurements only partly foreshadow left a small underestimate, ",
                      "\u22123 % to \u22125 % (survivors: \u221218 % to \u221221 %). When many animals die for ",
                      "reasons their growth does not show, read the Survival analysis alongside the AUC."
                    ),
                    div(class = "alert alert-warning", style = "font-size:0.9em; margin-bottom:0;",
                      icon("triangle-exclamation"), " ",
                      strong("Why not each animal's own AUC? "), "An animal euthanised on day 14 has a smaller ",
                      "trapezoidal AUC than one followed to day 28 because it was measured for less time, not because ",
                      "its tumour was smaller. Controls leave first, so comparing per-animal AUCs made an effective ",
                      "drug look worse than control in about half of simulated studies. The per-animal values are ",
                      "shown under \u2018AUC Individual Data\u2019 for description, not compared.")
                  )
                ),

                # ── Bayesian LMM Mode ─────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-success", "Bayesian Mode"), " \u2014 Method & Outputs")),
                  div(class = "panel-body",

                    h5(strong("Method")),
                    p("The Bayesian LMM fits the same formula as LME4 using full Bayesian inference via ",
                      tags$a("brms", href = "https://paul-buerkner.github.io/brms/", target = "_blank"),
                      " / Stan. Instead of point estimates and p-values, every parameter has a ",
                      strong("posterior distribution"), ": a probability distribution over plausible values ",
                      "given the data and priors."),
                    tags$blockquote(tags$code(
                      "log(Volume) ~ Treatment \u00d7 Day + (Day | Animal)"
                    )),
                    p("Results are summarised as ", strong("posterior means"),
                      " with 95 % credible intervals. A 95 % CrI has a direct probability ",
                      "interpretation: there is a 95 % posterior probability that the true ",
                      "parameter lies within the interval \u2014 unlike a frequentist confidence ",
                      "interval, which is a statement about the procedure rather than the parameter."),
                    p(tags$small(class = "text-muted",
                      "R14.7: two different interval types appear under the same ",
                      code("Lower_95_CrI"), " / ", code("Upper_95_CrI"), " headings, and it is ",
                      "worth knowing which you are reading. The ", strong("Summary"),
                      " tab reports raw model parameters straight from brms: an ",
                      em("equal-tailed"), " interval, cutting 2.5 % from each tail. The ",
                      strong("Treatment Effects"), " and ", strong("Comparisons"),
                      " tabs come from emmeans, which returns a genuine ", em("highest-density"),
                      " interval \u2014 the shortest interval containing 95 % of the posterior. ",
                      "The two coincide when a posterior is symmetric and diverge when it is ",
                      "skewed, which is exactly the case for variance components and for the ",
                      "back-transformed time ratios in the Bayesian survival model."
                    )),

                    hr(),
                    h5(strong("Prior Strength")),
                    p("Priors express what you believe about parameter values ", em("before"), " seeing the data. ",
                      "They regularise estimates and prevent the model from chasing noise in small samples. ",
                      "All presets use Normal priors on fixed effects and Exponential priors on variance components; only the widths differ."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Skeptical"),
                      tags$dd(tags$strong("Default."), " N(0, 0.25) on fixed effects; Exponential(2) on SDs. ",
                              "Encodes the reasonable belief that treatment effects are small \u2014 effects > \u00b10.5 log-units per day ",
                              "are considered very unlikely a priori. Recommended for exploratory preclinical studies where ",
                              "large effects have not been previously established."),
                      tags$dt("Weakly informative"),
                      tags$dd("N(0, 1) on fixed effects; Exponential(1) on SDs. ",
                              "Broader than Skeptical \u2014 appropriate when prior work suggests the treatment ",
                              "category can produce moderate-to-large effects."),
                      tags$dt("Informative"),
                      tags$dd("N(0, 0.5) on fixed effects; Exponential(2) on SDs. ",
                              "Between Skeptical and Weakly informative. Useful when effect sizes in the ",
                              "0.1\u20130.5 log-unit range are plausible based on prior data."),
                      tags$dt("Diffuse"),
                      tags$dd("N(0, 2.5) on fixed effects; Exponential(0.5) on SDs. ",
                              "Near-flat priors; use with caution on small datasets as estimates can be ",
                              "driven entirely by noisy data."),
                      tags$dt("Manual"),
                      tags$dd("Specify the prior family and parameters for each model class independently. ",
                              "Fixed effects (b): Normal or Student-t. Intercept: same. ",
                              "Random-effect SD and residual sigma: Exponential or half-Normal.")
                    ),

                    hr(),
                    h5(strong("Result Tabs")),
                    tags$ul(
                      tags$li(strong("Summary"), " \u2014 posterior median and 95 % CrI for every model parameter."),
                      tags$li(strong("Treatment Effects"), " \u2014 group-level estimated marginal means (posterior) at the mean study day."),
                      tags$li(strong("Comparisons"), " \u2014 posterior contrasts of each group vs the reference."),
                      tags$li(strong("Posterior Summary"), " \u2014 every model parameter's posterior estimate, 95 % credible interval, R-hat and effective sample sizes. It replaces the ANOVA tab for a Bayesian fit."),
                      tags$li(strong("Diagnostics: posterior predictive check"), " \u2014 overlay of observed data density vs draws from the posterior predictive distribution. Systematic mismatch indicates model mis-specification. A coverage plot reports whether the stated credible intervals actually contain the observed values at the claimed rate."),
                      tags$li(strong("Diagnostics"), " (Bayesian view) \u2014 Rhat and effective sample size (ESS) per parameter, with trace and rank plots. Rhat > 1.01 indicates the chains have not mixed; low ESS means the posterior summaries are themselves noisy. Rank plots detect the mixing failures that trace plots hide once chains overlap."),
                      tags$li(strong("Influence"), " \u2014 Pareto-k values from PSIS leave-one-out. Values above the threshold (0.7 at most; lower with fewer posterior draws, e.g. 0.667 at 1,000) mean the importance sampling is unreliable for that observation, which usually flags an animal the model fits badly rather than a coding error."),
                      tags$li(strong("Growth Rates"), " \u2014 per-animal log-linear regression rates (same as LME4, not MCMC-based).")
                    ),

                    div(class = "alert alert-info", style = "margin-top:10px;",
                      icon("clock"), " ",
                      strong("Computation time:"), " the first fit of a model compiles its Stan code (about two minutes on mouseexperiment.com, less on a fast computer); ",
                      "a later fit whose model has the same structure (same number of arms, random effects and prior preset) ",
                      "reuses the compiled model and starts sampling at once. ",
                      "Fits run one at a time on a background worker, so the app stays responsive while a model runs; ",
                      "the progress window shows how many fits are ahead of yours."
                    )
                  )
                ),

                # ── Necrotic Tumor Handling ───────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title",
                       tags$span(class = "label label-danger", "Necrotic Tumors"),
                       " \u2014 Handling Options")),
                  div(class = "panel-body",
                    p("When a ", strong("Necrotic Flag"), " column is mapped in Data Upload, a handling
                      selector appears in the sidebar. Necrotic tumors have a non-viable core that inflates
                      caliper-based volume estimates without contributing to proliferating tissue."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Exclude (recommended)"),
                      tags$dd("Removes flagged observations from the model input. Conservative and publication-defensible.
                               The Necrosis tab reports how many observations were flagged and which groups were affected,
                               regardless of handling mode."),
                      tags$dt("Covariate"),
                      tags$dd("Adds a 0/1 ", code(".necrotic_flag"), " fixed effect to the LMM formula. Estimates the
                               average volume bias attributable to necrosis. Useful when necrosis is approximately
                               evenly distributed across groups and you want to retain all timepoints."),
                      tags$dt("No adjustment"),
                      tags$dd("All observations are used as-is (current default when no necrotic column is mapped).
                               The Necrosis summary table is still shown for documentation purposes.")
                    ),
                    div(class = "alert alert-info", style = "font-size:0.85em; margin-bottom:0;",
                      icon("circle-info"),
                      " The Necrosis tab in the results panel is inserted automatically when a necrotic column
                        is mapped and an analysis has been run."
                    )
                  )
                ),

                # ── Shared Parameters ─────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Shared Parameters")),
                  div(class = "panel-body",
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Reference Group"),
                      tags$dd("The study's control group, chosen under Study Design on Data Upload. With a dose column mapped, arms are labelled 'Treatment Dose' and the control maps to its own dose label."),
                      tags$dt("Post-hoc Comparisons"),
                      tags$dd("Which family of comparisons to make \u2014 vs Reference Group (default) or All Comparisons. Paired with the separate Correction setting."),
                      tags$dt("Correction"),
                      tags$dd("How the chosen family is corrected: Dunnett (exact multivariate-t, default), Tukey HSD, Bonferroni, Holm, Benjamini-Hochberg FDR, or None."),
                      tags$dt("Extrapolation Points"),
                      tags$dd(strong("Affects the growth-trajectory plots only"), " (Plots tab, under Show extrapolated curves/points). Number of recent timepoints used to estimate a slope when projecting the curve for animals that left the study early; projected segments are drawn dashed. Set to 0 to disable. This no longer feeds any statistical model \u2014 treating an imputed point as measured data made standard errors ignore the imputation."),
                      tags$dt("Time Column"),
                      tags$dd("Numeric study day, or calculated from a date column mapped in Data Upload (days elapsed from first observation).")
                    )
                  )
                ),

                br()
              )
            )
          ),

          tabPanel("Summary",
            uiOutput(ns("tg_summary_tab_info")),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("tg_summary"))
            ),
            # CODE_REVIEW.md Round 2 F.1 — reproducible report export.
            br(),
            div(class = "text-right",
              downloadButton(ns("download_tg_report"),
                             "Download Full Report (HTML)",
                             class = "btn-sm btn-primary")
            )
          ),
          
          tabPanel("ANOVA",
            uiOutput(ns("tg_anova_tab_info")),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("tg_anova_table"))
            ),
            downloadButton(ns("download_tg_anova"), "Download ANOVA"),

            # R20.D30 / R20.D61: the Welch variance panel is gone; the
            # model-based AUC (backend 0.26.0) has no variance test.

            # Backend H.3 -- randomisation test of the interaction. Reported here
            # because it answers the same question as the ANOVA row above but
            # without the denominator-df or normality approximations.
            uiOutput(ns("tg_perm_test_ui"))
          ),

          # R23.4/R23.5: a Bayesian fit has no ANOVA; the server hides that tab
          # and shows this one in its place. The posterior summary sat in a
          # "Posterior Check" tab appended after Download Results.
          tabPanel("Posterior Summary",
            posterior_summary_help("log tumour volume"),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("tg_posterior_summary"))
            )
          ),
          
          tabPanel("Comparisons",
            uiOutput(ns("tg_comparisons_tab_info")),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("tg_comparisons_table"))
            ),
            downloadButton(ns("download_tg_comparisons"), "Download Comparisons")
          ),
          
          tabPanel("Treatment Effects",
            uiOutput(ns("tg_treatment_effects_tab_info")),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("tg_effects_table"))
            ),
            downloadButton(ns("download_tg_effects"), "Download Effects")
          ),
          
          tabPanel("Growth Rates",
            tab_info(
              strong("Per-animal exponential growth rates "), "estimated from each subject\u2019s individual tumor volume trajectory using linear regression on log\u2081\u209a(volume) vs time. Requires \u22653 time points per animal. ",
              strong("All rate, fold change, and percent change values are on a per-day basis."),
              tags$ul(
                tags$li(strong("Treatment"), " \u2014 treatment group."),
                tags$li(strong("ID"), " \u2014 individual animal identifier."),
                tags$li(strong("Cage"), " \u2014 cage identifier."),
                tags$li(strong("Growth Rate"), " \u2014 slope of the log\u2081\u209a(volume)\u2013time regression (units: log-volume per day). Represents the instantaneous exponential growth constant: a value of 0.10 means ~10% volume increase per day on a continuous compounding basis."),
                tags$li(strong("GR 95% CI Lower/Upper Bound"), " \u2014 95% confidence interval for the log-scale growth rate."),
                tags$li(strong("Fold Change"), " \u2014 exp(Growth Rate): multiplicative volume change per day. A value of 1.10 means the tumor is 10% larger each day; values below 1.0 indicate regression."),
                tags$li(strong("FC 95% CI Lower/Upper Bound"), " \u2014 exp(GR CI bounds): back-transformed confidence interval for the daily fold change."),
                tags$li(strong("Percent Change"), " \u2014 (Fold Change \u2212 1) \u00d7 100: daily volume change expressed as a percentage (e.g. +10.5%/day or \u22123.2%/day)."),
                tags$li(strong("PC 95% CI Lower/Upper Bound"), " \u2014 (exp(GR CI) \u2212 1) \u00d7 100: back-transformed CI bounds for the percent change."),
                tags$li(strong("R\u00b2"), " \u2014 coefficient of determination from the per-animal log-linear regression. Values close to 1.0 indicate a good exponential fit; lower values suggest irregular growth or too few time points.")
              ),
              em("These are descriptive per-animal values used to visualise within-group spread and for sample size estimation. For group-level inference, use the Treatment Effects tab.")
            ),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("tg_growth_rates_table"))
            ),
            downloadButton(ns("download_tg_rates"), "Download Growth Rates"),
            br(), br(),
            # Per-mouse drill-down (CODE_REVIEW.md Round 2 F.4)
            h5(strong(icon("magnifying-glass"), " Per-mouse drill-down"),
               tags$small(style = "color:#6c757d; margin-left:8px;",
                          "Click rows in the table above to see those mice's trajectories; click again to remove one.")),
            shinycssloaders::withSpinner(
              plotOutput(ns("tg_drilldown_plot"), height = "auto")
            )
          ),
          
          tg_endpoints_tab(ns),   # F1/F2

          tabPanel("Data Summary",
            tab_info(
              strong("Descriptive statistics for tumor volume "), "by treatment group and timepoint, on the volume scale in the declared units, whatever transform the model used.",
              tags$ul(
                tags$li(strong("Treatment"), " \u2014 treatment group."),
                tags$li(strong("Day"), " \u2014 study day (measurement timepoint)."),
                tags$li(strong("N"), " \u2014 number of animals with a volume measured on that day."),
                tags$li(strong("Mean"), " \u2014 arithmetic mean tumor volume for the group at that timepoint."),
                tags$li(strong("SD"), " \u2014 standard deviation of volume across animals in the group."),
                tags$li(strong("SEM"), " \u2014 standard error of the mean (SD / \u221aN)."),
                tags$li(strong("Geometric mean"), " (log transform only) \u2014 the mean on the log scale the model uses, back-transformed; zeros count at the floor stated in the Summary. Tumour volumes are right-skewed, so it is usually below the arithmetic mean.")
              ),
              em("Useful for inspecting raw descriptive trends across time and confirming data entry. Statistical comparisons should use the ANOVA, Comparisons, or AUC Summary tabs as appropriate.")
            ),
            shinycssloaders::withSpinner(
              DT::DTOutput(ns("tg_data_summary_table"))
            ),
            downloadButton(ns("download_tg_data_summary"), "Download Data Summary")
          ),
          
          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  selectInput(ns("tg_plot_selector"), "Plot Type",
                    choices = c(
                      "Tumor Growth Curves"      = "growth_curves",
                      "Growth Rates Plot"        = "growth_rate",
                      "Area Under the Curve Plot" = "auc",
                      "Waterfall Plot"           = "waterfall",
                      "Effect Size Forest Plot"  = "effect_forest"
                    )
                  ),
                  ctrl_section("group-order", "sort", "Groups",
                    uiOutput(ns("tg_group_order_color_ui"))
                  ),
                  ctrl_section("elements", "layer-group", "Elements",
                    conditionalPanel(ns = ns,
                      "['growth_rate','auc'].indexOf(input.tg_plot_selector) >= 0",
                      selectInput(ns("tg_opt_center_bar"), "Center line / bar",
                        choices = c("Mean" = "mean", "Median" = "median",
                                    "Geometric Mean" = "geometric_mean",
                                    "None" = "none"),
                        selected = "mean"
                      ),
                      selectInput(ns("tg_opt_error_bar"), "Error bars",
                        choices  = c("None" = "none", "SEM" = "SEM",
                                     "SD" = "SD", "95% CI" = "CI"),
                        selected = "SEM"
                      )
                    ),
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'growth_curves'",
                      checkboxInput(ns("tg_opt_group_line"),
                                    "Show group average line", value = TRUE),
                      checkboxInput(ns("tg_opt_show_extrap"),
                                    "Show extrapolated curves/points", value = TRUE),
                      # R23.9: from the sidebar, where it read as a model setting.
                      conditionalPanel(ns = ns, "input.tg_opt_show_extrap == true",
                        numericInput(ns("tg_extrapolation_points"),
                                    tagList("Extrapolation points", help_tip(paste(
                                      "The number of recent timepoints used to project the curve of an",
                                      "animal that left the study early, drawn dashed. 0 turns it off.",
                                      "The statistics never use projected values."))),
                                    value = 3, min = 0, max = 10, step = 1),
                        uiOutput(ns("tg_extrap_note"))
                      )
                    ),
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'effect_forest'",
                      checkboxInput(ns("tg_opt_forest_show_pval"),
                                    "Show p-value",     value = TRUE),
                      checkboxInput(ns("tg_opt_forest_show_effect"),
                                    "Show effect size", value = FALSE),
                      checkboxInput(ns("tg_opt_forest_show_ci"),
                                    "Show 95% CI",      value = FALSE)
                    ),
                    conditionalPanel(ns = ns,
                      "['growth_rate','auc','effect_forest'].indexOf(input.tg_plot_selector) >= 0",
                      checkboxInput(ns("tg_opt_backtransform"),
                                    "Display on original scale",
                                    value = FALSE),
                      helpText("Back-transforms log-scale model output to original units. ",
                               "Growth rates and contrasts become daily fold changes (exp(value)); ",
                               "forest plot uses a log10 axis centred on 1. ",
                               "Inference is still performed on the log scale.")
                    ),
                    selectInput(ns("tg_opt_legend_pos"), "Legend position",
                      choices = c("None" = "none",
                                  "Top" = "top", "Bottom" = "bottom",
                                  "Left" = "left", "Right" = "right",
                                  "In plot - top left"  = "inplot_tl",
                                  "In plot - top right" = "inplot_tr",
                                  "In plot - bottom left" = "inplot_bl",
                                  "In plot - bottom right" = "inplot_br"),
                      selected = "right"
                    ),
                    textInput(ns("tg_plot_title"), "Title",
                              placeholder = "Leave blank for default"),
                    selectInput(ns("tg_plot_title_align"), "Title alignment",
                      choices = c("Left" = "left", "Center" = "center", "Right" = "right"),
                      selected = "center"),
                    textInput(ns("tg_plot_subtitle"), "Subtitle",
                              placeholder = "Optional subtitle"),
                    textInput(ns("tg_plot_xaxis_label"), "X-axis label",
                              placeholder = "Leave blank for default"),
                    textInput(ns("tg_plot_yaxis_label"), "Y-axis label",
                              placeholder = "Leave blank for default"),
                    textInput(ns("tg_plot_legend_label"), "Legend label",
                              placeholder = "Leave blank for default")
                  ),
                  ctrl_section("scaling", "expand", "Scaling",
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'growth_curves'",
                      selectInput(ns("tg_opt_yscale"), "Y-axis scale",
                        choices = c("Linear" = "linear",
                                    "Log10"  = "log10"),
                        selected = "linear"
                      ),
                      numericInput(ns("tg_opt_xstart"), "X-axis start (Day)",
                                   value = 0, min = 0, step = 1)
                    ),
                    sliderInput(ns("tg_opt_text_size"), "Text size (pt)",
                                min = 10, max = 24, value = 16, step = 1),
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'growth_curves'",
                      sliderInput(ns("tg_opt_growth_pt_size"), "Point size",
                                  min = 0.5, max = 10, value = 4, step = 0.5)
                    ),
                    conditionalPanel(ns = ns,
                      "['growth_rate','auc'].indexOf(input.tg_plot_selector) >= 0",
                      sliderInput(ns("tg_opt_point_size"),   "Point size",
                                  0.5, 10, 4, step = 0.5),
                      sliderInput(ns("tg_opt_jitter_width"), "Jitter width",
                                  0, 0.5, 0.1, step = 0.02)
                    ),
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'effect_forest'",
                      sliderInput(ns("tg_opt_forest_pt_size"), "Point size",
                                  min = 1, max = 10, value = 3, step = 0.5),
                      sliderInput(ns("tg_opt_forest_ann_size"), "Annotation text size",
                                  min = 8, max = 22, value = 12, step = 1)
                    )
                  ),
                  ctrl_section("export", "arrows-alt", "Plot Size / Export",
                    fluidRow(
                      column(6, numericInput(ns("tg_plot_width"),  "Width (px)",  800, min = 200, max = 5000, step = 50)),
                      column(6, numericInput(ns("tg_plot_height"), "Height (px)", 550, min = 200, max = 5000, step = 50))
                    ),
                    fluidRow(
                      column(6, numericInput(ns("tg_plot_dpi"), "DPI", 100, min = 72, max = 300, step = 1)),
                      column(6, tags$div(style = "margin-top: 25px;", uiOutput(ns("tg_download_plot_btn"))))
                    )
                  )
                )
              ),
              column(9,
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'growth_curves'",
                  h4("Tumor Growth Curves"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_growth"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'growth_rate'",
                  h4("Growth Rates Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_rates"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'auc'",
                  h4("Area Under the Curve Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_auc"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'waterfall'",
                  h4("Waterfall Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_waterfall"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'effect_forest'",
                  h4("Effect Size Forest Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_effect_forest"), width = "auto", height = "auto")
                  )
                ),

              )
            )
          ),
          
          tabPanel("Diagnostics",
            tab_info(
              strong("Model-fit diagnostics. "),
              "Frequentist fits: residual diagnostics. Bayesian fits: MCMC ",
              "convergence + posterior diagnostics. The view switches based ",
              "on the selected model type."
            ),
            # ── Freq diagnostics (LME4 / AUC) ──────────────────────
            conditionalPanel(
              ns = ns,
              condition = "input.tg_model_type != 'bayes'",
            fluidRow(
              column(6,
                h4("Q-Q Plot of Residuals"),
                plot_help(
                    "The residuals are standardised (centred and divided by their SD), so ",
                    "for a well-fitting model the points fall along the red dashed line. ",
                    strong("S-shape: "), "residuals have heavier tails than normal \u2014 ",
                    "consider a robust likelihood or check for outliers. ",
                    strong("Bow: "), "skewed residuals \u2014 check the volume transform ",
                    "(log vs sqrt vs none). ",
                    strong("A few outliers at the ends: "), "common; not a problem ",
                    "unless they're clinically obvious errors."),
                div(style = "margin-bottom: 24px;",
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_diag_qq"), height = "380px")
                  )
                )
              ),
              column(6,
                h4("Residuals vs Fitted"),
                plot_help(
                    "Points should scatter randomly around the red zero line. ",
                    strong("Curved blue trend: "), "the model misses non-linearity ",
                    "\u2014 growth is not log-linear over the whole study; try a different volume ",
                    "transform. ",
                    strong("Fan shape (wider on the right): "), "heteroscedasticity ",
                    "\u2014 variance grows with the mean. Switch to log transform or a ",
                    "model with variance structure."),
                div(style = "margin-bottom: 24px;",
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_diag_resid_fitted"), height = "380px")
                  )
                )
              )
            ),
            br(),
            fluidRow(
              column(6,
                h4("Scale-Location"),
                plot_help(
                    "Plots sqrt(|residuals|) vs fitted to magnify variance ",
                    "trends. ",
                    strong("Flat blue line: "), "constant variance \u2014 LMM assumption ",
                    "satisfied. ",
                    strong("Upward slope: "), "heteroscedasticity \u2014 see ",
                    "Residuals vs Fitted notes above."),
                div(style = "margin-bottom: 24px;",
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_diag_scale_location"), height = "380px")
                  )
                )
              ),
              column(6,
                h4("Random-effects Q-Q"),
                plot_help(
                    "Per-animal random intercepts (and slopes, if used) should ",
                    "be normally distributed. The values are standardised, so points ",
                    "along the red line = assumption met. A component with no spread ",
                    "(estimated variance zero) is left out. Heavy tails / clear S-shape: the LMM ",
                    "may be giving a few animals too much / too little ",
                    "shrinkage. Rare and usually safe to ignore at N \u2264 30."),
                div(style = "margin-bottom: 24px;",
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_diag_re_qq"), height = "380px")
                  )
                )
              )
            ),
            # CODE_REVIEW.md DIAGNOSTICS gap (13) — LMM influence
            # diagnostics. Only shown when diag_cooks_distance was
            # populated (LME4 path with include_diagnostics = TRUE).
            conditionalPanel(
              ns = ns,
              condition = "output.tg_show_influence === true",
              hr(),
              h4("Influence diagnostics (case-deletion)"),
              tab_info(
                strong("Cook's distance"), " measures the overall effect of ",
                "leaving out each animal, one row per animal (arm, ID, cage). ",
                "Values above the threshold (4 / number of animals) flag ",
                "animals that move the fit. Inspect them individually \u2014 ",
                "they are not automatically wrong, but they are driving the ",
                "result. ",
                strong("DFBETAS"), " breaks the influence down per fixed ",
                "effect: each column is one coefficient, each row one ",
                "animal. Values | > 1 | are conventionally flagged."
              ),
              fluidRow(
                column(6,
                  h5("Cook's distance"),
                  shinycssloaders::withSpinner(
                    DT::DTOutput(ns("tg_cooks_distance"))
                  )
                ),
                column(6,
                  h5("DFBETAS"),
                  shinycssloaders::withSpinner(
                    DT::DTOutput(ns("tg_dfbetas"))
                  )
                )
              )
            )
            ),  # /freq conditionalPanel
            # ── Bayes MCMC diagnostics ────────────────────────────────────
            conditionalPanel(
              ns = ns,
              condition = "input.tg_model_type == 'bayes'",
              tab_info(
                strong("MCMC convergence diagnostics"),
                " for all fixed-effect parameters. All chains should have ",
                "converged (R\u0302 \u2248 1.00) and sufficient effective ",
                "samples before interpreting results.",
                tags$ul(
                  tags$li(strong("Parameter"), " \u2014 model parameter name."),
                  tags$li(strong("Rhat"), " \u2014 potential scale-reduction factor; values > 1.01 indicate non-convergence."),
                  tags$li(strong("Bulk_ESS"), " \u2014 bulk effective sample size; aim for > 400 per chain."),
                  tags$li(strong("Tail_ESS"), " \u2014 tail effective sample size; important for credible interval accuracy."),
                  tags$li(strong("Converged"), " \u2014 TRUE when R\u0302 \u2264 1.01.")
                ),
                em("If Rhat > 1.01, try increasing iterations or chains before reporting results.")
              ),
              uiOutput(ns("tg_bayes_diagnostics")),
              bayes_diagnostics_help(),   # R23.6
              shinycssloaders::withSpinner(DT::DTOutput(ns("tg_mcmc_diag_table"))),
              br(),
              h5(strong("Influential observations (Pareto-k > 0.7)")),
              # R23.6: the rows now name the animal and the measurement.
              em(class = "text-muted",
                 "Measurements that move the fit too much for the leave-one-out estimate, with their animal. ",
                 "Check them against the raw data; 'What these diagnostics mean' above explains Pareto-k."),
              shinycssloaders::withSpinner(DT::DTOutput(ns("tg_bayes_influential"))),
              br(),
              # R23.5: was its own "Posterior Check" tab.
              h5(strong("Posterior predictive check")),
              plot_help(
                "The model simulates new datasets from the posterior and overlays their densities ",
                "(y_rep, 50 draws) on the observed, transformed volumes (y). ",
                strong("The observed density inside the simulated envelope: "), "the model reproduces the data. ",
                strong("Outside it: "), "consider a different transformation or random-effects specification."),
              shinycssloaders::withSpinner(plotOutput(ns("tg_pp_check_plot"), height = "450px")),
              downloadButton(ns("download_tg_pp_check"), "Download Plot"),
              br(), br(),
              fluidRow(
                column(6,
                  h5(strong("Pareto-k distribution (all observations)")),
                  plot_help(
                      "Histogram of LOO Pareto-k values across all observations. ",
                      strong("All bars left of 0.5 (orange dotted): "), "LOO is reliable. ",
                      strong("Some bars between 0.5 and 0.7: "), "moderate influence; estimates usable but inspect. ",
                      strong("Bars right of 0.7 (red dashed): "), "those observations are too influential for PSIS-LOO to estimate reliably \u2014 see the table above for which."),
                  shinycssloaders::withSpinner(plotOutput(ns("tg_pareto_k_plot"), height = "300px"))
                ),
                column(6,
                  h5(strong("Posterior predictive coverage")),
                  plot_help(
                      "Each bar = empirical coverage of the corresponding nominal posterior-predictive interval. ",
                      strong("Bars near red dashed line: "), "well-calibrated. ",
                      strong("Bars below: "), "under-coverage \u2014 model intervals are too narrow. ",
                      strong("Bars above: "), "over-coverage \u2014 too wide; model is being conservative."),
                  shinycssloaders::withSpinner(plotOutput(ns("tg_ppc_coverage_plot"), height = "300px"))
                )
              ),
              br(),
              h5(strong("MCMC trace plots")),
              plot_help(
                "Each line is one chain's draws of a parameter, after warmup. Converged chains ",
                "overlap and look like stationary noise around one level (a \u2018fuzzy ",
                "caterpillar\u2019). ", strong("A drift, a trend, or a chain at a different level: "),
                "the chains have not converged and the estimates are unreliable; run more warmup ",
                "iterations or simplify the model."),
              shinycssloaders::withSpinner(plotOutput(ns("tg_mcmc_trace_plot"), height = "500px")),
              br(),
              h5(strong("MCMC rank plots")),
              plot_help(
                  "Modern convergence diagnostic. Each chain's draws are rank-normalised; chains that have converged look like uniform draws (flat histograms) and resemble each other. ",
                  strong("Distinctly non-uniform shapes or chain-specific patterns: "), "the chain hasn't mixed and the posterior estimate is unreliable.",
                  summary = "How to read these plots"),
              shinycssloaders::withSpinner(plotOutput(ns("tg_mcmc_rank_plot"), height = "500px")),
              br(),
              fluidRow(
                column(6,
                  h5(strong("Treatment-parameter posterior densities")),
                  plot_help(
                      "Posterior density (filled area) for each treatment coefficient. ",
                      "Look for: (a) modes shifted from the prior (the data informed the parameter), (b) credible-interval bars that don't include zero (significant posterior direction), (c) heavy tails / multiple modes (potential identifiability issue)."),
                  shinycssloaders::withSpinner(plotOutput(ns("tg_posterior_dist_plot"), height = "400px"))
                ),
                column(6,
                  h5(strong("Bayesian residuals vs day")),
                  plot_help(
                      "Posterior-mean residuals plotted against study day, faceted by treatment, with a loess smoother. Flat horizontal line at zero per group = the linear-time assumption is satisfied. ",
                      strong("Systematic curvature: "), "the log-linear time assumption is violated; read growth-rate comparisons with care."),
                  shinycssloaders::withSpinner(plotOutput(ns("tg_bayes_residuals_plot"), height = "400px"))
                )
              )
            )  # /bayes conditionalPanel
          ),

          # R20-D.X: replaced "Download All"; the user chooses what to export.
          tabPanel(
            title = tags$span(icon("file-zipper"), " Download Results"),
            value = "__download_results__",
            wellPanel(download_results_ui(ns, "tg"))
          )
        )
      )
    )
  )
}
