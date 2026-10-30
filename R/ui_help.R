# The Help and About pages (split from ui_modules.R, CODE_REVIEW.md R22-K).

#' Help/About UI
#' @keywords internal
helpUI <- function() {
  fluidPage(
    h2("Help & Documentation"),

    div(
      class = "alert alert-warning",
      style = "border-left: 5px solid #e2a000;",
      h4(icon("triangle-exclamation"), " Disclaimer"),
      p(
        strong("mouseExperiment"), " is provided \"as is\" for research purposes only,",
        " without warranty of any kind, express or implied. The authors and contributors",
        " accept no liability for errors, inaccuracies, or any consequences \u2014 including",
        " but not limited to incorrect results, data loss, or decisions made on the basis",
        " of outputs from this software.",
        " All results must be independently validated before use in publications, reports,",
        " regulatory submissions, clinical decisions, or any other application."
      )
    ),

    hr(),

    # \u2500\u2500 Workflow overview \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    # Quick start (26.10.009): the shortest path to a first result.
    div(
      class = "alert alert-info",
      h4(icon("rocket"), " Quick start"),
      tags$ol(style = "margin-bottom:0;",
        tags$li("On the ", strong("Data"), " tab, under ", strong("Try Demo Data"), ", choose ",
                strong("All Features"), " and click ", strong("Load Demo Data"), "."),
        tags$li("Below the column mapping, check ", strong("Study Design"), ": volume units mm\u00b3 and the ",
                "control group Vehicle."),
        tags$li("Open ", strong("Tumor Growth"), " and click ", strong("Run Analysis"),
                "; the results open on the Summary tab."),
        tags$li("Look through the other result tabs (Plots, Comparisons, Diagnostics), then save what you ",
                "need from ", strong("Download Results"), ". Every other analysis tab works the same way.")
      )
    ),

    h3(icon("arrows-turn-right"), " Basic Workflow"),
    tags$ol(
      tags$li(strong("Data tab"), " \u2014 upload a CSV file and map your columns (ID, time, volume, treatment). Optional columns: cage, dose, body weight, survival/censor indicator. The Data check below the summary then lists likely errors in the file (duplicates, negative volumes, large jumps, missed days) with their line numbers."),
      tags$li(strong("Study Design"), " (on the Data tab) \u2014 declare the tumour volume units (mm\u00b3 or cm\u00b3), the control group, and, when a cage column is mapped, whether cage labels are unique across the study or numbered within each arm (C1 and C2 in every arm are then different cages) and whether individual animals or whole cages were randomised. With a removal-reason column mapped, say what each reason means (death, weight loss, the volume limit, a planned end); Survival and Toxicity use it. Every analysis tab uses these settings; none has its own reference selector."),
      tags$li(strong("Analysis tabs"), " \u2014 each tab activates once the required columns are mapped:"),
      tags$ul(
        tags$li(strong("Tumor Growth"), " \u2014 always available after upload"),
        tags$li(strong("Survival"), " \u2014 requires a survival/censor column (or use the volume-threshold derivation)"),
        tags$li(strong("Toxicity"), " \u2014 requires a body weight column"),
        tags$li(strong("Dose Response"), " \u2014 requires a dose column"),
        tags$li(strong("Drug Synergy"), " \u2014 requires control, single-agent A, single-agent B, and combination groups"),
        tags$li(strong("Power Analysis"), " \u2014 always available; no data upload needed for a priori modes")
      ),
      tags$li("Select the statistical ", strong("method"), " in the sidebar, configure parameters, and click ", strong("Run Analysis"), ". The ", icon("circle-question"), " beside an option explains it: hover over it, or click or tap it."),
      tags$li("Explore results across sub-tabs (summary, plots, pairwise comparisons, model output). Each Diagnostics tab has a ", strong("How to read this plot"), " note under every plot."),
      tags$li("Download one table or plot from the button beside it, or choose any mix of tables, plots, diagnostics and data from a tab's ", strong("Download Results"), " tab: one ZIP, with plots in PNG, PDF or SVG and a README listing every file, the settings of the run and the package versions. Every tab also offers an ", strong("R script"), " that repeats the run with the mouseExperiment package, from the analysed data that comes with it; tick ", strong("Include my original file"), " and it starts from the file you uploaded, replaying your exclusions, the study days, the volumes and the arm labels. Tumor Growth also offers an HTML report."),
      tags$li("Long analyses (Bayesian fits, randomisation tests, the power simulation, the body-weight model) run on a background worker, so the app stays usable. You can close the browser and come back: a banner at the top lets you resume or cancel."),
      tags$li("To reuse the column mapping and settings with another file of the same layout, save them with ", strong("Save Preset"), " at the bottom of the Data tab and load them next time.")
    ),

    hr(),

    # \u2500\u2500 Data format \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    h3(icon("table"), " Data Format"),
    p("Upload data in ", strong("long format"), ": one row per measurement per animal, as a CSV file (.csv) of up to 30 MB. Commas, semicolons or tabs can separate the columns (with a semicolon or tab, a decimal comma is read correctly), and the file can be UTF-8, Latin-1 or UTF-16."),

    h4("Demo data"),
    p("No data of your own? The Data tab ships four demo datasets:"),
    tags$ul(
      tags$li(tags$strong("Combo Treatment"), " \u2014 four arms with survival data, for the Drug Synergy tab."),
      tags$li(tags$strong("Dose Levels"), " \u2014 graded doses with survival, for Dose Response."),
      tags$li(tags$strong("Weight"), " \u2014 four arms with body weight, for Toxicity."),
      tags$li(tags$strong("All Features"), " \u2014 the most complete, and the one that reaches every module: ",
              "six arms (vehicle, a three-dose series, a second drug and a combination), length/width/height, ",
              "body weight, dose, two cages per arm, necrosis flags, survival censoring and a removal-reason ",
              "column (tumour volume limit, weight loss, found dead, end of study).")
    ),

    h4("Column auto-detection"),
    p("The dashboard detects columns by name. Recognised patterns (case-insensitive):"),
    tags$table(
      class = "table table-condensed table-bordered",
      style = "max-width:680px;",
      tags$thead(
        tags$tr(tags$th("Role"), tags$th("Recognised names"), tags$th("Required?"))
      ),
      tags$tbody(
        tags$tr(tags$td("Animal ID"),         tags$td("ID, Mouse_ID, Animal_ID, Subject, Tag, Ear Tag, Mouse Tag"), tags$td("Yes")),
        tags$tr(tags$td("Time / Day"),        tags$td("Day, Time, Timepoint, Day_Number"),  tags$td("Yes \u2021")),
        tags$tr(tags$td("Date"),              tags$td("Date, Measurement_Date"),            tags$td("Yes \u2021")),
        tags$tr(tags$td("Tumor volume"),      tags$td("Volume, Tumor_Volume, TV"),          tags$td("Yes \u2020")),
        tags$tr(tags$td("Tumor length"),      tags$td("Length, Tumor_Length"),              tags$td("Yes \u2020")),
        tags$tr(tags$td("Tumor width"),       tags$td("Width, Tumor_Width"),                tags$td("Yes \u2020")),
        tags$tr(tags$td("Tumor height"),      tags$td("mapped manually in Data Upload"),    tags$td("No")),
        tags$tr(tags$td("Treatment / Group"), tags$td("Treatment, Group, Treatment_Group"), tags$td("Yes")),
        tags$tr(tags$td("Cage"),              tags$td("Cage, Cage_ID"),                     tags$td("No")),
        tags$tr(tags$td("Dose"),              tags$td("Dose, Dose_Level, Concentration, Conc"), tags$td("No")),
        tags$tr(tags$td("Body weight"),       tags$td("Weight, Mass, Body_Weight, Body Mass, BW"), tags$td("No")),
        tags$tr(tags$td("Sex"),               tags$td("Sex, Gender, M/F"),                  tags$td("No")),
        tags$tr(tags$td("Necrotic flag"),     tags$td("Necrotic, Necrosis, Necr, Is_Necrotic, Necrotic_Flag"), tags$td("No")),
        tags$tr(tags$td("Removal reason"),    tags$td("Removal_Reason, Reason, Cause_Of_Death, Fate"), tags$td("No")),
        tags$tr(tags$td("Survival / Censor"), tags$td("Survival_Censor, Censor, Event, Survival_Event"), tags$td("No"))
      )
    ),
    p(tags$em("\u2020 Provide either a pre-calculated Volume column, or Length + Width columns."),
      " When volume is calculated from calipers, six formulas are selectable in Data Upload: ",
      "ellipsoid (default, ", tags$code("\u03c0\u00b7L\u00b7W\u00b2/6"), "), modified ellipsoid (",
      tags$code("L\u00b7W\u00b2/2"), ") \u2014 the common preclinical convention \u2014 plus 3-axis ellipsoid, cylinder, ",
      "sphere and box. The 3-axis and box formulas use a Height column if one is mapped; the others ",
      "approximate height by width. Volumes may be supplied in mm\u00b3 or cm\u00b3."),
    p(tags$em("\u2021 Supply either a numeric Day column or a Date column; elapsed days are derived from the date and a reference start date.")),
    p("Survival/Censor must be numeric: ", strong("1 = endpoint reached"), " (e.g. tumour exceeded IACUC volume limit), ",
      strong("0 = censored"), " (alive or study-end censored). If absent, use the volume-threshold derivation in the Survival tab sidebar. An animal's time to event is its last measurement day, or with the volume threshold, the first day its tumour reaches the limit; no separate time column is needed."),

    hr(),

    # \u2500\u2500 Analysis modules \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    h3(icon("book"), " Key Concepts"),

    h4("Evaluable days"),
    p(evaluable_rule_text()),
    p("Comparisons on one day, or over a window, use only days on which every arm compared is ",
      "evaluable. The AUC's window ends at the last such day. Dose Response and the therapeutic ",
      "window use the last evaluable day unless you choose another, and Drug Synergy offers every ",
      "evaluable day or one of them. Later than that, an arm is represented by too few animals, ",
      "and by the ones whose tumours grew slowest, for a fair comparison."),

    h4("Glossary"),
    tags$dl(
      tags$dt("Reference group"),
      tags$dd("The control arm chosen under Study Design. Comparisons are made against it, or between ",
              "all pairs of arms where a tab offers that."),
      tags$dt("Geometric mean"),
      tags$dd("The average taken on the log scale and converted back: the typical tumour volume. A few ",
              "very large tumours move it far less than they move the ordinary (arithmetic) mean. ",
              "Model-based estimates of tumour volume are geometric means."),
      tags$dt("TGI (tumour growth inhibition)"),
      tags$dd("(1 \u2212 treated volume / control volume) \u00d7 100 at the evaluation day, from each arm's ",
              "mean volume: the geometric mean under the model-based estimand, the arithmetic mean under ",
              "the observed ones. 100 % means no tumour left, 0 % the same as control, and a negative ",
              "value a larger tumour than control."),
      tags$dt("Estimand"),
      tags$dd("What a number is meant to estimate, decided before looking at the data. For example, ",
              "each arm's mean volume on day 21 had every animal stayed on study (model-based), or the ",
              "mean of the animals still present on day 21 (survivors). They differ when animals leave ",
              "early, and the Toxicity and Drug Synergy sidebars let you choose."),
      tags$dt("Random effects (random slopes)"),
      tags$dd("Each animal's own starting size and growth rate, around its arm's. They let a model use ",
              "repeated measurements of the same animal correctly, instead of treating them as ",
              "independent animals."),
      tags$dt("Adjusted and unadjusted p-values"),
      tags$dd("Comparing several arms with the control makes several tests, and some p-values will be ",
              "small by chance. An adjusted p-value (Dunnett, Tukey, Holm, Bonferroni or FDR) accounts ",
              "for the whole family of comparisons shown, and adjusted intervals do the same for ",
              "intervals. Report adjusted values when several arms are compared; an unadjusted value ",
              "suits only a single planned comparison."),
      tags$dt("Confidence and credible intervals"),
      tags$dd("A 95 % confidence interval (frequentist) comes from a procedure that covers the true ",
              "value in 95 % of repeated studies. A 95 % credible interval (Bayesian) contains the ",
              "value with 95 % probability, given the data, the model and its priors."),
      tags$dt("Hazard ratio and time ratio"),
      tags$dd("A hazard ratio below 1 means the endpoint is reached at a lower rate than in the ",
              "reference arm. A time ratio above 1 (Bayesian survival) means it takes longer to reach."),
      tags$dt("Bliss difference"),
      tags$dd("The combination's observed effect minus the effect expected if the two drugs acted ",
              "independently, in percentage points. Positive suggests synergy, negative antagonism; ",
              "its 95 % interval decides the verdict.")
    ),

    hr(),

    h3(icon("flask"), " Analysis Modules"),

    h4("Tumor Growth"),
    tags$dl(
      tags$dt("LME4 (linear mixed-effects)"),
      tags$dd("Log-transforms tumor volume and fits a linear mixed model with a Treatment \u00d7 Day interaction. Each animal has its own intercept and growth rate (random slopes), falling back to simpler structures only if the fit fails; the ANOVA reports F-tests with Satterthwaite degrees of freedom. Cage handling is chosen to match the housing design rather than applied uniformly \u2014 see the tab's Info panel. Pairwise contrasts via emmeans; the comparison family and its correction are two separate settings, defaulting to vs-reference with an exact Dunnett correction."),
      tags$dt("AUC (area under the curve)"),
      tags$dd("The area under each arm's fitted growth curve over the window in which every arm is evaluable (at least 50 % of its enrolled animals, and at least 3, on study). Animals that die or are removed early count through all their measurements, so early deaths of large tumours do not pull the AUC down. Arms are compared by the ratio of their AUCs, with intervals and Wald p-values from the model and the selected correction. Each animal's own trapezoidal AUC is shown for description only: an animal removed early accumulates less area from being measured for less time."),
      tags$dt("Bayesian LMM (brms)"),
      tags$dd("The LME4 model fitted via MCMC (Stan). Returns posterior distributions for treatment effects with credible intervals, plus convergence diagnostics (Rhat, ESS, trace and rank plots) and PSIS-LOO influence measures. Requires the brms package."),
      tags$dt("Randomisation test (optional, any frequentist mode)"),
      tags$dd("Permutes treatment labels and refits, giving an exact p-value for the Treatment \u00d7 Day interaction that assumes neither a denominator-df approximation nor normal random effects. The unit permuted, individual animals or whole cages, is the unit of randomisation declared under Study Design on the Data tab \u2014 the test is only valid if that matches how the experiment was actually run."),
      tags$dt("Growth Delay"),
      tags$dd("A per-animal endpoint from the data the run analysed, without the model: T \u2212 C, ",
              "the median time for an arm's tumours to reach a volume you set, minus the reference ",
              "arm's, with a bootstrap interval, a log-rank test against the reference, and the ",
              "specific growth delay (T \u2212 C in the reference arm's doubling times). A table ",
              "gives each animal's time.")
    ),

    h4("Survival Analysis"),
    tags$dl(
      tags$dt("Kaplan-Meier"),
      tags$dd("Non-parametric survival estimates with 95% confidence bands. Log-rank test for group differences."),
      tags$dt("Cox proportional hazards"),
      tags$dd("Semi-parametric Cox PH model. Hazard ratios with 95% CI for each treatment vs. reference. When whole cages were randomised (Study Design), each arm is instead compared with the reference by a cage-level permutation log-rank, and the hazard ratios have no interval. Firth's bias-reduced estimation is applied automatically when any arm has zero events, where the standard partial likelihood would return infinite coefficients. The proportional-hazards assumption is tested and reported \u2014 if it fails, the hazard ratio is an average over a period in which the true ratio changed, and should not be quoted as a single number."),
      tags$dt("Accelerated failure time (AFT)"),
      tags$dd("A parametric model of the time to the endpoint (Weibull, log-normal or log-logistic, the best by AIC unless you choose). Each arm is compared with the reference as a time ratio: 1.8 means its animals take 1.8 times as long to reach the endpoint. It does not assume proportional hazards, so use it when the Cox model's check rejects them. Diagnostics compares the distributions and draws the fit against Kaplan-Meier."),
      tags$dt("Bayesian AFT (brms)"),
      tags$dd("Accelerated failure time model with a Weibull or log-normal distribution. Returns time ratios and derived hazard ratios with credible intervals.")
    ),

    h4("Toxicity (Body Weight)"),
    tags$dl(
      tags$dt("LME4 body weight model"),
      tags$dd("Linear mixed model of body weight over time, with each animal's own intercept and slope. Each arm's rate of weight change (per day) is compared with the reference arm's; weight at the mean study day is reported as well. Tumour mass can be subtracted from body weight first (on by default)."),
      tags$dt("Therapeutic window"),
      tags$dd("Sets each arm's tumor growth inhibition beside its mean worst weight loss (each animal's largest drop below its own first weighing), each with a 95% interval, and flags the weight loss as Tolerated, Not tolerated or Unclear against a threshold you set (default 20%). The two axes are read separately; a single efficacy-toxicity ratio was reported before version 26.09.005 and was withdrawn."),
      tags$dt("Weight-loss time-to-event"),
      tags$dd("Kaplan-Meier for time to crossing a weight-loss threshold. With a removal-reason column mapped on Data Upload, removals for weight loss count as events, planned ends as censoring, and other removals (such as tumour burden) as competing events, with an Aalen-Johansen cumulative incidence estimate. Without one, an animal whose record ends early is censored, and the tab says so. A per-arm breakdown of how animals left is reported."),
      tags$dt("Endpoint estimand"),
      tags$dd("Shared by the therapeutic window and the Drug Synergy tab. ", strong("Model-based"), " (recommended) estimates each arm's geometric mean at the endpoint day from a mixed model of every observation, so animals removed earlier still count. ", strong("Each animal's last observation"), " uses every animal's final measurement, whatever day it was taken. ", strong("Survivors"), " uses only animals still present at the endpoint day; because the fastest-growing tumours leave first, it flatters the control arm and understates TGI, and is kept only to reproduce older figures.")
    ),

    h4("Drug Synergy"),
    tags$dl(
      tags$dt("Bliss Independence"),
      tags$dd("The expected combined effect if the two drugs act independently: with FE = TGI / 100, Expected FE = FE_A + FE_B \u2212 FE_A \u00d7 FE_B. The observed combination is compared with it at one day or across every evaluable day, and the difference is reported in percentage points. Positive = synergy; negative = antagonism."),
      tags$dt("Why there is no Combination Index"),
      tags$dd("Loewe additivity, which a Combination Index measures, is ",
              "dose-equivalence: it needs a dose-response curve for each agent, ",
              "which single-dose designs do not provide. The usual single-dose ",
              "stand-in is ", tags$em("response additivity"),
              ", which fails the sham-combination test: combining a drug with ",
              "itself it predicts twice the effect, so an agent at 50 % inhibition ",
              "should reach 100 %, and when it does not the method calls the drug ",
              "antagonistic with itself. Measured across the effect grid with the ",
              "combination set to exactly Bliss-additive, 42 % of cells were ",
              "labelled antagonistic and none synergistic."),
      tags$dt("Interval on the verdict"),
      tags$dd("The synergy verdict comes from the 95 % interval for the Bliss excess: Synergy when it lies above zero, Antagonism when it lies below, and Additive (no departure from Bliss detected) when it contains zero. Without intervals the label says so.")
    ),

    h4("Dose-Response"),
    tags$dl(
      tags$dt("One agent's dose series"),
      tags$dd("A curve describes one agent at several doses against the study's control. Choose the agent's arms in the sidebar; other agents, combinations, and arms sharing a dose cannot join the same curve."),
      tags$dt("Hill / Emax"),
      tags$dd("Fits a four-parameter log-logistic (4PL) dose-response curve via nonlinear least squares; with six or more doses a five-parameter (asymmetric) curve is also tried and kept when AIC prefers it. The lower asymptote cannot fall below zero. Reports the EC50 (the dose giving half the fitted response range) with a 95% interval on log dose, the Hill slope, and the fitted curve; an EC50 outside the tested doses is flagged."),
      tags$dt("Trend tests"),
      tags$dd("Alongside the curve fit, the tab reports linear and quadratic trend p-values and a two-sided Jonckheere-Terpstra test. The last is the one that respects dose ordering: it asks whether response is monotone in dose, which is the question a dose-response study is usually posing, rather than merely whether the groups differ.")
    ),

    h4("Power Analysis"),
    tags$dl(
      tags$dt("Analytic (a priori)"),
      tags$dd("Classical power for each treated-vs-control comparison (two-sample t-test), with the alpha split across the comparisons a multi-arm study makes (Dunnett or Bonferroni) and an optional attrition allowance, across a range of effect sizes and sample sizes."),
      tags$dt("Survival endpoint (a priori)"),
      tags$dd("Under the analytic mode, choose \"Hazard ratio\" to plan the time to an event, such as the tumour volume limit, compared by the log-rank test. Schoenfeld's approximation gives the events one comparison needs from the hazard ratio, and the animals per arm from the share of controls expected to have the event by the end. It is asymptotic, so with few events treat the N as a minimum."),
      tags$dt("LMM simulation"),
      tags$dd("Simulates repeated-measures data under a linear mixed model and estimates power empirically across a grid of sample sizes. An optional euthanasia volume limit models attrition; leaving it blank assumes every animal is measured at every timepoint, which overstates power by roughly 6 points at 30 % attrition. Note that baseline SD and random-intercept SD do not affect the answer \u2014 the test is on growth rates, which are orthogonal to per-animal baselines. Random-slope SD is the parameter that matters.")
    ),

    # R23.15
    h4("Results"),
    tags$dl(
      tags$dt("Every run, kept"),
      tags$dd("Each analysis you run is listed on the Results tab, newest first, until the session ends. Open one to show it again in its tab, with its tables, plots and Download Results (the sidebar keeps its current settings). Select several to download them together: one folder each, with what the tab's Download Results offers, the plots drawn with the tab's current settings."),
      tags$dt("No refit for the same fit"),
      tags$dd("A Bayesian fit with the same data and settings as one already run in the session is shown again rather than refitted; the seed is fixed, so a refit would give the same results.")
    ),

    hr(),

    # \u2500\u2500 Troubleshooting \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    h3(icon("wrench"), " Troubleshooting"),
    tags$dl(
      tags$dt("File won't upload"),
      tags$dd("Save it as CSV (.csv), up to 30 MB. If a row is reported with the wrong number of fields, a value probably contains the separator, such as a comma inside a treatment name: quote that value or remove the comma."),
      tags$dt("Columns not auto-detected"),
      tags$dd("Use the column mapping dropdowns in the Data tab to assign columns manually."),
      tags$dt("Analysis tab is greyed out"),
      tags$dd("Upload and configure data first. Survival requires a censor column (or enable volume-threshold derivation). Toxicity requires a body weight column. Dose Response requires a dose column. Drug Synergy requires at least four treatment groups."),
      tags$dt("Analysis fails or produces no output"),
      tags$dd("Verify that the selected columns contain valid, non-missing data. Tumor volumes must be positive. Time points must be numeric. At least two groups must be present for any comparison."),
      tags$dt("Bayesian model is slow to start"),
      tags$dd("The first fit of a model compiles it, which takes about two minutes on mouseexperiment.com and less on a fast computer; later fits with the same structure reuse it and usually take under a minute. On mouseexperiment.com compiled models are kept across restarts and updates of the server, so each structure is compiled once; on your own computer, see ?clear_compiled_model_cache to keep them between R sessions. Fits run one at a time on a background worker, so other users and tabs stay responsive, and yours may wait for one ahead of it."),
      tags$dt("Plot download produces a blank file"),
      tags$dd("Run the analysis first, then click the download button. Download Results draws every plot of the run with the current settings (colours, titles, sizes), including plot types you have not opened.")
    )
  )
}

#' About UI
#' @keywords internal
citationUI <- function() {
  pkg_version   <- tryCatch(as.character(utils::packageVersion("mouseExperiment")),      error = function(e) "(not installed)")
  # Neither package has a Date field any more (it went stale; R20.D26), and
  # the fallback here was a fixed 2026-04-03. The build date is always there.
  built_date <- function(pkg) tryCatch({
    b <- utils::packageDescription(pkg)$Built
    d <- if (is.null(b)) NA_character_ else regmatches(b, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", b))
    if (length(d) && !is.na(d)) d else "unknown"
  }, error = function(e) "unknown")
  pkg_date  <- built_date("mouseExperiment")

  fluidPage(
    h2("About"),

    fluidRow(
      column(12,
        wellPanel(
          h4(icon("box-open"), " mouseExperiment"),
          tags$table(
            class = "table table-condensed",
            style = "margin-bottom:0;",
            tags$tbody(
              tags$tr(tags$td(strong("Version")), tags$td(pkg_version)),
              tags$tr(tags$td(strong("Built")), tags$td(pkg_date)),
              tags$tr(
                tags$td(strong("Source")),
                tags$td(tags$a("github.com/sciOmics/mouseExperiment",
                               href = "https://github.com/sciOmics/mouseExperiment",
                               target = "_blank"), " (MIT licence)")
              ),
              tags$tr(
                tags$td(strong("Run it yourself")),
                tags$td(tags$code('remotes::install_github("sciOmics/mouseExperiment")'),
                        ", then ", tags$code("mouseExperiment::run_dashboard()"),
                        ". Your data then never leaves your computer.")
              )
            )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        wellPanel(
          h4(icon("envelope"), " Contact"),
          p("Report bugs, request features and ask questions on ",
            tags$a("GitHub Issues", href = "https://github.com/sciOmics/mouseExperiment/issues",
                   target = "_blank"),
            ". Issues are public, so please don't attach unpublished study data: a few ",
            "made-up rows, or a demo dataset, that show the problem are enough.")
        )
      )
    ),

    # R20.D12: say what is counted.
    fluidRow(
      column(12,
        wellPanel(
          h4(icon("chart-simple"), " Usage analytics"),
          p("mouseexperiment.com counts page views and which features are used, with ",
            "self-hosted Umami analytics: no cookies, and nothing sent to a third party. ",
            "Events carry only the analysis type, the method, counts such as the number ",
            "of rows, and yes/no flags. They never include your data, file names, column ",
            "names, treatment names or error messages. Uploaded data stays in your ",
            "session on this server.")
        )
      )
    )
  )
}
