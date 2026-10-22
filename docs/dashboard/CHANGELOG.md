# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)

## [26.10.008] - 2026-10-22

### Changed
- **The navbar puts Power Analysis after Drug Synergy,** last among the analyses: it plans the
  next study, while the others analyse the loaded one.

## [26.10.007] - 2026-10-22

### Changed
- **Tumour Growth sidebar: explanations moved into `?` tooltips.** The paragraphs under each
  option now sit in a tooltip beside its label, cut to one to three sentences. The sidebar is
  half as tall (with LME4 and the randomisation test: 2,445 px to 1,199 px), and Run Analysis
  moves up from 2,274 px to 1,028 px.
  - **Still visible:** the reference group, the unit permuted, and the note that no curve is
    extrapolated without an Event/Censor column.
  - **Removed:** an internal review reference ("CODE_REVIEW.md R3.9 / G.4") in the
    extrapolation help.
  - The full method descriptions are unchanged in the Info tab.
- **About: a contact address,** david@sciomics.bio, for questions, bug reports and feature
  requests.
- **Survival: the model menu says "Cox Proportional Hazards"** instead of "Frequentist (Cox PH)".
- **Toxicity: "Estimation method" has a tooltip** explaining when to use REML and when ML.
- **AUC: how animals that die early are handled** is explained in the sidebar tooltip, the
  Info tab, the Help glossary and the AUC results.
  - Every measurement an animal had counts, so early deaths of large tumours do not pull an
    arm's AUC down.
  - In simulations, removal at the volume limit moved the AUC by −0.2 %, where an average of
    the survivors moved by −13 %.
  - Deaths that earlier measurements only partly foreshadow leave a −3 % to −5 % underestimate.
- **Diagnostics: the "How to read this plot" notes are easy to find.** They were a collapsed
  line of small grey text below each plot, often a screen further down. Each is now a clearly
  styled link with an info icon, under the plot's title. The MCMC trace plots (Tumour Growth
  and Survival), which had none, have one.

### Fixed
- **The group colour and shape controls overflowed the plot-options panel** (Tumour Growth,
  Survival, Toxicity). At a 1400 px window the panel is 209 px wide. Each group's row needed
  about 275 px, so the shape menus ran out of the box by 66 px, and the names were cut to
  "Drug_...".
  - Each group now takes two lines: the full name, then the colour and a shape menu that
    shows the symbol only.

### Removed
- **Tumour Growth's Compare tab.** It showed the previous and current runs' treatment effects
  side by side, but only after a change of model type. So it compared LME4 growth-rate
  differences with AUC ratios, and never two runs of the same model with different settings.
  Each session also kept a full copy of the previous result just for this tab.
- **The first Tumour Growth run ignored "Exclude" for necrotic observations** (R20-D.V). The
  Necrotic Observations selector appeared only after a run, and until then the run used "No
  adjustment". So the first analysis kept the necrotic observations (42 of 604 in the All
  Features demo). The selector then appeared showing "Exclude from analysis", and the next run
  excluded them, though nothing had been changed.
  - It now appears as soon as a necrotic column is mapped, and a run uses the default it shows.
  - Each run records its necrotic handling with its settings.

## [26.10.006] - 2026-10-11

### Fixed
- **The help icons (?) showed nothing on a click or a tap** (user report). They used the
  browser's title attribute, which appears only after the pointer rests on the icon.
  - **The fix:** they are Bootstrap tooltips now, shown on hover, on keyboard focus (Tab), and
    on a click or tap. They are also wider, for their sentence or two.
  - **The menu:** each icon sits inside its input's label, so a click on it was passed on to
    the input: the Control group menu opened. Clicks on an icon now stay on the icon. This
    covers all 17 help icons (Study Design, column mapping, Survival, Dose Response).

### Tests
- A browser test clicks the Control group icon with real mouse events. It checks that the
  tooltip shows, that no click or focus reaches the select box, and that the tooltip hides
  once the pointer and focus leave. It fails on 26.10.005.
- Tests that read the sources skip when run from an installed copy, as the VPS staging gate
  now runs them, instead of failing (`CODE_REVIEW.md` R20-D.T).

## [26.10.005] - 2026-10-10

Requires mouseExperiment ≥ 0.33.0. The Round 21 Minors, the grey synergy
bars and the slow plot labels (`CODE_REVIEW.md` R20-D.S).

### Fixed

- **Uploads (R20.D23).** These now load:
  - semicolon-separated files with decimal commas, as Excel writes in much of Europe;
  - tab-separated files, which loaded silently as a single column;
  - Latin-1 and UTF-16 files.

  A row that does not match the header is reported with its line number,
  where a stray comma used to shift the columns silently. Files up to 30 MB
  are accepted.
- **Typing a plot title froze the app (R20.D53).** Every keystroke rebuilt
  the growth plot. Labels are now applied on top of it after a short pause.
- **Synergy bars were all grey (R20.D20).**
- **The growth plot dropped the necrotic markers,** and cut off pre-treatment
  days (R20.D55). It now explains when extrapolation can't apply.
- **Labels and help that misstated results.**
  - Survival and toxicity KM plots label their log-rank p as unadjusted (R20.D29).
  - The survival help no longer calls an interval that excludes 1 significant, and no longer promises features that do not exist (R20.D54).
  - The cage and random-effects help match what the model does (R20.D30).
  - The Power summary uses the target power you chose (R20.D60).
  - A dropped covariate is explained (R20.D58).
  - The growth-rate mean bar no longer shows negative rates as zero, and an all-pairs forest plot has no reference row (R20.D57).
- **Bayesian results (R20.D59, R20.D32).**
  - Missing sampler diagnostics show as unavailable, not as clean.
  - The posterior summary is shown.
  - Bayesian runs get a forest plot.
  - The sampler advice names settings you can change.
- **Synergy, dose-response and the therapeutic window** pass the declared
  volume unit, for the backend's zero floor (0.33.0).

### Removed

- Dead code (R20.D24, R20.D25, R20.D61):
  - two broken entry-point files;
  - an unread log buffer;
  - nine unused functions;
  - two plots built on every run and never shown;
  - a variance panel that never rendered.

## [26.10.004] - 2026-10-09

Long analyses no longer freeze every user (`CODE_REVIEW.md` R20-D.R).

### Fixed

- **The power simulation, the randomisation test and the body-weight model
  ran in the Shiny process (R20.D6),** and every session on the server froze
  while they did: two minutes for the default simulation, which anyone could
  start without uploading data. They now run on the background worker, like
  the Bayesian fits, with a progress window, Cancel, and Resume from the
  banner. The tumour-growth fit appears at once and the randomisation test
  fills its panel when it finishes. Backend warnings are still shown.

## [26.10.003] - 2026-10-09

Requires mouseExperiment ≥ 0.31.0. Plots and tables that mislabelled their
results, and the rest of the analytics leak (`CODE_REVIEW.md` R20-D.Q).

### Fixed

- **Q-Q plots always looked non-normal (R20.D42).** The Tumor Growth
  residual Q-Q drew y = x over raw residuals. It now shows the backend's
  standardised plot (backend R20.84), as the other Q-Q plots already did.
- **Forest plots paired unadjusted intervals with adjusted p-values
  (R20.D43).** They now use the backend's adjusted intervals (R20.85). A
  caption says how the intervals are adjusted.
- **The Comparisons tab (R20.D44)** showed the dropdown's family rather than
  the run's, described Dunnett or Tukey whatever correction was used, and
  never showed its title. The table now follows the run, says when the
  dropdown has changed since, and has its title as a caption. The help
  describes the correction actually used.
- **The Data Summary showed log-scale means labelled mm³ (R20.D45).** It is
  now on the volume scale (backend R20.57), with the unit in the headers and
  a geometric mean for log-scale runs.
- **The HTML report recorded the sidebar, not the run (R20.D46).** Its
  settings now come from the run.
- **Power's "Seed SD from loaded data" pooled every day (R20.D51),** which
  inflated the required N up to 5-fold on the demo (68 instead of 14 per
  group). It now takes each animal's last volume up to one endpoint day,
  which you can set.
- **Analytics (R20.D12).** No event carries a column name, a treatment or
  reference name, or an error message. Staging is not counted under the
  production site. The About page says what is counted.

## [26.10.002] - 2026-10-09

Requires mouseExperiment ≥ 0.30.0. The dashboard side of backend 0.30.0
(package `CODE_REVIEW.md` R20-T), which changed what four results hold.

### Changed

- **Tumour growth passes the declared volume unit (backend R20.32).** On the
  log scale a zero volume now takes a fixed floor, 2 mm³ or 0.002 cm³, instead
  of half the smallest volume in the data. The Summary states the floor and
  how many values it applied to. Undeclared units are treated as mm³.
- **The ANOVA tab describes the new table (backend R20.31).** The table leads
  with the test of whether growth rates differ, then the arms at the mean
  study day, then average growth; a `Tests` column states each question. The
  help no longer describes a Type III table, and the LME4 Info tab no longer
  describes Wald chi-square tests, which the backend stopped using.
- **Body weight shows rates of weight change (backend R20.44).** The table
  that was "Estimated Marginal Means" is each arm's rate (g/day), and a new
  table, with its own download, compares each arm's rate with the
  reference's. The dashboard never showed the comparisons before.
- **Power simulation (backend R20.33).** The Summary says that power is for one
  treated arm against the control at α/(k−1), and names the new
  `Power_Any`, `Power_All` and `Power_Omnibus` columns. The curve's band is
  the Wilson interval, and the method text no longer describes a
  likelihood-ratio test.

## [26.10.001] - 2026-10-09

Requires mouseExperiment ≥ 0.28.0.

### Changed

- **Demo datasets cleaned up.** "Try Demo Data" now offers four datasets with
  short names, in this order: Combo Treatment, Dose Levels, Weight, All Features.
  - "Simple" (random values) and "Necrotic Data" are removed, with their files.
    The necrosis feature stays; All Features carries a necrosis flag.
  - "Master" is renamed **All Features** (`all_features_synthetic_data.csv`). It
    gains a `Removal_Reason` column, so the reason-dependent options have data:
    tumour volume limit (19 animals), end of study (26), weight loss >20% (two
    Drug_A High animals, whose weights now fall past 20 % before removal) and
    found dead (one combination animal, day 14). Both new removals are censored
    for the tumour endpoint, so Survival's "count deaths as events" and the
    Toxicity tab's weight-loss separation each change a result.

## [26.09.006] - 2026-10-01

Requires mouseExperiment ≥ 0.28.0 and future ≥ 1.49.0. Sixth implementation
step of the Round 21 review (`CODE_REVIEW.md` R20-D.N): Bayesian robustness.

### Fixed

- **Bayesian fits ran in the Shiny process and froze every session
  (R20.D62).** `plan(multisession, workers = 1L)` is sequential in future; the
  dashboard now uses one real background worker
  (`plan(multisession, workers = I(1))`). An LME4 fit in another session now
  finishes in seconds while a Bayesian fit runs.
- **Cancel did not stop a fit, and repeated clicks queued several (R20.D7).**
  The job registry now queues fits itself, runs one per browser at a time,
  and interrupts a running fit on Cancel, which is also in the progress
  window. The window shows the queue position or the elapsed time.
- **Resume lost the result it resumed (R20.D63),** because the job was retired
  before its poller could deliver it.
- **Resume from an earlier browser session did nothing (R20.D64):** the newest
  session owned every resume handler. Handlers are now per session.
- **"Dismiss all" stranded a waiting tab's modal and dropped its result
  (R20.D5).** It is now "Dismiss finished", and a tab whose job disappears
  closes its modal with a message.
- **Each session kept its Bayesian model alive, about 31 MB (R20.D18).** The
  diagnostics now come from the result's plots and posterior draws (backend
  R20.80), and the model is not returned.
- The influence diagnostics use the backend's Pareto-k threshold, which
  depends on the number of draws (backend R20.35).
- Bayesian survival from volume-derived events works (R20.D8, fixed in the
  backend's column handling).
- The Bayesian survival launch no longer sends the reference arm's name to
  analytics (R20.D12).

### Changed

- A second fit of the same model structure reuses the compiled model: 11 s
  against 40 s for the first fit on the Master demo (backend R20.78). The
  timing text in the Tumour Growth and Survival tabs says so.

## [26.09.005] - 2026-09-30

Requires mouseExperiment ≥ 0.27.0. Fifth implementation step of the Round 21
review (`CODE_REVIEW.md` R20-D.M): the two-axis therapeutic window, the
interval-based synergy verdict, one agent's dose series, and the optional
removal-reason column.

### Added

- **An optional removal-reason column (R20-D.H).** Map it under Optional
  Columns on Data Upload (it is auto-detected from names such as
  Removal_Reason, Reason or Cause_Of_Death). The Study Design panel then lists
  each of its values with a category: tumour volume limit, death or moribund,
  weight loss or body condition, planned end or data cut, or other. The
  categories are guessed from the wording and can be changed, and the preset
  saves them.
- **Survival counts deaths from other causes as events (R20.D13).** A
  checkbox, "Count deaths from other causes as events (recommended)",
  appears when a reason column is mapped. With it ticked, an animal whose
  reason is a death is an event at its last measurement, with either event
  source. The sidebar states the assumption in force, and the summary counts
  the deaths added.
- **The Dose Response tab shows the Jonckheere–Terpstra test and the EC50
  interval (R20.D49),** in the Summary, Model Results and the results CSV. The
  test is two-sided, with the direction of the trend stated. The EC50 has a
  95 % interval on log dose, and an EC50 outside the tested doses is flagged.

### Changed

- **The therapeutic window is two axes (R20-K).** The Toxicity tab sets each
  arm's TGI beside its mean worst weight loss, each with a 95 % interval, and
  flags the weight loss as Tolerated, Not tolerated or Unclear against a
  threshold set in the sidebar (default 20 %). The plot shows TGI against
  worst weight loss, with both intervals and the threshold. A note calls out a
  tolerated arm that still has animals over the threshold. The TWM ratio, its
  ranked bar chart and its ranking note are gone.
- **The synergy verdict comes from the interval (backend R20.19).** It reads
  Synergy, Antagonism, or Additive (no departure from Bliss detected), and the
  rule behind it is shown under it. A verdict without an interval says so and
  is not shown in green. "Not evaluable" becomes "Bliss does not apply", since
  "evaluable" now refers to evaluable days.
- **Dose Response analyses one agent's dose series (R20.D47).** A "Dose
  series" selector lists the arms; the largest set of single-agent arms that
  share a name and differ in dose is suggested. Other agents, combinations,
  and arms that share a dose are refused with the backend's reason. When one
  treatment name is used at every dose, its dose-0 rows are the control. The
  day selector lists the evaluable days of the chosen series.
- **The weight-loss threshold uses the removal reasons (backend R20.23).**
  Weight-loss removals are events, planned ends censoring, and other removals
  competing risks. The Kaplan-Meier panel states which assumption applied;
  without a reason column an early end is censoring.

### Fixed

- **The therapeutic-window plot download produced no file (R20.D50, in
  part).** The plot is now built directly for the download. The trajectory
  and KM downloads still read an output and remain open.
- The Dose Response backend no longer prints its report to the server log.

## [26.09.004] - 2026-09-30

Requires mouseExperiment ≥ 0.26.0. Fourth implementation step of the Round 21
review (`CODE_REVIEW.md` R20-D.L): the evaluable-day rule, stated to users; one
source for TGI; the model-based AUC; random slopes by default.

### Added

- **The evaluable-day rule is stated to users (R20-K).** An arm is evaluable
  on a day when at least 50 % of its enrolled animals, and at least 3, are
  still on study.
  - The Study Design panel and the Synergy, Toxicity (therapeutic window) and
    Dose Response Info tabs state it.
  - Every result that uses an endpoint day names the day, and lists the days
    and arms the rule excluded, e.g. "Not evaluated: day 30 (Control 3/8)".
  - Day selectors list only evaluable days. Synergy's "All dates" becomes
    "All evaluable days"; Dose Response's becomes "Last evaluable day".

### Changed

- **Tumour Growth:**
  - The default random effects are random slopes. An uncorrelated-slopes
    option is added, and intercept-only is labelled with its assumption. The
    summary states the structure used and any fallback.
  - The ANOVA help describes Satterthwaite F-tests; it called car's Wald
    chi-square a likelihood-ratio test (R20.D30, part).
- **The AUC is the backend's model-based AUC:** the area under each arm's
  fitted curve over the common evaluable window.
  - The comparisons table shows AUC ratios and differences with intervals.
  - The AUC plot shows each arm's AUC with its interval.
  - The summary gives the window.
  - The per-animal trapezoid tables are relabelled as descriptive only.
- **Intervals match their estimates (R20.D35).** Synergy's intervals, and the
  therapeutic window's TGI intervals, now describe the reported estimates,
  with the method stated. The TWM ratio has no interval under the model
  estimand, and the ranking note says so and points to the TGI and weight-loss
  intervals. The sidebar control is renamed "Interval draws".
- **Influence tables have one row per animal (backend R20.76),** labelled by
  arm, ID and cage. A Tumour Growth run on the Master demo takes about 2 s,
  down from 13 s with random slopes.

### Fixed

- **The Synergy TGI plot recomputed TGI from the animals present each day
  (R20.D37),** the survivor estimand the tab warns against, and contradicted
  the table. It now plots the backend's TGIs and their 95 % intervals.
- **The Dose Response TGI plot re-derived TGI per animal against the live
  Reference dropdown (R20.D48).**
  - It now plots `tgi_table`, TGI by dose from the model at the analysis day
    (the same TGI as Synergy and Toxicity), and the table appears under
    Dose-Response Table.
  - The day is the last evaluable day, not each animal's own last observation.
- **The Synergy Info tab documented response additivity, "min(FE_A + FE_B, 1)"
  (R20.D38).** It now gives the Bliss formula and a correct ceiling example
  (50 % and 50 % give 75 %). A test guards the formula.
- **The over-time synergy plot's `ylim(0, 100)` dropped negative TGIs
  (R20.D29, part).**
- **The Synergy diagnostics show the checks for the estimate reported:** the
  endpoint model's residuals by arm and each arm's fitted curve over its data,
  or the per-animal plots for the observed estimands.

## [26.09.003] - 2026-09-30

Requires mouseExperiment ≥ 0.25.2. Third implementation step of the Round 21
review (`CODE_REVIEW.md` R20-D.K): the Study Design panel, one animal key,
declared volume units, and dates read day-first or month-first from the data.

### Added

- **Study Design panel on Data Upload (R20-D.H).** It sits below the column
  mapping, is set once, and every tab uses it:
  - **Tumour volume units (mm³ / cm³).** The demos, and uploads whose volume is
    calculated from caliper measurements, start at mm³. An upload with a volume
    column starts undeclared.
  - **Control group.** It replaces the reference selectors on Tumour Growth,
    Survival, Toxicity, Dose-Response and Synergy. Each tab now shows a
    read-only note naming the reference it uses, and says so when the control
    is excluded by the subset filter.
  - **Unit of randomisation (individual animals / whole cages).** Shown when a
    cage column is mapped. It replaces the Tumour Growth permutation-unit
    selector and drives the new survival analysis by cage.
- **Date order selector** under the date column (Detect from the data / day
  first / month first).
- **Settings preset format 2.** The preset now also saves the study design and
  the subset filter. Loading one before a dataset is uploaded says so; before,
  its settings were silently lost.

### Fixed — Critical

- **DD/MM/YYYY dates were read as MM/DD/YYYY (R20.D1).** The layout was taken
  from the first row, and any d/d/yyyy was assumed month-first. For a study
  starting 1 March, days 0–17 came out as −57, 33, 155, 247, NA, NA. The whole
  column is now read:
  - a value such as 15/03/2025 decides the order;
  - a column no value decides uses the declared order, or asks for one;
  - data that contradict the declared order are an error;
  - year-less dates take the year nearest the reference date, so a study
    that crosses New Year no longer goes negative.
  - Data Upload says how the dates were read, and warns about unreadable
    values and days before the reference date.
- **Survival merged ear tags reused across arms (R20.D2).** Animals were keyed
  by ID (plus cage). With tags restarting in each arm and no cage column, 20
  animals became 10 records. They are now keyed by arm + ID + cage. The event
  time is the first day flagged, not the last day measured; an animal whose
  every row is flagged has an animal-level status, and its last day is used.
- **Toxicity's tumour-mass correction assumed its units (R20.D34).**
  - The declared units are passed to all three backend calls. With units
    undeclared, the correction does not run, and the tab says why.
  - The trajectory plot shows the body-weight model's own net weight instead of
    re-deriving it with a hard-coded `/1000`.
  - On the Weight demo, 0 of 20 animals now cross 20 % weight loss (20 of 20
    before). Declaring cm³ for these mm³ data is refused with the backend's
    explanation.
- **Backend warnings were dropped (R20.D34).** `safe_analysis()` had only an
  error handler, so a 1000× unit error ran to "Analysis completed
  successfully!". Warnings are now collected and logged, and the status box
  reports "completed with N warning(s)" and lists up to five of them.

### Fixed

- **Survival and "All dates" synergy needed an ID column named `ID` (R20.D14).**
  - Both now pass the mapped ID and cage columns.
  - "All dates" synergy passes the estimand and bootstrap settings, which it
    used to ignore.
  - Survival no longer prints its report to the server console.
- **Survival by unit of randomisation (package R20.6).** With whole cages
  randomised, each arm is compared with the reference by a cage-level
  permutation log-rank. Model Results then shows, per arm:
  - the cages in the arm and in the reference, and the number of cage
    assignments;
  - the permutation p-value;
  - the smallest attainable p-value, with the backend's warning when it
    exceeds 0.05.
  - With individual animals and one arm per cage, the backend's caveat about
    cage-mates is shown.
- **The "Treat natural death as event" checkbox was never read (R20.D13).** It
  is removed. The sidebar states the assumption the chosen event source makes.
  A working option returns with the removal-reason column (step 5).
- **Synergy failed after a Tumour Growth run with a dose column (R20.D36).** The
  selectors listed raw treatment names while the analysis used "Treatment Dose"
  labels. Both now use the arm labels, whether or not Tumour Growth has run.
- **Therapeutic window without the cage key (R20.D31).** `cage_column` is passed,
  so cage-mates sharing an ear tag are no longer one animal.
- **The Tumour Growth drill-down plotted cage-mates that share an ear tag as one
  mouse (R20.D57, key part).**
- **The subset filter excluded an ear tag from every arm.** Animals are now
  listed and excluded by arm, ID and cage, e.g. "3 (Vehicle, cage C1)".
- **Stale preprocessing after a mapping change.** Changing the reference date,
  date order, volume formula, height, cage or dose column left the cached Day
  and Volume in place. Tabs that read the cache (Survival, Toxicity,
  Dose-Response, Synergy) then disagreed with Tumour Growth, which recomputes
  on every run.
- **Toxicity and Dose-Response read the cache directly (R20.D29, part).** So a
  subset-filter change made after the first preprocessing did not apply. They
  now use `ensure_preprocessed()`.
- **Post-run displays read live widgets.**
  - The Tumour Growth report and forest-plot labels, the Dose-Response TGI plot
    and comparison note, and the Synergy summary now use the reference and
    roles of the displayed run (R20.D46 and R20.D48, part).
  - The Tumour Growth and Survival analytics events no longer send the
    reference group's name (R20.D12, part).

### Removed

- The per-tab reference selectors (`tg_reference`, `surv_reference_group`,
  `tox_reference_group`, `dr_control_group`, `syn_control`), the Tumour Growth
  permutation-unit selector, and the unread natural-death checkbox.
- The Data tab's Reset handler. It had no button, and it replaced the
  column-mapping, summary and preview renderers with permanent NULLs.

## [26.09.002] - 2026-09-30

Requires mouseExperiment ≥ 0.24.0. Second implementation step of the Round 21
review (`CODE_REVIEW.md` R20-D.J): the one-line and near-one-line bugs.

### Fixed

- **Bayesian fits in sessions without the identity cookie (R20.D3).** They took
  a "legacy inline" path whose promise operators were never imported. The fit
  never started, and the "fit running" modal stayed open until reload.
  - That path is gone. Such sessions get a registry key that no browser can
    present (`local-` + random hex), so their jobs run through the registry and
    stay private to the session.
  - This affected any browser that blocks the cookie. Verified in headless
    Chrome with cookies disabled: the fit now runs and reports its result.
    (Browsers do set the cookie on `http://localhost`, which counts as a secure
    context, so the README's older claim that local runs lack it was wrong.)
- **Reference date (R20.D4).** A date outside the accepted formats raised
  "object 'NA_Date_' not found".
  - Tumour Growth froze on "Running analysis…", Survival showed the raw
    message, and the Data tab hid the error.
  - The intended fallback, today's date, would have counted every Day from
    today without a word, so the date is now reported instead: under the
    reference-date box on the Data tab, and in the Tumour Growth and Survival
    status messages.
- **Browser tests never ran (R20.D11).**
  - The skip call used an argument `skip_if_not_installed()` does not have.
  - Behind that, the test apps used pre-module signatures, their test hooks
    were missing from the page, and they were passed as objects that a fresh R
    process cannot load.
  - The apps are now app directories that load the package themselves. All 7
    tests run, and the upload-to-fit workflow test now requires a successful
    fit, not just a result object.
- **Power tab (R20.D17), the dashboard half of package R20.7.** The analytic
  path passed no multiplicity or attrition settings. It now has three controls,
  and the summary states what the power refers to:
  - a correction (Bonferroni by default);
  - the number of treated-vs-control comparisons (follows k − 1);
  - expected attrition, which gives Enroll_N.
  - For 4 groups at d = 1 the tab recommended 7 per group; it now recommends
    23, or 29 to enrol at 20 % attrition.
- **Survival forest-plot download (R20.D19)** returned HTTP 500: its plot builder
  was a local inside the survival server, invisible to the download handler.
  The builder now lives in `helpers_plot.R`.
- **AUC failed on first click (R20.D39).** The default correction, Dunnett, is
  undefined for AUC's separate t-tests, and so is Tukey.
  - Neither is offered for AUC any more. Holm is selected on the way in, and
    Dunnett comes back on the way out.
  - A stale Dunnett or Tukey still maps to Holm on the server.
- **Synergy (package R20.8).** When an agent did not inhibit growth, the summary
  printed "Synergistic: TRUE" and a peak day. It now says the case is not
  evaluable, and the Info tab describes the new `Evaluable` column.
- **Console noise:**
  - Info alerts logged "does not correspond to a known icon" on every use
    (R20.D21).
  - Demo loads threw 20 shinyFeedback errors (R20.D33). The cause was not the
    empty feedback text the review suggested: the feedback was sent before the
    mapping selectors existed.

- **Forest-plot axis label.** "Hazard Ratio (log₁₀ scale)" rendered as "log" plus
  two boxes, both on screen and in the download, because the PNG device font has
  no subscript digits. It now reads "log10".
- **A test-suite trap.** The TTL sweeper's test reset only a flag and left its
  callback on the global `later` loop. Shiny's `testServer()` waits for that loop
  to empty on every output read, so any later test file that read an output
  hung. The existing files happened to sort before it; the new one did not.
  `.reset_ttl_sweeper()` now cancels the callback.

### Changed

- `runDashboard()` now builds its app with `dashboard_app()`, so the full app can
  run under test. A new browser test starts it and checks that a session works.
  `testServer()` could not have caught the start-up error found during this
  step: it runs the server body inside a reactive context.
- The Bayesian survival modal no longer says "Weibull" when the family is
  log-normal.
- The README, DEPLOYMENT, ARCHITECTURE and CONTRIBUTING docs no longer describe
  the legacy inline path. The Power tab's Info panel and Help describe the
  per-comparison method.

## [26.09.001] - 2026-09-29

Requires mouseExperiment ≥ 0.23.0. First implementation step of the Round 21
review (`CODE_REVIEW.md` R20-D.H): the dashboard is narrowed to the analyses the
maintainer decided to keep.

### Removed

- **Bayesian fits other than tumor growth and survival.** Bayesian synergy,
  Bayesian dose-response, the Bayesian body-weight model and the Bayesian
  therapeutic window go, with their settings, tabs, MCMC diagnostics, downloads
  and resume handlers (`bw_resume.R`, `syn_resume.R`, `dr_resume.R`,
  `twm_resume.R`). The resume banner now knows two job kinds, `bayes_tg` and
  `bayes_surv`.
- **GAMM** from Tumor Growth: the model option, the Bayesian GAMM variant, the
  GAMM Info panel and glossary entries, and the k-check / concurvity diagnostics
  (closes R20.D40, R20.D41). The Bayesian option always fits the LMM.
- **Survival:** the exponential and gamma Bayesian AFT families. Weibull and
  log-normal remain.
- **Toxicity:**
  - Body Weight Analysis loses body-weight AUC, together with its nadir table,
    comparison table, forest plot and the post-hoc comparison selector that only
    it used.
  - The efficacy-toxicity mode loses the safety-efficacy scatter, total benefit
    score, weight-corrected TGI and the Combined Report (which always failed,
    R20.D16), along with their sidebar inputs (efficacy metric, λ, safety
    threshold).
  - What remains of that mode is the Therapeutic Window Metric, so the mode is
    now labelled **Therapeutic Window**.
- **`.rda` / `.RData` upload** (R20.D9). Upload is CSV only.
- **Session save / restore** (R20.D10). The JSON settings preset stays.
- `R/mod_efficacy_toxicity.R`, whose stubs nothing called.

### Changed

- Help and Info text no longer describe the removed analyses, and the Toxicity
  Diagnostics tab no longer switches between frequentist and Bayesian views. The
  residual-plot advice no longer suggests a GAM.
- Therapeutic-window plot groups follow the data's arm order, limited to the arms
  in the TWM table. They previously came from the bivariate result, which no
  longer exists.
- README, `docs/ARCHITECTURE.md`, `DATA_FORMAT.md`, `SECURITY.md` and
  `DEPLOYMENT.md` updated to match.

### Verification

- Suite: 184 tests, 564 expectations, 0 failures. The difference from 574 is
  the four removed resume kinds and the two GAMM labels in the docs test.
  `test-shiny_app.R` still errors at file level wherever shinytest2 is installed
  (R20.D11, unchanged).
- Every module server was driven with `testServer()` on the master demo: body
  weight, therapeutic window, tumor growth (LME4 and AUC), survival, synergy
  (single day and over time) and dose-response. Results, plots and downloads
  were checked: 37 of 37 checks passed.
- Installed with backend 0.23.0; `runDashboard()` serves HTTP 200, and the page
  has no duplicate element IDs.
- In headless Chrome, the new build and 26.07.015 were run through the same
  script:
  - Toxicity switches between modes and inserts the right tabs.
  - The TWM table and plot render.
  - The selectors offer only the kept options.
  - Both builds log the same 20 console errors, all shinyFeedback on demo load
    (R20.D33), so none are new.
- R CMD check: 6 WARNINGs, 5 NOTEs, the same as 26.07.015. The detailed output
  differs only by lines that disappeared with the removed files.

## [26.07.015] - 2026-08-05

### Fixed — stale NAMESPACE imports broke the build (production 503)

v0.21.0 removed `plot_combination_index()` and `plot_synergy_combined()` from the
backend; two `importFrom()` lines for them stayed in NAMESPACE. `R CMD INSTALL`
resolves imports at lazy-load and fails with "object 'plot_combination_index' is
not exported by 'namespace:mouseExperiment'", so the Docker build died at
`install_github` and NPM had no healthy upstream.

Two gaps let it through. The cross-repo guard added in 26.07.014 scanned
`mouseExperiment::` **call syntax** and never looked at NAMESPACE imports — a
different mechanism. And `devtools::load_all()`, which every test run in this
review used, does not resolve imports the way installation does: the suite was
green at 569 against a package that could not be installed.

Verified by the path that matters: `R CMD INSTALL` succeeds, `library()` loads,
and `runDashboard()` serves HTTP 200 on a real port — a check that had never been
run here either.

The guard now scans NAMESPACE imports too, with a companion test proving it can
fail.

## [26.07.014] - 2026-07-31

### Removed — the Combination Index from the Synergy tab (R17.2)

Backend v0.21.0 dropped it: with the combination set to exactly Bliss-additive,
42 % of the effect grid disagreed with Bliss and every disagreement was
antagonistic, because the implementation was response additivity rather than
Loewe additivity.

Gone here: the Loewe Info panel, the CI plot and its download, the CI entry in
the plot selector, the shaded-region checkbox only that plot used, and the
associated metric descriptions.

The verdict box no longer calls itself a "consensus" requiring two models to
agree. With one model there is nothing to reconcile, so it names Bliss, states
the formula, and says plainly that a threshold on three group means should be
read against the interval below it. The Help entry now explains *why* there is no
CI, sham-combination argument included — a user expecting one deserves the reason.

### Fixed — three calls to backend functions that do not exist (R17.3)

`plot_survival` (removed from the backend in `ace5b94`), `prepare_dose_data` and
`generate_summary_statistics` (both internal, never exported). All three sat
inside `tryCatch`, so every survival and dose-response run threw, logged, and
fell back — output appeared, so nothing looked wrong.

The dose-response consequence was not cosmetic: `result$analysis_data` was never
set, three renderers hit a silent `req(ad)`, and **the scatter, box and trend
plots all rendered blank**. On the survival side, the comment says "use the
exported plot_survival() function (matches demo notebook approach)" — that plot
has never once been drawn.

### Added — a guard against cross-repo staleness (R17.4)

`test-backend-calls-resolve.R` scans every `mouseExperiment::` call and fails if
one does not resolve to an exported function. It strips comments first (comments
naming a removed function are wanted; only live calls matter), asserts the
matcher can actually fail, and checks that dose-response results carry
`analysis_data` so the blank-plot failure cannot return.

Suite: 558 → 569 passing, no skips.

## [26.07.013] - 2026-07-30

### Fixed — the forest plot marked significance from unadjusted p-values on AUC (R14.6)

`p_col <- intersect(c("p.value", "p_value", "p.adj", "adj.p.value", "P_Value"), ...)`
never listed `p_adjusted`. On the AUC path `p_value` is the **raw** t-test
p-value, so the Treatment Effects forest plot marked significance without any
multiplicity correction — while the GAM path, where `p_value` holds the adjusted
value, used corrected ones. Same plot, same legend, two different multiplicity
regimes depending on a model selection invisible in the figure.

Backend v0.18.0 now supplies `P_Value_Adjusted` on every path; the sniff prefers
it and keeps the old list only as a fallback for older backends.

### Added — the consensus synergy verdict is now shown (R14.5)

`overall_assessment` requires Bliss **and** Loewe to agree before claiming
synergy. It was computed by the backend, returned, and never displayed: the tab
showed only the Combination Index verdict, which reads "antagonistic" for a
genuinely Bliss-additive combination across roughly FE 0.25–0.65 — the range most
active single agents occupy.

It now appears above the CI verdict, styled by outcome, with a note that the CI
verdict reflects Loewe alone. When the two models point opposite ways the panel
says so explicitly, rather than leaving the reader to reconcile two contradictory
labels: they encode different null hypotheses, so disagreement is information
rather than an error.

### Fixed — Bayesian intervals were mislabelled (R14.7)

The Info tab said results are "posterior medians with 95 % highest posterior
density (HPD) credible intervals". Checked against brms directly: the point
estimate is the posterior **mean** (3.89831 vs a median of 3.89553) and the
interval is **equal-tailed**, matching `quantile(draws, c(.025, .975))` exactly.

The subtler part is that both interval types ship under the same column names.
`posterior_summary` comes from `summary(brmsfit)` and is equal-tailed; the
Treatment Effects and Comparisons tables come from emmeans and are genuine
highest-density intervals — and both are relabelled `Lower_95_CrI`/`Upper_95_CrI`.
They coincide for symmetric posteriors and diverge for skewed ones, which is
exactly the case for variance components and the back-transformed time ratios in
the Bayesian survival model. The tab now says which table is which.

## [26.07.012] - 2026-07-30

### Fixed — the Tumor Growth tab claimed a metric it does not produce (R14.3)

The Info panel stated "All modes compute Tumor Growth Inhibition (TGI) for each
treatment group relative to the reference." `tumor_growth_statistics()` returns no
TGI field and the tab renders none. TGI is produced by the Toxicity path (as the
efficacy half of the therapeutic window) and computed inside Drug Synergy (as the
fractional effect the models act on) — nowhere else.

The claim predates this review and Round 11 propagated it without checking, which
is the exact failure that round existed to fix. The replacement says where TGI
actually lives and why a trajectory model does not report it: collapsing a
trajectory to a one-day ratio of group means discards most of what was measured,
and a ratio of endpoint means is precisely where early removals bias the answer.

Pinned by `test-docs-current.R` so it cannot return.

### Known — the synergy verdict shown is the one most prone to bias (R14.5)

Backend R14.5 establishes that the linear-Loewe index labels a combination
"antagonistic" across roughly FE 0.25–0.65 even when it is exactly additive under
Bliss — the range most active single agents occupy. The Synergy tab renders only
`combination_index$interpretation`. The backend also returns `overall_assessment`,
which requires both models to agree before claiming synergy, and the dashboard
never displays it.

Not changed here: which model to headline is a scientific choice for the package
owner rather than a defect to patch silently. Recorded so the decision is
explicit.

## [26.07.011] - 2026-07-30

### Removed — the dosing-schedule annotation added in 26.07.010

The toggle, the dosing-days field, the style selector and the warning banner are
gone from the Tumor Growth, Toxicity and Survival sidebars, along with the
annotation on the growth, body-weight and Kaplan-Meier plots and the Info-tab
documentation for it. The feature did not work in practice and has been withdrawn
rather than patched.

The two treatment-schedule CSVs removed in 26.07.010 stay removed — nothing in
the app ever reached them. Every file in `inst/sample_data/` remains loadable
from the demo selector.

Suite: 621 → 551 passing, no skips. The drop is the 70 assertions that covered
the withdrawn feature.

## [26.07.010] - 2026-07-30

### Added — dosing schedules drawn onto the plots (D.15)

Replaces the backend's removed `plot_treatments()`, which needed a schedule
uploaded as its own file and a second panel kept aligned with the growth curve by
hand. Its two schedule CSVs shipped here for months with nothing in the app able
to open them.

Dosing days are now entered as text in the sidebar and drawn in a strip below the
plot, behind a toggle that is off by default. Available on the three time-axis
plots: growth trajectories, body-weight trajectories (both the spaghetti and
mean ± error variants), and Kaplan-Meier curves.

- **Input.** A bare list — `1, 5, 9, 13` — marks every treated arm; one line per
  arm (`DrugA: 1,5,9,13`) when schedules differ, in which case rows stack and
  take their arm's colour. Dosing timing cannot be read from uploaded data: the
  Dose column records how much, not when, and dosing days are usually not
  measurement days. This is the one piece of study design that has to be typed.
- **Marks.** Down-triangles for an intermittent schedule, a shaded window for a
  continuous one, switching automatically above eight doses or at daily dosing,
  and overridable. They sit below the panel rather than crossing it, so they
  never obscure a curve.
- **Feedback.** A mistyped arm name produces a warning banner beside the field
  rather than a silently missing row — the §K.2 failure this review has found
  five times. Annotation failure degrades to the unannotated plot; it never costs
  the user their analysis.

The Info tab documents the feature and one caveat worth stating plainly:
separation appearing shortly after a dose mark is not evidence the dose caused
it. With four doses and eight measurement days some separation will land near
some dose by chance, and growth curves diverge steadily regardless. Onset timing
needs a model with time-varying terms, not a visual coincidence.

### Removed — the two treatment-schedule CSVs

They ship no longer, and the demo-reachability exemption list is now empty — every
file in `inst/sample_data/` is loadable from the selector.

Suite: 557 → 621 passing, no skips.

## [26.07.009] - 2026-07-30

Demo data that could not be loaded, and Info/Help tabs describing an app that no
longer exists. Both are the same failure: prose and packaged assets have no
failure mode, so they stop being true without anything noticing.

### Added — two demo datasets that shipped but were unreachable (D.12)

`master_synthetic_data` (619 × 13) and `synthetic_data` (120 × 7) were shipped by
the backend with no CSV here and no entry in the demo loader, so neither could
ever be selected. `master` carries length/width/**height**, body weight, dose,
cage, necrosis flags and survival censoring across six arms — the only demo that
reaches every module without the user supplying their own data.

Both are now exported and wired into the loader and selector.
`test-demo-data-reachable.R` (32 assertions) closes the loop both ways: every
shipped file must be reachable or explicitly classified, every loader entry must
resolve to a file that parses, and the master demo must retain every optional
column.

### Fixed — Info and Help tabs (D.13)

Corrected, all user-visible:

- "two statistical modeling approaches" → **five** (LME4, GAMM, AUC, Bayesian
  LMM, Bayesian GAMM)
- "Type III ANOVA F-test" → Wald **chi-square**, with why it is a chi-square and
  when that matters
- AUC "Kruskal-Wallis / pairwise Wilcoxon" → **Welch one-way ANOVA** plus an
  equal-variance test
- "emmeans with Tukey adjustment" → comparison family and correction are separate
  settings; the default is an exact Dunnett
- "LOCF can extrapolate AUC to a common endpoint" → extrapolation was removed
  from the statistical path and now affects plots only
- cage "included as a covariate if not collinear" → structure-aware handling
  (crossed / nested-replicated / nested-confounded)
- Power's "Bayesian simulation" mode → removed; the selector offers two modes

Added: the randomisation test and its exchangeability requirement, the five
non-Dunnett corrections, GAMM with its k-check and concurvity diagnostics,
Bayesian GAMM, scale-location and random-effects Q–Q plots, influence measures
and Pareto-k, the endpoint estimand, bootstrap intervals, TWM ranking bounds,
per-arm attrition, competing risks and Aalen-Johansen, dose-response trend tests,
four auto-detected column roles (Date, Sex, Necrotic, Survival_Time), the Height
column, and five of the six volume formulas.

Limitations that could not be fixed are now stated rather than omitted: AUC
penalising early-removed animals, the observed estimand flattering the control
arm, and the chi-square being anti-conservative in small samples.

The Survival Info tab was audited and found already accurate.

### Added — a guard against documentation drift (D.14)

`test-docs-current.R`, 47 assertions pinning the checkable claims to the source of
truth: the documented model count must equal the number of selector choices,
every selectable model must be named, every auto-detected column role must appear
in the Help table, and the removed terms must not reappear. It cannot verify the
prose is good — only that its falsifiable claims are still true, which is the part
that rots silently.

Suite: 478 → 557 passing, no skips.

## [26.07.008] - 2026-07-30

Addresses D.1 (files over 1,000 LOC), reframed after looking at what was in them.

### Changed — `helpers.R` split into four themed files

`helpers.R` (1,227 LOC) became `helpers_data.R` (423), `helpers_format.R` (216),
`helpers_plot.R` (501) and `helpers_ui.R` (126), following the `helpers_*.R`
convention already set by `helpers_bayes.R` / `helpers_report.R`. The partition
was verified lossless before the original was deleted, all 39 functions confirmed
to still resolve, and the suite was unchanged at 314 across the move — which is
what makes it a refactor rather than a rewrite.

### Added — 164 assertions for helpers that had none (D.9)

The split is the means; this is the point. 29 of the 39 functions in `helpers.R`
had no test referencing them, and two of those decide statistical results rather
than presentation:

- **`back_transform()`** inverts the modelling transform. The backend fits
  `log(Volume)`, not `log1p`, so `expm1()` here would be wrong by one unit at
  every point — invisible on a 2,000 mm³ tumour, badly wrong on a 20 mm³ one. Now
  pinned against `exp()` with round-trip tests, and its silent fall-through for an
  unrecognised transform documented rather than left to be found.
- **`synthesize_censor_column()`** fabricates the event indicator the survival
  analysis consumes, by the rule *last observation before study end ⇒ event*. It
  is correct for the mechanism it models, and it has a limitation that cannot be
  fixed in this function: a survivor whose final measurement was skipped is coded
  identically to one that was euthanised. The data frame does not carry what would
  separate them. It does report `n_flagged`, and that contract is now pinned so
  the count stays available to display.

Also covered for the first time: `pick_reference_group()` (its substring matching
promotes an arm called "Uncontrolled" over a genuine "Vehicle" — asserted rather
than silently changed, since requiring exact equality would break the common
"VEHICLE-1" style), `apply_subset_filter()` (exclusions fail safe: a mis-mapped
column drops nothing rather than everything), `detect_date_format()`,
`format_pval()`, `translate_error()` (unrecognised messages pass through verbatim
instead of being replaced by a generic shrug) and `js_string()`.

### Declined — mechanical splits of `ui_modules.R` and the server modules (D.10)

`ui_modules.R` is 1,563 lines of static UI with no reactive coupling; nothing can
hide there. The server modules are dominated by render blocks, and moving a
renderer between files makes it neither more correct nor more testable — what
blocked I.2 for four rounds was `rv$x <- NULL`, not 1,877 lines. The move that
would pay is extracting the decision logic out of `renderUI` into pure functions;
recorded as D.11 rather than done here.

Suite: 314 → 478 passing, no skips.

## [26.07.007] - 2026-07-30

Closes I.2. Doing so turned up why it had resisted four rounds of attempts.

### Fixed — modules wiped their own result slot at init (D.7)

Every module server opened with an unconditional `rv$<x>_results <- NULL`. Since
`rv` is created empty in `app.R` and modules mount once, this never cleared
anything in production. What it did do is make the modules impossible to drive
through `testServer()` with a pre-populated result: the slot was wiped the instant
the module initialised, so every attempt to assert render behaviour got the "run
an analysis first" placeholder back.

That is why I.2 stayed open. The existing tests assert
`inherits(session, "session_proxy")` and that `session$ns()` returns a non-empty
string — they verify testServer built a mock session, and never touch the module.
They were not lazy; the harness made anything better impossible and nobody traced
why. Removed in all six modules; behaviour is unchanged, since a missing key on a
`reactiveValues` reads `NULL` either way. The deliberate reset on new data load is
a different thing and is kept.

### Fixed — the data-load reset named a slot that does not exist (D.8)

`reset_analyses()` cleared `rv$pw_results`, which appears in exactly one place in
the repository: that line. The Power module uses `ap_results`. `tox_results` was
missing from the reset altogether.

So loading a new dataset reset every analysis tab except Power and Toxicity, which
kept showing results from the previous dataset. For Toxicity that is plainly
wrong. For Power it is subtler — the a priori calculation is mostly
data-independent, except that "Seed SD from loaded data" ties the pooled SD to
whatever was loaded, so a stale curve can be calibrated against data the user has
since replaced.

Fifth instance in this review of the same root cause: a name written and never
asserted on fails silently and looks like it works.

### Added — render-logic tests (I.2)

14 assertions that set `rv$<X>_results` to a realistic backend result and read the
rendered output back: the Power summary's dropout statement in both directions, a
pre-K.15 result still rendering, the randomisation panel distinguishing "not run"
from a result, and the TWM ranking note separating a real ranking from arms whose
intervals overlap — identical point estimates in both fixtures, so only the
interval logic can produce the difference.

Suite: 300 → 314 passing, no skips.

## [26.07.006] - 2026-07-30

### Fixed — scale metadata records what the backend did, not what was asked (R3.D6)

R3.D6 was left open pending backend R3.29, which has landed: `me_result_meta()`
makes every analysis function declare `transform_used` and `model_type_used`,
with `"none"` meaning none was applied rather than the field being absent.

The predicted divergence is real. An AUC run submitted with the default
`transform = "log"` returns `transform_used == "none"`, because that path computes
on raw volumes and ignores the argument. The module recorded `input$tg_transform`
and so labelled that output `"log"`.

Resolution is now `result$meta$…` → top-level field → input, so an older backend
degrades to the previous behaviour rather than erroring. Two contract tests pin
it: one asserts the backend's own AUC/lme4 divergence, so the fix cannot rot if
the backend changes; the other asserts the module's resolution order.

Benign today — both consumers are gated on `model_type_used == "lme4"` — and fixed
because that gate is what the next feature forgets.

## [26.07.005] - 2026-07-30

Power module, reached from the backend side. The backend's first working coverage
baseline named `apriori_power_simulation()` at 0.00 %; testing it turned up two
findings that live here.

### Fixed — fallbacks encoded the parameter scale the backend abandoned (D.4)

`mod_power.R` passed `control_growth_rate %||% 15`, `treatment_effect %||% 10`,
`random_slope_sd %||% 3` and two more — values from the percent-scale
parameterisation the backend dropped when its data-generating model moved to the
log scale. All five input IDs exist with correct UI defaults, so these only fire
when an input reads `NULL`, transiently before the client syncs. Had one fired,
`control_growth_rate = 15` makes `exp(15 · t)` overflow to `Inf` by day 14 and the
fit dies pointing nowhere near the cause. Corrected, and `test-power-params.R` now
asserts each fallback equals the UI default beside it so they cannot drift apart
again.

### Fixed — the Method panel named the wrong power driver (D.5)

It described random intercept SD as *"Higher values require larger N."* That is
false: power is identical across a 30-fold range of it (0.05 → 1.50), because the
Treatment × Day slope contrast is orthogonal to per-mouse intercepts. Baseline SD
is inert for the same reason and aliased with it.

Meanwhile `random_slope_sd`, which moves power from 1.00 to 0.15 over 0.01 → 0.15
on the same design, was described only as "between-animal variability in growth
rate". A user calibrating carefully against historical data would have seen the
answer refuse to move. Both entries rewritten with the worked numbers.

### Added — plan against dropout (D.6)

New "Euthanasia volume limit" control, mapping blank/`NA`/non-positive to `Inf` so
the backend's validation never sees a malformed value. The Summary tab reports
realised attrition and states the assumption **even when no dropout is modelled**,
so a complete-data power figure is never read as unconditional. Ignoring dropout
costs roughly 6 points of power at 30 % attrition, 12 at 40 %.

Suite: 251 → 295 passing, no skips.

## [26.07.004] - 2026-07-30

### Added — the Round 3/5 backend outputs are now visible

Six result fields the backend added were being ignored by the UI. Five are now
shown, one partially.

- **Synergy Metrics** — bootstrap interval table under the verdict, with an explicit
  callout saying whether the Bliss-excess 95 % interval excludes zero. That callout
  is the point: the verdict is a categorical label derived from three group means,
  and a "Strong Synergy" whose interval spans zero was previously indistinguishable
  on screen from one that does not. Plus per-arm attrition and an endpoint-estimand
  selector.
- **Therapeutic Window** — the table is sorted by TWM, so it reads as a ranking. A
  note now compares the top two arms' intervals and says outright when they overlap,
  because at typical group sizes they usually do and the order is not resolvable.
  Plus attrition, endpoint day, and the estimand used.
- **Weight-loss KM** — a competing-risks warning when animals left before the study
  ended, stating that 1 − KM *overstates* weight-loss incidence in that case, plus a
  censoring breakdown by type. (The Aalen-Johansen curve is in the result object but
  not yet plotted.)
- **ANOVA tab** — the Brown-Forsythe variance check, framed as "was the
  unequal-variance choice needed" rather than a gate, since the AUC path uses Welch's
  tests either way.
- **ANOVA tab** — the H.3 randomisation test, opt-in via a sidebar checkbox that
  requires you to state the unit randomised. Three interpretive aids matter more than
  the p-value: a warning when the asymptotic p-value is >100× smaller than the
  randomisation one (the norm at these group sizes), a warning when the design cannot
  attain p < 0.05 at all, and the unit and exchangeable-unit count stated plainly.
- **Endpoint estimand selectors** in both the Synergy and Toxicity sidebars,
  defaulting to model-based, with the survivor-bias magnitude quoted (~23 points
  versus ~5 on a simulated study) so the default is justified rather than asserted.

`cage_column` is now passed to `weight_loss_threshold()` so reused IDs across cages
are not collapsed into one subject (backend R3.25).

### Added — backend contract tests

`test-backend-contract.R`: 39 assertions that call the backend the way the dashboard
does and verify every field the new renderers read exists with the expected name.
This is the first test here that checks the dashboard/backend *interface* rather than
module initialisation.

It closes a specific gap. The new renderers read ~20 field names off result objects
and nothing checked them; a rename upstream would have silently blanked a panel —
the same failure class as §R3.D4 (the GAM tab that could never populate) and §R4.D2
(a suite passing against a six-release-stale backend).

### Tests

251 passing (was 212), 0 failing.

**Still the substantive gap:** §I.1 and §I.2. The contract tests verify the
interface, not the Shiny reactives, so the render logic itself — the interval
comparison, the overlap check, the floor warning — is still only verified by reading
it. Driving those through `shiny::testServer` is the next step.

## [26.07.003] - 2026-07-30

### Fixed — R4.D2: tests could silently run against a stale backend

`devtools::load_all()` resolves `Imports: mouseExperiment` against the **installed**
library, not `../mouseExperiment`. During the Round 3 review the installed backend
was 0.4.10 while the source under review was 0.9.0, so this suite reported 210
passing against a backend six releases behind the code being reviewed.

New `test-backend-version.R` compares the installed backend against the sibling
source's `DESCRIPTION` and **fails** when the installed copy is older, naming the fix
(`devtools::install("../mouseExperiment", upgrade = FALSE)`).

It fails rather than skips deliberately: skipping is the failure mode Round 3 spent
its effort removing, and a green run against the wrong backend is worse than no run.
Verified in both directions — passes when aligned, fails when the sibling is bumped
ahead. Succeeds quietly with no sibling checkout (a CI job installing from GitHub),
where `load_all()` has already enforced the `DESCRIPTION` floor.

`CONTRIBUTING.md` now opens its Testing section with the install step.

**Not fixed by this:** the coverage problem underneath. The same 210 tests pass
against backend 0.4.10 and 0.11.0, so the suite still cannot tell a working backend
from a stale one — only that the versions differ. §I.1 and §I.2 remain the
substantive gap.

### Tests

212 passing (210 + the 2 new skew assertions), 0 failing.

## [26.07.002] - 2026-07-29

Dependency declarations (CODE_REVIEW.md Round 4).

### BREAKING — installation now requires the Stan stack

`brms`, `bayesplot`, `future` and `promises` moved from `Suggests` to `Imports`
(R4.D1). All four were optional in `DESCRIPTION` and mandatory in practice: every
analysis module offers a Bayesian model type in its UI, the F.7 job registry exists
solely to keep `future` promises alive across reconnects, and the VPS image installs
all four. Five dead "Bayesian mode requires 'brms', 'future', and 'promises'" banner
blocks removed.

`mouseExperiment` floor raised to `>= 0.10.0`.

### Known issue — the test suite has been running against a stale backend (R4.D2)

`devtools::load_all()` resolves `Imports: mouseExperiment` against the **installed**
library, not the sibling source tree. The installed backend was **0.4.10** while the
source under review was 0.9.0 — six releases and the whole of Round 3 apart. So the
previously-reported "210 passing" was passing against a backend from before any of
that work, and every dashboard test exercising backend behaviour was testing code
that had not changed.

This is the dashboard's version of the brms skip in the backend: a green suite that
is not testing what it appears to test, for a structural reason rather than a coding
mistake. **Not yet fixed** — needs a `Remotes:` entry, a test preamble that installs
the sibling source, or at minimum a helper that fails loudly on version skew. Raising
the floor to `>= 0.10.0` at least makes the skew *visible* (load now errors rather
than silently using 0.4.10) instead of silent.

Re-running against the freshly-installed 0.10.0 gives **210 passing / 0 failing —
the same 210 as against 0.4.10**. No break was hidden, but that identical result is
the real finding: a suite whose outcome does not move across six releases of its
primary dependency is not exercising it, which corroborates §I.1 and §I.2 on thin
coverage from a different direction.

## [26.07.001] - 2026-07-29

Round 3 code review (CODE_REVIEW.md §R3.D). Closes every dashboard finding; four
of the five were blocked on backend work that landed in mouseExperiment v0.6.0.

### BREAKING — displayed p-values will change

- **The dashboard no longer recomputes pairwise contrasts (R3.D1).** It used to
  discard `result$pairwise_comparisons` and rebuild them from `result$model` with
  a hardcoded Tukey/Dunnett adjustment. Because the backend applied *no*
  adjustment at the time (backend R3.1), the same contrast carried one p-value in
  the UI and a different one in an R session — a product reporting two answers to
  one question. The backend now honours `comparison_family` and
  `p_adjust_method`, so the recomputation is deleted and the UI renders the
  package's own table. There is one implementation and the two agree by
  construction.

### Added

- **Correction selector (R3.D3).** `p_adjust_method` was documented by the
  backend with four options and exposed by the UI with none, so users got
  whichever correction each code path happened to hardcode. New
  "Multiple-comparison correction" control offering Dunnett (exact),
  Tukey HSD, Bonferroni, Holm, Benjamini-Hochberg FDR, and none. Invalid
  family/correction pairings are silently substituted with the closest valid
  default rather than surfacing a backend error the user cannot act on.
- **Cage structure and ICC in the model summary.** New CAGE STRUCTURE block
  showing the detected design (crossed / nested-with-replication /
  confounded), where cage ended up (fixed effect / random intercept / omitted),
  the reason, and the cage ICC. The ICC is what tells a reader how much the
  clustering mattered — the design effect on effective sample size is
  1 + (mice per cage − 1) × ICC.
- **COMPARISONS provenance block** stating the family and correction actually
  applied, so the panel reports rather than assumes (partial R3.D7).

### Fixed

- **`adjust = "dunnett"` was silently downgraded (R3.D2).** That is not a value
  emmeans recognises; it fell back to the `dunnettx` approximation while the help
  text claimed exact Dunnett. The backend now maps "dunnett" to the exact
  multivariate-t, and the table header states what was applied instead of
  asserting a default.
- **The GAM "All Comparisons" tab could never populate (R3.D4).** emmeans cannot
  dispatch on a `gamm4` object, so the recomputation errored into
  `log_caught_error()` and the tab rendered "Comparisons not available for
  tukey". Resolved as a side effect of R3.D1 — rendering the backend's table
  removes the failing call entirely.
- **The survival reduction took the last row, not the last day (R3.D5).**
  Reducing to one row per animal is what makes the backend Cox fit valid, but it
  assumed the upload was sorted by day within animal. A file exported grouped by
  measurement session, sorted descending, or with a re-measure appended out of
  order silently used an arbitrary observation's survival time and event
  indicator. Now sorts by day first, and takes the animal's maximum event
  indicator so an event flagged at an intermediate visit is not dropped.

### Changed

- **Cage Effect Handling** gains "Auto (match the design)" as the default and
  "As random effect" as an option; "Include if not collinear" is relabelled
  legacy, since it drops cage whenever it is nested and so discards variance the
  design can estimate.
- **Extrapolation Points** is relabelled "(plots only)". It no longer feeds the
  statistical model — the backend deprecated that in v0.9.0 because appending an
  imputed observation and treating it as measured data made standard errors
  ignore the imputation. Plot-side projection is unchanged and still drawn
  dashed.

### Tests

- 210 passing, 0 failing.

## [v26.06.013] 2026 June 18th (staging)

### Changed

- **One "Diagnostics" tab per Bayesian module instead of a separate
  "MCMC Diagnostics" tab.** Previously the dashboard appended a
  "MCMC Diagnostics" tab when Bayes mode was selected, while the
  frequentist "Diagnostics" tab sat next to it (and in survival's
  case the freq tab was called "PH Assumption"). Now each Bayesian
  module has a single "Diagnostics" tab whose content switches via
  `conditionalPanel` keyed off the module's `_model_type` input:
    - Tumor Growth: freq residual plots ↔ MCMC traces / rank /
      Pareto-k / posterior densities.
    - Dose-Response: drc lack-of-fit + residuals ↔ MCMC convergence
      table + traces.
    - Synergy: per-group Q-Q + boxplot ↔ MCMC convergence + traces.
    - Toxicity (BW): LMM residual plots ↔ MCMC convergence + traces.
    - Survival: Cox PH Schoenfeld checks (previously "PH Assumption"
      tab — now folded in) ↔ MCMC convergence + traces.

  The `appendTab("MCMC Diagnostics", ...)` blocks were dropped from
  all five module servers. The freq Schoenfeld content in Survival
  was renamed from "PH Assumption" to "Diagnostics" so that switching
  from Cox to Bayes doesn't change which tab to look at.

## [v26.06.012] 2026 June 18th (staging)

### Added

- **Computation-time line on every Bayesian Summary tab.** Wall-clock
  duration of the fit (formatted "X min Y s") is appended to the
  Summary text in the Tumor Growth, Dose-Response, Synergy, Toxicity,
  and Survival modules. Sourced from the job registry's `started_at` /
  `completed_at` timestamps, stamped onto the result by
  `attach_bayes_poller()`. No per-module backend changes needed.

### Changed

- **Stan backend selector shows just the engine names** ("rstan",
  "cmdstanr"). Removed the "(default)" / "(faster, requires install)"
  parenthetical labels and the helpText below the selectInput.
- **Survival module: the IACUC-endpoint clarifying note moved from
  above the analysis pane into the top of the Info tab.** Same wording,
  better location — first-time users now see it where they go to read
  about the analysis, instead of every time they open the tab.

## [v26.06.011] 2026 June 18th (staging)

### Fixed

- **MCMC rank plot + trace plot now render on the Bayesian "MCMC
  Diagnostics" tab.** Used `brms::as.array(brmsfit)` to extract the
  3-D draws array, but brms 2.23 (current VPS image) no longer
  exports `as.array`. Replaced with `posterior::as_draws_array()` and
  updated the dimnames key from `$parameters` (old brms convention)
  to `$variable` (posterior convention) in
  `build_mcmc_rank_plot()`. Same swap applied to the 5 inline
  `mcmc_trace(brms::as.array(rv$bayes_models$…))` render blocks in
  the tumor-growth / dose-response / synergy / toxicity / survival
  modules.

### Changed

- **"App stays usable" copy on Bayesian fit alerts is now accurate.**
  The alert previously claimed the dashboard remained interactive
  while the model ran; in practice a progress modal blocks all other
  tabs. Updated the alert in 5 modules (tumor-growth, dose-response,
  synergy, toxicity, survival) to explain that the modal blocks the
  UI but the user can safely close the browser and reconnect later —
  the fit continues on the server and results return on reconnect
  via the F.7 session-resume registry.

## [v26.06.010] 2026 June 17th

### Changed

- **Subset Filter is now applied explicitly via an "Apply Filtering"
  button** on the Data tab. Previously the selectizeInputs mirrored
  directly into `rv$exclude_*` on every keystroke; users only learned
  the filter was "active" when an analysis re-fit picked it up. The
  new flow: edit the three exclusion lists → click Apply → the
  Dataset Summary and Data Preview both immediately reflect the
  filtered dataset, and a status line shows
  `N of M rows shown (K excluded)`. A Clear button empties all three
  lists at once. The static "click Run Analysis to re-fit" hint is
  replaced by a status row that tells the user, at a glance, what is
  and isn't currently filtered.

- **Subset Filter panel moved up** to sit between the Dataset Summary
  and the Data Preview (was below the preview). The natural flow is
  now top-down: see summary → choose what to exclude → click Apply →
  see the result in the preview directly below.

### Internals

- `setup_subset_filter()` no longer auto-syncs `input$exclude_*` → ``rv$exclude_*`. Pending choices now stage in the inputs; only the
  Apply / Clear handlers commit to `rv` and bump a new
  `rv$filter_applied_at` tick.
- `preview_data()` in the data-upload server takes a reactive
  dependency on `rv$filter_applied_at` and passes its output through
  `apply_subset_filter()` so the preview reflects the current applied
  state without depending on `rv$exclude_*` directly.

### Tests

- 210 / 210 pass, no regressions.

## [v26.06.009] 2026 June 17th

### Fixed

- **Diagnostic plot overlap on TG and Toxicity Diagnostics tabs.**
  Plot images were rendering taller than their declared
  `plotOutput(height = "350px")` containers, overlapping the
  "How to read this plot" details element below. Two-part fix:
  - Backend v0.4.10 ships compact-theme ggplots (smaller title,
    tighter margins) so the plot fits its canvas
  - Dashboard wraps each plot in a `div(style = "margin-bottom:
    24px;")` and bumps `plotOutput` height from 350 → 380 px on
    each diagnostic plot (TG Diagnostics × 4 plots; Toxicity / BW
    Diagnostics × 4 plots) so the details element sits visually
    apart from the plot

Requires backend ≥ 0.4.10 to get the compact ggplots.

## [v26.06.008] 2026 June 2nd

Closes the remaining DIAGNOSTICS audit items (2, 10, 11, 12) plus the
dashboard half of (6) and (13). Requires backend ≥ 0.4.9 for the new
result-list fields. After this release all 13 items in
`docs/DIAGNOSTICS.md` are addressed.

### Added — features

- **Pareto-k distribution histogram** (DIAGNOSTICS.md gap 11). In the
  TG MCMC Diagnostics tab, alongside the existing high-k table. Shows
  the full Pareto-k distribution with red (0.7) and orange (0.5) cut-off
  lines so users can see whether high-k observations are isolated or
  part of a broader pattern. New helper `helpers_bayes.R::
  build_pareto_k_plot()`.

- **PPC coverage bar chart** (gap 10). In the TG MCMC Diagnostics tab.
  Replaces the dl-row text with a quick-scan bar chart showing
  empirical vs nominal coverage at 50/80/95% with the nominal target
  drawn as a red dashed reference. New helper
  `helpers_bayes.R::build_ppc_coverage_plot()`.

- **MCMC rank plots** (gap 12). Modern convergence diagnostic that
  complements Rhat. Each chain's draws are rank-normalised; converged
  chains look like uniform histograms and resemble each other. Wraps
  `bayesplot::mcmc_rank_overlay()`. New helper
  `helpers_bayes.R::build_mcmc_rank_plot()`.

- **Bayesian residuals + posterior-density plots rendered** (gap 2).
  The backend computes `residuals_plot` (posterior-mean residuals vs
  day, faceted by treatment, with loess smoother) and
  `posterior_dist_plot` (`bayesplot::mcmc_areas` for treatment
  parameters) — both were computed and silently discarded. Now
  rendered in the TG MCMC Diagnostics tab. Flipped TG's F.7 path to
  `plots = TRUE` so the backend actually generates them.

- **TG Influence sub-tab** (gap 13). Conditional section on the TG
  Diagnostics tab showing Cook's distance + DFBETAS as DT tables.
  Shows only when the backend populated the fields (gated via the
  `tg_show_influence` reactive).

- **BW Influence sub-tab** (gap 13). Same on Toxicity → BW Diagnostics
  tab. `tox_show_influence` reactive gates visibility.

- **Synergy Diagnostics tab** (gap 6). New tab on Drug Synergy with
  the per-group Q-Q and boxplot from `analyze_drug_synergy()`.
  Available only in the frequentist path; falls back to a graceful
  "not available" message in Bayesian / over-time mode.

### Interpretation help — pattern

Every new diagnostic plot has a collapsible `tags$details(tags$summary
("How to read this plot"))` block describing:
- What pattern indicates the assumption is met
- What pattern indicates a problem
- What to do about it

Collapsed by default; expand inline without leaving the tab.

### Tests

- 210 / 210 pass, no regressions.

## [v26.06.007] 2026 June 2nd

Fixes the diagnostic-rendering bugs and adds the high-priority gaps
surfaced by the backend `docs/DIAGNOSTICS.md` audit. Requires
backend ≥ 0.4.8 for the new diagnostic fields.

### Fixed

- **TG Diagnostics tab no longer broken**. The dashboard was
  rendering `rv$tg_results$diag_qq_plot` and
  `rv$tg_results$diag_resid_fitted_plot` — fields the backend's
  `tumor_growth_statistics()` never produced for the LME4 path. Every
  LME4 fit silently showed "not available". Fixed in tandem with the
  backend's v0.4.8 update which now produces the ggplot diagnostics
  on the same field names.

### Added — features

- **Expanded TG Diagnostics tab**: Q-Q + Residuals vs Fitted now
  joined by Scale-Location and Random-effects Q-Q. Each plot has
  collapsible "How to read this plot" help text explaining what
  patterns to look for (S-shape vs bow vs heavy tails vs fan shape
  etc.). Conditional GAM-specific sub-section displays `k.check` and
  `concurvity` tables when `model_type = "gam"`.

- **Schoenfeld residuals plot** in the Survival → PH Assumption tab.
  Plots the smoothed Schoenfeld residuals via the `cox.zph` object's
  built-in `plot()` method. Surfaces the SHAPE of any PH violation
  (linear drift, sinusoidal, sudden change) that the `cox.zph` test
  table can't show. Includes interpretation help.

- **New DR Diagnostics tab** (Dose-Response). Surfaces the
  backend's new `dr_lack_of_fit` table (`drc::modelFit` output) and
  `dr_residuals_plot` (residuals vs dose with loess smoother). Each
  with collapsible interpretation help.

- **Expanded Toxicity / BW Diagnostics tab**: Scale-Location and
  Random-effects Q-Q joining the existing Q-Q and Residuals vs
  Fitted. Plots now use the ggplot2 versions from the backend
  (consistent with TG) instead of base-R `qqnorm` / `plot`. Fallback
  to base-R is kept for older cached results.

- **Inline interpretation help** under every diagnostic plot via
  collapsible `tags$details` blocks. Each describes:
  - What pattern indicates the assumption is met
  - What pattern indicates a problem (S-shape, bow, fan, curved
    trend, etc.)
  - What to do about it (try transform X, switch to model Y, etc.)
  Help is collapsed by default so it doesn't clutter the surface
  for experienced users.

### Tests

- 210 / 210 pass — same as v26.06.006; no regressions.

## [v26.06.006] 2026 June 2nd

Closes CODE_REVIEW.md Round 2 **F.7** — session resume for long-running
Bayesian fits. When the user's tab disconnects mid-fit (NPM proxy
timeout, browser sleep, network blip), the future continues running
in its worker; on reconnect, a banner surfaces the outstanding job and
the user can resume into the right module tab. Six commits implement
the full pattern.

### Added — features

- **Job registry** (`R/job_registry.R`). Process-level
  environment keyed by sid → job_id → entry. Every read API takes
  the sid as the first argument so cross-user enumeration is
  impossible. Exposes `launch_job()`, `get_job()`, `list_jobs()`,
  `cancel_job()`, `retire_job()`, `evict_stale_jobs()`.

- **Session-identity cookie** (`R/session_identity.R` +
  `inst/www/js/session-id.js`). 256-bit CSPRNG token generated by
  the browser on first connect; cookie attributes `SameSite=Strict;
  Secure; Path=/; max-age=86400`. The server reads the cookie from
  `session$request$HTTP_COOKIE` — never from any Shiny input — so
  forging an `input$mE_session_id` message gets the attacker
  nowhere.

- **Resume banner** (`R/resume_banner.R`). Renders nothing
  when the sid has no outstanding jobs. When jobs exist, shows one
  row per job with kind label, elapsed time, status, and Resume /
  Cancel buttons. Click delegation: a single inline JS handler per
  action sets a single module-scoped input with the data-jid
  attribute. Dispatch by `kind` via the `register_resume_handler()`
  registry decouples the banner from any specific module.

- **TTL sweeper** (`R/ttl_sweeper.R`). `later::later()`
  recurring task runs `evict_stale_jobs()` every 60 s. Defaults:
  completed jobs retained 30 min, abandoned running jobs 4 h.
  Idempotent — repeat calls don't stack sweepers.

- **Per-module resume helpers** for all 6 Bayesian fits:
  `bayes_tg`, `bayes_bw`, `bayes_surv`, `bayes_syn`, `bayes_dr`,
  `bayes_twm`. Each module's Bayesian branch (50-80 LOC of inline
  `future_promise` + `%...>%` / `%...!%` chains) replaced with one
  call to `launch_<X>_bayes_job()`. Each module registers its
  resume handler at server startup via `register_<X>_resume()`.

- **Shared runner** (`R/bayes_job_runner.R`).
  `attach_bayes_poller()` and `register_bayes_resume()` extract the
  common poll-and-dispatch logic so per-module files stay small
  (110-150 LOC each instead of duplicating 50+ LOC of polling
  boilerplate).

### Security guarantees (multi-user isolation)

- Registry keyed by sid: `get_job(sid_B, jid_A)` returns NULL.
  Job IDs are never access tokens; they're indices within an
  sid's namespace.
- Server reads sid from the HTTP Cookie header server-side; never
  trusts `input$mE_session_id`. Anti-spoofing test in
  `test-session_isolation.R` greps every R/ file to assert no
  exported code reads that input.
- 256-bit cookie entropy is unguessable. Leak audit greps every
  R/ file to assert no `log_*` / `message` / `stop` / `warning` /
  `showNotification` call references `session_id`. Banner HTML
  also asserted to never contain the sid value.
- Cookie attributes: `SameSite=Strict` blocks cross-site
  submission; `Secure` blocks plain-HTTP transmission.
- Workers don't `eval()` user-supplied strings, so worker memory
  isn't a cross-job leak vector. Tested via the assertion that no
  `eval`/`parse` of user-uploaded data occurs in the Bayesian fit
  paths.

### Tests

- `test-session_isolation.R` (39 assertions) — security-critical
  layers 1-4 audited above.
- `test-resume_banner.R` (10 assertions) — banner empty when no
  jobs, renders kinds correctly, sid not in HTML, dispatch
  routing works.
- `test-tg_resume.R` (11 assertions) — TG-specific success /
  failure / context round-trip.
- `test-resume_e2e.R` (26 assertions) — full disconnect/reconnect
  cycle, cross-sid isolation under dispatch, every Bayesian kind
  has a `register_*_resume`, cancel-race.

Full test suite: 210/210 pass.

### Known caveats

- Registry is in-memory. If the Shiny container restarts (deploy,
  OOM, server reboot), all in-flight fits are lost. This was an
  explicit design choice — `callr::r_bg` or a SQLite-backed
  job table would survive restart but add a child-process
  supervisor / DB dependency. Tracked for a future iteration if
  the need arises.
- Cookie is not HttpOnly because JS sets it client-side. XSS
  anywhere in the dashboard could exfiltrate sids. Mitigation:
  audit user-supplied strings flowing into HTML; set CSP header
  at the NPM layer in a follow-up.

## [v26.06.005] 2026 June 2nd

Closes the last remaining CODE_REVIEW.md item — **2.4 Dose-Response plot
customization controls**.

### Added — features

- **Dose-Response: per-dose centre / error overlay** (Round 1 item 2.4).
  Two new selectors on the Plots tab (centre: None / Mean / Median;
  error bar: None / SEM / SD / 95 % CI) add an overlay layer on both
  the Tumor Volume vs Dose scatter and the TGI subplot, matching the
  established pattern used by the Tumor Growth and Toxicity plot tabs.
  The overlay is computed per dose level and drawn on top of the raw
  scatter / fitted trend line.
- **Dose-Response: jitter width slider** for overlapping points at the
  same dose level. Useful for replicate-rich studies where multiple
  animals share a dose. Default 0 preserves the previous behaviour.

## [v26.06.004] 2026 June 2nd

Closes CODE_REVIEW.md Round 2 C.1 / C.3 for the remaining four mega-files
(`mod_toxicity.R`, `mod_survival.R`, `mod_synergy.R`, `helpers.R`) using the
same by-concern split established for `mod_tumor_growth.R` in v26.06.003.

### Changed

- `mod_toxicity.R`: 2,175 → 1,595 LOC (−27%). Extracted:
  - `R/mod_toxicity_downloads.R` (236 LOC) — all 17 download handlers
    bundled into `add_tox_download_handlers(input, output, session, rv)`.
    The active-plot download still reaches into `output[[plot_id]]()` —
    that cross-binding works because every handler shares the same
    `output` env.
  - `R/mod_toxicity_plots.R` (399 LOC) — 10 frequentist plot renderers
    (Q-Q, residuals, BW forest, BW trajectory, BW AUC, KM, ET TWM, ET
    scatter, ET benefit, ET weight-corrected TGI) plus the BW and ET
    group-order/palette UIs. Bayesian renderers stay in main because
    they're interleaved with the Bayesian tab-insertion observer.
    Exposed as `add_tox_plots(input, output, session, rv,
    ordered_tox_grps, ordered_et_grps)`.
- `mod_survival.R`: 1,554 → 1,467 LOC (−6%). Extracted:
  - `R/mod_survival_downloads.R` (108 LOC) — all 6 download handlers
    (2 Bayesian + 4 frequentist). `km_plot_r` is passed in explicitly
    so the KM-PNG handler's reactive dependency stays visible at the
    call site. Exposed as `add_surv_download_handlers(input, output,
    session, rv, km_plot_r)`.
- `mod_synergy.R`: 1,242 → 1,171 LOC (−6%). Extracted:
  - `R/mod_synergy_downloads.R` (103 LOC) — all 6 download handlers
    (1 Bayesian + 1 CSV + 4 plot PNGs). The four plot PNG handlers
    each pull from a `current_syn_*_plot` reactiveVal that the
    `renderPlot` bindings update; those reactiveVals are passed in
    explicitly. Exposed as `add_syn_download_handlers(input, output,
    session, rv, current_syn_main_plot, current_syn_trend_plot,
    current_syn_ci_plot, current_syn_tgi_plot)`.
- `helpers.R`: 1,417 → 1,153 LOC (−19%). Extracted:
  - `R/helpers_bayes.R` (162 LOC) — `bayes_diagnostics_panel()` +
    `bayes_influential_obs()`. Used by every Bayesian module.
  - `R/helpers_report.R` (112 LOC) — `build_html_report()`. Composes a
    standalone HTML audit document.

No behaviour change anywhere — implementation is moved verbatim. Full
test suite still passes 124/124.

## [v26.06.003] 2026 June 1st

Closes CODE_REVIEW.md Round 2 C.1 / C.3 — the tumor-growth mega-file refactor.
`R/mod_tumor_growth.R` shrinks from 2,659 LOC to 1,727 LOC (−932 LOC, −35%);
the extracted code lives in two new files split by concern, not by tab.

### Changed

- `R/mod_tumor_growth_downloads.R` (new, 253 LOC). Holds all 16
  download handlers (anova / comparisons / effects / rates / AUC /
  data summary / 5 plot PNGs / posterior predictive check /
  necrosis CSV / "Download All" ZIP) plus the necrosis-summary
  `renderDT` that sits adjacent to its download. Exposed as
  `add_tg_download_handlers(input, output, session, rv,
  tg_growth_plot, tg_rates_plot, tg_auc_plot, tg_forest_plot,
  tg_wf_plot)` — the plot reactives are passed in explicitly rather
  than captured by closure so the cross-file contract is visible.
- `R/mod_tumor_growth_plots.R` (new, 746 LOC). Holds all five plot
  reactives (`tg_growth_plot`, `tg_rates_plot`, `tg_auc_plot`,
  `tg_wf_plot`, `tg_forest_plot`), their associated `renderPlot`
  bindings, the group-order / colour `renderUI`, the palette-preset
  observer, and the two frequentist diagnostic plots (`tg_diag_qq`,
  `tg_diag_resid_fitted`). Exposed as `add_tg_plots(input, output,
  session, rv, ordered_tg_grps)` and returns a named list of the
  five reactives so the caller can hand them to the downloads helper.
- `R/mod_tumor_growth.R` now owns lifecycle only — observers, fit
  dispatch (frequentist + Bayesian), status messages, tab insertion,
  the Bayesian renderers (`tg_posterior_summary`,
  `tg_pp_check_plot`, `tg_mcmc_diag_table`, `tg_bayes_diagnostics`,
  `tg_bayes_influential`, `tg_mcmc_trace_plot`), and the table /
  drilldown / compare-panel renderers that remain interleaved with
  lifecycle code.

No behaviour change: implementation is moved verbatim, and the
existing `test-module_servers.R::tumor_growth_server initialises and
exposes session ns` test still exercises the full module startup
across all three files.

## [v26.06.002] 2026 June 1st

Closes 11 more CODE_REVIEW.md Round 2 items. After this release only the
mega-file refactor (C.1 / C.3) and the cross-session resume infrastructure
(F.7) remain open in CODE_REVIEW.md Round 2; both need a dedicated
architectural pass.

### Added — features

- **Subset Filter** on the Data Upload tab (F.5). Three selectizeInputs
  (Exclude IDs / Treatments / Cages) wired through to every analysis
  module via the new `apply_subset_filter()` helper. Re-running an
  analysis with a different subset is now a one-click operation — no
  CSV editing or re-upload required. Filter state survives session
  save/restore.
- **Per-mouse drill-down** on the Tumor Growth Growth Rates tab (F.4).
  Selecting any row in the per-mouse table now renders that mouse's
  volume trajectory below, overlaid with the dashed group-mean
  trajectory and a subtitle reporting its growth rate + back-transformed
  daily fold change.
- **Side-by-side compare** of consecutive Tumor Growth runs (F.3).
  Re-running TG with a different model type (LME4 → Bayesian, AUC →
  GAMM, …) preserves the previous result. A new Compare tab renders
  both treatment_effects tables side-by-side, each tagged with its
  `model_type_used`.
- **Reproducible HTML report** from the Tumor Growth Summary tab (F.1).
  Single self-contained file with timestamp + settings + Treatment
  Effects + Pairwise + ANOVA + Growth Rates + Data Summary. Printable
  to PDF from the browser; no new package dependencies (composed via
  the new `build_html_report()` base-R helper).
- **Settings preset save / load** in the Session Save / Restore panel
  (F.2). Small JSON preset (column mappings + time / volume / formula
  settings) is portable across datasets with the same headers.
- **Dose-Response plot customisation** controls (Round 1 item 2.4).
  Adds legend-position selector + X-axis / Y-axis / Legend label text
  inputs to match the Tumor Growth / Survival / Synergy plot tabs.

### Added — tests

- New `tests/testthat/test-modal-and-cachebust.R` covering
  `show_analysis_modal()`, the `updateTabsetPanel(..., selected = "Summary")`
  auto-switch behaviour, and the `?v=<version>` cache-busting on static
  assets (I.3).
- Extended `tests/testthat/test-shiny_app.R` with an end-to-end workflow
  test: load combo demo → switch to TG → click Run → assert
  `rv$tg_results` populates within 30s on the LME4 path (I.5).
- Init tests for `tumor_growth_server`, `toxicity_server`, and
  `data_upload_server` in `test-module_servers.R` (I.2). Test count for
  that file is now 30 (was 24).

### Changed

- **NPM proxy-timeout startup warning** when running under shiny-server
  (E.3). `runDashboard()` detects `SHINY_PORT` or `SHINY_SERVER_VERSION`
  env vars and emits an 8-line reminder banner about the
  `proxy_read_timeout >= 3600s` + `allow_websocket_upgrade = true` NPM
  requirements documented in Round 2 A.3. No effect for local
  development.
- **Module file naming standardised** (C.2). All `*_modules.R` files
  renamed to `mod_*.R` to match the existing `mod_data_upload.R` /
  `mod_tumor_growth.R` convention. Two internal references updated.

### Removed
- `.waylog/` (Claude Code's local session log directory) and `.claude/`
  added to `.gitignore`; one stale waylog markdown that had been
  accidentally committed in v26.05.023 is removed.

### Open — still deferred
- C.1 + C.3: split mega-files into per-tab UI / server / renders.
  `mod_tumor_growth.R` is still 2,400+ LOC; the split needs paired
  design with the analysis-module pattern and is a session of its own.
- F.7: cross-session resume of in-flight Bayesian futures (promise
  handle persistence keyed by a session ID / cookie). Substantial
  infrastructure piece requiring a server-side design decision.

## [v26.06.001] 2026 June 1st

Closes 11 items from CODE_REVIEW.md Round 2, primarily the D-series
"surface diagnostics that the mouseExperiment v0.4.5 / v0.4.6
backend just shipped". Requires mouseExperiment >= 0.4.6.

### Added — Bayesian diagnostics across every brms tab
- **NUTS sampler warning banner** (Round 2 D). New uiOutput at the
  top of every MCMC Diagnostics tab. Warning alert when any of
  n_divergent > 0 / n_max_treedepth > 0 / E-BFMI < 0.3 is detected,
  with the standard Stan remediation pointers (raise adapt_delta,
  reparameterise, increase max_treedepth). Success alert otherwise.
  Fed by mouseExperiment v0.4.5's nuts_diagnostics field.
- **Model fit & diagnostics panel** in every MCMC Diagnostics tab.
  Single block summarising Bayes R² (Estimate + 95 % CrI),
  Posterior Predictive Check coverage (empirical 50 / 80 / 95 %
  interval coverage), and PSIS-LOO (elpd_loo ± SE plus n_high_k
  count of observations with Pareto-k > 0.7). Pulled from
  bayes_R2 / ppc_coverage / loo_diagnostics result fields.
- **Influential observations table** below the per-parameter Rhat /
  ESS table in every MCMC Diagnostics tab. Lists the per-mouse
  Pareto-k > 0.7 values from PSIS-LOO; empty-state shows a clean
  "No observations exceed the threshold" message.
- Two new internal helpers in `R/helpers.R`:
  - `bayes_diagnostics_panel(result)` — rendering for the warning
    banner + stats panel (used by all five Bayesian modules).
  - `bayes_influential_obs(result, threshold = 0.7)` — extracts the
    per-mouse Pareto-k table from `loo_diagnostics$pareto_k`.

### Added — Survival diagnostic surfaces
- **PH Assumption tab** in the Survival results panel. Renders the
  full Schoenfeld-residual cox.zph table (per-covariate + GLOBAL row,
  with chisq / df / p columns). Info banner explains unavailability
  for the Firth (coxphf) and log-rank fallback methods.
- **Concordance / C-index** rendered above the Cox HR table with a
  95 % CI from sqrt(var). Frequentist Cox path only.
- **PH check global summary banner** above the Cox table — coloured
  alert (warning when global cox.zph p < 0.05, success otherwise),
  with a pointer to the new PH Assumption tab for details.

### Added — Bayesian P_direction column in pairwise tables
- Tumor Growth pairwise table now puts `Lower_CrI`, `Upper_CrI`, and
  `P_direction` in their preferred positions next to estimate / p.value.
  New roxygen bullets in the tab_info explain `P_direction` as the
  directional posterior probability replacing the awkward
  "does the 95 % CrI exclude zero?" check.

### Added — cmdstanr backend toggle in every Bayesian sidebar
- Every Bayesian module (TG / Survival / Synergy / Dose-Response /
  Toxicity BW) gains a "Stan backend" selectInput choosing between
  RStan (default) and cmdstanr. The choice is threaded through to the
  corresponding bayesian_*() call's new `backend = ...` argument
  (mE v0.4.6 E.1). RStan default keeps out-of-the-box use unchanged;
  cmdstanr compiles 3-10× faster on the VPS but needs a one-time
  CmdStan toolchain install via `cmdstanr::install_cmdstan()`.

### Added — Bayesian TWM diagnostic plots
- The TWM result section gains a "TWM Diagnostics" tab with three
  side-by-side panels (Posterior Predictive Check, Prior vs Posterior,
  MCMC Trace) showing both input models (TG on left, BW on right).
  Backed by mouseExperiment v0.4.6 G.5's aliased
  pp_check_tg/bw, prior_posterior_tg/bw, mcmc_trace_tg/bw fields.

### Added — Power Analysis Summary tab
- First results tab is now "Summary" (was "Info"). Displays method
  (Analytic submode or LMM Simulation), scenario count, and the
  recommended N per group (smallest N reaching target_power, derived
  from `power_curve_data` for LMM mode or `Required_N` for analytic
  mode). Both completion branches now call
  `updateTabsetPanel(..., selected = "Summary")` matching every other
  module's auto-switch behaviour.

### Changed — Tests
- `tests/testthat/test-module_servers.R`:
  - The four `expect_true(TRUE)` placeholders are replaced with real
    assertions on the mock session + namespace wiring; the new
    assertions catch the "module crashes silently on init" bug class.
  - `make_test_rv()` gains a `with_results = character(0)` argument
    that pre-populates `rv$<X>_results` with small dummy result lists
    for any subset of `{tg, surv, syn, pw, dr, tox, ap}`. Lets tests
    of post-completion UI behaviour (auto-switch, downloads enabled)
    set up the precondition without running the analysis.
  - The `analysis_log_server` test block is gone alongside the
    function it tested.

### Removed
- `analysis_log_server()` and its two tests. The Log tab was removed
  in v26.05.020 but the unused module + tests outlived it. Structured
  stdout logging via the `logger` package continues unchanged.

### Docs
- Backfilled the v26.05.020 CHANGELOG entry (Log tab removal +
  Help tab rewrite + Umami repoint).

### Requires
- **mouseExperiment >= 0.4.6** for the new diagnostic fields
  (nuts_diagnostics, bayes_R2, loo_diagnostics, ppc_coverage,
  P_direction, ph_test, c_index, TWM aliased plots, backend arg).

## [v26.05.023] 2026 May 30th

### Added
- **Tumor Growth — GAMM model option** — new "GAMM (smooth time × treatment)" and
  "Bayesian GAMM (brms + smoothers)" choices in the Tumor Growth Statistical Model
  selector. The frequentist GAMM path calls `mouseExperiment::tumor_growth_statistics(model_type = "gam")`
  (gamm4 backend); the Bayesian GAMM path passes `model_type = "gam"` to
  `bayesian_tumor_growth()`. Both fit a group-specific smoother on Day with
  auto-chosen basis dimension. Use these whenever growth is materially non-linear —
  the GAMM smoother is the principled replacement for fitting higher-order
  polynomials in time and extrapolates more honestly.
- **Tumor Growth — "Display on original scale" plot toggle** — new checkbox in
  the Plots-tab Elements panel for Growth Rates, Treatment Effects, and Effect
  Size Forest plots. Back-transforms log-scale model output to interpretable
  original units: growth rates become daily fold changes (`exp(rate)`); forest
  contrasts become fold changes vs. reference rendered on a log10 axis centred
  on 1. CIs are exponentiated from log-scale bounds (asymmetric in fold-change
  space, as expected). Inference is still performed on the log scale; the
  toggle only changes display. For sqrt-transformed fits the toggle is a no-op
  on the forest plot since back-transforming a difference of √volumes has no
  clean interpretation.

### Changed
- **Tumor Growth — model summary text** — the "Model summary" tooltip now
  describes the four model families (LME4 / GAMM / AUC / Bayesian) and reports
  the smoother basis dimension `k` when GAMM is selected.
- **Tumor Growth — backend call** — the dashboard no longer passes
  `polynomial_degree = 1L` to `tumor_growth_statistics()` (the backend argument
  was removed in mouseExperiment v0.4.4 in favour of `model_type = "gam"`).
- **Tumor Growth — analysis-scale axis labels** — plots that display values on
  the model scale now declare it in the axis title so readers don't mistake
  log/sqrt-scale numbers for raw volumes. Growth Rates and Treatment Effects
  default y-axis reads "Growth Rate (log-volume / day)" (or "Daily Fold Change"
  when the new back-transform toggle is on); the Forest plot x-axis appends
  " (log scale)" / " (sqrt scale)" when a non-trivial transform was used; the
  Q–Q and residuals-vs-fitted diagnostic axes carry the same suffix.

### Fixed
- **Tumor Growth — `back_transform()` helper** — the helper claimed to invert
  `log1p` (via `expm1`), but `tumor_growth_statistics()`, `bayesian_tumor_growth()`,
  and the GAMM path all transform with plain `log()` (volumes are strictly
  positive). The helper now uses `exp()` so it matches what the rest of the
  pipeline does. It had no callers in the dashboard before this release, so no
  downstream numbers change.

### Requires
- **mouseExperiment >= 0.4.4** (for the new `model_type = "gam"` argument on
  `tumor_growth_statistics()` and `bayesian_tumor_growth()`).

## [v26.05.022] 2026 May 29th

### Changed
- **All modules — auto-switch to Summary on completion** — when an analysis finishes (frequentist or Bayesian, success or error), the results tabset now jumps to the **Summary** tab automatically via `updateTabsetPanel()`. Applied to Tumor Growth, Survival, Drug Synergy, Dose Response, and Toxicity (TWM Bayesian, BW Bayesian, and frequentist paths). Power Analysis has no Summary tab and was left alone.

## [v26.05.020] 2026 May 20th

### Removed
- **Analysis Log tab** — the user-facing **Log** panel is gone. Structured
  logging via `logger` continues to emit to stdout (visible in
  `journalctl -u shiny` / `docker logs`), so the operational signal is
  preserved; only the in-app display tab was removed. `analysisLogUI` is no
  longer registered in `app.R`.

### Changed
- **Help tab — rewritten end-to-end.** Now covers every analysis module
  with the same level of detail: LME4 / AUC / Bayesian Tumor Growth,
  Survival (KM / Cox / AFT), Body Weight + Therapeutic Window Metric,
  Bliss / Loewe / Bayesian Synergy, Hill / Emax Dose-Response, and Power
  Analysis. Also adds: a column auto-detection table (which header names
  the dashboard recognises for ID / Cage / Day / Volume / Weight /
  Treatment / Survival_Censor), a step-by-step workflow guide
  (upload → map → configure → run → export), and an expanded
  troubleshooting section that explicitly addresses the long Bayesian
  fit times (3-12 min, expected during Stan compilation).

## [v26.05.021] 2026 May 29th

### Changed
- **All modules — analysis progress indicator** — the bottom-right `Progress$new()` notification widget has been replaced with a centered modal dialog (`show_analysis_modal()`), shown for the duration of every analysis (frequentist and Bayesian) for consistent, prominent feedback.
- **Static assets — cache busting** — CSS/JS includes now carry a `?v=<package-version>` query string so browsers fetch fresh assets after every release instead of serving stale cached `custom.css`/`overlay.js` (which previously kept old styling/behaviour live in clients after a deploy).
- **Analytics** — Umami tracking repointed to the new `analytics.mouseexperiment.com` instance.

### Fixed
- **All modules — spurious "required column" error / stuck overlay** — removed `progress$set()` calls left orphaned after the progress-widget-to-modal migration. These threw `object 'progress' not found`, which surfaced as a misleading "A required column was not found" error in Tumor Growth (LME4/AUC) and, in Survival, Drug Synergy, Power Analysis, and Dose-Response, crashed the analysis observer so the modal never cleared and no results appeared.
- **Bayesian analyses — full-screen overlay** — removed the redundant in-page `shinycssloaders::withSpinner()` status block from every Bayesian path; an empty, never-completing spinner there produced a loading overlay that covered the whole screen. The centered modal now conveys the running state alone.
- **Analysis modal — backdrop dimming & scroll lock** — the Bootstrap modal backdrop is now rendered transparent and non-interactive (CSS plus a `MutationObserver` that catches the element the moment Bootstrap inserts it), so the modal no longer dims the page; also fixed a case where `modal-open` could persist on `<body>` and leave the page unscrollable after a modal closed.
- **Toxicity — Bayesian body-weight modal** — `on.exit(removeModal)` was tearing down the Bayesian modal the instant the branch returned; the `future_promise` callbacks now own the modal lifecycle.

## [v26.05.019] 2026 May 20th

### Fixed
- **Toxicity — cage column fallback** — `cage_column = "Cage"` hardcoded fallback replaced with `NULL`; no longer causes a missing-column error when no cage column is mapped.
- **All modules — plot download validation** — download handlers now check that the specific plot object is non-NULL before writing the file; a user-visible warning notification is shown if the plot is not available, preventing silent corrupt-file downloads.
- **Survival — volume-threshold event derivation** — added early guard when `event_src == "volume"` and no volume column is mapped; shows a clear error notification instead of a cryptic internal crash.
- **All modules — date format detection** — `detect_date_format()` returning `NULL` now emits a warning notification at every call site before analysis proceeds, preventing silently corrupt `Day` values.
- **Survival — preprocessing consistency** — analysis loop now reads from `pp$processed` throughout rather than mixing `rv$uploaded_data` and preprocessed data.
- **Synergy — colour vector** — inline colour construction in the synergy module replaced with `build_color_vector()`; user-selected group colours now apply correctly to all synergy plots.

### Changed
- **All modules — tab naming** — result tab names standardised: primary output tab → "Summary", model detail → "Model Results", pairwise output → "Pairwise Comparisons", across Survival, Toxicity, Synergy, and Dose-Response.
- **All modules — download button placement** — all download buttons moved to the top-right of their respective tab panels for consistent discoverability.
- **All modules — reference group label** — standardised to "Reference Group" across Survival, Toxicity, and Synergy sidebar selectors.
- **Survival / Synergy / Dose-Response — progress steps** — added minimum three progress checkpoints (data validation ~10 %, model fitting ~60 %, results preparation ~95 %) to match Toxicity; progress bar no longer appears stuck at 20 % for most of the run.
- **All modules — download filenames** — renamed from generic (`kaplan_meier_plot.png`) to date-stamped and reference-group-qualified (`km_plot_Control_20260520.png`) across all download handlers.
- **All modules — error messages** — raw R error strings passed through `translate_error()` before display; common patterns (missing column, singular fit, convergence failure, insufficient data) produce actionable user-friendly messages.
- **All modules — download buttons** — buttons are `shinyjs::disable()`d when required analysis results are not yet available, with tooltip copy explaining what must run first; re-enabled on analysis completion.
- **All modules — JS analytics escaping** — error message strings embedded in `trackAnalysis()` / `trackError()` JS calls now serialised via `jsonlite::toJSON()` (`js_string()` helper), preventing broken JS from multi-line or unicode error messages.
- **Helpers — preprocessing cache** — `ensure_preprocessed()` cache check replaced by eager preprocessing at upload time in `mod_data_upload.R`; no longer re-runs `calculate_volume()` / `calculate_dates()` on first call from each analysis module.
- **All modules — group ordering memoization** — `get_ordered_groups()` result wrapped in `reactive()` and called via the reactive reference in all plot renders; no longer re-parses all Shiny input IDs on every render tick.
- **All modules — brmsfit memory handling** — `brmsfit` objects are now stored exclusively in `rv$bayes_models` (keyed `$tg`, `$surv`, `$bw`, `$syn`, `$dr`); the `$model` field is nulled from result lists and `gc()` called after extraction. This prevents 50–200 MB brmsfit objects accumulating in session reactive values.
- **Survival — KM plot reactive** — replaced duplicated ~170-line KM build block shared between `renderPlot` and `downloadHandler` with a single `km_plot_r <- reactive({...})`; both outputs now share one code path.

## [v26.05.018] 2026 May 20th

### Changed
- **Compatibility with mouseExperiment v0.4.3 API changes:**
  - Bayesian result tables now use `Lower_CrI`/`Upper_CrI` column names instead
    of `Lower_CL`/`Upper_CL`. The non-`lme4` display path in `mod_tumor_growth.R`
    renders columns generically — no dashboard code changes required.
  - `bayesian_survival()` parameter `include_frailty` renamed to
    `include_cage_effect`. No dashboard module calls this parameter directly.
  - `model_type_used` for `bayesian_tumor_growth()` changed from `"bayes"` to
    `"bayes_tg"`. Update any string comparisons on this value accordingly.

## [v26.05.017] 2026 May 20th

### Added
- **Synergy — Bayesian Combination Synergy (brms)** — new `"Bayesian (brms)"` option in the Drug Synergy analysis mode (uses `bayesian_synergy()` v0.4.2). Sidebar controls: Transform (log / sqrt / none), Prior Strength preset, Chains, Iterations. Result tabs: TGI Summary (per-group posterior median ± 95 % CrI), Bliss Summary (excess posterior distribution + P(synergy)), Loewe Summary (CI posterior distribution + P(CI < 1) + interpretation string), Synergy Table (all six rows including Bliss-expected and Loewe-expected), MCMC Diagnostics, Synergy Plot (bar chart with Bliss/Loewe reference lines), Posterior Distribution (density overlays for Bliss excess and Loewe CI).
- **Toxicity — Bayesian Therapeutic Window (brms)** — new `"Bayesian TWM (brms)"` option in the Therapeutic Window Metric mode (uses `bayesian_therapeutic_window()` v0.4.1). Requires both a Bayesian tumor growth result (`return_model = TRUE`) and a Bayesian body weight result (`return_model = TRUE`) to already be available in the session. Result tabs: TWM Table (per-group TWM median + 95 % CrI), TGI Summary, Weight-Loss Summary, TWM Plot (forest plot with break-even reference line at TWM = 1), Efficacy vs Safety Scatter (TGI vs weight-loss % with TWM = 1 isoline).
- **Tumor Growth — Bayesian Dose-Response (brms)** — new `"Bayesian Hill/Emax (brms)"` option in the Dose-Response analysis mode (uses `bayesian_dose_response()` v0.4.0). Sidebar controls: Prior Strength, Chains, Iterations. Result tabs: DR Parameters (EC50/Emax/Hill posterior median + 95 % CrI), Dose-Response Summary, MCMC Diagnostics, Dose-Response Curve (posterior median ± CrI ribbon + EC50 line), Prior vs Posterior.

## [v26.05.016] 2026 May 19th

### Changed
- **Tumor Growth — Bayesian LMM: cage random-intercept and prior controls** — updated to reflect `bayesian_tumor_growth()` v0.3.7 changes:
  - Prior Strength selector now includes `"Skeptical"` as a new option and sets it as the default (replaces `"Weakly informative"`). A `"Manual"` option is reserved for future sidebar expansion.
  - New **Cage Effects** checkbox (default: on) in the Bayesian sidebar, wired to the `include_cage_effect` argument. Displayed only when a cage column is mapped in Data Upload.
  - **Prior vs Posterior** plot tab dynamically populated from `res$prior_posterior_plot` (available when `plots = TRUE`); shows density overlay of prior (grey) and posterior (blue) per treatment coefficient.
  - **Residuals vs Day** plot tab dynamically populated from `res$residuals_plot`; faceted by treatment group with loess smoother.

### Added
- **Survival — Bayesian AFT (brms)** — new `"Bayesian AFT (brms)"` option in the Survival Statistical Model selector (uses `bayesian_survival()` v0.3.8). Sidebar controls: Family (Weibull / Log-normal / Exponential / Gamma), Prior Strength, Chains, Iterations, Include Cage Frailty checkbox. Result tabs: Treatment Effects (Time Ratio + HR), Posterior Summary, MCMC Diagnostics, Survival Curves (posterior credible bands + KM), Prior vs Posterior.
- **Toxicity — Bayesian Body Weight (brms)** — new `"Bayesian Body Weight (brms)"` option under Toxicity → Body Weight mode (uses `bayesian_body_weight()` v0.3.9). Result tabs: Treatment Effects, Weight Loss Summary (posterior % change first-to-last day with 95 % CrI), Weight Trajectory (posterior median ± CrI ribbon + observed points), MCMC Diagnostics, Prior vs Posterior.

## [v26.05.015] 2026 May 14th

### Added
- **Tumor Growth — Bayesian LMM (brms) analysis mode** — new "Bayesian LMM (brms)" option in the Statistical Model selector. When selected, the analysis runs asynchronously in a background worker via `future`/`promises` so the app remains fully responsive during the 3–12 minute Stan compilation and MCMC sampling. A spinner notification is shown while the model runs.
  - Sidebar controls: Prior Strength preset (Weakly informative / Informative / Diffuse), MCMC chains (default 4), Iterations per chain (default 2 000), and an amber timing notice.
  - Random Effects selector is now visible for both `lme4` and `bayes` modes (previously lme4-only).
  - Two Bayesian-specific result tabs are dynamically inserted when Bayesian mode is active and removed when switching away: **Posterior Check** (posterior predictive density overlay via `brms::pp_check`) and **MCMC Diagnostics** (Rhat, Bulk_ESS, Tail_ESS, Converged flag table + trace plots).
  - The **ANOVA** tab displays an informative "not applicable" message in Bayesian mode (ANOVA has no meaning for posterior inference).
  - The **Diagnostics** tab (Q–Q and Residuals vs Fitted plots) redirects users to the Posterior Check tab in Bayesian mode.
  - Treatment Effects and Pairwise Comparisons tabs populate from the same column structure as the lme4 path (95 % CrI for Bayes, 95 % CI for lme4).
  - **Bayesian LMM Mode** section added to the Tumor Growth Info tab describing priors, outputs, and interpretation guidance.
- **Async execution infrastructure** — `future::plan(multisession, workers = 1L)` initialised at app startup in `app.R`; cleaned up on exit via `on.exit`. Prevents blocking the Shiny event loop during long-running Stan fits.
- `future (>= 1.33.0)`, `promises (>= 1.2.0)`, `brms (>= 2.19)`, `bayesplot (>= 1.10)` added to `Suggests` in DESCRIPTION.

## [v26.05.014] 2026 May 14th

### Changed
- **Power Analysis — LMM Simulation defaults and labels** — growth-rate and variance-component inputs now use log-scale defaults (`Control log-growth rate = 0.15`, `Treatment reduction = 0.10`, `Random intercept SD = 0.20`, `Random slope SD = 0.05`, `Residual SD = 0.10`) matching the corrected `apriori_power_simulation()` which generates log-normal data and fits `log(Volume)`.
- **Toxicity — Therapeutic Window Metric formula** — TWM denominator now reads "Mean Body Weight Loss" throughout the UI (info panel, tab description, plot subtitle) and the summary text correctly references `Mean_Pct_Weight_Loss` and `TGI` columns from the updated package.

### Fixed
- **Toxicity — TWM summary text** — the combined-report text block referenced non-existent columns `TGI_Pct` and `Max_Weight_Loss_Pct`; corrected to `TGI` and `Mean_Pct_Weight_Loss`.

## [v26.05.013] 2026 May 14th

### Removed
- **Power Analysis — Post-hoc (observed data) mode** — eliminated in line with the removal of `post_power_analysis()` from the `mouseExperiment` package. Post-hoc power is statistically uninformative (a monotone function of the p-value) and should not be offered as an analysis option. Removed: the "Post-hoc" radio option, post-hoc sidebar controls (method selector, simulation count, custom effect sizes), post-hoc results panel (Summary, Effect Sizes, Power Analysis, Sample Size Recommendations, Power Curves tabs), and all associated server outputs and download handlers (`pw_summary`, `pw_effect_sizes_table`, `pw_power_table`, `pw_sample_sizes_table`, `pw_plot_curves`, `download_pw_effect_sizes`, `download_pw_power_table`, `download_pw_sample_sizes`, `download_pw_plot`).

### Changed
- **Power Analysis — default mode** changed from "Post-hoc (observed data)" to "A priori — Analytic".
- **Power Analysis — Target Power** panel now shown only for "A priori — Analytic" mode (previously also shown for post-hoc mode).
- **Power Analysis — results area** is now a single unconditional tabset showing Scenario Table, Power Curves, and SD Sensitivity (previously two conditional panels switching between post-hoc and a priori tabs).

## [v26.05.012] 2026 May 5th

### Added
- **Power Analysis — A priori Analytic mode** — new sidebar mode alongside the existing post-hoc. Users enter Cohen's d values (one or several, comma-separated) or a raw mean difference + pooled SD. Supports two-sample t-test (two groups) or one-way ANOVA omnibus (k ≥ 3 groups via `pwr::pwr.anova.test`). Sub-modes: "Find N" (required N per group for each target power × alpha combination) and "Find Power" (achieved power at user-specified N values). A "Seed SD from data" button fills the pooled SD from the currently loaded dataset. For mean-difference inputs, an SD sensitivity table shows how required N changes as the assumed SD is varied ±20% and ±40%, communicating uncertainty in the pilot SD estimate.
- **Power Analysis — A priori LMM Simulation mode** — simulates longitudinal tumour-growth data from user-supplied growth parameters (baseline volume, control rate, treatment rate-reduction, random effects, residual SD, timepoint schedule). Fits `Volume ~ Treatment * Day + (Day|ID)` via `lme4::lmer` and uses a likelihood-ratio test for the Treatment × Day interaction; power = proportion of simulations where p < α. Evaluates power across a user-specified range of N per group. Progress bar updates during each (N, alpha) scenario.
- **Power Analysis — Analysis Mode selector** — `radioButtons` at the top of the sidebar switches between "Post-hoc (observed data)", "A priori — Analytic", and "A priori — LMM Simulation". Result panels are conditionally displayed; post-hoc tabs only show for post-hoc mode, a priori tabs (Scenario Table, Power Curves, SD Sensitivity) only show for a priori modes.

## [v26.05.011] 2026 May 5th

### Added
- **Power Analysis tab** — activated as a full tab; previously rendered as a greyed-out "Coming soon" placeholder even though the module was fully implemented. The tab is now enabled whenever data is loaded (same condition as Tumor Growth). The analysis itself requires preprocessed data from a prior Tumor Growth run; a warning is shown if that has not been done yet.

## [v26.05.010] 2026 May 5th

### Removed
- **Toxicity — Log Cell Kill efficacy metric** — removed from the ET mode sidebar radio buttons, the Info tab Efficacy Metrics list, the scatter-plot label switch, and the Shared Parameters description; matches the removal in mouseExperiment 0.3.2

## [v26.05.009] 2026 May 5th

### Added
- **Toxicity — BW Trajectory combined view** — new display mode "Individual + Group Mean" overlays individual mouse trajectories (faint) with bold group-mean lines; default display changed to "Individual mice"
- **Toxicity — BW AUC plot controls** — replaced "Show box plot" checkbox with TG-style centre-bar selector (Mean / Median / None), error-bar selector (SEM / SD / 95% CI / None), and jitter-width slider; y-axis label updated to include units (%·days)
- **Toxicity — ET plot group controls** — group order, colour, and palette-preset controls now available for all ET mode plots (Therapeutic Window, Safety-Efficacy Scatter, Total Benefit Score, Weight-Corrected TGI)
- **Toxicity — Safety-Efficacy Scatter controls** — added centroid size slider, ellipse line-width slider, ellipse line-style selector (solid / dashed / dotted), and "Shape by group" checkbox

### Fixed
- **Toxicity — Safety-Efficacy Scatter only showing reference group** — `efficacy_toxicity_bivariate()` now uses composite `Treatment⁠ID` keys when iterating per-mouse efficacy, preventing shared IDs across groups from collapsing to a single group

### Removed
- **Toxicity — Correlation Heatmap** — removed from ET plot selector, UI, server, and download handler

## [v26.05.008] 2026 May 5th

### Added
- **Toxicity Analysis — Results tabs restructured** — the single scrollable "Results" tab has been replaced with separate, focused tabs for each analysis section:
  - *Body Weight mode*: **Mixed Model** (fixed effects + estimated marginal means), **AUC** (AUC summary, nadir, pairwise comparisons), **Weight Loss KM** (time-to-threshold event data)
  - *Efficacy-Toxicity mode*: **Therapeutic Window**, **Safety-Efficacy**, **Benefit Score**, **Corrected TGI**, **Combined Report** — dynamically inserted/removed when the Statistical Model selector is changed
- **Toxicity Analysis — Post-hoc comparison type** — new sidebar dropdown in Body Weight mode: "vs Reference Group" (BH-corrected Welch t-tests restricted to comparisons involving the control) or "All Pairwise (BH-corrected)" (all group pairs); the AUC Comparisons table, Comparisons title, and AUC Comparison Forest plot all react to this selector
- **Toxicity Analysis — BH p-value correction** — Benjamini-Hochberg correction is now actually applied to the AUC comparison p-values (previously the Info tab stated BH correction was used but the displayed p-values were unadjusted); a `P_Adjusted_BH` column is now shown alongside the raw `P_Value`
- **Toxicity Analysis — AUC Comparison Forest Plot** — new plot type in the BW Plots selector: effect-size forest plot of mean AUC differences with 95% CI, coloured by BH-corrected significance; respects the comparison type selector
- **Toxicity Analysis — Diagnostics tab** — new tab showing Q-Q plot of mixed-model residuals and residuals-vs-fitted scatter, available in Body Weight mode; displays a notice in ET mode
- **Toxicity Analysis — Download All tab** — new tab with a single button that exports all available result tables as a ZIP file of CSVs
- **Toxicity Analysis — Additional download buttons** — new per-table download buttons: EMMs CSV (Mixed Model tab), Nadir CSV and Comparisons CSV (AUC tab)
- **Toxicity Analysis — Reactive tab info** — AUC tab description and comparisons title now update dynamically to reflect the current comparison type

## [v26.05.007] 2026 May 5th

### Added
- **Survival Analysis — Endpoint source always visible** — the "Use Event/Censor column" / "Derive from tumour volume threshold" radio group is now shown unconditionally in the sidebar; previously it was hidden behind a static text message when no direct Volume column was mapped

### Fixed
- **Survival — volume-threshold endpoint with calculated volume** — when the Data Upload tab is configured to calculate tumour volume from length × width measurements (rather than mapping a direct Volume column), the volume-threshold event derivation now correctly uses the calculated `Volume` column produced by `ensure_preprocessed()`; previously `rv$settings$col_volume` was empty in this case, causing a validation error and blocking the analysis

## [v26.05.006] 2026 May 5th

### Added
- **Survival Analysis — Endpoint source always visible** — the "Use Event/Censor column" / "Derive from tumour volume threshold" radio group is now shown unconditionally in the Analysis Parameters sidebar; previously it was replaced by a static text note when no Volume column was configured in Data Upload. When "Derive from tumour volume threshold" is selected and no volume column is available, an inline amber warning directs the user to Data Upload instead of silently hiding the control.

## [v26.05.005] 2026 May 5th

### Added
- **Body Weight AUC plot** — new "Show box plot" checkbox (default: on) to toggle the boxplot layer, leaving only individual data points when unchecked; new "Point size" slider (0.5–6, default: 2); dot colours now follow the Groups palette when the boxplot is hidden
- **Toxicity Analysis — Info tab** — new first tab in the results panel with structured explanations of all analyses, metrics, and parameters, organised into four sections:
  - *Overview*: when to use Body Weight Analysis vs Efficacy-Toxicity Analysis
  - *Body Weight Analysis*: mixed-effects model formula and covariate effects, REML vs ML, AUC derivation, nadir analysis, KM weight-loss curve, and threshold definition
  - *Efficacy-Toxicity Analysis*: Final TGI, Tumor AUC, and Log Cell Kill formulas; Therapeutic Window Metric; Safety-Efficacy Scatter layers; Total Benefit Score formula and λ guidance; Weight-Corrected TGI rationale; Correlation Heatmap interpretation
  - *Shared Parameters*: reference group auto-detection, tumor weight adjustment, tumor density, and baseline day

### Fixed
- **Safety-Efficacy Scatter plot** — corrected column name mismatches between `efficacy_toxicity_bivariate()` return value and the `renderPlot` code (`Max_Pct_Weight_Loss`, `Efficacy_Mean`, `Toxicity_Mean`); previously returned `[object Object]`
- **Total Benefit Score plot** — corrected column reference from `Benefit` to `Benefit_Score`; previously returned `[object Object]`
- **Weight-Corrected vs Standard TGI plot** — corrected column references from `Uncorrected_TGI`/`Corrected_TGI` to `TGI_Uncorrected`/`TGI_Corrected` (merge-suffix order from `weight_corrected_tgi()`); previously returned "arguments imply differing number of rows"
- **Correlation Heatmap plot** — corrected AUC column references from `ea$AUC`/`ta$AUC` to `ea$Efficacy_AUC`/`ta$Toxicity_AUC`; previously returned "undefined columns selected"

## [v26.05.004] 2026 May 5th

### Fixed
- **Body Weight Trajectory "Individual mice" plot** — spaghetti lines now correctly group by `paste(Treatment, ID)` composite key; previously, animals sharing the same ID across treatment groups were joined into a single crossing line
- **Body Weight Trajectory, AUC, and Weight Loss KM plots** — added group recoloring and reordering controls (drag-and-drop order list + per-group colour pickers with palette presets) matching the Tumor Growth and Survival tab functionality

## [v26.05.003] 2026 May 4th

### Added
- **Survival tab — IACUC endpoint context and volume-threshold event derivation**
  - Page header now shows subtitle "(Time-to-Endpoint)" and a warning banner explaining that in most preclinical mouse studies the recorded endpoint is time-until-tumour-volume-threshold (not true death), with a note on competing risks and the recommended pragmatic approach
  - New **Endpoint** section in the Analysis Parameters sidebar:
    - When a Tumour Volume column is configured in Data Upload, a radio group lets users choose between "Use Event/Censor column" (existing behaviour) and "Derive from tumour volume threshold" (new)
    - Volume mode: `numericInput` for the IACUC volume limit (default 2000 mm³); automatically computes per-animal Time = first day volume ≥ threshold (Event = 1) or last observation day (Event = 0) if never exceeded
  - New **Natural Death** checkbox in the sidebar (default: checked); text explains the competing-risks tradeoff and the pragmatic recommendation
  - Summary output now shows "Event Source: Derived from tumour volume threshold (>= X mm³)" when in volume mode, replacing the "Event Column:" line
  - KM plot title and axis labels update dynamically when volume mode is active ("Time-to-Tumour-Volume-Threshold by Treatment", "Time to Endpoint (Days)", "Proportion Below Volume Threshold")
- **Data Upload — updated Event/Censor Status column descriptions**
  - `help_tip` now reads: "1 = endpoint reached (tumour exceeded IACUC volume limit, or animal euthanised/died), 0 = censored"
  - Below-input text notes that volume-threshold derivation is available in the Survival sidebar
  - Help page bullet updated to mention tumour volume threshold as the typical endpoint

### Changed
- Summary tab header renamed from "Survival summary" to "Time-to-Endpoint summary"
- Cox Model "Events" bullet updated to mention tumour volume threshold; Median_Survival bullet updated to say "reached the endpoint" rather than "experienced the event"
- Risk Table N_Censored bullet updated to say "without reaching the endpoint (e.g. still within tumour volume limit at study end)"

## [v26.05.002] 2026 May 4th

### Fixed
- **Body Weight AUC and Weight Loss KM plots showing only the control group** — `body_weight_auc()` and `weight_loss_threshold()` in the `mouseExperiment` package now use a composite `Treatment + ID` key internally when grouping observations per mouse. Previously, if the same ID value (e.g. an ear tag) appeared in multiple treatment groups, all mice were incorrectly attributed to the first group encountered (typically the vehicle control). The fix ensures each treatment-group × ID combination is processed as a separate animal.

## [v26.05.001] 2026 May 4th

### Changed
- **Toxicity tab** — merged the former "Efficacy-Toxicity" tab into the Toxicity tab
  - New **Statistical Model** dropdown in Analysis Parameters sidebar (matches Tumor Growth pattern):
    - "Body Weight Analysis" — mixed-effects model, AUC, and time-to-threshold KM curve
    - "Efficacy-Toxicity Analysis" — therapeutic window metric, safety-efficacy scatter, total benefit score, weight-corrected TGI, correlation heatmap, combined report
  - Replaced four separate sub-tabs (Trajectory, Mixed Model, AUC, KM) with a unified **Plots** tab containing a plot-type selector, plus a **Results** tab with all data tables
  - Plot controls (type selector, display style, scaling, export) consolidated into a single left-hand sidebar panel within the Plots tab, consistent with Tumor Growth module layout
  - Sidebar parameters are now model-specific: Body Weight options (covariates, REML/ML, threshold, baseline day) shown only when "Body Weight Analysis" is selected; Efficacy-Toxicity options (efficacy metric, lambda, safety threshold) shown only when "Efficacy-Toxicity Analysis" is selected
- Removed standalone **Efficacy-Toxicity** top-level tab from navigation bar

## [v26.04.005] 2026 April 20th

### Added
- **Toxicity tab** — new top-level analysis tab for body weight toxicity
  - Body Weight Trajectory: spaghetti or mean±SEM/SD with shaded ribbons
  - Mixed Model: longitudinal LME with tumor-weight adjustment, covariates (volume, sex, initial mass), auto-fallback random effects
  - Body Weight AUC: trapezoidal AUC of % weight change, boxplot + jitter, nadir analysis, pairwise comparisons
  - Weight Loss Survival: KM curves for time-to-threshold weight loss, log-rank, optional Cox PH, risk table
- **Efficacy-Toxicity tab** — new top-level analysis tab combining tumor and weight data
  - Therapeutic Window Metric (TWM): TGI/MaxWeightLoss% bar chart with noise-floor handling
  - Safety-Efficacy Plot: per-mouse bivariate scatter with group ellipses/centroids, 3 efficacy metrics (TGI, tumor AUC, log cell kill)
  - Total Benefit Score: B = AUC(efficacy) − λ×AUC(toxicity) with interactive λ slider
  - Weight-Corrected TGI: side-by-side comparison after excluding mice exceeding safety threshold
  - Correlation Heatmap: efficacy vs toxicity metrics
  - Combined Report: merged summary table with CSV export
- Weight and Sex column auto-detection in data upload
- Tab enable/disable: Toxicity requires weight column, Efficacy-Toxicity requires both weight and volume

## [v26.04.004] 2026 April 18th

### Added
- "Tumor Growth Inhibition vs Dose" plot in Dose Response module
  - Per-observation %TGI computed as `(1 - Volume / mean_control_volume) * 100`
  - Dashed reference line at 0% inhibition
- Trend line fit selector (shared by both Dose Response plots): Linear or Non-linear (LOESS)
- Point size slider for Dose Response scatter plots
- Confidence interval shading toggle for Dose Response trend lines

### Changed
- Renamed Dose Response "Scatter Plot" to "Tumor Volume vs Dose"
- Removed "Boxplot by Dose" and "Mean by Dose (Bar)" plot options from Dose Response module
- Appearance controls (point colour, line colour, point size, CI toggle) now shared across both Dose Response plots

### Fixed
- Dose Response module no longer requires running Tumor Growth Analysis first — `ensure_preprocessed()` now resolves dose, ID, and cage column mappings from settings

## [v26.04.003] 2026 April 17th

### Added
- Enabled Dose Response and Drug Synergy analysis tabs (previously disabled)
- Tumor Growth Inhibition (TGI) plot in Drug Synergy module
  - Days on x-axis, percent inhibition on y-axis
  - Per-mouse TGI computed from raw volume data against control group mean
  - Combination group displayed as dashed line; monotherapies as solid
  - Configurable error bars (SEM, SD, 95% CI, or None)
- Toggleable region shading and labels on Combination Index plot (Synergistic / Additive / Antagonistic bands)
- CI plot now rebuilt dynamically to support interactive option controls
- Mobile incompatibility overlay — on screens ≤ 768 px a full-screen message recommends using a desktop, with a "Continue anyway" dismiss button
- Viewport meta tag for proper mobile scaling

### Changed
- Updated Drug Synergy analysis to use proper Loewe Additivity model (Berenbaum, 1989)
  - Expected fractional effect: min(FE_A + FE_B, 1.0)
  - Combination Index: min(FE_A + FE_B, 1) / FE_combo
- Updated synergy info box, labels, and analytics tracking to reflect Loewe terminology

### Fixed
- Corrected synergy metrics info box text to match actual package output
- Fixed Combination Index formula: numerator now uses capped Loewe expected FE — previously the uncapped sum could inflate CI and falsely indicate antagonism when individual drug effects sum to > 1.0

## [v26.04.002] 2026 April 8th 

### Added
- Comprehensive Umami analytics event tracking across all modules via a new `umami-events.js` helper library
  - `trackAnalysisStart` / `trackAnalysisComplete` for tumor growth, survival, synergy, power, and dose-response modules
  - `trackDownload` for all CSV, PNG, and ZIP download handlers
  - `trackDataUpload` and `trackDemoDataLoaded` on file/demo load
  - `trackColumnMapping` when all required column mappings are set
  - `trackReferenceGroupSelected` and `trackError` for validation paths
  - `trackTabNavigation` for main navbar; `trackFeatureAccess` when a coming-soon tab is clicked
  - `model_type_selected` and `plot_type_selected` events for tumor growth model/plot switching
  - `session_saved`, `session_restored`, and `analysis_log` download events
  - `dark_mode_toggled` event fired directly from `darkmode.js`

### Fixed
- Loading new data now resets all downstream analysis results — tumor growth, survival, synergy, power, and dose-response results are cleared whenever a new file or demo dataset is loaded, preventing stale results from a previous dataset from remaining visible
