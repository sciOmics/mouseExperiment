mouseExperiment:

2026-10-26

- [] Get rid of mRECIST response analysis
- [] Expand growth delay to show statistics in summary table and an individual mouse table
- [] Do we need a menu for the different Stan backends? What is different between the two?
- [] After running a Bayesian analysis, the "ANOVA" tab should not be visible. 
- [] The posterior check tab from a Bayesian analysis should be merged into the "Diagnostics" tab.
- [] Need to provide users info on what "Influential observations (Pareto-k > 0.7)" means and how to interpret. Similarly for max-treedepth hit(s), E-BFMI, Bayes R², PPC coverage, PSIS-LOO, etc. 
- [] Users need info on what the posterior summary metrics (e.g. rhat, Estimate, ESS) mean and how to intepret
- [] The "Growth Rates" tab has a "Per-mouse drill-down" that graphs the growth rate for an individual mouse. Add functionality to select multiple mice. 
- [] The "Extrapolation Points (plots only) " control on Analysis parameters side menu should be moved to the plotting menu under "Show extrapolated curves/points"
- [] Remove the " The first fit of a model compiles it (about 20 s); later fits with the same structure reuse it. Fits run one at a time in the background: you can close the browser and resume or cancel from the banner when you return." text. This is not useful to users. 
- [] Add AFT analysis? 
- [] The tumor Growth Bayes analysis doesn't appear to have a prior vs posterior plot. Perhaps the posterior (with or without prior as an option) plots should be included in the plotting types options. Include the ability to shade for HDI (defined by user: e.g. 99%)
- [] The arrangement of the tabs in Survival Analysis and Tumor Growth are different. We need to standardize. For example, "Download Results" should always be the last tab. 
- [] When a user downloads results, how does this work for plots? Ideally the users adjustments to plots in the plots tab should be maintained and downloaded as such. 
- [] When using "Download Results" the user should be able to download the results for multiple analyses they have conducted (e.g. Tumor Growth LME4, Tumor Growth Bayes, Survival Bayes, etc.)
- [] The tooltip in Drug Synergy Endpoint estimand gives numbers (e.g. "underestimated TGI by about 23 points"). Remove specifics like this.
- [] The synergy plot doesn't show in the legend what colors are for the shading of Synergy vs. Anatagonism
- [] In the toxicity analysis parameters, the Covariates has "Tumor Volume" checked by default, but it indicates that tumor mass is already subtracted from body weight (also by default). 
- [] Given we don't have ways built in to compare models, do we need the option for "ML" in the estimation method for Toxicity? 
- [] Remove any language in the model fitting progress box overlays that talks about only fitting a Bayesian model once. This is documentation info, not needed on the user interface. 
- [] In the info tabs, it might be useful to have info on fixed vs random effects
- [] Is it appropriate to add p-values to the rate of weight change summary in the Model Results tab of Toxicity analysis?
- [] In some places lower.CL and upper.CL is reported without indicating the percentage (e.g. 95% CL)
- [] The demo weight dataset doesn't appear to have any mice/groups that drop below the 20% threshold. We should modify these synthetic data so that the combo group has severe toxicity and the individual drug groups have mild. 
- [] Should we add the ability to assess drug synergy in regards to toxicity? 
- [] The dose response analysis: does this assume a linear relationship? Can it capture others (e.g. polynomial) or describe the most parsimonous fit? The plotting function allows Linear and LOESS trend line fits. Perhaps we should include other common fits? 
- [] Do we want to move the "staging" environment so that it is no longer on a staging subdomain and is instead run local? This would free resources on the VPS. If so, we will want to change the name of "staging" to "development"







