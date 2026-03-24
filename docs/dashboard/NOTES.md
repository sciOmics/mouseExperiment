# Dashboard

## Tasks in Progress

- Polynomial transformation is missing from Tumor Growth Analysis. 
- Add option for plots determining what the horizontal bar shows: "none", "mean", "median", "geometric_mean"
- Controls for text size on plots
- Add option for "Legend Position": "None" "Above" "Below", "Left", "Right"...
- Add option for plot dimensions.
- AUC Comparison plot should be beeswarm plot. 
- Treatment schedule plots are not working. 
- The legend/lines/coloring in the Trend Over Time plot is mixing up groups. The dotted line that should be the Bliss Expected is being labeled "HDACi", while the line that is the combination treatment (HDACi + PD1) is being labeled Bliss Expected. The line that is one durg is being labeled as the combo (HDACi + PD1). 
- The Trend Over Time plot should always have the y-axis scaled from 0 to 100. 
- In the Trend Over Time plot, the Synergy and Antagonism legend needs to be elsewhere. It is currently ontop of the plot.
- Synergy Plot and Trend Over Time Plot are the same. Eliminate the Synergy Plot and rename Trend over Time to Synergy Plot
- Add y-axis scaling options (linear, log2, log10) on Tumor Growth Curves. 
- Extraploated points on the Tumor Growth Curves should be hollow circles. Include a legend with this info. 
- In the Drug Synergy Analysis, the Synergy Metrics and TGI Table have the same information. Eliminate the TGI Table.

## Tasks to be Done

- Compare results against manual analysis
- Fix "explanations" in the analysis pages
- Add "Methods" section
- Have each analysis section have it's own color scheme:
    - Tumor Growth:
    - Survival:
    - Dose Response:
    - Drug Synergy: 
    - Power Analysis: 
- Make the left side menus on the Analyses pages collapsable to give more horizontal room for the results/plots. 
- Smart detect/recommend Statistical Model and Transformation? 
- Add "in plot - top left" & "in plot - top right" options to legend locations.
- Remove "log2" option for Y-axis scale. 
- Add back the Antagonism/Synergy shading to the Synergy Plot


## Completed Tasks
- Merge mouseExperiment and mouseExperimentDemo (vignette) 
- Make plots interactive/customizable. Instead of showing all plots at once, have it selectable what plot to show. 
- Add mouse outline png as background.
- Generate demo data sets and code to use for automated testing of dashboard and statistics/analysis/plots. In particular, we need to ensure that the statistics being reported are correct. 
- For reference group choosing, the default should look for a group that contains similar wording to "vehicle" "control" "IgG" or "DMSO". 
- The forest plot x-axis scaling should be log(2) based and balanced such that whatever the max value one direction from the value of "0", the other direction should have the same. For example, if the plot min is -2 as determined by the data, the max should be 2. Also, do not color by significance. 
- Post-hoc Power Analysis is not working. "Error: At least two treatment groups are required for power analysis"
- Need to incorporate dose response analysis. Will need to have an optional "Dose" column in the data upload. 
- Drug Synergy Analysis is not working. 
- Dose levels not being understood by treatment summary or analysis. The different doses are being lumped together since they all share the same "Treatment" value. 
- Survival Analysis checkbox for "Show Risk Table" is pointless. Remove it.
- Use the same criteria for plotting the Effect Size Forest Plot as we did with the hazard ratio Forest Plot (i.e balancing the x-axis)
- The AUC Comparison Plot has numbers for the group labels instead of the actual treatment (& dose if applicable).
- The QQ Plot should be removed from the Plots tab since it now resides in Diagnostics. 
- Rework the main/top menu tabs.Tab Bar: have mouseExperiment left aligned, other tabs are centered except "Help" and "Citation" Break Analysis subtabs into invididual tabs: Tabs: Data, Tumor Growth, Survival, Dose Response, Drug Synergy, Power Analysis, Help, Citation. Tabs should be greyed out if not applicable. For example, all will be greyed out and unclickable without data uploaded. If the data does not contain survival data, the survival tab will be greyed out. If no dose column is present/specified, the Dose Resposne tab is inactive. Etc. Help should have to links, one for the GitHub repo for users to get documentation/other info and the other for the repo's Issues page for users to submit bugs or feature requests. "https://github.com/sciOmics/mouseExperiment" "https://github.com/sciOmics/mouseExperiment/issues"
- Include access to the demo data in mouseExperiment for users to use as demo data in the dashboard. 

## Notes