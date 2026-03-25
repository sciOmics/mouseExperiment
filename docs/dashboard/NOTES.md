# Dashboard

## Notes

## Tasks to be Done

### Features:
- Improve "explanations" in the analysis pages
- Add "Methods" section
- Change coloring of boxes, selections, etc. 
- Fix mouse icon (fill with white)
- Change loading icon to a mouse scratching
- Revise the Citation page

 
### Bugs:
- Remove the "Plot Export Settings" from the Analysis Parameters side menu. 
- Add to the plot controls the ability to resize the graph displayed (e.g. height = 6, width = 8) 
- The Effect Size Forest Plot is showing comparisons between groups despite having run Reference Group = Control. 
- When "Show extrapolated curves/points" is selected, the extrapolated points should not be filled, whereas non-extrapolated points are filled. There is a note below the graph indicating this is done, but it is not. 
- The Tumor Growth Curves plot is getting cut off on the right hand side of the x-axis. 

### Testing:
- Compare results against manual analysis
- Code review

## Future Improvements
- Add CHANGELOG.md
- Versioning

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
- Add "in plot - top left" & "in plot - top right" options to legend locations.
- Remove "log2" option for Y-axis scale. 
- Add back the Antagonism/Synergy shading to the Synergy Plot; have as controllable (show/do not show) from the side menu
- Add controls for the text size and legend location to the Synergy Plot. 
- The Synergy Plot legend drug names should use the real drug names selected in the data menu (e.g. Drug A (single agent) is aPD1; use aPD1 as name)
- Add the "X-axis start (Day)" feature from the KM plots to the Tuomr Growth Plots.
- Add controls for what to show on the HR Forest Plot: HR, 95% CI, p-value. 
- The x-axis labels on the HR Forest Plot are crowded. Try switching to a log10 scale. 
- Make the left side menus on the Analyses pages collapsable to give more horizontal room for the results/plots. 
- Add the control group to the forest plots to help orient users that the control group exists and the values are relative to it. 
- In the Tumor Growth Analysis section, the Analysis Parameters menu needs to be reordered. It should be Statistical Model, Reference Group, Post-hoc Comparisons, Data Transformation, Extrapolation Points, Cage Effect Handling, Random Effects. For now, we are going to not show the Polynomial Degree. Also do not have the "Advanced Options" text and drop down. Everything will be shown. Finally, add a brief description below each menu giving information (similar to what is done currently with Extrapolation Points). 
- In the growth rates tabe of the Tumor Growth Analysis, the growth_rate shown is presumably the log transformed value. Please include a "back-transformed" value in the original units as well. This applies to anywhere where the transformed values are shown.
- For plots, add the ability to change the height, width, and resolution.
- Dark Mode: Have clicking the mouse logo beside "mouseExperiment" control dark mode. Also, when dark mode is enabled, the mouse png background images need to be inverted. 
- In the Tumor Growth Curves plots, add the ability to change the size of the points. Also, change the points from squares to circles. 
- The group order adjustment on plots does not change the order of the groups on the plot. It should. 
- The AUC Comparison plot on Tumor Growth Analysis is not plotting correctly. All the individual dots (samples) are at the value 0. There are group mean dots at diffrent values with vertical error bars, but the horizontal ends of the error bars are missing. 
- The Effect Size Forest Plot is wrong. The values plotted for groups are not correct. For example, it is showing HDACi vs Control as above Control and HDACi + PD1 vs control as below. Something is wrong.
- Treatment schedule plotting is not working. Let's remove all Treatment Schedule plotting for now. 
Let's amend the text explanations/info for the Tumor Growth Analysis tabs:
In the ANOVA tab, give info on what the Chisq, Df, and Pr(>Chisq) columns represent and how to interpret.
In the Comparisons tab, give info (what they represent and how to interpret) on all the columns (estimate, SE, df, t.ratio, p.value)
In the Treatment Effects tab, give info all the columns.
Same for Growth Rates tab. 
- The plot sizing controls should control the plot shown too. Move all the plot sizing controls and the download plot button to the plotting control panel.
- Similar to the Tumor Growth tabs, have info boxes show what columns represent and how to interpret for all the tabs in all the other analyses (Survival, Dose Response, Drug Synergy, Power Analysis)
- Add option to KM plot on Survival Analysis for showing/not showing the risk table on the plot. 
- The Effect Size Forest Plot is not showing anything. Also the Forest plot for Survival Analysis is not working. Error: dim(X) must have a positive length. 
- The group order is still not changing the order of the groups on the graphs.
- Add statement of no liability
- Cannot return to other tabs (e.g. Data) by clicking once on the Help or Citation tabs. 
- Need to indicate for survival event notion, 1 is an event (e.g. death) and 0 is censored. Add this info to the Help & Documentation Column Naming section. Also add it to the Event/Censor Status Column subtext on the Data Upload & Configuration.



 mouseExperiment is provided "as is" for research purposes only, without warranty of any kind, express or implied. The authors and contributors accept no liability for errors, inaccuracies, or any consequences — including but not limited to incorrect results, data loss, or decisions made on the basis of outputs from this software.

All results must be independently validated before use in publications, reports, regulatory submissions, clinical decisions, or any other application. Use of this tool is entirely at your own risk.

By clicking I Understand , you acknowledge that you have read this disclaimer and accept full responsibility for validating any results produced by this software. 