# Dashboard

## Tasks in Progress



## Tasks to be Done

- Make plots interactive/customizable. Instead of showing all plots at once, have it selectable what plot to show. 
- Add mouse outline png as background (lower right side of page)

## Complete Tasks

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