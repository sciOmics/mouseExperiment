# Dashboard

## Notes

## Tasks to be Done

### Features:
- Add "Methods" section to documentation
- Change coloring of boxes, selections, etc. 
- Change loading icon to a mouse scratching
- Add ability to change "units" on plots. For example, the Growth Rate Distribution plot shows log transformed units. Add ability to change to original units. 
- Add a "Learning" section to the Help that gives info on the who, what, when, where and how of the analysis. Best practices in experiment design and analysis. Etc. 
    - Random intercept and analyzing the growth rate rather than comparing the volumes on a particular day controls for differences in when initial tumors develop (resulting in part to variance in tumor injection)
- Add ability in graph to add marker to day that treatment began. 
- Add "Height" column and Ellipsoid 3D volume calculation for tumor measurement
 
### Bugs:
- The AUC analysis in Tumor Growth Analysis: the Comparisons tab info box is showing info for columns that are not shown in the table (t.ratio, p.value, Fold_Change). 
- The AUC analysis Treatment Effects info box does not match the table columns at all. 

### Testing:
- Compare results against manual analysis
- Code review

### Other:
- Update documentation. 
- Update DESCRIPTION file.

## Future Improvements
- Add CHANGELOG.md (once initial version is published)
- Versioning
