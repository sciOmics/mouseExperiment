# Plan to Fix Mouse Identification in Additional Functions

## Problem Description

After fixing the `tumor_growth_statistics` function to properly account for mice with the same ID in different cages, we've identified two additional functions that have the same issue:

1. `tumor_auc_analysis` in R/tumor_auc_analysis.R
2. `calculate_auc_values` in R/post_power_analysis.R

Both functions process subjects based solely on their ID value without considering cage information. This means that if there are mice with the same ID in different cages (but in the same treatment group), they would be incorrectly treated as a single mouse. This could lead to inaccurate AUC calculations and statistical analyses.

## Proposed Solution

We need to modify both functions to create unique subject identifiers that include cage information, similar to the fix applied to `tumor_growth_statistics`. 

### For `tumor_auc_analysis` function:

1. Modify the function to create a composite subject identifier using ID, Treatment, and Cage values.
2. Use this composite identifier throughout the function instead of just the ID value.

### For `calculate_auc_values` function:

1. Add a `cage_column` parameter with a default value of "Cage".
2. Modify the function to create a composite subject identifier using ID, Treatment, and Cage values.
3. Use this composite identifier for AUC calculations.

## Implementation Steps

### For `tumor_auc_analysis`:

1. Add a `cage_column` parameter with a default value of "Cage".
2. Check for the existence of the cage column at the start of the function.
3. Create a composite subject identifier combining ID, Treatment, and Cage.
4. Use this composite identifier for subject-level operations.
5. Preserve the original ID in the output for reference.

### For `calculate_auc_values`:

1. Add a `cage_column` parameter with a default value of "Cage".
2. Check if the cage column exists in the data.
3. Create a composite subject identifier including cage information when available.
4. Use this composite identifier for AUC calculations.
5. Return the original ID in the output for reference.

## Testing Plan

1. Create a test script that explicitly demonstrates the issue with duplicate IDs across cages.
2. Verify that the original functions fail to correctly identify all unique mice.
3. Test the modified functions to confirm they correctly identify all unique mice.
4. Ensure the modifications do not break existing functionality for datasets without duplicate IDs.

## Documentation Updates

1. Update function documentation to explain the use of the composite ID for unique subject identification.
2. Add a note to the CHANGELOG.md about these fixes.

## Expected Results

After implementing this fix, the `tumor_growth_statistics` function will correctly:

1. Identify all unique mice in each treatment group, even when mice share the same ID across different cages
2. Calculate AUC values for all mice
3. Generate statistical analyses based on the complete dataset

This will ensure that users get accurate results from their experiments, especially in cases where mice are housed in multiple cages per treatment group. 