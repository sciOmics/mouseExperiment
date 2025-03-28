# Plan to Fix Mouse Identification in tumor_growth_statistics

## Problem Description

The `tumor_growth_statistics` function is not properly accounting for all mice in each treatment group when mice have the same ID but are in different cages. The issue specifically appears when using the `model_type="auc"` option.

The current function creates a composite ID using only the ID and Treatment columns:
```R
composite_id <- paste(auc_df[[id_column]], auc_df[[treatment_column]], sep = "_")
```

This approach cannot distinguish between mice that have the same ID but are in different cages within the same treatment group. While our test with the synthetic dataset showed all mice were properly identified, the user's real-world data likely has mice with duplicate IDs across different cages.

## Proposed Solution

We will modify the `composite_id` creation in the `tumor_growth_statistics` function to include the cage column, ensuring that mice with the same ID in different cages are uniquely identified:

```R
composite_id <- paste(auc_df[[id_column]], auc_df[[treatment_column]], auc_df[[cage_column]], sep = "_")
```

This change will maintain backward compatibility while resolving the issue of identifying all unique mice correctly.

## Implementation Steps

1. Edit the `R/tumor_growth_statistics.R` file at approximately line 491.

2. Locate the following code in the AUC calculation section:
   ```R
   # Calculate AUC for each individual (Composite ID ensures unique subject-treatment combinations)
   composite_id <- paste(auc_df[[id_column]], auc_df[[treatment_column]], sep = "_")
   auc_data <- data.frame()
   ```

3. Replace it with:
   ```R
   # Calculate AUC for each individual using ID + Treatment + Cage for unique identification
   composite_id <- paste(auc_df[[id_column]], auc_df[[treatment_column]], auc_df[[cage_column]], sep = "_")
   auc_data <- data.frame()
   ```

4. Update the comment above the code to clearly document the change and reasoning.

5. Create a test script that specifically validates this fix using a dataset that contains mice with the same ID but in different cages.

6. Update the CHANGELOG.md to document this fix.

## Testing Plan

1. Create a synthetic dataset that explicitly reproduces the issue:
   - Multiple mice with the same ID in different cages within the same treatment group

2. Test the function before and after the fix to verify:
   - The original function fails to identify all unique mice
   - The fixed function correctly identifies all unique mice

3. Verify the fix doesn't break existing functionality for datasets that don't have this specific issue.

4. Ensure the modified function produces correct AUC values for all mice.

## Documentation Updates

1. Update the function documentation to clarify how unique mice are identified.

2. Add a note to the CHANGELOG.md about this fix, explaining the issue and solution.

3. Consider updating examples to illustrate the proper structure for datasets with multiple cages.

## Expected Results

After implementing this fix, the `tumor_growth_statistics` function will correctly:

1. Identify all unique mice in each treatment group, even when mice share the same ID across different cages
2. Calculate AUC values for all mice
3. Generate statistical analyses based on the complete dataset

This will ensure that users get accurate results from their experiments, especially in cases where mice are housed in multiple cages per treatment group. 