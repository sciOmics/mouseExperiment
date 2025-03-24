# Test script for date parsing with different formats
# This script tests the fixes to the calculate_dates function

# Load the package
library(mouseExperiment)

# Create a simple test data frame
test_data <- data.frame(
  Date = c("17-Mar", "18-Mar", "19-Mar", "20-Mar", "21-Mar"),
  Cage = 1,
  ID = c("A", "B", "C", "D", "E"),
  Treatment = "Test",
  Length = c(1, 1.5, 2, 2.5, 3),
  Width = c(0.5, 0.7, 1, 1.2, 1.5)
)

# Print the test data
cat("Test data:\n")
print(test_data)

# Try with just the start date (should auto-detect and use current year)
cat("\nTest 1: With just start date '17-Mar':\n")
result1 <- tryCatch({
  calculate_dates(test_data, start_date = "17-Mar")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(result1)) {
  cat("SUCCESS: Parsed correctly\n")
  print(head(result1))
}

# Try with explicit date format and year
cat("\nTest 2: With explicit date format and year:\n")
result2 <- tryCatch({
  calculate_dates(test_data, start_date = "17-Mar", 
                 date_format = "%d-%b", year = 2023)
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(result2)) {
  cat("SUCCESS: Parsed correctly\n")
  print(head(result2))
}

# Try with MM/DD/YYYY format
cat("\nTest 3: With MM/DD/YYYY format:\n")
mm_dd_yyyy_data <- data.frame(
  Date = c("03/17/2023", "03/18/2023", "03/19/2023", "03/20/2023", "03/21/2023"),
  Cage = 1,
  ID = c("A", "B", "C", "D", "E"),
  Treatment = "Test",
  Length = c(1, 1.5, 2, 2.5, 3),
  Width = c(0.5, 0.7, 1, 1.2, 1.5)
)

result3 <- tryCatch({
  calculate_dates(mm_dd_yyyy_data, start_date = "03/17/2023")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(result3)) {
  cat("SUCCESS: Parsed correctly\n")
  print(head(result3))
}

cat("\nAll tests completed\n")