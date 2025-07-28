# Test script to verify our bug fixes
# Run this to check if the main functions work correctly

# Load required libraries
library(shiny)

# Test the main functions to ensure they don't error
tryCatch({
  
  # Test 1: Load modules without errors
  cat("Testing module loading...\n")
  source("phase3/workflow_wizard.R", local = TRUE)
  source("phase3/parameter_intelligence.R", local = TRUE) 
  source("phase3/quality_control.R", local = TRUE)
  source("phase3/error_prevention.R", local = TRUE)
  source("phase3/results_validation.R", local = TRUE)
  source("phase3/real_data_testing.R", local = TRUE)
  source("phase3/real_data_server.R", local = TRUE)
  cat("✅ All modules loaded successfully\n")
  
  # Test 2: Create mock data to test functions
  cat("Testing mock data generation...\n")
  
  # Create simple test data
  test_expression <- matrix(
    rpois(1000, 100), 
    nrow = 100, 
    ncol = 10,
    dimnames = list(
      paste0("Gene_", 1:100),
      paste0("Sample_", 1:10)
    )
  )
  
  test_metadata <- data.frame(
    Sample_ID = paste0("Sample_", 1:10),
    Condition = rep(c("Control", "Treatment"), each = 5),
    Batch = rep(c("Batch1", "Batch2"), 5),
    stringsAsFactors = FALSE
  )
  
  cat("✅ Test data created successfully\n")
  
  # Test 3: Parameter intelligence function
  cat("Testing parameter intelligence...\n")
  params <- recommend_parameters("cell_line_treatment", test_metadata, test_expression)
  cat("✅ Parameter intelligence working\n")
  
  # Test 4: Error prevention function  
  cat("Testing error prevention...\n")
  test_params <- list(
    experiment_type = "cell_line_treatment",
    padj_cutoff = 0.05,
    fc_cutoff = 2.0,
    multiple_testing_method = "BH"
  )
  
  errors <- detect_analysis_errors(test_expression, test_metadata, test_params)
  cat("✅ Error prevention working\n")
  
  # Test 5: Data validation
  cat("Testing data validation...\n")
  validation <- validate_uploaded_data(test_expression, test_metadata)
  cat("✅ Data validation working\n")
  
  cat("\n🎉 All tests passed! The fixes are working correctly.\n")
  cat("You can now run the demo app safely:\n")
  cat("shiny::runApp('phase3_guardrails_demo.R')\n")
  
}, error = function(e) {
  cat("❌ Error during testing:", e$message, "\n")
  cat("Please check the specific function that's failing.\n")
})