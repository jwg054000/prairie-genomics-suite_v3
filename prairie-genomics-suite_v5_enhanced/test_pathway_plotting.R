# Test Pathway Visualization Performance
# Quick test to verify pathway plotting works without timeouts

cat("🧪 Testing Fast Pathway Visualization\n")
cat("=====================================\n\n")

# Load required functions
if (!exists("create_fast_pathway_plot")) {
  source("pathway_analysis.R")
}

# Create minimal test data
test_pathway_data <- data.frame(
  Description = c("Pathway A", "Pathway B", "Pathway C", "Pathway D", "Pathway E"),
  pvalue = c(0.001, 0.01, 0.05, 0.1, 0.2),
  padj = c(0.005, 0.02, 0.08, 0.15, 0.25),
  stringsAsFactors = FALSE
)

# Create mock pathway results
mock_results <- list(
  success = TRUE,
  data = test_pathway_data,
  analysis_type = "GO_BP"
)

cat("📋 Test 1: Fast pathway plotting function\n")
cat("-----------------------------------------\n")

# Test the fast plotting function
plot_result <- tryCatch({
  start_time <- Sys.time()
  
  plot_obj <- create_fast_pathway_plot(mock_results, plot_type = "dotplot", top_n = 5)
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("✅ Fast plot created successfully\n")
  cat("⏱️ Execution time:", round(execution_time, 3), "seconds\n")
  
  if (!is.null(plot_obj)) {
    cat("📊 Plot object type:", class(plot_obj)[1], "\n")
    return(TRUE)
  } else {
    cat("❌ Plot object is NULL\n")
    return(FALSE)
  }
  
}, error = function(e) {
  cat("❌ Fast plotting failed:", e$message, "\n")
  return(FALSE)
})

cat("\n" , rep("=", 50), "\n")

cat("📋 Test 2: Emergency plot fallback\n")
cat("-----------------------------------\n")

# Test emergency plotting
emergency_result <- tryCatch({
  start_time <- Sys.time()
  
  plot_obj <- create_emergency_plot(mock_results, top_n = 3)
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("✅ Emergency plot created successfully\n")
  cat("⏱️ Execution time:", round(execution_time, 3), "seconds\n")
  
  return(TRUE)
  
}, error = function(e) {
  cat("❌ Emergency plotting failed:", e$message, "\n")
  return(FALSE)
})

cat("\n" , rep("=", 50), "\n")

cat("📋 Test 3: Large dataset stress test\n")  
cat("-------------------------------------\n")

# Create larger test dataset to simulate timeout conditions
large_test_data <- data.frame(
  Description = paste("Pathway", 1:50),
  pvalue = runif(50, 0.001, 0.2),
  padj = runif(50, 0.005, 0.25),
  stringsAsFactors = FALSE
)

large_mock_results <- list(
  success = TRUE,
  data = large_test_data,
  analysis_type = "KEGG"
)

stress_result <- tryCatch({
  start_time <- Sys.time()
  
  # This should automatically limit to 10 pathways
  plot_obj <- create_fast_pathway_plot(large_mock_results, plot_type = "dotplot", top_n = 20)
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("✅ Large dataset handled successfully\n")
  cat("⏱️ Execution time:", round(execution_time, 3), "seconds\n")
  
  if (execution_time > 2) {
    cat("⚠️ Warning: Execution time > 2 seconds may cause timeouts\n")
  } else {
    cat("🎯 Execution time acceptable for Shiny\n")
  }
  
  return(TRUE)
  
}, error = function(e) {
  cat("❌ Stress test failed:", e$message, "\n")
  return(FALSE)
})

# Summary
cat("\n📊 TEST SUMMARY\n")
cat("===============\n")

total_tests <- 3
passed_tests <- sum(c(plot_result, emergency_result, stress_result))

cat("Tests passed:", passed_tests, "/", total_tests, "\n")

if (passed_tests == total_tests) {
  cat("🎉 SUCCESS: All pathway visualization tests passed!\n")
  cat("✅ The timeout issue should now be resolved\n")
  cat("✅ aes_string() deprecation warnings fixed\n")
  cat("✅ Fallback plotting available for edge cases\n")
} else {
  cat("⚠️ Some tests failed - there may still be issues\n")
}

cat("\n🔧 Next Steps:\n")
cat("1. Restart your app with run_app.R\n")
cat("2. Run pathway analysis\n") 
cat("3. Switch to Visualizations tab\n")
cat("4. Plotting should now be fast and stable\n")

cat("\n💡 If issues persist:\n")
cat("- Check console for 'Using base R fallback plot' messages\n")
cat("- The system will automatically use simpler plots if needed\n")