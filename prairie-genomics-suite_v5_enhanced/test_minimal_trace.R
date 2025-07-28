# Minimal test with enhanced error tracing
# Temporarily override the stop function to get better error location
original_stop <- base::stop
error_locations <- list()

traced_stop <- function(...) {
  # Capture the call stack
  calls <- sys.calls()
  error_locations[[length(error_locations) + 1]] <<- calls
  
  # Get the error message
  error_msg <- paste0(...)
  
  cat("\n🔥 TRACED ERROR:", error_msg, "\n")
  cat("📍 Error occurred at depth:", length(calls), "\n")
  
  # Show relevant calls
  if (length(calls) > 0) {
    cat("📍 Recent calls:\n")
    start_idx <- max(1, length(calls) - 5)
    for (i in start_idx:length(calls)) {
      cat("  ", i, ":", deparse(calls[[i]])[1], "\n")
    }
  }
  
  # Call original stop
  original_stop(...)
}

# Override stop temporarily
assignInNamespace("stop", traced_stop, ns = "base")

# Now run the test
tryCatch({
  source('pathway_analysis.R')
  
  test_data <- data.frame(
    Gene = paste0('Gene', 1:5),
    baseMean = rep(1000, 5),
    log2FoldChange = c(-2, 2, -3, 3, -1.5),
    lfcSE = rep(0.3, 5),
    pvalue = rep(1e-5, 5),
    padj = rep(1e-4, 5)
  )
  
  cat('\n🧪 Running test with error tracing...\n')
  result <- run_pathway_analysis(
    deseq2_results = test_data,
    analysis_type = 'GO',
    species = 'mouse'
  )
  
  if (!result$success) {
    cat('\nAnalysis failed with error:', result$error, '\n')
  }
  
}, finally = {
  # Restore original stop function
  assignInNamespace("stop", original_stop, ns = "base")
  cat('\n✅ Error tracing complete\n')
})