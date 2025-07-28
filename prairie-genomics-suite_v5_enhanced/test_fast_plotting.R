# Test Ultra-Fast Plotting System for Pathway Visualization Issues
# 
# This script tests the new ultra-fast plotting system designed to resolve
# visualization issues where plots get stuck with large datasets (538+ pathways)
#
# Author: Prairie Genomics Team
# Date: January 27, 2025

cat("⚡ TESTING ULTRA-FAST PLOTTING SYSTEM\n")
cat("=====================================\n\n")

# Load required modules
source("pathway_analysis.R")

# Test fast plotting with mock GSEA results (similar to user's 538 pathways)
test_fast_gsea_plotting <- function() {
  cat("📊 TEST 1: FAST GSEA PLOTTING (Large Dataset)\n")
  cat("=============================================\n")
  
  # Create mock GSEA results similar to user's data (538 pathways)
  set.seed(123)
  n_pathways <- 538  # Same size as user's stuck analysis
  
  mock_gsea_data <- data.frame(
    pathway = paste0("PATHWAY_", 1:n_pathways),
    pval = runif(n_pathways, 0, 0.1),
    padj = runif(n_pathways, 0, 0.2),
    NES = rnorm(n_pathways, 0, 2),
    size = sample(20:300, n_pathways, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Sort by significance like real GSEA results
  mock_gsea_data <- mock_gsea_data[order(mock_gsea_data$padj), ]
  
  # Create mock pathway result object
  mock_gsea_result <- list(
    success = TRUE,
    data = mock_gsea_data,
    analysis_type = "GSEA",
    method = "fgseaMultilevel",
    species = "mouse",
    n_pathways = n_pathways,
    n_significant = sum(mock_gsea_data$padj < 0.05)
  )
  
  cat("✅ Mock GSEA data created:\n")
  cat("  - Total pathways:", n_pathways, "\n")
  cat("  - Significant pathways:", mock_gsea_result$n_significant, "\n")
  
  # Test different plot types with timing
  plot_types <- c("dotplot", "gsea", "enrichment")
  plot_results <- list()
  
  for (plot_type in plot_types) {
    cat(paste0("\n⚡ Testing ", plot_type, " plot...\n"))
    
    start_time <- Sys.time()
    
    plot_obj <- tryCatch({
      create_fast_pathway_plot(mock_gsea_result, plot_type = plot_type, top_n = 20)
    }, error = function(e) {
      cat("  ❌ Error:", e$message, "\n")
      NULL
    })
    
    end_time <- Sys.time()
    duration <- as.numeric(end_time - start_time, units = "secs")
    
    if (!is.null(plot_obj)) {
      cat(paste0("  ✅ ", plot_type, " plot created in ", round(duration, 2), " seconds\n"))
      plot_results[[plot_type]] <- list(success = TRUE, duration = duration)
    } else {
      cat(paste0("  ❌ ", plot_type, " plot failed in ", round(duration, 2), " seconds\n"))
      plot_results[[plot_type]] <- list(success = FALSE, duration = duration)
    }
  }
  
  return(plot_results)
}

# Test fast plotting with mock ORA results (GO/KEGG/MSigDB/Reactome)
test_fast_ora_plotting <- function() {
  cat("\n📊 TEST 2: FAST ORA PLOTTING (GO/KEGG/MSigDB/Reactome)\n")
  cat("========================================================\n")
  
  analysis_types <- c("GO", "KEGG", "MSigDB", "Reactome")
  all_results <- list()
  
  for (analysis_type in analysis_types) {
    cat(paste0("\n🔄 Testing ", analysis_type, " plotting...\n"))
    
    # Create mock ORA results
    set.seed(456)
    n_pathways <- sample(50:200, 1)  # Variable pathway count
    
    mock_ora_data <- data.frame(
      ID = paste0(analysis_type, "_", sprintf("%04d", 1:n_pathways)),
      Description = paste0(analysis_type, " Pathway ", 1:n_pathways),
      pvalue = runif(n_pathways, 0, 0.1),
      p.adjust = runif(n_pathways, 0, 0.2),
      Count = sample(5:50, n_pathways, replace = TRUE),
      stringsAsFactors = FALSE
    )
    
    # Sort by significance
    mock_ora_data <- mock_ora_data[order(mock_ora_data$p.adjust), ]
    
    mock_ora_result <- list(
      success = TRUE,
      data = mock_ora_data,
      analysis_type = analysis_type,
      species = "mouse",
      n_pathways = n_pathways
    )
    
    # Test plotting with timing
    start_time <- Sys.time()
    
    plot_obj <- tryCatch({
      create_fast_pathway_plot(mock_ora_result, plot_type = "dotplot", top_n = 20)
    }, error = function(e) {
      cat("  ❌ Error:", e$message, "\n")
      NULL
    })
    
    end_time <- Sys.time()
    duration <- as.numeric(end_time - start_time, units = "secs")
    
    if (!is.null(plot_obj)) {
      cat(paste0("  ✅ ", analysis_type, " plot created in ", round(duration, 2), " seconds\n"))
      all_results[[analysis_type]] <- list(success = TRUE, duration = duration, pathways = n_pathways)
    } else {
      cat(paste0("  ❌ ", analysis_type, " plot failed in ", round(duration, 2), " seconds\n"))
      all_results[[analysis_type]] <- list(success = FALSE, duration = duration, pathways = n_pathways)
    }
  }
  
  return(all_results)
}

# Test emergency fallback system
test_emergency_plotting <- function() {
  cat("\n🚨 TEST 3: EMERGENCY FALLBACK PLOTTING\n")
  cat("======================================\n")
  
  # Create problematic mock data that might cause issues
  mock_problematic_result <- list(
    success = TRUE,
    data = data.frame(
      pathway = c("BROKEN_PATHWAY_1", "BROKEN_PATHWAY_2"),
      weird_column = c(NA, Inf),
      another_col = c("text", "more_text")
    ),
    analysis_type = "PROBLEMATIC",
    n_pathways = 2
  )
  
  cat("🔄 Testing emergency fallback with problematic data...\n")
  
  start_time <- Sys.time()
  
  emergency_plot <- tryCatch({
    create_emergency_plot(mock_problematic_result, top_n = 20)
  }, error = function(e) {
    cat("  ❌ Even emergency failed:", e$message, "\n")
    NULL
  })
  
  end_time <- Sys.time()
  duration <- as.numeric(end_time - start_time, units = "secs")
  
  if (!is.null(emergency_plot)) {
    cat(paste0("  ✅ Emergency plot created in ", round(duration, 2), " seconds\n"))
    return(TRUE)
  } else {
    cat(paste0("  ❌ Emergency plot failed in ", round(duration, 2), " seconds\n"))
    return(FALSE)
  }
}

# Run all tests
cat("🚀 STARTING ULTRA-FAST PLOTTING TESTS\n")
cat("======================================\n")

# Test 1: GSEA plotting (large datasets like user's 538 pathways)
gsea_results <- test_fast_gsea_plotting()

# Test 2: ORA plotting (GO, KEGG, MSigDB, Reactome)
ora_results <- test_fast_ora_plotting()

# Test 3: Emergency fallback
emergency_success <- test_emergency_plotting()

# Performance Summary
cat("\n📊 PERFORMANCE SUMMARY\n")
cat("=====================\n")

cat("1. GSEA Plotting Performance:\n")
for (plot_type in names(gsea_results)) {
  result <- gsea_results[[plot_type]]
  status <- if (result$success) "✅ FAST" else "❌ FAILED"
  cat(paste0("   ", plot_type, ": ", status, " (", round(result$duration, 2), "s)\n"))
}

cat("\n2. ORA Plotting Performance:\n")
total_ora_time <- 0
successful_ora <- 0
for (analysis_type in names(ora_results)) {
  result <- ora_results[[analysis_type]]
  status <- if (result$success) "✅ FAST" else "❌ FAILED"
  cat(paste0("   ", analysis_type, ": ", status, " (", round(result$duration, 2), "s, ", result$pathways, " pathways)\n"))
  if (result$success) {
    total_ora_time <- total_ora_time + result$duration
    successful_ora <- successful_ora + 1
  }
}

cat(paste0("\n3. Emergency Fallback: ", if (emergency_success) "✅ WORKING" else "❌ FAILED", "\n"))

# Overall Assessment
gsea_success_rate <- sum(sapply(gsea_results, function(x) x$success)) / length(gsea_results) * 100
ora_success_rate <- successful_ora / length(ora_results) * 100
avg_ora_time <- if (successful_ora > 0) total_ora_time / successful_ora else 0

cat("\n🎯 OVERALL ASSESSMENT\n")
cat("====================\n")
cat(paste0("GSEA Plot Success Rate: ", round(gsea_success_rate, 1), "%\n"))
cat(paste0("ORA Plot Success Rate: ", round(ora_success_rate, 1), "%\n"))
cat(paste0("Average ORA Plot Time: ", round(avg_ora_time, 2), " seconds\n"))
cat(paste0("Emergency Fallback: ", if (emergency_success) "WORKING" else "BROKEN", "\n"))

if (gsea_success_rate >= 80 && ora_success_rate >= 80 && emergency_success) {
  cat("\n✅ ULTRA-FAST PLOTTING SYSTEM WORKING!\n")
  cat("🚀 Ready to resolve user's visualization issues\n")
} else {
  cat("\n⚠️ SOME PLOTTING ISSUES REMAIN\n")
  cat("🔧 Additional debugging may be needed\n")
}

cat("\n💡 NEXT STEPS FOR USER:\n")
cat("=======================\n")
cat("1. Restart the Shiny app to load the ultra-fast plotting system\n")
cat("2. Run pathway analysis (GSEA should complete as before)\n")
cat("3. Plots should now appear within 1-2 seconds instead of getting stuck\n")
cat("4. If plots still don't appear, emergency fallback will show summary\n")
cat("5. Results table will always work regardless of plot issues\n")

cat("\n⚡ SPEED OPTIMIZATIONS IMPLEMENTED:\n")
cat("==================================\n")
cat("✅ 10-second timeout protection on all plotting\n")
cat("✅ Immediate data limiting to top N pathways only\n")
cat("✅ Simplified ggplot2 code with minimal styling\n")
cat("✅ Analysis-type-specific optimized plotting functions\n")
cat("✅ Emergency fallback with text-based summaries\n")
cat("✅ Error handling at every level with graceful degradation\n")