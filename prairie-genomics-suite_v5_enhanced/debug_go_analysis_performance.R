# Debug GO Analysis Performance and Gene Selection Issues
# Investigate why GO analysis is slow and verify correct gene filtering
#
# Author: Prairie Genomics Team
# Date: January 27, 2025

cat("🔍 DEBUGGING GO ANALYSIS PERFORMANCE ISSUES\n")
cat("===========================================\n\n")

# Function to trace exactly what genes are being used in GO analysis
trace_go_analysis_genes <- function(deseq2_results, padj_cutoff = 0.05, fc_cutoff = 1.0) {
  cat("📊 STEP-BY-STEP GO ANALYSIS GENE TRACING\n")
  cat("========================================\n\n")
  
  # Step 1: Check input data
  cat("1. INPUT DATA ANALYSIS\n")
  cat("----------------------\n")
  cat("Total genes in DESeq2 results:", nrow(deseq2_results), "\n")
  cat("Data structure:", class(deseq2_results), "\n")
  cat("Columns:", paste(colnames(deseq2_results), collapse = ", "), "\n\n")
  
  # Step 2: Manual filtering to see what SHOULD be used
  cat("2. MANUAL GENE FILTERING\n")
  cat("------------------------\n")
  
  # Apply the exact same filters as prepare_gene_list_ora
  valid_data <- deseq2_results[
    !is.na(deseq2_results$padj) & 
    !is.na(deseq2_results$log2FoldChange), 
  ]
  cat("After removing NAs:", nrow(valid_data), "genes\n")
  
  significant_genes <- valid_data[
    valid_data$padj < padj_cutoff & 
    abs(valid_data$log2FoldChange) > fc_cutoff, 
  ]
  cat("After significance filters:", nrow(significant_genes), "genes\n")
  
  if (nrow(significant_genes) > 0) {
    cat("Sample significant genes:", paste(head(rownames(significant_genes), 5), collapse = ", "), "\n")
    cat("Top 3 by significance:\n")
    top_genes <- significant_genes[order(significant_genes$padj), ][1:min(3, nrow(significant_genes)), ]
    for (i in 1:nrow(top_genes)) {
      cat(sprintf("  %s: padj=%.2e, FC=%.2f\n", 
                  rownames(top_genes)[i], 
                  top_genes$padj[i], 
                  top_genes$log2FoldChange[i]))
    }
  }
  cat("\n")
  
  # Step 3: Test the actual prepare_gene_list_ora function
  cat("3. TESTING prepare_gene_list_ora FUNCTION\n")
  cat("-----------------------------------------\n")
  
  start_time <- Sys.time()
  gene_list <- prepare_gene_list_ora(deseq2_results, padj_cutoff, fc_cutoff, "human")
  prep_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  cat("Function execution time:", round(prep_time, 2), "seconds\n")
  
  if (!is.null(gene_list) && length(gene_list) > 0) {
    cat("Genes returned by function:", length(gene_list), "\n")
    cat("Sample genes from function:", paste(head(gene_list, 5), collapse = ", "), "\n")
    
    # Check if these are the same genes we manually filtered
    if (nrow(significant_genes) > 0) {
      manual_count <- nrow(significant_genes)
      function_count <- length(gene_list)
      cat("Manual filtering found:", manual_count, "genes\n")
      cat("Function returned:", function_count, "genes\n")
      
      if (function_count < manual_count / 10) {
        cat("⚠️ WARNING: Function returned far fewer genes than expected!\n")
        cat("This suggests gene ID conversion is the bottleneck.\n")
      }
    }
  } else {
    cat("❌ Function returned NULL or empty list\n")
  }
  cat("\n")
  
  # Step 4: Test GO analysis timing
  cat("4. TESTING GO ANALYSIS TIMING\n")
  cat("------------------------------\n")
  
  if (!is.null(gene_list) && length(gene_list) > 0) {
    cat("Testing GO analysis with", length(gene_list), "genes...\n")
    
    start_time <- Sys.time()
    go_result <- run_go_analysis(gene_list, "human", "BP")
    go_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    cat("GO analysis execution time:", round(go_time, 2), "seconds\n")
    
    if (go_result$success) {
      cat("✅ GO analysis successful\n")
      cat("Terms found:", nrow(go_result$data), "\n")
    } else {
      cat("❌ GO analysis failed:", go_result$error, "\n")
    }
    
    # Performance assessment
    if (go_time > 60) {
      cat("🚨 PERFORMANCE ISSUE: GO analysis took >60 seconds\n")
      cat("This suggests the gene list might still be too large.\n")
    } else if (go_time > 30) {
      cat("⚠️ SLOW: GO analysis took >30 seconds\n")
      cat("Consider further optimization.\n")
    } else {
      cat("✅ GOOD: GO analysis completed in reasonable time\n")
    }
  }
  cat("\n")
  
  return(list(
    total_genes = nrow(deseq2_results),
    significant_genes = if(exists("significant_genes")) nrow(significant_genes) else 0,
    converted_genes = if(!is.null(gene_list)) length(gene_list) else 0,
    prep_time = if(exists("prep_time")) prep_time else NA,
    go_time = if(exists("go_time")) go_time else NA
  ))
}

# Function to check for performance bottlenecks in the GO analysis pipeline
check_go_performance_bottlenecks <- function() {
  cat("5. CHECKING FOR PERFORMANCE BOTTLENECKS\n")
  cat("========================================\n\n")
  
  # Test with different gene list sizes to identify bottleneck
  test_sizes <- c(10, 50, 100, 500, 1000, 2000)
  
  for (size in test_sizes) {
    cat("Testing with", size, "genes...\n")
    
    # Create test gene list
    test_genes <- as.character(1:size)  # Use simple Entrez IDs
    
    start_time <- Sys.time()
    go_result <- tryCatch({
      run_go_analysis(test_genes, "human", "BP")
    }, error = function(e) {
      list(success = FALSE, error = e$message)
    })
    end_time <- Sys.time()
    
    exec_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    if (go_result$success) {
      cat("  ✅", size, "genes:", round(exec_time, 2), "seconds\n")
    } else {
      cat("  ❌", size, "genes: FAILED after", round(exec_time, 2), "seconds\n")
    }
    
    # Stop if we hit a performance wall
    if (exec_time > 120) {
      cat("  🚨 Performance wall hit at", size, "genes\n")
      break
    }
  }
  cat("\n")
}

# Create a lightweight test to verify gene filtering
quick_gene_filtering_test <- function() {
  cat("6. QUICK GENE FILTERING VERIFICATION\n")
  cat("====================================\n\n")
  
  # Create minimal test data
  set.seed(123)
  test_data <- data.frame(
    baseMean = runif(100, 10, 1000),
    log2FoldChange = rnorm(100, 0, 2),
    lfcSE = runif(100, 0.1, 0.5),
    stat = rnorm(100, 0, 3),
    pvalue = runif(100, 0, 1),
    padj = runif(100, 0, 1)
  )
  
  # Make 10 genes significant
  test_data$padj[1:10] <- runif(10, 0.001, 0.049)
  test_data$log2FoldChange[1:10] <- rnorm(10, 0, 3)
  
  rownames(test_data) <- paste0("Gene_", 1:100)
  
  cat("Created test data: 100 genes, 10 significant\n")
  
  # Test filtering
  manual_sig <- test_data[
    !is.na(test_data$padj) & 
    !is.na(test_data$log2FoldChange) &
    test_data$padj < 0.05 & 
    abs(test_data$log2FoldChange) > 1.0,
  ]
  
  cat("Manual filtering result:", nrow(manual_sig), "genes\n")
  
  # Test function
  function_result <- prepare_gene_list_ora(test_data, 0.05, 1.0, "human")
  
  if (!is.null(function_result)) {
    cat("Function result:", length(function_result), "genes\n")
    
    if (nrow(manual_sig) > 0 && length(function_result) > 0) {
      cat("✅ Both manual and function filtering found genes\n")
    } else if (nrow(manual_sig) == 0) {
      cat("⚠️ No genes meet significance criteria in test data\n")
    } else {
      cat("❌ Function filtering failed despite manual success\n")
    }
  } else {
    cat("❌ Function returned NULL\n")
  }
  cat("\n")
}

cat("🔧 GO Analysis Debugging Functions Loaded\n")
cat("==========================================\n\n")
cat("Available functions:\n")
cat("• trace_go_analysis_genes(deseq2_results) - Trace gene selection process\n")
cat("• check_go_performance_bottlenecks() - Find performance issues\n") 
cat("• quick_gene_filtering_test() - Verify gene filtering works\n\n")
cat("💡 Run these functions to diagnose GO analysis issues!\n")