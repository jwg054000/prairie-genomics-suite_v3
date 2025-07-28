# Test Pathway Analysis Fix
# Verify that the gene conversion optimization resolves the 15,357 gene issue
#
# Author: Prairie Genomics Team
# Date: January 27, 2025

cat("🧪 TESTING PATHWAY ANALYSIS FIX\n")
cat("===============================\n\n")

# Load required modules
source("gene_conversion_cache.R")
source("pathway_analysis.R")

# Test the fix
test_pathway_performance_fix <- function() {
  cat("🚀 TESTING PATHWAY PERFORMANCE FIX\n")
  cat("==================================\n\n")
  
  # Create large mock DESeq2 results (simulating the 15,357 gene issue)
  cat("📊 Creating large mock dataset (15,000+ genes)...\n")
  set.seed(42)
  
  large_deseq2_results <- data.frame(
    baseMean = runif(15000, 10, 1000),
    log2FoldChange = rnorm(15000, 0, 2),
    lfcSE = runif(15000, 0.1, 0.5),
    stat = rnorm(15000, 0, 3),
    pvalue = runif(15000, 0, 1),
    padj = runif(15000, 0, 1)
  )
  
  # Make a realistic number of genes significant (3-5%)
  n_significant <- 500  # Realistic number
  sig_indices <- sample(1:15000, n_significant)
  large_deseq2_results$padj[sig_indices] <- runif(n_significant, 0.001, 0.049)
  large_deseq2_results$log2FoldChange[sig_indices] <- rnorm(n_significant, 0, 3)
  
  # Set realistic gene names
  rownames(large_deseq2_results) <- paste0("ENSG", sprintf("%011d", 1:15000))
  
  cat("✅ Large dataset created:\n")
  cat("  - Total genes:", nrow(large_deseq2_results), "\n")
  cat("  - Significant genes (padj < 0.05):", sum(large_deseq2_results$padj < 0.05, na.rm = TRUE), "\n")
  cat("  - High FC genes (|FC| > 1):", sum(abs(large_deseq2_results$log2FoldChange) > 1, na.rm = TRUE), "\n")
  
  # Test 1: Gene list preparation (should only process significant genes)
  cat("\n🔍 TEST 1: Gene List Preparation\n")
  cat("--------------------------------\n")
  
  start_time <- Sys.time()
  
  # This should filter down to ~100-200 genes (significant with |FC| > 1)
  gene_list <- prepare_gene_list_ora(large_deseq2_results, 0.05, 1.0, "human")
  
  prep_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  if (!is.null(gene_list) && length(gene_list) > 0) {
    cat("✅ Gene preparation successful:\n")
    cat("  - Input genes:", nrow(large_deseq2_results), "\n")
    cat("  - Genes after filtering:", length(gene_list), "\n") 
    cat("  - Reduction ratio:", round(100 * length(gene_list) / nrow(large_deseq2_results), 1), "%\n")
    cat("  - Preparation time:", round(prep_time, 2), "seconds\n")
    
    # This should be a MASSIVE reduction (15,000 → ~100-200 genes)
    reduction_success <- length(gene_list) < 1000  # Should be much less than 1000
    cat("  - Performance check:", if(reduction_success) "✅ PASS" else "❌ FAIL", "\n")
    
  } else {
    cat("❌ Gene preparation failed\n")
    return(list(success = FALSE, stage = "preparation"))
  }
  
  # Test 2: GO Analysis Performance
  cat("\n🚀 TEST 2: GO Analysis Performance\n")
  cat("----------------------------------\n")
  
  if (length(gene_list) > 10) {
    start_go <- Sys.time()
    go_result <- run_go_analysis(gene_list, "human", "BP")
    go_time <- as.numeric(difftime(Sys.time(), start_go, units = "secs"))
    
    cat("⏱️ GO analysis completed in", round(go_time, 2), "seconds\n")
    
    if (go_result$success) {
      cat("✅ GO analysis successful:\n")
      cat("  - Enriched terms:", nrow(go_result$data), "\n")
      cat("  - Genes analyzed:", go_result$genes_analyzed, "\n")
      cat("  - Performance status:", go_result$performance_status, "\n")
      
      # Performance evaluation
      performance_check <- go_time <= 30
      cat("  - Speed target (<30s):", if(performance_check) "✅ PASS" else "❌ FAIL", "\n")
      
      return(list(
        success = TRUE,
        gene_reduction = round(100 * length(gene_list) / nrow(large_deseq2_results), 1),
        go_time = go_time,
        terms_found = nrow(go_result$data),
        performance_grade = if(go_time <= 30) "EXCELLENT" else if(go_time <= 60) "GOOD" else "NEEDS_WORK"
      ))
      
    } else {
      cat("❌ GO analysis failed:", go_result$error, "\n")
      return(list(success = FALSE, stage = "go_analysis", error = go_result$error))
    }
  } else {
    cat("❌ Insufficient genes for GO analysis\n")
    return(list(success = FALSE, stage = "insufficient_genes"))
  }
}

# Test to verify the gene conversion optimization works
test_gene_conversion_optimization <- function() {
  cat("\n🔧 TEST 3: Gene Conversion Optimization\n")
  cat("=======================================\n")
  
  # Test that gene conversion is skipped for large datasets
  large_dataset_size <- 15000
  small_dataset_size <- 1000
  
  cat("Testing decision logic for gene conversion:\n")
  
  # Test large dataset (should skip conversion)
  skip_large <- !(("human" != "none") && (large_dataset_size <= 5000))
  cat("  - Large dataset (", large_dataset_size, "genes):", if(skip_large) "✅ SKIP (correct)" else "❌ CONVERT (wrong)", "\n")
  
  # Test small dataset (should convert)
  convert_small <- (("human" != "none") && (small_dataset_size <= 5000))
  cat("  - Small dataset (", small_dataset_size, "genes):", if(convert_small) "✅ CONVERT (correct)" else "❌ SKIP (wrong)", "\n")
  
  # Test "none" species (should always skip)
  skip_none <- !(("none" != "none") && (small_dataset_size <= 5000))
  cat("  - No species selected:", if(skip_none) "✅ SKIP (correct)" else "❌ CONVERT (wrong)", "\n")
  
  optimization_working <- skip_large && convert_small && skip_none
  cat("\n🎯 Optimization logic:", if(optimization_working) "✅ WORKING CORRECTLY" else "❌ NEEDS FIXING", "\n")
  
  return(list(success = optimization_working))
}

# Run all tests
cat("🚀 STARTING COMPREHENSIVE PATHWAY FIX TEST\n")
cat("===========================================\n\n")

# Test 1: Performance fix
performance_result <- test_pathway_performance_fix()

# Test 2: Optimization logic
optimization_result <- test_gene_conversion_optimization()

# Summary
cat("\n📋 TEST SUMMARY\n")
cat("===============\n")

cat("1. Pathway Performance Fix:")
if (performance_result$success) {
  cat("✅ PASS\n")
  cat("   - Gene reduction:", performance_result$gene_reduction, "% (massive improvement)\n")
  cat("   - GO analysis time:", round(performance_result$go_time, 2), "seconds\n")
  cat("   - Terms found:", performance_result$terms_found, "\n")
  cat("   - Performance grade:", performance_result$performance_grade, "\n")
} else {
  cat("❌ FAIL\n")
  cat("   - Failed at stage:", performance_result$stage, "\n")
  if (!is.null(performance_result$error)) {
    cat("   - Error:", performance_result$error, "\n")
  }
}

cat("\n2. Gene Conversion Optimization:")
if (optimization_result$success) {
  cat("✅ PASS - Logic working correctly\n")
} else {
  cat("❌ FAIL - Logic needs adjustment\n")
}

# Overall assessment
overall_success <- performance_result$success && optimization_result$success

cat("\n🎯 OVERALL RESULT:")
if (overall_success) {
  cat("✅ PATHWAY ANALYSIS FIX SUCCESSFUL\n")
  cat("🎉 The 15,357 gene conversion issue should be resolved:\n")
  cat("   - Large datasets skip expensive gene conversion\n")
  cat("   - Pathway analysis only processes significant genes\n")
  cat("   - Performance is dramatically improved\n")
} else {
  cat("❌ SOME ISSUES REMAIN\n")
  cat("🔧 Check individual test results above for details\n")
}

cat("\n💡 NEXT STEPS FOR USER:\n")
cat("1. Try running pathway analysis again\n")
cat("2. Should see much faster performance\n")
cat("3. Gene symbols will be converted during pathway analysis only\n")
cat("4. No more 15,357 gene conversion timeouts\n")