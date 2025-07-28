# Debug Pathway Results Display Issue
# Check why pathway analysis results are not showing in the UI
#
# Author: Prairie Genomics Team
# Date: January 27, 2025

cat("🔍 DEBUGGING PATHWAY RESULTS DISPLAY ISSUE\n")
cat("==========================================\n\n")

# Source required modules
source("pathway_analysis.R")
source("deseq2_analysis.R")

# Test function to check pathway analysis results structure
test_pathway_results_structure <- function() {
  cat("🧪 TESTING PATHWAY RESULTS STRUCTURE\n")
  cat("====================================\n\n")
  
  # Create mock DESeq2 results (similar to what would come from the app)
  set.seed(123)
  n_genes <- 200
  
  mock_deseq2_results <- data.frame(
    baseMean = runif(n_genes, 10, 1000),
    log2FoldChange = rnorm(n_genes, 0, 2),
    lfcSE = runif(n_genes, 0.1, 0.5),
    stat = rnorm(n_genes, 0, 3),
    pvalue = runif(n_genes, 0, 1),
    padj = runif(n_genes, 0, 1),
    gene = paste0("GENE_", 1:n_genes),
    gene_symbol = paste0("SYMBOL_", 1:n_genes),
    stringsAsFactors = FALSE
  )
  
  # Make some genes significant
  n_sig <- 30
  sig_idx <- sample(1:n_genes, n_sig)
  mock_deseq2_results$padj[sig_idx] <- runif(n_sig, 0.001, 0.049)
  mock_deseq2_results$log2FoldChange[sig_idx] <- rnorm(n_sig, 0, 2)
  
  cat("✅ Mock DESeq2 results created:\n")
  cat("  - Total genes:", n_genes, "\n")
  cat("  - Significant genes:", sum(mock_deseq2_results$padj < 0.05, na.rm = TRUE), "\n")
  
  # Test pathway analysis call
  cat("\n🚀 Testing pathway analysis function...\n")
  
  pathway_result <- tryCatch({
    run_pathway_analysis(
      deseq2_results = mock_deseq2_results,
      analysis_type = "GO",
      species = "human",
      ontology = "BP",
      padj_cutoff = 0.05,
      fc_cutoff = 1.0
    )
  }, error = function(e) {
    list(
      success = FALSE,
      error = e$message,
      message = paste("Pathway analysis failed:", e$message)
    )
  })
  
  cat("⏱️ Pathway analysis completed\n")
  
  # Examine the result structure
  cat("\n🔍 EXAMINING PATHWAY RESULT STRUCTURE:\n")
  cat("=====================================\n")
  
  cat("Result class:", class(pathway_result), "\n")
  cat("Result names:", paste(names(pathway_result), collapse = ", "), "\n")
  
  # Check critical fields
  if ("success" %in% names(pathway_result)) {
    cat("✅ Has 'success' field:", pathway_result$success, "\n")
  } else {
    cat("❌ Missing 'success' field\n")
  }
  
  if ("data" %in% names(pathway_result)) {
    if (!is.null(pathway_result$data)) {
      cat("✅ Has 'data' field with", nrow(pathway_result$data), "rows\n")
      cat("  - Data columns:", paste(colnames(pathway_result$data), collapse = ", "), "\n")
    } else {
      cat("⚠️ Has 'data' field but it's NULL\n")
    }
  } else {
    cat("❌ Missing 'data' field\n")
  }
  
  if ("message" %in% names(pathway_result)) {
    cat("📝 Message:", pathway_result$message, "\n")
  }
  
  if ("error" %in% names(pathway_result)) {
    cat("❌ Error:", pathway_result$error, "\n")
  }
  
  # Check what the UI condition would evaluate to
  ui_condition <- !is.null(pathway_result) && 
                  ("success" %in% names(pathway_result)) && 
                  pathway_result$success
  
  cat("\n🎯 UI DISPLAY CONDITION:\n")
  cat("=======================\n")
  cat("!is.null(pathway_result):", !is.null(pathway_result), "\n")
  if (!is.null(pathway_result) && "success" %in% names(pathway_result)) {
    cat("pathway_result$success:", pathway_result$success, "\n")
  }
  cat("Overall UI condition (show_pathway_results):", ui_condition, "\n")
  
  if (!ui_condition) {
    cat("\n❌ PROBLEM IDENTIFIED: UI condition is FALSE\n")
    cat("🔧 This is why pathway results are not displaying in the UI\n")
    
    # Provide specific diagnosis
    if (is.null(pathway_result)) {
      cat("📋 Issue: pathway_result is NULL\n")
    } else if (!("success" %in% names(pathway_result))) {
      cat("📋 Issue: pathway_result missing 'success' field\n")
    } else if (!pathway_result$success) {
      cat("📋 Issue: pathway_result$success is FALSE\n")
      if ("message" %in% names(pathway_result)) {
        cat("   Reason:", pathway_result$message, "\n")
      }
      if ("error" %in% names(pathway_result)) {
        cat("   Error:", pathway_result$error, "\n")
      }
    }
  } else {
    cat("\n✅ UI condition is TRUE - results should display\n")
  }
  
  return(pathway_result)
}

# Test the specific pathway analysis functions
test_specific_pathway_functions <- function() {
  cat("\n🧬 TESTING SPECIFIC PATHWAY FUNCTIONS\n")
  cat("====================================\n")
  
  # Test GO analysis specifically
  cat("Testing GO analysis function directly...\n")
  
  # Create simple gene list
  gene_list <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "AKT1", "MTOR", "RB1", "CDK4")
  
  go_result <- tryCatch({
    run_go_analysis(gene_list, "human", "BP")
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  cat("GO analysis result:\n")
  cat("  - Success:", go_result$success, "\n")
  if (go_result$success && !is.null(go_result$data)) {
    cat("  - GO terms found:", nrow(go_result$data), "\n")
  }
  if (!go_result$success && "error" %in% names(go_result)) {
    cat("  - Error:", go_result$error, "\n")
  }
  
  return(go_result)
}

# Check if the pathway analysis functions exist
check_function_availability <- function() {
  cat("\n📋 CHECKING FUNCTION AVAILABILITY\n")
  cat("=================================\n")
  
  functions_to_check <- c(
    "run_pathway_analysis",
    "run_go_analysis", 
    "run_kegg_analysis",
    "prepare_gene_list_ora"
  )
  
  for (func_name in functions_to_check) {
    if (exists(func_name)) {
      cat("✅", func_name, "- Available\n")
    } else {
      cat("❌", func_name, "- Missing\n")
    }
  }
}

# Run all diagnostic tests
cat("🚀 STARTING PATHWAY DISPLAY DIAGNOSTIC\n")
cat("======================================\n\n")

# Test 1: Check function availability
check_function_availability()

# Test 2: Test pathway result structure
pathway_result <- test_pathway_results_structure()

# Test 3: Test specific functions
go_result <- test_specific_pathway_functions()

# Summary and recommendations
cat("\n📋 DIAGNOSTIC SUMMARY\n")
cat("====================\n")

if (exists("run_pathway_analysis")) {
  cat("✅ Core pathway analysis function exists\n")
} else {
  cat("❌ Missing run_pathway_analysis function\n")
}

if (!is.null(pathway_result) && "success" %in% names(pathway_result)) {
  if (pathway_result$success) {
    cat("✅ Pathway analysis completes successfully\n")
  } else {
    cat("❌ Pathway analysis fails\n")
    cat("   This is likely why results don't show in UI\n")
  }
} else {
  cat("❌ Pathway analysis returns invalid structure\n")
}

cat("\n💡 RECOMMENDATIONS:\n")
cat("===================\n")

if (!exists("run_pathway_analysis")) {
  cat("1. ✅ Ensure pathway_analysis.R is properly loaded\n")
  cat("2. ✅ Check that all required packages are installed\n")
} else {
  if (!is.null(pathway_result) && "success" %in% names(pathway_result) && !pathway_result$success) {
    cat("1. ✅ Fix pathway analysis errors (check package installation)\n")
    cat("2. ✅ Ensure gene list preparation is working correctly\n")
    cat("3. ✅ Check internet connection for KEGG analysis\n")
  } else {
    cat("1. ✅ Check app.R server logic for pathway_results assignment\n")
    cat("2. ✅ Verify reactive values are being updated correctly\n")
    cat("3. ✅ Check browser console for JavaScript errors\n")
  }
}

cat("\n🎯 NEXT STEPS:\n")
cat("==============\n")
cat("1. Run this diagnostic script to identify the specific issue\n")
cat("2. Address the root cause (function missing, analysis failing, or UI logic)\n")
cat("3. Test with a simple pathway analysis to verify the fix\n")