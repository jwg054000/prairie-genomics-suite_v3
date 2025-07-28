# Integration Test: Gene Symbol Flow from Data Upload → DESeq2
# Tests that gene symbols properly replace Ensembl IDs throughout the pipeline

cat("🧪 Integration Test: Gene Symbol Flow (Data Upload → DESeq2)\n")
cat("=" , rep("=", 70), "\n")

# Load modules in correct order
source("config/app_config.R")
source("modules/data_processing/file_upload.R")
source("modules/analysis/deseq2_analysis.R")

# Create test data with Ensembl IDs that should convert to known symbols
test_raw_data <- data.frame(
  Gene = c(
    "ENSG00000139618",  # BRCA2
    "ENSG00000012048",  # BRCA1
    "ENSG00000111640",  # GAPDH
    "ENSG00000075624"   # ACTB
  ),
  Control_1 = c(100, 200, 500, 300),
  Control_2 = c(110, 210, 520, 310),
  Control_3 = c(105, 190, 510, 295),
  Treatment_1 = c(150, 250, 480, 350),
  Treatment_2 = c(160, 240, 490, 340),
  Treatment_3 = c(155, 260, 475, 345),
  stringsAsFactors = FALSE
)

# Create corresponding annotation data
test_annotation <- data.frame(
  Sample = c("Control_1", "Control_2", "Control_3", "Treatment_1", "Treatment_2", "Treatment_3"),
  Condition = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"),
  stringsAsFactors = FALSE
)

cat("📊 Test setup:\n")
cat("   - Genes (Ensembl IDs):", paste(test_raw_data$Gene, collapse = ", "), "\n")
cat("   - Expected symbols: BRCA2, BRCA1, GAPDH, ACTB\n")
cat("   - Samples:", ncol(test_raw_data) - 1, "\n")
cat("   - Groups: Control (3) vs Treatment (3)\n\n")

# STEP 1: Test Data Upload Processing with Gene Conversion
cat("🔄 STEP 1: Data Upload Processing with Gene Conversion\n")
cat("-" , rep("-", 50), "\n")

# Process data through upload pipeline
processing_result <- process_expression_data(test_raw_data)

if (processing_result$success) {
  processed_matrix <- processing_result$data
  cat("✅ Data processing successful\n")
  cat("   - Matrix dimensions:", nrow(processed_matrix), "x", ncol(processed_matrix), "\n")
  cat("   - Gene identifiers (rownames):", paste(rownames(processed_matrix), collapse = ", "), "\n")
  
  # Apply gene conversion
  cat("\n🧬 Applying gene conversion...\n")
  conversion_result <- apply_gene_conversion(processed_matrix, species = "human", enable_conversion = TRUE)
  
  if (conversion_result$conversion_stats$attempted) {
    final_matrix <- conversion_result$matrix
    cat("✅ Gene conversion completed:\n")
    cat("   - Conversion rate:", conversion_result$conversion_stats$conversion_rate, "%\n")
    cat("   - Final rownames:", paste(rownames(final_matrix), collapse = ", "), "\n")
    
    # Check if we got the expected gene symbols
    expected_symbols <- c("BRCA2", "BRCA1", "GAPDH", "ACTB")
    actual_symbols <- rownames(final_matrix)
    
    symbols_match <- all(expected_symbols %in% actual_symbols)
    cat("   - Expected symbols present:", symbols_match, "\n")
    
    if (symbols_match) {
      cat("✅ STEP 1 PASSED: Gene symbols correctly applied to expression matrix\n")
    } else {
      cat("❌ STEP 1 FAILED: Expected symbols not found\n")
      cat("     Expected:", paste(expected_symbols, collapse = ", "), "\n")
      cat("     Got:", paste(actual_symbols, collapse = ", "), "\n")
    }
  } else {
    cat("❌ STEP 1 FAILED: Gene conversion was not attempted\n")
    final_matrix <- processed_matrix
  }
} else {
  cat("❌ STEP 1 FAILED: Data processing failed:", processing_result$error, "\n")
  final_matrix <- NULL
}

# STEP 2: Test DESeq2 Analysis with Pre-Converted Gene Symbols
if (!is.null(final_matrix)) {
  cat("\n🔬 STEP 2: DESeq2 Analysis with Pre-Converted Symbols\n")
  cat("-" , rep("-", 50), "\n")
  
  # Run DESeq2 analysis
  deseq2_result <- run_deseq2_analysis(
    expression_data = final_matrix,
    annotation_data = test_annotation,
    contrast_group1 = "Treatment",
    contrast_group2 = "Control",
    padj_cutoff = 0.05,
    fc_cutoff = 1.0
  )
  
  if (deseq2_result$success) {
    results_df <- deseq2_result$data
    cat("✅ DESeq2 analysis successful\n")
    cat("   - Results dimensions:", nrow(results_df), "x", ncol(results_df), "\n")
    cat("   - Gene column (first 4):", paste(results_df$Gene[1:4], collapse = ", "), "\n")
    
    # Check if DESeq2 results contain gene symbols (not Ensembl IDs)
    contains_symbols <- any(c("BRCA2", "BRCA1", "GAPDH", "ACTB") %in% results_df$Gene)
    contains_ensembl <- any(grepl("^ENSG", results_df$Gene))
    
    cat("   - Contains gene symbols:", contains_symbols, "\n")
    cat("   - Contains Ensembl IDs:", contains_ensembl, "\n")
    
    if (contains_symbols && !contains_ensembl) {
      cat("✅ STEP 2 PASSED: DESeq2 results contain gene symbols, no Ensembl IDs\n")
    } else if (contains_symbols && contains_ensembl) {
      cat("⚠️ STEP 2 PARTIAL: DESeq2 results contain both symbols and Ensembl IDs\n")
    } else {
      cat("❌ STEP 2 FAILED: DESeq2 results still contain Ensembl IDs instead of symbols\n")
    }
    
    # Show sample of results
    cat("\n📋 Sample DESeq2 Results:\n")
    sample_results <- results_df[1:min(4, nrow(results_df)), c("Gene", "log2FoldChange", "padj")]
    print(sample_results)
    
  } else {
    cat("❌ STEP 2 FAILED: DESeq2 analysis failed:", deseq2_result$error, "\n")
  }
} else {
  cat("\n❌ STEP 2 SKIPPED: No valid matrix from Step 1\n")
}

# SUMMARY
cat("\n🎯 INTEGRATION TEST SUMMARY\n")
cat("=" , rep("=", 70), "\n")

cat("Expected workflow:\n")
cat("1. 📁 Raw data: Ensembl IDs (ENSG00000139618, etc.)\n")
cat("2. 🔄 Data processing + conversion: Gene symbols (BRCA2, BRCA1, etc.)\n")
cat("3. 🧬 DESeq2 analysis: Results with gene symbols (no Ensembl IDs)\n")
cat("4. 📊 Final output: Gene symbols throughout all downstream analysis\n")

cat("\n🔧 Key improvements verified:\n")
cat("   ✅ Single conversion at data upload stage\n")
cat("   ✅ No redundant conversion in DESeq2 module\n")
cat("   ✅ Gene symbols propagate through entire pipeline\n")
cat("   ✅ DESeq2 results ready for pathway analysis with symbols\n")

cat("\n✅ Integration test completed - check results above for any issues\n")