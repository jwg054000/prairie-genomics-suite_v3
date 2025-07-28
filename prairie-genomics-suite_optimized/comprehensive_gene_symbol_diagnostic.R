# COMPREHENSIVE GENE SYMBOL DIAGNOSTIC
# Trace complete data flow: Raw Upload → Aggregation → Conversion → Display
# Identify exactly where gene symbols are lost in the pipeline

cat("🔍 COMPREHENSIVE GENE SYMBOL DIAGNOSTIC\n")
cat("=" , rep("=", 80), "\n")

# Load all modules in correct order (same as run_app.R)
source("config/app_config.R")
source("modules/data_processing/file_upload.R")

# Create test data with Ensembl IDs that should convert to known symbols
create_test_data <- function() {
  test_raw_data <- data.frame(
    Gene = c(
      "ENSG00000139618",  # BRCA2
      "ENSG00000012048",  # BRCA1
      "ENSG00000111640",  # GAPDH
      "ENSG00000075624",  # ACTB
      "ENSG00000139618",  # BRCA2 (duplicate for aggregation test)
      "ENSG00000111640"   # GAPDH (duplicate for aggregation test)
    ),
    Control_1 = c(100, 200, 500, 300, 50, 250),
    Control_2 = c(110, 210, 520, 310, 60, 260),
    Control_3 = c(105, 190, 510, 295, 55, 255),
    Treatment_1 = c(150, 250, 480, 350, 75, 240),
    Treatment_2 = c(160, 240, 490, 340, 80, 245),
    Treatment_3 = c(155, 260, 475, 345, 70, 235),
    stringsAsFactors = FALSE
  )
  
  cat("📊 Test data created:\n")
  cat("   - Rows:", nrow(test_raw_data), "\n")
  cat("   - Columns:", ncol(test_raw_data), "\n")
  cat("   - Gene column:", paste(test_raw_data$Gene, collapse = ", "), "\n")
  cat("   - Duplicates present: BRCA2 and GAPDH (should be aggregated)\n")
  cat("   - Expected symbols after conversion: BRCA2, BRCA1, GAPDH, ACTB\n\n")
  
  return(test_raw_data)
}

# STEP 1: RAW DATA PROCESSING (same as file upload)
cat("🔄 STEP 1: RAW DATA PROCESSING\n")
cat("-", rep("-", 50), "\n")

test_data <- create_test_data()

# Simulate file upload processing
processing_result <- process_expression_data(test_data)

if (processing_result$success) {
  processed_matrix <- processing_result$data
  cat("✅ Data processing successful\n")
  cat("   - Matrix dimensions:", nrow(processed_matrix), "x", ncol(processed_matrix), "\n")
  cat("   - Gene identifiers (rownames):", paste(rownames(processed_matrix), collapse = ", "), "\n")
  cat("   - Aggregation stats:\n")
  if (!is.null(processing_result$aggregation_stats)) {
    agg_stats <- processing_result$aggregation_stats
    cat("     * Duplicates found:", agg_stats$duplicates_found, "\n")
    if (agg_stats$duplicates_found) {
      cat("     * Original genes:", agg_stats$original_genes, "\n")
      cat("     * Final genes:", agg_stats$final_genes, "\n")
      cat("     * Genes aggregated:", agg_stats$genes_aggregated, "\n")
      cat("     * Method:", agg_stats$aggregation_method %||% "standard", "\n")
    }
  }
  
  cat("\n📋 Processed matrix sample:\n")
  print(processed_matrix[1:min(4, nrow(processed_matrix)), 1:min(3, ncol(processed_matrix))])
  
} else {
  cat("❌ STEP 1 FAILED: Data processing failed:", processing_result$error, "\n")
  stop("Cannot continue without successful data processing")
}

# STEP 2: GENE CONVERSION
cat("\n🧬 STEP 2: GENE CONVERSION\n")
cat("-", rep("-", 50), "\n")

# Apply gene conversion exactly as in the module
conversion_result <- apply_gene_conversion(processed_matrix, species = "human", enable_conversion = TRUE)

if (conversion_result$conversion_stats$attempted) {
  final_matrix <- conversion_result$matrix
  conversion_stats <- conversion_result$conversion_stats
  conversion_table <- conversion_result$conversion_table
  
  cat("✅ Gene conversion completed:\n")
  cat("   - Conversion attempted:", conversion_stats$attempted, "\n")
  cat("   - Total genes:", conversion_stats$total_genes, "\n")
  cat("   - Converted count:", conversion_stats$converted_count, "\n")
  cat("   - Conversion rate:", conversion_stats$conversion_rate, "%\n")
  
  cat("\n📋 Final matrix rownames after conversion:\n")
  cat("   - All rownames:", paste(rownames(final_matrix), collapse = ", "), "\n")
  
  cat("\n📊 Final matrix sample (with gene symbols as rownames):\n")
  print(final_matrix[1:min(4, nrow(final_matrix)), 1:min(3, ncol(final_matrix))])
  
  # Check if we got the expected gene symbols
  expected_symbols <- c("BRCA2", "BRCA1", "GAPDH", "ACTB")
  actual_symbols <- rownames(final_matrix)
  
  symbols_match <- any(expected_symbols %in% actual_symbols)
  ensembl_still_present <- any(grepl("^ENSG", actual_symbols))
  
  cat("\n🎯 Symbol conversion verification:\n")
  cat("   - Expected symbols present:", symbols_match, "\n")
  cat("   - Ensembl IDs still present:", ensembl_still_present, "\n")
  
  if (symbols_match && !ensembl_still_present) {
    cat("✅ STEP 2 PASSED: Gene symbols successfully replaced Ensembl IDs\n")
  } else if (symbols_match && ensembl_still_present) {
    cat("⚠️ STEP 2 PARTIAL: Both symbols and Ensembl IDs present\n")
  } else {
    cat("❌ STEP 2 FAILED: Gene symbols not properly applied\n")
  }
  
  # Show conversion table if available
  if (!is.null(conversion_table)) {
    cat("\n📋 Conversion table sample:\n")
    print(head(conversion_table, 6))
  }
  
} else {
  cat("❌ STEP 2 FAILED: Gene conversion was not attempted\n")
  final_matrix <- processed_matrix
}

# STEP 3: SIMULATE STORAGE IN REACTIVE VALUES (as done in dataUploadServer)
cat("\n💾 STEP 3: SIMULATING STORAGE IN REACTIVE VALUES\n")
cat("-", rep("-", 50), "\n")

# This simulates what happens in the dataUploadServer function
simulated_values <- list()
simulated_values$expression_data <- final_matrix  # This is line 1441 in file_upload.R

cat("✅ Data stored in simulated reactive values\n")
cat("   - Storage dimensions:", nrow(simulated_values$expression_data), "x", ncol(simulated_values$expression_data), "\n")
cat("   - Stored rownames:", paste(rownames(simulated_values$expression_data), collapse = ", "), "\n")

# STEP 4: SIMULATE DATA PREVIEW (as rendered in UI)
cat("\n📊 STEP 4: SIMULATING DATA PREVIEW RENDERING\n")
cat("-", rep("-", 50), "\n")

# This simulates the processed_preview output (lines 1537-1570 in file_upload.R)
preview_data <- as.data.frame(simulated_values$expression_data)
preview_data <- preview_data[1:min(4, nrow(preview_data)), 1:min(3, ncol(preview_data))]

# Add Gene_Symbol column as first column (exactly as in the module)
preview_data_with_symbols <- data.frame(
  Gene_Symbol = rownames(preview_data),
  preview_data,
  stringsAsFactors = FALSE
)

cat("📋 Simulated UI preview data:\n")
print(preview_data_with_symbols)

# STEP 5: SIMULATE DESEQ2 ANALYSIS INPUT
cat("\n🔬 STEP 5: SIMULATING DESEQ2 ANALYSIS INPUT\n")
cat("-", rep("-", 50), "\n")

# Check what DESeq2 would receive as input
deseq2_input_matrix <- simulated_values$expression_data

cat("📊 DESeq2 would receive:\n")
cat("   - Matrix dimensions:", nrow(deseq2_input_matrix), "x", ncol(deseq2_input_matrix), "\n")
cat("   - Matrix rownames:", paste(rownames(deseq2_input_matrix), collapse = ", "), "\n")
cat("   - Contains gene symbols:", any(c("BRCA2", "BRCA1", "GAPDH", "ACTB") %in% rownames(deseq2_input_matrix)), "\n")
cat("   - Contains Ensembl IDs:", any(grepl("^ENSG", rownames(deseq2_input_matrix))), "\n")

cat("\n📋 DESeq2 input matrix sample:\n")
print(deseq2_input_matrix[1:min(4, nrow(deseq2_input_matrix)), 1:min(3, ncol(deseq2_input_matrix))])

# FINAL SUMMARY
cat("\n🎯 COMPREHENSIVE DIAGNOSTIC SUMMARY\n")
cat("=", rep("=", 80), "\n")

cat("Data flow trace:\n")
cat("1. 📁 Raw data: 6 genes (including duplicates) with Ensembl IDs\n")
cat("2. 🔄 Processing: Applied exact user aggregation pattern\n")
cat("3. 🧬 Conversion: Converted Ensembl IDs to gene symbols\n")
cat("4. 💾 Storage: Stored in reactive values with symbols as rownames\n")
cat("5. 📊 Preview: UI should display Gene_Symbol column with converted names\n")
cat("6. 🔬 DESeq2: Analysis should receive matrix with gene symbol rownames\n")

# Key diagnostic questions
cat("\nKEY DIAGNOSTIC RESULTS:\n")
cat("❓ Are duplicates properly aggregated?", 
    if(!is.null(processing_result$aggregation_stats) && processing_result$aggregation_stats$duplicates_found) "YES" else "NO", "\n")

cat("❓ Is conversion working?", 
    if(conversion_result$conversion_stats$attempted && conversion_result$conversion_stats$converted_count > 0) "YES" else "NO", "\n")

cat("❓ Are gene symbols in final matrix rownames?", 
    if(any(c("BRCA2", "BRCA1", "GAPDH", "ACTB") %in% rownames(final_matrix))) "YES" else "NO", "\n")

cat("❓ Are gene symbols preserved in reactive values?", 
    if(any(c("BRCA2", "BRCA1", "GAPDH", "ACTB") %in% rownames(simulated_values$expression_data))) "YES" else "NO", "\n")

cat("❓ Will UI preview show gene symbols?", 
    if(any(c("BRCA2", "BRCA1", "GAPDH", "ACTB") %in% preview_data_with_symbols$Gene_Symbol)) "YES" else "NO", "\n")

cat("❓ Will DESeq2 receive gene symbols?", 
    if(any(c("BRCA2", "BRCA1", "GAPDH", "ACTB") %in% rownames(deseq2_input_matrix))) "YES" else "NO", "\n")

# If everything looks good in this test but user reports issues, the problem is likely:
cat("\n🔧 TROUBLESHOOTING GUIDANCE:\n")
cat("If this diagnostic shows gene symbols working correctly but the running app shows Ensembl IDs:\n")
cat("1. ✅ Check that run_app.R is being used as entry point\n")
cat("2. ✅ Check that the exact user aggregation pattern is being used\n") 
cat("3. ✅ Verify no caching issues in the browser\n")
cat("4. ✅ Ensure gene conversion is enabled in the UI\n")
cat("5. ✅ Check for any errors in the R console during processing\n")

cat("\n✅ Comprehensive diagnostic completed - check results above for any issues\n")