# Test Improved Gene Conversion with v5 Implementation
# Tests the updated biomaRt conversion with multiple mirrors and retry logic

cat("🧬 Testing Improved Gene Conversion (v5 Implementation)\n")
cat("=" , rep("=", 60), "\n")

# Load modules
source("config/app_config.R")
source("modules/data_processing/file_upload.R")

# Test with real Ensembl IDs from test data
cat("📊 Loading test data...\n")
if (file.exists("test_data_with_patterns.csv")) {
  test_data <- read.csv("test_data_with_patterns.csv", stringsAsFactors = FALSE)
  test_genes <- test_data$Gene[1:5]  # Test with first 5 genes
  
  cat("✅ Test genes loaded:", paste(test_genes, collapse = ", "), "\n")
  
  # Test the improved conversion function
  cat("\n🔬 Testing improved gene conversion...\n")
  cat("Key improvements from v5:\n")
  cat("  - useEnsembl() instead of useMart()\n") 
  cat("  - Multiple mirror fallback (www, useast, asia)\n")
  cat("  - Smaller batch sizes (200 vs 1000)\n")
  cat("  - Progressive retry with backoff\n")
  cat("  - Rate limiting between batches\n")
  cat("  - Proper symbol attributes (hgnc_symbol, mgi_symbol)\n\n")
  
  # Test human conversion
  cat("🧬 Testing human gene conversion...\n")
  human_result <- tryCatch({
    convert_gene_ids_to_symbols(test_genes, species = "human")
  }, error = function(e) {
    cat("❌ Human conversion failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(human_result)) {
    cat("✅ Human conversion completed:\n")
    for (i in 1:nrow(human_result)) {
      cat("   -", human_result$original_id[i], "→", human_result$gene_symbol[i], "\n")
    }
    
    # Check conversion success
    converted <- sum(human_result$gene_symbol != human_result$original_id)
    cat("📊 Conversion rate:", converted, "out of", nrow(human_result), "genes\n")
  }
  
  # Test mouse conversion
  cat("\n🐭 Testing mouse gene conversion...\n")
  mouse_result <- tryCatch({
    convert_gene_ids_to_symbols(test_genes[1:2], species = "mouse")  # Just 2 genes for mouse
  }, error = function(e) {
    cat("❌ Mouse conversion failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(mouse_result)) {
    cat("✅ Mouse conversion completed:\n")
    for (i in 1:nrow(mouse_result)) {
      cat("   -", mouse_result$original_id[i], "→", mouse_result$gene_symbol[i], "\n")
    }
  }
  
} else {
  cat("❌ Test data file not found\n")
}

# Test the matrix conversion function
cat("\n📊 Testing matrix conversion with improved function...\n")
if (exists("test_genes")) {
  # Create small test matrix
  test_matrix <- matrix(
    c(100, 200, 150,
      50, 75, 60,
      300, 250, 280),
    nrow = 3,
    dimnames = list(test_genes[1:3], c("Sample1", "Sample2", "Sample3"))
  )
  
  cat("📋 Test matrix:\n")
  print(test_matrix)
  
  cat("\n🧬 Applying matrix conversion...\n")
  matrix_result <- apply_gene_conversion(test_matrix, species = "human", enable_conversion = TRUE)
  
  cat("✅ Matrix conversion results:\n")
  cat("   - Attempted:", matrix_result$conversion_stats$attempted, "\n")
  cat("   - Total genes:", matrix_result$conversion_stats$total_genes, "\n")
  cat("   - Converted count:", matrix_result$conversion_stats$converted_count, "\n")
  cat("   - Conversion rate:", matrix_result$conversion_stats$conversion_rate, "%\n")
  
  if (matrix_result$conversion_stats$attempted) {
    cat("   - New row names:", paste(rownames(matrix_result$matrix), collapse = ", "), "\n")
  }
}

cat("\n🎉 IMPROVED GENE CONVERSION TEST COMPLETED!\n")
cat("=" , rep("=", 60), "\n")

if (biomart_available) {
  cat("✅ biomaRt is available - testing with improved v5 implementation\n")
  cat("🔧 Key improvements implemented:\n")
  cat("   ✅ Multiple Ensembl mirrors (www, useast, asia)\n")
  cat("   ✅ Progressive retry logic with backoff\n")
  cat("   ✅ Smaller batch sizes for reliability (200 genes)\n")
  cat("   ✅ Rate limiting between batches\n")
  cat("   ✅ Proper symbol attributes per species\n")
  cat("   ✅ Legacy fallback connection method\n")
} else {
  cat("⚠️ biomaRt not available - conversion will use fallback mode\n")
}

cat("\n🚀 Ready to test in the Shiny application!\n\n")

# Cleanup
if (exists("test_matrix")) rm(test_matrix)
if (exists("human_result")) rm(human_result)
if (exists("mouse_result")) rm(mouse_result)
smart_gc()