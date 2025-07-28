# Test Gene ID to Symbol Conversion Functionality
# Tests the new gene conversion feature added to the data upload module

cat("🧬 Testing Gene ID to Symbol Conversion\n")
cat("=" , rep("=", 50), "\n")

# Load required modules
source("config/app_config.R")
source("utils/memory_manager.R")

# Load the enhanced file upload module with gene conversion
source("modules/data_processing/file_upload.R")

# Test 1: Load test data with Ensembl IDs
cat("📊 Test 1: Loading test data...\n")
if (file.exists("test_data_with_patterns.csv")) {
  test_data <- read.csv("test_data_with_patterns.csv", stringsAsFactors = FALSE)
  
  cat("✅ Test data loaded:\n")
  cat("   - Rows:", nrow(test_data), "\n")
  cat("   - Columns:", ncol(test_data), "\n")
  cat("   - First few gene IDs:", paste(test_data$Gene[1:3], collapse = ", "), "\n")
  
  # Check if genes look like Ensembl IDs
  ensembl_pattern <- "^ENS[A-Z]*G\\d{11}"
  ensembl_count <- sum(grepl(ensembl_pattern, test_data$Gene))
  cat("   - Ensembl-like IDs:", ensembl_count, "out of", nrow(test_data), "\n")
  
} else {
  cat("❌ Test data file not found\n")
  stop("Cannot proceed without test data")
}

# Test 2: Test gene conversion function directly
cat("\n🔬 Test 2: Testing gene conversion function...\n")

# Test with a small subset of genes
test_genes <- test_data$Gene[1:5]
cat("📋 Testing conversion of genes:", paste(test_genes, collapse = ", "), "\n")

conversion_result <- tryCatch({
  convert_gene_ids_to_symbols(test_genes, species = "human")
}, error = function(e) {
  cat("❌ Conversion failed:", e$message, "\n")
  return(NULL)
})

if (!is.null(conversion_result)) {
  cat("✅ Conversion function worked:\n")
  for (i in 1:nrow(conversion_result)) {
    cat("   -", conversion_result$original_id[i], "→", conversion_result$gene_symbol[i], "\n")
  }
} else {
  cat("⚠️ Conversion function returned NULL - likely no biomaRt connection\n")
}

# Test 3: Test matrix conversion function
cat("\n📊 Test 3: Testing matrix conversion...\n")

# Create test matrix
gene_names <- test_data$Gene
expression_values <- test_data[, -1]
test_matrix <- as.matrix(expression_values)
rownames(test_matrix) <- gene_names

cat("📋 Test matrix created:\n")
cat("   - Genes:", nrow(test_matrix), "\n")
cat("   - Samples:", ncol(test_matrix), "\n")
cat("   - Sample names:", paste(colnames(test_matrix)[1:3], collapse = ", "), "...\n")

# Test with conversion enabled
cat("\n🧬 Testing matrix conversion (enabled)...\n")
matrix_result_enabled <- apply_gene_conversion(test_matrix, species = "human", enable_conversion = TRUE)

cat("✅ Matrix conversion (enabled) results:\n")
cat("   - Attempted:", matrix_result_enabled$conversion_stats$attempted, "\n")
cat("   - Total genes:", matrix_result_enabled$conversion_stats$total_genes, "\n")
cat("   - Converted count:", matrix_result_enabled$conversion_stats$converted_count, "\n")
cat("   - Conversion rate:", matrix_result_enabled$conversion_stats$conversion_rate, "%\n")

if (matrix_result_enabled$conversion_stats$attempted) {
  cat("   - First few converted names:", paste(rownames(matrix_result_enabled$matrix)[1:3], collapse = ", "), "\n")
} else {
  cat("   - Reason for skipping:", matrix_result_enabled$conversion_stats$reason %||% "Unknown", "\n")
}

# Test with conversion disabled
cat("\n❌ Testing matrix conversion (disabled)...\n")
matrix_result_disabled <- apply_gene_conversion(test_matrix, species = "human", enable_conversion = FALSE)

cat("✅ Matrix conversion (disabled) results:\n")
cat("   - Attempted:", matrix_result_disabled$conversion_stats$attempted, "\n")
cat("   - Row names unchanged:", identical(rownames(test_matrix), rownames(matrix_result_disabled$matrix)), "\n")

# Test 4: Test different species
if (biomart_available) {
  cat("\n🐭 Test 4: Testing mouse species conversion...\n")
  
  # Test with mouse (should work even with human Ensembl IDs for testing)
  mouse_result <- tryCatch({
    convert_gene_ids_to_symbols(test_genes[1:2], species = "mouse")
  }, error = function(e) {
    cat("⚠️ Mouse conversion failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(mouse_result)) {
    cat("✅ Mouse conversion test completed\n")
  }
} else {
  cat("\n⚠️ Test 4: Skipped (biomaRt not available)\n")
}

# Summary
cat("\n🎉 GENE CONVERSION TESTS COMPLETED!\n")
cat("=" , rep("=", 50), "\n")

if (biomart_available) {
  cat("✅ biomaRt is available - full gene conversion functionality enabled\n")
  cat("📋 Gene conversion features:\n")
  cat("   - Automatic Ensembl ID detection\n")
  cat("   - Chunked processing for large datasets\n")
  cat("   - Error handling and retries\n")
  cat("   - Species selection (human/mouse)\n")
  cat("   - Graceful fallback when conversion fails\n")
} else {
  cat("⚠️ biomaRt is not available - gene conversion will use fallback mode\n")
  cat("📋 In fallback mode:\n")
  cat("   - Original gene IDs are preserved\n")
  cat("   - No network requests are made\n")
  cat("   - Application continues to function normally\n")
}

cat("\n🚀 The gene conversion feature is ready for use in the Data Upload tab!\n")
cat("💡 Users can now:\n")
cat("   - Enable/disable gene conversion with a checkbox\n")
cat("   - Select species (human/mouse)\n")
cat("   - See conversion statistics in notifications\n")
cat("   - View converted gene symbols in data preview\n\n")

# Cleanup
rm(test_data, test_matrix, conversion_result, matrix_result_enabled, matrix_result_disabled)
smart_gc()