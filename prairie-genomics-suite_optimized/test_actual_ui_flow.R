# TEST ACTUAL UI DATA FLOW
# Simulate exactly what happens in the Shiny app to identify where gene symbols are lost

cat("🔍 TESTING ACTUAL UI DATA FLOW - STEP BY STEP\n")
cat("=" , rep("=", 80), "\n")

# Load modules exactly as the app does
source("config/app_config.R")
source("modules/data_processing/file_upload.R")

# Simulate reactive values structure (like Shiny does)
values <- list()
local_values <- list()

# Create test data that should convert
test_data <- data.frame(
  Gene = c("ENSG00000139618", "ENSG00000012048", "ENSG00000111640", "ENSG00000075624"),
  Control_1 = c(100, 200, 500, 300),
  Control_2 = c(110, 210, 520, 310),
  Treatment_1 = c(150, 250, 480, 350),
  Treatment_2 = c(160, 240, 490, 340),
  stringsAsFactors = FALSE
)

cat("📊 STEP 1: RAW DATA\n")
cat("Raw gene IDs:", paste(test_data$Gene, collapse = ", "), "\n\n")

# STEP 1: Process expression data (same as clicking "Process Data" button)
cat("🔄 STEP 2: PROCESSING DATA (simulating process_file button click)\n")
cat("-", rep("-", 50), "\n")

# This simulates lines 1397-1443 in file_upload.R
processing_result <- process_expression_data(test_data)

if (processing_result$success) {
  processed_matrix <- processing_result$data
  cat("✅ Processing successful\n")
  cat("Processed matrix rownames:", paste(rownames(processed_matrix), collapse = ", "), "\n")
  
  # Store aggregation stats (line 1407 in file_upload.R)
  local_values$aggregation_stats <- processing_result$aggregation_stats
  
  # STEP 2: Apply gene conversion (lines 1415-1437 in file_upload.R)
  cat("\n🧬 STEP 3: APPLYING GENE CONVERSION\n")
  cat("-", rep("-", 50), "\n")
  
  enable_conversion <- TRUE  # Simulating checkbox checked
  species <- "human"  # Simulating dropdown selection
  
  if (enable_conversion) {
    cat("🧬 Starting gene ID conversion...\n")
    
    conversion_result <- apply_gene_conversion(
      processed_matrix, 
      species = species, 
      enable_conversion = TRUE
    )
    
    processed_matrix <- conversion_result$matrix
    local_values$conversion_stats <- conversion_result$conversion_stats
    local_values$conversion_table <- conversion_result$conversion_table
    
    cat("Conversion result matrix rownames:", paste(rownames(processed_matrix), collapse = ", "), "\n")
  }
  
  # STEP 3: Store in reactive values (lines 1440-1441 in file_upload.R)
  cat("\n💾 STEP 4: STORING IN REACTIVE VALUES\n")
  cat("-", rep("-", 50), "\n")
  
  local_values$processed_data <- processed_matrix
  values$expression_data <- processed_matrix
  
  cat("local_values$processed_data rownames:", paste(rownames(local_values$processed_data), collapse = ", "), "\n")
  cat("values$expression_data rownames:", paste(rownames(values$expression_data), collapse = ", "), "\n")
  
  # STEP 4: Simulate UI preview rendering (lines 1537-1570 in file_upload.R)
  cat("\n📊 STEP 5: SIMULATING UI PREVIEW RENDERING\n")
  cat("-", rep("-", 50), "\n")
  
  # This simulates the processed_preview output
  if (!is.null(local_values$processed_data)) {
    preview_data <- as.data.frame(local_values$processed_data)
    preview_data <- preview_data[1:min(4, nrow(preview_data)), 1:min(3, ncol(preview_data))]
    
    # Add Gene_Symbol column as first column
    preview_data_with_symbols <- data.frame(
      Gene_Symbol = rownames(preview_data),
      preview_data,
      stringsAsFactors = FALSE
    )
    
    cat("📋 UI Preview Data (what user should see):\n")
    print(preview_data_with_symbols)
    
    # Check if gene symbols are actually present
    contains_symbols <- any(c("BRCA2", "BRCA1", "GAPDH", "ACTB") %in% preview_data_with_symbols$Gene_Symbol)
    contains_ensembl <- any(grepl("^ENSG", preview_data_with_symbols$Gene_Symbol))
    
    cat("\n🎯 PREVIEW ANALYSIS:\n")
    cat("Contains expected gene symbols (BRCA2, BRCA1, etc.):", contains_symbols, "\n")
    cat("Still contains Ensembl IDs (ENSG...):", contains_ensembl, "\n")
    
    if (contains_symbols && !contains_ensembl) {
      cat("✅ SUCCESS: UI preview shows gene symbols correctly\n")
    } else if (!contains_symbols && contains_ensembl) {
      cat("❌ FAILURE: UI preview still shows Ensembl IDs - conversion failed\n")
    } else {
      cat("⚠️ MIXED: UI preview shows both symbols and Ensembl IDs\n")
    }
  }
  
  # STEP 5: Check conversion stats
  cat("\n📈 STEP 6: CONVERSION STATISTICS\n")
  cat("-", rep("-", 50), "\n")
  
  if (!is.null(local_values$conversion_stats)) {
    stats <- local_values$conversion_stats
    cat("Conversion attempted:", stats$attempted, "\n")
    cat("Total genes:", stats$total_genes, "\n") 
    cat("Converted count:", stats$converted_count, "\n")
    cat("Conversion rate:", stats$conversion_rate, "%\n")
    if (!is.null(stats$species)) {
      cat("Species used:", stats$species, "\n")
    }
  } else {
    cat("❌ No conversion stats available\n")
  }
  
} else {
  cat("❌ CRITICAL: Data processing failed:", processing_result$error, "\n")
}

cat("\n🎯 ROOT CAUSE ANALYSIS\n")
cat("=" , rep("=", 80), "\n")

cat("If the UI preview still shows Ensembl IDs, the issue is at one of these steps:\n")
cat("1. 🔄 Data processing: Aggregation worked but didn't preserve row structure\n")
cat("2. 🧬 Gene conversion: Conversion failed despite reporting success\n")
cat("3. 💾 Storage: Converted data not properly stored in reactive values\n")
cat("4. 📊 UI rendering: Preview pulling from wrong data source\n")
cat("5. 🐛 Code bugs: Recent edits introduced scope/variable issues\n")

cat("\n✅ Actual UI flow test completed - check results above\n")