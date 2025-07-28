# ULTRA DEBUG: Gene Conversion Pipeline Analysis
# This will diagnose exactly where the gene conversion is failing

cat("🔬 ULTRA DEBUG: Gene Conversion Pipeline Analysis\n")
cat("=" , rep("=", 70), "\n")

# Load modules
source("config/app_config.R")
source("modules/data_processing/file_upload.R")

# Test with various real-world gene ID formats
test_gene_ids <- c(
  # Standard Ensembl human gene IDs
  "ENSG00000139618",      # BRCA2
  "ENSG00000012048",      # BRCA1  
  "ENSG00000111640",      # GAPDH
  "ENSG00000075624",      # ACTB
  
  # Versioned Ensembl IDs (common in real data)
  "ENSG00000139618.13",   # BRCA2 with version
  "ENSG00000012048.22",   # BRCA1 with version
  
  # Different species
  "ENSMUSG00000041147",   # Mouse gene
  
  # Non-Ensembl IDs (should be skipped)
  "GAPDH",
  "ACTB", 
  "NM_005228"             # RefSeq ID
)

cat("🧪 Testing with real gene ID formats:\n")
for (i in 1:length(test_gene_ids)) {
  cat("   ", i, ":", test_gene_ids[i], "\n")
}

# DIAGNOSTIC 1: Pattern Matching Analysis
cat("\n🔍 DIAGNOSTIC 1: Ensembl Pattern Matching\n")
cat("-" , rep("-", 50), "\n")

current_pattern <- "^ENS[A-Z]*G\\d{11}"
cat("Current pattern:", current_pattern, "\n")

for (gene_id in test_gene_ids) {
  matches <- grepl(current_pattern, gene_id)
  cat("   ", gene_id, "→ Matches:", matches, "\n")
}

# Test with more flexible pattern
flexible_pattern <- "^ENS[A-Z]*G\\d{11}(\\.\\d+)?$"  # Allow version numbers
cat("\nTesting flexible pattern:", flexible_pattern, "\n")

for (gene_id in test_gene_ids) {
  matches <- grepl(flexible_pattern, gene_id)
  cat("   ", gene_id, "→ Matches:", matches, "\n")
}

# DIAGNOSTIC 2: Biomart Availability Check
cat("\n🔍 DIAGNOSTIC 2: Gene Conversion Method Availability\n")
cat("-" , rep("-", 50), "\n")

cat("biomart_available:", exists("biomart_available") && biomart_available, "\n")
cat("orgdb_available:", exists("orgdb_available") && orgdb_available, "\n")

if (exists("biomart_available") && biomart_available) {
  cat("✅ biomaRt is loaded and available\n")
} else {
  cat("❌ biomaRt is NOT available\n")
}

if (exists("orgdb_available") && orgdb_available) {
  cat("✅ org.Hs.eg.db is loaded and available\n")
} else {
  cat("❌ org.Hs.eg.db is NOT available\n")  
}

# DIAGNOSTIC 3: Test Conversion Functions Directly
cat("\n🔍 DIAGNOSTIC 3: Direct Conversion Function Testing\n")
cat("-" , rep("-", 50), "\n")

# Test just the human Ensembl IDs
human_ensembl_ids <- c("ENSG00000139618", "ENSG00000012048", "ENSG00000111640", "ENSG00000075624")

cat("Testing convert_gene_ids_to_symbols directly...\n")
tryCatch({
  conversion_result <- convert_gene_ids_to_symbols(human_ensembl_ids, species = "human")
  cat("✅ Conversion function executed successfully\n")
  cat("Results:\n")
  print(conversion_result)
}, error = function(e) {
  cat("❌ Conversion function failed:", e$message, "\n")
})

# DIAGNOSTIC 4: Test apply_gene_conversion Function
cat("\n🔍 DIAGNOSTIC 4: Testing apply_gene_conversion Function\n")
cat("-" , rep("-", 50), "\n")

# Create test matrix with human Ensembl IDs
test_matrix <- matrix(
  data = rnorm(20),
  nrow = 4,
  ncol = 5,
  dimnames = list(
    human_ensembl_ids,
    c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
  )
)

cat("Test matrix created with gene IDs:", paste(rownames(test_matrix), collapse = ", "), "\n")

# Test the apply_gene_conversion function
cat("Testing apply_gene_conversion...\n")
tryCatch({
  gene_conversion_result <- apply_gene_conversion(test_matrix, species = "human", enable_conversion = TRUE)
  
  cat("✅ apply_gene_conversion executed\n")
  cat("Conversion attempted:", gene_conversion_result$conversion_stats$attempted, "\n")
  
  if (gene_conversion_result$conversion_stats$attempted) {
    cat("Conversion rate:", gene_conversion_result$conversion_stats$conversion_rate, "%\n")
    cat("Final matrix rownames:", paste(rownames(gene_conversion_result$matrix), collapse = ", "), "\n")
    
    if (!is.null(gene_conversion_result$conversion_table)) {
      cat("Conversion table:\n")
      print(gene_conversion_result$conversion_table)
    } else {
      cat("❌ No conversion table created\n")
    }
  } else {
    cat("❌ Conversion was not attempted\n")
    if (!is.null(gene_conversion_result$conversion_stats$reason)) {
      cat("Reason:", gene_conversion_result$conversion_stats$reason, "\n")
    }
  }
  
}, error = function(e) {
  cat("❌ apply_gene_conversion failed:", e$message, "\n")
})

# DIAGNOSTIC 5: Test Full Pipeline
cat("\n🔍 DIAGNOSTIC 5: Full Pipeline Test\n")
cat("-" , rep("-", 50), "\n")

# Create full test data as would come from file upload
full_test_data <- data.frame(
  Gene = human_ensembl_ids,
  Control_1 = c(100, 200, 300, 400),
  Control_2 = c(110, 210, 310, 410),
  Treatment_1 = c(120, 220, 320, 420),
  Treatment_2 = c(130, 230, 330, 430),
  stringsAsFactors = FALSE
)

cat("Testing full pipeline from raw data...\n")
tryCatch({
  # Test processing
  processing_result <- process_expression_data(full_test_data)
  
  if (processing_result$success) {
    cat("✅ Data processing successful\n")
    processed_matrix := processing_result$data
    
    # Test conversion on processed data
    final_conversion_result := apply_gene_conversion(processed_matrix, species = "human", enable_conversion = TRUE)
    
    cat("Final conversion stats:\n")
    cat("  - Attempted:", final_conversion_result$conversion_stats$attempted, "\n")
    cat("  - Rate:", final_conversion_result$conversion_stats$conversion_rate, "%\n")
    cat("  - Final rownames:", paste(rownames(final_conversion_result$matrix), collapse = ", "), "\n")
    
  } else {
    cat("❌ Data processing failed:", processing_result$error, "\n")
  }
  
}, error = function(e) {
  cat("❌ Full pipeline failed:", e$message, "\n")
})

# DIAGNOSTIC 6: Check Current Configuration
cat("\n🔍 DIAGNOSTIC 6: Current Configuration Check\n")
cat("-" , rep("-", 50), "\n")

cat("Current working directory:", getwd(), "\n")
cat("R version:", R.version.string, "\n")

# Check if required packages are actually loaded
required_packages <- c("biomaRt", "org.Hs.eg.db", "AnnotationDbi")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✅", pkg, "is available\n")
  } else {
    cat("❌", pkg, "is NOT available\n")
  }
}

cat("\n🎯 SUMMARY OF FINDINGS:\n")
cat("=" , rep("=", 70), "\n")
cat("This diagnostic should reveal:\n")
cat("1. Whether gene ID pattern matching is working\n")
cat("2. Whether conversion functions are available\n") 
cat("3. Whether conversion is happening but not displaying\n")
cat("4. Where exactly in the pipeline the failure occurs\n")

cat("\n✅ Diagnostic completed - check results above for issues\n")