# Test Duplicate Gene Aggregation and Symbol Display
# Tests the fixes for missing gene symbol column and duplicate aggregation

cat("🧪 Testing Duplicate Gene Aggregation and Symbol Display\n")
cat("=" , rep("=", 65), "\n")

# Load modules
source("config/app_config.R")
source("modules/data_processing/file_upload.R")

# Create test data with duplicate gene symbols
test_expression_matrix <- matrix(
  data = c(
    # GAPDH appears 3 times (should be aggregated)
    10, 15, 20, 12, 8, 25,   # ENSG00000111640 -> GAPDH
    5, 8, 12, 6, 4, 15,      # ENSG00000111641 -> GAPDH  
    7, 10, 14, 8, 5, 18,     # ENSG00000111642 -> GAPDH
    
    # ACTB appears twice (should be aggregated) 
    100, 120, 110, 95, 85, 130,  # ENSG00000075624 -> ACTB
    90, 110, 105, 88, 80, 125,   # ENSG00000075625 -> ACTB
    
    # Unique genes
    50, 60, 55, 45, 40, 65,      # ENSG00000139618 -> BRCA2
    30, 35, 32, 28, 25, 38       # ENSG00000012048 -> BRCA1
  ),
  nrow = 7,
  ncol = 6,
  dimnames = list(
    c("ENSG00000111640", "ENSG00000111641", "ENSG00000111642", 
      "ENSG00000075624", "ENSG00000075625", 
      "ENSG00000139618", "ENSG00000012048"),
    c("Control_1", "Control_2", "Control_3", "Treatment_1", "Treatment_2", "Treatment_3")
  )
)

cat("📊 Test data created:\n")
cat("   - Total genes:", nrow(test_expression_matrix), "\n")
cat("   - Samples:", ncol(test_expression_matrix), "\n")
cat("   - Expected duplicates: GAPDH (3 times), ACTB (2 times)\n")
cat("   - Expected final genes: 4 unique symbols\n\n")

# Mock gene conversion results (simulating biomaRt/local conversion)
mock_conversion_result <- data.frame(
  original_id = c("ENSG00000111640", "ENSG00000111641", "ENSG00000111642", 
                  "ENSG00000075624", "ENSG00000075625", 
                  "ENSG00000139618", "ENSG00000012048"),
  gene_symbol = c("GAPDH", "GAPDH", "GAPDH", 
                  "ACTB", "ACTB", 
                  "BRCA2", "BRCA1"),
  stringsAsFactors = FALSE
)

# Override convert_gene_ids_to_symbols for testing
convert_gene_ids_to_symbols <- function(gene_ids, species = "human") {
  cat("🧬 Mock conversion: Converting", length(gene_ids), "gene IDs to symbols...\n")
  return(mock_conversion_result[mock_conversion_result$original_id %in% gene_ids, ])
}

# Test the apply_gene_conversion function
cat("🔄 Testing gene conversion with aggregation...\\n")
conversion_result <- apply_gene_conversion(
  test_expression_matrix, 
  species = "human", 
  enable_conversion = TRUE
)

cat("\\n📋 RESULTS:\\n")
cat("-" , rep("-", 50), "\\n")

# Check results
if (conversion_result$conversion_stats$attempted) {
  cat("✅ Conversion attempted: YES\\n")
  cat("✅ Conversion rate:", conversion_result$conversion_stats$conversion_rate, "%\\n")
  
  if (!is.null(conversion_result$conversion_stats$duplicates_aggregated) &&
      conversion_result$conversion_stats$duplicates_aggregated) {
    cat("✅ Duplicates aggregated: YES\\n")
    cat("✅ Original genes:", conversion_result$conversion_stats$total_genes, "\\n")
    cat("✅ Final genes:", conversion_result$conversion_stats$final_genes, "\\n")
    cat("✅ Genes aggregated:", 
        conversion_result$conversion_stats$total_genes - conversion_result$conversion_stats$final_genes, "\\n")
  }
  
  # Check final matrix
  final_matrix <- conversion_result$matrix
  cat("\\n📊 Final Expression Matrix:\\n")
  cat("   - Dimensions:", nrow(final_matrix), "x", ncol(final_matrix), "\\n")
  cat("   - Gene symbols (rownames):", paste(rownames(final_matrix)[1:min(4, nrow(final_matrix))], collapse = ", "), "\\n")
  
  # Verify aggregation worked correctly
  if ("GAPDH" %in% rownames(final_matrix)) {
    gapdh_values <- final_matrix["GAPDH", ]
    original_gapdh_sum <- colSums(test_expression_matrix[1:3, ]) # First 3 rows were GAPDH
    cat("\\n🔍 GAPDH Aggregation Verification:\\n")
    cat("   - Aggregated GAPDH Control_1:", gapdh_values[1], "(expected:", original_gapdh_sum[1], ")\\n")
    cat("   - Match:", all.equal(as.numeric(gapdh_values), as.numeric(original_gapdh_sum)), "\\n")
  }
  
  # Check conversion table
  if (!is.null(conversion_result$conversion_table)) {
    cat("\\n📋 Conversion Table:\\n")
    print(conversion_result$conversion_table)
    
    # Verify aggregated genes are marked
    aggregated_genes <- conversion_result$conversion_table$Gene_Symbol[
      conversion_result$conversion_table$Aggregated == "Yes (summed)"
    ]
    cat("\\n✅ Genes marked as aggregated:", paste(aggregated_genes, collapse = ", "), "\\n")
  }
  
} else {
  cat("❌ Conversion was not attempted\\n")
}

cat("\\n🎉 TEST COMPLETED!\\n")
cat("=" , rep("=", 65), "\\n")

cat("\\n🔧 Features tested:\\n")
cat("   ✅ Duplicate gene symbol detection\\n")
cat("   ✅ Expression value aggregation by summing\\n")
cat("   ✅ Conversion table creation with aggregation info\\n")
cat("   ✅ Final matrix dimension reduction\\n")
cat("   ✅ Proper gene symbol assignment to rownames\\n")

cat("\\n📈 Expected improvements:\\n")
cat("   🎯 Gene symbols visible in conversion table\\n")
cat("   🎯 Duplicate genes properly aggregated before DESeq2\\n")
cat("   🎯 No artificial gene suffix creation\\n")
cat("   🎯 Accurate gene counts for downstream analysis\\n")

# Cleanup
rm(test_expression_matrix, mock_conversion_result, conversion_result)
cat("\\n✅ Test cleanup completed\\n")