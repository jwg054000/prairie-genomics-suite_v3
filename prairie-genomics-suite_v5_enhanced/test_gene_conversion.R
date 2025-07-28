# Test Script for Simple Gene Conversion
# Verifies that local gene conversion is working properly

cat("🧪 Testing Simple Gene Conversion System\n")
cat("========================================\n\n")

# Load the simple conversion functions
source("simple_gene_conversion.R")

# Test with sample mouse genes (from your data)
cat("📋 Test 1: Mouse Gene Conversion\n")
cat("--------------------------------\n")

test_mouse_genes <- c(
  "ENSMUSG00000037742",
  "ENSMUSG00000064351", 
  "ENSMUSG00000119584",
  "ENSMUSG00000025902",
  "ENSMUSG00000033845"
)

cat("Testing with mouse genes:", paste(test_mouse_genes, collapse = ", "), "\n\n")

mouse_result <- convert_genes_simple(test_mouse_genes, species = "mouse")

if (!is.null(mouse_result)) {
  cat("✅ Mouse conversion result:\n")
  print(mouse_result)
  
  converted_count <- sum(!is.na(mouse_result$gene_symbol) & mouse_result$gene_symbol != "")
  cat("Conversion rate:", round(100 * converted_count / length(test_mouse_genes), 1), "%\n")
} else {
  cat("❌ Mouse conversion failed\n")
}

cat("\n", rep("=", 50), "\n")

# Test with sample human genes
cat("📋 Test 2: Human Gene Conversion\n")
cat("--------------------------------\n")

test_human_genes <- c(
  "ENSG00000141510",  # TP53
  "ENSG00000012048",  # BRCA1
  "ENSG00000136997",  # MYC
  "ENSG00000133703",  # KRAS
  "ENSG00000171862"   # PTEN
)

cat("Testing with human genes:", paste(test_human_genes, collapse = ", "), "\n\n")

human_result <- convert_genes_simple(test_human_genes, species = "human")

if (!is.null(human_result)) {
  cat("✅ Human conversion result:\n")
  print(human_result)
  
  converted_count <- sum(!is.na(human_result$gene_symbol) & human_result$gene_symbol != "")
  cat("Conversion rate:", round(100 * converted_count / length(test_human_genes), 1), "%\n")
} else {
  cat("❌ Human conversion failed\n")
}

cat("\n", rep("=", 50), "\n")

# Test auto-detection
cat("📋 Test 3: Species Auto-Detection\n")
cat("----------------------------------\n")

mixed_genes <- c(test_mouse_genes[1:2], test_human_genes[1:2])
cat("Testing mixed genes:", paste(mixed_genes, collapse = ", "), "\n")

detected_species <- detect_species_from_genes(mixed_genes)
cat("Auto-detected species:", detected_species, "\n")

mouse_detected <- detect_species_from_genes(test_mouse_genes)
human_detected <- detect_species_from_genes(test_human_genes)

if (mouse_detected == "mouse" && human_detected == "human") {
  cat("✅ Species auto-detection working correctly\n")
} else {
  cat("❌ Species auto-detection failed\n")
  cat("   Mouse genes detected as:", mouse_detected, "\n")
  cat("   Human genes detected as:", human_detected, "\n")
}

cat("\n📊 SUMMARY\n")
cat("==========\n")

if (exists("mouse_result") && exists("human_result")) {
  mouse_success <- !is.null(mouse_result) && sum(!is.na(mouse_result$gene_symbol)) > 0
  human_success <- !is.null(human_result) && sum(!is.na(human_result$gene_symbol)) > 0
  
  if (mouse_success && human_success) {
    cat("🎉 SUCCESS: Both mouse and human gene conversion working!\n")
    cat("✅ The app should now convert your mouse genes properly\n")
  } else if (mouse_success) {
    cat("✅ Mouse conversion working, human needs org.Hs.eg.db package\n")
  } else if (human_success) {
    cat("✅ Human conversion working, mouse needs org.Mm.eg.db package\n")
  } else {
    cat("❌ Both conversions failed - missing Bioconductor packages\n")
    cat("💡 Run: install_gene_conversion_packages()\n")
  }
} else {
  cat("❌ Test failed to complete\n")
}

cat("\n🔧 Next Steps:\n")
cat("1. If packages are missing, run: install_gene_conversion_packages()\n")
cat("2. If successful, restart your app with run_app.R\n")
cat("3. Upload your data and enable gene conversion\n")
cat("4. The system will now use fast local conversion instead of BioMart\n")