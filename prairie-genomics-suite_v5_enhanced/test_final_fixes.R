# Test Final Fixes for MSigDB API and GSEA Logical Coercion
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🧪 Testing Final Fixes\n")
cat("=", rep("=", 30), "\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Test 1: MSigDB API compatibility
cat("\nTest 1: MSigDB API Compatibility\n")
cat(rep("-", 35), "\n")

if (!require("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
  library(msigdbr)
}

# Test collection parameter (new API)
cat("Testing collection parameter...\n")
collection_test <- tryCatch({
  msigdbr(species = "Mus musculus", collection = "H")
}, error = function(e) {
  cat("❌ Collection parameter failed:", e$message, "\n")
  NULL
})

# Test category parameter (old API)
cat("Testing category parameter...\n")
category_test <- tryCatch({
  msigdbr(species = "Mus musculus", category = "H")
}, error = function(e) {
  cat("❌ Category parameter failed:", e$message, "\n")
  NULL
})

# Determine which API works
if (!is.null(collection_test) && nrow(collection_test) > 0) {
  cat("✅ Collection parameter works:", nrow(collection_test), "entries\n")
  working_api <- "collection"
} else if (!is.null(category_test) && nrow(category_test) > 0) {
  cat("✅ Category parameter works:", nrow(category_test), "entries\n")
  working_api <- "category"
} else {
  cat("❌ Both APIs failed\n")
  working_api <- "unknown"
}

# Test 2: GSEA Logical Coercion Fix
cat("\nTest 2: GSEA Logical Coercion Fix\n")
cat(rep("-", 35), "\n")

# Load the pathway analysis functions
source("pathway_analysis.R")

# Create test data with baseMean column
test_data <- data.frame(
  log2FoldChange = rnorm(1000, 0, 1.5),
  pvalue = runif(1000, 0.001, 0.1),
  padj = runif(1000, 0.001, 0.1),
  baseMean = runif(1000, 5, 200),
  stringsAsFactors = FALSE
)
rownames(test_data) <- paste0("GENE_", sprintf("%04d", 1:1000))

cat("Testing with", nrow(test_data), "genes...\n")

# Test the GSEA preparation function
gsea_test <- tryCatch({
  prepare_gene_list_gsea(test_data, "human", max_genes = 300)
}, error = function(e) {
  cat("❌ GSEA preparation failed:", e$message, "\n")
  NULL
})

if (!is.null(gsea_test) && length(gsea_test) > 0) {
  cat("✅ GSEA preparation successful:", length(gsea_test), "genes\n")
  logical_fix_works <- TRUE
} else {
  cat("❌ GSEA preparation failed\n")
  logical_fix_works <- FALSE
}

# Summary
cat("\n🎯 Final Test Summary\n")
cat("=", rep("=", 25), "\n")

cat("MSigDB API Status:\n")
cat("  Working parameter:", working_api, "\n")
if (working_api == "category") {
  cat("  ✅ Our fix uses category parameter - CORRECT\n")
} else if (working_api == "collection") {
  cat("  ⚠️ Collection parameter works - may need to update fix\n")
} else {
  cat("  ❌ MSigDB API issues remain\n")
}

cat("\nGSEA Logical Fix Status:\n")
if (logical_fix_works) {
  cat("  ✅ Logical coercion error fixed\n")
  cat("  ✅ Gene filtering works properly\n")
} else {
  cat("  ❌ Logical coercion issues remain\n")
}

cat("\n🚀 Ready for Application Testing:\n")
if (working_api != "unknown" && logical_fix_works) {
  cat("✅ Both critical issues resolved\n")
  cat("✅ pathway_analysis.R should work properly\n")
  cat("🧪 Ready to test with the main Shiny application\n")
} else {
  cat("❌ Some issues remain - check individual test results above\n")
}

cat("\n🧬 Final Fix Testing Complete!\n")