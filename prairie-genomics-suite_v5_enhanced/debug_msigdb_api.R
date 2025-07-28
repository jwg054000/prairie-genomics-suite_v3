# Debug MSigDB API Issues
# Test different parameter combinations to find what works
# Author: Prairie Genomics Suite Development Team

cat("🔍 Debugging MSigDB API Issues\n")
cat("=", rep("=", 40), "\n")

# Load msigdbr
if (!require("msigdbr", quietly = TRUE)) {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  install.packages("msigdbr")
  library(msigdbr)
}

# Check available collections and species
cat("\n📋 Available Collections:\n")
collections <- msigdbr_collections()
print(collections)

cat("\n🐭 Available Species (mouse):\n")
species_info <- msigdbr_species()
mouse_info <- species_info[species_info$species_name == "Mus musculus", ]
print(mouse_info)

# Test different API approaches
cat("\n🧪 Testing Different API Approaches\n")
cat(rep("-", 45), "\n")

test_results <- list()

# Test 1: Modern API with collection only
cat("Test 1: collection='H', species='Mus musculus'\n")
test_results$test1 <- tryCatch({
  result <- msigdbr(species = "Mus musculus", collection = "H")
  cat("✅ Success:", nrow(result), "entries\n")
  list(success = TRUE, entries = nrow(result))
}, error = function(e) {
  cat("❌ Failed:", e$message, "\n")
  list(success = FALSE, error = e$message)
})

# Test 2: Modern API with collection and db_species
cat("\nTest 2: collection='H', species='Mus musculus', db_species='MM'\n")
test_results$test2 <- tryCatch({
  result <- msigdbr(species = "Mus musculus", collection = "H", db_species = "MM")
  cat("✅ Success:", nrow(result), "entries\n")
  list(success = TRUE, entries = nrow(result))
}, error = function(e) {
  cat("❌ Failed:", e$message, "\n")
  list(success = FALSE, error = e$message)
})

# Test 3: Old API with category only
cat("\nTest 3: category='H', species='Mus musculus'\n")
test_results$test3 <- tryCatch({
  suppressWarnings({
    result <- msigdbr(species = "Mus musculus", category = "H")
  })
  cat("✅ Success:", nrow(result), "entries\n")
  list(success = TRUE, entries = nrow(result))
}, error = function(e) {
  cat("❌ Failed:", e$message, "\n")
  list(success = FALSE, error = e$message)
})

# Test 4: Old API with category and db_species
cat("\nTest 4: category='H', species='Mus musculus', db_species='MM'\n")
test_results$test4 <- tryCatch({
  suppressWarnings({
    result <- msigdbr(species = "Mus musculus", category = "H", db_species = "MM")
  })
  cat("✅ Success:", nrow(result), "entries\n")
  list(success = TRUE, entries = nrow(result))
}, error = function(e) {
  cat("❌ Failed:", e$message, "\n")
  list(success = FALSE, error = e$message)
})

# Test 5: Check if db_species parameter is causing issues
cat("\nTest 5: collection='H', species='Mus musculus' (no db_species)\n")
test_results$test5 <- tryCatch({
  result <- msigdbr(species = "Mus musculus", collection = "H")
  cat("✅ Success:", nrow(result), "entries\n")
  
  # Check if this gives native mouse or human orthologs
  sample_genes <- head(unique(result$gene_symbol), 10)
  cat("📋 Sample genes:", paste(sample_genes, collapse = ", "), "\n")
  
  list(success = TRUE, entries = nrow(result), sample_genes = sample_genes)
}, error = function(e) {
  cat("❌ Failed:", e$message, "\n")
  list(success = FALSE, error = e$message)
})

# Summary and recommendation
cat("\n🎯 Test Summary and Recommendation\n")
cat("=", rep("=", 35), "\n")

working_tests <- names(test_results)[sapply(test_results, function(x) x$success)]
cat("✅ Working approaches:", paste(working_tests, collapse = ", "), "\n")

if ("test2" %in% working_tests) {
  cat("🏆 RECOMMENDED: Use collection='H' with db_species='MM' for native mouse\n")
  recommended_code <- '
msigdbr(species = "Mus musculus", collection = "H", db_species = "MM")
'
} else if ("test1" %in% working_tests || "test5" %in% working_tests) {
  cat("🏆 RECOMMENDED: Use collection='H' without db_species (may use orthologs)\n")
  recommended_code <- '
msigdbr(species = "Mus musculus", collection = "H")
'
} else if ("test3" %in% working_tests) {
  cat("🏆 FALLBACK: Use deprecated category='H' (shows warnings)\n")
  recommended_code <- '
msigdbr(species = "Mus musculus", category = "H")
'
} else {
  cat("❌ NO WORKING APPROACH FOUND - MSigDB may have issues\n")
  recommended_code <- "# No working approach found"
}

cat("\n💡 Recommended Code:\n")
cat(recommended_code, "\n")

cat("\n🔧 Ready to Update pathway_analysis.R\n")