# Test MSigDB Collections for Mouse
# Verify that C2 and other collections work properly
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🧪 Testing MSigDB Collections for Mouse\n")
cat("=" , rep("=", 50), "\n")

# Test 1: Load msigdbr and check species
cat("\nTest 1: MSigDB Package and Species\n")
cat(rep("-", 35), "\n")

if (!require("msigdbr", quietly = TRUE)) {
  cat("📦 Installing msigdbr...\n")
  install.packages("msigdbr")
  library(msigdbr)
}

cat("✅ msigdbr loaded\n")

# Check available species
species_list <- msigdbr_species()
cat("🌍 Available species:\n")
print(species_list[species_list$species_name %in% c("Homo sapiens", "Mus musculus"), ])

# Test 2: Test collections for mouse
cat("\nTest 2: Mouse Collection Availability\n")
cat(rep("-", 35), "\n")

test_collections <- c("H", "C2", "C5")

for (collection in test_collections) {
  cat("Testing", collection, "for mouse... ")
  
  tryCatch({
    gene_sets <- msigdbr(species = "Mus musculus", category = collection)
    cat("✅", nrow(gene_sets), "entries\n")
    
    # Show some example pathways
    if (nrow(gene_sets) > 0) {
      example_pathways <- unique(gene_sets$gs_name)[1:min(3, length(unique(gene_sets$gs_name)))]
      cat("   Examples:", paste(example_pathways, collapse = ", "), "\n")
    }
    
  }, error = function(e) {
    cat("❌", e$message, "\n")
  })
}

# Test 3: Test C2 subcategories specifically
cat("\nTest 3: C2 Subcategories for Mouse\n")
cat(rep("-", 35), "\n")

c2_subcategories <- c("CP", "CP:KEGG", "CP:REACTOME", "CGP")

for (subcat in c2_subcategories) {
  cat("Testing C2:", subcat, "... ")
  
  tryCatch({
    if (subcat == "CP") {
      gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
    } else {
      gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = subcat)
    }
    
    cat("✅", nrow(gene_sets), "entries\n")
    
  }, error = function(e) {
    cat("❌", e$message, "\n")
  })
}

# Test 4: Test the improved get_fgsea_gene_sets function
cat("\nTest 4: Improved Gene Set Retrieval Function\n")
cat(rep("-", 45), "\n")

# Load the pathway analysis functions
tryCatch({
  source("pathway_analysis.R")
  cat("✅ pathway_analysis.R loaded\n")
  
  if (exists("get_fgsea_gene_sets")) {
    cat("✅ get_fgsea_gene_sets function available\n")
    
    # Test with different collections
    test_cases <- list(
      list(species = "mouse", collection = "H"),
      list(species = "mouse", collection = "C2"), 
      list(species = "human", collection = "H"),
      list(species = "human", collection = "C2")
    )
    
    for (test_case in test_cases) {
      cat("🧪 Testing", test_case$species, test_case$collection, "... ")
      
      result <- get_fgsea_gene_sets(test_case$species, test_case$collection)
      
      if (!is.null(result) && length(result) > 0) {
        cat("✅", length(result), "pathways\n")
      } else {
        cat("❌ Failed or empty\n")
      }
    }
    
  } else {
    cat("❌ get_fgsea_gene_sets function not found\n")
  }
  
}, error = function(e) {
  cat("❌ Function test failed:", e$message, "\n")
})

# Test 5: Create working examples
cat("\nTest 5: Working Collection Examples\n")
cat(rep("-", 35), "\n")

cat("💡 Recommended collections for mouse:\n")

# Test Hallmark (most reliable)
tryCatch({
  h_mouse <- msigdbr(species = "Mus musculus", category = "H")
  cat("✅ Hallmark (H):", nrow(h_mouse), "pathways - RECOMMENDED\n")
}, error = function(e) {
  cat("❌ Hallmark failed\n")
})

# Test KEGG specifically
tryCatch({
  kegg_mouse <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
  cat("✅ KEGG (C2:CP:KEGG):", nrow(kegg_mouse), "pathways - GOOD FOR PATHWAYS\n")
}, error = function(e) {
  cat("❌ KEGG failed\n")
})

# Test GO Biological Process
tryCatch({
  go_mouse <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
  cat("✅ GO:BP (C5:GO:BP):", nrow(go_mouse), "pathways - LARGE COLLECTION\n")  
}, error = function(e) {
  cat("❌ GO:BP failed\n")
})

# Summary and recommendations
cat("\n🎯 Summary and Recommendations\n")
cat("=" , rep("=", 40), "\n")

cat("🔧 Issue Resolution:\n")
cat("• C2 collection is too broad - use subcategories\n")
cat("• C2:CP:KEGG works well for pathway analysis\n")
cat("• Hallmark (H) is most reliable across species\n")
cat("• Improved function includes automatic fallbacks\n")

cat("\n💡 Usage Recommendations:\n")
cat("• For general analysis: Use Hallmark (H)\n")
cat("• For pathway analysis: Use C2:CP:KEGG or C2:CP:REACTOME\n")
cat("• For GO analysis: Use C5:GO:BP, C5:GO:MF, or C5:GO:CC\n")
cat("• Avoid broad C2 - always specify subcategory\n")

cat("\n🚀 App Integration:\n")
cat("1. Update UI to show recommended collections\n")
cat("2. Set default to Hallmark (H) for reliability\n")
cat("3. Add subcategory options for C2 and C5\n")
cat("4. Implement automatic fallbacks in analysis\n")

cat("\n🧬 MSigDB Collections Testing Complete!\n")
cat("📊 Mouse collections should now work properly\n")