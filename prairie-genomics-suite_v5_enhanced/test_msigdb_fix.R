# Test MSigDB Fix
# Verify the updated get_fgsea_gene_sets function works

cat("🧪 Testing MSigDB Fix\n")
cat("=", rep("=", 25), "\n")

# Load the updated pathway analysis
source("pathway_analysis.R")

cat("\n🔬 Testing get_fgsea_gene_sets function...\n")

# Test with mouse Hallmark gene sets
mouse_hallmark <- tryCatch({
  get_fgsea_gene_sets("mouse", "H")
}, error = function(e) {
  cat("❌ Test failed:", e$message, "\n")
  NULL
})

if (!is.null(mouse_hallmark) && length(mouse_hallmark) > 0) {
  cat("✅ MSigDB function works!\n")
  cat("📊 Retrieved", length(mouse_hallmark), "Hallmark pathways\n")
  
  # Show some pathway names
  pathway_names <- names(mouse_hallmark)[1:min(3, length(mouse_hallmark))]
  cat("🛤️ Example pathways:", paste(pathway_names, collapse = ", "), "\n")
  
  # Show some genes from first pathway
  first_pathway_genes <- mouse_hallmark[[1]][1:min(5, length(mouse_hallmark[[1]]))]
  cat("🧬 Genes in first pathway:", paste(first_pathway_genes, collapse = ", "), "\n")
  
} else {
  cat("❌ MSigDB function still has issues\n")
}

cat("\n🎯 MSigDB Fix Test Complete!\n")