
# Quick test of both fixes
cat("🧪 Testing Both Fixes\n")

# Load fixes
source("complete_fixes_both_errors.R")

# Test 1: MSigDB function with category parameter
cat("\n1. Testing MSigDB function...\n")
if (exists("get_fgsea_gene_sets")) {
  mouse_pathways <- get_fgsea_gene_sets("mouse", "H")
  if (!is.null(mouse_pathways) && length(mouse_pathways) > 0) {
    cat("✅ MSigDB function works:", length(mouse_pathways), "pathways\n")
  } else {
    cat("❌ MSigDB function failed\n")
  }
}

# Test 2: GSEA preparation with mock data
cat("\n2. Testing GSEA preparation...\n")
if (exists("prepare_gene_list_gsea")) {
  # Create mock data without baseMean column
  mock_data <- data.frame(
    log2FoldChange = rnorm(100, 0, 1),
    pvalue = runif(100, 0.001, 0.1),
    padj = runif(100, 0.001, 0.1),
    stringsAsFactors = FALSE
  )
  rownames(mock_data) <- paste0("GENE_", 1:100)
  
  gene_list <- prepare_gene_list_gsea(mock_data, "human", max_genes = 50)
  if (!is.null(gene_list) && length(gene_list) > 0) {
    cat("✅ GSEA preparation works:", length(gene_list), "genes\n")
  } else {
    cat("❌ GSEA preparation failed\n")
  }
}

cat("\n🎉 Both fixes tested!\n")

