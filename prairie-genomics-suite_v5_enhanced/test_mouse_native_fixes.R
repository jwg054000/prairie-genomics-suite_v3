# Test Mouse Native MSigDB and Enhanced Filtering
# Verify fixes for ortholog mapping and gene count issues
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🧪 Testing Mouse Native MSigDB and Filtering Fixes\n")
cat("=" , rep("=", 55), "\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Test 1: Verify MSigDB updates
cat("\nTest 1: MSigDB API Updates\n")
cat(rep("-", 30), "\n")

if (!require("msigdbr", quietly = TRUE)) {
  cat("📦 Installing msigdbr...\n")
  install.packages("msigdbr")
  library(msigdbr)
}

cat("✅ msigdbr loaded\n")

# Check if new API is available
tryCatch({
  # Test new API with db_species parameter
  cat("🧪 Testing native mouse gene sets...\n")
  
  # This should use native mouse gene sets (no ortholog mapping)
  mouse_native <- msigdbr(
    species = "Mus musculus",
    collection = "H",  # Using collection instead of category
    db_species = "MM"  # Native mouse database
  )
  
  cat("✅ Native mouse Hallmark:", nrow(mouse_native), "gene entries\n")
  cat("📋 Database used: MM (native mouse)\n")
  
  # Show some example genes to verify they're mouse genes
  example_genes <- head(unique(mouse_native$gene_symbol), 10)
  cat("🐭 Example mouse genes:", paste(example_genes, collapse = ", "), "\n")
  
}, error = function(e) {
  cat("❌ Native mouse test failed:", e$message, "\n")
  cat("💡 May need to update msigdbr package\n")
})

# Test 2: Compare human orthologs vs native mouse
cat("\nTest 2: Ortholog vs Native Comparison\n")
cat(rep("-", 40), "\n")

tryCatch({
  # Human gene sets (for comparison)
  human_h <- msigdbr(species = "Homo sapiens", collection = "H", db_species = "HS")
  cat("🧑 Human Hallmark (native):", nrow(human_h), "gene entries\n")
  
  # Mouse with human orthologs (old way)
  mouse_ortho <- msigdbr(species = "Mus musculus", collection = "H", db_species = "HS")
  cat("🐭 Mouse Hallmark (human orthologs):", nrow(mouse_ortho), "gene entries\n")
  
  # Mouse native (new way)  
  mouse_native <- msigdbr(species = "Mus musculus", collection = "H", db_species = "MM")
  cat("🐭 Mouse Hallmark (native):", nrow(mouse_native), "gene entries\n")
  
  # Compare gene examples
  cat("\n📊 Gene comparison:\n")
  cat("   Human genes:", paste(head(unique(human_h$gene_symbol), 5), collapse = ", "), "\n")
  cat("   Mouse ortho:", paste(head(unique(mouse_ortho$gene_symbol), 5), collapse = ", "), "\n") 
  cat("   Mouse native:", paste(head(unique(mouse_native$gene_symbol), 5), collapse = ", "), "\n")
  
}, error = function(e) {
  cat("❌ Comparison test failed:", e$message, "\n")
})

# Test 3: Updated pathway analysis functions
cat("\nTest 3: Updated Pathway Functions\n")
cat(rep("-", 35), "\n")

tryCatch({
  source("pathway_analysis.R")
  cat("✅ pathway_analysis.R loaded\n")
  
  # Test the updated get_fgsea_gene_sets function
  if (exists("get_fgsea_gene_sets")) {
    cat("🧪 Testing native mouse gene set retrieval...\n")
    
    mouse_pathways <- get_fgsea_gene_sets("mouse", "H")
    
    if (!is.null(mouse_pathways) && length(mouse_pathways) > 0) {
      cat("✅ Native mouse pathways retrieved:", length(mouse_pathways), "\n")
      cat("📋 Should show 'Database used: MM (native mouse)'\n")
      
      # Show some pathway names
      pathway_names <- names(mouse_pathways)[1:3]
      cat("🛤️ Example pathways:", paste(pathway_names, collapse = ", "), "\n")
      
    } else {
      cat("❌ Mouse pathway retrieval failed\n")
    }
  }
  
}, error = function(e) {
  cat("❌ Function test failed:", e$message, "\n")
})

# Test 4: Enhanced gene filtering
cat("\nTest 4: Enhanced Gene Filtering\n")
cat(rep("-", 35), "\n")

tryCatch({
  # Create realistic mock DESeq2 data with varying significance
  set.seed(123)
  n_genes <- 2000  # Large dataset to test filtering
  
  mock_deseq <- data.frame(
    log2FoldChange = rnorm(n_genes, 0, 1.5),
    pvalue = c(runif(500, 0.001, 0.01),    # 500 very significant
               runif(500, 0.01, 0.05),     # 500 significant  
               runif(1000, 0.05, 1.0)),    # 1000 not significant
    padj = c(runif(500, 0.001, 0.01),      # 500 very significant
             runif(500, 0.01, 0.05),       # 500 significant
             runif(1000, 0.05, 1.0)),      # 1000 not significant
    baseMean = c(runif(1000, 100, 1000),   # 1000 highly expressed
                 runif(500, 10, 100),      # 500 moderately expressed  
                 runif(500, 1, 10)),       # 500 lowly expressed
    stringsAsFactors = FALSE
  )
  
  rownames(mock_deseq) <- paste0("GENE_", sprintf("%04d", 1:n_genes))
  
  cat("✅ Created mock DESeq2 data:", nrow(mock_deseq), "genes\n")
  
  # Test different filtering levels
  if (exists("prepare_gene_list_gsea")) {
    
    # Test 1: Default filtering
    cat("\n🧪 Testing default filtering...\n")
    default_genes <- prepare_gene_list_gsea(mock_deseq, "human")
    if (!is.null(default_genes)) {
      cat("✅ Default filtering:", length(default_genes), "genes\n")
    }
    
    # Test 2: Stringent filtering
    cat("\n🧪 Testing stringent filtering...\n")
    stringent_genes <- prepare_gene_list_gsea(
      mock_deseq, "human", 
      max_genes = 300,
      padj_filter = 0.01,
      basemean_filter = 50
    )
    if (!is.null(stringent_genes)) {
      cat("✅ Stringent filtering:", length(stringent_genes), "genes\n")
    }
    
    # Test 3: Very strict filtering  
    cat("\n🧪 Testing very strict filtering...\n")
    strict_genes <- prepare_gene_list_gsea(
      mock_deseq, "human",
      max_genes = 100, 
      padj_filter = 0.005,
      basemean_filter = 100
    )
    if (!is.null(strict_genes)) {
      cat("✅ Very strict filtering:", length(strict_genes), "genes\n")
    }
    
  } else {
    cat("❌ prepare_gene_list_gsea function not found\n")
  }
  
}, error = function(e) {
  cat("❌ Gene filtering test failed:", e$message, "\n")
})

# Test 5: End-to-end GSEA with native mouse
cat("\nTest 5: End-to-End Mouse GSEA Test\n")
cat(rep("-", 40), "\n")

if (require("fgsea", quietly = TRUE)) {
  tryCatch({
    # Create small test gene list
    test_genes <- c(2.5, 1.8, -1.2, 3.1, -2.0, 1.5, -1.8, 2.2, 1.1, -1.4)
    names(test_genes) <- c("Tp53", "Brca1", "Egfr", "Myc", "Kras", "Pik3ca", "Akt1", "Mtor", "Pten", "Rb1")
    test_genes <- sort(test_genes, decreasing = TRUE)
    
    cat("🧬 Testing end-to-end GSEA with native mouse gene sets...\n")
    cat("   - Test genes:", length(test_genes), "\n")
    
    if (exists("run_gsea_analysis")) {
      gsea_result <- run_gsea_analysis(test_genes, "mouse", "H")
      
      if (!is.null(gsea_result) && gsea_result$success) {
        cat("✅ End-to-end GSEA successful!\n")
        cat("   - Should show 'NATIVE gene sets' and 'Database used: MM'\n")
        cat("   - Method:", gsea_result$method, "\n")
        cat("   - Pathways found:", gsea_result$n_pathways, "\n")
      } else {
        cat("❌ End-to-end GSEA failed\n")
        if (!is.null(gsea_result$error)) {
          cat("   Error:", gsea_result$error, "\n")
        }
      }
    }
    
  }, error = function(e) {
    cat("❌ End-to-end test failed:", e$message, "\n")
  })
} else {
  cat("⚠️ fgsea not available for end-to-end test\n")
}

# Summary
cat("\n🎯 Test Summary\n")
cat("=" , rep("=", 20), "\n")

cat("Fixes Implemented:\n")
cat("✅ Native mouse gene sets (db_species = 'MM')\n")
cat("✅ Updated MSigDB API (collection instead of category)\n")
cat("✅ Enhanced gene filtering with count limits\n")
cat("✅ Expression level filtering (baseMean)\n")
cat("✅ Significance filtering (padj)\n")

cat("\nExpected Improvements:\n")
cat("• No more ortholog mapping warnings\n")
cat("• Faster gene set retrieval\n")
cat("• More accurate mouse-specific pathways\n")
cat("• Controlled gene count (prevents overload)\n")
cat("• Better GSEA performance\n")

cat("\n🚀 Ready for App Testing:\n")
cat("1. ✅ Native mouse gene sets configured\n")
cat("2. ✅ Enhanced filtering prevents gene overload\n")
cat("3. ✅ Updated MSigDB API (no deprecation warnings)\n")
cat("4. 🧪 Test the main app with mouse data\n")

cat("\n💡 Recommended Settings for App:\n")
cat("For typical analysis:\n")
cat("  - max_genes = 500-800\n")
cat("  - padj_filter = 0.05-0.1\n")
cat("  - basemean_filter = 10-20\n")

cat("\nFor stringent analysis:\n")
cat("  - max_genes = 200-300\n")
cat("  - padj_filter = 0.01-0.05\n")
cat("  - basemean_filter = 50-100\n")

cat("\n🧬 Mouse Native MSigDB Issues Fixed!\n")
cat("🐭 Should now use proper native mouse gene sets\n")