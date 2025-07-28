# Debug script to check actual column names in clusterProfiler results
suppressMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

cat("🔍 DEBUGGING: Checking actual clusterProfiler column names\n")
cat("=======================================================\n")

# Test with known mouse genes
test_genes <- c('11545', '74778', '17869', '16590', '12043', '13982', '16193', '11593', '20926', '75560')

cat("📊 Testing GO enrichment...\n")
tryCatch({
  go_result <- enrichGO(gene = test_genes, 
                       OrgDb = org.Mm.eg.db, 
                       ont = 'BP', 
                       pvalueCutoff = 0.1, 
                       qvalueCutoff = 0.2)
  
  if (!is.null(go_result) && nrow(go_result@result) > 0) {
    go_df <- as.data.frame(go_result@result)
    cat("✅ GO results found:", nrow(go_df), "pathways\n")
    cat("📋 GO result columns:\n")
    for (i in 1:length(colnames(go_df))) {
      cat("  ", i, ":", colnames(go_df)[i], "\n")
    }
    
    # Check specific columns I'm trying to filter on
    cat("\n🔍 Column existence check:\n")
    cat("  - 'p.adjust' exists:", "p.adjust" %in% colnames(go_df), "\n")
    cat("  - 'Count' exists:", "Count" %in% colnames(go_df), "\n")
    cat("  - 'pvalue' exists:", "pvalue" %in% colnames(go_df), "\n")
    
    # Show sample data
    cat("\n📊 Sample data (first 2 rows, key columns):\n")
    if ("p.adjust" %in% colnames(go_df) && "Count" %in% colnames(go_df)) {
      sample_df <- go_df[1:min(2, nrow(go_df)), c("ID", "Description", "pvalue", "p.adjust", "Count")]
      print(sample_df)
    }
    
  } else {
    cat("❌ No GO results found\n")
  }
}, error = function(e) {
  cat("❌ GO enrichment failed:", e$message, "\n")
})

cat("\n📊 Testing KEGG enrichment...\n")
tryCatch({
  kegg_result <- enrichKEGG(gene = test_genes, 
                           organism = 'mmu', 
                           pvalueCutoff = 0.1, 
                           qvalueCutoff = 0.2)
  
  if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
    kegg_df <- as.data.frame(kegg_result@result)
    cat("✅ KEGG results found:", nrow(kegg_df), "pathways\n")
    cat("📋 KEGG result columns:\n")
    for (i in 1:length(colnames(kegg_df))) {
      cat("  ", i, ":", colnames(kegg_df)[i], "\n")
    }
    
    # Check specific columns 
    cat("\n🔍 Column existence check:\n")
    cat("  - 'p.adjust' exists:", "p.adjust" %in% colnames(kegg_df), "\n")
    cat("  - 'Count' exists:", "Count" %in% colnames(kegg_df), "\n")
    cat("  - 'pvalue' exists:", "pvalue" %in% colnames(kegg_df), "\n")
    
  } else {
    cat("❌ No KEGG results found\n")
  }
}, error = function(e) {
  cat("❌ KEGG enrichment failed:", e$message, "\n")
})

cat("\n✅ Column name debugging complete\n")