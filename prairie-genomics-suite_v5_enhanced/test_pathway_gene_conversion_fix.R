# Test Pathway Gene Conversion Fix
# Verify that the gene ID detection and conversion works correctly
#
# Author: Prairie Genomics Team
# Date: January 27, 2025

cat("🧪 TESTING PATHWAY GENE CONVERSION FIX\n")
cat("======================================\n\n")

# Load required modules
source("pathway_analysis.R")

# Test the fixed prepare_gene_list_ora function
test_gene_id_detection_and_conversion <- function() {
  cat("🔍 TESTING GENE ID DETECTION AND CONVERSION\n")
  cat("==========================================\n\n")
  
  # Test 1: Mock DESeq2 results with gene symbols (like from our pre-conversion)
  cat("📊 TEST 1: Gene symbols (post pre-DESeq2 conversion)\n")
  cat("---------------------------------------------------\n")
  
  set.seed(123)
  n_genes <- 100
  
  # Create mock results with gene symbols as rownames
  gene_symbols <- c("Tp53", "Brca1", "Myc", "Kras", "Pik3ca", "Akt1", "Mtor", "Rb1", "Cdk4", "Ccnd1",
                   "Egfr", "Vegfa", "Hif1a", "Nf1", "Apc", "Ctnnb1", "Wnt1", "Tgfb1", "Smad4", "Cdkn2a",
                   paste0("Gene", 21:n_genes))  # Fill with generic names
  
  mock_results_symbols <- data.frame(
    baseMean = runif(n_genes, 10, 1000),
    log2FoldChange = rnorm(n_genes, 0, 2),
    lfcSE = runif(n_genes, 0.1, 0.5),
    stat = rnorm(n_genes, 0, 3),
    pvalue = runif(n_genes, 0, 1),
    padj = runif(n_genes, 0, 1),
    row.names = gene_symbols
  )
  
  # Make some genes significant
  n_sig <- 20
  sig_idx <- sample(1:n_genes, n_sig)
  mock_results_symbols$padj[sig_idx] <- runif(n_sig, 0.001, 0.049)
  mock_results_symbols$log2FoldChange[sig_idx] <- rnorm(n_sig, 0, 2)
  
  cat("✅ Mock results created with gene symbols:\n")
  cat("  - Total genes:", n_genes, "\n")
  cat("  - Significant genes:", sum(mock_results_symbols$padj < 0.05, na.rm = TRUE), "\n")
  cat("  - Sample gene names:", paste(head(rownames(mock_results_symbols), 5), collapse = ", "), "\n")
  
  # Test the gene list preparation
  cat("\n🔄 Testing prepare_gene_list_ora with gene symbols...\n")
  
  gene_list_symbols <- prepare_gene_list_ora(
    deseq2_results = mock_results_symbols,
    padj_cutoff = 0.05,
    fc_cutoff = 1.0,
    species = "mouse"
  )
  
  if (!is.null(gene_list_symbols)) {
    cat("✅ Symbol conversion test PASSED\n")
    cat("  - Converted genes:", length(gene_list_symbols), "\n")
    cat("  - Sample Entrez IDs:", paste(head(gene_list_symbols, 5), collapse = ", "), "\n")
  } else {
    cat("❌ Symbol conversion test FAILED\n")
  }
  
  # Test 2: Mock DESeq2 results with Ensembl IDs (original format)
  cat("\n📊 TEST 2: Ensembl IDs (original format)\n")
  cat("----------------------------------------\n")
  
  ensembl_ids <- paste0("ENSG", sprintf("%011d", 1:n_genes))
  
  mock_results_ensembl <- data.frame(
    baseMean = runif(n_genes, 10, 1000),
    log2FoldChange = rnorm(n_genes, 0, 2),
    lfcSE = runif(n_genes, 0.1, 0.5),
    stat = rnorm(n_genes, 0, 3),
    pvalue = runif(n_genes, 0, 1),
    padj = runif(n_genes, 0, 1),
    row.names = ensembl_ids
  )
  
  # Make some genes significant
  mock_results_ensembl$padj[sig_idx] <- runif(n_sig, 0.001, 0.049)
  mock_results_ensembl$log2FoldChange[sig_idx] <- rnorm(n_sig, 0, 2)
  
  cat("✅ Mock results created with Ensembl IDs:\n")
  cat("  - Total genes:", n_genes, "\n")
  cat("  - Significant genes:", sum(mock_results_ensembl$padj < 0.05, na.rm = TRUE), "\n")
  cat("  - Sample gene names:", paste(head(rownames(mock_results_ensembl), 3), collapse = ", "), "\n")
  
  # Test the gene list preparation
  cat("\n🔄 Testing prepare_gene_list_ora with Ensembl IDs...\n")
  
  gene_list_ensembl <- prepare_gene_list_ora(
    deseq2_results = mock_results_ensembl,
    padj_cutoff = 0.05,
    fc_cutoff = 1.0,
    species = "human"
  )
  
  if (!is.null(gene_list_ensembl)) {
    cat("✅ Ensembl conversion test PASSED\n")
    cat("  - Converted genes:", length(gene_list_ensembl), "\n")
    cat("  - Sample Entrez IDs:", paste(head(gene_list_ensembl, 5), collapse = ", "), "\n")
  } else {
    cat("❌ Ensembl conversion test FAILED\n")
  }
  
  return(list(
    symbols_test = !is.null(gene_list_symbols),
    ensembl_test = !is.null(gene_list_ensembl),
    symbols_count = if(!is.null(gene_list_symbols)) length(gene_list_symbols) else 0,
    ensembl_count = if(!is.null(gene_list_ensembl)) length(gene_list_ensembl) else 0
  ))
}

# Test pathway analysis end-to-end
test_full_pathway_analysis <- function() {
  cat("\n🚀 TESTING FULL PATHWAY ANALYSIS\n")
  cat("================================\n")
  
  # Create realistic test data with gene symbols
  cat("📊 Creating realistic test data with gene symbols...\n")
  
  # Use real mouse gene symbols
  real_mouse_genes <- c("Tp53", "Brca1", "Myc", "Kras", "Pik3ca", "Akt1", "Mtor", "Rb1", "Cdk4", "Ccnd1",
                       "Egfr", "Vegfa", "Hif1a", "Nf1", "Apc", "Ctnnb1", "Tgfb1", "Smad4", "Cdkn2a", "Pten",
                       "Braf", "Mapk1", "Jun", "Fos", "Stat3", "Stat1", "Ifng", "Il6", "Tnf", "Nfkb1")
  
  set.seed(42)
  n_genes <- length(real_mouse_genes)
  
  realistic_results <- data.frame(
    baseMean = runif(n_genes, 50, 2000),
    log2FoldChange = rnorm(n_genes, 0, 1.5),
    lfcSE = runif(n_genes, 0.1, 0.4),
    stat = rnorm(n_genes, 0, 2.5),
    pvalue = runif(n_genes, 0, 1),
    padj = runif(n_genes, 0, 1),
    row.names = real_mouse_genes
  )
  
  # Make most genes significant (for testing)
  n_sig <- 25
  sig_idx <- sample(1:n_genes, n_sig)
  realistic_results$padj[sig_idx] <- runif(n_sig, 0.001, 0.049)
  realistic_results$log2FoldChange[sig_idx] <- rnorm(n_sig, 0, 2)
  
  cat("✅ Realistic test data created:\n")
  cat("  - Real mouse genes:", n_genes, "\n")
  cat("  - Significant genes:", sum(realistic_results$padj < 0.05, na.rm = TRUE), "\n")
  
  # Test GO pathway analysis
  cat("\n🧬 Testing GO pathway analysis...\n")
  
  pathway_result <- tryCatch({
    run_pathway_analysis(
      deseq2_results = realistic_results,
      analysis_type = "GO",
      species = "mouse",
      ontology = "BP",
      padj_cutoff = 0.05,
      fc_cutoff = 1.0
    )
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  if (pathway_result$success) {
    cat("🎉 FULL PATHWAY ANALYSIS SUCCESS!\n")
    cat("  - Analysis type:", pathway_result$analysis_type, "\n")
    cat("  - Species:", pathway_result$species, "\n")
    cat("  - Enriched terms:", pathway_result$n_terms %||% nrow(pathway_result$data), "\n")
    cat("  - Execution time:", round(pathway_result$execution_time, 2), "seconds\n")
  } else {
    cat("❌ Full pathway analysis FAILED:", pathway_result$error, "\n")
  }
  
  return(pathway_result)
}

# Run all tests
cat("🚀 STARTING GENE CONVERSION FIX TESTS\n")
cat("=====================================\n\n")

# Test 1: Gene ID detection and conversion
conversion_tests <- test_gene_id_detection_and_conversion()

# Test 2: Full pathway analysis
pathway_test <- test_full_pathway_analysis()

# Summary
cat("\n📋 TEST SUMMARY\n")
cat("===============\n")

cat("1. Gene Symbol Detection & Conversion:")
if (conversion_tests$symbols_test) {
  cat("✅ PASS (", conversion_tests$symbols_count, "genes converted)\n")
} else {
  cat("❌ FAIL\n")
}

cat("2. Ensembl ID Detection & Conversion:")
if (conversion_tests$ensembl_test) {
  cat("✅ PASS (", conversion_tests$ensembl_count, "genes converted)\n")
} else {
  cat("❌ FAIL\n")
}

cat("3. Full Pathway Analysis:")
if (pathway_test$success) {
  cat("✅ PASS - Found", pathway_test$n_terms %||% 0, "enriched pathways\n")
} else {
  cat("❌ FAIL -", pathway_test$error %||% "Unknown error", "\n")
}

# Overall assessment
overall_success <- conversion_tests$symbols_test && pathway_test$success

cat("\n🎯 OVERALL RESULT:")
if (overall_success) {
  cat("✅ GENE CONVERSION FIX WORKING\n")
  cat("🎉 The pathway analysis should now work correctly with gene symbols!\n")
} else {
  cat("❌ SOME ISSUES REMAIN\n")
  cat("🔧 Check individual test results above\n")
}

cat("\n💡 NEXT STEPS:\n")
cat("===============\n")
cat("1. Restart the Shiny app\n")
cat("2. Run DESeq2 analysis (with pre-conversion gene symbols)\n")
cat("3. Try pathway analysis - should now detect gene symbols correctly\n")
cat("4. Gene conversion rate should be much higher (50-80%+)\n")
cat("5. Pathway analysis should find enriched terms successfully\n")