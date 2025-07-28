# Test pathway analysis fixes
source('pathway_analysis.R')

# Create test data
test_data <- data.frame(
  Gene = c('Gapdh', 'Actb', 'Ppia', 'B2m', 'Rplp0', 'Hprt', 'Tbp', 'Pgk1', 'Ldha', 'Sdha'),
  baseMean = c(1000, 800, 600, 500, 400, 300, 200, 150, 100, 80),
  log2FoldChange = c(-2.5, 2.3, -1.8, 1.5, -1.2, 2.8, -3.1, 1.9, -2.2, 2.0),
  lfcSE = rep(0.3, 10),
  stat = c(-8.3, 7.7, -6.0, 5.0, -4.0, 9.3, -10.3, 6.3, -7.3, 6.7),
  pvalue = c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-11, 1e-12, 1e-8, 1e-9, 1e-8),
  padj = c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-10, 1e-11, 1e-7, 1e-8, 1e-7)
)

cat('========================================\n')
cat('TESTING PATHWAY ANALYSIS FIXES\n')
cat('========================================\n\n')

# Test GO analysis
cat('1. Testing GO analysis...\n')
go_result <- run_pathway_analysis(
  deseq2_results = test_data,
  analysis_type = 'GO',
  species = 'mouse',
  padj_cutoff = 0.05,
  fc_cutoff = 1.0
)

if (is.list(go_result) && 'success' %in% names(go_result)) {
  if (go_result$success) {
    cat('   ✅ GO analysis SUCCEEDED!\n')
    if ('go_results' %in% names(go_result) && is.data.frame(go_result$go_results)) {
      cat('   ✅ Found', nrow(go_result$go_results), 'GO pathways\n')
    }
  } else {
    cat('   ❌ GO analysis FAILED\n')
    cat('   Error:', go_result$error, '\n')
  }
}

# Test KEGG analysis
cat('\n2. Testing KEGG analysis...\n')
kegg_result <- run_pathway_analysis(
  deseq2_results = test_data,
  analysis_type = 'KEGG',
  species = 'mouse',
  padj_cutoff = 0.05,
  fc_cutoff = 1.0
)

if (is.list(kegg_result) && 'success' %in% names(kegg_result)) {
  if (kegg_result$success) {
    cat('   ✅ KEGG analysis SUCCEEDED!\n')
    if ('kegg_results' %in% names(kegg_result) && is.data.frame(kegg_result$kegg_results)) {
      cat('   ✅ Found', nrow(kegg_result$kegg_results), 'KEGG pathways\n')
    }
  } else {
    cat('   ❌ KEGG analysis FAILED\n')
    cat('   Error:', kegg_result$error, '\n')
  }
}

cat('\n========================================\n')
cat('TEST COMPLETE\n')
cat('========================================\n')