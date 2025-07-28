# Test with traceback to find exact error location
options(error = function() {
  cat("\n🔥 ERROR TRACEBACK:\n")
  print(traceback())
})

source('pathway_analysis.R')

# Create simple test data
test_data <- data.frame(
  Gene = paste0('Gene', 1:5),
  baseMean = rep(1000, 5),
  log2FoldChange = c(-2, 2, -3, 3, -1.5),
  lfcSE = rep(0.3, 5),
  pvalue = rep(1e-5, 5),
  padj = rep(1e-4, 5)
)

cat('Running GO test with traceback enabled...\n')
result <- run_pathway_analysis(
  deseq2_results = test_data,
  analysis_type = 'GO',
  species = 'mouse'
)

if (is.list(result) && 'success' %in% names(result) && !result$success) {
  cat('\nError occurred:', result$error, '\n')
}