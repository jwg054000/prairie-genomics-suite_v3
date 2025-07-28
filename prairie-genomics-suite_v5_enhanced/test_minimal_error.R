# Minimal test to isolate error
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

cat('Running minimal GO test...\n')
result <- run_pathway_analysis(
  deseq2_results = test_data,
  analysis_type = 'GO',
  species = 'mouse'
)

# Check detailed error
if (is.list(result) && 'success' %in% names(result)) {
  if (!result$success) {
    cat('\nDETAILED ERROR INFO:\n')
    cat('Error message:', result$error, '\n')
    cat('Full message:', result$message, '\n')
    if ('suggestion' %in% names(result)) {
      cat('Suggestion:', result$suggestion, '\n')
    }
  } else {
    cat('\n✅ GO analysis succeeded!\n')
    if ('go_results' %in% names(result)) {
      cat('Found', nrow(result$go_results), 'GO pathways\n')
    }
  }
}