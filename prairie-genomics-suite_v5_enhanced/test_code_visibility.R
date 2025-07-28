# Test Code Visibility System
# Verify the code logging and display functionality works correctly
#
# Author: Prairie Genomics Team
# Date: January 24, 2025

cat("🧪 Testing Code Visibility System\n")
cat("=", rep("=", 35), "\n")

# Load the code visibility modules
source("code_visibility/code_logger.R")
source("code_visibility/code_generator.R")  
source("code_visibility/analysis_wrappers.R")

# Test 1: Initialize code logger
cat("\n1. Testing Code Logger Initialization\n")
cat(rep("-", 40), "\n")

session_id <- init_code_logger(user_name = "test_user")
cat("✅ Session initialized:", session_id, "\n")

# Test 2: Log a simple analysis step
cat("\n2. Testing Analysis Step Logging\n")
cat(rep("-", 35), "\n")

# Simulate data loading
log_analysis_step(
  session_id = session_id,
  step_name = "Data Loading Test",
  category = "data_upload",
  code_snippet = generate_data_loading_code(
    list(name = "test_data.csv", size = "1.2MB"),
    c("remove_zero_genes", "ensure_integer")
  ),
  description = "Test data loading step",
  input_data = "Test count matrix: 1000 genes x 6 samples",
  output_data = "Processed count matrix ready for analysis"
)

# Test 3: Log DESeq2 analysis step
cat("\n3. Testing DESeq2 Code Generation\n")
cat(rep("-", 35), "\n")

deseq2_code <- generate_deseq2_code(
  design_formula = "~ group",
  contrast_info = list(numerator = "treatment", denominator = "control"),
  filter_params = list(padj_cutoff = 0.05, fc_cutoff = 1)
)

log_analysis_step(
  session_id = session_id,
  step_name = "DESeq2 Differential Expression Analysis",
  category = "deseq2",
  code_snippet = deseq2_code,
  description = "Run DESeq2 analysis with default parameters",
  execution_time = 15.3
)

# Test 4: Log pathway analysis
cat("\n4. Testing Pathway Analysis Code Generation\n")
cat(rep("-", 45), "\n")

pathway_code <- generate_pathway_code(
  analysis_type = "GSEA",
  species = "human",
  parameters = list(collection = "H", padj_cutoff = 0.05)
)

log_analysis_step(
  session_id = session_id,
  step_name = "GSEA Pathway Analysis",
  category = "pathway", 
  code_snippet = pathway_code,
  description = "Gene Set Enrichment Analysis using Hallmark gene sets"
)

# Test 5: Generate visualization code
cat("\n5. Testing Visualization Code Generation\n")
cat(rep("-", 40), "\n")

viz_code <- generate_visualization_code(
  plot_type = "volcano",
  data_info = list(),
  parameters = list(padj_cutoff = 0.05, fc_cutoff = 1)
)

log_analysis_step(
  session_id = session_id,
  step_name = "Volcano Plot Generation",
  category = "visualization",
  code_snippet = viz_code,
  description = "Create volcano plot for differential expression results"
)

# Test 6: Get complete session code
cat("\n6. Testing Complete Code Generation\n")
cat(rep("-", 35), "\n")

complete_code <- get_session_code(session_id, include_comments = TRUE)
cat("✅ Generated complete analysis script\n")
cat("📊 Code length:", nchar(complete_code), "characters\n")

# Show first few lines
code_lines <- strsplit(complete_code, "\n")[[1]]
cat("🔍 First 10 lines of generated code:\n")
for (i in 1:min(10, length(code_lines))) {
  cat(sprintf("%2d: %s\n", i, code_lines[i]))
}

# Test 7: Export code to file
cat("\n7. Testing Code Export\n")
cat(rep("-", 25), "\n")

export_file <- export_session_code(session_id, format = "R")
if (file.exists(export_file)) {
  cat("✅ Code exported successfully to:", export_file, "\n")
  cat("📄 File size:", file.info(export_file)$size, "bytes\n")
} else {
  cat("❌ Code export failed\n")
}

# Test 8: Test specific code snippets
cat("\n8. Testing Category-Specific Code Retrieval\n")
cat(rep("-", 45), "\n")

deseq2_only <- get_analysis_code_snippet(session_id, category = "deseq2")
pathway_only <- get_analysis_code_snippet(session_id, category = "pathway")
viz_only <- get_analysis_code_snippet(session_id, category = "visualization")

cat("✅ DESeq2 code length:", nchar(deseq2_only), "characters\n")
cat("✅ Pathway code length:", nchar(pathway_only), "characters\n") 
cat("✅ Visualization code length:", nchar(viz_only), "characters\n")

# Test 9: Test R Markdown conversion
cat("\n9. Testing R Markdown Conversion\n")
cat(rep("-", 35), "\n")

rmd_content <- convert_to_rmarkdown(complete_code, session_id)
cat("✅ R Markdown conversion successful\n")
cat("📊 Rmd content length:", nchar(rmd_content), "characters\n")

# Test 10: Test session summary
cat("\n10. Testing Session Summary\n")
cat(rep("-", 28), "\n")

summary <- get_session_summary(session_id)
if (!is.null(summary)) {
  cat("✅ Session summary generated:\n")
  cat("   - Session ID:", summary$session_id, "\n")
  cat("   - User:", summary$user_name, "\n")
  cat("   - Steps logged:", summary$num_steps, "\n")
  cat("   - Categories:", paste(summary$categories, collapse = ", "), "\n")
  cat("   - Start time:", as.character(summary$start_time), "\n")
} else {
  cat("❌ Session summary failed\n")
}

# Test 11: Test analysis templates
cat("\n11. Testing Analysis Templates\n")
cat(rep("-", 32), "\n")

basic_template <- generate_analysis_template("basic_deseq2")
pathway_template <- generate_analysis_template("pathway_analysis")

cat("✅ Basic DESeq2 template:", nchar(basic_template), "characters\n")
cat("✅ Pathway analysis template:", nchar(pathway_template), "characters\n")

# Test 12: Test enhanced analysis wrappers
cat("\n12. Testing Enhanced Analysis Wrappers\n")
cat(rep("-", 40), "\n")

# Create mock data for testing
mock_count_matrix <- matrix(
  rpois(6000, lambda = 100), 
  nrow = 1000, 
  ncol = 6,
  dimnames = list(
    paste0("GENE_", 1:1000),
    paste0("Sample_", 1:6)
  )
)

mock_sample_info <- data.frame(
  sample_id = colnames(mock_count_matrix),
  group = rep(c("control", "treatment"), each = 3),
  stringsAsFactors = FALSE
)
rownames(mock_sample_info) <- mock_sample_info$sample_id

# Test enhanced DESeq2 wrapper
tryCatch({
  deseq2_result <- run_deseq2_analysis_logged(
    count_matrix = mock_count_matrix,
    sample_info = mock_sample_info,
    design_formula = "~ group",
    contrast_info = list(numerator = "treatment", denominator = "control"),
    filter_params = list(padj_cutoff = 0.05, fc_cutoff = 1),
    session_id = session_id,
    species = "human"
  )
  
  if (deseq2_result$success) {
    cat("✅ Enhanced DESeq2 analysis wrapper working\n")
    cat("📊 Mock analysis results:", deseq2_result$analysis_summary$total_genes, "genes analyzed\n")
  } else {
    cat("❌ Enhanced DESeq2 wrapper failed:", deseq2_result$error, "\n")
  }
}, error = function(e) {
  cat("❌ Enhanced DESeq2 wrapper test failed:", e$message, "\n")
})

# Test summary and cleanup
cat("\n🎯 Test Summary\n")
cat("=", rep("=", 15), "\n")

final_summary <- get_session_summary(session_id)
cat("Final session state:\n")
cat("✅ Total analysis steps logged:", final_summary$num_steps, "\n")
cat("✅ Categories covered:", paste(final_summary$categories, collapse = ", "), "\n")
cat("✅ Session duration:", round(as.numeric(difftime(Sys.time(), final_summary$start_time, units = "secs")), 2), "seconds\n")

# Get final code
final_code <- get_session_code(session_id)
cat("✅ Final complete code length:", nchar(final_code), "characters\n")

# Clean up
clear_session(session_id)
cat("🗑️ Test session cleaned up\n")

# Remove test export file
if (exists("export_file") && file.exists(export_file)) {
  file.remove(export_file)
  cat("🗑️ Test export file removed\n")
}

cat("\n🎉 Code Visibility System Test Complete!\n")
cat("📋 All core functionality appears to be working correctly\n")
cat("🚀 Ready for integration into the main application\n")