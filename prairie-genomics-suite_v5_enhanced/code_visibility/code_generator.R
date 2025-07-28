# Code Generator Module for Prairie Genomics Suite
# Specialized functions for generating analysis-specific R code
#
# Author: Prairie Genomics Team
# Date: January 24, 2025
# Purpose: Generate clean, documented R code for each analysis type

# Source the code logger
source("code_visibility/code_logger.R")

# Generate visualization code
generate_visualization_code <- function(plot_type, data_info, parameters = list()) {
  code_lines <- c(
    paste0("# ", stringr::str_to_title(plot_type), " Plot Generation"),
    paste0("# ", paste(rep("=", nchar(plot_type) + 16), collapse = "")),
    "",
    "library(ggplot2)"
  )
  
  if (plot_type == "volcano") {
    code_lines <- c(code_lines,
      "",
      "# Create volcano plot",
      "volcano_data <- results_df",
      "volcano_data$significance <- ifelse(",
      paste0("  !is.na(volcano_data$padj) & volcano_data$padj < ", parameters$padj_cutoff %||% 0.05, " &"),
      paste0("  abs(volcano_data$log2FoldChange) > ", parameters$fc_cutoff %||% 1),
      "  'Significant', 'Not Significant'",
      ")",
      "",
      "# Generate volcano plot",
      "volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +",
      "  geom_point(aes(color = significance), alpha = 0.6) +",
      "  scale_color_manual(values = c('Significant' = 'red', 'Not Significant' = 'gray')) +",
      "  geom_vline(xintercept = c(-1, 1), linetype = 'dashed', alpha = 0.5) +",
      "  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', alpha = 0.5) +",
      "  labs(",
      "    title = 'Volcano Plot',",
      "    x = 'Log2 Fold Change',",
      "    y = '-Log10 Adjusted P-value'",
      "  ) +",
      "  theme_minimal()",
      "",
      "print(volcano_plot)"
    )
  } else if (plot_type == "ma") {
    code_lines <- c(code_lines,
      "",
      "# Create MA plot",
      "ma_data <- results_df[!is.na(results_df$baseMean) & !is.na(results_df$log2FoldChange), ]",
      "ma_data$significance <- ifelse(",
      paste0("  !is.na(ma_data$padj) & ma_data$padj < ", parameters$padj_cutoff %||% 0.05),
      "  'Significant', 'Not Significant'",
      ")",
      "",
      "# Generate MA plot",
      "ma_plot <- ggplot(ma_data, aes(x = log10(baseMean), y = log2FoldChange)) +",
      "  geom_point(aes(color = significance), alpha = 0.6) +",
      "  scale_color_manual(values = c('Significant' = 'red', 'Not Significant' = 'gray')) +",
      "  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +",
      "  labs(",
      "    title = 'MA Plot',",
      "    x = 'Log10 Mean Expression',",
      "    y = 'Log2 Fold Change'",
      "  ) +",
      "  theme_minimal()",
      "",
      "print(ma_plot)"
    )
  } else if (plot_type == "heatmap") {
    code_lines <- c(code_lines,
      "library(pheatmap)",
      "",
      "# Prepare data for heatmap",
      paste0("top_genes <- head(significant_genes[order(significant_genes$padj), ], ", parameters$top_genes %||% 50, ")"),
      "",
      "# Get normalized counts",
      "norm_counts <- counts(dds, normalized = TRUE)",
      "heatmap_data <- norm_counts[rownames(top_genes), ]",
      "",
      "# Log transform for better visualization",
      "heatmap_data <- log2(heatmap_data + 1)",
      "",
      "# Generate heatmap",
      "pheatmap(",
      "  heatmap_data,",
      "  cluster_rows = TRUE,",
      "  cluster_cols = TRUE,",
      "  show_rownames = TRUE,",
      "  show_colnames = TRUE,",
      "  scale = 'row',",
      "  main = 'Top Significant Genes Heatmap'",
      ")"
    )
  } else if (plot_type == "pca") {
    code_lines <- c(code_lines,
      "",
      "# PCA analysis",
      "# Use variance stabilizing transformation",
      "vsd <- vst(dds, blind = FALSE)",
      "",
      "# Calculate PCA",
      "pca_data <- plotPCA(vsd, intgroup = 'group', returnData = TRUE)",
      "percent_var <- round(100 * attr(pca_data, 'percentVar'))",
      "",
      "# Generate PCA plot",
      "pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = group)) +",
      "  geom_point(size = 3) +",
      "  labs(",
      "    title = 'Principal Component Analysis',",
      "    x = paste0('PC1: ', percent_var[1], '% variance'),",
      "    y = paste0('PC2: ', percent_var[2], '% variance')",
      "  ) +",
      "  theme_minimal()",
      "",
      "print(pca_plot)"
    )
  }
  
  return(paste(code_lines, collapse = "\n"))
}

# Generate gene conversion code
generate_gene_conversion_code <- function(species, conversion_type = "ensembl_to_symbol") {
  code_lines <- c(
    "# Gene ID Conversion",
    "# ==================",
    ""
  )
  
  if (species == "human") {
    code_lines <- c(code_lines,
      "library(org.Hs.eg.db)",
      ""
    )
    
    if (conversion_type == "ensembl_to_symbol") {
      code_lines <- c(code_lines,
        "# Convert Ensembl IDs to Gene Symbols",
        "gene_symbols <- mapIds(",
        "  org.Hs.eg.db,",
        "  keys = rownames(results_df),",
        "  column = 'SYMBOL',",
        "  keytype = 'ENSEMBL',",
        "  multiVals = 'first'",
        ")",
        "",
        "# Add gene symbols to results",
        "results_df$gene_symbol <- gene_symbols[rownames(results_df)]",
        ""
      )
    } else if (conversion_type == "ensembl_to_entrez") {
      code_lines <- c(code_lines,
        "# Convert Ensembl IDs to Entrez IDs",
        "entrez_ids <- mapIds(",
        "  org.Hs.eg.db,",
        "  keys = rownames(results_df),",
        "  column = 'ENTREZID',",
        "  keytype = 'ENSEMBL',",
        "  multiVals = 'first'",
        ")",
        "",
        "# Add Entrez IDs to results",
        "results_df$entrez_id <- entrez_ids[rownames(results_df)]",
        ""
      )
    }
  } else if (species == "mouse") {
    code_lines <- c(code_lines,
      "library(org.Mm.eg.db)",
      ""
    )
    
    if (conversion_type == "ensembl_to_symbol") {
      code_lines <- c(code_lines,
        "# Convert Ensembl IDs to Gene Symbols",
        "gene_symbols <- mapIds(",
        "  org.Mm.eg.db,",
        "  keys = rownames(results_df),",
        "  column = 'SYMBOL',",
        "  keytype = 'ENSEMBL',",
        "  multiVals = 'first'",
        ")",
        "",
        "# Add gene symbols to results",
        "results_df$gene_symbol <- gene_symbols[rownames(results_df)]",
        ""
      )
    }
  }
  
  return(paste(code_lines, collapse = "\n"))
}

# Generate export code
generate_export_code <- function(export_items = c("results", "significant_genes", "plots")) {
  code_lines <- c(
    "# Export Results",
    "# ==============",
    ""
  )
  
  if ("results" %in% export_items) {
    code_lines <- c(code_lines,
      "# Export complete results",
      "write.csv(results_df, 'deseq2_results.csv', row.names = TRUE)",
      ""
    )
  }
  
  if ("significant_genes" %in% export_items) {
    code_lines <- c(code_lines,
      "# Export significant genes only",
      "write.csv(significant_genes, 'significant_genes.csv', row.names = TRUE)",
      ""
    )
  }
  
  if ("plots" %in% export_items) {
    code_lines <- c(code_lines,
      "# Save plots",
      "ggsave('volcano_plot.png', volcano_plot, width = 8, height = 6, dpi = 300)",
      "ggsave('ma_plot.png', ma_plot, width = 8, height = 6, dpi = 300)",
      ""
    )
  }
  
  return(paste(code_lines, collapse = "\n"))
}

# Generate complete analysis template
generate_analysis_template <- function(analysis_type = "basic_deseq2") {
  if (analysis_type == "basic_deseq2") {
    template <- c(
      "# Prairie Genomics Suite - Basic DESeq2 Analysis Template",
      "# ======================================================",
      "",
      "# This template provides a complete workflow for differential",
      "# expression analysis using DESeq2",
      "",
      "# Load required libraries",
      "library(DESeq2)",
      "library(ggplot2)",
      "library(org.Hs.eg.db)  # Use org.Mm.eg.db for mouse",
      "",
      "# 1. Load data",
      "# ============",
      "# Replace 'your_count_file.csv' with your actual file",
      "count_matrix <- read.csv('your_count_file.csv', row.names = 1, check.names = FALSE)",
      "",
      "# Create sample information",
      "# Modify this to match your experimental design",
      "sample_info <- data.frame(",
      "  sample_id = colnames(count_matrix),",
      "  group = rep(c('control', 'treatment'), each = 3),  # Adjust as needed",
      "  stringsAsFactors = FALSE",
      ")",
      "rownames(sample_info) <- sample_info$sample_id",
      "",
      "# 2. Quality control and preprocessing",
      "# ===================================",
      "# Remove genes with zero counts",
      "count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]",
      "",
      "# Ensure integer counts",
      "count_matrix <- round(count_matrix)",
      "",
      "# Check data",
      "cat('Count matrix dimensions:', dim(count_matrix), '\\n')",
      "head(count_matrix[,1:min(5, ncol(count_matrix))])",
      "",
      "# 3. DESeq2 analysis",
      "# ==================",
      "# Create DESeq2 dataset",
      "dds <- DESeqDataSetFromMatrix(",
      "  countData = count_matrix,",
      "  colData = sample_info,",
      "  design = ~ group",
      ")",
      "",
      "# Run analysis",
      "dds <- DESeq(dds)",
      "",
      "# Extract results",
      "results <- results(dds, contrast = c('group', 'treatment', 'control'))",
      "results_df <- as.data.frame(results)",
      "results_df$gene_id <- rownames(results_df)",
      "",
      "# 4. Filter significant genes",
      "# ===========================",
      "significant_genes <- results_df[",
      "  !is.na(results_df$padj) & ",
      "  results_df$padj < 0.05 & ",
      "  abs(results_df$log2FoldChange) > 1,",
      "]",
      "",
      "cat('Found', nrow(significant_genes), 'significant genes\\n')",
      "",
      "# 5. Add gene symbols",
      "# ===================",
      "gene_symbols <- mapIds(",
      "  org.Hs.eg.db,",
      "  keys = results_df$gene_id,",
      "  column = 'SYMBOL',",
      "  keytype = 'ENSEMBL',",
      "  multiVals = 'first'",
      ")",
      "results_df$gene_symbol <- gene_symbols",
      "",
      "# 6. Visualization",
      "# ================",
      "# Volcano plot",
      "volcano_data <- results_df",
      "volcano_data$significance <- ifelse(",
      "  !is.na(volcano_data$padj) & volcano_data$padj < 0.05 &",
      "  abs(volcano_data$log2FoldChange) > 1,",
      "  'Significant', 'Not Significant'",
      ")",
      "",
      "volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +",
      "  geom_point(aes(color = significance), alpha = 0.6) +",
      "  scale_color_manual(values = c('Significant' = 'red', 'Not Significant' = 'gray')) +",
      "  geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +",
      "  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +",
      "  labs(title = 'Volcano Plot', x = 'Log2 Fold Change', y = '-Log10 Adjusted P-value') +",
      "  theme_minimal()",
      "",
      "print(volcano_plot)",
      "",
      "# 7. Export results",
      "# =================",
      "write.csv(results_df, 'deseq2_results.csv', row.names = TRUE)",
      "write.csv(significant_genes, 'significant_genes.csv', row.names = TRUE)",
      "ggsave('volcano_plot.png', volcano_plot, width = 8, height = 6, dpi = 300)",
      "",
      "# 8. Session information",
      "# ======================",
      "sessionInfo()",
      ""
    )
  } else if (analysis_type == "pathway_analysis") {
    template <- c(
      "# Prairie Genomics Suite - Pathway Analysis Template",
      "# ==================================================",
      "",
      "# This template provides pathway analysis following DESeq2",
      "",
      "# Load required libraries",
      "library(clusterProfiler)",
      "library(fgsea)",
      "library(msigdbr)",
      "library(org.Hs.eg.db)  # Use org.Mm.eg.db for mouse",
      "library(ggplot2)",
      "",
      "# Assuming you have results_df from DESeq2 analysis",
      "# If not, run the basic DESeq2 template first",
      "",
      "# 1. Gene Ontology (GO) Analysis",
      "# ==============================",
      "# Prepare significant genes",
      "sig_genes <- results_df[",
      "  !is.na(results_df$padj) & results_df$padj < 0.05 &",
      "  abs(results_df$log2FoldChange) > 1,",
      "]",
      "",
      "# Convert to Entrez IDs",
      "entrez_ids <- mapIds(org.Hs.eg.db, keys = sig_genes$gene_id,",
      "                     column = 'ENTREZID', keytype = 'ENSEMBL')",
      "entrez_ids <- entrez_ids[!is.na(entrez_ids)]",
      "",
      "# Run GO enrichment",
      "go_results <- enrichGO(",
      "  gene = entrez_ids,",
      "  OrgDb = org.Hs.eg.db,",
      "  ont = 'BP',",
      "  pAdjustMethod = 'BH',",
      "  pvalueCutoff = 0.05,",
      "  readable = TRUE",
      ")",
      "",
      "# 2. GSEA Analysis",
      "# ================",
      "# Prepare ranked gene list",
      "clean_results <- results_df[!is.na(results_df$log2FoldChange) & !is.na(results_df$pvalue), ]",
      "gene_ranks <- sign(clean_results$log2FoldChange) * (-log10(clean_results$pvalue))",
      "names(gene_ranks) <- clean_results$gene_id",
      "gene_ranks <- sort(gene_ranks, decreasing = TRUE)",
      "",
      "# Get Hallmark gene sets",
      "gene_sets <- msigdbr(species = 'Homo sapiens', collection = 'H')",
      "pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)",
      "",
      "# Run GSEA",
      "gsea_results <- fgsea(",
      "  pathways = pathways,",
      "  stats = gene_ranks,",
      "  minSize = 15,",
      "  maxSize = 500",
      ")",
      "",
      "# 3. Visualization",
      "# ================",
      "# GO dotplot",
      "if (!is.null(go_results) && nrow(go_results@result) > 0) {",
      "  dotplot(go_results, showCategory = 20)",
      "}",
      "",
      "# GSEA enrichment plot",
      "if (!is.null(gsea_results) && nrow(gsea_results) > 0) {",
      "  plotEnrichment(pathways[[1]], gene_ranks)",
      "}",
      "",
      "# Export results",
      "if (!is.null(go_results)) write.csv(go_results@result, 'go_results.csv')",
      "if (!is.null(gsea_results)) write.csv(gsea_results, 'gsea_results.csv')",
      ""
    )
  }
  
  return(paste(template, collapse = "\n"))
}

# Helper function for NULL coalescing
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

# Convert R script to R Markdown format
convert_to_rmarkdown <- function(r_code, session_id = NULL) {
  # Split code into lines
  code_lines <- strsplit(r_code, "\n")[[1]]
  
  # Initialize Rmd content
  rmd_lines <- c(
    "---",
    "title: \"Prairie Genomics Suite Analysis Report\"",
    if (!is.null(session_id)) paste0("subtitle: \"Session: ", session_id, "\"") else NULL,
    paste0("date: \"", Sys.Date(), "\""),
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    code_folding: show",
    "    theme: flatly",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)",
    "```",
    "",
    "## Analysis Overview",
    "",
    "This report contains the complete R code and results from your genomics analysis",
    "performed using the Prairie Genomics Suite web application.",
    "",
    "## Analysis Code",
    ""
  )
  
  # Process code lines to create appropriate chunks
  in_chunk <- FALSE
  current_chunk <- character()
  chunk_title <- ""
  
  for (line in code_lines) {
    # Check if this is a section header
    if (grepl("^# [A-Z].*[=]{3,}$|^# [A-Z].*Analysis$|^# [A-Z].*Setup$", line)) {
      # Close previous chunk if open
      if (in_chunk) {
        rmd_lines <- c(rmd_lines, "```", "")
        in_chunk <- FALSE
      }
      
      # Extract title
      chunk_title <- gsub("^# (.*?)( [=]+)?$", "\\1", line)
      rmd_lines <- c(rmd_lines, paste0("### ", chunk_title), "")
      
      # Start new chunk
      rmd_lines <- c(rmd_lines, "```{r}", line)
      in_chunk <- TRUE
    } else if (line != "" && !grepl("^#", line)) {
      # Regular code line
      if (!in_chunk) {
        rmd_lines <- c(rmd_lines, "```{r}")
        in_chunk <- TRUE
      }
      rmd_lines <- c(rmd_lines, line)
    } else if (grepl("^# ", line) && !grepl("[=]+", line)) {
      # Comment line
      if (!in_chunk) {
        rmd_lines <- c(rmd_lines, "```{r}")
        in_chunk <- TRUE
      }
      rmd_lines <- c(rmd_lines, line)
    } else if (line == "") {
      # Empty line
      if (in_chunk) {
        rmd_lines <- c(rmd_lines, line)
      }
    }
  }
  
  # Close final chunk if open
  if (in_chunk) {
    rmd_lines <- c(rmd_lines, "```")
  }
  
  # Add session info section
  rmd_lines <- c(rmd_lines,
    "",
    "## Session Information",
    "",
    "```{r sessioninfo}",
    "sessionInfo()",
    "```"
  )
  
  return(paste(rmd_lines, collapse = "\n"))
}

cat("✅ Code generator module loaded\n")
cat("📋 Functions available: generate_visualization_code(), generate_gene_conversion_code(), generate_analysis_template()\n")