#!/usr/bin/env Rscript

# Prairie Genomics Suite - Detailed Analysis with Gene Symbol Conversion
# MC9 vs MLM Analysis with Complete Transparency

cat("🧬 Prairie Genomics Suite - Detailed MC9 vs MLM Analysis\n")
cat("🔍 Complete Code Transparency + Gene Symbol Conversion\n")
cat("=" , rep("=", 70), "\n", sep="")

# Load required libraries
suppressPackageStartupMessages({
  library(readr)
  library(DESeq2)
  
  # For gene annotation
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", quiet = TRUE)
    }
    BiocManager::install("biomaRt", quiet = TRUE)
  }
  library(biomaRt)
})

cat("📦 Libraries loaded: readr, DESeq2, biomaRt\n\n")

# =============================================================================
# STEP 1: DATA LOADING AND PREPARATION
# =============================================================================

cat("📊 STEP 1: DATA LOADING AND PREPARATION\n")
cat("-" , rep("-", 50), "\n", sep="")

# Load expression data
expr_file <- "/Users/joshuagarton/Desktop/MC9.raw.counts.test.csv"
meta_file <- "/Users/joshuagarton/Desktop/MC9_sample_metadata.csv"

cat("📈 Loading expression matrix...\n")
expr_data <- read_csv(expr_file, show_col_types = FALSE)

# Convert to matrix format  
gene_names <- expr_data[[1]]
count_matrix <- as.matrix(expr_data[, -1])
rownames(count_matrix) <- gene_names

cat("📋 Loading sample metadata...\n")
meta_data <- read_csv(meta_file, show_col_types = FALSE)

cat("✅ Data loaded successfully\n")
cat("   - Expression matrix:", nrow(count_matrix), "genes ×", ncol(count_matrix), "samples\n")
cat("   - Metadata:", nrow(meta_data), "samples\n\n")

# Filter for MC9 vs MLM comparison
mc9_samples <- meta_data$Sample_ID[meta_data$Condition == "MC9"]
mlm_samples <- meta_data$Sample_ID[meta_data$Condition == "MLM"]

comparison_samples <- c(mc9_samples, mlm_samples)
comparison_matrix <- count_matrix[, comparison_samples]
comparison_metadata <- meta_data[meta_data$Sample_ID %in% comparison_samples, ]
comparison_metadata <- comparison_metadata[match(colnames(comparison_matrix), comparison_metadata$Sample_ID), ]

cat("🔬 Prepared MC9 vs MLM comparison:\n")
cat("   - MC9 samples:", length(mc9_samples), "(", paste(mc9_samples, collapse = ", "), ")\n")
cat("   - MLM samples:", length(mlm_samples), "(", paste(mlm_samples, collapse = ", "), ")\n")
cat("   - Total genes:", nrow(comparison_matrix), "\n\n")

# =============================================================================
# STEP 2: DESEQ2 ANALYSIS WITH TRANSPARENT CODE
# =============================================================================

cat("🧮 STEP 2: DESEQ2 ANALYSIS - SHOWING ALL CODE\n")
cat("-" , rep("-", 50), "\n", sep="")

cat("📝 CODE: Creating DESeq2 dataset\n")
cat("```R\n")
cat("# Ensure condition is a factor with MC9 as reference\n")
cat("comparison_metadata$Condition <- factor(comparison_metadata$Condition, levels = c('MC9', 'MLM'))\n")
cat("\n")
cat("# Create DESeq2 dataset\n")
cat("dds <- DESeqDataSetFromMatrix(\n")
cat("  countData = comparison_matrix,\n")
cat("  colData = comparison_metadata,\n")
cat("  design = ~ Condition\n")
cat(")\n")
cat("```\n\n")

# Execute the code
comparison_metadata$Condition <- factor(comparison_metadata$Condition, levels = c("MC9", "MLM"))

dds <- DESeqDataSetFromMatrix(
  countData = comparison_matrix,
  colData = comparison_metadata,
  design = ~ Condition
)

cat("📝 CODE: Filtering low-count genes\n")
cat("```R\n")
cat("# Filter genes with total counts < 10 across all samples\n")
cat("keep <- rowSums(counts(dds)) >= 10\n")
cat("dds <- dds[keep, ]\n")
cat("```\n\n")

# Execute filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("✅ Filtered from", nrow(comparison_matrix), "to", nrow(dds), "genes\n\n")

cat("📝 CODE: Running DESeq2 differential expression\n")
cat("```R\n")
cat("# Run DESeq2 analysis (normalization + dispersion + testing)\n")
cat("dds <- DESeq(dds, quiet = TRUE)\n")
cat("\n")
cat("# Extract results: MLM vs MC9 (MLM is numerator)\n")
cat("res <- results(dds, contrast = c('Condition', 'MLM', 'MC9'))\n")
cat("```\n\n")

# Execute DESeq2
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("Condition", "MLM", "MC9"))

cat("✅ DESeq2 analysis completed\n")
cat("   - Genes tested:", sum(!is.na(res$padj)), "\n")
cat("   - Results structure: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj\n\n")

# =============================================================================
# STEP 3: SIGNIFICANCE DETERMINATION - EXACT CODE
# =============================================================================

cat("🎯 STEP 3: SIGNIFICANCE DETERMINATION - EXACT CODE\n")
cat("-" , rep("-", 50), "\n", sep="")

cat("📝 CODE: Defining significance criteria\n")
cat("```R\n")
cat("# Expert-validated parameters (from Joshua's approval)\n")
cat("padj_threshold <- 0.05    # Adjusted p-value (FDR)\n")
cat("fc_threshold <- 1.5       # Fold change threshold\n")
cat("\n")
cat("# Convert fold change to log2 scale for filtering\n")
cat("log2fc_threshold <- log2(fc_threshold)  # log2(1.5) = 0.585\n")
cat("\n")
cat("# Apply significance filters\n")
cat("significant_genes <- subset(res, \n")
cat("                           padj < padj_threshold & \n")
cat("                           abs(log2FoldChange) >= log2fc_threshold)\n")
cat("\n")
cat("# Sort by adjusted p-value (most significant first)\n")
cat("significant_genes <- significant_genes[order(significant_genes$padj), ]\n")
cat("```\n\n")

# Execute significance filtering
padj_threshold <- 0.05
fc_threshold <- 1.5
log2fc_threshold <- log2(fc_threshold)

cat("🔍 Applied thresholds:\n")
cat("   - Adjusted p-value < ", padj_threshold, "\n")
cat("   - Absolute fold change ≥ ", fc_threshold, "x (log2FC ≥ ", round(log2fc_threshold, 3), ")\n\n")

significant_genes <- subset(res, 
                           padj < padj_threshold & 
                           abs(log2FoldChange) >= log2fc_threshold)

significant_genes <- significant_genes[order(significant_genes$padj), ]

cat("✅ Significance filtering results:\n")
cat("   - Total genes tested:", sum(!is.na(res$padj)), "\n")
cat("   - Significant genes:", nrow(significant_genes), "\n")
cat("   - Percentage significant:", round(nrow(significant_genes) / sum(!is.na(res$padj)) * 100, 1), "%\n\n")

# =============================================================================
# STEP 4: GENE SYMBOL CONVERSION
# =============================================================================

cat("🧬 STEP 4: CONVERTING ENSEMBL IDs TO GENE SYMBOLS\n")
cat("-" , rep("-", 50), "\n", sep="")

cat("📝 CODE: Setting up biomaRt for mouse gene annotation\n")
cat("```R\n")
cat("# Connect to Ensembl mouse database\n")
cat("mart <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')\n")
cat("\n")
cat("# Get gene symbols for significant genes\n")
cat("ensembl_ids <- rownames(significant_genes)\n")
cat("gene_info <- getBM(\n")
cat("  attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),\n")
cat("  filters = 'ensembl_gene_id',\n")
cat("  values = ensembl_ids,\n")
cat("  mart = mart\n")
cat(")\n")
cat("```\n\n")

# Execute gene annotation
tryCatch({
  cat("🔍 Connecting to Ensembl database...\n")
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  cat("🧬 Retrieving gene symbols for", nrow(significant_genes), "significant genes...\n")
  ensembl_ids <- rownames(significant_genes)
  
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id", 
    values = ensembl_ids,
    mart = mart
  )
  
  cat("✅ Retrieved gene information for", nrow(gene_info), "genes\n\n")
  
  # Merge with results
  significant_with_symbols <- merge(
    data.frame(ensembl_gene_id = rownames(significant_genes), significant_genes, stringsAsFactors = FALSE),
    gene_info,
    by = "ensembl_gene_id",
    all.x = TRUE
  )
  
  # Sort by p-value again
  significant_with_symbols <- significant_with_symbols[order(significant_with_symbols$padj), ]
  
}, error = function(e) {
  cat("❌ Error connecting to Ensembl:", e$message, "\n")
  cat("📝 Creating results without gene symbols...\n")
  
  significant_with_symbols <- data.frame(
    ensembl_gene_id = rownames(significant_genes),
    significant_genes,
    external_gene_name = "Symbol_unavailable",
    description = "Description_unavailable",
    stringsAsFactors = FALSE
  )
})

# =============================================================================
# STEP 5: TOP 10 SIGNIFICANT GENES WITH SYMBOLS
# =============================================================================

cat("🏆 STEP 5: TOP 10 MOST SIGNIFICANT GENES\n")
cat("-" , rep("-", 50), "\n", sep="")

cat("📊 MC9 vs MLM Differential Expression Results\n")
cat("🔍 Significance criteria: padj < 0.05, |FC| ≥ 1.5x\n")
cat("📈 Direction: Positive = Higher in MLM, Negative = Higher in MC9\n\n")

# Display top 10 with detailed information
top_10 <- head(significant_with_symbols, 10)

cat("Rank | Ensembl ID        | Gene Symbol | Fold Change | Direction | Adj P-value   | Description\n")
cat(rep("-", 100), "\n", sep="")

for (i in 1:nrow(top_10)) {
  ensembl <- top_10$ensembl_gene_id[i]
  symbol <- ifelse(is.na(top_10$external_gene_name[i]) | top_10$external_gene_name[i] == "", 
                   "Unknown", top_10$external_gene_name[i])
  log2fc <- top_10$log2FoldChange[i]
  fold_change <- round(2^abs(log2fc), 2)
  direction <- ifelse(log2fc > 0, "MLM ↑", "MC9 ↑")
  padj <- formatC(top_10$padj[i], format = "e", digits = 2)
  description <- ifelse(is.na(top_10$description[i]), "No description", 
                       substr(top_10$description[i], 1, 40))
  
  cat(sprintf("%4d | %-17s | %-11s | %8.2fx | %6s | %11s | %s\n", 
              i, ensembl, symbol, fold_change, direction, padj, description))
}

cat("\n")

# =============================================================================
# STEP 6: DETAILED STATISTICS AND VALIDATION
# =============================================================================

cat("📊 STEP 6: DETAILED STATISTICAL SUMMARY\n")
cat("-" , rep("-", 50), "\n", sep="")

# Direction analysis
up_in_mlm <- sum(significant_genes$log2FoldChange > 0, na.rm = TRUE)
up_in_mc9 <- sum(significant_genes$log2FoldChange < 0, na.rm = TRUE)

cat("🔢 Statistical breakdown:\n")
cat("   - Total significant genes:", nrow(significant_genes), "\n")
cat("   - Higher in MLM:", up_in_mlm, "(", round(up_in_mlm/nrow(significant_genes)*100, 1), "%)\n")
cat("   - Higher in MC9:", up_in_mc9, "(", round(up_in_mc9/nrow(significant_genes)*100, 1), "%)\n\n")

# Fold change distribution
fc_ranges <- c(
  "1.5-2x" = sum(abs(significant_genes$log2FoldChange) >= log2(1.5) & abs(significant_genes$log2FoldChange) < log2(2)),
  "2-3x" = sum(abs(significant_genes$log2FoldChange) >= log2(2) & abs(significant_genes$log2FoldChange) < log2(3)),
  "3-5x" = sum(abs(significant_genes$log2FoldChange) >= log2(3) & abs(significant_genes$log2FoldChange) < log2(5)),
  ">5x" = sum(abs(significant_genes$log2FoldChange) >= log2(5))
)

cat("📈 Fold change distribution:\n")
for (i in 1:length(fc_ranges)) {
  cat("   -", names(fc_ranges)[i], ":", fc_ranges[i], "genes\n")
}

cat("\n")

# P-value ranges
pval_ranges <- c(
  "< 1e-50" = sum(significant_genes$padj < 1e-50, na.rm = TRUE),
  "1e-50 to 1e-20" = sum(significant_genes$padj >= 1e-50 & significant_genes$padj < 1e-20, na.rm = TRUE),
  "1e-20 to 1e-10" = sum(significant_genes$padj >= 1e-20 & significant_genes$padj < 1e-10, na.rm = TRUE),
  "1e-10 to 0.001" = sum(significant_genes$padj >= 1e-10 & significant_genes$padj < 0.001, na.rm = TRUE),
  "0.001 to 0.05" = sum(significant_genes$padj >= 0.001 & significant_genes$padj < 0.05, na.rm = TRUE)
)

cat("🎯 P-value distribution:\n")
for (i in 1:length(pval_ranges)) {
  cat("   -", names(pval_ranges)[i], ":", pval_ranges[i], "genes\n")
}

cat("\n")

# =============================================================================
# STEP 7: SAVE RESULTS
# =============================================================================

cat("💾 STEP 7: SAVING DETAILED RESULTS\n")
cat("-" , rep("-", 50), "\n", sep="")

# Save comprehensive results
output_file <- "MC9_vs_MLM_detailed_results.csv"
write.csv(significant_with_symbols, output_file, row.names = FALSE)

cat("✅ Detailed results saved to:", output_file, "\n")
cat("   - Includes: Ensembl IDs, gene symbols, fold changes, p-values, descriptions\n")
cat("   - Total records:", nrow(significant_with_symbols), "\n\n")

# Summary for expert validation
cat("🧠 EXPERT VALIDATION SUMMARY:\n")
cat("-" , rep("-", 50), "\n", sep="")

cat("📋 Key findings for expert review:\n")
cat("1. **", nrow(significant_genes), "significant genes** identified (", 
    round(nrow(significant_genes) / sum(!is.na(res$padj)) * 100, 1), "% of tested)\n")
cat("2. **Balanced regulation**: ", up_in_mlm, " up in MLM, ", up_in_mc9, " up in MC9\n")
cat("3. **Strong effects**: ", fc_ranges["3-5x"] + fc_ranges[">5x"], " genes with >3x fold change\n")
cat("4. **High confidence**: ", pval_ranges["< 1e-50"], " genes with p.adj < 1e-50\n")
cat("5. **Top gene**: ", ifelse(exists("top_10") && nrow(top_10) > 0, 
                                paste(top_10$external_gene_name[1], "(", round(2^abs(top_10$log2FoldChange[1]), 2), "x change)"), 
                                "See results above"), "\n\n")

cat("=" , rep("=", 70), "\n", sep="")
cat("🎉 COMPLETE ANALYSIS FINISHED WITH FULL TRANSPARENCY!\n")
cat("🔍 All code shown, gene symbols retrieved, top genes identified\n")
cat("📊 Ready for expert biological validation!\n")
cat("=" , rep("=", 70), "\n", sep="")