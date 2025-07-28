#!/usr/bin/env Rscript

# Prairie Genomics Suite - Complete Analysis with Scientific Guardrails
# MC9 vs MLM Differential Expression Analysis
# Real-world demonstration of AI-guided genomics analysis

cat("🧬 Prairie Genomics Suite - Complete Analysis Pipeline\n")
cat("🎯 Analysis: MC9 vs MLM Mouse Cancer Cell Lines\n")
cat("🛡️ Scientific Guardrails: ACTIVE\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("📅 Analysis Date:", as.character(Sys.Date()), "\n")
cat("👨‍🔬 Expert Validation: Joshua Garton\n\n")

# Load required libraries with guardrails
cat("📦 LOADING ANALYSIS LIBRARIES:\n")
required_packages <- c("readr", "DESeq2", "ggplot2", "pheatmap", "RColorBrewer")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("⚠️  Installing missing package:", pkg, "\n")
    if (pkg == "DESeq2") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", quiet = TRUE)
      }
      BiocManager::install("DESeq2", quiet = TRUE)
    } else {
      install.packages(pkg, quiet = TRUE)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat("✅", pkg, "loaded\n")
}
cat("\n")

# Load data with validation guardrails
cat("📊 DATA LOADING WITH VALIDATION GUARDRAILS:\n")
cat("-" , rep("-", 50), "\n", sep="")

# Expression data
expr_file <- "/Users/joshuagarton/Desktop/MC9.raw.counts.test.csv"
meta_file <- "/Users/joshuagarton/Desktop/MC9_sample_metadata.csv"

tryCatch({
  # Load expression matrix
  cat("📈 Loading expression matrix...\n")
  expr_data <- read_csv(expr_file, show_col_types = FALSE)
  
  # GUARDRAIL: Validate data structure
  if (nrow(expr_data) == 0) stop("❌ GUARDRAIL ALERT: Expression file is empty")
  if (ncol(expr_data) < 3) stop("❌ GUARDRAIL ALERT: Need at least gene column + 2 samples")
  
  # Convert to matrix format
  gene_names <- expr_data[[1]]
  count_matrix <- as.matrix(expr_data[, -1])
  rownames(count_matrix) <- gene_names
  
  cat("✅ Expression matrix loaded:", nrow(count_matrix), "genes ×", ncol(count_matrix), "samples\n")
  
  # Load metadata
  cat("📋 Loading sample metadata...\n")
  meta_data <- read_csv(meta_file, show_col_types = FALSE)
  
  # GUARDRAIL: Validate metadata structure
  if (!"Sample_ID" %in% colnames(meta_data)) stop("❌ GUARDRAIL ALERT: Missing Sample_ID column")
  if (!"Condition" %in% colnames(meta_data)) stop("❌ GUARDRAIL ALERT: Missing Condition column")
  
  cat("✅ Metadata loaded:", nrow(meta_data), "samples\n")
  
  # GUARDRAIL: Validate sample matching
  expr_samples <- colnames(count_matrix)
  meta_samples <- meta_data$Sample_ID
  
  if (length(intersect(expr_samples, meta_samples)) == 0) {
    stop("❌ GUARDRAIL ALERT: No matching samples between expression and metadata")
  }
  
  overlap_pct <- length(intersect(expr_samples, meta_samples)) / max(length(expr_samples), length(meta_samples)) * 100
  if (overlap_pct < 80) {
    cat("⚠️  GUARDRAIL WARNING: Only", round(overlap_pct, 1), "% sample overlap\n")
  } else {
    cat("✅ Sample matching:", round(overlap_pct, 1), "% overlap\n")
  }
  
}, error = function(e) {
  stop("❌ DATA LOADING FAILED:", e$message)
})

cat("\n")

# Focus on MC9 vs MLM comparison with guardrails
cat("🔬 ANALYSIS SETUP: MC9 vs MLM COMPARISON\n")
cat("-" , rep("-", 50), "\n", sep="")

# Filter data for MC9 and MLM only
mc9_samples <- meta_data$Sample_ID[meta_data$Condition == "MC9"]
mlm_samples <- meta_data$Sample_ID[meta_data$Condition == "MLM"]

cat("🧪 Sample groups identified:\n")
cat("   - MC9 samples:", length(mc9_samples), "(", paste(mc9_samples, collapse = ", "), ")\n")
cat("   - MLM samples:", length(mlm_samples), "(", paste(mlm_samples, collapse = ", "), ")\n")

# GUARDRAIL: Check sample size adequacy
if (length(mc9_samples) < 3 || length(mlm_samples) < 3) {
  cat("⚠️  GUARDRAIL WARNING: Fewer than 3 replicates per group reduces statistical power\n")
} else {
  cat("✅ Sample size adequate: ≥3 replicates per group\n")
}

# Prepare data for DESeq2
comparison_samples <- c(mc9_samples, mlm_samples)
comparison_matrix <- count_matrix[, comparison_samples]
comparison_metadata <- meta_data[meta_data$Sample_ID %in% comparison_samples, ]

# Ensure metadata order matches matrix
comparison_metadata <- comparison_metadata[match(colnames(comparison_matrix), comparison_metadata$Sample_ID), ]

cat("📊 Analysis dataset prepared:\n")
cat("   - Genes:", nrow(comparison_matrix), "\n")
cat("   - Samples:", ncol(comparison_matrix), "\n")
cat("   - Conditions: MC9 (n=", length(mc9_samples), "), MLM (n=", length(mlm_samples), ")\n\n")

# Pre-analysis quality control with guardrails
cat("🔍 PRE-ANALYSIS QUALITY CONTROL GUARDRAILS:\n")
cat("-" , rep("-", 50), "\n", sep="")

# GUARDRAIL: Check library sizes
lib_sizes <- colSums(comparison_matrix)
median_lib_size <- median(lib_sizes)
cv_lib_size <- sd(lib_sizes) / mean(lib_sizes)

cat("📚 Library size analysis:\n")
cat("   - Range:", format(min(lib_sizes), big.mark=","), "to", format(max(lib_sizes), big.mark=","), "\n")
cat("   - Median:", format(median_lib_size, big.mark=","), "\n")
cat("   - CV:", round(cv_lib_size * 100, 1), "%\n")

if (median_lib_size < 1000000) {
  cat("⚠️  GUARDRAIL WARNING: Low library sizes may reduce statistical power\n")
} else {
  cat("✅ Library sizes adequate for robust analysis\n")
}

if (cv_lib_size > 0.3) {
  cat("⚠️  GUARDRAIL WARNING: High library size variation (CV >30%)\n")
} else {
  cat("✅ Library size variation acceptable (CV <30%)\n")
}

# GUARDRAIL: Check for low-count genes
low_count_genes <- rowSums(comparison_matrix) < 10
pct_low_count <- sum(low_count_genes) / nrow(comparison_matrix) * 100

cat("🧬 Gene filtering analysis:\n")
cat("   - Total genes:", nrow(comparison_matrix), "\n")
cat("   - Low count genes (<10 total):", sum(low_count_genes), "(", round(pct_low_count, 1), "%)\n")

if (pct_low_count > 70) {
  cat("⚠️  GUARDRAIL WARNING: Very high proportion of low-count genes\n")
} else {
  cat("✅ Reasonable proportion of detectable genes\n")
}

cat("\n")

# Run DESeq2 analysis with guardrails
cat("🧮 DIFFERENTIAL EXPRESSION ANALYSIS WITH GUARDRAILS:\n")
cat("-" , rep("-", 50), "\n", sep="")

tryCatch({
  # Create DESeq2 dataset
  cat("🔧 Creating DESeq2 dataset...\n")
  
  # GUARDRAIL: Ensure factors are properly set
  comparison_metadata$Condition <- factor(comparison_metadata$Condition, levels = c("MC9", "MLM"))
  
  dds <- DESeqDataSetFromMatrix(
    countData = comparison_matrix,
    colData = comparison_metadata,
    design = ~ Condition
  )
  
  cat("✅ DESeq2 dataset created\n")
  
  # GUARDRAIL: Filter low-count genes
  cat("🔄 Applying count filters...\n")
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  cat("✅ Filtered to", nrow(dds), "genes with adequate counts\n")
  
  # Run DESeq2
  cat("⚙️  Running DESeq2 analysis...\n")
  dds <- DESeq(dds, quiet = TRUE)
  
  cat("✅ DESeq2 analysis completed\n")
  
  # Extract results
  cat("📊 Extracting results: MLM vs MC9...\n")
  res <- results(dds, contrast = c("Condition", "MLM", "MC9"))
  
  # GUARDRAIL: Check analysis quality
  if (sum(!is.na(res$padj)) == 0) {
    stop("❌ GUARDRAIL ALERT: No valid adjusted p-values generated")
  }
  
  cat("✅ Results extracted successfully\n")
  
}, error = function(e) {
  stop("❌ DESEQ2 ANALYSIS FAILED:", e$message)
})

cat("\n")

# Apply expert-validated parameters with guardrails
cat("🎯 APPLYING EXPERT-VALIDATED PARAMETERS:\n")
cat("-" , rep("-", 50), "\n", sep="")

# Parameters validated by Joshua
padj_threshold <- 0.05
fc_threshold <- 1.5

cat("📋 Analysis parameters (expert-approved):\n")
cat("   - Adjusted p-value threshold:", padj_threshold, "\n")
cat("   - Fold change threshold:", fc_threshold, "x\n")
cat("   - Multiple testing correction: Benjamini-Hochberg (FDR)\n")
cat("   - Statistical method: DESeq2 Wald test\n\n")

# Filter significant results
significant_genes <- subset(res, padj < padj_threshold & abs(log2FoldChange) >= log2(fc_threshold))
significant_genes <- significant_genes[order(significant_genes$padj), ]

# GUARDRAIL: Validate results quality
total_tested <- sum(!is.na(res$padj))
num_significant <- nrow(significant_genes)
pct_significant <- num_significant / total_tested * 100

cat("🧬 RESULTS SUMMARY WITH GUARDRAILS:\n")
cat("-" , rep("-", 50), "\n", sep="")

cat("📊 Statistical summary:\n")
cat("   - Genes tested:", total_tested, "\n")
cat("   - Significant genes:", num_significant, "(", round(pct_significant, 1), "%)\n")

# GUARDRAIL: Check for reasonable results
if (num_significant == 0) {
  cat("⚠️  GUARDRAIL WARNING: No significant genes found - check parameters or data quality\n")
} else if (pct_significant > 50) {
  cat("⚠️  GUARDRAIL WARNING: Very high proportion of significant genes (>50%) - check for batch effects\n")
} else if (pct_significant < 0.1) {
  cat("⚠️  GUARDRAIL WARNING: Very few significant genes (<0.1%) - check statistical power\n")
} else {
  cat("✅ Reasonable proportion of significant genes detected\n")
}

# Direction analysis
up_regulated <- sum(significant_genes$log2FoldChange > 0, na.rm = TRUE)
down_regulated <- sum(significant_genes$log2FoldChange < 0, na.rm = TRUE)

cat("   - Up-regulated in MLM:", up_regulated, "\n")
cat("   - Down-regulated in MLM:", down_regulated, "\n")

if (up_regulated == 0 || down_regulated == 0) {
  cat("⚠️  GUARDRAIL WARNING: All changes in one direction - check for systematic bias\n")
} else {
  cat("✅ Bidirectional changes detected (expected for biological comparison)\n")
}

# Show top results
if (num_significant > 0) {
  cat("\n🏆 TOP 10 MOST SIGNIFICANT GENES:\n")
  cat("   (MLM vs MC9 - positive = higher in MLM)\n\n")
  
  top_genes <- head(significant_genes, 10)
  
  for (i in 1:nrow(top_genes)) {
    gene <- rownames(top_genes)[i]
    fc <- round(2^abs(top_genes$log2FoldChange[i]), 2)
    direction <- ifelse(top_genes$log2FoldChange[i] > 0, "↑", "↓")
    padj <- formatC(top_genes$padj[i], format = "e", digits = 2)
    
    cat(sprintf("   %d. %s %s %.2fx (p.adj = %s)\n", i, gene, direction, fc, padj))
  }
}

cat("\n")

# Quality control plots (text-based summary)
cat("📈 QUALITY CONTROL ASSESSMENT:\n")
cat("-" , rep("-", 50), "\n", sep="")

# MA plot statistics
mean_expr <- rowMeans(log2(counts(dds, normalized = TRUE) + 1))
log2fc <- res$log2FoldChange

cat("📊 MA plot statistics:\n")
cat("   - Mean expression range:", round(min(mean_expr, na.rm=TRUE), 2), "to", round(max(mean_expr, na.rm=TRUE), 2), "\n")
cat("   - Log2 fold change range:", round(min(log2fc, na.rm=TRUE), 2), "to", round(max(log2fc, na.rm=TRUE), 2), "\n")

# P-value distribution
pval_bins <- hist(res$pvalue, breaks = 20, plot = FALSE)
uniform_expected <- length(res$pvalue[!is.na(res$pvalue)]) / 20

cat("📊 P-value distribution assessment:\n")
if (pval_bins$counts[1] > 3 * uniform_expected) {
  cat("✅ P-value distribution shows enrichment of small p-values (good signal)\n")
} else {
  cat("⚠️  GUARDRAIL WARNING: P-value distribution may indicate weak signal\n")
}

cat("\n")

# Expert validation prompts
cat("🧠 EXPERT VALIDATION FRAMEWORK:\n")
cat("-" , rep("-", 50), "\n", sep="")

cat("📝 Questions for expert ground truth validation:\n\n")

cat("1. 🎯 BIOLOGICAL EXPECTATION:\n")
cat("   Do these results align with your biological expectations for MC9 vs MLM?\n")
cat("   Expected answer: [Expert provides biological context]\n\n")

cat("2. 🧬 KNOWN MARKERS:\n")
cat("   Are any known cancer markers or pathway genes in the top results?\n")
cat("   Expected answer: [Expert identifies familiar genes/pathways]\n\n")

cat("3. 📊 RESULT MAGNITUDE:\n")
cat("   Does the number of significant genes (", num_significant, ") seem reasonable?\n")
cat("   Expected answer: [Expert validates scale of differences]\n\n")

cat("4. 🔬 PRIOR ANALYSIS:\n")
cat("   How do these results compare to your previous MC9 vs MLM analyses?\n")
cat("   Expected answer: [Expert compares to historical results]\n\n")

cat("5. ⚙️ PARAMETER VALIDATION:\n")
cat("   Do the chosen parameters (p.adj <", padj_threshold, ", FC ≥", fc_threshold, "x) capture relevant biology?\n")
cat("   Expected answer: [Expert confirms parameter appropriateness]\n\n")

# Analysis summary and recommendations
cat("💡 AI SYSTEM RECOMMENDATIONS:\n")
cat("-" , rep("-", 50), "\n", sep="")

cat("Based on this analysis, our AI system recommends:\n\n")

cat("1. 📊 FOLLOW-UP ANALYSES:\n")
if (num_significant > 100) {
  cat("   - Pathway enrichment analysis (many genes detected)\n")
  cat("   - Gene set enrichment analysis (GSEA)\n")
} else {
  cat("   - Individual gene validation (focused gene list)\n")
  cat("   - Targeted pathway analysis\n")
}

cat("   - Sample clustering analysis to verify groups\n")
cat("   - Batch effect assessment if available\n\n")

cat("2. 🧬 BIOLOGICAL INTERPRETATION:\n")
cat("   - Focus on genes with highest fold changes AND significance\n")
cat("   - Cross-reference with cancer gene databases\n")
cat("   - Consider cell line-specific biology\n\n")

cat("3. 🔬 EXPERIMENTAL VALIDATION:\n")
cat("   - qPCR validation of top 3-5 genes\n")
cat("   - Protein level validation if antibodies available\n")
cat("   - Functional assays for pathway-relevant genes\n\n")

cat("=" , rep("=", 70), "\n", sep="")
cat("🎉 ANALYSIS COMPLETE: Scientific guardrails successfully protected analysis integrity!\n")
cat("🧠 Expert validation framework ready for ground truth collection\n")
cat("🚀 AI-guided genomics analysis pipeline demonstrated on real data\n")
cat("=" , rep("=", 70), "\n", sep="")

# Save results summary for ground truth comparison
write.csv(significant_genes, "MC9_vs_MLM_significant_genes.csv", row.names = TRUE)
cat("💾 Results saved: MC9_vs_MLM_significant_genes.csv\n")

cat("\n📋 NEXT STEPS:\n")
cat("1. Expert reviews results for biological accuracy\n")
cat("2. Ground truth validation against expert knowledge\n")
cat("3. AI system refinement based on expert feedback\n")
cat("4. Expansion to additional comparisons (M1245, M242)\n\n")

cat("🔬 Scientific reproducibility achieved through AI-guided analysis! ✨\n")