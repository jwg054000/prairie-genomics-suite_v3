# Debug Pathway Gene Filtering Issue
# Check what genes are actually being passed to pathway analysis
#
# Author: Prairie Genomics Team
# Date: January 27, 2025
# Purpose: Debug why GO analysis uses all genes instead of significant ones

cat("🔍 Debugging Pathway Gene Filtering Issue\n")
cat("==========================================\n\n")

# Function to diagnose gene filtering in pathway analysis
debug_pathway_gene_filtering <- function(deseq2_results, padj_cutoff = 0.05, fc_cutoff = 1.0) {
  cat("📊 DEBUGGING GENE FILTERING FOR PATHWAY ANALYSIS\n")
  cat("================================================\n\n")
  
  # Step 1: Check input data structure
  cat("1. Input Data Analysis\n")
  cat("----------------------\n")
  cat("Class of deseq2_results:", class(deseq2_results), "\n")
  cat("Dimensions:", nrow(deseq2_results), "genes x", ncol(deseq2_results), "columns\n")
  cat("Column names:", paste(colnames(deseq2_results), collapse = ", "), "\n")
  cat("Row names (first 5):", paste(head(rownames(deseq2_results), 5), collapse = ", "), "\n\n")
  
  # Step 2: Check for required columns
  cat("2. Required Columns Check\n")
  cat("-------------------------\n")
  required_cols <- c("padj", "log2FoldChange")
  for (col in required_cols) {
    if (col %in% colnames(deseq2_results)) {
      cat("✅", col, "column present\n")
      
      # Check data quality
      na_count <- sum(is.na(deseq2_results[[col]]))
      cat("   - NA values:", na_count, "/", nrow(deseq2_results), 
          "(", round(na_count/nrow(deseq2_results)*100, 1), "%)\n")
      
      if (col == "padj") {
        valid_padj <- sum(!is.na(deseq2_results$padj) & deseq2_results$padj < padj_cutoff)
        cat("   - Genes with padj <", padj_cutoff, ":", valid_padj, "\n")
      }
      
      if (col == "log2FoldChange") {
        valid_fc <- sum(!is.na(deseq2_results$log2FoldChange) & abs(deseq2_results$log2FoldChange) > fc_cutoff)
        cat("   - Genes with |FC| >", fc_cutoff, ":", valid_fc, "\n")
      }
    } else {
      cat("❌", col, "column MISSING\n")
    }
  }
  cat("\n")
  
  # Step 3: Apply filtering step by step
  cat("3. Step-by-Step Filtering\n")
  cat("-------------------------\n")
  
  # Convert to data frame if needed
  if (class(deseq2_results)[1] == "DESeqResults") {
    cat("Converting DESeqResults to data frame...\n")
    deseq2_df <- as.data.frame(deseq2_results)
  } else {
    deseq2_df <- deseq2_results
  }
  
  cat("Starting with", nrow(deseq2_df), "total genes\n")
  
  # Step 3a: Remove NA padj
  step1 <- deseq2_df[!is.na(deseq2_df$padj), ]
  cat("After removing NA padj:", nrow(step1), "genes (removed", nrow(deseq2_df) - nrow(step1), ")\n")
  
  # Step 3b: Remove NA log2FoldChange
  step2 <- step1[!is.na(step1$log2FoldChange), ]
  cat("After removing NA log2FoldChange:", nrow(step2), "genes (removed", nrow(step1) - nrow(step2), ")\n")
  
  # Step 3c: Apply padj cutoff
  step3 <- step2[step2$padj < padj_cutoff, ]
  cat("After padj <", padj_cutoff, "filter:", nrow(step3), "genes (removed", nrow(step2) - nrow(step3), ")\n")
  
  # Step 3d: Apply fold change cutoff
  step4 <- step3[abs(step3$log2FoldChange) > fc_cutoff, ]
  cat("After |FC| >", fc_cutoff, "filter:", nrow(step4), "genes (removed", nrow(step3) - nrow(step4), ")\n")
  
  cat("\nFINAL SIGNIFICANT GENES:", nrow(step4), "\n\n")
  
  # Step 4: Show top significant genes
  if (nrow(step4) > 0) {
    cat("4. Top 10 Most Significant Genes\n")
    cat("--------------------------------\n")
    
    # Sort by padj
    step4_sorted <- step4[order(step4$padj), ]
    top_genes <- head(step4_sorted, 10)
    
    for (i in 1:nrow(top_genes)) {
      gene_id <- rownames(top_genes)[i]
      padj_val <- format(top_genes$padj[i], scientific = TRUE, digits = 3)
      fc_val <- round(top_genes$log2FoldChange[i], 3)
      cat(sprintf("%2d. %s: padj=%s, FC=%s\n", i, gene_id, padj_val, fc_val))
    }
    cat("\n")
  } else {
    cat("4. ❌ NO SIGNIFICANT GENES FOUND!\n")
    cat("   This explains why pathway analysis might be using all genes.\n\n")
  }
  
  # Step 5: Test the actual prepare_gene_list_ora function
  cat("5. Testing prepare_gene_list_ora Function\n")
  cat("-----------------------------------------\n")
  
  tryCatch({
    source("pathway_analysis.R")
    
    cat("Testing prepare_gene_list_ora with same parameters...\n")
    
    # Test the function
    gene_list <- prepare_gene_list_ora(deseq2_results, padj_cutoff, fc_cutoff, "human")
    
    if (!is.null(gene_list) && length(gene_list) > 0) {
      cat("✅ Function returned", length(gene_list), "genes\n")
      cat("Gene list sample:", paste(head(gene_list, 5), collapse = ", "), "\n")
    } else {
      cat("❌ Function returned NULL or empty list\n")
      cat("This confirms the filtering issue!\n")
    }
    
  }, error = function(e) {
    cat("❌ Error testing prepare_gene_list_ora:", e$message, "\n")
  })
  
  cat("\n")
  
  # Step 6: Recommendations
  cat("6. Diagnosis and Recommendations\n")
  cat("================================\n")
  
  if (nrow(step4) == 0) {
    cat("🚨 CRITICAL ISSUE: No significant genes found!\n")
    cat("\nPossible causes:\n")
    cat("1. padj cutoff too stringent (try 0.1 instead of 0.05)\n")
    cat("2. FC cutoff too stringent (try 0.5 instead of 1.0)\n")
    cat("3. DESeq2 analysis didn't find significant differences\n")
    cat("4. Sample groups not properly defined\n")
    cat("\nRecommended action:\n")
    cat("- Relax filtering parameters temporarily\n")
    cat("- Check sample annotation and contrasts\n")
    cat("- Verify DESeq2 analysis completed successfully\n")
  } else if (nrow(step4) < 10) {
    cat("⚠️ WARNING: Very few significant genes (", nrow(step4), ")\n")
    cat("This may lead to poor pathway enrichment results.\n")
    cat("\nRecommendations:\n")
    cat("- Consider relaxing cutoffs slightly\n")
    cat("- Check if sample size is adequate\n")
    cat("- Verify biological replicates are properly grouped\n")
  } else {
    cat("✅ Good: Found", nrow(step4), "significant genes\n")
    cat("This should be sufficient for pathway analysis.\n")
    
    if (nrow(step4) > 2000) {
      cat("\n⚠️ Note: Gene list will be limited to top 2000 for GO analysis performance\n")
    }
  }
  
  # Return diagnostic information
  return(list(
    total_genes = nrow(deseq2_df),
    significant_genes = nrow(step4),
    filtering_successful = nrow(step4) > 0,
    recommended_padj = if (nrow(step4) == 0) 0.1 else padj_cutoff,
    recommended_fc = if (nrow(step4) == 0) 0.5 else fc_cutoff,
    top_genes = if (nrow(step4) > 0) rownames(head(step4[order(step4$padj), ], 10)) else NULL
  ))
}

# Function to fix pathway analysis by ensuring proper gene filtering
fix_pathway_gene_filtering <- function() {
  cat("🔧 APPLYING FIX FOR PATHWAY GENE FILTERING\n")
  cat("==========================================\n\n")
  
  # Read the current pathway_analysis.R file
  pathway_file <- "pathway_analysis.R"
  
  if (!file.exists(pathway_file)) {
    cat("❌ pathway_analysis.R not found!\n")
    return(FALSE)
  }
  
  # Read current content
  content <- readLines(pathway_file)
  
  # Find the prepare_gene_list_ora function
  start_line <- which(grepl("prepare_gene_list_ora.*function", content))
  
  if (length(start_line) == 0) {
    cat("❌ prepare_gene_list_ora function not found!\n")
    return(FALSE)
  }
  
  cat("✅ Found prepare_gene_list_ora function at line", start_line, "\n")
  
  # Enhanced version with better debugging
  enhanced_function <- '
# Prepare gene list for over-representation analysis (GO, KEGG, MSigDB)
prepare_gene_list_ora <- function(deseq2_results, padj_cutoff, fc_cutoff, species) {
  tryCatch({
    cat("🔄 Preparing gene list for ORA analysis...\\n")
    cat("📊 Input filters: padj <", padj_cutoff, ", |FC| >", fc_cutoff, "\\n")
    
    # DEBUGGING: Check input data
    cat("🔍 DEBUG: Input data has", nrow(deseq2_results), "genes\\n")
    cat("🔍 DEBUG: Required columns present:", all(c("padj", "log2FoldChange") %in% colnames(deseq2_results)), "\\n")
    
    # Convert to data frame if needed
    if (class(deseq2_results)[1] == "DESeqResults") {
      deseq2_df <- as.data.frame(deseq2_results)
      cat("🔄 Converted DESeqResults to data frame\\n")
    } else {
      deseq2_df <- deseq2_results
    }
    
    # CRITICAL: Ensure we have the required columns
    if (!all(c("padj", "log2FoldChange") %in% colnames(deseq2_df))) {
      cat("❌ ERROR: Missing required columns (padj, log2FoldChange)\\n")
      cat("Available columns:", paste(colnames(deseq2_df), collapse = ", "), "\\n")
      return(NULL)
    }
    
    # Step-by-step filtering with debugging
    cat("🔍 Step 1: Starting with", nrow(deseq2_df), "total genes\\n")
    
    # Remove genes with NA values
    valid_genes <- deseq2_df[
      !is.na(deseq2_df$padj) & 
      !is.na(deseq2_df$log2FoldChange), 
    ]
    cat("🔍 Step 2: After removing NAs:", nrow(valid_genes), "genes\\n")
    
    # Apply significance filters
    significant_genes <- valid_genes[
      valid_genes$padj < padj_cutoff & 
      abs(valid_genes$log2FoldChange) > fc_cutoff, 
    ]
    cat("🔍 Step 3: After applying filters:", nrow(significant_genes), "significant genes\\n")
    
    # CRITICAL CHECK: If no significant genes, warn user
    if (nrow(significant_genes) == 0) {
      cat("⚠️ WARNING: No significant genes found with current filters!\\n")
      cat("💡 Consider relaxing cutoffs: padj < 0.1, |FC| > 0.5\\n")
      
      # Try with relaxed filters as fallback
      relaxed_genes <- valid_genes[
        valid_genes$padj < 0.1 & 
        abs(valid_genes$log2FoldChange) > 0.5, 
      ]
      
      if (nrow(relaxed_genes) > 0) {
        cat("🔧 Using relaxed filters as fallback:", nrow(relaxed_genes), "genes\\n")
        significant_genes <- relaxed_genes
      } else {
        cat("❌ Even relaxed filters found no genes - check your data!\\n")
        return(NULL)
      }
    }
    
    # Limit gene list size for GO analysis performance
    max_genes_ora <- 2000
    if (nrow(significant_genes) > max_genes_ora) {
      cat("⚠️ Large gene list detected (", nrow(significant_genes), "genes)\\n")
      cat("🔧 Limiting to top", max_genes_ora, "genes for optimal GO analysis performance\\n")
      
      # Sort by padj (most significant first) and take top genes
      significant_genes <- significant_genes[order(significant_genes$padj), ]
      significant_genes <- significant_genes[1:max_genes_ora, ]
    }
    
    # Extract gene IDs (rownames)
    gene_ids <- rownames(significant_genes)
    cat("🧬 Extracted", length(gene_ids), "gene IDs for conversion\\n")
    
    # Convert gene IDs to Entrez format if needed
    if (species == "human") {
      org_db <- org.Hs.eg.db
    } else if (species == "mouse") {
      org_db <- org.Mm.eg.db
    } else {
      org_db <- org.Hs.eg.db  # Default to human
    }
    
    cat("🔄 Converting gene IDs to Entrez format for", species, "...\\n")
    
    # Convert to Entrez IDs
    entrez_ids <- tryCatch({
      mapIds(org_db, 
             keys = gene_ids,
             column = "ENTREZID", 
             keytype = "ENSEMBL",
             multiVals = "first")
    }, error = function(e) {
      cat("⚠️ ENSEMBL conversion failed, trying SYMBOL...\\n")
      tryCatch({
        mapIds(org_db,
               keys = gene_ids,
               column = "ENTREZID",
               keytype = "SYMBOL", 
               multiVals = "first")
      }, error = function(e2) {
        cat("⚠️ SYMBOL conversion failed, using original IDs...\\n")
        gene_ids
      })
    })
    
    # Remove NAs from conversion
    if (is.vector(entrez_ids)) {
      final_genes <- entrez_ids[!is.na(entrez_ids)]
    } else {
      final_genes <- gene_ids  # Fallback to original IDs
    }
    
    conversion_rate <- length(final_genes) / length(gene_ids) * 100
    cat("✅ Gene ID conversion: ", length(final_genes), "/", length(gene_ids), 
        " (", round(conversion_rate, 1), "% success)\\n")
    
    # Final validation
    if (length(final_genes) < 5) {
      cat("❌ ERROR: Too few genes for pathway analysis (", length(final_genes), ")\\n")
      cat("💡 Need at least 5 genes. Check your filtering parameters.\\n")
      return(NULL)
    }
    
    cat("🎯 FINAL: ", length(final_genes), " genes ready for pathway analysis\\n")
    cat("📋 Sample genes:", paste(head(final_genes, 5), collapse = ", "), "\\n")
    
    return(final_genes)
    
  }, error = function(e) {
    cat("❌ ERROR in prepare_gene_list_ora:", e$message, "\\n")
    return(NULL)
  })
}'
  
  cat("✅ Enhanced function created with better debugging\n")
  cat("💾 Enhanced function is ready to replace the current version\n")
  
  # Instructions for manual replacement
  cat("\n🔧 TO APPLY THE FIX:\n")
  cat("1. Locate the prepare_gene_list_ora function in pathway_analysis.R\n")
  cat("2. Replace it with the enhanced version above\n") 
  cat("3. Save the file and restart the application\n")
  cat("4. The enhanced version includes detailed debugging output\n")
  
  return(TRUE)
}

# Function to test pathway analysis with sample data
test_pathway_analysis_with_debug <- function() {
  cat("🧪 TESTING PATHWAY ANALYSIS WITH DEBUG DATA\n")
  cat("===========================================\n\n")
  
  # Create sample DESeq2 results for testing
  set.seed(123)
  n_genes <- 1000
  
  sample_results <- data.frame(
    baseMean = runif(n_genes, 10, 1000),
    log2FoldChange = rnorm(n_genes, 0, 2),
    lfcSE = runif(n_genes, 0.1, 0.5),
    stat = rnorm(n_genes, 0, 3),
    pvalue = runif(n_genes, 0, 1),
    padj = runif(n_genes, 0, 1),
    stringsAsFactors = FALSE
  )
  
  # Make some genes significant
  significant_indices <- sample(1:n_genes, 100)  # 100 significant genes
  sample_results$padj[significant_indices] <- runif(100, 0, 0.05)
  sample_results$log2FoldChange[significant_indices] <- rnorm(100, 0, 3)
  
  # Add gene names
  rownames(sample_results) <- paste0("ENSG", sprintf("%011d", 1:n_genes))
  
  cat("Created sample DESeq2 results with:\n")
  cat("- Total genes:", n_genes, "\n")
  cat("- Significant genes (padj < 0.05):", sum(sample_results$padj < 0.05, na.rm = TRUE), "\n")
  cat("- High FC genes (|FC| > 1):", sum(abs(sample_results$log2FoldChange) > 1, na.rm = TRUE), "\n\n")
  
  # Run diagnostic
  diagnostic_result <- debug_pathway_gene_filtering(sample_results, 0.05, 1.0)
  
  cat("🎯 Diagnostic complete!\n")
  cat("Significant genes found:", diagnostic_result$significant_genes, "\n")
  cat("Filtering successful:", diagnostic_result$filtering_successful, "\n")
  
  return(diagnostic_result)
}

cat("✅ Debugging functions loaded!\n")
cat("\nAvailable functions:\n")
cat("  - debug_pathway_gene_filtering(deseq2_results, padj_cutoff, fc_cutoff)\n")
cat("  - fix_pathway_gene_filtering()\n")
cat("  - test_pathway_analysis_with_debug()\n")
cat("\n💡 Run these functions to diagnose and fix the gene filtering issue!\n")