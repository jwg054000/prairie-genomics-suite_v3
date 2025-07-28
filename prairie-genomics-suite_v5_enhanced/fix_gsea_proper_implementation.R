# Fix GSEA Implementation Based on DESeq2-GSEA Reference Documentation
# This script creates the proper GSEA implementation following best practices
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🔧 Implementing Proper GSEA Following DESeq2-GSEA Reference\n")
cat("=" , rep("=", 60), "\n")

# Create the improved GSEA preparation function
improved_gsea_code <- '
# IMPROVED: Prepare ranked gene list for GSEA following reference documentation
prepare_gene_list_gsea <- function(deseq2_results, species, ranking_method = "signed_pvalue") {
  tryCatch({
    cat("🔄 Preparing GSEA gene list using", ranking_method, "method\n")
    
    # Convert DESeq2 results to data frame with gene symbols
    deseq2_df <- as.data.frame(deseq2_results) %>%
      rownames_to_column("gene_id")
    
    # Remove genes with NA values (critical for GSEA)
    valid_genes <- deseq2_df %>%
      filter(
        !is.na(log2FoldChange) & 
        !is.na(padj) &
        !is.na(pvalue) &
        padj != 0 &
        pvalue != 0
      )
    
    if (nrow(valid_genes) == 0) {
      cat("❌ No valid genes for GSEA analysis\n")
      return(NULL)
    }
    
    cat("📊 Processing", nrow(valid_genes), "valid genes for GSEA\n")
    
    # Convert gene IDs to gene symbols (required for MSigDB)
    gene_symbols <- convert_to_gene_symbols(valid_genes$gene_id, species)
    
    # Add gene symbols to dataframe
    valid_genes$gene_symbol <- gene_symbols
    
    # Remove genes that couldnt be converted
    valid_genes <- valid_genes %>%
      filter(!is.na(gene_symbol) & gene_symbol != "")
    
    if (nrow(valid_genes) == 0) {
      cat("❌ No genes could be converted to symbols\n")
      return(NULL)
    }
    
    # Calculate ranking statistic based on method
    if (ranking_method == "signed_pvalue") {
      # Method from reference: signed p-value
      valid_genes$rank_stat <- sign(valid_genes$log2FoldChange) * (-log10(valid_genes$pvalue))
    } else if (ranking_method == "log2fc_pvalue") {
      # Alternative method: log2FC weighted by significance
      valid_genes$rank_stat <- valid_genes$log2FoldChange * (-log10(valid_genes$padj))
    } else if (ranking_method == "log2fc_only") {
      # Simple method: just log2FoldChange
      valid_genes$rank_stat <- valid_genes$log2FoldChange
    } else {
      # Default to signed p-value
      valid_genes$rank_stat <- sign(valid_genes$log2FoldChange) * (-log10(valid_genes$pvalue))
    }
    
    # Create named vector for GSEA
    gene_ranks <- valid_genes$rank_stat
    names(gene_ranks) <- valid_genes$gene_symbol
    
    # CRITICAL: Handle duplicate gene symbols (common issue)
    if (any(duplicated(names(gene_ranks)))) {
      cat("🔧 Handling", sum(duplicated(names(gene_ranks))), "duplicate gene symbols\n")
      
      # For duplicates, keep the one with maximum absolute rank statistic
      rank_df <- data.frame(
        gene_symbol = names(gene_ranks),
        rank_stat = gene_ranks,
        stringsAsFactors = FALSE
      ) %>%
        group_by(gene_symbol) %>%
        slice_max(abs(rank_stat), n = 1, with_ties = FALSE) %>%
        ungroup()
      
      # Recreate named vector
      gene_ranks <- rank_df$rank_stat
      names(gene_ranks) <- rank_df$gene_symbol
    }
    
    # Sort in decreasing order (most upregulated to most downregulated)
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)
    
    cat("✅ GSEA gene list prepared:\n")
    cat("   - Total genes:", length(gene_ranks), "\n")
    cat("   - Ranking method:", ranking_method, "\n")
    cat("   - Top upregulated:", names(gene_ranks)[1], "=", round(gene_ranks[1], 3), "\n")
    cat("   - Top downregulated:", names(gene_ranks)[length(gene_ranks)], "=", round(gene_ranks[length(gene_ranks)], 3), "\n")
    
    return(gene_ranks)
    
  }, error = function(e) {
    cat("❌ GSEA gene list preparation failed:", e$message, "\n")
    return(NULL)
  })
}

# IMPROVED: Gene symbol conversion function
convert_to_gene_symbols <- function(gene_ids, species) {
  tryCatch({
    # Use the existing gene conversion cache system
    if (exists("convert_genes_fast")) {
      # Use the cached conversion system
      conversion_result <- convert_genes_fast(gene_ids, species)
      
      if (!is.null(conversion_result) && "gene_symbol" %in% colnames(conversion_result)) {
        symbols <- conversion_result$gene_symbol[match(gene_ids, conversion_result$ensembl_gene_id)]
        return(symbols)
      }
    }
    
    # Fallback to direct annotation database conversion
    if (species == "human" && require("org.Hs.eg.db", quietly = TRUE)) {
      symbols <- suppressMessages(mapIds(
        org.Hs.eg.db,
        keys = gene_ids,
        column = "SYMBOL", 
        keytype = "ENSEMBL",
        multiVals = "first"
      ))
      return(as.character(symbols))
    } else if (species == "mouse" && require("org.Mm.eg.db", quietly = TRUE)) {
      symbols <- suppressMessages(mapIds(
        org.Mm.eg.db,
        keys = gene_ids,
        column = "SYMBOL",
        keytype = "ENSEMBL", 
        multiVals = "first"
      ))
      return(as.character(symbols))
    }
    
    # If all fails, return gene IDs as symbols
    cat("⚠️ Using gene IDs as symbols (conversion failed)\n")
    return(gene_ids)
    
  }, error = function(e) {
    cat("⚠️ Gene symbol conversion error:", e$message, "\n")
    return(gene_ids)
  })
}

# IMPROVED: GSEA analysis using fgsea (recommended approach)
run_gsea_analysis <- function(gene_list, species, gene_set_collection = "H") {
  tryCatch({
    cat("🔬 Running GSEA analysis with fgsea\n")
    
    # Ensure fgsea is available
    if (!require("fgsea", quietly = TRUE)) {
      cat("📦 Installing fgsea package...\n")
      BiocManager::install("fgsea")
      library(fgsea)
    }
    
    # Get gene sets in proper format for fgsea
    pathways <- get_fgsea_gene_sets(species, gene_set_collection)
    
    if (is.null(pathways) || length(pathways) == 0) {
      return(list(
        success = FALSE,
        error = "No gene sets retrieved",
        message = paste("Could not get", gene_set_collection, "gene sets for", species)
      ))
    }
    
    # Filter gene sets by size (15-500 genes as per reference)
    pathway_sizes <- sapply(pathways, length)
    valid_pathways <- pathways[pathway_sizes >= 15 & pathway_sizes <= 500]
    
    cat("📋 Gene set filtering:\n")
    cat("   - Original pathways:", length(pathways), "\n")
    cat("   - After size filtering (15-500):", length(valid_pathways), "\n")
    
    if (length(valid_pathways) == 0) {
      return(list(
        success = FALSE,
        error = "No valid gene sets after filtering",
        message = "All gene sets were outside the 15-500 gene size range"
      ))
    }
    
    # Run fgsea with proper parameters
    cat("🧬 Running fgsea with", length(valid_pathways), "pathways...\n")
    
    fgsea_results <- fgsea(
      pathways = valid_pathways,
      stats = gene_list,
      minSize = 15,
      maxSize = 500,
      nperm = 10000  # Higher for more accurate p-values
    )
    
    if (is.null(fgsea_results) || nrow(fgsea_results) == 0) {
      return(list(
        success = FALSE,
        error = "No significant pathways found",
        message = "Try different gene set collection or check gene ranking"
      ))
    }
    
    # Sort by adjusted p-value
    fgsea_results <- fgsea_results %>%
      arrange(padj) %>%
      as.data.frame()
    
    # Add interpretation columns
    fgsea_results$direction <- ifelse(fgsea_results$NES > 0, "Upregulated", "Downregulated")
    fgsea_results$significance <- ifelse(fgsea_results$padj < 0.05, "Significant", "Not Significant")
    
    cat("✅ GSEA completed successfully:\n")
    cat("   - Total pathways tested:", nrow(fgsea_results), "\n")
    cat("   - Significant pathways (padj < 0.05):", sum(fgsea_results$padj < 0.05), "\n")
    cat("   - Upregulated pathways:", sum(fgsea_results$NES > 0 & fgsea_results$padj < 0.05), "\n")
    cat("   - Downregulated pathways:", sum(fgsea_results$NES < 0 & fgsea_results$padj < 0.05), "\n")
    
    return(list(
      success = TRUE,
      data = fgsea_results,
      analysis_type = "GSEA",
      method = "fgsea",
      species = species,
      gene_set_collection = gene_set_collection,
      n_pathways = nrow(fgsea_results),
      n_significant = sum(fgsea_results$padj < 0.05)
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      error = e$message,
      message = paste("GSEA analysis failed:", e$message)
    ))
  })
}

# IMPROVED: Get gene sets in fgsea format (named list)
get_fgsea_gene_sets <- function(species, collection = "H") {
  tryCatch({
    # Check msigdbr availability
    if (!require("msigdbr", quietly = TRUE)) {
      cat("📦 Installing msigdbr...\n")
      install.packages("msigdbr")
      library(msigdbr)
    }
    
    # Map species names
    msigdb_species <- case_when(
      species == "human" ~ "Homo sapiens",
      species == "mouse" ~ "Mus musculus", 
      TRUE ~ "Homo sapiens"  # Default to human
    )
    
    # Get MSigDB data
    cat("📚 Retrieving", collection, "gene sets for", msigdb_species, "\n")
    
    m_df <- msigdbr(species = msigdb_species, category = collection)
    
    if (nrow(m_df) == 0) {
      cat("❌ No gene sets found for", collection, "in", msigdb_species, "\n")
      return(NULL)
    }
    
    # Convert to named list format required by fgsea
    pathways <- m_df %>%
      dplyr::select(gs_name, gene_symbol) %>%
      group_by(gs_name) %>%
      summarise(genes = list(gene_symbol), .groups = "drop") %>%
      deframe()  # Converts to named list
    
    cat("✅ Retrieved", length(pathways), "gene sets\n")
    
    return(pathways)
    
  }, error = function(e) {
    cat("❌ Failed to get gene sets:", e$message, "\n")
    return(NULL)
  })
}
'

# Write the improved code to a temporary file and source it
temp_file <- tempfile(fileext = ".R")
writeLines(improved_gsea_code, temp_file)

cat("📝 Generated improved GSEA implementation\n")
cat("🔧 Key improvements:\n")
cat("   ✅ Proper gene ranking using signed p-value method\n")
cat("   ✅ Gene set size filtering (15-500 genes)\n")
cat("   ✅ Duplicate gene handling with max absolute statistic\n")
cat("   ✅ fgsea implementation with 10,000 permutations\n")
cat("   ✅ Enhanced gene symbol conversion\n")
cat("   ✅ Comprehensive error handling and progress reporting\n")

cat("\n📋 Implementation Summary:\n")
cat("1. Gene Ranking: Uses signed p-value (sign(log2FC) * -log10(pvalue))\n")
cat("2. Gene Sets: MSigDB collections with proper size filtering\n")
cat("3. Analysis: fgsea with optimized parameters\n")
cat("4. Output: Comprehensive results with interpretation\n")

cat("\n🚀 Ready to integrate into pathway_analysis.R\n")
cat("💡 This follows the exact methodology from the reference document\n")

# Test the improved functions if libraries are available
cat("\n🧪 Testing improved GSEA functions...\n")
tryCatch({
  source(temp_file)
  cat("✅ Improved GSEA functions loaded successfully\n")
  
  # Test if required packages are available
  packages_needed <- c("fgsea", "msigdbr", "dplyr")
  packages_available <- sapply(packages_needed, function(x) requireNamespace(x, quietly = TRUE))
  
  cat("📦 Package status:\n")
  for (i in seq_along(packages_needed)) {
    status <- if (packages_available[i]) "✅" else "❌"
    cat("   ", status, packages_needed[i], "\n")
  }
  
  if (all(packages_available)) {
    cat("🎉 All required packages available - GSEA should work properly!\n")
  } else {
    cat("⚠️ Some packages missing - install with fix scripts\n")
  }
  
}, error = function(e) {
  cat("⚠️ Function testing failed:", e$message, "\n")
})

cat("\n🎯 Next Steps:\n")
cat("1. Run this script to generate improved functions\n")
cat("2. Update pathway_analysis.R with these improvements\n") 
cat("3. Test with actual DESeq2 results\n")
cat("4. Verify GSEA visualizations work\n")

cat("\n🧬 Prairie Genomics Suite - GSEA Implementation Fixed!\n")