# Fix GSEA and KEGG Performance Issues
# 1. Fix missing rownames_to_column function
# 2. Optimize KEGG analysis performance
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🔧 Fixing GSEA and KEGG Performance Issues\n")
cat("=" , rep("=", 60), "\n")

# Issue 1: Fix missing rownames_to_column function
cat("\n📝 Issue 1: Missing rownames_to_column Function\n")
cat(rep("-", 45), "\n")

cat("The error occurs because rownames_to_column is from tibble package\n")
cat("Solution: Add tibble dependency or use base R alternative\n")

# Create improved GSEA function without tibble dependency
improved_gsea_prep <- '
# FIXED: Prepare ranked gene list for GSEA (no tibble dependency)
prepare_gene_list_gsea <- function(deseq2_results, species, ranking_method = "signed_pvalue") {
  tryCatch({
    cat("🔄 Preparing GSEA gene list using", ranking_method, "method\\n")
    
    # Convert DESeq2 results to data frame WITHOUT rownames_to_column
    deseq2_df <- as.data.frame(deseq2_results)
    deseq2_df$gene_id <- rownames(deseq2_df)  # Add gene_id column manually
    
    # Remove genes with NA values (critical for GSEA)
    valid_genes <- deseq2_df[
      !is.na(deseq2_df$log2FoldChange) & 
      !is.na(deseq2_df$padj) &
      !is.na(deseq2_df$pvalue) &
      deseq2_df$padj != 0 &
      deseq2_df$pvalue != 0,
    ]
    
    if (nrow(valid_genes) == 0) {
      cat("❌ No valid genes for GSEA analysis\\n")
      return(NULL)
    }
    
    cat("📊 Processing", nrow(valid_genes), "valid genes for GSEA\\n")
    
    # Convert gene IDs to gene symbols (required for MSigDB)
    gene_symbols <- convert_to_gene_symbols(valid_genes$gene_id, species)
    
    # Add gene symbols to dataframe
    valid_genes$gene_symbol <- gene_symbols
    
    # Remove genes that couldnt be converted
    valid_genes <- valid_genes[
      !is.na(valid_genes$gene_symbol) & valid_genes$gene_symbol != "",
    ]
    
    if (nrow(valid_genes) == 0) {
      cat("❌ No genes could be converted to symbols\\n")
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
      cat("🔧 Handling", sum(duplicated(names(gene_ranks))), "duplicate gene symbols\\n")
      
      # For duplicates, keep the one with maximum absolute rank statistic
      # Using base R approach without dplyr
      dup_names <- names(gene_ranks)[duplicated(names(gene_ranks))]
      unique_names <- unique(names(gene_ranks))
      
      clean_ranks <- numeric(length(unique_names))
      names(clean_ranks) <- unique_names
      
      for (name in unique_names) {
        matching_indices <- which(names(gene_ranks) == name)
        if (length(matching_indices) > 1) {
          # Take the one with maximum absolute value
          max_idx <- matching_indices[which.max(abs(gene_ranks[matching_indices]))]
          clean_ranks[name] <- gene_ranks[max_idx]
        } else {
          clean_ranks[name] <- gene_ranks[matching_indices]
        }
      }
      
      gene_ranks <- clean_ranks
    }
    
    # Sort in decreasing order (most upregulated to most downregulated)
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)
    
    cat("✅ GSEA gene list prepared:\\n")
    cat("   - Total genes:", length(gene_ranks), "\\n")
    cat("   - Ranking method:", ranking_method, "\\n")
    cat("   - Top upregulated:", names(gene_ranks)[1], "=", round(gene_ranks[1], 3), "\\n")
    cat("   - Top downregulated:", names(gene_ranks)[length(gene_ranks)], "=", round(gene_ranks[length(gene_ranks)], 3), "\\n")
    
    return(gene_ranks)
    
  }, error = function(e) {
    cat("❌ GSEA gene list preparation failed:", e$message, "\\n")
    return(NULL)
  })
}
'

cat("✅ Created improved GSEA function without tibble dependency\n")

# Issue 2: KEGG Performance Optimization
cat("\n⚡ Issue 2: KEGG Analysis Performance\n")
cat(rep("-", 45), "\n")

cat("Problems identified:\n")
cat("• Multiple simultaneous KEGG queries\n")
cat("• Large gene lists (895 genes) causing timeouts\n")
cat("• No caching of KEGG annotation data\n")
cat("• Duplicate analysis calls\n")

# Create optimized KEGG function
optimized_kegg <- '
# OPTIMIZED: KEGG analysis with performance improvements
run_kegg_analysis <- function(gene_list, species) {
  tryCatch({
    cat("🔍 Running optimized KEGG analysis...\\n")
    
    # Performance check - warn for large gene lists
    if (length(gene_list) > 500) {
      cat("⚠️ Large gene list (", length(gene_list), "genes) - KEGG analysis may be slow\\n")
      cat("💡 Consider using stricter filtering (lower padj, higher |FC|) for faster results\\n")
      
      # Option to limit genes for performance
      if (length(gene_list) > 1000) {
        cat("🔧 Limiting to top 1000 genes for performance\\n")
        gene_list <- gene_list[1:1000]
      }
    }
    
    # Map species to KEGG organism codes
    kegg_organism <- case_when(
      species == "human" ~ "hsa",
      species == "mouse" ~ "mmu",
      species == "rat" ~ "rno", 
      species == "fly" ~ "dme",
      species == "worm" ~ "cel",
      species == "yeast" ~ "sce",
      species == "zebrafish" ~ "dre",
      TRUE ~ "hsa"  # Default to human
    )
    
    cat("🌐 Querying KEGG database for", kegg_organism, "with", length(gene_list), "genes...\\n")
    cat("⏱️ This may take 30-60 seconds for large gene lists...\\n")
    
    # Run KEGG enrichment with timeout handling
    kegg_result <- tryCatch({
      # Set timeout for KEGG query (2 minutes max)
      setTimeLimit(cpu = 120, elapsed = 120)
      
      enrichKEGG(
        gene = gene_list,
        organism = kegg_organism,
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2,
        minGSSize = 10,  # Minimum pathway size
        maxGSSize = 500  # Maximum pathway size
      )
      
    }, error = function(e) {
      if (grepl("timeout|time limit", e$message, ignore.case = TRUE)) {
        cat("⏰ KEGG query timed out - try with fewer genes or check internet connection\\n")
      } else {
        cat("❌ KEGG query failed:", e$message, "\\n")
      }
      return(NULL)
    }, finally = {
      # Reset time limit
      setTimeLimit(cpu = Inf, elapsed = Inf)
    })
    
    if (is.null(kegg_result)) {
      return(list(
        success = FALSE,
        error = "KEGG query failed or timed out",
        message = "Try with fewer genes or check internet connection. KEGG can be slow with large gene lists."
      ))
    }
    
    if (nrow(kegg_result@result) == 0) {
      return(list(
        success = FALSE,
        error = "No significant KEGG pathways found",
        message = "Try relaxing thresholds or using a different gene set"
      ))
    }
    
    kegg_df <- as.data.frame(kegg_result)
    
    cat("✅ KEGG analysis completed:\\n")
    cat("   - Pathways found:", nrow(kegg_df), "\\n")
    cat("   - Significant (padj < 0.05):", sum(kegg_df$p.adjust < 0.05), "\\n")
    
    return(list(
      success = TRUE,
      data = kegg_df,
      enrichment_object = kegg_result,
      analysis_type = "KEGG",
      species = species,
      organism_code = kegg_organism,
      n_pathways = nrow(kegg_df)
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      error = e$message,
      message = paste("KEGG analysis failed:", e$message)
    ))
  })
}
'

cat("✅ Created optimized KEGG function with timeout handling\n")

# Additional performance improvements
cat("\n🚀 Additional Performance Improvements\n")
cat(rep("-", 45), "\n")

performance_improvements <- '
# PERFORMANCE: Add analysis queue to prevent duplicate calls
analysis_queue <- new.env()
analysis_queue$running <- list()
analysis_queue$completed <- list()

# Queue management for pathway analysis
queue_analysis <- function(analysis_type, gene_list, species, ...) {
  # Create unique key for this analysis
  analysis_key <- paste(analysis_type, species, length(gene_list), digest::digest(gene_list), sep = "_")
  
  # Check if already running
  if (analysis_key %in% names(analysis_queue$running)) {
    cat("⏳ Analysis already in progress, please wait...\\n")
    return(list(success = FALSE, error = "Analysis in progress"))
  }
  
  # Check if recently completed (cache for 5 minutes)
  if (analysis_key %in% names(analysis_queue$completed)) {
    cache_time <- analysis_queue$completed[[analysis_key]]$timestamp
    if (difftime(Sys.time(), cache_time, units = "mins") < 5) {
      cat("💾 Returning cached result from", round(difftime(Sys.time(), cache_time, units = "mins"), 1), "minutes ago\\n")
      return(analysis_queue$completed[[analysis_key]]$result)
    }
  }
  
  # Mark as running
  analysis_queue$running[[analysis_key]] <- Sys.time()
  
  # Run analysis
  result <- tryCatch({
    if (analysis_type == "KEGG") {
      run_kegg_analysis(gene_list, species, ...)
    } else if (analysis_type == "GSEA") {
      run_gsea_analysis(gene_list, species, ...)
    } else if (analysis_type == "GO") {
      run_go_analysis(gene_list, species, ...)
    }
  }, finally = {
    # Remove from running queue
    if (analysis_key %in% names(analysis_queue$running)) {
      analysis_queue$running[[analysis_key]] <- NULL
    }
  })
  
  # Cache successful results
  if (!is.null(result) && result$success) {
    analysis_queue$completed[[analysis_key]] <- list(
      result = result,
      timestamp = Sys.time()
    )
  }
  
  return(result)
}

# Clear analysis cache (for debugging)
clear_analysis_cache <- function() {
  analysis_queue$running <- list()
  analysis_queue$completed <- list()
  cat("🗑️ Analysis cache cleared\\n")
}
'

cat("✅ Created analysis queue system to prevent duplicates\n")

# Create complete fix script
cat("\n📝 Creating Complete Fix\n")
cat(rep("-", 30), "\n")

complete_fix <- paste(
  "# Complete GSEA and KEGG Performance Fix",
  "# Generated:", Sys.time(),
  "",
  "# Fix 1: GSEA without tibble dependency", 
  improved_gsea_prep,
  "",
  "# Fix 2: Optimized KEGG analysis",
  optimized_kegg, 
  "",
  "# Fix 3: Analysis queue system",
  performance_improvements,
  "",
  "cat('✅ All fixes loaded successfully\\n')",
  sep = "\n"
)

# Write complete fix to file
writeLines(complete_fix, "pathway_analysis_performance_fixes.R")

cat("✅ Complete fix saved to pathway_analysis_performance_fixes.R\n")

# Summary and recommendations
cat("\n🎯 Fix Summary\n")
cat("=" , rep("=", 30), "\n")

cat("GSEA Issues Fixed:\n")
cat("✅ Removed tibble::rownames_to_column dependency\n") 
cat("✅ Used base R rownames() instead\n")
cat("✅ Improved duplicate handling with base R\n")

cat("\nKEGG Performance Issues Fixed:\n")
cat("✅ Added timeout handling (2 minutes max)\n")
cat("✅ Gene list size warnings and limits\n") 
cat("✅ Analysis queue to prevent duplicates\n")
cat("✅ Result caching for 5 minutes\n")

cat("\n💡 Usage Recommendations:\n")
cat("• For KEGG: Use <500 genes for best performance\n")
cat("• Apply stricter filtering: padj < 0.01, |FC| > 1.5\n")
cat("• Try GO analysis first (much faster)\n")
cat("• Use GSEA/MSigDB as good speed/quality compromise\n")

cat("\n🚀 Integration Steps:\n")
cat("1. Source the performance fixes file\n")
cat("2. Update pathway_analysis.R with new functions\n")
cat("3. Test with smaller gene lists first\n")
cat("4. Monitor for duplicate analysis calls\n")

cat("\n⚠️ Important Notes:\n")
cat("• KEGG queries external databases - slow is normal\n")
cat("• Internet connection required for KEGG\n") 
cat("• Large gene lists (>1000) will be automatically limited\n")
cat("• Analysis results cached for 5 minutes\n")

cat("\n🧬 Prairie Genomics Suite - Performance Issues Fixed!\n")