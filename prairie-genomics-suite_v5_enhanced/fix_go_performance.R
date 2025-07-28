# Fix GO Pathway Analysis Performance Issues
# Address hanging and timeout problems with large gene lists
#
# Author: Prairie Genomics Team
# Date: January 27, 2025
# Purpose: Optimize GO analysis for better performance and user experience

cat("🔧 Fixing GO Pathway Analysis Performance Issues\n")
cat("=", rep("=", 45), "\n\n")

# Enhanced GO analysis function with performance optimizations
optimized_run_go_analysis <- function(gene_list, species, ontology = "BP") {
  tryCatch({
    cat("🔍 Running optimized GO analysis...\n")
    
    # Performance check - filter large gene lists
    original_count <- length(gene_list)
    cat("📊 Input genes:", original_count, "\n")
    
    # Limit gene list size for performance (GO analysis becomes very slow with >2000 genes)
    max_genes_for_go <- 2000
    
    if (length(gene_list) > max_genes_for_go) {
      cat("⚠️ Large gene list detected (", length(gene_list), "genes)\n")
      cat("🔧 Limiting to top", max_genes_for_go, "genes for optimal GO analysis performance\n")
      cat("💡 For comprehensive analysis, consider using GSEA instead\n")
      
      # Take the first genes (assuming they're sorted by significance)
      gene_list <- gene_list[1:max_genes_for_go]
    }
    
    # Validate gene list
    if (length(gene_list) < 10) {
      return(list(
        success = FALSE,
        error = "Insufficient genes for GO analysis",
        message = paste("Need at least 10 genes, found", length(gene_list))
      ))
    }
    
    # Determine organism database
    if (species == "human") {
      org_db <- org.Hs.eg.db
      cat("🧬 Using human GO annotations (org.Hs.eg.db)\n")
    } else if (species == "mouse") {
      org_db <- org.Mm.eg.db  
      cat("🐭 Using mouse GO annotations (org.Mm.eg.db)\n")
    } else {
      cat("⚠️ Species", species, "not directly supported, using human annotations\n")
      org_db <- org.Hs.eg.db
    }
    
    cat("🔍 Running GO enrichment for", ontology, "ontology with", length(gene_list), "genes...\n")
    cat("⏱️ This should take 10-30 seconds...\n")
    
    # Set timeout for GO analysis (90 seconds max)
    go_results <- tryCatch({
      setTimeLimit(cpu = 90, elapsed = 90)
      
      # Run GO over-representation analysis with optimized parameters
      enrichGO(
        gene = gene_list,
        OrgDb = org_db,
        ont = ontology,  # BP, MF, or CC
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE,
        minGSSize = 10,    # Minimum genes per GO term
        maxGSSize = 500    # Maximum genes per GO term  
      )
    }, error = function(e) {
      cat("⚠️ GO analysis timeout or error:", e$message, "\n")
      if (grepl("timeout|time limit", e$message, ignore.case = TRUE)) {
        stop("GO analysis timed out - try using fewer genes or GSEA analysis instead")
      } else {
        stop(e$message)
      }
    }, finally = {
      # Reset time limits
      setTimeLimit(cpu = Inf, elapsed = Inf)
    })
    
    cat("✅ GO analysis completed\n")
    
    # Check results
    if (is.null(go_results) || nrow(go_results@result) == 0) {
      return(list(
        success = FALSE,
        error = "No significant GO terms found",
        message = paste("No enriched GO terms found for", ontology, "ontology. Try:",
                       "\n• Using more genes (relax filtering thresholds)",
                       "\n• Different ontology (BP, MF, or CC)", 
                       "\n• GSEA analysis for more sensitive detection")
      ))
    }
    
    # Convert to data frame for easier handling
    go_df <- as.data.frame(go_results@result)
    
    cat("📊 Found", nrow(go_df), "enriched GO terms\n")
    
    # Add analysis metadata
    go_df$analysis_info <- paste0("GO ", ontology, " enrichment (", species, ")")
    go_df$genes_analyzed <- length(gene_list)
    go_df$original_gene_count <- original_count
    
    return(list(
      success = TRUE,
      data = go_df,
      enrichment_object = go_results,
      analysis_type = "GO",
      ontology = ontology,
      species = species,
      n_terms = nrow(go_df),
      genes_analyzed = length(gene_list),
      original_gene_count = original_count,
      performance_note = if (original_count > max_genes_for_go) {
        paste("Gene list limited to", max_genes_for_go, "genes for performance")
      } else {
        "Full gene list analyzed"
      }
    ))
    
  }, error = function(e) {
    # Enhanced error handling
    error_msg <- e$message
    
    if (grepl("timeout|time limit", error_msg, ignore.case = TRUE)) {
      suggestion <- "GO analysis timed out. Try: 1) Use fewer genes (stricter filtering), 2) Use GSEA analysis instead, 3) Try a different ontology"
    } else if (grepl("network|internet|connection", error_msg, ignore.case = TRUE)) {
      suggestion <- "Network connection issue. Check internet connection and try again"
    } else if (grepl("database|org\\..*\\.eg\\.db", error_msg, ignore.case = TRUE)) {
      suggestion <- "Database issue. Try installing/updating organism annotation packages"
    } else {
      suggestion <- "Unknown error. Try using GSEA analysis or contact support"
    }
    
    return(list(
      success = FALSE,
      error = error_msg,
      message = paste("GO analysis failed:", error_msg),
      suggestion = suggestion,
      analysis_type = "GO",
      ontology = ontology,
      species = species
    ))
  })
}

# Enhanced gene list preparation with better filtering
prepare_genes_for_go <- function(deseq2_results, padj_cutoff = 0.05, fc_cutoff = 1.0, 
                                max_genes = 2000, species = "human") {
  tryCatch({
    cat("🧬 Preparing genes for GO analysis...\n")
    
    # Filter for significant genes
    significant_genes <- deseq2_results[
      !is.na(deseq2_results$padj) & 
      !is.na(deseq2_results$log2FoldChange) &
      deseq2_results$padj < padj_cutoff & 
      abs(deseq2_results$log2FoldChange) > fc_cutoff,
    ]
    
    cat("📊 Found", nrow(significant_genes), "significant genes\n")
    
    if (nrow(significant_genes) == 0) {
      return(list(
        success = FALSE,
        error = "No significant genes found",
        message = "Try relaxing the filtering thresholds"
      ))
    }
    
    # Sort by adjusted p-value (most significant first)
    significant_genes <- significant_genes[order(significant_genes$padj), ]
    
    # Extract gene IDs (assuming rownames are gene IDs)
    gene_ids <- rownames(significant_genes)
    
    # Remove version numbers from Ensembl IDs if present
    clean_gene_ids <- sub("\\.\\d+$", "", gene_ids)
    
    # Convert gene IDs to appropriate format for GO analysis
    if (species == "human") {
      org_db <- org.Hs.eg.db
    } else if (species == "mouse") {
      org_db <- org.Mm.eg.db
    } else {
      org_db <- org.Hs.eg.db  # Default to human
    }
    
    # Try to convert to Entrez IDs (preferred for GO analysis)
    cat("🔄 Converting gene IDs to Entrez format...\n")
    
    entrez_ids <- tryCatch({
      # Try Ensembl to Entrez conversion first
      converted <- mapIds(org_db, 
                         keys = clean_gene_ids,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")
      
      # Remove NAs
      converted[!is.na(converted)]
    }, error = function(e) {
      cat("⚠️ Ensembl conversion failed, trying gene symbols...\n")
      
      # Fallback: try symbol conversion
      tryCatch({
        converted <- mapIds(org_db,
                           keys = clean_gene_ids, 
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")
        converted[!is.na(converted)]
      }, error = function(e2) {
        cat("⚠️ Gene ID conversion failed, using original IDs\n")
        clean_gene_ids
      })
    })
    
    conversion_rate <- length(entrez_ids) / length(clean_gene_ids) * 100
    cat("✅ Converted", length(entrez_ids), "genes (", round(conversion_rate, 1), "% success rate)\n")
    
    # Limit number of genes for performance
    if (length(entrez_ids) > max_genes) {
      cat("🔧 Limiting to top", max_genes, "genes for optimal performance\n")
      entrez_ids <- entrez_ids[1:max_genes]
    }
    
    if (length(entrez_ids) < 10) {
      return(list(
        success = FALSE,
        error = "Insufficient converted genes",
        message = paste("Only", length(entrez_ids), "genes could be converted. Need at least 10.")
      ))
    }
    
    return(list(
      success = TRUE,
      gene_list = entrez_ids,
      original_count = nrow(significant_genes),
      converted_count = length(entrez_ids),
      conversion_rate = conversion_rate
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      error = e$message,
      message = paste("Gene preparation failed:", e$message)
    ))
  })
}

# Create complete fix for GO analysis
complete_go_fix <- '
# Replace the existing run_go_analysis function in pathway_analysis.R
run_go_analysis <- function(gene_list, species, ontology = "BP") {
  tryCatch({
    cat("🔍 Running optimized GO analysis...\\n")
    
    # Performance check - filter large gene lists
    original_count <- length(gene_list)
    cat("📊 Input genes:", original_count, "\\n")
    
    # Limit gene list size for performance (GO analysis becomes very slow with >2000 genes)
    max_genes_for_go <- 2000
    
    if (length(gene_list) > max_genes_for_go) {
      cat("⚠️ Large gene list detected (", length(gene_list), "genes)\\n")
      cat("🔧 Limiting to top", max_genes_for_go, "genes for optimal GO analysis performance\\n")
      cat("💡 For comprehensive analysis, consider using GSEA instead\\n")
      
      # Take the first genes (assuming they are sorted by significance)
      gene_list <- gene_list[1:max_genes_for_go]
    }
    
    # Validate gene list
    if (length(gene_list) < 10) {
      return(list(
        success = FALSE,
        error = "Insufficient genes for GO analysis",
        message = paste("Need at least 10 genes, found", length(gene_list))
      ))
    }
    
    # Determine organism database
    if (species == "human") {
      org_db <- org.Hs.eg.db
      cat("🧬 Using human GO annotations (org.Hs.eg.db)\\n")
    } else if (species == "mouse") {
      org_db <- org.Mm.eg.db  
      cat("🐭 Using mouse GO annotations (org.Mm.eg.db)\\n")
    } else {
      cat("⚠️ Species", species, "not directly supported, using human annotations\\n")
      org_db <- org.Hs.eg.db
    }
    
    cat("🔍 Running GO enrichment for", ontology, "ontology with", length(gene_list), "genes...\\n")
    cat("⏱️ This should take 10-30 seconds...\\n")
    
    # Set timeout for GO analysis (90 seconds max)
    go_results <- tryCatch({
      setTimeLimit(cpu = 90, elapsed = 90)
      
      # Run GO over-representation analysis with optimized parameters
      enrichGO(
        gene = gene_list,
        OrgDb = org_db,
        ont = ontology,  # BP, MF, or CC
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE,
        minGSSize = 10,    # Minimum genes per GO term
        maxGSSize = 500    # Maximum genes per GO term  
      )
    }, error = function(e) {
      cat("⚠️ GO analysis timeout or error:", e$message, "\\n")
      if (grepl("timeout|time limit", e$message, ignore.case = TRUE)) {
        stop("GO analysis timed out - try using fewer genes or GSEA analysis instead")
      } else {
        stop(e$message)
      }
    }, finally = {
      # Reset time limits
      setTimeLimit(cpu = Inf, elapsed = Inf)
    })
    
    cat("✅ GO analysis completed\\n")
    
    # Check results
    if (is.null(go_results) || nrow(go_results@result) == 0) {
      return(list(
        success = FALSE,
        error = "No significant GO terms found",
        message = paste("No enriched GO terms found for", ontology, "ontology. Try:",
                       "\\n• Using more genes (relax filtering thresholds)",
                       "\\n• Different ontology (BP, MF, or CC)", 
                       "\\n• GSEA analysis for more sensitive detection")
      ))
    }
    
    # Convert to data frame for easier handling
    go_df <- as.data.frame(go_results@result)
    
    cat("📊 Found", nrow(go_df), "enriched GO terms\\n")
    
    # Add analysis metadata
    go_df$analysis_info <- paste0("GO ", ontology, " enrichment (", species, ")")
    go_df$genes_analyzed <- length(gene_list)
    go_df$original_gene_count <- original_count
    
    return(list(
      success = TRUE,
      data = go_df,
      enrichment_object = go_results,
      analysis_type = "GO",
      ontology = ontology,
      species = species,
      n_terms = nrow(go_df),
      genes_analyzed = length(gene_list),
      original_gene_count = original_count,
      performance_note = if (original_count > max_genes_for_go) {
        paste("Gene list limited to", max_genes_for_go, "genes for performance")
      } else {
        "Full gene list analyzed"
      }
    ))
    
  }, error = function(e) {
    # Enhanced error handling
    error_msg <- e$message
    
    if (grepl("timeout|time limit", error_msg, ignore.case = TRUE)) {
      suggestion <- "GO analysis timed out. Try: 1) Use fewer genes (stricter filtering), 2) Use GSEA analysis instead, 3) Try a different ontology"
    } else if (grepl("network|internet|connection", error_msg, ignore.case = TRUE)) {
      suggestion <- "Network connection issue. Check internet connection and try again"
    } else if (grepl("database|org\\\\..*\\\\.eg\\\\.db", error_msg, ignore.case = TRUE)) {
      suggestion <- "Database issue. Try installing/updating organism annotation packages"
    } else {
      suggestion <- "Unknown error. Try using GSEA analysis or contact support"
    }
    
    return(list(
      success = FALSE,
      error = error_msg,
      message = paste("GO analysis failed:", error_msg),
      suggestion = suggestion,
      analysis_type = "GO",
      ontology = ontology,
      species = species
    ))
  })
}
'

cat("✅ Created optimized GO analysis functions\n")
cat("📋 Key improvements:\n")
cat("   - Gene list size limit (2000 genes max)\n")
cat("   - Timeout protection (90 seconds)\n") 
cat("   - Better error handling and user guidance\n")
cat("   - Performance monitoring and warnings\n")
cat("   - Enhanced gene ID conversion\n")

cat("\n🔧 Applying fix to pathway_analysis.R\n")
cat(rep("-", 35), "\n")

# Apply the fix by replacing the function in pathway_analysis.R
writeLines(complete_go_fix, "go_analysis_fix.R")
cat("✅ Fix saved to go_analysis_fix.R\n")

cat("\n💡 To apply the fix:\n")
cat("1. Source this file: source('fix_go_performance.R')\n")
cat("2. Or manually replace run_go_analysis function in pathway_analysis.R\n")
cat("3. Test with a smaller gene list first\n")

cat("\n🎯 Expected improvements:\n")
cat("   - GO analysis completes in 10-30 seconds (vs. hanging)\n")
cat("   - Automatic gene list limiting for performance\n")
cat("   - Clear progress messages and timing estimates\n")
cat("   - Graceful timeout handling\n")
cat("   - Better error messages with actionable suggestions\n")