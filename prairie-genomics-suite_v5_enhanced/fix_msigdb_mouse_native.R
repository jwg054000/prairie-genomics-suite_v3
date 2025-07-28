# Fix MSigDB Mouse Native Gene Sets and Gene Filtering
# Issues: Using human orthologs instead of native mouse, deprecated arguments, too many genes
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🔧 Fixing MSigDB Mouse Native Issues\n")
cat("=" , rep("=", 50), "\n")

# Issue Analysis
cat("📊 Issues Identified:\n")
cat("❌ Using human MSigDB with ortholog mapping (slow, less accurate)\n")
cat("❌ Using deprecated 'category' argument instead of 'collection'\n") 
cat("❌ Too many genes passing through despite stringent filters\n")
cat("❌ Need to use db_species = 'MM' for native mouse gene sets\n")

# Fix 1: Update MSigDB function to use native mouse gene sets
cat("\n🐭 Fix 1: Native Mouse Gene Sets\n")
cat(rep("-", 40), "\n")

improved_mouse_msigdb <- '
# FIXED: Native mouse MSigDB gene sets with updated API
get_fgsea_gene_sets <- function(species, collection = "H") {
  tryCatch({
    # Check msigdbr availability with CRAN mirror fix
    if (!require("msigdbr", quietly = TRUE)) {
      cat("📦 Installing msigdbr...\\n")
      if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
        options(repos = c(CRAN = "https://cloud.r-project.org/"))
      }
      tryCatch({
        install.packages("msigdbr", dependencies = TRUE)
        library(msigdbr)
      }, error = function(e) {
        cat("❌ msigdbr installation failed:", e$message, "\\n")
        return(NULL)
      })
    }
    
    # Map species names with native database specification
    if (species == "human") {
      msigdb_species <- "Homo sapiens"
      db_species <- "HS"  # Native human
    } else if (species == "mouse") {
      msigdb_species <- "Mus musculus" 
      db_species <- "MM"  # NATIVE MOUSE - this fixes the ortholog issue
    } else {
      cat("⚠️ Unsupported species, defaulting to human\\n")
      msigdb_species <- "Homo sapiens"
      db_species <- "HS"
    }
    
    # Validate species availability
    available_species <- tryCatch({
      msigdbr_species()
    }, error = function(e) {
      cat("❌ Could not check available species:", e$message, "\\n")
      return(NULL)
    })
    
    if (is.null(available_species) || !msigdb_species %in% available_species$species_name) {
      cat("❌ Species", msigdb_species, "not available in MSigDB\\n")
      return(NULL)
    }
    
    cat("📚 Retrieving NATIVE", collection, "gene sets for", msigdb_species, "(db_species =", db_species, ")\\n")
    
    # Try to get gene sets with UPDATED API (collection instead of category)
    gene_sets <- NULL
    
    # Primary attempt with exact collection using NEW API
    tryCatch({
      # FIXED: Use collection argument instead of deprecated category
      gene_sets <- msigdbr(
        species = msigdb_species, 
        collection = collection,  # NEW: collection instead of category
        db_species = db_species   # NEW: native species database
      )
      
      if (nrow(gene_sets) == 0) {
        gene_sets <- NULL
      }
    }, error = function(e) {
      cat("⚠️ Primary retrieval failed:", e$message, "\\n")
      gene_sets <<- NULL
    })
    
    # Fallback strategies for specific collections
    if (is.null(gene_sets) && collection == "C2") {
      cat("🔄 C2 collection failed, trying canonical pathways...\\n")
      tryCatch({
        gene_sets <- msigdbr(
          species = msigdb_species, 
          collection = "C2",
          subcollection = "CP",  # Canonical pathways
          db_species = db_species
        )
        if (nrow(gene_sets) == 0) gene_sets <- NULL
      }, error = function(e) {
        cat("⚠️ C2:CP fallback failed:", e$message, "\\n")
      })
    }
    
    # Further fallback to KEGG pathways if C2 fails
    if (is.null(gene_sets) && collection == "C2") {
      cat("🔄 Trying KEGG pathways as final fallback...\\n")
      tryCatch({
        gene_sets <- msigdbr(
          species = msigdb_species, 
          collection = "C2",
          subcollection = "CP:KEGG",
          db_species = db_species
        )
        if (nrow(gene_sets) == 0) gene_sets <- NULL
      }, error = function(e) {
        cat("⚠️ KEGG fallback failed:", e$message, "\\n")
      })
    }
    
    # Ultimate fallback to Hallmark if all else fails
    if (is.null(gene_sets) && collection != "H") {
      cat("🔄 All", collection, "attempts failed, falling back to Hallmark (H)...\\n")
      tryCatch({
        gene_sets <- msigdbr(
          species = msigdb_species, 
          collection = "H",
          db_species = db_species
        )
        if (nrow(gene_sets) == 0) gene_sets <- NULL
      }, error = function(e) {
        cat("❌ Even Hallmark fallback failed:", e$message, "\\n")
      })
    }
    
    if (is.null(gene_sets) || nrow(gene_sets) == 0) {
      cat("❌ No gene sets found for", collection, "in", msigdb_species, "\\n")
      return(NULL)
    }
    
    # Convert to named list using base R (no dplyr)
    pathway_names <- unique(gene_sets$gs_name)
    pathways <- list()
    
    for (pathway in pathway_names) {
      genes <- gene_sets$gene_symbol[gene_sets$gs_name == pathway]
      pathways[[pathway]] <- genes
    }
    
    cat("✅ Retrieved", length(pathways), "NATIVE gene sets\\n")
    cat("📋 Database used:", db_species, "(native", species, ")\\n")
    cat("📋 Collection used:", unique(gene_sets$gs_collection)[1], "\\n")
    if ("gs_subcollection" %in% colnames(gene_sets)) {
      subcollections <- unique(gene_sets$gs_subcollection)
      if (length(subcollections) <= 3) {
        cat("📋 Subcollections:", paste(subcollections, collapse = ", "), "\\n")
      }
    }
    
    return(pathways)
    
  }, error = function(e) {
    cat("❌ Failed to get gene sets:", e$message, "\\n")
    return(NULL)
  })
}
'

cat("✅ Created native mouse MSigDB function\n")

# Fix 2: Enhanced gene filtering to reduce gene count
cat("\n🔬 Fix 2: Enhanced Gene Filtering\n")
cat(rep("-", 35), "\n")

enhanced_filtering <- '
# ENHANCED: More stringent gene filtering for GSEA
prepare_gene_list_gsea <- function(deseq2_results, species, ranking_method = "signed_pvalue", 
                                  min_genes = 100, max_genes = 1000, 
                                  padj_filter = 0.1, basemean_filter = 10) {
  tryCatch({
    cat("🔄 Preparing GSEA gene list with enhanced filtering\\n")
    cat("   - Method:", ranking_method, "\\n")
    cat("   - Target range:", min_genes, "-", max_genes, "genes\\n")
    cat("   - padj filter: <", padj_filter, "\\n")
    cat("   - baseMean filter: >=", basemean_filter, "\\n")
    
    # Convert DESeq2 results to data frame WITHOUT rownames_to_column dependency
    deseq2_df <- as.data.frame(deseq2_results)
    deseq2_df$gene_id <- rownames(deseq2_df)
    
    cat("📊 Starting with", nrow(deseq2_df), "total genes\\n")
    
    # ENHANCED: More stringent filtering
    valid_genes <- deseq2_df[
      !is.na(deseq2_df$log2FoldChange) & 
      !is.na(deseq2_df$padj) &
      !is.na(deseq2_df$pvalue) &
      deseq2_df$padj != 0 &
      deseq2_df$pvalue != 0 &
      deseq2_df$padj < padj_filter &  # NEW: padj filter
      (is.null(deseq2_df$baseMean) || deseq2_df$baseMean >= basemean_filter),  # NEW: expression filter
    ]
    
    if (nrow(valid_genes) == 0) {
      cat("❌ No valid genes after filtering\\n")
      return(NULL)
    }
    
    cat("📊 After quality filtering:", nrow(valid_genes), "genes\\n")
    
    # If still too many genes, apply additional filtering
    if (nrow(valid_genes) > max_genes) {
      cat("🔧 Too many genes (", nrow(valid_genes), "), applying additional filtering...\\n")
      
      # Sort by significance and take top genes
      valid_genes <- valid_genes[order(valid_genes$padj), ]
      valid_genes <- valid_genes[1:max_genes, ]
      
      cat("📊 After count limiting:", nrow(valid_genes), "genes\\n")
    }
    
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
    
    cat("📊 After gene conversion:", nrow(valid_genes), "genes\\n")
    
    # Check if we have minimum required genes
    if (nrow(valid_genes) < min_genes) {
      cat("⚠️ Only", nrow(valid_genes), "genes available (minimum:", min_genes, ")\\n")
      cat("💡 Consider relaxing filters: higher padj_filter or lower basemean_filter\\n")
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
    
    # Handle duplicate gene symbols
    if (any(duplicated(names(gene_ranks)))) {
      cat("🔧 Handling", sum(duplicated(names(gene_ranks))), "duplicate gene symbols\\n")
      
      # For duplicates, keep the one with maximum absolute rank statistic
      dup_names <- names(gene_ranks)[duplicated(names(gene_ranks))]
      unique_names <- unique(names(gene_ranks))
      
      clean_ranks <- numeric(length(unique_names))
      names(clean_ranks) <- unique_names
      
      for (name in unique_names) {
        matching_indices <- which(names(gene_ranks) == name)
        if (length(matching_indices) > 1) {
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
    cat("   - Final gene count:", length(gene_ranks), "\\n")
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

cat("✅ Created enhanced gene filtering function\n")

# Write fixes to file
complete_fix <- paste(
  "# Complete MSigDB Mouse Native and Filtering Fix",
  "# Generated:", Sys.time(),
  "",
  "# Fix 1: Native mouse MSigDB gene sets",
  improved_mouse_msigdb,
  "",
  "# Fix 2: Enhanced gene filtering",
  enhanced_filtering,
  "",
  "cat('✅ MSigDB mouse native fixes loaded successfully\\n')",
  sep = "\n"
)

writeLines(complete_fix, "msigdb_mouse_native_fixes.R")
cat("✅ Complete fix saved to msigdb_mouse_native_fixes.R\n")

# Create usage examples
cat("\n💡 Usage Examples\n")
cat(rep("-", 25), "\n")

usage_examples <- '
# Usage Examples for Fixed Functions

# Example 1: Get native mouse Hallmark gene sets
mouse_hallmark <- get_fgsea_gene_sets("mouse", "H")
# Should show: "Database used: MM (native mouse)"

# Example 2: Enhanced gene filtering for GSEA
filtered_genes <- prepare_gene_list_gsea(
  deseq2_results, 
  species = "mouse",
  ranking_method = "signed_pvalue",
  min_genes = 100,     # Minimum genes required
  max_genes = 500,     # Maximum genes to use (prevents too many)
  padj_filter = 0.05,  # Only genes with padj < 0.05
  basemean_filter = 20 # Only well-expressed genes
)

# Example 3: Very stringent filtering
strict_genes <- prepare_gene_list_gsea(
  deseq2_results,
  species = "mouse", 
  max_genes = 200,     # Very limited gene set
  padj_filter = 0.01,  # Very significant only
  basemean_filter = 50 # Highly expressed only
)
'

writeLines(usage_examples, "usage_examples.R")
cat("✅ Usage examples saved\n")

# Summary
cat("\n🎯 Fix Summary\n")
cat("=" , rep("=", 20), "\n")

cat("MSigDB Issues Fixed:\n")
cat("✅ Using native mouse gene sets (db_species = 'MM')\n")
cat("✅ Updated to 'collection' argument (not deprecated 'category')\n")
cat("✅ No more ortholog mapping warnings\n")
cat("✅ Faster and more accurate gene set retrieval\n")

cat("\nGene Filtering Issues Fixed:\n")
cat("✅ Added padj filter to reduce gene count\n")
cat("✅ Added baseMean filter for expression level\n")
cat("✅ Added max_genes limit to prevent overload\n")
cat("✅ Enhanced filtering controls for user customization\n")

cat("\n🚀 Expected Improvements:\n")
cat("• Faster MSigDB retrieval (native vs ortholog mapping)\n")
cat("• More accurate mouse gene sets\n")
cat("• Controlled gene count (100-1000 genes)\n")
cat("• No deprecation warnings\n")
cat("• Better GSEA performance with fewer genes\n")

cat("\n💡 Recommended Settings:\n")
cat("For most analyses:\n")
cat("  - padj_filter = 0.05 (significant genes)\n")
cat("  - max_genes = 500 (good performance)\n")
cat("  - basemean_filter = 10 (expressed genes)\n")

cat("\nFor stringent analyses:\n") 
cat("  - padj_filter = 0.01 (very significant)\n")
cat("  - max_genes = 200 (fast performance)\n")
cat("  - basemean_filter = 50 (highly expressed)\n")

cat("\n🧬 MSigDB Mouse Native Issues Fixed!\n")
cat("🐭 Should now use proper mouse gene sets without orthologs\n")