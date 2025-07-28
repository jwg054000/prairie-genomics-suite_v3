# Complete Fix for MSigDB API and GSEA Logical Errors
# Generated:
2025-07-24 15:55:23.682665

# Set CRAN mirror
options(repos = c(CRAN = 'https://cloud.r-project.org/'))

# Fix 1: MSigDB function using working API

# FIXED: MSigDB function using working API
get_fgsea_gene_sets <- function(species, collection = "H") {
  tryCatch({
    # Check msigdbr availability with CRAN mirror fix
    if (!require("msigdbr", quietly = TRUE)) {
      cat("📦 Installing msigdbr...\n")
      if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
        options(repos = c(CRAN = "https://cloud.r-project.org/"))
      }
      install.packages("msigdbr", dependencies = TRUE)
      library(msigdbr)
    }
    
    # Map species names
    msigdb_species <- if (species == "human") {
      "Homo sapiens"
    } else if (species == "mouse") {
      "Mus musculus"
    } else {
      "Homo sapiens"  # Default to human
    }
    
    cat("📚 Retrieving", collection, "gene sets for", msigdb_species, "\n")
    
    # Try to get gene sets with fallbacks - using CATEGORY (not collection)
    gene_sets <- NULL
    
    # Primary attempt using category parameter (older but working API)
    tryCatch({
      gene_sets <- msigdbr(species = msigdb_species, category = collection)
      if (nrow(gene_sets) == 0) {
        gene_sets <- NULL
      }
    }, error = function(e) {
      cat("⚠️ Primary retrieval failed:", e$message, "\n")
      gene_sets <<- NULL
    })
    
    # Fallback strategies for specific collections
    if (is.null(gene_sets) && collection == "C2") {
      cat("🔄 C2 collection failed, trying C2:CP...\n")
      tryCatch({
        gene_sets <- msigdbr(species = msigdb_species, category = "C2", subcategory = "CP")
        if (nrow(gene_sets) == 0) gene_sets <- NULL
      }, error = function(e) {
        cat("⚠️ C2:CP fallback failed:", e$message, "\n")
      })
    }
    
    # Further fallback to KEGG pathways
    if (is.null(gene_sets) && collection == "C2") {
      cat("🔄 Trying C2:CP:KEGG...\n")
      tryCatch({
        gene_sets <- msigdbr(species = msigdb_species, category = "C2", subcategory = "CP:KEGG")
        if (nrow(gene_sets) == 0) gene_sets <- NULL
      }, error = function(e) {
        cat("⚠️ KEGG fallback failed:", e$message, "\n")
      })
    }
    
    # Ultimate fallback to Hallmark
    if (is.null(gene_sets) && collection != "H") {
      cat("🔄 Falling back to Hallmark (H)...\n")
      tryCatch({
        gene_sets <- msigdbr(species = msigdb_species, category = "H")
        if (nrow(gene_sets) == 0) gene_sets <- NULL
      }, error = function(e) {
        cat("❌ Hallmark fallback failed:", e$message, "\n")
      })
    }
    
    if (is.null(gene_sets) || nrow(gene_sets) == 0) {
      cat("❌ No gene sets found for", collection, "\n")
      return(NULL)
    }
    
    # Convert to named list using base R
    pathway_names <- unique(gene_sets$gs_name)
    pathways <- list()
    
    for (pathway in pathway_names) {
      genes <- gene_sets$gene_symbol[gene_sets$gs_name == pathway]
      pathways[[pathway]] <- genes
    }
    
    cat("✅ Retrieved", length(pathways), "gene sets\n")
    cat("📋 Collection used:", unique(gene_sets$gs_cat)[1], "\n")
    
    return(pathways)
    
  }, error = function(e) {
    cat("❌ Failed to get gene sets:", e$message, "\n")
    return(NULL)
  })
}


# Fix 2: GSEA preparation with fixed logical operations

# FIXED: GSEA preparation with corrected logical operations
prepare_gene_list_gsea <- function(deseq2_results, species, ranking_method = "signed_pvalue", 
                                  max_genes = 800, padj_filter = 0.1, basemean_filter = 10) {
  tryCatch({
    cat("🔄 Preparing GSEA gene list with enhanced filtering\n")
    cat("   - Method:", ranking_method, "\n")
    cat("   - Max genes:", max_genes, "\n")
    cat("   - padj filter: <", padj_filter, "\n")
    cat("   - baseMean filter: >=", basemean_filter, "\n")
    
    # Convert DESeq2 results to data frame
    deseq2_df <- as.data.frame(deseq2_results)
    deseq2_df$gene_id <- rownames(deseq2_df)
    
    cat("📊 Starting with", nrow(deseq2_df), "total genes\n")
    
    # FIXED: Proper logical operations (fix for coercion error)
    # Check if baseMean column exists
    has_basemean <- "baseMean" %in% colnames(deseq2_df)
    
    if (has_basemean) {
      # Enhanced filtering WITH baseMean
      valid_genes <- deseq2_df[
        !is.na(deseq2_df$log2FoldChange) & 
        !is.na(deseq2_df$padj) &
        !is.na(deseq2_df$pvalue) &
        deseq2_df$padj != 0 &
        deseq2_df$pvalue != 0 &
        deseq2_df$padj < padj_filter &
        deseq2_df$baseMean >= basemean_filter,  # FIXED: Direct comparison
      ]
    } else {
      # Enhanced filtering WITHOUT baseMean
      cat("⚠️ No baseMean column found, skipping expression filter\n")
      valid_genes <- deseq2_df[
        !is.na(deseq2_df$log2FoldChange) & 
        !is.na(deseq2_df$padj) &
        !is.na(deseq2_df$pvalue) &
        deseq2_df$padj != 0 &
        deseq2_df$pvalue != 0 &
        deseq2_df$padj < padj_filter,
      ]
    }
    
    if (nrow(valid_genes) == 0) {
      cat("❌ No valid genes after filtering\n")
      return(NULL)
    }
    
    cat("📊 After quality filtering:", nrow(valid_genes), "genes\n")
    
    # If still too many genes, limit by significance
    if (nrow(valid_genes) > max_genes) {
      cat("🔧 Too many genes (", nrow(valid_genes), "), limiting to top", max_genes, "\n")
      
      # Sort by significance and take top genes
      valid_genes <- valid_genes[order(valid_genes$padj), ]
      valid_genes <- valid_genes[1:max_genes, ]
      
      cat("📊 After count limiting:", nrow(valid_genes), "genes\n")
    }
    
    # Convert gene IDs to gene symbols
    gene_symbols <- convert_to_gene_symbols(valid_genes$gene_id, species)
    valid_genes$gene_symbol <- gene_symbols
    
    # Remove genes that couldnt be converted
    valid_genes <- valid_genes[
      !is.na(valid_genes$gene_symbol) & valid_genes$gene_symbol != "",
    ]
    
    if (nrow(valid_genes) == 0) {
      cat("❌ No genes could be converted to symbols\n")
      return(NULL)
    }
    
    cat("📊 After gene conversion:", nrow(valid_genes), "genes\n")
    
    # Calculate ranking statistic
    if (ranking_method == "signed_pvalue") {
      valid_genes$rank_stat <- sign(valid_genes$log2FoldChange) * (-log10(valid_genes$pvalue))
    } else if (ranking_method == "log2fc_pvalue") {
      valid_genes$rank_stat <- valid_genes$log2FoldChange * (-log10(valid_genes$padj))
    } else {
      valid_genes$rank_stat <- valid_genes$log2FoldChange
    }
    
    # Create named vector for GSEA
    gene_ranks <- valid_genes$rank_stat
    names(gene_ranks) <- valid_genes$gene_symbol
    
    # Handle duplicates
    if (any(duplicated(names(gene_ranks)))) {
      cat("🔧 Handling", sum(duplicated(names(gene_ranks))), "duplicate gene symbols\n")
      
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
    
    # Sort in decreasing order
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)
    
    cat("✅ GSEA gene list prepared:\n")
    cat("   - Final gene count:", length(gene_ranks), "\n")
    cat("   - Top upregulated:", names(gene_ranks)[1], "=", round(gene_ranks[1], 3), "\n")
    cat("   - Top downregulated:", names(gene_ranks)[length(gene_ranks)], "=", round(gene_ranks[length(gene_ranks)], 3), "\n")
    
    return(gene_ranks)
    
  }, error = function(e) {
    cat("❌ GSEA gene list preparation failed:", e$message, "\n")
    return(NULL)
  })
}


cat('✅ Both fixes loaded successfully\n')
