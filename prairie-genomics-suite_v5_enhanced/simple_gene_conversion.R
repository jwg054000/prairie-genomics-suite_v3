# Simple Gene Conversion - Local Bioconductor Only
# Fast, reliable gene ID to symbol conversion without external dependencies
# 
# This avoids BioMart connection issues and works offline
# Uses org.Hs.eg.db and org.Mm.eg.db for fast local conversion

# Simple gene conversion function using only local packages
convert_genes_simple <- function(ensembl_ids, species = "human") {
  tryCatch({
    cat("🧬 Starting simple gene conversion for", length(ensembl_ids), "genes\n")
    cat("   - Species:", species, "\n")
    
    # Initialize result data frame
    result_df <- data.frame(
      ensembl_gene_id = ensembl_ids,
      gene_symbol = NA_character_,
      stringsAsFactors = FALSE
    )
    
    if (species == "human") {
      # Try human conversion with org.Hs.eg.db
      if (requireNamespace("org.Hs.eg.db", quietly = TRUE) && 
          requireNamespace("AnnotationDbi", quietly = TRUE)) {
        
        cat("✅ Using org.Hs.eg.db for human gene conversion\n")
        
        # Remove version numbers from Ensembl IDs
        clean_ids <- sub("\\.\\d+$", "", ensembl_ids)
        
        # Get mappings
        mapping <- tryCatch({
          AnnotationDbi::select(
            org.Hs.eg.db::org.Hs.eg.db,
            keys = clean_ids,
            columns = c("ENSEMBL", "SYMBOL"),
            keytype = "ENSEMBL"
          )
        }, error = function(e) {
          cat("⚠️ org.Hs.eg.db query failed:", e$message, "\n")
          NULL
        })
        
        if (!is.null(mapping) && nrow(mapping) > 0) {
          # Map symbols back to original IDs
          for (i in 1:nrow(result_df)) {
            clean_id <- sub("\\.\\d+$", "", ensembl_ids[i])
            match_idx <- which(mapping$ENSEMBL == clean_id)
            if (length(match_idx) > 0 && !is.na(mapping$SYMBOL[match_idx[1]])) {
              result_df$gene_symbol[i] <- mapping$SYMBOL[match_idx[1]]
            }
          }
        }
        
      } else {
        cat("⚠️ org.Hs.eg.db not available - symbols will be NA\n")
      }
      
    } else if (species == "mouse") {
      # Try mouse conversion with org.Mm.eg.db
      if (requireNamespace("org.Mm.eg.db", quietly = TRUE) && 
          requireNamespace("AnnotationDbi", quietly = TRUE)) {
        
        cat("✅ Using org.Mm.eg.db for mouse gene conversion\n")
        
        # Remove version numbers from Ensembl IDs
        clean_ids <- sub("\\.\\d+$", "", ensembl_ids)
        
        # Get mappings
        mapping <- tryCatch({
          AnnotationDbi::select(
            org.Mm.eg.db::org.Mm.eg.db,
            keys = clean_ids,
            columns = c("ENSEMBL", "SYMBOL"),
            keytype = "ENSEMBL"
          )
        }, error = function(e) {
          cat("⚠️ org.Mm.eg.db query failed:", e$message, "\n")
          NULL
        })
        
        if (!is.null(mapping) && nrow(mapping) > 0) {
          # Map symbols back to original IDs
          for (i in 1:nrow(result_df)) {
            clean_id <- sub("\\.\\d+$", "", ensembl_ids[i])
            match_idx <- which(mapping$ENSEMBL == clean_id)
            if (length(match_idx) > 0 && !is.na(mapping$SYMBOL[match_idx[1]])) {
              result_df$gene_symbol[i] <- mapping$SYMBOL[match_idx[1]]
            }
          }
        }
        
      } else {
        cat("⚠️ org.Mm.eg.db not available - symbols will be NA\n")
      }
    }
    
    # Report conversion statistics
    converted_count <- sum(!is.na(result_df$gene_symbol) & result_df$gene_symbol != "")
    conversion_rate <- round(100 * converted_count / length(ensembl_ids), 1)
    
    cat("📊 Simple conversion completed:\n")
    cat("   - Total genes:", length(ensembl_ids), "\n")
    cat("   - Successfully converted:", converted_count, "\n")
    cat("   - Conversion rate:", conversion_rate, "%\n")
    
    return(result_df)
    
  }, error = function(e) {
    cat("❌ Simple gene conversion failed:", e$message, "\n")
    
    # Return fallback result
    return(data.frame(
      ensembl_gene_id = ensembl_ids,
      gene_symbol = NA_character_,
      stringsAsFactors = FALSE
    ))
  })
}

# Test if required packages are available
test_gene_conversion_packages <- function() {
  cat("🧪 Testing gene conversion package availability:\n")
  
  # Test human packages
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) && 
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    cat("✅ Human conversion available (org.Hs.eg.db)\n")
  } else {
    cat("❌ Human conversion NOT available - install org.Hs.eg.db\n")
  }
  
  # Test mouse packages
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE) && 
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    cat("✅ Mouse conversion available (org.Mm.eg.db)\n")
  } else {
    cat("❌ Mouse conversion NOT available - install org.Mm.eg.db\n")
  }
}

# Auto-detect species from gene IDs
detect_species_from_genes <- function(gene_ids) {
  human_count <- sum(grepl("^ENSG[0-9]", gene_ids))
  mouse_count <- sum(grepl("^ENSMUSG[0-9]", gene_ids))
  
  if (mouse_count > human_count) {
    return("mouse")
  } else if (human_count > 0) {
    return("human") 
  } else {
    return("human")  # Default
  }
}

# Install required packages if missing
install_gene_conversion_packages <- function() {
  cat("📦 Installing required gene conversion packages...\n")
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager")
  }
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    cat("Installing org.Hs.eg.db (human gene annotations)...\n")
    BiocManager::install("org.Hs.eg.db")
  }
  
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    cat("Installing org.Mm.eg.db (mouse gene annotations)...\n")
    BiocManager::install("org.Mm.eg.db")
  }
  
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    cat("Installing AnnotationDbi...\n")
    BiocManager::install("AnnotationDbi")
  }
  
  cat("✅ Package installation completed!\n")
}

cat("✅ Simple gene conversion functions loaded\n")
cat("📋 Available functions: convert_genes_simple(), test_gene_conversion_packages(), detect_species_from_genes(), install_gene_conversion_packages()\n")

# Test packages on load
test_gene_conversion_packages()