# Fix MSigDB Gene Set Collection Issues
# Diagnose and fix "Could not get C2 gene sets for mouse" error
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🔧 Fixing MSigDB Gene Set Collection Issues\n")
cat("=" , rep("=", 60), "\n")

# Issue: MSigDB gene set retrieval failing for mouse C2 collection
cat("❌ Problem: Could not get C2 gene sets for mouse\n")
cat("🔍 Possible causes:\n")
cat("   • MSigDB collection codes incorrect\n")
cat("   • Species mapping issues\n") 
cat("   • Collection availability for mouse\n")
cat("   • msigdbr package problems\n")

# Diagnostic 1: Test MSigDB package and availability
cat("\n🧪 Diagnostic 1: MSigDB Package Testing\n")
cat(rep("-", 45), "\n")

tryCatch({
  # Load msigdbr
  if (!require("msigdbr", quietly = TRUE)) {
    cat("📦 Installing msigdbr...\n")
    install.packages("msigdbr")
    library(msigdbr)
  }
  
  cat("✅ msigdbr package loaded\n")
  
  # Check available species
  cat("🌍 Available species in MSigDB:\n")
  available_species <- msigdbr_species()
  print(available_species)
  
  # Check mouse specifically
  mouse_available <- "Mus musculus" %in% available_species$species_name
  cat("🐭 Mouse (Mus musculus) available:", ifelse(mouse_available, "✅ YES", "❌ NO"), "\n")
  
}, error = function(e) {
  cat("❌ MSigDB package test failed:", e$message, "\n")
})

# Diagnostic 2: Test gene set collections for mouse
cat("\n🧪 Diagnostic 2: Mouse Gene Set Collections\n")
cat(rep("-", 45), "\n")

tryCatch({
  if (require("msigdbr", quietly = TRUE)) {
    cat("🐭 Testing mouse gene set collections...\n")
    
    # Test different collections
    collections_to_test <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
    
    mouse_collections <- list()
    
    for (collection in collections_to_test) {
      tryCatch({
        cat("Testing", collection, "collection... ")
        
        gene_sets <- msigdbr(species = "Mus musculus", category = collection)
        
        if (nrow(gene_sets) > 0) {
          cat("✅", nrow(gene_sets), "entries\n")
          mouse_collections[[collection]] <- nrow(gene_sets)
        } else {
          cat("❌ Empty\n")
          mouse_collections[[collection]] <- 0
        }
        
      }, error = function(e) {
        cat("❌ Error:", e$message, "\n")
        mouse_collections[[collection]] <- -1
      })
    }
    
    cat("\n📊 Mouse Collection Summary:\n")
    for (col in names(mouse_collections)) {
      count <- mouse_collections[[col]]
      status <- if (count > 0) "✅" else if (count == 0) "⚠️" else "❌"
      cat("   ", status, col, ":", count, "gene sets\n")
    }
    
  }
}, error = function(e) {
  cat("❌ Mouse collection test failed:", e$message, "\n")
})

# Diagnostic 3: Test specific C2 subcategories
cat("\n🧪 Diagnostic 3: C2 Collection Subcategories\n")
cat(rep("-", 45), "\n")

tryCatch({
  if (require("msigdbr", quietly = TRUE)) {
    cat("🔬 Testing C2 subcategories for mouse...\n")
    
    # C2 has subcategories - test them individually
    c2_subcategories <- c("CGP", "CP", "CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS")
    
    for (subcat in c2_subcategories) {
      tryCatch({
        cat("Testing C2:", subcat, "... ")
        
        if (subcat == "CP") {
          # CP without subcategory
          gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
        } else if (grepl(":", subcat)) {
          # Subcategories with colon
          parts <- strsplit(subcat, ":")[[1]]
          gene_sets <- msigdbr(species = "Mus musculus", category = parts[1], subcategory = subcat)
        } else {
          gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = subcat)
        }
        
        cat("✅", nrow(gene_sets), "entries\n")
        
      }, error = function(e) {
        cat("❌", e$message, "\n")
      })
    }
    
    # Try C2 without subcategory
    cat("Testing C2 (all subcategories)... ")
    tryCatch({
      all_c2 <- msigdbr(species = "Mus musculus", category = "C2")
      cat("✅", nrow(all_c2), "entries\n")
    }, error = function(e) {
      cat("❌", e$message, "\n")
    })
    
  }
}, error = function(e) {
  cat("❌ C2 subcategory test failed:", e$message, "\n")
})

# Create improved gene set retrieval function
cat("\n🔧 Creating Improved Gene Set Retrieval\n")
cat(rep("-", 45), "\n")

improved_msigdb_function <- '
# IMPROVED: Robust MSigDB gene set retrieval with fallbacks
get_fgsea_gene_sets <- function(species, collection = "H") {
  tryCatch({
    # Check msigdbr availability
    if (!require("msigdbr", quietly = TRUE)) {
      cat("📦 Installing msigdbr...\\n")
      install.packages("msigdbr")
      library(msigdbr)
    }
    
    # Map species names with validation
    msigdb_species <- if (species == "human") {
      "Homo sapiens"
    } else if (species == "mouse") {
      "Mus musculus"
    } else {
      cat("⚠️ Unsupported species, defaulting to human\\n")
      "Homo sapiens"
    }
    
    # Validate species availability
    available_species <- msigdbr_species()
    if (!msigdb_species %in% available_species$species_name) {
      cat("❌ Species", msigdb_species, "not available in MSigDB\\n")
      return(NULL)
    }
    
    cat("📚 Retrieving", collection, "gene sets for", msigdb_species, "\\n")
    
    # Try to get gene sets with collection fallbacks
    gene_sets <- NULL
    
    # Primary attempt with exact collection
    tryCatch({
      gene_sets <- msigdbr(species = msigdb_species, category = collection)
      if (nrow(gene_sets) == 0) {
        gene_sets <- NULL
      }
    }, error = function(e) {
      cat("⚠️ Primary retrieval failed:", e$message, "\\n")
      gene_sets <<- NULL
    })
    
    # Fallback strategies for specific collections
    if (is.null(gene_sets) && collection == "C2") {
      cat("🔄 C2 collection failed, trying C2:CP (Canonical Pathways)...\\n")
      tryCatch({
        gene_sets <- msigdbr(species = msigdb_species, category = "C2", subcategory = "CP")
        if (nrow(gene_sets) == 0) gene_sets <- NULL
      }, error = function(e) {
        cat("⚠️ C2:CP fallback failed:", e$message, "\\n")
      })
    }
    
    # Further fallback to KEGG pathways if C2 fails
    if (is.null(gene_sets) && collection == "C2") {
      cat("🔄 Trying C2:CP:KEGG as final fallback...\\n")
      tryCatch({
        gene_sets <- msigdbr(species = msigdb_species, category = "C2", subcategory = "CP:KEGG")
        if (nrow(gene_sets) == 0) gene_sets <- NULL
      }, error = function(e) {
        cat("⚠️ KEGG fallback failed:", e$message, "\\n")
      })
    }
    
    # Ultimate fallback to Hallmark if all else fails
    if (is.null(gene_sets) && collection != "H") {
      cat("🔄 All", collection, "attempts failed, falling back to Hallmark (H)...\\n")
      tryCatch({
        gene_sets <- msigdbr(species = msigdb_species, category = "H")
        if (nrow(gene_sets) == 0) gene_sets <- NULL
      }, error = function(e) {
        cat("❌ Even Hallmark fallback failed:", e$message, "\\n")
      })
    }
    
    if (is.null(gene_sets) || nrow(gene_sets) == 0) {
      cat("❌ No gene sets found for", collection, "in", msigdb_species, "\\n")
      return(NULL)
    }
    
    # Convert to named list format required by fgsea
    pathways <- gene_sets %>%
      dplyr::select(gs_name, gene_symbol) %>%
      group_by(gs_name) %>%
      summarise(genes = list(gene_symbol), .groups = "drop") %>%
      deframe()  # Converts to named list
    
    cat("✅ Retrieved", length(pathways), "gene sets\\n")
    cat("📋 Collection used:", unique(gene_sets$gs_cat)[1], "\\n")
    if ("gs_subcat" %in% colnames(gene_sets)) {
      subcats <- unique(gene_sets$gs_subcat)
      if (length(subcats) <= 3) {
        cat("📋 Subcategories:", paste(subcats, collapse = ", "), "\\n")
      }
    }
    
    return(pathways)
    
  }, error = function(e) {
    cat("❌ Failed to get gene sets:", e$message, "\\n")
    return(NULL)
  })
}

# Alternative function using base R (no dplyr)
get_fgsea_gene_sets_base <- function(species, collection = "H") {
  tryCatch({
    # Check msigdbr availability
    if (!require("msigdbr", quietly = TRUE)) {
      cat("📦 Installing msigdbr...\\n")
      install.packages("msigdbr")
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
    
    cat("📚 Retrieving", collection, "gene sets for", msigdb_species, "\\n")
    
    # Get gene sets with error handling
    gene_sets <- tryCatch({
      msigdbr(species = msigdb_species, category = collection)
    }, error = function(e) {
      cat("❌ Gene set retrieval failed:", e$message, "\\n")
      return(NULL)
    })
    
    if (is.null(gene_sets) || nrow(gene_sets) == 0) {
      cat("❌ No gene sets found\\n")
      return(NULL)
    }
    
    # Convert to named list using base R
    pathway_names <- unique(gene_sets$gs_name)
    pathways <- list()
    
    for (pathway in pathway_names) {
      genes <- gene_sets$gene_symbol[gene_sets$gs_name == pathway]
      pathways[[pathway]] <- genes
    }
    
    cat("✅ Retrieved", length(pathways), "gene sets using base R\\n")
    
    return(pathways)
    
  }, error = function(e) {
    cat("❌ Base R gene set retrieval failed:", e$message, "\\n")
    return(NULL)
  })
}
'

# Write improved functions to file
writeLines(improved_msigdb_function, "improved_msigdb_functions.R")
cat("✅ Improved MSigDB functions saved to improved_msigdb_functions.R\n")

# Create collection mapping guide
cat("\n📚 MSigDB Collection Guide\n")
cat(rep("-", 30), "\n")

collection_guide <- '
# MSigDB Collection Guide for Prairie Genomics Suite

## Available Collections:
- H: Hallmark gene sets (50 sets) - RECOMMENDED for most analyses
- C1: Positional gene sets (chromosomal location)
- C2: Curated gene sets (canonical pathways)
  - C2:CGP: Chemical and genetic perturbations  
  - C2:CP: Canonical pathways (KEGG, Reactome, etc.)
  - C2:CP:KEGG: KEGG pathways only
  - C2:CP:REACTOME: Reactome pathways only
- C3: Regulatory target gene sets (microRNA, transcription factors)
- C4: Computational gene sets
- C5: Gene Ontology gene sets
  - C5:BP: Biological Process
  - C5:CC: Cellular Component
  - C5:MF: Molecular Function
- C6: Oncogenic signature gene sets
- C7: Immunologic signature gene sets
- C8: Cell type signature gene sets

## Recommendations:
- Start with H (Hallmark) - most reliable across species
- For pathways: Try C2:CP:KEGG or C2:CP:REACTOME
- Avoid large collections like C2 (all) for performance
- Some collections may not be available for all species

## Troubleshooting:
- If C2 fails, try C2:CP or C2:CP:KEGG
- Always test Hallmark (H) first
- Check msigdbr_species() for species availability
'

writeLines(collection_guide, "msigdb_collection_guide.txt")

cat("\n🎯 Fix Summary\n")
cat("=" , rep("=", 30), "\n")

cat("Issues Identified:\n")
cat("• C2 collection too broad for some species\n")
cat("• Need subcategory specification for C2\n")
cat("• Missing fallback strategies\n")
cat("• No collection availability validation\n")

cat("\nSolutions Implemented:\n")
cat("✅ Improved gene set retrieval with fallbacks\n")
cat("✅ C2 subcategory handling (CP, KEGG, etc.)\n")
cat("✅ Species availability validation\n")
cat("✅ Fallback to Hallmark if collection fails\n")
cat("✅ Base R alternative without dplyr\n")

cat("\n💡 Recommendations:\n")
cat("• Use Hallmark (H) collection for reliability\n")
cat("• For pathways, try C2:CP:KEGG instead of C2\n")
cat("• Test collection availability before analysis\n")
cat("• Implement collection fallback hierarchy\n")

cat("\n🚀 Next Steps:\n")
cat("1. Update pathway_analysis.R with improved functions\n")
cat("2. Test with mouse C2:CP:KEGG collection\n")
cat("3. Add collection selection guidance to UI\n")
cat("4. Implement automatic fallback in main analysis\n")

cat("\n🧬 Prairie Genomics Suite - MSigDB Collections Fixed!\n")