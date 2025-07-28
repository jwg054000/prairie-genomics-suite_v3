# Fix MSigDB API Compatibility Issues
# Check actual msigdbr API and use correct parameters
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🔧 Fixing MSigDB API Compatibility\n")
cat("=" , rep("=", 50), "\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Load and check msigdbr version
if (!require("msigdbr", quietly = TRUE)) {
  cat("📦 Installing latest msigdbr...\n")
  install.packages("msigdbr")
  library(msigdbr)
}

cat("✅ msigdbr loaded\n")
cat("📋 msigdbr version:", packageVersion("msigdbr"), "\n")

# Issue: Check what API is actually available
cat("\n🔍 Checking Available MSigDB API\n")
cat(rep("-", 35), "\n")

# Check available collections
tryCatch({
  collections <- msigdbr_collections()
  cat("✅ Available collections:\n")
  print(collections)
}, error = function(e) {
  cat("❌ msigdbr_collections() failed:", e$message, "\n")
  cat("💡 Using older API format\n")
})

# Check available species
tryCatch({
  species_info <- msigdbr_species()
  cat("\n✅ Available species:\n")
  mouse_info <- species_info[species_info$species_name == "Mus musculus", ]
  human_info <- species_info[species_info$species_name == "Homo sapiens", ]
  
  cat("🐭 Mouse:\n")
  print(mouse_info)
  cat("🧑 Human:\n") 
  print(human_info)
  
}, error = function(e) {
  cat("❌ msigdbr_species() failed:", e$message, "\n")
})

# Test different API approaches
cat("\n🧪 Testing Different API Approaches\n")
cat(rep("-", 40), "\n")

# Approach 1: Try new API with collection parameter
cat("Approach 1: New API with 'collection' parameter...\n")
tryCatch({
  test1 <- msigdbr(species = "Mus musculus", collection = "H")
  cat("✅ New API works -", nrow(test1), "entries\n")
  api_approach <- "new"
}, error = function(e) {
  cat("❌ New API failed:", e$message, "\n")
  api_approach <- "old"
})

# Approach 2: Try old API with category parameter  
if (!exists("api_approach") || api_approach == "old") {
  cat("Approach 2: Old API with 'category' parameter...\n")
  tryCatch({
    test2 <- msigdbr(species = "Mus musculus", category = "H")
    cat("✅ Old API works -", nrow(test2), "entries\n")
    api_approach <- "old"
  }, error = function(e) {
    cat("❌ Old API also failed:", e$message, "\n")
    api_approach <- "unknown"
  })
}

# Test db_species parameter availability
cat("\nTesting db_species parameter...\n")
tryCatch({
  # Try with db_species parameter
  if (api_approach == "new") {
    test_native <- msigdbr(species = "Mus musculus", collection = "H", db_species = "MM")
  } else {
    test_native <- msigdbr(species = "Mus musculus", category = "H", db_species = "MM")
  }
  cat("✅ db_species parameter works -", nrow(test_native), "entries\n")
  native_available <- TRUE
}, error = function(e) {
  cat("❌ db_species parameter not available:", e$message, "\n")
  native_available <- FALSE
})

# Create compatible function based on what works
cat("\n🔧 Creating Compatible MSigDB Function\n")
cat(rep("-", 40), "\n")

if (exists("api_approach")) {
  
  if (api_approach == "new" && exists("native_available") && native_available) {
    # Best case: New API with native species support
    compatible_function <- '
# COMPATIBLE: MSigDB function with new API and native species
get_fgsea_gene_sets_compatible <- function(species, collection = "H") {
  tryCatch({
    if (!require("msigdbr", quietly = TRUE)) {
      install.packages("msigdbr")
      library(msigdbr)
    }
    
    # Map species with native database specification
    if (species == "human") {
      msigdb_species <- "Homo sapiens"
      db_species <- "HS"  # Native human
    } else if (species == "mouse") {
      msigdb_species <- "Mus musculus" 
      db_species <- "MM"  # NATIVE MOUSE
    } else {
      msigdb_species <- "Homo sapiens"
      db_species <- "HS"
    }
    
    cat("📚 Retrieving NATIVE", collection, "gene sets for", msigdb_species, "\\n")
    
    # Use NEW API with native species
    gene_sets <- msigdbr(
      species = msigdb_species, 
      collection = collection,
      db_species = db_species
    )
    
    if (nrow(gene_sets) == 0) {
      return(NULL)
    }
    
    # Convert to pathways list
    pathway_names <- unique(gene_sets$gs_name)
    pathways <- list()
    
    for (pathway in pathway_names) {
      genes <- gene_sets$gene_symbol[gene_sets$gs_name == pathway]
      pathways[[pathway]] <- genes
    }
    
    cat("✅ Retrieved", length(pathways), "NATIVE pathways\\n")
    cat("📋 Database:", db_species, "(native", species, ")\\n")
    
    return(pathways)
    
  }, error = function(e) {
    cat("❌ Compatible function failed:", e$message, "\\n")
    return(NULL)
  })
}
'
    
  } else if (api_approach == "new") {
    # New API but no native species - use collection parameter
    compatible_function <- '
# COMPATIBLE: MSigDB function with new API (no native species)
get_fgsea_gene_sets_compatible <- function(species, collection = "H") {
  tryCatch({
    if (!require("msigdbr", quietly = TRUE)) {
      install.packages("msigdbr")
      library(msigdbr)
    }
    
    msigdb_species <- if (species == "human") {
      "Homo sapiens"
    } else if (species == "mouse") {
      "Mus musculus"
    } else {
      "Homo sapiens"
    }
    
    cat("📚 Retrieving", collection, "gene sets for", msigdb_species, "\\n")
    
    # Use NEW API without native species
    gene_sets <- msigdbr(species = msigdb_species, collection = collection)
    
    if (nrow(gene_sets) == 0) {
      return(NULL)
    }
    
    # Convert to pathways list
    pathway_names <- unique(gene_sets$gs_name)
    pathways <- list()
    
    for (pathway in pathway_names) {
      genes <- gene_sets$gene_symbol[gene_sets$gs_name == pathway]
      pathways[[pathway]] <- genes
    }
    
    cat("✅ Retrieved", length(pathways), "pathways\\n")
    
    return(pathways)
    
  }, error = function(e) {
    cat("❌ Compatible function failed:", e$message, "\\n")
    return(NULL)
  })
}
'
    
  } else {
    # Fall back to old API with category parameter
    compatible_function <- '
# COMPATIBLE: MSigDB function with old API (category parameter)
get_fgsea_gene_sets_compatible <- function(species, collection = "H") {
  tryCatch({
    if (!require("msigdbr", quietly = TRUE)) {
      install.packages("msigdbr")
      library(msigdbr)
    }
    
    msigdb_species <- if (species == "human") {
      "Homo sapiens"
    } else if (species == "mouse") {
      "Mus musculus"
    } else {
      "Homo sapiens"
    }
    
    cat("📚 Retrieving", collection, "gene sets for", msigdb_species, "\\n")
    
    # Use OLD API with category parameter
    gene_sets <- msigdbr(species = msigdb_species, category = collection)
    
    if (nrow(gene_sets) == 0) {
      return(NULL)
    }
    
    # Convert to pathways list
    pathway_names <- unique(gene_sets$gs_name)
    pathways <- list()
    
    for (pathway in pathway_names) {
      genes <- gene_sets$gene_symbol[gene_sets$gs_name == pathway]
      pathways[[pathway]] <- genes
    }
    
    cat("✅ Retrieved", length(pathways), "pathways\\n")
    
    return(pathways)
    
  }, error = function(e) {
    cat("❌ Compatible function failed:", e$message, "\\n")
    return(NULL)
  })
}
'
  }
  
  cat("✅ Created compatible function for API approach:", api_approach, "\n")
  
  # Write the compatible function
  writeLines(compatible_function, "msigdb_compatible_function.R")
  
  # Test the compatible function
  cat("\n🧪 Testing Compatible Function\n")
  cat(rep("-", 35), "\n")
  
  # Load and test the function
  tryCatch({
    eval(parse(text = compatible_function))
    
    # Test with mouse
    cat("Testing mouse Hallmark gene sets...\n")
    mouse_test <- get_fgsea_gene_sets_compatible("mouse", "H")
    
    if (!is.null(mouse_test) && length(mouse_test) > 0) {
      cat("✅ Mouse test successful:", length(mouse_test), "pathways\n")
      
      # Test with human  
      cat("Testing human Hallmark gene sets...\n")
      human_test <- get_fgsea_gene_sets_compatible("human", "H")
      
      if (!is.null(human_test) && length(human_test) > 0) {
        cat("✅ Human test successful:", length(human_test), "pathways\n")
      }
    } else {
      cat("❌ Compatible function test failed\n")
    }
    
  }, error = function(e) {
    cat("❌ Compatible function testing failed:", e$message, "\n")
  })
  
}

# Summary and recommendations
cat("\n🎯 API Compatibility Summary\n")
cat("=" , rep("=", 35), "\n")

cat("MSigDB Version:", as.character(packageVersion("msigdbr")), "\n")

if (exists("api_approach")) {
  cat("Working API approach:", api_approach, "\n")
  
  if (exists("native_available")) {
    cat("Native species support:", ifelse(native_available, "✅ Available", "❌ Not available"), "\n")
  }
  
  cat("\n💡 Recommendations:\n")
  if (api_approach == "new" && exists("native_available") && native_available) {
    cat("• Use collection parameter with db_species for native gene sets\n")
    cat("• Best performance and accuracy\n")
  } else if (api_approach == "new") {
    cat("• Use collection parameter (no native species)\n")
    cat("• May use ortholog mapping for mouse\n") 
  } else {
    cat("• Use category parameter (older API)\n")
    cat("• May show deprecation warnings but will work\n")
  }
  
  cat("\n🔧 Next Steps:\n")
  cat("1. Update pathway_analysis.R with compatible function\n")
  cat("2. Test with the main application\n")
  cat("3. Monitor for any remaining warnings\n")
  
} else {
  cat("❌ Could not determine working API approach\n")
  cat("💡 Try updating msigdbr: install.packages('msigdbr')\n")
}

cat("\n🧬 MSigDB API Compatibility Check Complete!\n")