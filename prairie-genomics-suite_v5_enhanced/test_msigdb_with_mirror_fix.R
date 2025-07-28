# Test MSigDB with CRAN Mirror Fix
# Fixed version that handles CRAN mirror configuration
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🧪 Testing MSigDB Collections with CRAN Mirror Fix\n")
cat("=" , rep("=", 55), "\n")

# CRITICAL FIX: Set CRAN mirror FIRST
cat("🔧 Setting CRAN Mirror\n")
cat(rep("-", 25), "\n")
options(repos = c(CRAN = "https://cloud.r-project.org/"))
cat("✅ CRAN mirror set to:", getOption("repos")["CRAN"], "\n")

# Test CRAN accessibility
tryCatch({
  available_packages <- available.packages()
  cat("✅ CRAN accessible -", nrow(available_packages), "packages available\n")
}, error = function(e) {
  cat("❌ CRAN access failed:", e$message, "\n")
  cat("🔧 Trying alternative mirror...\n")
  options(repos = c(CRAN = "https://cran.rstudio.com/"))
})

# Test 1: Install and load msigdbr
cat("\nTest 1: MSigDB Package Installation\n")
cat(rep("-", 40), "\n")

msigdbr_loaded <- FALSE

tryCatch({
  if (!require("msigdbr", quietly = TRUE)) {
    cat("📦 Installing msigdbr...\n")
    install.packages("msigdbr", dependencies = TRUE)
    library(msigdbr)
  }
  msigdbr_loaded <- TRUE
  cat("✅ msigdbr loaded successfully\n")
}, error = function(e) {
  cat("❌ msigdbr installation failed:", e$message, "\n")
  cat("💡 Try manually: install.packages('msigdbr')\n")
})

# Test 2: Check species availability (only if msigdbr loaded)
if (msigdbr_loaded) {
  cat("\nTest 2: Species Availability\n")
  cat(rep("-", 30), "\n")
  
  tryCatch({
    species_list <- msigdbr_species()
    cat("🌍 Available species:\n")
    
    # Show key species
    key_species <- c("Homo sapiens", "Mus musculus")
    for (species in key_species) {
      available <- species %in% species_list$species_name
      status <- ifelse(available, "✅", "❌")
      cat("   ", status, species, "\n")
    }
    
  }, error = function(e) {
    cat("❌ Species check failed:", e$message, "\n")
  })
  
  # Test 3: Quick collection test
  cat("\nTest 3: Quick Collection Test\n") 
  cat(rep("-", 35), "\n")
  
  # Test Hallmark for mouse (most reliable)
  tryCatch({
    cat("Testing Hallmark (H) for mouse... ")
    h_mouse <- msigdbr(species = "Mus musculus", category = "H")
    cat("✅", nrow(h_mouse), "pathways\n")
    
    # Show a few examples
    if (nrow(h_mouse) > 0) {
      examples <- unique(h_mouse$gs_name)[1:3]
      cat("   Examples:", paste(examples, collapse = ", "), "\n")
    }
    
  }, error = function(e) {
    cat("❌", e$message, "\n")
  })
  
  # Test C2:CP:KEGG for mouse  
  tryCatch({
    cat("Testing C2:CP:KEGG for mouse... ")
    kegg_mouse <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
    cat("✅", nrow(kegg_mouse), "pathways\n")
    
  }, error = function(e) {
    cat("❌", e$message, "\n")
  })
  
  # Test 4: Test improved get_fgsea_gene_sets function
  cat("\nTest 4: Updated Gene Set Function\n")
  cat(rep("-", 40), "\n")
  
  # Create a simple version of the function for testing
  get_fgsea_gene_sets_test <- function(species, collection = "H") {
    tryCatch({
      # Map species
      msigdb_species <- if (species == "human") {
        "Homo sapiens"
      } else if (species == "mouse") {
        "Mus musculus"
      } else {
        "Homo sapiens"
      }
      
      cat("📚 Retrieving", collection, "for", msigdb_species, "... ")
      
      # Try primary collection
      gene_sets <- msigdbr(species = msigdb_species, category = collection)
      
      # Fallback for C2
      if (nrow(gene_sets) == 0 && collection == "C2") {
        cat("fallback to C2:CP... ")
        gene_sets <- msigdbr(species = msigdb_species, category = "C2", subcategory = "CP")
      }
      
      if (nrow(gene_sets) == 0) {
        cat("❌ No gene sets found\n")
        return(NULL)
      }
      
      # Convert to named list
      pathway_names <- unique(gene_sets$gs_name)
      pathways <- list()
      
      for (pathway in pathway_names) {
        genes <- gene_sets$gene_symbol[gene_sets$gs_name == pathway]
        pathways[[pathway]] <- genes
      }
      
      cat("✅", length(pathways), "pathways\n")
      return(pathways)
      
    }, error = function(e) {
      cat("❌", e$message, "\n")
      return(NULL)
    })
  }
  
  # Test different scenarios
  test_cases <- list(
    list(species = "mouse", collection = "H", name = "Mouse Hallmark"),
    list(species = "mouse", collection = "C2", name = "Mouse C2 (with fallback)"),
    list(species = "human", collection = "H", name = "Human Hallmark")
  )
  
  for (test_case in test_cases) {
    cat("🧪", test_case$name, ": ")
    result <- get_fgsea_gene_sets_test(test_case$species, test_case$collection)
    
    if (!is.null(result) && length(result) > 0) {
      cat("   SUCCESS -", length(result), "pathways retrieved\n")
    } else {
      cat("   FAILED\n")
    }
  }
  
} else {
  cat("\n⚠️ Skipping further tests - msigdbr not available\n")
}

# Summary and next steps
cat("\n🎯 Test Summary\n")
cat("=" , rep("=", 20), "\n")

cat("CRAN Mirror Fix:\n")
if (getOption("repos")["CRAN"] != "@CRAN@") {
  cat("✅ CRAN mirror properly configured\n")
} else {
  cat("❌ CRAN mirror still not set\n")
}

cat("\nmsigdbr Package:\n")
if (msigdbr_loaded) {
  cat("✅ msigdbr successfully loaded\n")
  cat("✅ Mouse collections should now work\n")
} else {
  cat("❌ msigdbr installation failed\n")
  cat("💡 Manual fix: options(repos = c(CRAN = 'https://cloud.r-project.org/'))\n")
  cat("💡 Then run: install.packages('msigdbr')\n")
}

cat("\n🚀 Next Steps:\n")
if (msigdbr_loaded) {
  cat("1. ✅ CRAN mirror is working\n")
  cat("2. ✅ msigdbr is available\n") 
  cat("3. 🧪 Test the main app with mouse C2 collection\n")
  cat("4. 📊 Mouse pathway analysis should now work\n")
} else {
  cat("1. 🔧 Fix CRAN mirror: options(repos = c(CRAN = 'https://cloud.r-project.org/'))\n")
  cat("2. 📦 Install msigdbr: install.packages('msigdbr')\n")
  cat("3. 🔄 Re-run this test script\n")
  cat("4. 🧪 Then test the main app\n")
}

cat("\n🧬 MSigDB Testing with Mirror Fix Complete!\n")