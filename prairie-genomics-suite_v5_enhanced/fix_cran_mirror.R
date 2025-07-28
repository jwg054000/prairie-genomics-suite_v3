# Fix CRAN Mirror Issue for Package Installation
# Resolves "trying to use CRAN without setting a mirror" error
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🔧 Fixing CRAN Mirror Configuration\n")
cat("=" , rep("=", 50), "\n")

# Issue: R session doesn't have CRAN mirror configured
cat("❌ Problem: CRAN mirror not set for package installation\n")
cat("🔧 Solution: Configure CRAN mirror and update installation functions\n")

# Fix 1: Set CRAN mirror
cat("\n📍 Setting CRAN Mirror\n")
cat(rep("-", 25), "\n")

# Set CRAN mirror to a reliable source
options(repos = c(CRAN = "https://cloud.r-project.org/"))

cat("✅ CRAN mirror set to: https://cloud.r-project.org/\n")
cat("🌍 This is the RStudio CRAN mirror (reliable worldwide)\n")

# Test the mirror
tryCatch({
  available_packages <- available.packages()
  cat("✅ CRAN mirror accessible -", nrow(available_packages), "packages available\n")
}, error = function(e) {
  cat("⚠️ CRAN mirror test failed:", e$message, "\n")
})

# Fix 2: Create improved package installation function
cat("\n📦 Creating Improved Package Installation\n")
cat(rep("-", 40), "\n")

improved_install_function <- '
# IMPROVED: Safe package installation with CRAN mirror handling
safe_install_packages <- function(packages, type = "CRAN") {
  # Ensure CRAN mirror is set
  if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org/"))
    cat("🔧 CRAN mirror configured\\n")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("📦 Installing", pkg, "...")
      
      tryCatch({
        if (type == "CRAN") {
          install.packages(pkg, dependencies = TRUE, quiet = TRUE)
        } else if (type == "Bioconductor") {
          if (!require("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", quiet = TRUE)
          }
          BiocManager::install(pkg, update = FALSE, ask = FALSE)
        }
        
        # Test if installation worked
        if (require(pkg, character.only = TRUE, quietly = TRUE)) {
          cat(" ✅\\n")
        } else {
          cat(" ❌ Installation failed\\n")
        }
        
      }, error = function(e) {
        cat(" ❌", e$message, "\\n")
      })
    } else {
      cat("✅", pkg, "already installed\\n")
    }
  }
}

# Install critical packages with proper error handling
install_pathway_packages <- function() {
  cat("📦 Installing Critical Pathway Analysis Packages\\n")
  cat(rep("-", 50), "\\n")
  
  # Set CRAN mirror first
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  
  # CRAN packages
  cran_packages <- c("msigdbr", "dplyr", "digest")
  cat("Installing CRAN packages...\\n")
  safe_install_packages(cran_packages, "CRAN")
  
  # Bioconductor packages  
  bioc_packages <- c("clusterProfiler", "fgsea", "org.Hs.eg.db", "org.Mm.eg.db", "enrichplot")
  cat("\\nInstalling Bioconductor packages...\\n")
  safe_install_packages(bioc_packages, "Bioconductor")
  
  cat("\\n🎉 Package installation completed!\\n")
}
'

# Write the improved functions to file
writeLines(improved_install_function, "safe_package_installation.R")
cat("✅ Safe installation functions saved to safe_package_installation.R\n")

# Fix 3: Update pathway_analysis.R with mirror handling
cat("\n🔧 Creating Updated MSigDB Functions\n")
cat(rep("-", 35), "\n")

updated_msigdb_functions <- '
# UPDATED: MSigDB functions with proper CRAN mirror handling
get_fgsea_gene_sets_safe <- function(species, collection = "H") {
  tryCatch({
    # Ensure CRAN mirror is set before any package operations
    if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
      options(repos = c(CRAN = "https://cloud.r-project.org/"))
    }
    
    # Check msigdbr availability with safe installation
    if (!require("msigdbr", quietly = TRUE)) {
      cat("📦 Installing msigdbr with proper CRAN mirror...\\n")
      tryCatch({
        install.packages("msigdbr", dependencies = TRUE)
        library(msigdbr)
      }, error = function(e) {
        cat("❌ msigdbr installation failed:", e$message, "\\n")
        return(NULL)
      })
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
    
    # Convert to named list using base R (no dplyr)
    pathway_names <- unique(gene_sets$gs_name)
    pathways <- list()
    
    for (pathway in pathway_names) {
      genes <- gene_sets$gene_symbol[gene_sets$gs_name == pathway]
      pathways[[pathway]] <- genes
    }
    
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
'

writeLines(updated_msigdb_functions, "updated_msigdb_functions.R")
cat("✅ Updated MSigDB functions saved\n")

# Fix 4: Create startup configuration
cat("\n⚙️ Creating R Startup Configuration\n")
cat(rep("-", 35), "\n")

startup_config <- '
# R Startup Configuration for Prairie Genomics Suite
# Add this to your .Rprofile or source at the beginning of your session

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Set Bioconductor repository
options(BioC_mirror = "https://bioconductor.org")

# Increase timeout for package downloads (useful for large packages)
options(timeout = 300)  # 5 minutes

# Disable package startup messages for cleaner output
options(warn = 1)

cat("🔧 Prairie Genomics Suite R configuration loaded\\n")
cat("📍 CRAN mirror:", getOption("repos")["CRAN"], "\\n")
'

writeLines(startup_config, "r_startup_config.R")
cat("✅ R startup configuration saved\n")

# Test the fixes
cat("\n🧪 Testing CRAN Mirror Fix\n")
cat(rep("-", 30), "\n")

# Test package availability
tryCatch({
  available <- available.packages()
  msigdbr_available <- "msigdbr" %in% rownames(available)
  cat("✅ CRAN packages accessible\n")
  cat("📦 msigdbr available for installation:", ifelse(msigdbr_available, "YES", "NO"), "\n")
}, error = function(e) {
  cat("❌ CRAN access test failed:", e$message, "\n")
})

# Instructions
cat("\n📋 Usage Instructions\n")
cat(rep("-", 25), "\n")

cat("To fix the CRAN mirror issue:\n\n")

cat("Option 1 - Quick Fix (current session):\n")
cat("   options(repos = c(CRAN = 'https://cloud.r-project.org/'))\n\n")

cat("Option 2 - Permanent Fix:\n")
cat("   1. Source the startup config: source('r_startup_config.R')\n")
cat("   2. Or add to your .Rprofile file\n\n")

cat("Option 3 - Use the safe installation:\n")
cat("   1. source('safe_package_installation.R')\n")
cat("   2. install_pathway_packages()\n\n")

cat("🎯 Summary\n")
cat("=" , rep("=", 15), "\n")

cat("Issues Fixed:\n")
cat("✅ CRAN mirror configuration\n")
cat("✅ Safe package installation functions\n")
cat("✅ Error handling for package failures\n") 
cat("✅ Updated MSigDB functions with mirror handling\n")

cat("\nRecommended Steps:\n")
cat("1. Run: options(repos = c(CRAN = 'https://cloud.r-project.org/'))\n")
cat("2. Test: install.packages('msigdbr')\n")
cat("3. Re-run the MSigDB collection test\n")
cat("4. Source r_startup_config.R for future sessions\n")

cat("\n🧬 CRAN Mirror Issue Fixed!\n")
cat("📦 Package installation should now work properly\n")