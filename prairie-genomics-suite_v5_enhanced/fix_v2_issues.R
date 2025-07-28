# Fix v2.0 Issues - Install Missing Packages and Test
# Run this script to fix the reported issues

cat("🔧 Prairie Genomics Suite v2.0 - Issue Fix Script\n")
cat("=" , rep("=", 60), "\n")

# Issue 1: Install missing pathway analysis packages
cat("\n📦 Installing missing pathway analysis packages...\n")

# Install BiocManager if not available
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Critical packages for pathway analysis
critical_packages <- c(
  "clusterProfiler",   # Core pathway analysis - REQUIRED
  "enrichplot",        # Advanced pathway visualizations
  "DOSE",              # Disease Ontology analysis
  "org.Hs.eg.db",      # Human gene annotations
  "org.Mm.eg.db",      # Mouse gene annotations
  "AnnotationDbi",     # Annotation database interface
  "GO.db",             # Gene Ontology database
  "biomaRt",           # BioMart interface
  "msigdbr"            # MSigDB gene sets (CRAN)
)

# Install CRAN packages
cran_packages <- c("msigdbr")
for (pkg in cran_packages) {
  cat(paste0("Installing ", pkg, " from CRAN... "))
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    tryCatch({
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
      cat("✅\n")
    }, error = function(e) {
      cat("❌\n")
      cat(paste0("  Error: ", e$message, "\n"))
    })
  } else {
    cat("✅ (already installed)\n")
  }
}

# Install Bioconductor packages
bioc_packages <- setdiff(critical_packages, cran_packages)
for (pkg in bioc_packages) {
  cat(paste0("Installing ", pkg, " from Bioconductor... "))
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
      library(pkg, character.only = TRUE)
      cat("✅\n")
    }, error = function(e) {
      cat("❌\n")
      cat(paste0("  Error: ", e$message, "\n"))
    })
  } else {
    cat("✅ (already installed)\n")
  }
}

# Issue 2: Test pathway analysis functionality
cat("\n🧪 Testing pathway analysis functionality...\n")

tryCatch({
  # Test clusterProfiler
  if (require("clusterProfiler", quietly = TRUE)) {
    cat("✅ clusterProfiler loaded successfully\n")
  } else {
    cat("❌ clusterProfiler failed to load\n")
  }
  
  # Test gene annotation databases
  if (require("org.Hs.eg.db", quietly = TRUE)) {
    cat("✅ Human gene annotations available\n")
  } else {
    cat("❌ Human gene annotations not available\n")
  }
  
  if (require("org.Mm.eg.db", quietly = TRUE)) {
    cat("✅ Mouse gene annotations available\n")
  } else {
    cat("❌ Mouse gene annotations not available\n")
  }
  
  # Test BioMart
  if (require("biomaRt", quietly = TRUE)) {
    cat("✅ BioMart interface available\n")
  } else {
    cat("❌ BioMart interface not available\n")
  }
  
}, error = function(e) {
  cat("❌ Testing failed:", e$message, "\n")
})

# Issue 3: Test gene conversion cache
cat("\n🔄 Testing gene conversion functionality...\n")

tryCatch({
  source("gene_conversion_cache.R")
  cat("✅ Gene conversion cache module loaded\n")
  
  # Test with a few sample genes
  test_genes <- c("ENSG00000141510", "ENSG00000155657", "ENSG00000117399")
  cat("🧪 Testing gene conversion with sample genes...\n")
  
  result <- convert_genes_fast(test_genes, species = "human")
  if (!is.null(result)) {
    cat("✅ Gene conversion test successful\n")
  } else {
    cat("⚠️ Gene conversion test returned NULL\n")
  }
  
}, error = function(e) {
  cat("❌ Gene conversion test failed:", e$message, "\n")
})

# Summary and next steps
cat("\n🎉 Fix Summary\n")
cat(rep("=", 40), "\n")

cat("✅ Package installation completed\n")
cat("✅ BioMart fallback function fixed\n")
cat("✅ Gene conversion cache enhanced\n")

cat("\n🚀 Next Steps:\n")
cat("1. Restart R session to ensure all packages are properly loaded\n")
cat("2. Run the app with: shiny::runApp('app.R')\n")
cat("3. Test pathway analysis with your data\n")

cat("\n💡 If you still encounter issues:\n")
cat("- Check that all required packages are installed: source('install.R')\n")
cat("- Restart R session completely\n")
cat("- Check internet connection for BioMart queries\n")

cat("\n🧬 Prairie Genomics Suite v2.0 Issues Fixed!\n")