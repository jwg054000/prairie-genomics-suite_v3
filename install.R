# Install Required Packages for Prairie Genomics Suite - R Shiny
# Run this script to install all necessary R packages

cat("🧬 Installing Prairie Genomics Suite R Shiny Dependencies\n")
cat("=" , rep("=", 60), "\n")

# List of required packages
required_packages <- c(
  # Shiny framework
  "shiny",
  "shinydashboard", 
  "shinyWidgets",
  "shinyjs",
  
  # Data manipulation
  "dplyr",
  "readr",
  "readxl",
  "DT",
  
  # Visualization
  "ggplot2",
  "plotly",
  "RColorBrewer",
  "pheatmap",
  "ggrepel",
  "viridis",
  
  # Statistical analysis
  "DESeq2",
  
  # Bioconductor dependencies
  "BiocManager",
  
  # Utilities
  "stringr",
  "future",
  "promises"
)

# Bioconductor packages
bioc_packages <- c(
  "DESeq2",
  "EnhancedVolcano"
)

# Function to install packages
install_packages <- function(packages, source = "CRAN") {
  cat(paste0("\n📦 Installing ", source, " packages...\n"))
  
  for (pkg in packages) {
    cat(paste0("  - ", pkg, "... "))
    
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      tryCatch({
        if (source == "CRAN") {
          install.packages(pkg, dependencies = TRUE, quiet = TRUE)
        } else if (source == "Bioconductor") {
          BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
        }
        
        # Test if package loads
        library(pkg, character.only = TRUE, quietly = TRUE)
        cat("✅\n")
        
      }, error = function(e) {
        cat("❌\n")
        cat(paste0("    Error: ", e$message, "\n"))
      })
    } else {
      cat("✅ (already installed)\n")
    }
  }
}

# Install BiocManager first if needed
if (!require("BiocManager", quietly = TRUE)) {
  cat("📦 Installing BiocManager...\n")
  install.packages("BiocManager")
}

# Install CRAN packages
cran_packages <- setdiff(required_packages, bioc_packages)
install_packages(cran_packages, "CRAN")

# Install Bioconductor packages
install_packages(bioc_packages, "Bioconductor")

# Test installation
cat("\n🧪 Testing Installation...\n")
cat(rep("=", 40), "\n")

all_packages <- c(required_packages, bioc_packages)
failed_packages <- character(0)

for (pkg in all_packages) {
  cat(paste0("Testing ", pkg, "... "))
  
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✅\n")
  } else {
    cat("❌\n")
    failed_packages <- c(failed_packages, pkg)
  }
}

# Summary
cat("\n🎉 Installation Summary\n")
cat(rep("=", 40), "\n")

if (length(failed_packages) == 0) {
  cat("✅ All packages installed successfully!\n")
  cat("\n🚀 You can now run the Prairie Genomics Suite with:\n")
  cat("   R -e \"shiny::runApp('app.R')\"\n")
  cat("\n📖 Or in RStudio:\n")
  cat("   1. Open app.R\n")
  cat("   2. Click 'Run App' button\n")
} else {
  cat("❌ Failed to install the following packages:\n")
  for (pkg in failed_packages) {
    cat(paste0("   - ", pkg, "\n"))
  }
  cat("\n💡 Try installing these manually:\n")
  cat("   install.packages(c(", paste0("'", failed_packages, "'", collapse = ", "), "))\n")
}

cat("\n📋 Package Versions:\n")
for (pkg in setdiff(all_packages, failed_packages)) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    version <- packageVersion(pkg)
    cat(paste0("   ", pkg, ": ", version, "\n"))
  }
}

cat("\n🧬 Prairie Genomics Suite R Shiny Setup Complete!\n")