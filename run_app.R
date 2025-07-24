# Prairie Genomics Suite - Startup Script
# Simplified launcher that checks dependencies and provides helpful guidance

cat("🧬 Prairie Genomics Suite R Shiny Launcher\n")
cat("=========================================\n")

# Check R version
r_version <- R.version.string
cat(paste0("R Version: ", r_version, "\n"))

if (getRversion() < "4.0.0") {
  cat("⚠️  Warning: R version 4.0.0 or higher recommended\n")
}

# Check core packages
core_packages <- c("shiny", "shinydashboard", "DT", "ggplot2", "plotly", "dplyr")
missing_core <- c()

cat("\n📦 Checking core packages...\n")
for (pkg in core_packages) {
  cat(paste0("  ", pkg, "... "))
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✅\n")
  } else {
    cat("❌\n")
    missing_core <- c(missing_core, pkg)
  }
}

# Check optional packages  
optional_packages <- c("DESeq2", "readxl", "RColorBrewer", "pheatmap")
missing_optional <- c()

cat("\n📦 Checking optional packages...\n")
for (pkg in optional_packages) {
  cat(paste0("  ", pkg, "... "))
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✅\n")
  } else {
    cat("⚠️\n")
    missing_optional <- c(missing_optional, pkg)
  }
}

# Report results
if (length(missing_core) > 0) {
  cat("\n❌ CANNOT START - Missing core packages:\n")
  for (pkg in missing_core) {
    cat(paste0("   - ", pkg, "\n"))
  }
  cat("\n💡 Install with:\n")
  cat(paste0('   install.packages(c("', paste(missing_core, collapse = '", "'), '"))\n'))
  cat("\n🔧 Or run the install script:\n")
  cat("   source('install.R')\n")
  quit(status = 1)
}

if (length(missing_optional) > 0) {
  cat("\n⚠️  Some optional packages missing (features will be limited):\n")
  for (pkg in missing_optional) {
    cat(paste0("   - ", pkg, "\n"))
  }
  cat("\n💡 For full functionality, install with:\n")
  cat(paste0('   install.packages(c("', paste(missing_optional, collapse = '", "'), '"))\n'))
  
  if ("DESeq2" %in% missing_optional) {
    cat("\n🧬 For DESeq2 (Bioconductor package):\n")
    cat('   if (!require("BiocManager")) install.packages("BiocManager")\n')
    cat('   BiocManager::install("DESeq2")\n')
  }
}

cat("\n✅ Core packages available - starting application...\n")
cat("🌐 Your app will open in your default web browser\n")
cat("🛑 Press Ctrl+C in this terminal to stop the app\n\n")

# Check if app.R exists
if (!file.exists("app.R")) {
  cat("❌ Error: app.R not found in current directory\n")
  cat("💡 Make sure you're in the prairie-genomics-suite_shiny directory\n")
  quit(status = 1)
}

# Start the app
tryCatch({
  shiny::runApp("app.R", launch.browser = TRUE)
}, error = function(e) {
  cat(paste0("\n❌ Error starting app: ", e$message, "\n"))
  cat("\n🔧 Troubleshooting tips:\n")
  cat("1. Make sure all packages are installed\n")
  cat("2. Check that app.R and modules/ directory exist\n")
  cat("3. Try running: source('simple_test.R') first\n")
})